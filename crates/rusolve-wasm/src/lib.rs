//! Agent-written MIP hooks as fuel-metered WebAssembly with no host access. Design and ABI:
//! docs/design/wasm-hooks.md.
//!
//! [`WasmHooks`] implements `rusolve_ekk::mip::hooks::MipHooks` by calling exported functions of a
//! Wasm module run by the `wasmi` interpreter (fuel metering: `Config::consume_fuel`, see the
//! wasmi documentation of `Store::set_fuel`). The module gets no imports (no WASI, no host
//! functions, no clock, no randomness), a hard memory cap, and a fixed fuel budget per call, so
//! one call costs the same on every run. A trap, fuel exhaustion, a missing export or an invalid
//! answer makes that single call fall back to the default (HiGHS) behaviour and is counted in
//! [`WasmStats`]; the solve is never interrupted.
//!
//! ABI v1 (all numbers little-endian, as Wasm linear memory):
//! * optional export `memory`: the linear memory (needed by `branch_score` and `select_node`);
//! * optional export `branch_score(depth: i32, n: i32, in_ptr: i32, out_ptr: i32) -> i32`:
//!   the host writes `n` records of four `f64` at `in_ptr` (`col`, `value`, `local_lower`,
//!   `local_upper`; 32 bytes each) and the module writes `n` `f64` scores at `out_ptr`
//!   (`in_ptr = 0`, `out_ptr = 32 * n`; highest score wins). Return 0 for "scores written",
//!   anything else for "no opinion";
//! * optional export `cut_score(violation: f64, num_nonzeros: i32, num_active_nonzeros: i32,
//!   active_norm: f64, default_score: f64) -> f64`: the replacement score of a pool cut; a NaN
//!   or infinite result means "no opinion" (the default score is used);
//! * optional export `select_node(n: i32, in_ptr: i32) -> i32`: the host writes `n` records of
//!   four `f64` at `in_ptr = 0` (`lower_bound`, `estimate`, `depth`, `num_domain_changes`;
//!   32 bytes each) and the module returns the index of the open node to dive into next, in
//!   `[0, n)`, or -1 for "no opinion" (the HiGHS best-bound pop is kept). An out-of-range
//!   answer is invalid. If the `n` records do not fit in the component's memory the call is
//!   declined. Needs an exported `memory`;
//! * optional export `separate(n: i32, in_ptr: i32, out_ptr: i32, out_cap: i32) -> i32`: the user
//!   separator with a variable-size result. The host writes the node LP solution, `n` `f64`
//!   (`in_ptr = 0`), and offers `out_cap` `f64` slots at `out_ptr = 8 * n` (all remaining
//!   memory). The module writes a sequence of cut records, each `nnz, rhs, (col, val) * nnz` as
//!   `f64` (cut `sum val * x[col] <= rhs`, assumed globally valid), and returns the number of
//!   `f64` slots it wrote (0 = no cuts); -1 means "no opinion". A count above `out_cap`, a
//!   negative count other than -1, a truncated record or a non-integer `nnz`/`col` is invalid
//!   and drops every cut of that call; a well-formed cut that is not valid for the model
//!   (index out of range, non-finite, cuts off the debug solution) is rejected by the host gate
//!   like any user cut. Needs an exported `memory`.
//!
//! An exported function that is absent is simply not hooked (no stats counted). A function with
//! the wrong signature is rejected when the hooks are built.

use rusolve_ekk::mip::hooks::{
    BranchingInfo, CutScoreInfo, MipHooks, OpenNodeInfo, SeparationInfo, UserCut,
};
use std::sync::{Arc, Mutex};
use wasmi::{
    Config, Engine, Linker, Memory, Module, Store, StoreLimits, StoreLimitsBuilder, TrapCode,
    TypedFunc,
};

/// Budget and caps for one component.
#[derive(Clone, Copy, Debug)]
pub struct WasmLimits {
    /// Wasmi fuel granted before every hook call (deterministic: independent of the machine).
    pub fuel_per_call: u64,
    /// Linear-memory cap in bytes (also applies to the declared initial size).
    pub max_memory_bytes: usize,
}

impl Default for WasmLimits {
    fn default() -> Self {
        WasmLimits {
            fuel_per_call: 1_000_000,
            max_memory_bytes: 1 << 20,
        }
    }
}

/// What the component did. `calls` counts invocations; each one is exactly one of `ok`,
/// `declined` (valid "no opinion"), `traps`, `out_of_fuel` or `invalid` (wrong shape or NaN).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct WasmStats {
    pub calls: u64,
    pub ok: u64,
    pub declined: u64,
    pub traps: u64,
    pub out_of_fuel: u64,
    pub invalid: u64,
    /// Total fuel consumed over all calls.
    pub fuel_used: u64,
}

pub type WasmStatsHandle = Arc<Mutex<WasmStats>>;

type BranchFn = TypedFunc<(i32, i32, i32, i32), i32>;
type CutFn = TypedFunc<(f64, i32, i32, f64, f64), f64>;
type SelectFn = TypedFunc<(i32, i32), i32>;
type SepFn = TypedFunc<(i32, i32, i32, i32), i32>;

/// The sandboxed component (`MipHooks` implementation).
pub struct WasmHooks {
    store: Store<StoreLimits>,
    memory: Option<Memory>,
    branch: Option<BranchFn>,
    cut: Option<CutFn>,
    select: Option<SelectFn>,
    sep: Option<SepFn>,
    limits: WasmLimits,
    stats: WasmStatsHandle,
    name: String,
}

impl WasmHooks {
    /// Build from a binary module or WAT text (text is detected by a leading `(`).
    /// The module must have no imports.
    pub fn new(name: &str, source: &[u8], limits: WasmLimits) -> Result<Self, String> {
        let is_text = source.iter().find(|b| !b.is_ascii_whitespace()) == Some(&b'(');
        let bytes: Vec<u8> = if is_text {
            wat::parse_bytes(source)
                .map_err(|e| format!("wat: {e}"))?
                .into_owned()
        } else {
            source.to_vec()
        };
        let mut config = Config::default();
        config.consume_fuel(true);
        let engine = Engine::new(&config);
        let module = Module::new(&engine, &bytes[..]).map_err(|e| format!("module: {e}"))?;
        if let Some(imp) = module.imports().next() {
            return Err(format!(
                "module imports '{}.{}': components get no host access",
                imp.module(),
                imp.name()
            ));
        }
        let mut store = Store::new(
            &engine,
            StoreLimitsBuilder::new()
                .memory_size(limits.max_memory_bytes)
                .memories(1)
                .tables(1)
                .table_elements(1024)
                .instances(1)
                .trap_on_grow_failure(false)
                .build(),
        );
        store.limiter(|l| l);
        // the start function, if any, runs under the same budget
        store
            .set_fuel(limits.fuel_per_call)
            .map_err(|e| e.to_string())?;
        let instance = Linker::<StoreLimits>::new(&engine)
            .instantiate_and_start(&mut store, &module)
            .map_err(|e| format!("instantiate: {e}"))?;
        let memory = instance.get_memory(&store, "memory");
        let branch = if instance.get_export(&store, "branch_score").is_some() {
            Some(
                instance
                    .get_typed_func::<(i32, i32, i32, i32), i32>(&store, "branch_score")
                    .map_err(|e| format!("branch_score signature: {e}"))?,
            )
        } else {
            None
        };
        let cut = if instance.get_export(&store, "cut_score").is_some() {
            Some(
                instance
                    .get_typed_func::<(f64, i32, i32, f64, f64), f64>(&store, "cut_score")
                    .map_err(|e| format!("cut_score signature: {e}"))?,
            )
        } else {
            None
        };
        let select = if instance.get_export(&store, "select_node").is_some() {
            Some(
                instance
                    .get_typed_func::<(i32, i32), i32>(&store, "select_node")
                    .map_err(|e| format!("select_node signature: {e}"))?,
            )
        } else {
            None
        };
        let sep = if instance.get_export(&store, "separate").is_some() {
            Some(
                instance
                    .get_typed_func::<(i32, i32, i32, i32), i32>(&store, "separate")
                    .map_err(|e| format!("separate signature: {e}"))?,
            )
        } else {
            None
        };
        if sep.is_some() && memory.is_none() {
            return Err("separate needs an exported 'memory'".into());
        }
        if branch.is_some() && memory.is_none() {
            return Err("branch_score needs an exported 'memory'".into());
        }
        if select.is_some() && memory.is_none() {
            return Err("select_node needs an exported 'memory'".into());
        }
        Ok(WasmHooks {
            store,
            memory,
            branch,
            cut,
            select,
            sep,
            limits,
            stats: Arc::new(Mutex::new(WasmStats::default())),
            name: name.to_string(),
        })
    }

    /// Counters, shared: keep the handle before boxing the hooks into the solver.
    pub fn stats(&self) -> WasmStatsHandle {
        Arc::clone(&self.stats)
    }

    fn bump(&self, f: impl FnOnce(&mut WasmStats)) {
        f(&mut self.stats.lock().unwrap_or_else(|e| e.into_inner()));
    }

    /// Fuel accounting and failure classification for one finished call.
    fn finish<T>(&mut self, r: Result<T, wasmi::Error>) -> Option<T> {
        let left = self.store.get_fuel().unwrap_or(0);
        let used = self.limits.fuel_per_call.saturating_sub(left);
        let (v, trap, fuel) = match r {
            Ok(v) => (Some(v), false, false),
            Err(e) => {
                let fuel = e.as_trap_code() == Some(TrapCode::OutOfFuel);
                (None, !fuel, fuel)
            }
        };
        self.bump(|s| {
            s.calls += 1;
            s.fuel_used += used;
            s.traps += trap as u64;
            s.out_of_fuel += fuel as u64;
        });
        v
    }
}

impl MipHooks for WasmHooks {
    fn name(&self) -> &str {
        &self.name
    }

    fn branching_score(&mut self, info: &BranchingInfo<'_>) -> Option<Vec<f64>> {
        let f = self.branch?;
        let mem = self.memory?;
        let n = info.candidates.len();
        let out_ptr = 32 * n;
        if out_ptr + 8 * n > mem.data_size(&self.store) {
            self.bump(|s| {
                s.calls += 1;
                s.invalid += 1;
            });
            return None;
        }
        let mut buf = Vec::with_capacity(32 * n);
        for c in info.candidates {
            let j = c.col as usize;
            for v in [c.col as f64, c.value, info.col_lower[j], info.col_upper[j]] {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        if mem.write(&mut self.store, 0, &buf).is_err() {
            self.bump(|s| {
                s.calls += 1;
                s.invalid += 1;
            });
            return None;
        }
        let _ = self.store.set_fuel(self.limits.fuel_per_call);
        let r = f.call(&mut self.store, (info.depth, n as i32, 0, out_ptr as i32));
        let status = self.finish(r)?;
        if status != 0 {
            self.bump(|s| s.declined += 1);
            return None;
        }
        let mut out = vec![0u8; 8 * n];
        if mem.read(&self.store, out_ptr, &mut out).is_err() {
            self.bump(|s| s.invalid += 1);
            return None;
        }
        let scores: Vec<f64> = out
            .as_chunks::<8>()
            .0
            .iter()
            .map(|b| f64::from_le_bytes(*b))
            .collect();
        if scores.iter().any(|s| s.is_nan()) {
            self.bump(|s| s.invalid += 1);
            return None;
        }
        self.bump(|s| s.ok += 1);
        Some(scores)
    }

    fn select_node(&mut self, open: &[OpenNodeInfo]) -> Option<usize> {
        let f = self.select?;
        let mem = self.memory?;
        let n = open.len();
        let need = 32 * n;
        if need > mem.data_size(&self.store) {
            // the n records do not fit in the component's (capped) memory: no call is made
            self.bump(|s| {
                s.calls += 1;
                s.declined += 1;
            });
            return None;
        }
        let mut buf = Vec::with_capacity(need);
        for o in open {
            for v in [
                o.lower_bound,
                o.estimate,
                o.depth as f64,
                o.num_domain_changes as f64,
            ] {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        if mem.write(&mut self.store, 0, &buf).is_err() {
            self.bump(|s| {
                s.calls += 1;
                s.invalid += 1;
            });
            return None;
        }
        let _ = self.store.set_fuel(self.limits.fuel_per_call);
        let r = f.call(&mut self.store, (n as i32, 0));
        let pick = self.finish(r)?;
        if pick == -1 {
            self.bump(|s| s.declined += 1);
            return None;
        }
        if pick < 0 || pick as usize >= n {
            self.bump(|s| s.invalid += 1);
            return None;
        }
        self.bump(|s| s.ok += 1);
        Some(pick as usize)
    }

    fn separate(&mut self, info: &SeparationInfo<'_>) -> Vec<UserCut> {
        let (Some(f), Some(mem)) = (self.sep, self.memory) else {
            return Vec::new();
        };
        let n = info.solution.len();
        let size = mem.data_size(&self.store);
        let out_ptr = 8 * n;
        // the solution must fit with room for at least one record; else no call is made
        if out_ptr + 16 > size || out_ptr + 16 > i32::MAX as usize {
            self.bump(|s| {
                s.calls += 1;
                s.declined += 1;
            });
            return Vec::new();
        }
        let out_cap = (size - out_ptr) / 8;
        let mut buf = Vec::with_capacity(out_ptr);
        for v in info.solution {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        if mem.write(&mut self.store, 0, &buf).is_err() {
            self.bump(|s| {
                s.calls += 1;
                s.invalid += 1;
            });
            return Vec::new();
        }
        let _ = self.store.set_fuel(self.limits.fuel_per_call);
        let r = f.call(
            &mut self.store,
            (n as i32, 0, out_ptr as i32, out_cap as i32),
        );
        let Some(count) = self.finish(r) else {
            return Vec::new();
        };
        if count == -1 {
            self.bump(|s| s.declined += 1);
            return Vec::new();
        }
        if count < 0 || count as usize > out_cap {
            self.bump(|s| s.invalid += 1);
            return Vec::new();
        }
        let mut raw = vec![0u8; 8 * count as usize];
        if mem.read(&self.store, out_ptr, &mut raw).is_err() {
            self.bump(|s| s.invalid += 1);
            return Vec::new();
        }
        let vals: Vec<f64> = raw
            .as_chunks::<8>()
            .0
            .iter()
            .map(|b| f64::from_le_bytes(*b))
            .collect();
        match parse_cuts(&vals) {
            Some(cuts) => {
                self.bump(|s| s.ok += 1);
                cuts
            }
            None => {
                self.bump(|s| s.invalid += 1);
                Vec::new()
            }
        }
    }

    fn cut_score(&mut self, info: &CutScoreInfo) -> Option<f64> {
        let f = self.cut?;
        let _ = self.store.set_fuel(self.limits.fuel_per_call);
        let r = f.call(
            &mut self.store,
            (
                info.violation,
                info.num_nonzeros,
                info.num_active_nonzeros,
                info.active_norm,
                info.default_score,
            ),
        );
        let v = self.finish(r)?;
        if v.is_finite() {
            self.bump(|s| s.ok += 1);
            Some(v)
        } else {
            self.bump(|s| s.invalid += 1);
            None
        }
    }
}

/// Cut records `nnz, rhs, (col, val) * nnz` -> cuts, or `None` when a record is truncated or a
/// count or column is not a non-negative integer.
fn parse_cuts(vals: &[f64]) -> Option<Vec<UserCut>> {
    let int =
        |v: f64| (v.is_finite() && v >= 0.0 && v.fract() == 0.0 && v < 2e9).then_some(v as usize);
    let mut cuts = Vec::new();
    let mut at = 0;
    while at < vals.len() {
        let nnz = int(*vals.get(at)?)?;
        let rhs = *vals.get(at + 1)?;
        let end = at.checked_add(2)?.checked_add(nnz.checked_mul(2)?)?;
        if end > vals.len() {
            return None;
        }
        let mut indices = Vec::with_capacity(nnz);
        let mut values = Vec::with_capacity(nnz);
        for k in 0..nnz {
            indices.push(int(vals[at + 2 + 2 * k])? as i32);
            values.push(vals[at + 3 + 2 * k]);
        }
        cuts.push(UserCut {
            indices,
            values,
            rhs,
        });
        at = end;
    }
    Some(cuts)
}
