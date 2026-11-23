//! Rust mirror of the oracle IPX trace: `tools/highs-oracle/highs-oracle-ipx-hooks.patch`
//! (`highs/ipm/ipx/HighsOracleIpxDump.h`). Not part of upstream HiGHS.
//!
//! The patched C++ appends one JSON line after every `IPM::PrintOutput`
//! (`start`, `iter`, `centring`), one at the end of every `IPM::Driver`
//! (`driver_end`) and one at the end of `LpSolver::Solve` (`solve_end`).
//! The Rust writes the identical records at the same points, so a fixture
//! comparison is a byte diff. Doubles are the 16 hex digits of their bit
//! pattern; vectors are FNV-1a 64 hashes of the little-endian bit patterns.
//!
//! Sinks: [`capture_begin`] / [`capture_end`] install an in-memory buffer on
//! the current thread (tests). Otherwise, when the environment variable
//! [`kIpxTraceEnv`] (`RUSOLVE_IPX_TRACE`) names a file, it is truncated by the
//! first record of the thread and appended to afterwards, like the C++
//! `HIGHS_ORACLE_IPX_DUMP`. Without either nothing is written and nothing is
//! evaluated.

use super::iterate::Iterate;
use super::kkt_solver::KKTSolver;
use super::{Info, Vector};
use std::cell::RefCell;
use std::fmt::Write as _;
use std::io::Write as _;

/// Environment variable naming the trace file (C++ `HIGHS_ORACLE_IPX_DUMP`).
#[allow(non_upper_case_globals)]
pub const kIpxTraceEnv: &str = "RUSOLVE_IPX_TRACE";

enum Sink {
    Unopened,
    Off,
    File(std::fs::File),
    Buffer(String),
}

struct Trace {
    sink: Sink,
    call: i64,
}

thread_local! {
    static TRACE: RefCell<Trace> = const { RefCell::new(Trace { sink: Sink::Unopened, call: 0 }) };
}

/// Starts capturing the trace of this thread in memory (call counter reset).
pub fn capture_begin() {
    TRACE.with(|t| {
        *t.borrow_mut() = Trace {
            sink: Sink::Buffer(String::new()),
            call: 0,
        }
    });
}

/// Ends an in-memory capture and returns the records written since
/// [`capture_begin`]; the thread falls back to the environment variable.
pub fn capture_end() -> String {
    TRACE.with(|t| {
        let old = std::mem::replace(
            &mut *t.borrow_mut(),
            Trace {
                sink: Sink::Unopened,
                call: 0,
            },
        );
        match old.sink {
            Sink::Buffer(s) => s,
            _ => String::new(),
        }
    })
}

/// True when a record would be written (opens the env file on first use).
fn enabled() -> bool {
    TRACE.with(|t| {
        let mut t = t.borrow_mut();
        if let Sink::Unopened = t.sink {
            t.sink = match std::env::var_os(kIpxTraceEnv) {
                Some(p) if !p.is_empty() => match std::fs::File::create(&p) {
                    Ok(f) => Sink::File(f),
                    Err(_) => Sink::Off,
                },
                _ => Sink::Off,
            };
        }
        !matches!(t.sink, Sink::Off)
    })
}

/// Writes one record; `body` gets the call number.
fn emit(head: &str, kind: &str, body: impl FnOnce(&mut String)) {
    if !enabled() {
        return;
    }
    TRACE.with(|t| {
        let mut t = t.borrow_mut();
        let call = t.call;
        t.call += 1;
        let mut s = String::new();
        let _ = write!(s, "{{\"t\":\"{head}\",\"call\":{call},\"kind\":\"{kind}\"");
        body(&mut s);
        s.push_str("}\n");
        match &mut t.sink {
            Sink::Buffer(b) => b.push_str(&s),
            Sink::File(f) => {
                let _ = f.write_all(s.as_bytes());
                let _ = f.flush();
            }
            _ => {}
        }
    });
}

fn put_f(s: &mut String, key: &str, v: f64) {
    let _ = write!(s, ",\"{key}\":\"{:016x}\"", v.to_bits());
}

fn put_i(s: &mut String, key: &str, v: i64) {
    let _ = write!(s, ",\"{key}\":{v}");
}

/// FNV-1a 64 over the little-endian bit patterns (FORMAT.md "Hashes").
fn put_h(s: &mut String, key: &str, v: &Vector) {
    let mut h: u64 = 0xcbf29ce484222325;
    for x in v.iter() {
        let bits = x.to_bits();
        for b in 0..8 {
            h ^= (bits >> (8 * b)) & 0xff;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    let _ = write!(s, ",\"{key}\":\"{h:016x}\"");
}

/// `highsOracleIpxIter(kind, *iterate_, *info_, *kkt_, step_primal_, step_dual_)`.
pub fn ipx_iter(
    kind: &str,
    it: &Iterate,
    info: &Info,
    kkt: &dyn KKTSolver,
    step_primal: f64,
    step_dual: f64,
) {
    emit("ipm", kind, |s| {
        put_i(s, "iter", info.iter as i64);
        put_f(s, "mu", it.mu());
        put_f(s, "step_primal", step_primal);
        put_f(s, "step_dual", step_dual);
        put_f(s, "presidual", it.presidual());
        put_f(s, "dresidual", it.dresidual());
        put_f(s, "pobjective", it.pobjective());
        put_f(s, "dobjective", it.dobjective());
        put_f(
            s,
            "pobjective_after_postproc",
            it.pobjective_after_postproc(),
        );
        put_f(
            s,
            "dobjective_after_postproc",
            it.dobjective_after_postproc(),
        );
        put_f(s, "complementarity", it.complementarity());
        put_f(s, "mu_min", it.mu_min());
        put_f(s, "mu_max", it.mu_max());
        put_i(s, "kkt_iter", kkt.iter() as i64);
        put_i(s, "basis_changes", kkt.basis_changes() as i64);
        put_i(s, "kktiter1", info.kktiter1 as i64);
        put_i(s, "kktiter2", info.kktiter2 as i64);
        put_i(s, "updates_ipm", info.updates_ipm as i64);
        put_i(s, "primal_dropped", info.primal_dropped as i64);
        put_i(s, "dual_dropped", info.dual_dropped as i64);
        put_i(s, "errflag", info.errflag as i64);
        put_h(s, "h_x", it.x());
        put_h(s, "h_xl", it.xl());
        put_h(s, "h_xu", it.xu());
        put_h(s, "h_y", it.y());
        put_h(s, "h_zl", it.zl());
        put_h(s, "h_zu", it.zu());
    });
}

/// `highsOracleIpxInfo(kind, info)`: `ipx::Info` without the `time_*` fields.
pub fn ipx_info(kind: &str, info: &Info) {
    emit("info", kind, |s| {
        macro_rules! i {
            ($($f:ident)*) => { $( put_i(s, stringify!($f), info.$f as i64); )* };
        }
        macro_rules! d {
            ($($f:ident)*) => { $( put_f(s, stringify!($f), info.$f); )* };
        }
        i!(status status_ipm status_crossover errflag num_var num_constr num_entries
           num_rows_solver num_cols_solver num_entries_solver dualized dense_cols
           centring_tried centring_success dependent_rows dependent_cols
           rows_inconsistent cols_inconsistent primal_dropped dual_dropped);
        d!(abs_presidual abs_dresidual rel_presidual rel_dresidual pobjval dobjval
           rel_objgap complementarity normx normy normz objval primal_infeas dual_infeas);
        i!(iter kktiter1 kktiter2 basis_repairs updates_start updates_ipm updates_crossover);
        d!(ftran_sparse btran_sparse mean_fill max_fill);
        i!(maxvol_updates maxvol_skipped maxvol_passes tbl_nnz);
        d!(tbl_max frobnorm_squared lambdamax volume_increase);
    });
}
