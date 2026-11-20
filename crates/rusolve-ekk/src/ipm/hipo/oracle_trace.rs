//! HiGHS source: none (oracle trace; mirrors `HighsOracleW4Dump.h`)
//!
//! Rust mirror of the oracle HiPO trace: `tools/highs-oracle/highs-oracle-w4-hooks.patch`
//! (`highs/ipm/hipo/HighsOracleW4Dump.h`). Not part of upstream HiGHS.
//!
//! The patched C++ appends one JSON line at each hook site:
//!
//! | record      | site |
//! |-------------|------|
//! | `hipo_ord`  | `FactorHiGHSSolver.cpp::chooseOrdering`, inside `run_ordering_and_analyse`, right after the `logger_.printInfo("Finished %s for %s\n", ...)` line (before the `failure[i]` test) |
//! | `hipo_fact` | `FactorHiGHSSolver.cpp::setup`, right after `setParallel()` |
//! | `hipo_iter` | `Solver.cpp::runIpm`, right after `makeStep()` |
//! | `hipo_end`  | `IpxWrapper.cpp::solveHipo`, right after `const hipo::Info hipo_info = hipo.getInfo()` |
//!
//! The Rust writes the identical records at the same points, so a fixture
//! comparison is a byte diff. Doubles are the 16 hex digits of their bit
//! pattern; double vectors are FNV-1a 64 hashes of the little-endian bit
//! patterns; integer vectors (`ptr`, `rows`, `perm`) are written in full, so
//! that `tests/t_w4_orderings.rs` can feed the pattern straight into the
//! ported METIS / AMD / RCM (`tools/highs-oracle/FORMAT.md`, HiPO records).
//!
//! Sinks: [`capture_begin`] / [`capture_end`] install an in-memory buffer on
//! the current thread (tests). Otherwise, when the environment variable
//! [`kHipoTraceEnv`] (`RUSOLVE_HIPO_TRACE`) names a file, it is truncated by
//! the first record of the thread and appended to afterwards (the C++
//! `HIGHS_ORACLE_HIPO_DUMP` file is only appended to, since the oracle writes
//! the `w4_header` line into it before the run). Without either nothing is
//! written and nothing is evaluated. Same design as `pdlp/oracle_trace.rs`.

use super::Int;
use super::ipm::info::Info;
use std::cell::RefCell;
use std::fmt::Write as _;
use std::io::Write as _;

/// Environment variable naming the trace file (C++ `HIGHS_ORACLE_HIPO_DUMP`).
#[allow(non_upper_case_globals)]
pub const kHipoTraceEnv: &str = "RUSOLVE_HIPO_TRACE";

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
pub fn enabled() -> bool {
    TRACE.with(|t| {
        let mut t = t.borrow_mut();
        if let Sink::Unopened = t.sink {
            t.sink = match std::env::var_os(kHipoTraceEnv) {
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
fn emit(head: &str, body: impl FnOnce(&mut String)) {
    if !enabled() {
        return;
    }
    TRACE.with(|t| {
        let mut t = t.borrow_mut();
        let call = t.call;
        t.call += 1;
        let mut s = String::new();
        let _ = write!(s, "{{\"t\":\"{head}\",\"call\":{call}");
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

fn put_s(s: &mut String, key: &str, v: &str) {
    let _ = write!(s, ",\"{key}\":\"{v}\"");
}

/// An integer vector in full: `,"key":[1,2,3]` (`[]` when empty).
fn put_iv(s: &mut String, key: &str, v: &[Int]) {
    let _ = write!(s, ",\"{key}\":[");
    for (k, x) in v.iter().enumerate() {
        if k > 0 {
            s.push(',');
        }
        let _ = write!(s, "{x}");
    }
    s.push(']');
}

/// FNV-1a 64 over the little-endian bit patterns of `v` (same hash as
/// `pdlp::oracle_trace::hash_f64` over the whole slice).
pub fn hash_f64(v: &[f64]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for x in v {
        let bits = x.to_bits();
        for b in 0..8 {
            h ^= (bits >> (8 * b)) & 0xff;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    h
}

fn put_h(s: &mut String, key: &str, v: &[f64]) {
    let _ = write!(s, ",\"{key}\":\"{:016x}\"", hash_f64(v));
}

/// One fill-reducing ordering (record `hipo_ord`). `nla` is the C++ `nla`
/// argument of `chooseOrdering` (`"AS"` / `"NE"`), `ordering` is `orderings_to_try[i]`, `ptr` /
/// `rows` are `full_ptr` / `full_rows` (the full symmetric pattern passed
/// to the ordering, 0-based CSC), `perm` is `permutations[i]` after the
/// call and `failure` is `failure[i]` after the ordering (before
/// `FH_.analyse`).
#[allow(clippy::too_many_arguments)]
pub fn hipo_ord(
    nla: &str,
    ordering: &str,
    n: Int,
    ptr: &[Int],
    rows: &[Int],
    perm: &[Int],
    failure: bool,
) {
    emit("hipo_ord", |s| {
        put_s(s, "nla", nla);
        put_s(s, "ordering", ordering);
        put_i(s, "n", n as i64);
        put_i(s, "failure", failure as i64);
        put_iv(s, "ptr", ptr);
        put_iv(s, "rows", rows);
        put_iv(s, "perm", perm);
    });
}

/// The chosen factorisation (record `hipo_fact`): `nla` is `"AS"` for
/// `kOptionNlaAugmented` and `"NE"` for `kOptionNlaNormEq` (`options_.nla`
/// after `setNla`), `ordering` is `ordering_AS_` / `ordering_NE_` for that
/// system, `flops` / `spops` / `nz` are `S_.flops()`, `S_.spops()`,
/// `(double)S_.nz()`.
pub fn hipo_fact(nla: &str, ordering: &str, flops: f64, spops: f64, nz: f64) {
    emit("hipo_fact", |s| {
        put_s(s, "nla", nla);
        put_s(s, "ordering", ordering);
        put_f(s, "flops", flops);
        put_f(s, "spops", spops);
        put_f(s, "nz", nz);
    });
}

/// The scalars of one IPM iteration (record `hipo_iter`), read after
/// `makeStep()`: `iter` = `iter_`, `alpha_p` / `alpha_d` / `sigma` =
/// `alpha_primal_` / `alpha_dual_` / `sigma_`, `correctors` =
/// `info_.correctors`, the rest from `it_` (`pobj`, `dobj`, `pinf`, `dinf`,
/// `pdgap`, `mu`, the values of the last `indicators()` / `computeMu()`).
#[derive(Clone, Copy, Debug, Default)]
pub struct IterRecord {
    pub iter: Int,
    pub pobj: f64,
    pub dobj: f64,
    pub pinf: f64,
    pub dinf: f64,
    pub pdgap: f64,
    pub mu: f64,
    pub alpha_p: f64,
    pub alpha_d: f64,
    pub sigma: f64,
    pub correctors: Int,
}

/// One IPM iteration (record `hipo_iter`); the vectors are `it_->x`, `xl`,
/// `xu`, `y`, `zl`, `zu` (full length).
#[allow(clippy::too_many_arguments)]
pub fn hipo_iter(
    r: &IterRecord,
    x: &[f64],
    xl: &[f64],
    xu: &[f64],
    y: &[f64],
    zl: &[f64],
    zu: &[f64],
) {
    emit("hipo_iter", |s| {
        put_i(s, "iter", r.iter as i64);
        put_f(s, "pobj", r.pobj);
        put_f(s, "dobj", r.dobj);
        put_f(s, "pinf", r.pinf);
        put_f(s, "dinf", r.dinf);
        put_f(s, "pdgap", r.pdgap);
        put_f(s, "mu", r.mu);
        put_f(s, "alpha_p", r.alpha_p);
        put_f(s, "alpha_d", r.alpha_d);
        put_f(s, "sigma", r.sigma);
        put_i(s, "correctors", r.correctors as i64);
        put_h(s, "h_x", x);
        put_h(s, "h_xl", xl);
        put_h(s, "h_xu", xu);
        put_h(s, "h_y", y);
        put_h(s, "h_zl", zl);
        put_h(s, "h_zu", zu);
    });
}

/// End of the HiPO solve (record `hipo_end`): the `hipo::Info` returned by
/// `getInfo()` in `solveHipo`.
pub fn hipo_end(info: &Info) {
    emit("hipo_end", |s| {
        put_i(s, "status", info.status as i64);
        put_i(s, "ipm_iter", info.ipm_iter as i64);
        put_i(s, "ipx_used", info.ipx_used as i64);
        put_i(s, "ipx_iter", info.ipx_info.iter as i64);
        put_i(s, "status_crossover", info.ipx_info.status_crossover as i64);
        put_i(
            s,
            "updates_crossover",
            info.ipx_info.updates_crossover as i64,
        );
        put_i(s, "m_solver", info.m_solver as i64);
        put_i(s, "n_solver", info.n_solver as i64);
        put_f(s, "p_obj", info.p_obj);
        put_f(s, "d_obj", info.d_obj);
        put_f(s, "pd_gap", info.pd_gap);
        put_f(s, "p_res_abs", info.p_res_abs);
        put_f(s, "p_res_rel", info.p_res_rel);
        put_f(s, "d_res_abs", info.d_res_abs);
        put_f(s, "d_res_rel", info.d_res_rel);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capture_round_trip() {
        capture_begin();
        hipo_ord("NE", "amd", 2, &[0, 1, 2], &[1, 0], &[1, 0], false);
        hipo_fact("NE", "amd", 1.0, 2.0, 3.0);
        hipo_end(&Info::default());
        let out = capture_end();
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(
            lines[0],
            "{\"t\":\"hipo_ord\",\"call\":0,\"nla\":\"NE\",\"ordering\":\"amd\",\"n\":2,\"failure\":0,\"ptr\":[0,1,2],\"rows\":[1,0],\"perm\":[1,0]}"
        );
        assert!(lines[1].starts_with("{\"t\":\"hipo_fact\",\"call\":1,\"nla\":\"NE\",\"ordering\":\"amd\",\"flops\":\"3ff0000000000000\""));
        // `Info::status` defaults to kStatusNotRun (1).
        assert!(
            lines[2].starts_with("{\"t\":\"hipo_end\",\"call\":2,\"status\":1,\"ipm_iter\":0,")
        );
        assert!(capture_end().is_empty());
    }

    #[test]
    fn fnv_of_empty_and_one() {
        assert_eq!(hash_f64(&[]), 0xcbf29ce484222325);
        assert_ne!(hash_f64(&[1.0]), hash_f64(&[2.0]));
    }
}
