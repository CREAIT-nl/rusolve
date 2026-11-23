//! Rust mirror of the oracle PDLP trace: `tools/highs-oracle/highs-oracle-w3-hooks.patch`
//! (`highs/pdlp/HighsOracleW3Dump.h`). Not part of upstream HiGHS.
//!
//! The patched C++ appends one JSON line at each hook site:
//!
//! | record           | site |
//! |------------------|------|
//! | `pdlp_check`     | `cupdlp_solver.c::PDHG_Solve`, right after `PDHG_Compute_Infeas_Residuals` in the `bool_checking` block |
//! | `pdlp_restart`   | `cupdlp_proj.c::PDHG_Restart_Iterate_GPU`, last statement (only when a restart happened) |
//! | `pdlp_end`       | `cupdlp_solver.c::PDHG_Solve`, right after the `for` loop (every exit path) |
//! | `hipdlp_check`   | `pdhg.cc::PDLPSolver::solve`, right after each `runConvergenceCheck` call (initial and in-loop) |
//! | `hipdlp_restart` | `pdhg.cc::PDLPSolver::solve`, end of the `if (do_restart)` block |
//! | `hipdlp_end`     | `pdhg.cc::PDLPSolver::solveReturn`, after `results_.term_code = term_code` |
//!
//! The Rust writes the identical records at the same points, so a fixture
//! comparison is a byte diff. Doubles are the 16 hex digits of their bit
//! pattern; vectors are FNV-1a 64 hashes of the little-endian bit patterns
//! of the dense array (`tools/highs-oracle/FORMAT.md`, "Hashes" and the PDLP records).
//!
//! Sinks: [`capture_begin`] / [`capture_end`] install an in-memory buffer on
//! the current thread (tests). Otherwise, when the environment variable
//! [`kPdlpTraceEnv`] (`RUSOLVE_PDLP_TRACE`) names a file, it is truncated by
//! the first record of the thread and appended to afterwards (the C++
//! `HIGHS_ORACLE_PDLP_DUMP` file is only appended to, since the oracle
//! writes the `w3_header` line into it before the run). Without either
//! nothing is written and nothing is evaluated. Same design as
//! `qpsolver/oracle_trace.rs`.

use super::cupdlp::{CUPDLPwork, PDHG_restart_choice};
use super::hipdlp::PDLPSolver;
use std::cell::RefCell;
use std::fmt::Write as _;
use std::io::Write as _;

/// Environment variable naming the trace file (C++ `HIGHS_ORACLE_PDLP_DUMP`).
#[allow(non_upper_case_globals)]
pub const kPdlpTraceEnv: &str = "RUSOLVE_PDLP_TRACE";

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
            t.sink = match std::env::var_os(kPdlpTraceEnv) {
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

/// FNV-1a 64 over the little-endian bit patterns of `v[..n]` (`n` clipped
/// to `v.len()`; the C++ hashes exactly `n` entries of the raw array).
pub fn hash_f64(v: &[f64], n: usize) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for x in &v[..n.min(v.len())] {
        let bits = x.to_bits();
        for b in 0..8 {
            h ^= (bits >> (8 * b)) & 0xff;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    h
}

fn put_h(s: &mut String, key: &str, v: &[f64], n: usize) {
    let _ = write!(s, ",\"{key}\":\"{:016x}\"", hash_f64(v, n));
}

/// cuPDLP check iteration (record `pdlp_check`, see `tools/highs-oracle/FORMAT.md`).
/// `x` / `y` are `iterates->x[nIter % 2]` / `y[nIter % 2]`, hashed over
/// `nCols` / `nRows` entries.
pub fn pdlp_check(w: &CUPDLPwork) {
    emit("pdlp_check", |s| {
        let (it, st, r) = (&w.iterates, &w.stepsize, &w.resobj);
        let k = (w.timers.nIter % 2) as usize;
        let (n, m) = (
            w.problem.nCols.max(0) as usize,
            w.problem.nRows.max(0) as usize,
        );
        put_i(s, "iter", w.timers.nIter as i64);
        put_f(s, "p_step", st.dPrimalStep);
        put_f(s, "d_step", st.dDualStep);
        put_f(s, "beta", st.dBeta);
        put_i(s, "ss_iter", st.nStepSizeIter as i64);
        put_i(s, "last_restart", it.iLastRestartIter as i64);
        put_f(s, "pobj", r.dPrimalObj);
        put_f(s, "dobj", r.dDualObj);
        put_f(s, "gap", r.dDualityGap);
        put_f(s, "compl", r.dComplementarity);
        put_f(s, "pfeas", r.dPrimalFeas);
        put_f(s, "dfeas", r.dDualFeas);
        put_f(s, "relgap", r.dRelObjGap);
        put_f(s, "pobj_a", r.dPrimalObjAverage);
        put_f(s, "dobj_a", r.dDualObjAverage);
        put_f(s, "gap_a", r.dDualityGapAverage);
        put_f(s, "compl_a", r.dComplementarityAverage);
        put_f(s, "pfeas_a", r.dPrimalFeasAverage);
        put_f(s, "dfeas_a", r.dDualFeasAverage);
        put_f(s, "relgap_a", r.dRelObjGapAverage);
        put_f(s, "pinf_obj", r.dPrimalInfeasObj);
        put_f(s, "dinf_obj", r.dDualInfeasObj);
        put_f(s, "pinf_res", r.dPrimalInfeasRes);
        put_f(s, "dinf_res", r.dDualInfeasRes);
        put_h(s, "h_x", &it.x[k].data, n);
        put_h(s, "h_y", &it.y[k].data, m);
        put_h(s, "h_xa", &it.xAverage.data, n);
        put_h(s, "h_ya", &it.yAverage.data, m);
    });
}

/// cuPDLP restart (record `pdlp_restart`): `restart_choice` is
/// `PDHG_RESTART_TO_CURRENT` (1) or `PDHG_RESTART_TO_AVERAGE` (2).
pub fn pdlp_restart(w: &CUPDLPwork, restart_choice: PDHG_restart_choice) {
    emit("pdlp_restart", |s| {
        let (it, st, r) = (&w.iterates, &w.stepsize, &w.resobj);
        let k = (w.timers.nIter % 2) as usize;
        let (n, m) = (
            w.problem.nCols.max(0) as usize,
            w.problem.nRows.max(0) as usize,
        );
        put_i(s, "iter", w.timers.nIter as i64);
        put_i(s, "choice", restart_choice as i64);
        put_f(s, "p_step", st.dPrimalStep);
        put_f(s, "d_step", st.dDualStep);
        put_f(s, "beta", st.dBeta);
        put_f(s, "pfeas_lr", r.dPrimalFeasLastRestart);
        put_f(s, "dfeas_lr", r.dDualFeasLastRestart);
        put_f(s, "gap_lr", r.dDualityGapLastRestart);
        put_h(s, "h_x", &it.x[k].data, n);
        put_h(s, "h_y", &it.y[k].data, m);
    });
}

/// End of the cuPDLP main loop (record `pdlp_end`).
pub fn pdlp_end(w: &CUPDLPwork) {
    emit("pdlp_end", |s| {
        let (it, r) = (&w.iterates, &w.resobj);
        let k = (w.timers.nIter % 2) as usize;
        let (n, m) = (
            w.problem.nCols.max(0) as usize,
            w.problem.nRows.max(0) as usize,
        );
        put_i(s, "iter", w.timers.nIter as i64);
        put_i(s, "term_code", r.termCode as i64);
        put_i(s, "term_iterate", r.termIterate as i64);
        put_h(s, "h_x", &it.x[k].data, n);
        put_h(s, "h_y", &it.y[k].data, m);
        put_h(s, "h_xa", &it.xAverage.data, n);
        put_h(s, "h_ya", &it.yAverage.data, m);
    });
}

fn put_hipdlp_iterates(s: &mut String, p: &PDLPSolver) {
    put_h(s, "h_x", &p.x_current_, p.x_current_.len());
    put_h(s, "h_y", &p.y_current_, p.y_current_.len());
    put_h(s, "h_xn", &p.x_next_, p.x_next_.len());
    put_h(s, "h_yn", &p.y_next_, p.y_next_.len());
}

/// HiPDLP convergence check (record `hipdlp_check`); `stop` is the return
/// value of `runConvergenceCheck`.
pub fn hipdlp_check(p: &PDLPSolver, stop: bool) {
    emit("hipdlp_check", |s| {
        let r = &p.results_;
        put_i(s, "iter", p.final_iter_count_ as i64);
        put_i(s, "halpern", p.halpern_iteration_ as i64);
        put_i(s, "stop", stop as i64);
        put_f(s, "p_step", p.stepsize_.primal_step);
        put_f(s, "d_step", p.stepsize_.dual_step);
        put_f(s, "beta", p.stepsize_.beta);
        put_f(s, "omega", p.params_.omega);
        put_f(s, "eta", p.params_.eta);
        put_f(s, "pw", p.primal_weight_);
        put_f(s, "fpe", p.fpe_);
        put_f(s, "ifpe", p.initial_fpe_);
        put_f(s, "pobj", r.primal_obj);
        put_f(s, "dobj", r.dual_obj);
        put_f(s, "gap", r.duality_gap);
        put_f(s, "pfeas", r.primal_feasibility);
        put_f(s, "dfeas", r.dual_feasibility);
        put_f(s, "relgap", r.relative_obj_gap);
        put_hipdlp_iterates(s, p);
    });
}

/// HiPDLP restart (record `hipdlp_restart`).
pub fn hipdlp_restart(p: &PDLPSolver) {
    emit("hipdlp_restart", |s| {
        put_i(s, "iter", p.final_iter_count_ as i64);
        put_f(s, "p_step", p.stepsize_.primal_step);
        put_f(s, "d_step", p.stepsize_.dual_step);
        put_f(s, "beta", p.stepsize_.beta);
        put_f(s, "omega", p.params_.omega);
        put_f(s, "pw", p.primal_weight_);
        put_f(s, "best_pw", p.best_primal_weight_);
        put_f(s, "pw_err_sum", p.primal_weight_error_sum_);
        put_f(s, "pw_last_err", p.primal_weight_last_error_);
        put_hipdlp_iterates(s, p);
    });
}

/// HiPDLP end (record `hipdlp_end`): `term_code` is the `TerminationStatus`
/// integer value (`NOTSET` 0 ... `ERROR` 6).
pub fn hipdlp_end(p: &PDLPSolver) {
    emit("hipdlp_end", |s| {
        put_i(s, "iter", p.final_iter_count_ as i64);
        put_i(s, "term_code", p.results_.term_code as i64);
        put_hipdlp_iterates(s, p);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capture_round_trip() {
        capture_begin();
        // `SolverResults::term_code` defaults to TIMEOUT (5) as in the header.
        let p = PDLPSolver::default();
        hipdlp_end(&p);
        let w = CUPDLPwork::default();
        pdlp_end(&w);
        let out = capture_end();
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(
            lines[0].starts_with("{\"t\":\"hipdlp_end\",\"call\":0,\"iter\":0,\"term_code\":5,")
        );
        assert!(lines[1].starts_with("{\"t\":\"pdlp_end\",\"call\":1,\"iter\":0,"));
        assert!(out.ends_with("}\n"));
        assert!(capture_end().is_empty());
    }

    #[test]
    fn fnv_of_empty_and_one() {
        assert_eq!(hash_f64(&[], 0), 0xcbf29ce484222325);
        assert_eq!(hash_f64(&[1.0, 2.0], 0), 0xcbf29ce484222325);
        assert_ne!(hash_f64(&[1.0], 1), hash_f64(&[2.0], 1));
    }
}
