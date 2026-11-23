//! Parallelism replay helpers for `tests/t_w5_threads.rs`.
//!
//! Each thread count runs on a fresh `std::thread` (own executor),
//! replays the `w5/` oracle fixtures, and exposes the port's captured trace so the
//! cross-thread-count identity checks compare port to port. Record layout:
//! tools/highs-oracle/FORMAT.md.

use crate::common::*;
use crate::port::w1::{diff_trace, hook_lines, wave_fixtures};
use crate::port::w4::{check_w4_result, w4_header, w4_lp, w4_options};
use crate::port::{enters_of, full_run_limit, iters_of, run_happ, trace_gate};
use rusolve_ekk::ipm::hipo::oracle_trace;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::parallel::highs_parallel::initialize_scheduler;
use rusolve_ekk::parallel::highs_task_executor::HighsTaskExecutor;
use rusolve_ekk::prelude::*;
use serde_json::Value;

/// The HiPO oracle-trace record types a parallel HiPO fixture carries.
pub const HIPO_HOOKS: &[&str] = &["hipo_ord", "hipo_fact", "hipo_iter", "hipo_end"];

/// Parallel SIP fixtures: `w5/simplex/*_s2_t{1,2,4,8}.jsonl`.
pub fn w5_sip_fixtures() -> Vec<Fixture> {
    wave_fixtures("w5/simplex/", "parallel oracle commands")
        .into_iter()
        .filter(|fx| fx.file.contains("_s2_t"))
        .collect()
}

/// Parallel HiPO fixtures: `w5/hipo/<case>_t{1,2,4,8}.jsonl`.
pub fn w5_hipo_fixtures() -> Vec<Fixture> {
    wave_fixtures("w5/hipo/", "parallel oracle commands")
}

/// The thread count encoded in a parallel fixture file name: the `N` in `..._t<N>.jsonl`.
pub fn threads_of(file: &str) -> HighsInt {
    let (_, tail) = file
        .rsplit_once("_t")
        .expect("parallel fixture name has _t<N>");
    tail.trim_end_matches(".jsonl")
        .parse::<HighsInt>()
        .expect("parallel fixture _t<N> is a number")
}

/// The fixture name without its `_t<N>.jsonl` suffix: the group of counts.
pub fn group_of(file: &str) -> String {
    let (base, _) = file
        .rsplit_once("_t")
        .expect("parallel fixture name has _t<N>");
    base.to_string()
}

/// Runs `case` on a fresh thread (its own `thread_local` executor); the
/// executor is shut down afterwards.
pub fn fresh<T: Send + 'static>(case: impl FnOnce() -> T + Send + 'static) -> T {
    let handle = std::thread::spawn(move || {
        let out = case();
        HighsTaskExecutor::shutdown(true);
        out
    });
    match handle.join() {
        Ok(out) => out,
        Err(e) => std::panic::resume_unwind(e),
    }
}

/// The oracle's options for a parallel HiPO fixture: the `w4_header` overrides,
/// with the `threads` and `parallel` overrides (which `w4_options` does not
/// know) taken from the same string, as the oracle applies them.
pub fn w5_options(overrides: &str) -> HighsOptions {
    let mut threads: Option<HighsInt> = None;
    let mut parallel: Option<String> = None;
    let mut rest: Vec<&str> = Vec::new();
    for kv in overrides.split(',') {
        if let Some(v) = kv.strip_prefix("threads=") {
            threads = Some(v.parse::<HighsInt>().expect("threads override is a number"));
        } else if let Some(v) = kv.strip_prefix("parallel=") {
            parallel = Some(v.to_string());
        } else {
            rest.push(kv);
        }
    }
    let mut o = w4_options(&rest.join(","));
    o.threads = threads.expect("parallel HiPO overrides set threads");
    o.parallel = parallel.expect("parallel HiPO overrides set parallel");
    o.sync_log_options();
    o
}

/// Replays one parallel HiPO fixture at `threads` (as `replay_w4` does for HiPO): the
/// hook lines must equal the fixture's byte for byte, then `w4_result` bit for
/// bit. Also returns the port's captured hook trace for the cross-count checks.
pub fn replay_w5_hipo(fx: &Fixture, threads: HighsInt, hooks: &[&str]) -> (Checker, String) {
    let text = std::fs::read_to_string(&fx.path).unwrap();
    let header = w4_header(fx);
    let options = w5_options(header.str("overrides"));
    assert_eq!(options.threads, threads, "{}: overrides threads", fx.file);
    assert_eq!(
        header.int("threads") as HighsInt,
        threads,
        "{}: header threads",
        fx.file
    );
    assert_eq!(
        header.str("parallel"),
        options.parallel,
        "{}: header parallel",
        fx.file
    );
    let want = hook_lines(fx, &text, hooks);
    let result = fx.first("w4_result");
    let mut chk = Checker::new(fx.file.clone());
    let mut highs = Highs::with_options(options);
    chk.int(
        "pass_model",
        HighsStatus::kOk as i32,
        highs.pass_model(w4_lp(fx)) as i32,
    );
    oracle_trace::capture_begin();
    let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| highs.run()));
    let got = oracle_trace::capture_end();
    let run_status = match run {
        Ok(s) => s,
        Err(e) => std::panic::resume_unwind(e),
    };
    diff_trace(&mut chk, &want, &got);
    check_w4_result(&mut chk, &result, run_status, &highs);
    (chk, got)
}

/// HiPO replay of `file` at its own thread count, on a fresh thread with no
/// preset scheduler (`Highs::run` initialises it, `initializeMultiThreading`).
pub fn hipo_replay(file: &str) -> (Checker, String) {
    let file = file.to_string();
    fresh(move || {
        let fx = load_required(&file);
        replay_w5_hipo(&fx, threads_of(&file), HIPO_HOOKS)
    })
}

/// Two HiPO replays of `file` on one thread, `Highs::resetGlobalScheduler(true)`
/// between them (the second `Highs::run` initialises a new scheduler): the
/// port's two hook traces.
pub fn hipo_replay_twice(file: &str) -> (String, String) {
    let file = file.to_string();
    fresh(move || {
        let fx = load_required(&file);
        let threads = threads_of(&file);
        let first = replay_w5_hipo(&fx, threads, HIPO_HOOKS).1;
        Highs::reset_global_scheduler(true);
        let second = replay_w5_hipo(&fx, threads, HIPO_HOOKS).1;
        (first, second)
    })
}

/// The thread count of a parallel SIP fixture: its `trace_header`, which must agree
/// with the file name and record `simplex_strategy` 2 and `parallel` "on".
fn sip_threads(fx: &Fixture) -> HighsInt {
    let h = fx.rec(0);
    assert_eq!(h.t(), "trace_header", "{}: first record", fx.file);
    let threads = h.int("threads") as HighsInt;
    assert_eq!(threads, threads_of(&fx.file), "{}: header threads", fx.file);
    assert_eq!(h.int("simplex_strategy"), 2, "{}: SIP fixture", fx.file);
    assert_eq!(h.str("parallel"), "on", "{}: header parallel", fx.file);
    threads
}

/// The port's SIP trace records (`run_happ` with emit, as the oracle emits them).
/// `run_happ` starts below `Highs::run`, so the scheduler `initializeMultiThreading`
/// would set up there is initialised here, to the fixture's `threads`.
fn sip_run(fx: &Fixture) -> Vec<Value> {
    initialize_scheduler(sip_threads(fx));
    let enters = enters_of(fx);
    let limit = full_run_limit(&iters_of(fx));
    run_happ(&enters, &fx.file, limit, None, true).records
}

/// The SIP gate of `file` (`trace_gate`, linear scan) at its own thread count.
pub fn sip_trace_gate(file: &str) -> Option<String> {
    let file = file.to_string();
    fresh(move || {
        let fx = load_required(&file);
        initialize_scheduler(sip_threads(&fx));
        trace_gate(&fx, true)
    })
}

/// The port's SIP trace records for `file` at its own thread count.
pub fn sip_records(file: &str) -> Vec<Value> {
    let file = file.to_string();
    fresh(move || {
        let fx = load_required(&file);
        sip_run(&fx)
    })
}

/// Two SIP runs of `file` on one thread, `Highs::resetGlobalScheduler(true)`
/// between them: the port's two record lists.
pub fn sip_records_twice(file: &str) -> (Vec<Value>, Vec<Value>) {
    let file = file.to_string();
    fresh(move || {
        let fx = load_required(&file);
        let first = sip_run(&fx);
        Highs::reset_global_scheduler(true);
        let second = sip_run(&fx);
        (first, second)
    })
}
