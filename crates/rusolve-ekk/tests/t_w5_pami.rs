//! Golden test for PAMI (`simplex_strategy` = 3, `HEkkDualMulti.cpp`) at threads = 1, 2,
//! 4, 8: the `w5/simplex/*_s3_t<N>.jsonl` fixtures (`highs-oracle
//! w5-simplex`).
//!
//! PAMI depends on `num_concurrency` (HEkk.cpp 1720-1785), not on the thread count, so
//! the port's trace must equal the oracle's at every thread count and the port's own
//! traces must be identical across thread counts.
//!
//! All these fixtures are target-only; with none generated a test prints a SKIP banner
//! (`wave_fixtures`; RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a failure).

mod common;
mod port;

use std::collections::BTreeMap;

use common::*;
use port::w1::wave_fixtures;
use port::w5::{fresh, group_of, threads_of};
use port::{Failures, enters_of, full_run_limit, iters_of, run_happ, run_trace_files, trace_gate};
use rusolve_ekk::parallel::highs_parallel::initialize_scheduler;
use rusolve_ekk::prelude::*;
use serde_json::Value;

/// The PAMI fixture files `w5/simplex/*_s3_t<N>.jsonl` (`simplex_strategy` 3).
fn pami_files() -> Vec<String> {
    wave_fixtures("w5/simplex/", "parallel oracle commands")
        .into_iter()
        .filter(|fx| fx.file.contains("_s3_t"))
        .map(|fx| fx.file)
        .collect()
}

/// The thread count of a PAMI fixture: its `trace_header`, which must agree with the
/// file name and record `simplex_strategy` 3 and `parallel` "on".
fn pami_threads(fx: &Fixture) -> HighsInt {
    let h = fx.rec(0);
    assert_eq!(h.t(), "trace_header", "{}: first record", fx.file);
    let threads = h.int("threads") as HighsInt;
    assert_eq!(threads, threads_of(&fx.file), "{}: header threads", fx.file);
    assert_eq!(h.int("simplex_strategy"), 3, "{}: PAMI fixture", fx.file);
    assert_eq!(h.str("parallel"), "on", "{}: header parallel", fx.file);
    threads
}

/// The port's PAMI trace records (`run_happ` with emit, as the oracle emits them), at
/// the fixture's own thread count. `run_happ` starts below `Highs::run`, so the
/// scheduler is initialised here, as `initializeMultiThreading` would do.
fn pami_run(fx: &Fixture) -> Vec<Value> {
    initialize_scheduler(pami_threads(fx));
    let enters = enters_of(fx);
    let limit = full_run_limit(&iters_of(fx));
    run_happ(&enters, &fx.file, limit, None, true).records
}

/// The PAMI gate of `file` (`trace_gate`, linear scan) at its own thread count.
fn pami_trace_gate(file: &str) -> Option<String> {
    let file = file.to_string();
    fresh(move || {
        let fx = load_required(&file);
        initialize_scheduler(pami_threads(&fx));
        trace_gate(&fx, true)
    })
}

/// The port's PAMI trace records for `file` at its own thread count.
fn pami_records(file: &str) -> Vec<Value> {
    let file = file.to_string();
    fresh(move || {
        let fx = load_required(&file);
        pami_run(&fx)
    })
}

/// The port's PAMI records with the `options.threads` echo removed: that field is
/// the thread count itself, the only field that differs across counts by design.
fn sans_threads(mut recs: Vec<Value>) -> Value {
    for rec in recs.iter_mut() {
        if let Some(opts) = rec.get_mut("options").and_then(Value::as_object_mut) {
            opts.remove("threads");
        }
    }
    Value::Array(recs)
}

#[test]
#[ignore = "target-only parallel fixtures (highs-oracle w5-simplex, OR5): cargo test -p rusolve-ekk --release --test t_w5_pami -- --include-ignored"]
fn w5_pami_matches_highs() {
    let files = pami_files();
    if files.is_empty() {
        return;
    }
    run_trace_files(files, "t_w5_pami (PAMI, threads 1/2/4/8)", |fx| {
        pami_trace_gate(&fx.file)
    });
}

#[test]
#[ignore = "target-only parallel fixtures (highs-oracle w5-simplex, OR5): cargo test -p rusolve-ekk --release --test t_w5_pami -- --include-ignored"]
fn w5_pami_identical_across_thread_counts() {
    let files = pami_files();
    if files.is_empty() {
        return;
    }
    let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for file in files {
        groups.entry(group_of(&file)).or_default().push(file);
    }
    let mut fails = Failures::default();
    for (group, mut counts) in groups {
        counts.sort_by_key(|f| threads_of(f));
        let base_file = counts[0].clone();
        let base = sans_threads(pami_records(&base_file));
        for file in counts.iter().skip(1) {
            let got = sans_threads(pami_records(file));
            let mut chk = Checker::new(format!(
                "{group}: threads={} vs threads={}",
                threads_of(file),
                threads_of(&base_file)
            ));
            chk.value("pami_records", &base, &got);
            fails.take(&chk);
        }
    }
    fails.finish("t_w5_pami (PAMI across thread counts)");
}
