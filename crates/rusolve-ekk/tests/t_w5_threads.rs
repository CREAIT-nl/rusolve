//! Golden test for parallelism at threads = 1, 2, 4, 8: the
//! `w5/simplex/*_s2_*` and `w5/hipo/` fixtures (`highs-oracle w5-simplex` / `w5-hipo`).
//! PAMI is `t_w5_pami.rs`, MIP is `rusolve-cli/tests/mip_trace_threads.rs`.
//!
//! Per fixture, the port runs with the fixture's `threads` and must reproduce the
//! oracle trace byte for byte, as the HiPO replay does (`tests/port/w4.rs`). LP (SIP,
//! `simplex_strategy` = 2) and HiPO are also checked *across* thread counts: the port's
//! trace at threads 2, 4, 8 equals its trace at threads 1.
//!
//! All these fixtures are target-only; with none generated a test prints a SKIP banner
//! (RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a failure).

mod common;
mod port;

use std::collections::BTreeMap;

use common::*;
use port::w5::*;
use port::{Failures, run_trace_files};
use rusolve_ekk::prelude::*;
use serde_json::Value;

/// The port's SIP records with the `options.threads` echo removed: that field is
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
#[ignore = "target-only parallel fixtures (highs-oracle w5-simplex, OR5): cargo test -p rusolve-ekk --release --test t_w5_threads -- --include-ignored"]
fn w5_threads_sip_matches_highs() {
    let files: Vec<String> = w5_sip_fixtures().into_iter().map(|fx| fx.file).collect();
    if files.is_empty() {
        return;
    }
    run_trace_files(files, "t_w5_threads (SIP, threads 1/2/4/8)", |fx| {
        sip_trace_gate(&fx.file)
    });
}

#[test]
#[ignore = "target-only parallel fixtures (highs-oracle w5-simplex, OR5): cargo test -p rusolve-ekk --release --test t_w5_threads -- --include-ignored"]
fn w5_threads_sip_identical_across_thread_counts() {
    let files: Vec<String> = w5_sip_fixtures().into_iter().map(|fx| fx.file).collect();
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
        let base = sans_threads(sip_records(&base_file));
        for file in counts.iter().skip(1) {
            let got = sans_threads(sip_records(file));
            let mut chk = Checker::new(format!(
                "{group}: threads={} vs threads={}",
                threads_of(file),
                threads_of(&base_file)
            ));
            chk.value("sip_records", &base, &got);
            fails.take(&chk);
        }
    }
    fails.finish("t_w5_threads (SIP across thread counts)");
}

#[test]
#[ignore = "target-only parallel fixtures (highs-oracle w5-hipo, OR5): cargo test -p rusolve-ekk --release --test t_w5_threads -- --include-ignored"]
fn w5_threads_hipo_matches_highs() {
    let fxs = w5_hipo_fixtures();
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    let mut traces: BTreeMap<String, Vec<(HighsInt, String)>> = BTreeMap::new();
    for fx in &fxs {
        let file = fx.file.clone();
        let (chk, got) = hipo_replay(&file);
        traces
            .entry(group_of(&file))
            .or_default()
            .push((threads_of(&file), got));
        fails.take(&chk);
    }
    for (group, mut runs) in traces {
        runs.sort_by_key(|r| r.0);
        let (base_threads, base) = runs[0].clone();
        for (threads, got) in runs.iter().skip(1) {
            let mut chk = Checker::new(format!(
                "{group}: threads={threads} vs threads={base_threads}"
            ));
            chk.str("port trace", &base, got);
            fails.take(&chk);
        }
    }
    fails.finish("t_w5_threads (HiPO, oracle and across thread counts)");
}

#[test]
#[ignore = "target-only parallel fixtures (highs-oracle w5-simplex / w5-hipo, OR5): cargo test -p rusolve-ekk --release --test t_w5_threads -- --include-ignored"]
fn w5_threads_run_twice_is_identical() {
    let hipo: Vec<String> = w5_hipo_fixtures()
        .into_iter()
        .map(|fx| fx.file)
        .filter(|f| threads_of(f) == 8)
        .collect();
    let sip: Vec<String> = w5_sip_fixtures()
        .into_iter()
        .map(|fx| fx.file)
        .filter(|f| threads_of(f) == 8)
        .collect();
    if hipo.is_empty() && sip.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for file in &hipo {
        let (first, second) = hipo_replay_twice(file);
        let mut chk = Checker::new(format!("{file}: run 2 vs run 1"));
        chk.str("port trace", &first, &second);
        fails.take(&chk);
    }
    for file in &sip {
        let (first, second) = sip_records_twice(file);
        let (first, second) = (Value::Array(first), Value::Array(second));
        let mut chk = Checker::new(format!("{file}: run 2 vs run 1"));
        chk.value("sip_records", &first, &second);
        fails.take(&chk);
    }
    fails.finish("t_w5_threads (run twice at threads=8)");
}
