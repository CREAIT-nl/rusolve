//! Throughput counters for the local-search heuristics (FJ, Local-MIP, FPR).
//!
//! Not part of the mip-heuristics port: rusolve instrumentation to test
//! whether local search pays only when move evaluation is fast (Gardi et al.,
//! *Mathematical Programming Solver Based on Local Search*, Wiley 2014). Counters are process-wide and
//! include every sub-MIP call. They never feed back into the search: the
//! heuristics read none of them.
//!
//! `RUSOLVE_LS_STATS=1` prints one cumulative `ls-stats` line to stderr after
//! each call; without it the counters are still kept (a few relaxed atomic
//! adds per call) but nothing is printed.

use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// Cumulative counters of one heuristic.
pub(crate) struct HeurCounters {
    name: &'static str,
    calls: AtomicU64,
    nanos: AtomicU64,
    effort: AtomicU64,
    evals: AtomicU64,
    commits: AtomicU64,
}

impl HeurCounters {
    const fn new(name: &'static str) -> Self {
        HeurCounters {
            name,
            calls: AtomicU64::new(0),
            nanos: AtomicU64::new(0),
            effort: AtomicU64::new(0),
            evals: AtomicU64::new(0),
            commits: AtomicU64::new(0),
        }
    }

    /// Moves scored (`evals`) and moves applied (`commits`) during a call.
    pub(crate) fn add_moves(&self, evals: u64, commits: u64) {
        self.evals.fetch_add(evals, Ordering::Relaxed);
        self.commits.fetch_add(commits, Ordering::Relaxed);
    }

    /// One finished call that started at `t0` and charged `effort` units.
    pub(crate) fn finish_call(&self, t0: Instant, effort: u64) {
        let ns = t0.elapsed().as_nanos() as u64;
        self.calls.fetch_add(1, Ordering::Relaxed);
        self.nanos.fetch_add(ns, Ordering::Relaxed);
        self.effort.fetch_add(effort, Ordering::Relaxed);
        if enabled() {
            eprintln!("{}", self.line());
        }
    }

    fn line(&self) -> String {
        let calls = self.calls.load(Ordering::Relaxed);
        let ns = self.nanos.load(Ordering::Relaxed);
        let effort = self.effort.load(Ordering::Relaxed);
        let evals = self.evals.load(Ordering::Relaxed);
        let commits = self.commits.load(Ordering::Relaxed);
        let s = (ns as f64) * 1e-9;
        let per_s = |x: u64| if s > 0.0 { (x as f64) / s } else { 0.0 };
        format!(
            "ls-stats heur={} calls={} ms={:.3} ms_per_call={:.4} effort={} evals={} commits={} \
             evals_per_s={:.0} commits_per_s={:.0} effort_per_s={:.0}",
            self.name,
            calls,
            s * 1e3,
            if calls > 0 {
                s * 1e3 / (calls as f64)
            } else {
                0.0
            },
            effort,
            evals,
            commits,
            per_s(evals),
            per_s(commits),
            per_s(effort)
        )
    }
}

pub(crate) static FJ: HeurCounters = HeurCounters::new("fj");
pub(crate) static LOCAL_MIP: HeurCounters = HeurCounters::new("local_mip");
pub(crate) static FPR: HeurCounters = HeurCounters::new("fpr");

/// Whether `RUSOLVE_LS_STATS=1` is set (read once).
pub(crate) fn enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("RUSOLVE_LS_STATS").is_ok_and(|v| v == "1"))
}
