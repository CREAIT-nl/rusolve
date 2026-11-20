//! rusolve (not HiGHS): wall-clock profile of the MIP/LP presolve main loop behind env
//! `RUSOLVE_PRESOLVE_PROF`. Debug output only, never read by the
//! search: with the variable unset [`scope`] takes no clock reading.
//!
//! Phases are the calls of `HPresolve::h16_presolve`'s main loop; the times are inclusive
//! (the fast loop called from inside a phase's follow-up is counted under `fast_loop`).

use std::sync::atomic::{AtomicU64, Ordering};

pub const PHASES: [&str; 17] = [
    "initial",
    "fast_loop",
    "fourier_motzkin",
    "aggregator",
    "sparsify",
    "parallel",
    "strengthen",
    "dominated_cols",
    "enumerate",
    "probing",
    "dependent_eq",
    "donor",
    "probing_init",
    "probing_final",
    "clique_extract",
    "clique_merge",
    "domain_propagate",
];
pub const INITIAL: usize = 0;
pub const FAST_LOOP: usize = 1;
pub const FOURIER_MOTZKIN: usize = 2;
pub const AGGREGATOR: usize = 3;
pub const SPARSIFY: usize = 4;
pub const PARALLEL: usize = 5;
pub const STRENGTHEN: usize = 6;
pub const DOMINATED_COLS: usize = 7;
pub const ENUMERATE: usize = 8;
pub const PROBING: usize = 9;
pub const DEPENDENT_EQ: usize = 10;
pub const DONOR: usize = 11;
// Sub-phases of probing and enumeration (inclusive, also counted in their parent).
pub const PROBING_INIT: usize = 12;
pub const PROBING_FINAL: usize = 13;
pub const CLIQUE_EXTRACT: usize = 14;
pub const CLIQUE_MERGE: usize = 15;
pub const DOMAIN_PROPAGATE: usize = 16;

static NS: [AtomicU64; 17] = [const { AtomicU64::new(0) }; 17];
static CALLS: [AtomicU64; 17] = [const { AtomicU64::new(0) }; 17];

pub fn on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("RUSOLVE_PRESOLVE_PROF").is_some())
}

pub struct Guard(usize, Option<std::time::Instant>);

impl Drop for Guard {
    fn drop(&mut self) {
        if let Some(t) = self.1 {
            NS[self.0].fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
            CALLS[self.0].fetch_add(1, Ordering::Relaxed);
        }
    }
}

/// Times the enclosing scope under `phase` when profiling is on.
pub fn scope(phase: usize) -> Guard {
    Guard(phase, on().then(std::time::Instant::now))
}

/// One line on stderr: `presolve_prof <tag> phase=seconds/calls ...`, then the counters reset.
pub fn report(tag: &str) {
    if !on() {
        return;
    }
    let mut line = format!("presolve_prof {tag}");
    for (i, name) in PHASES.iter().enumerate() {
        let ns = NS[i].swap(0, Ordering::Relaxed);
        let calls = CALLS[i].swap(0, Ordering::Relaxed);
        if calls > 0 {
            line.push_str(&format!(" {name}={:.3}/{calls}", ns as f64 * 1e-9));
        }
    }
    eprintln!("{line}");
}
