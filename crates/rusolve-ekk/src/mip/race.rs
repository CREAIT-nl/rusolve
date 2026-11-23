//! rusolve, not a HiGHS file: concurrent racing of diversified
//! copies of the port's MIP search with incumbent sharing.
//!
//! Sources: Shinano et al., "FiberSCIP - a shared memory parallelization of SCIP", INFORMS J.
//! Computing 30(1), 2018 (racing ramp-up: the first racer to finish wins); Carvajal, Ahmed,
//! Nemhauser et al., "Using diversification, communication and parallelism to solve
//! mixed-integer linear programs", Oper. Res. Letters 42(2), 2014 (sharing incumbents between
//! differently-configured searches accounts for most of the gain).
//!
//! Each racer is a full `HighsMipSolver` on its own thread with its own one-thread scheduler.
//! Racers share only this store: a stop flag and, per presolved-model fingerprint, the best
//! incumbent found so far. A racer whose presolved model hashes differently (seed- or
//! cutoff-dependent presolve) neither imports nor is imported; its objective still races.

use crate::lp_data::highs_options::HighsOptions;
use crate::prelude::*;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

#[derive(Debug)]
struct RaceEntry {
    fp: u64,
    obj: f64,
    sol: Vec<f64>,
}

/// State shared by the racers of one solve.
#[derive(Debug, Default)]
pub struct RaceShared {
    /// Set by the first racer to finish with a terminal status; racers stop at `check_limits`.
    pub stop: AtomicBool,
    /// Bumped on every improving publish; racers import only when it moved.
    pub version: AtomicU64,
    best: Mutex<Vec<RaceEntry>>,
}

impl RaceShared {
    /// Offer an incumbent (`obj` and `sol` in the presolved space of fingerprint `fp`).
    pub fn publish(&self, fp: u64, obj: f64, sol: &[f64]) {
        let mut best = self.best.lock().unwrap();
        match best.iter_mut().find(|e| e.fp == fp) {
            Some(e) if obj < e.obj => {
                e.obj = obj;
                e.sol = sol.to_vec();
            }
            Some(_) => return,
            None => best.push(RaceEntry {
                fp,
                obj,
                sol: sol.to_vec(),
            }),
        }
        self.version.fetch_add(1, Ordering::Release);
    }

    /// The best shared incumbent for fingerprint `fp` that beats `below`.
    pub fn fetch(&self, fp: u64, below: f64) -> Option<(f64, Vec<f64>)> {
        let best = self.best.lock().unwrap();
        best.iter()
            .find(|e| e.fp == fp && e.obj < below)
            .map(|e| (e.obj, e.sol.clone()))
    }
}

/// Fingerprint of a presolved model: racers share incumbents only on equal fingerprints.
pub fn model_fingerprint(lp: &HighsLp) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut mix = |x: u64| {
        h ^= x;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
        h ^= h >> 29;
    };
    mix(lp.num_col_ as u64);
    mix(lp.num_row_ as u64);
    for v in lp
        .col_cost_
        .iter()
        .chain(&lp.col_lower_)
        .chain(&lp.col_upper_)
        .chain(&lp.row_lower_)
        .chain(&lp.row_upper_)
        .chain(&lp.a_matrix_.value_)
    {
        mix(v.to_bits());
    }
    for &i in &lp.a_matrix_.index_ {
        mix(i as u64);
    }
    h
}

/// The options of racer `k` (k = 0 is the caller's configuration, unchanged). Each racer moves
/// the search away from the caller's settings in one direction: primal heuristics, cutoff
/// big-M restarts, cut-heavy, or reliability branching depth; seeds differ for all of them.
pub fn racer_options(base: &HighsOptions, k: usize) -> HighsOptions {
    let mut o = base.clone();
    o.threads = 1;
    o.mip_racing_threads = 0;
    if k == 0 {
        return o;
    }
    o.random_seed = base.random_seed.wrapping_add(k as HighsInt);
    o.output_flag = false;
    o.log_to_console = false;
    o.log_options.output_flag = false;
    o.log_options.log_to_console = false;
    o.mip_improving_solution_file = String::new();
    match (k - 1) % 4 {
        0 => {
            o.mip_proximity = true;
            o.mip_alns = true;
        }
        1 => o.mip_cutoff_bigm = true,
        2 => {
            o.mip_cut_selection_hybrid = true;
            o.mip_pool_soft_limit = base.mip_pool_soft_limit.saturating_mul(2);
        }
        _ => {
            o.mip_pscost_minreliable = base.mip_pscost_minreliable.saturating_mul(2).max(4);
            o.mip_heuristic_effort = (base.mip_heuristic_effort * 2.0).min(0.5);
        }
    }
    o
}

/// Run `n` racers on `n` threads; return them in index order after the race is decided.
/// `run_one` builds and runs racer `k` and returns it. The first racer to return with a terminal
/// status sets the stop flag; racer 0 returning for any reason does too (user interrupt).
pub fn run_race<S, F>(n: usize, shared: &Arc<RaceShared>, run_one: F) -> Vec<S>
where
    S: Send,
    F: Fn(usize, &Arc<RaceShared>) -> (S, bool) + Sync,
{
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..n)
            .map(|k| {
                let run_one = &run_one;
                let shared = Arc::clone(shared);
                std::thread::Builder::new()
                    .stack_size(256 << 20)
                    .spawn_scoped(scope, move || {
                        crate::parallel::highs_parallel::initialize_scheduler(1);
                        let (s, terminal) = run_one(k, &shared);
                        if terminal || k == 0 {
                            shared.stop.store(true, Ordering::Release);
                        }
                        s
                    })
                    .expect("spawn racer")
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("racer panicked"))
            .collect()
    })
}
