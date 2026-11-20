//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! mode_dispatch.cpp `run_presolve` + `run_sequential`, the
//! sequential (`num_workers = 1`) presolve chain of the FPR and LocalMIP
//! entries. Reached only from G5 in mip_solver.rs, with a flag on.

use crate::lp_data::highs_options::HighsOptions;
use crate::mip::fpr_local_mip::context::{make_exec, make_problem};
use crate::mip::fpr_local_mip::rng::{
    heuristic_effort_budget, make_budget, patience_threshold, saturating_mul,
};
use crate::mip::fpr_local_mip::sink::{EffortLedger, S5Sink};
use crate::mip::fpr_local_mip::stats::{self, HeurCounters};
use crate::mip::fpr_local_mip::{
    DispatchOutcome, ExecutionContext, HeuristicBudget, RestartPoint, S5_FPR_EFFORT,
    S5_FPR_PATIENCE, S5_LOCAL_MIP_EFFORT, S5_LOCAL_MIP_PATIENCE, S5_LOCAL_MIP_RESTART,
    S5_LOCAL_MIP_RESTART_WEIGHTS, S5_SKIP_SUBMIP_CHAIN, S5Problem, S5Report, fpr, local_mip,
};
use crate::mip::mip_solver::HighsMipSolver;
use crate::mip::mip_solver_data::{
    kSolutionSourceFPR, kSolutionSourceHeuristic, kSolutionSourceLocalMIP,
};

/// The chain runner signature. FPR reaches the solver through the sink,
/// so both entries share LocalMIP's `without_solver` shape.
type ChainRun =
    fn(&S5Problem, &HeuristicBudget, &ExecutionContext, &mut S5Sink<'_>) -> DispatchOutcome;

/// One entry of the presolve chain. The option pointers of the C++ record
/// are the §0 constants and the bools.
///
/// mip-heuristics (36f7d9d): mode_dispatch.cpp::HeuristicConfig (lines 60-107)
struct HeuristicConfig {
    name: &'static str,
    source_tag: i32,
    effort: f64,
    enable_switch: fn(&HighsOptions) -> bool,
    budget_is_per_worker: bool,
    patience: f64,
    restart: Option<u8>,
    restart_weights: Option<bool>,
    run: ChainRun,
    /// rusolve throughput counters (stats.rs), not in the C++ record.
    stats: &'static HeurCounters,
}

/// mip-heuristics (36f7d9d): mode_dispatch.cpp::kChain (lines 120-134).
/// The "fj" entry (HiGHS keeps its own FJ, Q5) and the "scylla" entry are
/// out of scope (§0).
const K_CHAIN: [HeuristicConfig; 2] = [
    HeuristicConfig {
        name: "fpr",
        source_tag: kSolutionSourceFPR,
        effort: S5_FPR_EFFORT,
        enable_switch: |options| options.mip_heuristic_run_fpr,
        budget_is_per_worker: false,
        patience: S5_FPR_PATIENCE,
        restart: None,
        restart_weights: None,
        run: fpr::run,
        stats: &stats::FPR,
    },
    HeuristicConfig {
        name: "local_mip",
        source_tag: kSolutionSourceLocalMIP,
        effort: S5_LOCAL_MIP_EFFORT,
        enable_switch: |options| options.mip_heuristic_run_local_mip,
        budget_is_per_worker: false,
        patience: S5_LOCAL_MIP_PATIENCE,
        restart: Some(S5_LOCAL_MIP_RESTART),
        restart_weights: Some(S5_LOCAL_MIP_RESTART_WEIGHTS),
        run: local_mip::run,
        stats: &stats::LOCAL_MIP,
    },
];

/// Whether `h` runs at all under `options` (D1: the effort constant and the
/// bool, which takes the place of upstream's enable switch).
///
/// mip-heuristics (36f7d9d): mode_dispatch.cpp::entry_enabled (lines 155-160)
fn entry_enabled(options: &HighsOptions, h: &HeuristicConfig) -> bool {
    if h.effort <= 0.0 {
        return false;
    }
    (h.enable_switch)(options)
}

/// `static_cast<RestartPoint>(options.*h.restart)`; the record bounds the
/// option to the enum's range.
fn restart_point_from(v: u8) -> RestartPoint {
    match v {
        1 => RestartPoint::kBest,
        2 => RestartPoint::kHybrid,
        _ => {
            debug_assert!(v == 0);
            RestartPoint::kDefault
        }
    }
}

/// Each enabled heuristic runs in turn, with its own effort budget, through
/// one shared sink.
///
/// mip-heuristics (36f7d9d): mode_dispatch.cpp::run_sequential (lines 170-311)
fn run_sequential(ms: &mut HighsMipSolver) -> bool {
    // the report is created here if absent and never cleared.
    if ms.mipdata_.s5_report.is_none() {
        ms.mipdata_.s5_report = Some(Box::<S5Report>::default());
    }

    let any_chain_enabled = K_CHAIN.iter().any(|h| entry_enabled(&ms.options_mip_, h));
    if !any_chain_enabled {
        return false;
    }
    let exec = make_exec(ms);

    // exec.terminated() before the sink exists: the same predicate as
    // ExecutionContext::terminated, the terminator lambda of
    // highs_context.h:53 or the deadline.
    if ms.terminator_terminated() || exec.past_deadline() {
        return false;
    }

    let ledger = EffortLedger;

    let problem = make_problem(ms);

    let mut sink = S5Sink::new_incumbent(ms, kSolutionSourceHeuristic);

    for h in K_CHAIN.iter() {
        if !entry_enabled(&sink.mipsolver().unwrap().options_mip_, h) || exec.terminated(&sink) {
            continue;
        }
        let sized = heuristic_effort_budget(problem.nnz, h.effort);
        let total = if h.budget_is_per_worker {
            saturating_mul(sized, exec.num_workers)
        } else {
            sized
        };
        let patience_option = h.patience;
        let patience_per_base = if h.budget_is_per_worker {
            patience_option * (exec.num_workers as f64)
        } else {
            patience_option
        };
        let mut slice = make_budget(
            total,
            exec.num_workers,
            patience_threshold(problem.nnz, patience_per_base, total),
        );
        if let Some(restart) = h.restart {
            slice.restart = restart_point_from(restart);
        }
        if let Some(restart_weights) = h.restart_weights {
            slice.randomize_restart_weights = restart_weights;
        }
        sink.set_source(h.source_tag);

        // run_and_charge (mode_dispatch.cpp:232-256)
        let accepted_before = sink.accepted();
        let t0_s = ledger.now_s(sink.mipsolver().unwrap());
        let t0 = std::time::Instant::now();
        let outcome = (h.run)(&problem, &slice, &exec, &mut sink);
        h.stats.finish_call(t0, outcome.effort as u64);
        let found = sink.accepted() > accepted_before;
        let t1_s = ledger.now_s(sink.mipsolver().unwrap());
        ledger.charge_presolve(
            sink.mipsolver_mut().unwrap(),
            h.name,
            outcome.effort,
            found,
            t0_s,
            t1_s,
            outcome.abandoned_setup,
        );

        // rusolve (mip_ls_int_lp_repair, `ls_int_lp_repair.rs`). LP-polish the
        // pool's best point after an entry that had an offer accepted.
        if found && sink.mipsolver().unwrap().options_mip_.mip_ls_int_lp_repair {
            let mut best: Vec<f64> = Vec::new();
            if sink.copy_best(&mut best) {
                let _ = crate::mip::ls_int_lp_repair::repair(
                    sink.mipsolver_mut().unwrap(),
                    &best,
                    h.source_tag,
                    true,
                );
            }
        }
    }

    false
}

/// The presolve heuristic chain. Returns true when it proved the model
/// infeasible (never, on this path).
///
/// mip-heuristics (36f7d9d): mode_dispatch.cpp::run_presolve (lines 330-370).
/// The fj_effort warning (lines 338-346 comment, 355-367) concerns the
/// excluded FJ entry and is not ported.
pub(crate) fn s5_run_presolve(ms: &mut HighsMipSolver) -> bool {
    if ms.submip && S5_SKIP_SUBMIP_CHAIN {
        return false;
    }
    run_sequential(ms)
}
