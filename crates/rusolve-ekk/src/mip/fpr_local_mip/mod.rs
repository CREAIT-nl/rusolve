//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! The sequential (`num_workers = 1`), effort-budgeted FPR and LocalMIP
//! presolve chain of mip-heuristics (`heuristics::run_presolve`), gated by
//! `mip_heuristic_run_fpr` / `mip_heuristic_run_local_mip`.
//!
//! This file holds the constants and the types (fields,
//! derives and `Default` with the C++ initialisers; no logic). Methods live
//! in the submodules.
//!
//! Deadline/time_limit: `Deadline` is ported with no timer, so it never
//! expires; every poll stays in the code as a literal port. With a finite
//! `time_limit` the chain can therefore overrun it, where upstream would stop.
//!
//! Pointer fields of the C++ structs become borrowed slices with a lifetime
//! (`&'a [T]`), and `nullptr`-able pointers become `Option<&'a ...>`.

pub(crate) mod clique_cover;
pub(crate) mod context;
pub(crate) mod driver;
pub(crate) mod fpr;
pub(crate) mod fpr_core;
pub(crate) mod fpr_val_select;
pub(crate) mod fpr_var_order;
pub(crate) mod gub_swap;
pub(crate) mod lm_construction;
pub(crate) mod lm_core;
pub(crate) mod lm_search;
pub(crate) mod lm_worker;
pub(crate) mod local_mip;
pub(crate) mod pool;
pub(crate) mod prop_engine;
pub(crate) mod repair_search;
pub(crate) mod repair_walk;
pub(crate) mod rng;
pub(crate) mod sink;
pub(crate) mod stats;
#[cfg(test)]
mod tests_e2e;
pub(crate) mod walksat;

pub(crate) use prop_engine::{HeapEntry, IndexedMinHeap, PropEngine, PropResult};

use crate::prelude::*;

use std::rc::Rc;

// ---------------------------------------------------------------------------
// §0 constants
// ---------------------------------------------------------------------------

/// mip-heuristics apply_patch.cmake:502 (mip_heuristic_fpr_effort default)
pub(crate) const S5_FPR_EFFORT: f64 = 3.161;
/// mip-heuristics apply_patch.cmake:503 (mip_heuristic_local_mip_effort default)
pub(crate) const S5_LOCAL_MIP_EFFORT: f64 = 3.2865;
/// mip-heuristics apply_patch.cmake:507 (mip_heuristic_fpr_patience default)
pub(crate) const S5_FPR_PATIENCE: f64 = 0.3372;
/// mip-heuristics apply_patch.cmake:508 (mip_heuristic_local_mip_patience default)
pub(crate) const S5_LOCAL_MIP_PATIENCE: f64 = 3.1943;
/// mip-heuristics apply_patch.cmake:511 (mip_heuristic_skip_submip_chain default, Q6)
pub(crate) const S5_SKIP_SUBMIP_CHAIN: bool = false;
/// mip-heuristics apply_patch.cmake:513 (mip_heuristic_local_mip_restart default; 0 = RestartPoint::kDefault)
pub(crate) const S5_LOCAL_MIP_RESTART: u8 = 0;
/// mip-heuristics apply_patch.cmake:514 (mip_heuristic_local_mip_restart_weights default)
pub(crate) const S5_LOCAL_MIP_RESTART_WEIGHTS: bool = false;
/// mip-heuristics mode_dispatch.cpp:124 (kChain fpr entry, budget_is_per_worker)
pub(crate) const S5_FPR_BUDGET_IS_PER_WORKER: bool = false;
/// mip-heuristics mode_dispatch.cpp:126-130 (kChain local_mip entry, budget_is_per_worker)
pub(crate) const S5_LOCAL_MIP_BUDGET_IS_PER_WORKER: bool = false;

/// mip-heuristics heuristic_common.h:133 (kInfBoundShiftWindow)
pub(crate) const K_INF_BOUND_SHIFT_WINDOW: f64 = 64.0;
/// mip-heuristics heuristic_common.h:141 (kSafeInt64DoubleRange)
pub(crate) const K_SAFE_INT64_DOUBLE_RANGE: f64 = 1e18;
/// mip-heuristics heuristic_common.h:152 (kBaseSeedOffset)
pub(crate) const K_BASE_SEED_OFFSET: u32 = 42;
/// mip-heuristics heuristic_common.h:153 (kSeedStride)
pub(crate) const K_SEED_STRIDE: u32 = 997;
/// mip-heuristics heuristic_common.h:176 (kBudgetBaseShift)
pub(crate) const K_BUDGET_BASE_SHIFT: i32 = 10;
/// mip-heuristics heuristic_common.h:270 (kPatienceCeilingDivisor)
pub(crate) const K_PATIENCE_CEILING_DIVISOR: usize = 4;
/// mip-heuristics heuristic_context.h:666 (make_until_stopped_budget kUnbounded)
pub(crate) const kUnbounded: usize = usize::MAX;

/// mip-heuristics fpr_core.cpp:39
pub(crate) const kDeadlinePollNodes: HighsInt = 16;
/// mip-heuristics fpr_core.cpp:95 (boxed_domain)
pub(crate) const kBox: f64 = 1e5;
/// mip-heuristics prop_engine.cpp:34
pub(crate) const kPropagateBudgetPerNnz: usize = 100;
/// mip-heuristics prop_engine.cpp:61
pub(crate) const kPropagateDeadlinePollWork: usize = 256;
/// mip-heuristics repair_walk.cpp:20
pub(crate) const kTabuLength: HighsInt = 3;
/// mip-heuristics repair_walk.cpp:25
pub(crate) const kSoftRestartPeriod: HighsInt = 10;
/// mip-heuristics repair_walk.cpp:47
pub(crate) const kRepairWalkBudgetPerNnz: usize = 100;
/// mip-heuristics walksat.cpp:31
pub(crate) const kWalkSatBudgetPerNnz: usize = 100;

/// mip-heuristics solution_pool.h:12
pub(crate) const kPoolCapacity: i32 = 10;
/// mip-heuristics solution_pool.h:16
pub(crate) const kDiversityObjTolerance: f64 = 0.10;
/// mip-heuristics solution_pool.h:18
pub(crate) const kDiversityMinHammingFrac: f64 = 0.05;
/// mip-heuristics solution_pool.h:31
pub(crate) const kImprovementObjMargin: f64 = 1e-9;

/// mip-heuristics local_mip_caches.h:58
pub(crate) const kScoreTol: f64 = 5e-7;
/// mip-heuristics local_mip_caches.h:63
pub(crate) const kViolDeltaTol: f64 = 5e-7;
/// mip-heuristics local_mip_caches.h:64
pub(crate) const kRestartInterval: HighsInt = 200000;
/// mip-heuristics local_mip_caches.h:80
pub(crate) const kTermCheckWork: usize = 65536;
/// mip-heuristics local_mip_caches.h:81
pub(crate) const kActivityPeriod: HighsInt = 100000;
/// mip-heuristics local_mip_caches.h:82
pub(crate) const kSmoothProb: f64 = 3e-4;
/// mip-heuristics local_mip_caches.h:83
pub(crate) const kBmsConstraints: HighsInt = 12;
/// mip-heuristics local_mip_caches.h:84
pub(crate) const kBmsBudget: HighsInt = 2250;
/// mip-heuristics local_mip_caches.h:85
pub(crate) const kBmsSatCon: HighsInt = 1;
/// mip-heuristics local_mip_caches.h:86
pub(crate) const kBmsSatBudget: HighsInt = 80;
/// mip-heuristics local_mip_caches.h:87
pub(crate) const kBoolFlipBudget: HighsInt = 5000;
/// mip-heuristics local_mip_caches.h:88
pub(crate) const kEasyBudget: HighsInt = 5;
/// mip-heuristics local_mip_caches.h:89
pub(crate) const kTabuBase: HighsInt = 3;
/// mip-heuristics local_mip_caches.h:90
pub(crate) const kTabuVar: HighsInt = 10;
/// mip-heuristics local_mip_caches.h:91
pub(crate) const kFeasibleRecheckPeriod: HighsInt = 100;
/// mip-heuristics local_mip_caches.h:104
pub(crate) const kFeasiblePlateau: HighsInt = 5000;
/// mip-heuristics local_mip_caches.h:105
pub(crate) const kFeasibleMaxRandomWalks: HighsInt = 20;
/// mip-heuristics local_mip_caches.h:106
pub(crate) const kEpsZero: f64 = 1e-15;
/// mip-heuristics local_mip_worker.h:18
pub(crate) const kPerturbBinaryFraction: f64 = 0.2;
/// mip-heuristics local_mip_construction.h:53
pub(crate) const kConstructionEffortFraction: f64 = 0.10;
/// mip-heuristics continuous_loop.h:63 (ContinuousLoopState::kNoPoller)
pub(crate) const kNoPoller: i32 = -1;

// ---------------------------------------------------------------------------
// heuristic_common.h / deadline.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): heuristic_common.h::VarState (lines 16-19)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct VarState {
    pub(crate) lb: f64,
    pub(crate) ub: f64,
    pub(crate) val: f64,
    pub(crate) fixed: bool,
}

/// mip-heuristics (36f7d9d): heuristic_common.h::CscMatrix (lines 57-61)
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct CscMatrix {
    pub(crate) col_start: Vec<HighsInt>,
    pub(crate) col_row: Vec<HighsInt>,
    pub(crate) col_val: Vec<f64>,
}

/// mip-heuristics (36f7d9d): deadline.h::Deadline (lines 21-55)
///
/// `timer` is always `None`: `Infallible` makes a `Some` unrepresentable,
/// so `expired()` is always false and `remaining()` is always `kHighsInf`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Deadline {
    pub(crate) timer: Option<std::convert::Infallible>,
    pub(crate) limit: f64,
}

// ---------------------------------------------------------------------------
// heuristic_context.h / highs_context.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): heuristic_context.h::DispatchOutcome (lines 102-122)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct DispatchOutcome {
    pub(crate) effort: usize,
    pub(crate) abandoned_setup: bool,
}

/// mip-heuristics (36f7d9d): heuristic_context.h::ProblemView (lines 133-238),
/// as the owned snapshot built by highs_context.h::make_problem (lines 62-109),
/// with ProblemStorage's CSC held in `csc` (D3, R6).
#[derive(Clone, Debug, Default)]
pub(crate) struct S5Problem {
    pub(crate) col_lower: Vec<f64>,
    pub(crate) col_upper: Vec<f64>,
    pub(crate) row_lower: Vec<f64>,
    pub(crate) row_upper: Vec<f64>,
    pub(crate) col_cost: Vec<f64>,
    pub(crate) offset: f64,
    pub(crate) integrality: Vec<HighsVarType>,
    pub(crate) ar_start: Vec<HighsInt>,
    pub(crate) ar_index: Vec<HighsInt>,
    pub(crate) ar_value: Vec<f64>,
    pub(crate) csc: CscMatrix,
    pub(crate) uplocks: Vec<HighsInt>,
    pub(crate) downlocks: Vec<HighsInt>,
    pub(crate) feastol: f64,
    pub(crate) epsilon: f64,
    pub(crate) ncol: HighsInt,
    pub(crate) nrow: HighsInt,
    pub(crate) nnz: usize,
    pub(crate) incumbent: Vec<f64>,
    pub(crate) binary: Vec<u8>,
}

/// mip-heuristics (36f7d9d): heuristic_context.h::RestartPoint (line 319)
#[repr(u8)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)] // C++ k-prefixed enumerator names
pub(crate) enum RestartPoint {
    #[default]
    kDefault = 0,
    kBest = 1,
    kHybrid = 2,
}

/// mip-heuristics (36f7d9d): heuristic_context.h::HeuristicBudget (lines 324-393)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct HeuristicBudget {
    pub(crate) total: usize,
    pub(crate) per_worker: usize,
    pub(crate) attempt_cap: usize,
    pub(crate) stale: usize,
    pub(crate) worker_stale: usize,
    pub(crate) fj_stall_ends_run: bool,
    pub(crate) restart: RestartPoint,
    pub(crate) randomize_restart_weights: bool,
}

/// mip-heuristics (36f7d9d): heuristic_context.h::ExecutionContext (lines 409-512)
///
/// Only the fields the sequential HiGHS path reads: `log_options` and
/// `terminator` are reached through the sink, and `stop`, `least_violated`
/// and `restart_now` are null on that path (highs_context.h::make_exec).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ExecutionContext {
    pub(crate) num_workers: usize,
    pub(crate) base_seed: u32,
    pub(crate) deadline: Deadline,
    /// rusolve GUB swap (`gub_swap.rs`), not a mip-heuristics field:
    /// `mip_ls_gub_swap`, read by `LocalMipWorker::new`.
    pub(crate) gub_swap: bool,
}

impl Default for ExecutionContext {
    // heuristic_context.h:409-435 has no initialisers for num_workers/base_seed;
    // num_workers = 1 is the fixed value, deadline is Deadline{}.
    fn default() -> Self {
        ExecutionContext {
            num_workers: 1,
            base_seed: 0,
            deadline: Deadline::default(),
            gub_swap: false,
        }
    }
}

// ---------------------------------------------------------------------------
// worker_base.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): worker_base.h::AttemptResult (lines 14-17)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct AttemptResult {
    pub(crate) effort: usize,
    pub(crate) found_improvement: bool,
}

/// mip-heuristics (36f7d9d): worker_base.h::WorkerTrace (lines 38-49)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct WorkerTrace {
    pub(crate) worker: i32,
    pub(crate) effort_base: usize,
}

impl Default for WorkerTrace {
    // worker_base.h:42-44
    fn default() -> Self {
        WorkerTrace {
            worker: -1,
            effort_base: 0,
        }
    }
}

/// mip-heuristics (36f7d9d): worker_base.h::WorkerBudgetState (lines 80-131)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct WorkerBudgetState {
    pub(crate) total_budget: usize,
    pub(crate) stale_budget: usize,
    pub(crate) total_effort: usize,
    pub(crate) effort_since_improvement: usize,
    pub(crate) finished: bool,
}

// ---------------------------------------------------------------------------
// fpr_strategies.h / fpr_val_select.h / fpr.cpp
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr_strategies.h::VarStrategy (lines 19-28)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum VarStrategy {
    kLR,
    kType,
    kRandom,
    kLocks,
    kTypecl,
    kCliques,
    kCliques2,
    kDomainSize,
}

/// mip-heuristics (36f7d9d): fpr_val_select.h::ValStrategy (lines 10-20)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)] // C++ k-prefixed enumerator names
pub(crate) enum ValStrategy {
    kUp,
    kRandom,
    kGoodobj,
    kBadobj,
    kLoosedyn,
    kZerocore,
    kZerolp,
    kCore,
    kLp,
}

/// mip-heuristics (36f7d9d): fpr_strategies.h::FprStrategyConfig (lines 39-42)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct FprStrategyConfig {
    pub(crate) var_strategy: VarStrategy,
    pub(crate) val_strategy: ValStrategy,
}

impl Default for FprStrategyConfig {
    // fpr_strategies.h:40-41
    fn default() -> Self {
        FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        }
    }
}

/// mip-heuristics (36f7d9d): fpr_strategies.h::FrameworkMode (lines 44-50)
///
/// No `Default`: the C++ enum has none; FprConfig's member initialiser
/// (kDiveprop) is written out in `FprConfig::default`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)] // C++ k-prefixed enumerator names
pub(crate) enum FrameworkMode {
    kDfs,
    kDfsrep,
    kDive,
    kDiveprop,
    kRepairSearch,
}

/// mip-heuristics fpr_strategies.h:88
pub(crate) const kStratRandom: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kTypecl,
    val_strategy: ValStrategy::kRandom,
};
/// mip-heuristics fpr_strategies.h:89
pub(crate) const kStratRandom2: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kRandom,
    val_strategy: ValStrategy::kRandom,
};
/// mip-heuristics fpr_strategies.h:90
pub(crate) const kStratBadobj: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kType,
    val_strategy: ValStrategy::kBadobj,
};
/// mip-heuristics fpr_strategies.h:91
pub(crate) const kStratBadobjcl: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kTypecl,
    val_strategy: ValStrategy::kBadobj,
};
/// mip-heuristics fpr_strategies.h:92
pub(crate) const kStratGoodobj: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kType,
    val_strategy: ValStrategy::kGoodobj,
};
/// mip-heuristics fpr_strategies.h:93
pub(crate) const kStratGoodobjcl: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kTypecl,
    val_strategy: ValStrategy::kGoodobj,
};
/// mip-heuristics fpr_strategies.h:94
pub(crate) const kStratLocks: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kLR,
    val_strategy: ValStrategy::kLoosedyn,
};
/// mip-heuristics fpr_strategies.h:95
pub(crate) const kStratLocks2: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kLocks,
    val_strategy: ValStrategy::kLoosedyn,
};
/// mip-heuristics fpr_strategies.h:96
pub(crate) const kStratCliques: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kCliques,
    val_strategy: ValStrategy::kUp,
};
/// mip-heuristics fpr_strategies.h:97
pub(crate) const kStratCliques2: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kCliques2,
    val_strategy: ValStrategy::kUp,
};
/// mip-heuristics fpr_strategies.h:98
pub(crate) const kStratDomsize: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kDomainSize,
    val_strategy: ValStrategy::kLoosedyn,
};
/// mip-heuristics fpr_strategies.h:101
pub(crate) const kStratZerocore: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kTypecl,
    val_strategy: ValStrategy::kZerocore,
};
/// mip-heuristics fpr_strategies.h:102
pub(crate) const kStratZerolp: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kTypecl,
    val_strategy: ValStrategy::kZerolp,
};
/// mip-heuristics fpr_strategies.h:103
pub(crate) const kStratCore: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kTypecl,
    val_strategy: ValStrategy::kCore,
};
/// mip-heuristics fpr_strategies.h:104
pub(crate) const kStratLp: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kTypecl,
    val_strategy: ValStrategy::kLp,
};
/// mip-heuristics fpr_strategies.h:119
pub(crate) const kStratLocksLp: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kLocks,
    val_strategy: ValStrategy::kLp,
};
/// mip-heuristics fpr_strategies.h:120
pub(crate) const kStratLRLp: FprStrategyConfig = FprStrategyConfig {
    var_strategy: VarStrategy::kLR,
    val_strategy: ValStrategy::kLp,
};

/// mip-heuristics (36f7d9d): fpr_strategies.h::NamedConfig (lines 126-129)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct NamedConfig {
    pub(crate) strat: FprStrategyConfig,
    pub(crate) mode: FrameworkMode,
}

/// mip-heuristics fpr.cpp:158-169 (kFprStrategies)
pub(crate) const kFprStrategies: [FprStrategyConfig; 8] = [
    kStratBadobjcl,
    kStratLocks2,
    kStratLocks,
    kStratRandom,
    kStratDomsize,
    kStratRandom2,
    kStratBadobj,
    kStratGoodobj,
];
/// mip-heuristics fpr.cpp:170 (kNumFprStrategies)
pub(crate) const kNumFprStrategies: i32 = kFprStrategies.len() as i32;

/// mip-heuristics (36f7d9d): fpr.cpp::InitialFprConfig (lines 175-178)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct InitialFprConfig {
    pub(crate) strat_idx: i32,
    pub(crate) mode: FrameworkMode,
}

/// mip-heuristics fpr.cpp:179-188 (kInitialFprConfigs)
pub(crate) const kInitialFprConfigs: [InitialFprConfig; 8] = [
    InitialFprConfig {
        strat_idx: 0,
        mode: FrameworkMode::kDfs,
    },
    InitialFprConfig {
        strat_idx: 1,
        mode: FrameworkMode::kDfs,
    },
    InitialFprConfig {
        strat_idx: 1,
        mode: FrameworkMode::kDive,
    },
    InitialFprConfig {
        strat_idx: 2,
        mode: FrameworkMode::kDfsrep,
    },
    InitialFprConfig {
        strat_idx: 0,
        mode: FrameworkMode::kDfsrep,
    },
    InitialFprConfig {
        strat_idx: 3,
        mode: FrameworkMode::kDiveprop,
    },
    InitialFprConfig {
        strat_idx: 2,
        mode: FrameworkMode::kRepairSearch,
    },
    InitialFprConfig {
        strat_idx: 4,
        mode: FrameworkMode::kDfs,
    },
];
/// mip-heuristics fpr.cpp:189 (kNumInitialFprConfigs)
pub(crate) const kNumInitialFprConfigs: i32 = kInitialFprConfigs.len() as i32;

// ---------------------------------------------------------------------------
// walksat.h / repair_walk.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): walksat.h::WalkSatMove (lines 11-14)
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct WalkSatMove {
    pub(crate) var: HighsInt,
    pub(crate) val: f64,
}

impl Default for WalkSatMove {
    // walksat.h:12-13
    fn default() -> Self {
        WalkSatMove { var: -1, val: 0.0 }
    }
}

/// mip-heuristics (36f7d9d): walksat.h::WalkSatScratch::Candidate (lines 24-28)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct WalkSatCandidate {
    pub(crate) var: HighsInt,
    pub(crate) new_val: f64,
    pub(crate) damage: f64,
}

/// mip-heuristics (36f7d9d): walksat.h::WalkSatScratch::UndoEntry (lines 33-36)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct WalkSatUndoEntry {
    pub(crate) idx: HighsInt,
    pub(crate) old_val: f64,
}

/// mip-heuristics (36f7d9d): walksat.h::WalkSatScratch (lines 23-41)
#[derive(Clone, Debug, Default)]
pub(crate) struct WalkSatScratch {
    pub(crate) cand: Vec<WalkSatCandidate>,
    pub(crate) best_indices: Vec<HighsInt>,
    pub(crate) violated: Vec<HighsInt>,
    /// size = nrow, -1 = absent
    pub(crate) violated_pos: Vec<HighsInt>,
    pub(crate) sol_undo: Vec<WalkSatUndoEntry>,
    pub(crate) lhs_undo: Vec<WalkSatUndoEntry>,
}

/// mip-heuristics (36f7d9d): repair_walk.h::RepairWalkScratch::Candidate (lines 29-33)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct RepairWalkCandidate {
    pub(crate) var: HighsInt,
    pub(crate) shift: f64,
    pub(crate) damage: f64,
}

/// mip-heuristics (36f7d9d): repair_walk.h::RepairWalkScratch (lines 24-46)
#[derive(Clone, Debug, Default)]
pub(crate) struct RepairWalkScratch {
    pub(crate) cand: Vec<RepairWalkCandidate>,
    pub(crate) best_indices: Vec<HighsInt>,
    pub(crate) violated: Vec<HighsInt>,
    pub(crate) violated_pos: Vec<HighsInt>,
    pub(crate) tabu: Vec<HighsInt>,
}

// ---------------------------------------------------------------------------
// fpr_core.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr_core.h::FprDfsNode (lines 22-30)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct FprDfsNode {
    pub(crate) var: HighsInt,
    pub(crate) val: f64,
    pub(crate) vs_mark: HighsInt,
    pub(crate) sol_mark: HighsInt,
    pub(crate) act_mark: HighsInt,
    pub(crate) pq_mark: HighsInt,
    pub(crate) cursor_reset: HighsInt,
}

/// mip-heuristics (36f7d9d): fpr_core.h::RepairSearchNode (lines 35-58)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct RepairSearchNode {
    /// variable to branch on (-1 for root)
    pub(crate) var: HighsInt,
    /// fix value or bound value
    pub(crate) val: f64,
    /// true = fix(var, val), false = tighten bound
    pub(crate) is_fix: bool,
    /// if !is_fix: true = tighten_lb, false = tighten_ub
    pub(crate) is_lb: bool,
    pub(crate) e_vs_mark: HighsInt,
    pub(crate) e_sol_mark: HighsInt,
    /// E PQ mark (-1 if pq not initialized)
    pub(crate) e_pq_mark: HighsInt,
    pub(crate) r_vs_mark: HighsInt,
    pub(crate) r_sol_mark: HighsInt,
    pub(crate) sol_undo_mark: HighsInt,
    pub(crate) lhs_undo_mark: HighsInt,
    /// total violation at parent (for BacktrackBestOpen)
    pub(crate) violation: f64,
}

/// mip-heuristics (36f7d9d): fpr_core.h::FprScratch (lines 64-122)
#[derive(Clone, Debug, Default)]
pub(crate) struct FprScratch<'p> {
    pub(crate) var_order: Vec<HighsInt>,
    pub(crate) dfs_stack: Vec<FprDfsNode>,
    pub(crate) lhs_cache: Vec<f64>,
    pub(crate) solution: Vec<f64>,
    pub(crate) walksat: WalkSatScratch,
    pub(crate) repair_walk: RepairWalkScratch,
    pub(crate) repair_dfs_stack: Vec<RepairSearchNode>,
    pub(crate) repair_best_solution: Vec<f64>,
    pub(crate) repair_best_lhs: Vec<f64>,
    pub(crate) prop_engine: Option<PropEngine<'p>>,
    pub(crate) repair_prop_engine_r: Option<PropEngine<'p>>,
}

/// mip-heuristics (36f7d9d): fpr_core.h::FprConfig (lines 124-228)
///
/// `precomputed_var_order` carries its own length; `precomputed_var_order_size`
/// is kept as the C++ field.
#[derive(Debug)]
pub(crate) struct FprConfig<'a, 'p> {
    pub(crate) max_effort: usize,
    pub(crate) cont_fallback: Option<&'a [f64]>,
    pub(crate) mode: FrameworkMode,
    pub(crate) strategy: Option<&'a FprStrategyConfig>,
    pub(crate) lp_ref: Option<&'a [f64]>,
    pub(crate) precomputed_var_order: Option<&'a [HighsInt]>,
    pub(crate) precomputed_var_order_size: HighsInt,
    pub(crate) repair_noise: f64,
    pub(crate) repair_iterations: HighsInt,
    pub(crate) walksat_iterations: HighsInt,
    pub(crate) repair_track_best: bool,
    pub(crate) scratch: Option<&'a mut FprScratch<'p>>,
}

impl Default for FprConfig<'_, '_> {
    // fpr_core.h:153-227. max_effort and cont_fallback have no C++
    // initialiser; they are value-initialised here (0 / nullptr).
    fn default() -> Self {
        FprConfig {
            max_effort: 0,
            cont_fallback: None,
            mode: FrameworkMode::kDiveprop,
            strategy: None,
            lp_ref: None,
            precomputed_var_order: None,
            precomputed_var_order_size: 0,
            repair_noise: 0.75,
            repair_iterations: 50,
            walksat_iterations: 200,
            repair_track_best: true,
            scratch: None,
        }
    }
}

/// mip-heuristics (36f7d9d): fpr_core.h::FprAttemptState::Phase (lines 300-309)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)] // C++ k-prefixed enumerator names
pub(crate) enum FprPhase {
    #[default]
    kIdle,
    kDfs,
    kReadyToFinish,
}

/// mip-heuristics (36f7d9d): fpr_core.h::FprAttemptState (lines 272-311)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct FprAttemptState {
    pub(crate) ncol: HighsInt,
    pub(crate) nrow: HighsInt,
    pub(crate) attempt_idx: i32,
    pub(crate) dynamic_var: bool,
    pub(crate) do_propagate: bool,
    pub(crate) do_backtrack: bool,
    pub(crate) do_repair: bool,
    pub(crate) node_limit: HighsInt,
    pub(crate) var_order_size: HighsInt,
    pub(crate) var_order_cursor: HighsInt,
    pub(crate) nodes_visited: HighsInt,
    pub(crate) found_complete: bool,
    pub(crate) effort_consumed: usize,
    pub(crate) phase: FprPhase,
}

/// mip-heuristics (36f7d9d): fpr_core.h::FprStepResult (lines 313-321)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FprStepResult {
    kBudgetGate,
    kVerdictReady,
}

// ---------------------------------------------------------------------------
// solution_pool.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): solution_pool.h::SolutionPool::Snapshot (lines 47-50)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct PoolSnapshot {
    pub(crate) has_solution: bool,
    pub(crate) best_objective: f64,
}

/// mip-heuristics (36f7d9d): solution_pool.h::SolutionPool::Entry (lines 54-63)
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct PoolEntry {
    pub(crate) objective: f64,
    pub(crate) solution: Vec<f64>,
    pub(crate) source: i32,
}

/// mip-heuristics (36f7d9d): solution_pool.h::SolutionPool::EntryRef (lines 67-71);
/// `SolutionPtr` (shared_ptr<const vector<double>>) is `Rc<Vec<f64>>`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PoolEntryRef {
    pub(crate) objective: f64,
    pub(crate) source: i32,
    pub(crate) solution: Rc<Vec<f64>>,
}

/// mip-heuristics (36f7d9d): solution_pool.h::SolutionPool::AddResult (lines 101-104)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct PoolAddResult {
    pub(crate) accepted: bool,
    pub(crate) improved_best: bool,
}

/// mip-heuristics (36f7d9d): solution_pool.h::SolutionPool::RestartPolicy (line 154)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum RestartPolicy {
    #[default]
    kMixed,
    kBestOnly,
}

/// mip-heuristics (36f7d9d): solution_pool.h::SolutionPool::IntegerColumns (lines 203-206)
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct IntegerColumns {
    pub(crate) mask: Vec<bool>,
    pub(crate) index: Vec<i32>,
}

/// mip-heuristics (36f7d9d): solution_pool.h::SolutionPool (lines 45-248)
///
/// `mtx_` (HighsSpinMutex) is dropped: the path is single-threaded.
/// `MaskPtr` (shared_ptr<const IntegerColumns>) is `Option<Rc<IntegerColumns>>`.
pub(crate) struct SolutionPool {
    /// sorted best first
    pub(crate) entries: Vec<PoolEntryRef>,
    pub(crate) capacity: i32,
    pub(crate) integer_mask: Option<Rc<IntegerColumns>>,
    pub(crate) generation: u64,
    pub(crate) scan_hook: Option<Box<dyn FnMut()>>,
    pub(crate) best_seen: f64,
    pub(crate) has_best_seen: bool,
}

// ---------------------------------------------------------------------------
// local_mip_caches.h / local_mip_core.h
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): local_mip_caches.h::IndexedSet (lines 115-155)
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct IndexedSet {
    pub(crate) elements: Vec<HighsInt>,
    /// -1 = absent
    pub(crate) pos: Vec<HighsInt>,
}

// mip-heuristics (36f7d9d): local_mip_caches.h::ViolCache (lines 158-181)
// is not ported: compute_candidate_scores recomputes the pre-move
// violation instead (throughput).

/// The row-only part of `WorkerCtx::compute_tight_delta` (lm_core.rs).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct TightRow {
    pub(crate) gap: f64,
    pub(crate) row_violated: bool,
    pub(crate) is_equality: bool,
}

/// mip-heuristics (36f7d9d): local_mip_caches.h::Candidate (lines 184-189)
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Candidate {
    pub(crate) var_idx: HighsInt,
    pub(crate) new_val: f64,
    pub(crate) score: f64,
    pub(crate) bonus: f64,
}

impl Default for Candidate {
    // local_mip_caches.h:185-188
    fn default() -> Self {
        Candidate {
            var_idx: -1,
            new_val: 0.0,
            score: f64::NEG_INFINITY,
            bonus: 0.0,
        }
    }
}

/// mip-heuristics (36f7d9d): local_mip_caches.h::BatchCand (lines 191-194)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct BatchCand {
    pub(crate) var_idx: HighsInt,
    pub(crate) new_val: f64,
}

/// mip-heuristics (36f7d9d): local_mip_caches.h::WeightedCon (lines 196-199)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct WeightedCon {
    pub(crate) ci: HighsInt,
    pub(crate) w: u64,
}

/// mip-heuristics (36f7d9d): local_mip_caches.h::LiftCache (lines 205-237)
///
/// `costed_vars` points at the owning LocalMipWorker's own vector in C++
/// (local_mip_worker.cpp:102); it is shared as `Rc<Vec<HighsInt>>` here
/// (None = nullptr).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct LiftCache {
    pub(crate) lo: Vec<f64>,
    pub(crate) hi: Vec<f64>,
    pub(crate) score: Vec<f64>,
    pub(crate) dirty: Vec<bool>,
    pub(crate) dirty_list: Vec<HighsInt>,
    pub(crate) all_dirty: bool,
    pub(crate) positive_list: Vec<HighsInt>,
    pub(crate) in_positive: Vec<bool>,
    pub(crate) costed_vars: Option<Rc<Vec<HighsInt>>>,
}

impl Default for LiftCache {
    // local_mip_caches.h:206-212
    fn default() -> Self {
        LiftCache {
            lo: Vec::new(),
            hi: Vec::new(),
            score: Vec::new(),
            dirty: Vec::new(),
            dirty_list: Vec::new(),
            all_dirty: true,
            positive_list: Vec::new(),
            in_positive: Vec::new(),
            costed_vars: None,
        }
    }
}

/// mip-heuristics (36f7d9d): local_mip_core.h::WorkerCtx (lines 60-194)
///
/// The `const std::vector<T>&` members borrow the `S5Problem` as `&'a [T]`.
#[derive(Clone, Debug)]
pub(crate) struct WorkerCtx<'a> {
    pub(crate) ar_start: &'a [HighsInt],
    pub(crate) ar_index: &'a [HighsInt],
    pub(crate) ar_value: &'a [f64],
    pub(crate) col_lb: &'a [f64],
    pub(crate) col_ub: &'a [f64],
    pub(crate) col_cost: &'a [f64],
    pub(crate) offset: f64,
    pub(crate) row_lo: &'a [f64],
    pub(crate) row_hi: &'a [f64],
    pub(crate) integrality: &'a [HighsVarType],
    pub(crate) csc: &'a CscMatrix,
    pub(crate) feastol: f64,
    pub(crate) epsilon: f64,
    pub(crate) ncol: HighsInt,
    pub(crate) nrow: HighsInt,
    pub(crate) binary: &'a [u8],
    pub(crate) solution: Vec<f64>,
    pub(crate) lhs: Vec<f64>,
    pub(crate) weight: Vec<u64>,
    /// C++ initialiser: 1 (local_mip_core.h:88)
    pub(crate) obj_weight: u64,
    pub(crate) current_obj: f64,
    pub(crate) violated: IndexedSet,
    pub(crate) satisfied: IndexedSet,
    pub(crate) lift: LiftCache,
    pub(crate) tabu_inc_until: Vec<HighsInt>,
    pub(crate) tabu_dec_until: Vec<HighsInt>,
    pub(crate) batch: Vec<BatchCand>,
    pub(crate) sampled: Vec<WeightedCon>,
    /// C++ initialiser: true (local_mip_core.h:113)
    pub(crate) was_infeasible: bool,
    pub(crate) feasible_recheck_counter: HighsInt,
    pub(crate) effort: usize,
    pub(crate) track_violation: bool,
    pub(crate) total_violation: f64,
    /// rusolve (throughput): bit `j` is set
    /// when costed column `j` can yield a breakthrough candidate, i.e. it is
    /// not parked at the bound its cost favours. Kept by `apply_move` and
    /// `rebuild_state`; read by `infeasible_step` Phase 1b only while
    /// `bt_live_valid`.
    pub(crate) bt_live: Vec<u64>,
    pub(crate) bt_live_valid: bool,
    /// Number of costed columns when `bt_live` was last rebuilt.
    pub(crate) bt_costed: usize,
    /// rusolve throughput counter (stats.rs): moves scored. Never read by the search.
    pub(crate) ls_evals: u64,
    /// rusolve throughput counter (stats.rs): moves applied. Never read by the search.
    pub(crate) ls_commits: u64,
    /// rusolve GUB swap (`gub_swap.rs`), not a mip-heuristics field:
    /// the GUB compound-move operator; `None` unless `mip_ls_gub_swap` is on.
    pub(crate) gub: Option<Box<gub_swap::GubSwap>>,
}

// ---------------------------------------------------------------------------
// test-observable report
// ---------------------------------------------------------------------------

/// one accepted offer (the `[HeurSol]` fields of incumbent_sink.cpp:102-130).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct S5Offer {
    pub(crate) name: &'static str,
    pub(crate) worker: i32,
    pub(crate) effort_at: usize,
    pub(crate) obj: f64,
}

/// one `[Heur]` charge of effort_ledger.cpp:28-50.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct S5HeurCharge {
    pub(crate) name: &'static str,
    pub(crate) effort: usize,
    pub(crate) found: bool,
    pub(crate) abandoned_setup: bool,
}

/// accepted offers and per-heuristic charges of the presolve chain.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct S5Report {
    pub(crate) offers: Vec<S5Offer>,
    pub(crate) heur: Vec<S5HeurCharge>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s5_data_option_constants() {
        assert_eq!(S5_FPR_EFFORT, 3.161);
        assert_eq!(S5_LOCAL_MIP_EFFORT, 3.2865);
        assert_eq!(S5_FPR_PATIENCE, 0.3372);
        assert_eq!(S5_LOCAL_MIP_PATIENCE, 3.1943);
        const { assert!(!S5_SKIP_SUBMIP_CHAIN) };
        assert_eq!(S5_LOCAL_MIP_RESTART, RestartPoint::kDefault as u8);
        const { assert!(!S5_LOCAL_MIP_RESTART_WEIGHTS) };
        const { assert!(!S5_FPR_BUDGET_IS_PER_WORKER) };
        const { assert!(!S5_LOCAL_MIP_BUDGET_IS_PER_WORKER) };
    }

    #[test]
    fn s5_data_common_constants() {
        assert_eq!(K_INF_BOUND_SHIFT_WINDOW, 64.0);
        assert_eq!(K_SAFE_INT64_DOUBLE_RANGE, 1e18);
        assert_eq!(K_BASE_SEED_OFFSET, 42u32);
        assert_eq!(K_SEED_STRIDE, 997u32);
        assert_eq!(K_BUDGET_BASE_SHIFT, 10);
        assert_eq!(K_PATIENCE_CEILING_DIVISOR, 4usize);
        assert_eq!(kUnbounded, usize::MAX);
        assert_eq!(kNoPoller, -1);
    }

    #[test]
    fn s5_data_fpr_constants() {
        assert_eq!(kDeadlinePollNodes, 16);
        assert_eq!(kBox, 1e5);
        assert_eq!(kPropagateBudgetPerNnz, 100);
        assert_eq!(kPropagateDeadlinePollWork, 256);
        assert_eq!(IndexedMinHeap::kNotPresent, -1);
        assert_eq!(kTabuLength, 3);
        assert_eq!(kSoftRestartPeriod, 10);
        assert_eq!(kRepairWalkBudgetPerNnz, 100);
        assert_eq!(kWalkSatBudgetPerNnz, 100);
    }

    #[test]
    fn s5_data_fpr_strategy_tables() {
        use FrameworkMode::*;
        let s = |v, w| FprStrategyConfig {
            var_strategy: v,
            val_strategy: w,
        };
        assert_eq!(
            kFprStrategies,
            [
                s(VarStrategy::kTypecl, ValStrategy::kBadobj),
                s(VarStrategy::kLocks, ValStrategy::kLoosedyn),
                s(VarStrategy::kLR, ValStrategy::kLoosedyn),
                s(VarStrategy::kTypecl, ValStrategy::kRandom),
                s(VarStrategy::kDomainSize, ValStrategy::kLoosedyn),
                s(VarStrategy::kRandom, ValStrategy::kRandom),
                s(VarStrategy::kType, ValStrategy::kBadobj),
                s(VarStrategy::kType, ValStrategy::kGoodobj),
            ]
        );
        assert_eq!(kNumFprStrategies, 8);
        let got: Vec<(i32, FrameworkMode)> = kInitialFprConfigs
            .iter()
            .map(|cv| (cv.strat_idx, cv.mode))
            .collect();
        assert_eq!(
            got,
            vec![
                (0, kDfs),
                (1, kDfs),
                (1, kDive),
                (2, kDfsrep),
                (0, kDfsrep),
                (3, kDiveprop),
                (2, kRepairSearch),
                (4, kDfs),
            ]
        );
        assert_eq!(kNumInitialFprConfigs, 8);
        assert_eq!(
            FprStrategyConfig::default(),
            s(VarStrategy::kType, ValStrategy::kGoodobj)
        );
    }

    #[test]
    fn s5_data_pool_and_local_mip_constants() {
        assert_eq!(kPoolCapacity, 10);
        assert_eq!(kDiversityObjTolerance, 0.10);
        assert_eq!(kDiversityMinHammingFrac, 0.05);
        assert_eq!(kImprovementObjMargin, 1e-9);
        assert_eq!(kScoreTol, 5e-7);
        assert_eq!(kViolDeltaTol, 5e-7);
        assert_eq!(kRestartInterval, 200000);
        assert_eq!(kTermCheckWork, 65536usize);
        assert_eq!(kActivityPeriod, 100000);
        assert_eq!(kSmoothProb, 3e-4);
        assert_eq!(kBmsConstraints, 12);
        assert_eq!(kBmsBudget, 2250);
        assert_eq!(kBmsSatCon, 1);
        assert_eq!(kBmsSatBudget, 80);
        assert_eq!(kBoolFlipBudget, 5000);
        assert_eq!(kEasyBudget, 5);
        assert_eq!(kTabuBase, 3);
        assert_eq!(kTabuVar, 10);
        assert_eq!(kFeasibleRecheckPeriod, 100);
        assert_eq!(kFeasiblePlateau, 5000);
        assert_eq!(kFeasibleMaxRandomWalks, 20);
        assert_eq!(kEpsZero, 1e-15);
        assert_eq!(kPerturbBinaryFraction, 0.2);
        assert_eq!(kConstructionEffortFraction, 0.10);
    }

    #[test]
    fn s5_data_heuristic_budget_default() {
        let b = HeuristicBudget::default();
        assert_eq!(b.total, 0);
        assert_eq!(b.per_worker, 0);
        assert_eq!(b.attempt_cap, 0);
        assert_eq!(b.stale, 0);
        assert_eq!(b.worker_stale, 0);
        assert!(!b.fj_stall_ends_run);
        assert_eq!(b.restart, RestartPoint::kDefault);
        assert!(!b.randomize_restart_weights);
    }

    #[test]
    fn s5_data_fpr_config_default() {
        let cv = FprConfig::default();
        assert_eq!(cv.max_effort, 0);
        assert!(cv.cont_fallback.is_none());
        assert_eq!(cv.mode, FrameworkMode::kDiveprop);
        assert!(cv.strategy.is_none());
        assert!(cv.lp_ref.is_none());
        assert!(cv.precomputed_var_order.is_none());
        assert_eq!(cv.precomputed_var_order_size, 0);
        assert_eq!(cv.repair_noise, 0.75);
        assert_eq!(cv.repair_iterations, 50);
        assert_eq!(cv.walksat_iterations, 200);
        assert!(cv.repair_track_best);
        assert!(cv.scratch.is_none());
    }

    #[test]
    fn s5_data_execution_context_default() {
        let e = ExecutionContext::default();
        assert_eq!(e.num_workers, 1);
        assert_eq!(e.base_seed, 0);
        assert!(e.deadline.timer.is_none());
        assert_eq!(e.deadline.limit, 0.0);
    }

    #[test]
    fn s5_data_small_struct_defaults() {
        assert_eq!(
            DispatchOutcome::default(),
            DispatchOutcome {
                effort: 0,
                abandoned_setup: false
            }
        );
        assert_eq!(
            WorkerTrace::default(),
            WorkerTrace {
                worker: -1,
                effort_base: 0
            }
        );
        assert_eq!(AttemptResult::default().effort, 0);
        assert!(!AttemptResult::default().found_improvement);
        assert_eq!(WorkerBudgetState::default().total_budget, 0);
        assert!(!WorkerBudgetState::default().finished);
        let cv = Candidate::default();
        assert_eq!(cv.var_idx, -1);
        assert_eq!(cv.score, f64::NEG_INFINITY);
        assert!(LiftCache::default().all_dirty);
        assert_eq!(WalkSatMove::default().var, -1);
        assert_eq!(FprAttemptState::default().phase, FprPhase::kIdle);
        assert_eq!(PoolAddResult::default(), PoolAddResult::default());
        assert!(!PoolAddResult::default().accepted);
        assert!(S5Report::default().offers.is_empty());
    }
}
