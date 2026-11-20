//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! fpr_core.cpp + fpr_core.h inline code (boxed_domain, the lifecycle begin/step/finish).

use super::fpr_val_select::choose_value;
use super::prop_engine::PropEngine;
use super::repair_search::{K_REPAIR_PROGRESS_THRESHOLD, repair_search};
use super::repair_walk::{any_violated_row_in_column, repair_walk};
use super::rng::{is_integer, is_row_violated};
use super::walksat::{greedy_1opt, walksat_repair};
use super::{
    CscMatrix, Deadline, FprAttemptState, FprConfig, FprDfsNode, FprPhase, FprScratch,
    FprStepResult, FprStrategyConfig, FrameworkMode, PropResult, S5Problem, ValStrategy,
    VarStrategy, kBox, kDeadlinePollNodes,
};
use crate::HighsInt;
use crate::lp_data::hconst::kHighsInf;
use crate::mip::fpr_local_mip::rng::{S5Rng, s5_shuffle};
use crate::util::highs_utils::{cpp_max, cpp_min};

// ---------------------------------------------------------------------------
// Anonymous-namespace helpers (fpr_core.cpp lines 21-201)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr_core.cpp::AttemptCtx (lines 47-67)
///
/// The per-call references begin/step/finish rehydrate the lambdas from. The
/// C++ holds const refs; here they are `&'p` slices borrowed from the
/// `S5Problem` (which owns the vectors), so the struct borrows nothing else.
struct AttemptCtx<'p> {
    ar_start: &'p [HighsInt],
    ar_index: &'p [HighsInt],
    ar_value: &'p [f64],
    col_lb: &'p [f64],
    col_ub: &'p [f64],
    // The view's minimisation form (#170).
    col_cost: &'p [f64],
    row_lo: &'p [f64],
    row_hi: &'p [f64],
    integrality: &'p [crate::lp_data::hconst::HighsVarType],
    feastol: f64,
    ncol: HighsInt,
    nrow: HighsInt,
    // Dispatch-time `isBinary` snapshot; see `ProblemView::binary`.
    binary: &'p [u8],
    uplocks: &'p [HighsInt],
    downlocks: &'p [HighsInt],
}

impl<'p> AttemptCtx<'p> {
    #[inline]
    fn is_binary(&self, j: HighsInt) -> bool {
        self.binary[j as usize] != 0
    }
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::make_ctx (lines 69-76)
fn make_ctx<'p>(problem: &'p S5Problem) -> AttemptCtx<'p> {
    AttemptCtx {
        ar_start: &problem.ar_start,
        ar_index: &problem.ar_index,
        ar_value: &problem.ar_value,
        col_lb: &problem.col_lower,
        col_ub: &problem.col_upper,
        col_cost: &problem.col_cost,
        row_lo: &problem.row_lower,
        row_hi: &problem.row_upper,
        integrality: &problem.integrality,
        feastol: problem.feastol,
        ncol: problem.ncol,
        nrow: problem.nrow,
        binary: &problem.binary,
        uplocks: &problem.uplocks,
        downlocks: &problem.downlocks,
    }
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::Box (lines 90-93)
struct Box {
    lo: f64,
    hi: f64,
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::boxed_domain (lines 95-107)
///
/// Paper Sect. 6: an artificial bounding box of [-1e5, +1e5] for infinite
/// bounds. A box with one finite side is anchored at it.
fn boxed_domain(lo: f64, hi: f64) -> Box {
    let k_box: f64 = kBox;
    if lo > -kHighsInf && hi < kHighsInf {
        return Box { lo, hi };
    }
    if lo > -kHighsInf {
        return Box { lo, hi: lo + k_box };
    }
    if hi < kHighsInf {
        return Box { lo: hi - k_box, hi };
    }
    Box {
        lo: -k_box,
        hi: k_box,
    }
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::finite_clamp_helper (lines 109-112)
fn finite_clamp_helper(val: f64, lo: f64, hi: f64) -> f64 {
    let b = boxed_domain(lo, hi);
    cpp_max(b.lo, cpp_min(b.hi, val))
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::acquire_engine (lines 125-154)
///
/// Lazy-construct (or reset) the cached PropEngine inside the scratch. Called
/// exactly once per attempt (only from `fpr_attempt_begin`). The pointer-
/// identity check guards against problem-buffer reuse across attempts on a
/// stale scratch; the C++ compares raw pointers, ported here with
/// `std::ptr::eq` on `as_ptr()`.
fn acquire_engine<'p>(
    prop_engine: &mut Option<PropEngine<'p>>,
    c: &AttemptCtx<'p>,
    csc: &'p CscMatrix,
    deadline: &Deadline,
) {
    let engine_valid = {
        match prop_engine {
            Some(e) => {
                e.ncol() == c.ncol
                    && e.nrow() == c.nrow
                    && std::ptr::eq(e.ar_start().as_ptr(), c.ar_start.as_ptr())
                    && std::ptr::eq(e.ar_index().as_ptr(), c.ar_index.as_ptr())
                    && std::ptr::eq(e.ar_value().as_ptr(), c.ar_value.as_ptr())
                    && std::ptr::eq(e.csc_start().as_ptr(), csc.col_start.as_ptr())
                    && std::ptr::eq(e.csc_row().as_ptr(), csc.col_row.as_ptr())
                    && std::ptr::eq(e.csc_val().as_ptr(), csc.col_val.as_ptr())
                    && std::ptr::eq(e.col_lb().as_ptr(), c.col_lb.as_ptr())
                    && std::ptr::eq(e.col_ub().as_ptr(), c.col_ub.as_ptr())
                    && std::ptr::eq(e.row_lo().as_ptr(), c.row_lo.as_ptr())
                    && std::ptr::eq(e.row_hi().as_ptr(), c.row_hi.as_ptr())
                    && std::ptr::eq(e.integrality().as_ptr(), c.integrality.as_ptr())
                    && e.feastol() == c.feastol
            }
            None => false,
        }
    };
    if !engine_valid {
        *prop_engine = Some(PropEngine::new(
            c.ncol,
            c.nrow,
            c.ar_start,
            c.ar_index,
            c.ar_value,
            &csc.col_start,
            &csc.col_row,
            &csc.col_val,
            c.col_lb,
            c.col_ub,
            c.row_lo,
            c.row_hi,
            c.integrality,
            c.feastol,
        ));
    } else {
        prop_engine.as_mut().unwrap().reset();
    }
    // Arm the engine's own wall-clock poll (issue #151) on both paths: the
    // cached engine survives across attempts, so a `set_deadline` only on the
    // emplace path would leave every attempt after the first propagating
    // against a null deadline.
    prop_engine.as_mut().unwrap().set_deadline(deadline);
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::attempt_engine (lines 160-166)
///
/// The engine `fpr_attempt_begin` engaged for the attempt in flight. Reaching
/// step/finish without one is calling them out of order; the C++ aborts
/// unconditionally (an `assert` is a no-op under NDEBUG).
fn attempt_engine<'a, 'p>(prop_engine: &'a mut Option<PropEngine<'p>>) -> &'a mut PropEngine<'p> {
    if prop_engine.is_none() {
        eprintln!("fpr_core: attempt step/finish without fpr_attempt_begin.");
        std::process::abort();
    }
    prop_engine.as_mut().unwrap()
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::choose_fix_value (lines 174-183)
///
/// Strategy-aware value selection. `cfg.strategy` must be non-null (issue
/// #120 deleted the legacy null-strategy branch).
fn choose_fix_value<'p>(
    j: HighsInt,
    strategy: &FprStrategyConfig,
    lp_ref: Option<&[f64]>,
    c: &AttemptCtx<'p>,
    e: &mut PropEngine<'p>,
    csc: &CscMatrix,
    rng: &mut S5Rng,
) -> f64 {
    let b = boxed_domain(e.var(j).lb, e.var(j).ub);
    choose_value(
        j,
        b.lo,
        b.hi,
        is_integer(c.integrality, j),
        c.col_cost[j as usize],
        strategy.val_strategy,
        rng,
        lp_ref,
        Some(c.row_lo),
        Some(c.row_hi),
        if e.activities_initialized() {
            Some(e.min_activity_data())
        } else {
            None
        },
        if e.activities_initialized() {
            Some(e.max_activity_data())
        } else {
            None
        },
        Some(csc),
    )
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::compute_alt (lines 185-195)
fn compute_alt<'p>(j: HighsInt, preferred: f64, c: &AttemptCtx<'p>, e: &PropEngine<'p>) -> f64 {
    if c.is_binary(j) {
        return if preferred < 0.5 { 1.0 } else { 0.0 };
    }
    let b = boxed_domain(e.var(j).lb, e.var(j).ub);
    let mut alt = if (preferred - b.lo).abs() < c.feastol {
        b.hi
    } else {
        b.lo
    };
    if is_integer(c.integrality, j) {
        alt = alt.round();
    }
    alt
}

/// mip-heuristics (36f7d9d): fpr_core.cpp::is_row_violated_in_ctx (lines 197-199)
fn is_row_violated_in_ctx(i: HighsInt, lhs: f64, c: &AttemptCtx<'_>) -> bool {
    is_row_violated(lhs, c.row_lo[i as usize], c.row_hi[i as usize], c.feastol)
}

/// mip-heuristics (36f7d9d): heuristic_context.h::ProblemView::objective (lines 231-237)
///
/// The view's minimisation-form objective: obj = offset; obj += col_cost[j]*x[j]
/// in j order.
fn problem_objective(problem: &S5Problem, x: &[f64]) -> f64 {
    let mut obj = problem.offset;
    for j in 0..problem.ncol {
        obj += problem.col_cost[j as usize] * x[j as usize];
    }
    obj
}

// ---------------------------------------------------------------------------
// Mode predicates (fpr_strategies.h lines 31-61)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr_strategies.h::is_dynamic_var_strategy (lines 31-33)
fn is_dynamic_var_strategy(s: VarStrategy) -> bool {
    s == VarStrategy::kDomainSize
}

/// mip-heuristics (36f7d9d): fpr_strategies.h::mode_propagates (lines 53-55)
fn mode_propagates(m: FrameworkMode) -> bool {
    m != FrameworkMode::kDive
}

/// mip-heuristics (36f7d9d): fpr_strategies.h::mode_repairs (lines 57-61)
fn mode_repairs(m: FrameworkMode) -> bool {
    // RepairSearch has its own Phase 3 dispatch — not the WalkSAT path.
    m == FrameworkMode::kDfsrep || m == FrameworkMode::kDive || m == FrameworkMode::kDiveprop
}

/// mip-heuristics (36f7d9d): fpr_strategies.h::mode_backtracks (lines 63-66)
fn mode_backtracks(m: FrameworkMode) -> bool {
    m == FrameworkMode::kDfs || m == FrameworkMode::kDfsrep || m == FrameworkMode::kRepairSearch
}

// ---------------------------------------------------------------------------
// HeuristicResult (heuristic_common.h lines 34-54)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): heuristic_common.h::HeuristicResult (lines 34-54)
#[derive(Clone, Debug)]
pub(crate) struct HeuristicResult {
    pub(crate) found_feasible: bool,
    pub(crate) solution: Vec<f64>,
    pub(crate) objective: f64,
    pub(crate) effort: usize,
}

impl Default for HeuristicResult {
    fn default() -> Self {
        Self {
            found_feasible: false,
            solution: Vec::new(),
            objective: f64::INFINITY,
            effort: 0,
        }
    }
}

impl HeuristicResult {
    /// mip-heuristics (36f7d9d): heuristic_common.h::HeuristicResult::infeasible_point (lines 49-54)
    ///
    /// A complete integer assignment that violates at least one row.
    /// `found_feasible` stays false; `objective` is left at its default and is
    /// *not* the cost of `point`.
    pub(crate) fn infeasible_point(point: Vec<f64>, e: usize) -> Self {
        let mut r = HeuristicResult::default();
        r.solution = point;
        r.effort = e;
        r
    }
}

// ---------------------------------------------------------------------------
// fpr_attempt_begin (fpr_core.cpp lines 213-411)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr_core.cpp::fpr_attempt_begin (lines 213-411)
///
/// Phases 1-2 of the paper's Fix-Propagate-Repair (Fig. 4): variable order,
/// value selection, propagation and DFS seeding, as one resumable state
/// machine so an attempt can pause at the budget gate.
pub(crate) fn fpr_attempt_begin<'a, 'p>(
    state: &mut FprAttemptState,
    problem: &'p S5Problem,
    deadline: &Deadline,
    cfg: &mut FprConfig<'a, 'p>,
    rng: &mut S5Rng,
    attempt_idx: i32,
) {
    // Read cfg fields into locals before borrowing scratch (cfg is not passed
    // to helpers; rule 3). The strategy / var-order checks stay after the
    // degenerate return below, as in the C++ (:235-242).
    let mode = cfg.mode;
    let cfg_strategy = cfg.strategy;
    let lp_ref = cfg.lp_ref;
    let cfg_precomputed = cfg.precomputed_var_order;
    let var_order_size = cfg.precomputed_var_order_size;

    let scratch: &mut FprScratch<'p> = cfg
        .scratch
        .as_deref_mut()
        .expect("fpr_attempt_begin requires cfg.scratch");
    let c = make_ctx(problem);

    // Lifecycle reset.
    *state = FprAttemptState::default();
    state.ncol = c.ncol;
    state.nrow = c.nrow;
    state.attempt_idx = attempt_idx;

    if c.ncol == 0 || c.nrow == 0 {
        // Degenerate model — no DFS to do. Leave phase = kIdle so a subsequent
        // finish takes its `ncol == 0 || nrow == 0` guard (the only return that
        // hands back an empty solution, #155).
        state.phase = FprPhase::kIdle;
        return;
    }

    let csc = &problem.csc;
    let strategy = cfg_strategy.expect("FprConfig::strategy must be set (issue #120)");
    // Required, and not only in a debug build (C++ :236-242).
    let Some(precomputed) = cfg_precomputed else {
        eprintln!("fpr_core: fpr_attempt_begin without a precomputed variable order.");
        std::process::abort();
    };

    // Destructure the scratch into disjoint fields (rule 4).
    let FprScratch {
        prop_engine,
        dfs_stack,
        var_order,
        lhs_cache,
        ..
    } = scratch;

    // --- Phase 1: variable ranking -------------------------------------------
    var_order.clear();
    var_order.extend_from_slice(&precomputed[..var_order_size as usize]);
    state.var_order_size = var_order.len() as HighsInt;

    // Ensure scratch.lhs_cache has capacity for finish().
    lhs_cache.resize(c.nrow as usize, 0.0);

    // --- Acquire PropEngine (resets if cached engine is from a previous attempt) ---
    acquire_engine(prop_engine, &c, csc, deadline);

    // No E.sol(j) seeding here (issue #122): nothing between acquire_engine's
    // reset() and the final extraction reads a column's value before writing it.

    // Shuffle top 30% of ranking for diversity on later attempts.
    if attempt_idx > 0 {
        let shuffle_len = cpp_max(1 as HighsInt, c.ncol * 3 / 10);
        s5_shuffle(&mut var_order[..shuffle_len as usize], rng);
    }

    // Row activities are needed by two independent consumers, and either one
    // alone arms them: `val_loosedyn`'s dynamic lock counts, and (since #124)
    // `repair_walk`, whose whole violation measure *is* the activity range.
    let arm_activities = strategy.val_strategy == ValStrategy::kLoosedyn || mode_repairs(mode);
    if arm_activities {
        attempt_engine(prop_engine).init_activities();
    }

    // Trivially-roundable fixings (paper Section 6).
    if !c.uplocks.is_empty() {
        for j in 0..c.ncol {
            if !is_integer(c.integrality, j) || attempt_engine(prop_engine).var(j).fixed {
                continue;
            }
            let b = boxed_domain(
                attempt_engine(prop_engine).var(j).lb,
                attempt_engine(prop_engine).var(j).ub,
            );
            let up = c.uplocks[j as usize];
            let down = c.downlocks[j as usize];
            if up == 0 && down != 0 {
                attempt_engine(prop_engine).fix(j, b.hi);
            } else if down == 0 && up != 0 {
                attempt_engine(prop_engine).fix(j, b.lo);
            }
        }
    }

    // First round of constraint propagation to fixpoint.
    for j in 0..c.ncol {
        if attempt_engine(prop_engine).var(j).fixed {
            attempt_engine(prop_engine).seed_worklist(j);
        }
    }
    // Verdict deliberately discarded (pre-existing; #127). Nothing below reads
    // it; an undetected kInfeasible here is not silently wrong (see C++ :335).
    let _ = attempt_engine(prop_engine).propagate(-1);

    // --- Phase 2 setup -------------------------------------------------------
    state.dynamic_var = is_dynamic_var_strategy(strategy.var_strategy);
    state.do_propagate = mode_propagates(mode);
    state.do_backtrack = mode_backtracks(mode);
    // Fig. 1's `repair` parameter. `mode_repairs` is exactly the paper's three
    // repair-enabled presets (dfsrep / dive / diveprop); it excludes
    // `kRepairSearch`, whose repair procedure (Fig. 5) runs only at the leaf.
    state.do_repair = mode_repairs(mode);
    state.node_limit = c.ncol + 1;
    state.var_order_cursor = 0;
    state.nodes_visited = 0;
    state.found_complete = false;

    dfs_stack.clear();
    let dfs_reserve: usize = if state.do_backtrack {
        2 * (c.ncol as usize)
    } else {
        c.ncol as usize
    };
    if dfs_stack.capacity() < dfs_reserve {
        dfs_stack.reserve(dfs_reserve);
    }

    if state.dynamic_var {
        attempt_engine(prop_engine).init_domain_pq();
    }

    // Seed root DFS node.
    let mut first_var: HighsInt = -1;
    let mut first_idx: HighsInt = -1;
    if state.dynamic_var {
        first_var = attempt_engine(prop_engine).pq_top();
        first_idx = 0;
    } else {
        while state.var_order_cursor < state.var_order_size {
            let j = var_order[state.var_order_cursor as usize];
            if is_integer(c.integrality, j) && !attempt_engine(prop_engine).var(j).fixed {
                first_var = j;
                first_idx = state.var_order_cursor;
                break;
            }
            state.var_order_cursor += 1;
        }
    }

    if first_var < 0 {
        // All integers fixed by propagation; DFS is trivial — go straight to
        // finish.
        state.found_complete = true;
        state.phase = FprPhase::kReadyToFinish;
    } else {
        let pref = choose_fix_value(
            first_var,
            strategy,
            lp_ref,
            &c,
            attempt_engine(prop_engine),
            csc,
            rng,
        );
        let alt = compute_alt(first_var, pref, &c, attempt_engine(prop_engine));
        let vs_m = attempt_engine(prop_engine).vs_mark();
        let sol_m = attempt_engine(prop_engine).sol_mark();
        let act_m = attempt_engine(prop_engine).act_mark();
        let pq_m = if attempt_engine(prop_engine).pq_initialized() {
            attempt_engine(prop_engine).pq_mark()
        } else {
            -1
        };
        let cursor_pt = first_idx + 1;

        if state.do_backtrack {
            dfs_stack.push(FprDfsNode {
                var: first_var,
                val: alt,
                vs_mark: vs_m,
                sol_mark: sol_m,
                act_mark: act_m,
                pq_mark: pq_m,
                cursor_reset: cursor_pt,
            });
        }
        dfs_stack.push(FprDfsNode {
            var: first_var,
            val: pref,
            vs_mark: vs_m,
            sol_mark: sol_m,
            act_mark: act_m,
            pq_mark: pq_m,
            cursor_reset: cursor_pt,
        });
        state.phase = FprPhase::kDfs;
    }

    state.effort_consumed = attempt_engine(prop_engine).effort();
}

// ---------------------------------------------------------------------------
// fpr_attempt_step (fpr_core.cpp lines 422-685)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr_core.cpp::fpr_attempt_step (lines 422-685)
pub(crate) fn fpr_attempt_step<'a, 'p>(
    state: &mut FprAttemptState,
    problem: &'p S5Problem,
    deadline: &Deadline,
    cfg: &mut FprConfig<'a, 'p>,
    rng: &mut S5Rng,
    effort_remaining: usize,
) -> FprStepResult {
    debug_assert!(
        state.phase == FprPhase::kDfs,
        "fpr_attempt_step called outside kDfs phase"
    );
    // Read cfg fields into locals before borrowing scratch (rule 3).
    let strategy = cfg
        .strategy
        .expect("FprConfig::strategy must be set (issue #120)");
    let lp_ref = cfg.lp_ref;
    let walksat_iterations = cfg.walksat_iterations;
    let repair_noise = cfg.repair_noise;

    let scratch: &mut FprScratch<'p> = cfg.scratch.as_deref_mut().expect("cfg.scratch != nullptr");
    let c = make_ctx(problem);
    let csc = &problem.csc;

    // Destructure the scratch into disjoint fields (rule 4). The binding is
    // named `repair_walk_scratch` so it does not shadow the `repair_walk` fn.
    let FprScratch {
        prop_engine,
        dfs_stack,
        var_order,
        repair_walk: repair_walk_scratch,
        ..
    } = scratch;

    // Per-call DFS budget is `effort_remaining` (a slice of the worker's
    // `attempt_budget`). Gate on the *delta* `E.effort() - effort_at_call_start`,
    // not on absolute `E.effort()` (the bug that hangs `infeasible-mip0`).
    let effort_at_call_start = attempt_engine(prop_engine).effort();
    let effort_target_delta = effort_remaining;

    // The wall clock, alongside the effort gate above (issue #117).
    let mut nodes_since_poll: HighsInt = 0;
    // stats.rs throughput counter: DFS fixings (each one is a committed move).
    let ls_nodes_start = state.nodes_visited;

    while !dfs_stack.is_empty()
        && state.nodes_visited < state.node_limit
        && !state.found_complete
        && attempt_engine(prop_engine)
            .effort()
            .wrapping_sub(effort_at_call_start)
            < effort_target_delta
    {
        nodes_since_poll += 1;
        if nodes_since_poll >= kDeadlinePollNodes {
            nodes_since_poll = 0;
            if deadline.expired() {
                // Leave the attempt alive: the stack is non-empty, so the
                // verdict test below reports `kBudgetGate` and the DFS resumes
                // on the next call.
                break;
            }
        }
        let node = dfs_stack.pop().unwrap();
        state.nodes_visited += 1;

        attempt_engine(prop_engine).backtrack_to(
            node.vs_mark,
            node.sol_mark,
            node.act_mark,
            node.pq_mark,
        );
        state.var_order_cursor = node.cursor_reset;

        // Node processing, Fig. 1 lines 4-10, in the paper's own order.
        // `infeas` starts as the verdict of `Apply(fixing, P)` and is then
        // handed to propagation, to repair, and finally to the backtrack
        // decision. Reading it as "prune here" earlier is the defect #124 names.
        let mut infeas = !attempt_engine(prop_engine).fix(node.var, node.val);

        // The rest of `Apply(fixing, P)`: the fixing is also infeasible when it
        // leaves one of the fixed column's rows unsatisfiable by *any*
        // completion of the current domain (issue #124). Gated on `do_repair`.
        if !infeas && state.do_repair {
            let mut apply_effort: usize = 0;
            infeas = any_violated_row_in_column(
                attempt_engine(prop_engine),
                node.var,
                &mut apply_effort,
            );
            attempt_engine(prop_engine).add_effort(apply_effort);
        }

        if !infeas && state.do_propagate {
            // Only a proven inconsistency makes the node infeasible (issue
            // #127). Budget exhaustion is sound but incomplete propagation.
            let pr = attempt_engine(prop_engine).propagate(node.var);
            if pr == PropResult::kDeadlineExpired {
                // The fixpoint's own poll fired (issue #151). Stop the way the
                // pre-node `deadline.expired()` break above stops: that one
                // breaks *before* popping, so push `node` back to leave a
                // paused stack of exactly the same shape.
                dfs_stack.push(node);
                break;
            }
            infeas = pr == PropResult::kInfeasible;
        }

        // Fig. 1 lines 7-8: `if infeas and repair: infeas = RepairWalk(P)`. The
        // repair runs on the *partial* assignment this node's domain encodes
        // (issue #124).
        if infeas && state.do_repair {
            let mut walk_effort: usize = 0;
            let repaired = repair_walk(
                attempt_engine(prop_engine),
                walksat_iterations,
                repair_noise,
                rng,
                &mut walk_effort,
                repair_walk_scratch,
                deadline,
            );
            // Charge into the engine's own counter, which is what both the
            // loop's budget gate above and `state.effort_consumed` below read.
            attempt_engine(prop_engine).add_effort(walk_effort);
            infeas = !repaired;
        }

        // Fig. 1 lines 9-10: `if infeas and backtrackOnInfeas: Backtrack`.
        // Popping the next node *is* the backtrack here.
        if infeas && state.do_backtrack {
            continue;
        }

        // Fig. 1 line 11: `branches = Branch(P)`. Reached with `infeas` still
        // true in a non-backtracking mode: `dive`/`diveprop` construct "a
        // complete solution in a single big dive" and are not entitled to stop
        // early.
        let (next_var, next_idx) = {
            if state.dynamic_var {
                (attempt_engine(prop_engine).pq_top(), 0 as HighsInt)
            } else {
                let mut found = (-1 as HighsInt, -1 as HighsInt);
                while state.var_order_cursor < state.var_order_size {
                    let j = var_order[state.var_order_cursor as usize];
                    if is_integer(c.integrality, j) && !attempt_engine(prop_engine).var(j).fixed {
                        found = (j, state.var_order_cursor);
                        break;
                    }
                    state.var_order_cursor += 1;
                }
                found
            }
        };

        if next_var < 0 {
            // Fig. 1 lines 12-16: no branches left. Line 14's backtrack has
            // already run for every mode that has one.
            debug_assert!(
                !(infeas && state.do_backtrack),
                "Fig. 1 line 14 already backtracked this node"
            );
            // Non-backtracking mode. `found_complete` is set even on a
            // still-infeasible leaf; Phase 2.5/3 are a leaf post-process whose
            // precondition has always been `found_complete`, not Fig. 1's
            // running `infeas` (see C++ :615-656).
            state.found_complete = true;
            break;
        }

        let pref = choose_fix_value(
            next_var,
            strategy,
            lp_ref,
            &c,
            attempt_engine(prop_engine),
            csc,
            rng,
        );
        let alt = compute_alt(next_var, pref, &c, attempt_engine(prop_engine));
        let vs_m = attempt_engine(prop_engine).vs_mark();
        let sol_m = attempt_engine(prop_engine).sol_mark();
        let act_m = attempt_engine(prop_engine).act_mark();
        let pq_m = if attempt_engine(prop_engine).pq_initialized() {
            attempt_engine(prop_engine).pq_mark()
        } else {
            -1
        };
        let cursor_pt = next_idx + 1;

        if state.do_backtrack {
            dfs_stack.push(FprDfsNode {
                var: next_var,
                val: alt,
                vs_mark: vs_m,
                sol_mark: sol_m,
                act_mark: act_m,
                pq_mark: pq_m,
                cursor_reset: cursor_pt,
            });
        }
        dfs_stack.push(FprDfsNode {
            var: next_var,
            val: pref,
            vs_mark: vs_m,
            sol_mark: sol_m,
            act_mark: act_m,
            pq_mark: pq_m,
            cursor_reset: cursor_pt,
        });
    }

    state.effort_consumed = attempt_engine(prop_engine).effort();
    super::stats::FPR.add_moves(0, (state.nodes_visited - ls_nodes_start) as u64);

    // Verdict determined? Found a leaf or stack/node-limit exhausted.
    if state.found_complete || dfs_stack.is_empty() || state.nodes_visited >= state.node_limit {
        state.phase = FprPhase::kReadyToFinish;
        return FprStepResult::kVerdictReady;
    }

    // Budget gate hit; attempt is alive.
    FprStepResult::kBudgetGate
}

// ---------------------------------------------------------------------------
// fpr_attempt_finish (fpr_core.cpp lines 696-908)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): fpr_core.cpp::fpr_attempt_finish (lines 696-908)
pub(crate) fn fpr_attempt_finish<'a, 'p>(
    state: &mut FprAttemptState,
    problem: &'p S5Problem,
    deadline: &Deadline,
    cfg: &mut FprConfig<'a, 'p>,
    rng: &mut S5Rng,
) -> HeuristicResult {
    // Read cfg fields into locals before borrowing scratch (rule 3).
    let mode = cfg.mode;
    let cfg_strategy = cfg.strategy;
    let lp_ref = cfg.lp_ref;
    let cont_fallback = cfg.cont_fallback;
    let repair_iterations = cfg.repair_iterations;
    let repair_noise = cfg.repair_noise;
    let repair_track_best = cfg.repair_track_best;
    let walksat_iterations = cfg.walksat_iterations;

    let scratch: &mut FprScratch<'p> = cfg.scratch.as_deref_mut().expect("cfg.scratch != nullptr");
    let c = make_ctx(problem);

    // Degenerate model from begin() — short-circuit cleanly.
    if c.ncol == 0 || c.nrow == 0 {
        state.phase = FprPhase::kIdle;
        return HeuristicResult::default();
    }

    let csc = &problem.csc;

    // repair_search needs the whole scratch while the primary engine E must
    // also be passed: take the engine out (rule 6) and put it back before
    // every return (the C++ engine persists in scratch). A missing engine is
    // `attempt_engine`'s abort (C++ :160-166).
    let mut e: PropEngine<'p> = {
        attempt_engine(&mut scratch.prop_engine);
        scratch.prop_engine.take().unwrap()
    };
    // C++ asserts this inside choose_fix_value (:176); read after the
    // degenerate return so a degenerate model never reaches it.
    let strategy = cfg_strategy.expect("FprConfig::strategy must be set (issue #120)");

    // Phase 2.5: fix remaining unfixed variables (continuous + residual
    // integers). This runs whether or not the DFS reached a leaf (issue #155):
    // Sect. 2.3 of Mexi et al. is explicit that fix-and-propagate produces an
    // integer point even when it ran into an infeasibility.
    for j in 0..c.ncol {
        if e.var(j).fixed {
            continue;
        }
        let lo = e.var(j).lb;
        let hi = e.var(j).ub;

        if !is_integer(c.integrality, j) {
            if c.col_cost[j as usize].abs() > 1e-15 {
                let want_low = c.col_cost[j as usize] > 0.0;
                let v = finite_clamp_helper(if want_low { lo } else { hi }, lo, hi);
                *e.sol_mut(j) = v;
            } else {
                let fallback = if let Some(cf) = cont_fallback {
                    cf[j as usize]
                } else {
                    0.0
                };
                let v = finite_clamp_helper(fallback, lo, hi);
                *e.sol_mut(j) = v;
            }
        } else {
            let v = choose_fix_value(j, strategy, lp_ref, &c, &mut e, csc, rng);
            let v = cpp_max(lo, cpp_min(hi, v)).round();
            *e.sol_mut(j) = v;
        }
        let v = cpp_max(
            c.col_lb[j as usize],
            cpp_min(c.col_ub[j as usize], e.sol(j)),
        );
        *e.sol_mut(j) = v;
    }

    // `auto& solution = scratch.solution` (C++ :761): take it out so the
    // repair_search call (which also borrows scratch whole) does not alias it.
    let mut solution: Vec<f64> = std::mem::take(&mut scratch.solution);
    solution.clear();
    solution.extend_from_slice(&e.sol_data()[..c.ncol as usize]);
    let mut total_prop_work = e.effort();

    // `auto& lhs_cache = scratch.lhs_cache` (C++ :765): likewise taken out.
    let mut lhs_cache: Vec<f64> = std::mem::take(&mut scratch.lhs_cache);
    lhs_cache.resize(c.nrow as usize, 0.0);
    total_prop_work += c.ar_index.len();
    for i in 0..c.nrow {
        let mut lhs = 0.0;
        let mut k = c.ar_start[i as usize];
        while k < c.ar_start[(i + 1) as usize] {
            lhs += c.ar_value[k as usize] * solution[c.ar_index[k as usize] as usize];
            k += 1;
        }
        lhs_cache[i as usize] = lhs;
    }

    // The DFS never reached a leaf: no Phase 3, no verdict, but the point the
    // fill just completed goes back to the caller (issue #155). It is charged
    // the row rebuild rather than the untouched `E.effort()`.
    if !state.found_complete {
        state.phase = FprPhase::kIdle;
        state.effort_consumed = total_prop_work;
        scratch.prop_engine = Some(e);
        // C++ moves scratch.solution into the result here; it stays empty.
        scratch.lhs_cache = lhs_cache;
        return HeuristicResult::infeasible_point(solution, total_prop_work);
    }

    let mut feasible = true;
    for i in 0..c.nrow {
        if is_row_violated_in_ctx(i, lhs_cache[i as usize], &c) {
            feasible = false;
            break;
        }
    }

    // Phase 3: RepairSearch (Fig. 5) or WalkSAT. The deadline gates entry to
    // both (issue #117). Neither takes an effort cap from here any more (#156).
    if !feasible && !deadline.expired() && mode == FrameworkMode::kRepairSearch {
        let mut rs_effort: usize = 0;
        feasible = repair_search(
            &mut e,
            &mut solution,
            &mut lhs_cache,
            c.col_lb,
            c.col_ub,
            c.row_lo,
            c.row_hi,
            repair_iterations,
            repair_noise,
            repair_track_best,
            rng,
            &mut rs_effort,
            scratch,
            deadline,
            None,
            K_REPAIR_PROGRESS_THRESHOLD,
        );
        total_prop_work += rs_effort;
    } else if !feasible && !deadline.expired() && mode_repairs(mode) {
        let mut walk_effort: usize = 0;
        feasible = walksat_repair(
            &e,
            &mut solution,
            &mut lhs_cache,
            c.col_lb,
            c.col_ub,
            walksat_iterations,
            repair_noise,
            repair_track_best,
            rng,
            &mut walk_effort,
            &mut scratch.walksat,
        );
        total_prop_work += walk_effort;
    }

    // Phase 3's answer is not the verdict; this re-check is (see
    // `fpr_attempt_step`, which names this the one verdict site). Folded into
    // one `feasible` rather than into two returns (see C++ :869-880).
    if feasible {
        for i in 0..c.nrow {
            if is_row_violated_in_ctx(i, lhs_cache[i as usize], &c) {
                feasible = false;
                break;
            }
        }
    }

    if !feasible {
        state.phase = FprPhase::kIdle;
        state.effort_consumed = total_prop_work;
        scratch.prop_engine = Some(e);
        // C++ moves scratch.solution into the result here; it stays empty.
        scratch.lhs_cache = lhs_cache;
        return HeuristicResult::infeasible_point(solution, total_prop_work);
    }

    greedy_1opt(
        &e,
        &mut solution,
        &mut lhs_cache,
        c.col_cost,
        &mut total_prop_work,
    );

    let obj = problem_objective(problem, &solution);

    let mut result = HeuristicResult::default();
    result.found_feasible = true;
    result.solution = solution;
    result.objective = obj;
    result.effort = total_prop_work;
    state.phase = FprPhase::kIdle;
    state.effort_consumed = total_prop_work;
    scratch.prop_engine = Some(e);
    scratch.lhs_cache = lhs_cache;
    result
}

// fpr_attempt (fpr_core.cpp:919-955, the one-shot wrapper) is NOT ported: its
// callers (scylla, fpr_lp) are out of scope for S5.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;
    use crate::mip::fpr_local_mip::rng::S5Rng;
    use crate::mip::fpr_local_mip::{CscMatrix, FprScratch, S5Problem};

    fn make_problem_2x2() -> S5Problem {
        // CSR: row 0 has entries at col 0 (coeff 1) and col 1 (coeff 1)
        //       row 1 has entries at col 0 (coeff 1) and col 1 (coeff 1)
        // CSC:  col 0 has entries at row 0 (coeff 1) and row 1 (coeff 1)
        //       col 1 has entries at row 0 (coeff 1) and row 1 (coeff 1)
        S5Problem {
            col_lower: vec![0.0, 0.0],
            col_upper: vec![10.0, 10.0],
            row_lower: vec![-kHighsInf, 1.0],
            row_upper: vec![3.0, kHighsInf],
            col_cost: vec![1.0, 1.0],
            offset: 0.0,
            integrality: vec![HighsVarType::kInteger, HighsVarType::kInteger],
            ar_start: vec![0, 2, 4],
            ar_index: vec![0, 1, 0, 1],
            ar_value: vec![1.0, 1.0, 1.0, 1.0],
            csc: CscMatrix {
                col_start: vec![0, 2, 4],
                col_row: vec![0, 1, 0, 1],
                col_val: vec![1.0, 1.0, 1.0, 1.0],
            },
            uplocks: vec![0, 0],
            downlocks: vec![0, 0],
            feastol: 1e-6,
            epsilon: 1e-6,
            ncol: 2,
            nrow: 2,
            nnz: 4,
            incumbent: vec![0.0, 0.0],
            binary: vec![0, 0],
        }
    }

    fn run_lifecycle(
        problem: &S5Problem,
        mode: FrameworkMode,
        strategy: &FprStrategyConfig,
        order: &[HighsInt],
    ) -> HeuristicResult {
        let mut scratch = FprScratch::default();
        let mut rng = S5Rng::new(42);
        let deadline = Deadline {
            timer: None,
            limit: f64::INFINITY,
        };
        let mut cfg = FprConfig::default();
        cfg.mode = mode;
        cfg.strategy = Some(strategy);
        cfg.precomputed_var_order = Some(order);
        cfg.precomputed_var_order_size = order.len() as HighsInt;
        cfg.scratch = Some(&mut scratch);

        let mut state = FprAttemptState::default();
        fpr_attempt_begin(&mut state, problem, &deadline, &mut cfg, &mut rng, 0);
        // `step` only in kDfs (C++ asserts it, fpr_core.cpp:425).
        while state.phase == FprPhase::kDfs {
            let r = fpr_attempt_step(
                &mut state,
                problem,
                &deadline,
                &mut cfg,
                &mut rng,
                usize::MAX,
            );
            if r == FprStepResult::kVerdictReady {
                break;
            }
        }
        fpr_attempt_finish(&mut state, problem, &deadline, &mut cfg, &mut rng)
    }

    fn check_solution(problem: &S5Problem, result: &HeuristicResult) {
        assert_eq!(result.solution.len(), problem.ncol as usize);
        for j in 0..problem.ncol {
            let v = result.solution[j as usize];
            assert!(v.is_finite(), "col {j} value {v} is not finite");
            assert!(
                v >= problem.col_lower[j as usize] - problem.feastol
                    && v <= problem.col_upper[j as usize] + problem.feastol,
                "col {j} value {v} outside bounds"
            );
            if is_integer(&problem.integrality, j) {
                assert!(
                    (v - v.round()).abs() < 1e-6,
                    "col {j} value {v} not integral"
                );
            }
        }
    }

    #[test]
    fn s5_fprcore_dfs_mode() {
        let p = make_problem_2x2();
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        let r = run_lifecycle(&p, FrameworkMode::kDfs, &strat, &[0, 1]);
        check_solution(&p, &r);
        assert!(r.found_feasible);
    }

    #[test]
    fn s5_fprcore_dive_mode() {
        let p = make_problem_2x2();
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        let r = run_lifecycle(&p, FrameworkMode::kDive, &strat, &[0, 1]);
        check_solution(&p, &r);
        assert!(r.found_feasible);
    }

    #[test]
    fn s5_fprcore_dfsrep_mode() {
        let p = make_problem_2x2();
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        let r = run_lifecycle(&p, FrameworkMode::kDfsrep, &strat, &[0, 1]);
        check_solution(&p, &r);
        assert!(r.found_feasible);
    }

    #[test]
    fn s5_fprcore_diveprop_mode() {
        let p = make_problem_2x2();
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        let r = run_lifecycle(&p, FrameworkMode::kDiveprop, &strat, &[0, 1]);
        check_solution(&p, &r);
        assert!(r.found_feasible);
    }

    #[test]
    fn s5_fprcore_repair_search_mode() {
        let p = make_problem_2x2();
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        let r = run_lifecycle(&p, FrameworkMode::kRepairSearch, &strat, &[0, 1]);
        check_solution(&p, &r);
        assert!(r.found_feasible);
    }

    #[test]
    fn s5_fprcore_finish_walksat_greedy() {
        // Feasible 2x2 model (row 0: x0 + x1 <= 3, row 1: x0 + x1 >= 1). The
        // DFS leaf is feasible here, so this covers finish's feasible path and
        // greedy_1opt; the WalkSAT branch (finish's Phase 3 on an infeasible
        // leaf) is not reached by this model and has no test here.
        let p = make_problem_2x2();
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        let r = run_lifecycle(&p, FrameworkMode::kDfsrep, &strat, &[0, 1]);
        check_solution(&p, &r);
        assert!(r.found_feasible);
    }

    /// Row x0 + x1 >= 3 with x in [0,1]: no assignment satisfies it, so a
    /// non-backtracking dive reaches a leaf (`found_complete`) whose point is
    /// violated, and finish enters Phase 3.
    fn make_problem_infeasible() -> S5Problem {
        let mut p = make_problem_2x2();
        p.col_upper = vec![1.0, 1.0];
        p.row_lower = vec![3.0, -kHighsInf];
        p.row_upper = vec![kHighsInf, kHighsInf];
        p
    }

    #[test]
    fn s5_fprcore_finish_walksat_infeasible_leaf() {
        // kDive: mode_repairs, so finish's Phase 3 is walksat_repair, which
        // cannot repair this row; the result is an infeasible_point (#155).
        let p = make_problem_infeasible();
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        let r = run_lifecycle(&p, FrameworkMode::kDive, &strat, &[0, 1]);
        check_solution(&p, &r);
        assert!(!r.found_feasible);
        assert_eq!(r.objective, f64::INFINITY);
        assert!(r.effort >= p.ar_index.len());
    }

    #[test]
    fn s5_fprcore_finish_repair_search_infeasible_leaf() {
        // Phase 1's trivially-roundable fixings bypass `Apply` (C++ :650-656):
        // downlocks != 0 and uplocks == 0 fix both columns at their upper
        // bound 1 in begin, every integer is fixed, so begin goes straight to
        // kReadyToFinish with found_complete and a violated leaf (1 + 1 < 3).
        // kRepairSearch then enters repair_search, the repairing modes
        // walksat_repair, and kDfs no Phase 3; none can satisfy the row.
        let mut p = make_problem_infeasible();
        p.downlocks = vec![1, 1];
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        for mode in [
            FrameworkMode::kRepairSearch,
            FrameworkMode::kDiveprop,
            FrameworkMode::kDfsrep,
            FrameworkMode::kDfs,
        ] {
            let r = run_lifecycle(&p, mode, &strat, &[0, 1]);
            check_solution(&p, &r);
            assert!(!r.found_feasible, "{mode:?}");
            assert_eq!(r.objective, f64::INFINITY, "{mode:?}");
        }
        // kDfs on the same model without the locks backtracks and exhausts the
        // tree (found_complete false): the Phase 2.5 point still comes back.
        p.downlocks = vec![0, 0];
        let r = run_lifecycle(&p, FrameworkMode::kDfs, &strat, &[0, 1]);
        check_solution(&p, &r);
        assert!(!r.found_feasible);
    }

    #[test]
    fn s5_fprcore_attempt_idx_shuffle() {
        // 10 columns, one row sum(x) >= 1: shuffle_len = max(1, 10*3/10) = 3,
        // so attempt_idx > 0 shuffles exactly var_order[..3] with the
        // attempt's rng (no draw precedes it in begin) and keeps the rest.
        let n: HighsInt = 10;
        let p = S5Problem {
            col_lower: vec![0.0; 10],
            col_upper: vec![1.0; 10],
            row_lower: vec![1.0],
            row_upper: vec![kHighsInf],
            col_cost: vec![1.0; 10],
            offset: 0.0,
            integrality: vec![HighsVarType::kInteger; 10],
            ar_start: vec![0, n],
            ar_index: (0..n).collect(),
            ar_value: vec![1.0; 10],
            csc: CscMatrix {
                col_start: (0..=n).collect(),
                col_row: vec![0; 10],
                col_val: vec![1.0; 10],
            },
            uplocks: vec![0; 10],
            downlocks: vec![0; 10],
            feastol: 1e-6,
            epsilon: 1e-6,
            ncol: n,
            nrow: 1,
            nnz: 10,
            incumbent: vec![0.0; 10],
            binary: vec![1; 10],
        };
        let strat = FprStrategyConfig {
            var_strategy: VarStrategy::kType,
            val_strategy: ValStrategy::kGoodobj,
        };
        let order: Vec<HighsInt> = (0..n).rev().collect();
        let deadline = Deadline {
            timer: None,
            limit: f64::INFINITY,
        };
        for attempt_idx in [0, 1, 2] {
            let mut scratch = FprScratch::default();
            let mut rng = S5Rng::new(7);
            let mut cfg = FprConfig::default();
            cfg.mode = FrameworkMode::kDfs;
            cfg.strategy = Some(&strat);
            cfg.precomputed_var_order = Some(&order);
            cfg.precomputed_var_order_size = n;
            cfg.scratch = Some(&mut scratch);
            let mut state = FprAttemptState::default();
            fpr_attempt_begin(&mut state, &p, &deadline, &mut cfg, &mut rng, attempt_idx);
            assert_eq!(state.var_order_size, n);

            let mut expected = order.clone();
            if attempt_idx > 0 {
                let mut ref_rng = S5Rng::new(7);
                s5_shuffle(&mut expected[..3], &mut ref_rng);
            }
            assert_eq!(scratch.var_order, expected, "attempt_idx {attempt_idx}");
        }
    }

    #[test]
    fn s5_fprcore_boxed_domain() {
        let b = boxed_domain(-kHighsInf, kHighsInf);
        assert_eq!(b.lo, -1e5);
        assert_eq!(b.hi, 1e5);

        let b = boxed_domain(3.0, kHighsInf);
        assert_eq!(b.lo, 3.0);
        assert_eq!(b.hi, 3.0 + 1e5);

        let b = boxed_domain(-kHighsInf, 4.0);
        assert_eq!(b.lo, 4.0 - 1e5);
        assert_eq!(b.hi, 4.0);

        let b = boxed_domain(-2.0, 5.0);
        assert_eq!(b.lo, -2.0);
        assert_eq!(b.hi, 5.0);
    }
}
