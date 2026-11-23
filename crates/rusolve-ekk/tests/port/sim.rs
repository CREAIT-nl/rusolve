//! Simplex-trace drivers shared by `t_ekk_init`, `t_dual_*` and
//! `t_primal_trace`.
//!
//! The port has no trace hooks (signature note S1), so the per-record
//! stream of a fixture cannot be produced. What can be observed:
//!
//! 1. **HEkk initialisation**: the state after `initialise_for_solve()`
//!    reproduces the `ekk_solve_enter`, first `invert` and
//!    `ekk_solve_initialised` records exactly (nothing between the
//!    `kOracleEkkComputeFactor` hook and the end of `initialiseForSolve`
//!    touches the factor or the fields of the `invert` record).
//! 2. **Solve exit**: after `HEkk::solve`, the `ekk_solve_exit` record
//!    (`HEkk::returnFromEkkSolve` changes none of its fields after the hook).
//! 3. **State after iteration k** (iteration-limit method): with
//!    `simplex_iteration_limit = k`, HEkkDual (HEkkDual.cpp lines 651, 916)
//!    and HEkkPrimal (HEkkPrimal.cpp lines 458, 554) test `bailout()`
//!    immediately after `iterate()`, so the solve stops right after the
//!    iteration that the fixture's `iter` record with `iteration == k`
//!    describes. `HEkk::returnFromSolve` (kIterationLimit branch) then
//!    recomputes the work arrays (`initialiseBound`, `initialiseCost`,
//!    `computePrimal`, `computeDual`, `zeroBasicDuals`), so `h_workDual`,
//!    `h_baseValue`, `h_workValue` are NOT comparable; it leaves
//!    `basis_.basicIndex_`, `dual_edge_weight_`, `info_.update_count`,
//!    `info_.updated_*_objective_value`, `total_synthetic_tick_` and the
//!    `info_.*_density` fields alone, and those ARE compared with record k.
//!    `HEkkPrimal::edge_weight_` lives in the local primal solver, so the
//!    primal `h_edge_weight` / `edge_weight_in` are not observable.
//! 4. **Second solve**: when the scaled LP ends Infeasible / Unbounded /
//!    with unscaled infeasibilities, HApp::solveLpSimplex moves the LP
//!    back, unscales it and calls `HEkk::solve(true)` again on the same
//!    HEkk (solve_call 1). [`run_happ`] replays that refinement block, and
//!    the `ekk_solve_enter`, `ekk_solve_initialised` and `ekk_solve_exit`
//!    records of every solve call are compared.
#![allow(dead_code, unused_imports)]

use super::*;
use crate::common::events::*;
use crate::common::*;
use rusolve_ekk::prelude::*;
use serde_json::Value;

/// `MatrixRec` view of a port matrix (format 3 carries `p_end`).
pub fn matrix_rec_of(m: &HighsSparseMatrix) -> MatrixRec {
    MatrixRec {
        format: m.format_ as i32,
        num_col: m.num_col_,
        num_row: m.num_row_,
        start: m.start_.clone(),
        p_end: (m.format_ == MatrixFormat::kRowwisePartitioned).then(|| m.p_end_.clone()),
        index: m.index_.clone(),
        value: m.value_.clone(),
    }
}

/// Every `ekk_solve_enter` of the trace, in order. The drivers replay
/// HApp::solveLpSimplex, which calls `HEkk::solve` once, or twice when it
/// refines the unscaled LP ([`run_happ`]).
pub fn enters_of(fx: &Fixture) -> Vec<EkkSolveEnter> {
    let v: Vec<_> = fx
        .of_type("ekk_solve_enter")
        .map(|r| EkkSolveEnter::from_rec(&r))
        .collect();
    assert!(
        (1..=2).contains(&v.len()),
        "{}: the drivers replay one or two HEkk::solve calls, trace has {}",
        fx.file,
        v.len()
    );
    for (i, e) in v.iter().enumerate() {
        assert_eq!(e.solve_call, i as i32, "{}: solve_call order", fx.file);
    }
    v
}

/// The trace's first `ekk_solve_enter` (the test input).
pub fn enter_of(fx: &Fixture) -> EkkSolveEnter {
    enters_of(fx).swap_remove(0)
}

/// `ekk_solve_enter` as the oracle would emit it for `e` (before
/// `initialiseForSolve`). `options` is the input echoed back (the port's
/// HighsOptions has no by-name view).
pub fn enter_obj(e: &HEkk, input: &EkkSolveEnter) -> Value {
    enter_obj_call(e, &input.options, 0)
}

/// [`enter_obj`] for solve call `solve_call`, echoing `options` except for
/// the fields HApp changes between solves (`simplex_strategy`), which are
/// read from the port so that a wrong choice shows as a mismatch.
pub fn enter_obj_call(e: &HEkk, options: &Value, solve_call: i32) -> Value {
    let mut options = options.clone();
    options["simplex_strategy"] = Value::from(e.options_.simplex_strategy);
    EkkSolveEnter {
        solve_call,
        iteration_count: e.iteration_count_,
        status: status_of(&e.status_),
        options,
        cost_scale: e.cost_scale_,
        lp: lp_rec_of(&e.lp_),
        basis: e.status_.has_basis.then(|| basis_of(&e.basis_)),
    }
    .to_obj()
    .build()
}

/// The first `invert` record, from the state after `initialise_for_solve`.
pub fn invert_obj(e: &HEkk) -> Value {
    let f = &e.simplex_nla_.factor_;
    InvertEvent {
        solve_call: 0,
        iteration_count: e.iteration_count_,
        rank_deficiency: f.rank_deficiency,
        factor_build_synthetic_tick: f.build_synthetic_tick,
        build_synthetic_tick: e.build_synthetic_tick_,
        total_synthetic_tick: e.total_synthetic_tick_,
        basis_matrix_num_el: f.basis_matrix_num_el,
        invert_num_el: f.invert_num_el,
        kernel_dim: f.kernel_dim,
        kernel_num_el: f.kernel_num_el,
        has_invert: e.status_.has_invert,
        h_basic_index: fnv_i32(&e.basis_.basic_index_),
    }
    .to_obj()
    .build()
}

pub fn initialised_obj(e: &HEkk) -> Value {
    initialised_obj_call(e, 0)
}

pub fn initialised_obj_call(e: &HEkk, solve_call: i32) -> Value {
    let i = &e.info_;
    EkkSolveInitialised {
        solve_call,
        status: status_of(&e.status_),
        basis: basis_of(&e.basis_),
        scalars: scalars_of(e),
        num_tot_random_value: i.num_tot_random_value_.clone(),
        num_tot_permutation: i.num_tot_permutation_.clone(),
        num_col_permutation: i.num_col_permutation_.clone(),
        simplex_strategy: i.simplex_strategy,
        dual_edge_weight_strategy: i.dual_edge_weight_strategy,
        primal_edge_weight_strategy: i.primal_edge_weight_strategy,
        price_strategy: i.price_strategy,
        dual_simplex_cost_perturbation_multiplier: i.dual_simplex_cost_perturbation_multiplier,
        primal_simplex_bound_perturbation_multiplier: i
            .primal_simplex_bound_perturbation_multiplier,
        factor_pivot_threshold: i.factor_pivot_threshold,
        update_limit: i.update_limit,
        arrays: arrays_of(e),
        ar_matrix: matrix_rec_of(&e.ar_matrix_),
    }
    .to_obj()
    .build()
}

pub fn exit_obj(e: &HEkk, solve_call: i32) -> Value {
    EkkSolveExit {
        solve_call,
        model_status: e.model_status_ as i32,
        exit_algorithm: e.exit_algorithm_ as i32,
        solve_bailout: e.solve_bailout_,
        scalars: scalars_of(e),
        status: status_of(&e.status_),
        basis: basis_of(&e.basis_),
        arrays: arrays_of(e),
    }
    .to_obj()
    .build()
}

/// HApp prefix + `HEkk::solve(false)` with `simplex_iteration_limit = limit`
/// (first solve only).
pub fn run_solve(enter: &EkkSolveEnter, ctx: &str, limit: HighsInt) -> HEkk {
    let (mut ekk, _parts) = ekk_after_move_lp(enter, ctx, Some(limit));
    at(format!("HEkk::solve (simplex_iteration_limit = {limit})"));
    ekk.solve(false);
    ekk
}

/// The outcome of [`run_happ`].
pub struct HAppRun {
    pub ekk: HEkk,
    /// Number of `HEkk::solve` calls made.
    pub solves: usize,
    /// `ekk_solve_enter`, `ekk_solve_initialised` (only with `emit`) and
    /// `ekk_solve_exit` per solve call, as the oracle would emit them.
    pub records: Vec<Value>,
}

/// Emit the `ekk_solve_enter` and `ekk_solve_initialised` records of the
/// coming `HEkk::solve`: the statements before each hook are run on a
/// clone (`ekk_solve_prefix`, then `initialiseForSolve`), which is what
/// `HEkk::solve` does first (HEkk.cpp:1003-1018).
fn emit_solve_entry(ekk: &HEkk, options: &Value, solve_call: i32, out: &mut Vec<Value>) {
    let mut ekk = ekk.clone();
    ekk_solve_prefix(&mut ekk);
    out.push(enter_obj_call(&ekk, options, solve_call));
    ekk.initialise_for_solve();
    out.push(initialised_obj_call(&ekk, solve_call));
}

/// HApp::solveLpSimplex (HApp.h:94-488) on the path the traces take: a
/// fresh `Highs` with no basis, no dualize/permute, and
/// `simplex_unscaled_solution_strategy = kSimplexUnscaledSolutionStrategyRefine`.
/// An LP without scaling is solved once.
/// `HEkk::solve` runs once, then the refinement block decides whether to
/// solve the unscaled LP again (solve_call 1).
///
/// `limit0` is the `simplex_iteration_limit` of the first solve; `limit1`,
/// when given, overrides it for the second (the iteration-limit method on
/// a second-solve iteration). `enters` are the trace's `ekk_solve_enter`
/// records: `[0]` is the input, `[1]` (if any) only supplies the options
/// echo of the emitted second `ekk_solve_enter`.
///
/// Not replayed: work that only fills `Highs` members HEkk never reads
/// again (`unscaleSolution`, the HighsInfo copies after the last solve,
/// restoring the options), and the third solve HApp makes when the
/// unscaled LP first reaches kObjectiveBound (no objective bound is set in
/// the traces; asserted).
pub fn run_happ(
    enters: &[EkkSolveEnter],
    ctx: &str,
    limit0: HighsInt,
    limit1: Option<HighsInt>,
    emit: bool,
) -> HAppRun {
    use HighsModelStatus as M;
    let mut records = Vec::new();
    let (mut ekk, mut parts) = ekk_after_move_lp(&enters[0], ctx, Some(limit0));
    assert_eq!(
        parts.options.simplex_unscaled_solution_strategy, kSimplexUnscaledSolutionStrategyRefine,
        "{ctx}: test replays the refine path"
    );
    if emit {
        emit_solve_entry(&ekk, &enters[0].options, 0, &mut records);
    }
    at(format!(
        "HEkk::solve, solve_call 0 (simplex_iteration_limit = {limit0})"
    ));
    let return_status = ekk.solve(false);
    records.push(exit_obj(&ekk, 0));
    let done = |ekk: HEkk, records: Vec<Value>, solves: usize| HAppRun {
        ekk,
        solves,
        records,
    };
    // `incumbent_lp.scale_.has_scaling`: C++ reads it from the moved-from
    // incumbent, whose bool members are copied by the move, so it equals
    // the HEkk LP's (the port's mem::take reset the incumbent's).
    if !ekk.lp_.scale_.has_scaling {
        // Unscaled LP: one solve, no refinement (HApp.h:235-254).
        ekk.unpermute();
        ekk.undualize();
        assert!(!ekk.status_.is_permuted && !ekk.status_.is_dualized);
        assert_eq!(parts.options.cost_scale_factor, 0, "{ctx}");
        let _ = return_status;
        return done(ekk, records, 1);
    }
    at("HApp between solves (getSolution, moveBackLpAndUnapplyScaling, getUnscaledInfeasibilities, moveLp)".into());
    ekk.unpermute();
    ekk.undualize();
    assert!(!ekk.status_.is_permuted && !ekk.status_.is_dualized);
    assert_eq!(parts.options.cost_scale_factor, 0, "{ctx}");
    if return_status == HighsStatus::kError {
        return done(ekk, records, 1);
    }
    let scaled_model_status = ekk.model_status_;
    parts.info.objective_function_value = ekk.info_.primal_objective_value;
    parts.info.simplex_iteration_count = ekk.iteration_count_;
    parts.solution = ekk.get_solution();
    // getHighsBasis(ekk_lp): it only reads the LP, so pass a copy rather
    // than borrow `ekk.lp_` mutably next to `&ekk`.
    let mut ekk_lp = ekk.lp_.clone();
    parts.basis = ekk.get_highs_basis(&mut ekk_lp);
    assert!(parts.basis.valid);
    parts.info.basis_validity = kBasisValidityValid;
    parts.lp.move_back_lp_and_unapply_scaling(&mut ekk.lp_);
    ekk.set_nla_pointers_for_lp_and_scale(&parts.lp);
    // unscaleSolution(solution, incumbent_lp.scale_): Highs-only, skipped.
    get_unscaled_infeasibilities(
        &parts.options,
        &parts.lp.scale_,
        &ekk.basis_,
        &ekk.info_,
        &mut parts.info,
    );
    let num_unscaled_primal_infeasibilities = parts.info.num_primal_infeasibilities;
    let num_unscaled_dual_infeasibilities = parts.info.num_dual_infeasibilities;
    let scaled_optimality_but_unscaled_infeasibilities = scaled_model_status == M::kOptimal
        && (num_unscaled_primal_infeasibilities != 0 || num_unscaled_dual_infeasibilities != 0);
    let scaled_objective_target_but_unscaled_primal_infeasibilities =
        scaled_model_status == M::kObjectiveTarget && parts.info.num_primal_infeasibilities > 0;
    let scaled_objective_bound_but_unscaled_dual_infeasibilities =
        scaled_model_status == M::kObjectiveBound && parts.info.num_dual_infeasibilities > 0;
    let refine_solution = parts.options.simplex_unscaled_solution_strategy
        == kSimplexUnscaledSolutionStrategyRefine
        && (scaled_optimality_but_unscaled_infeasibilities
            || scaled_model_status == M::kInfeasible
            || scaled_model_status == M::kUnboundedOrInfeasible
            || scaled_model_status == M::kUnbounded
            || scaled_objective_bound_but_unscaled_dual_infeasibilities
            || scaled_objective_target_but_unscaled_primal_infeasibilities
            || scaled_model_status == M::kUnknown);
    if !refine_solution {
        return done(ekk, records, 1);
    }
    {
        let mut so = HighsLpSolverObject::new(
            &mut parts.lp,
            &mut parts.basis,
            &mut parts.solution,
            &mut parts.info,
            &mut parts.callback,
            &mut parts.options,
            &mut parts.timer,
        );
        ekk.move_lp(&mut so);
    }
    let mut solve_unscaled_lp = true;
    if scaled_model_status == M::kInfeasible && ekk.exit_algorithm_ == SimplexAlgorithm::kDual {
        assert!(ekk.dual_ray_record_.index != kNoRayIndex, "{ctx}");
    }
    if scaled_model_status == M::kInfeasible && ekk.dual_ray_record_.index != kNoRayIndex {
        at("HEkk::proofOfPrimalInfeasibility (unscaled LP, scaled ray)".into());
        ekk.set_nla_pointers_for_lp_and_scale(&parts.lp);
        if ekk.proof_of_primal_infeasibility() {
            solve_unscaled_lp = false;
        }
    }
    if !solve_unscaled_lp {
        return done(ekk, records, 1);
    }
    // HApp changes `options` after moveLp; in C++ HEkk sees it through its
    // pointer, the port's HEkk holds a copy (set_pointers), so set both.
    let strategy =
        if num_unscaled_primal_infeasibilities == 0 || scaled_model_status == M::kObjectiveTarget {
            kSimplexStrategyPrimal
        } else {
            assert!(num_unscaled_primal_infeasibilities > 0);
            if (ekk.status_.has_basis || parts.basis.valid)
                && !ekk.status_.has_dual_steepest_edge_weights
            {
                ekk.info_.dual_edge_weight_strategy = kSimplexEdgeWeightStrategyDevex;
            }
            kSimplexStrategyDual
        };
    parts.options.simplex_strategy = strategy;
    ekk.options_.simplex_strategy = strategy;
    if let Some(k) = limit1 {
        parts.options.simplex_iteration_limit = k;
        ekk.options_.simplex_iteration_limit = k;
    }
    let force_phase2 = parts.options.simplex_unscaled_solution_strategy
        != kSimplexUnscaledSolutionStrategyDirect
        && scaled_model_status != M::kObjectiveTarget;
    if emit {
        let options = enters.get(1).map_or(&enters[0].options, |e| &e.options);
        emit_solve_entry(&ekk, options, 1, &mut records);
    }
    at(format!(
        "HEkk::solve, solve_call 1 (force_phase2 = {force_phase2}, simplex_iteration_limit = {})",
        ekk.options_.simplex_iteration_limit
    ));
    ekk.solve(force_phase2);
    records.push(exit_obj(&ekk, 1));
    assert!(
        !(scaled_model_status != M::kObjectiveBound
            && ekk.model_status_ == M::kObjectiveBound
            && ekk.info_.num_dual_infeasibilities > 0),
        "{ctx}: third HEkk::solve (objective-bound refinement) is not replayed"
    );
    done(ekk, records, 2)
}

// ---------------------------------------------------------------------------
// Expected per-iteration state
// ---------------------------------------------------------------------------

/// One `iter` record of a trace (dual or primal).
#[derive(Clone, Debug)]
pub struct IterExp {
    /// 1-based fixture line.
    pub line: usize,
    /// The `HEkk::solve` call the iteration belongs to (0 or 1).
    pub solve_call: i32,
    pub iteration: i32,
    pub dual: bool,
    pub solve_phase: i32,
    pub row_out: i32,
    pub variable_in: i32,
    pub variable_out: i32,
    pub d: Option<DualIterDone>,
    pub p: Option<(PrimalIterP, PrimalIterDone)>,
    /// One-line summaries of the non-`iter` records since the previous
    /// `iter` (rebuilds, inverts, phase ends): what happened in between.
    pub between: Vec<String>,
}

pub fn iters_of(fx: &Fixture) -> Vec<IterExp> {
    let mut out = Vec::new();
    let mut between = Vec::new();
    let mut solve_call = -1;
    for (i, v) in fx.records.iter().enumerate() {
        let r = fx.rec(i);
        let t = r.t();
        if t == "ekk_solve_enter" {
            solve_call = r.int("solve_call");
        }
        if t != "iter" {
            if t != "iter_none" && t != "trace_header" && t != "result" {
                let mut s = format!("line {} {t}", i + 1);
                if r.has("reason") {
                    s.push_str(&format!(" reason={}", r.int("reason")));
                }
                if r.has("solve_phase") {
                    s.push_str(&format!(" solve_phase={}", r.int("solve_phase")));
                }
                if r.has("algorithm") {
                    s.push_str(&format!(" ({})", r.str("algorithm")));
                }
                between.push(s);
            }
            continue;
        }
        let _ = v;
        let dual = r.str("algorithm") == "dual";
        let e = if dual {
            let it = DualIter::from_rec(&r);
            let d = it.done.unwrap();
            IterExp {
                line: i + 1,
                solve_call,
                iteration: it.iteration,
                dual,
                solve_phase: it.solve_phase,
                row_out: it.row_out,
                variable_in: d.variable_in,
                variable_out: d.variable_out,
                d: Some(d),
                p: None,
                between: std::mem::take(&mut between),
            }
        } else {
            let it = PrimalIter::from_rec(&r);
            IterExp {
                line: i + 1,
                solve_call,
                iteration: it.iteration,
                dual,
                solve_phase: it.p.solve_phase,
                row_out: it.p.row_out,
                variable_in: it.p.variable_in,
                variable_out: it.p.variable_out.unwrap_or_else(|| {
                    panic!("{}:{} primal iter without variable_out", fx.file, i + 1)
                }),
                p: Some((it.p.clone(), it.done.unwrap())),
                d: None,
                between: std::mem::take(&mut between),
            }
        };
        out.push(e);
    }
    out
}

/// Basic set after each iteration, replayed from the initial logical basis
/// (`ekk_solve_initialised.basis.basicIndex`) by the records' in/out pairs
/// (bound flips: `row_out < 0`, no basis change).
pub fn replay_sets(fx: &Fixture, iters: &[IterExp]) -> Vec<Vec<i32>> {
    let init = fx.first("ekk_solve_initialised").basis("basis").basic_index;
    let mut cur: Vec<i32> = init;
    cur.sort();
    let mut out = vec![cur.clone()];
    for it in iters {
        if it.row_out >= 0 && it.variable_out >= 0 {
            let pos = cur.binary_search(&it.variable_out).unwrap_or_else(|_| {
                panic!("{}:{} variable_out not basic in replay", fx.file, it.line)
            });
            cur.remove(pos);
            let pos = cur.binary_search(&it.variable_in).unwrap_err();
            cur.insert(pos, it.variable_in);
        }
        out.push(cur.clone());
    }
    out
}

fn sorted(v: &[i32]) -> Vec<i32> {
    let mut v = v.to_vec();
    v.sort();
    v
}

fn minus(a: &[i32], b: &[i32]) -> Vec<i32> {
    a.iter()
        .filter(|x| b.binary_search(x).is_err())
        .copied()
        .collect()
}

/// Compare the port's state after a solve stopped at iteration k with
/// fixture record `iters[k-1]`; `sets[k]` / `sets[k-1]` from
/// [`replay_sets`]. Returns `None` when everything observable matches.
pub fn check_iteration(
    fx: &Fixture,
    iters: &[IterExp],
    sets: &[Vec<i32>],
    k: usize,
    e: &HEkk,
) -> Option<Checker> {
    let it = &iters[k - 1];
    let alg = if it.dual { "dual" } else { "primal" };
    let mut chk = Checker::new(format!(
        "{}:{} {alg} iter iteration={} (solve_phase={}, expected row_out={}, variable_in={}, variable_out={}), port stopped by simplex_iteration_limit={}",
        fx.file,
        it.line,
        it.iteration,
        it.solve_phase,
        it.row_out,
        it.variable_in,
        it.variable_out,
        it.iteration
    ));
    chk.int(
        "HEkk::iteration_count_ (control flow: the port must reach this iteration)",
        it.iteration,
        e.iteration_count_,
    );
    chk.int(
        "HEkk::model_status_ (14 = kIterationLimit: the port must still be iterating)",
        HighsModelStatus::kIterationLimit as i32,
        e.model_status_ as i32,
    );
    let port_set = sorted(&e.basis_.basic_index_);
    let exp_set = &sets[k];
    if &port_set != exp_set {
        let prev = &sets[k - 1];
        let entered = minus(&port_set, prev);
        let left = minus(prev, &port_set);
        let mut note = String::new();
        if prev == &port_set {
            note.push_str("port made no basis change (bound flip or no iteration)");
        } else {
            note.push_str(&format!(
                "port entered {entered:?} and removed {left:?} (relative to iteration {} which matched)",
                it.iteration - 1
            ));
        }
        let culprit = if entered != [it.variable_in] && left == [it.variable_out] {
            if it.dual {
                "CHUZC: HEkkDualRow::chooseMakepack/choosePossible/chooseFinal (HEkkDual::chooseColumn)"
            } else {
                "CHUZC: HEkkPrimal::chuzc / chooseColumn / hyperChooseColumn"
            }
        } else if left != [it.variable_out] && entered == [it.variable_in] {
            if it.dual {
                "CHUZR: HEkkDualRHS::chooseNormal (HEkkDual::chooseRow)"
            } else {
                "CHUZR: HEkkPrimal::chooseRow / phase1ChooseRow"
            }
        } else {
            "control flow / both pivots (HEkkDual::iterate, HEkkPrimal::iterate, rebuild)"
        };
        chk.push(
            "basic set after the iteration",
            format!("+{} -{}", it.variable_in, it.variable_out),
            note,
            format!("likely: {culprit}"),
        );
    }
    chk.hash(
        "h_basicIndex (order; HEkk::updatePivots, HFactor::build permutation)",
        match (&it.d, &it.p) {
            (Some(d), _) => d.h_basic_index,
            (_, Some((_, p))) => p.h_basic_index,
            _ => unreachable!(),
        },
        fnv_i32(&e.basis_.basic_index_),
    );
    let i = &e.info_;
    if let Some(d) = &it.d {
        chk.int("info_.update_count", d.update_count, i.update_count);
        chk.f64(
            "info_.updated_dual_objective_value (HEkkDual::updateDual, HEkkDualRow::updateFlip)",
            d.updated_dual_objective_value,
            i.updated_dual_objective_value,
        );
        chk.f64(
            "HEkk::total_synthetic_tick_ (tick contract: HFactor/HVector synthetic_tick)",
            d.total_synthetic_tick,
            e.total_synthetic_tick_,
        );
        chk.f64(
            "info_.row_ep_density (HEkk::updateOperationResultDensity)",
            d.row_ep_density,
            i.row_ep_density,
        );
        chk.f64("info_.row_ap_density", d.row_ap_density, i.row_ap_density);
        chk.f64("info_.col_aq_density", d.col_aq_density, i.col_aq_density);
        chk.f64(
            "info_.row_DSE_density",
            d.row_dse_density,
            i.row_dse_density,
        );
        if let Some(w) = d.edge_weight_row_out_after {
            let a = e
                .dual_edge_weight_
                .get(it.row_out as usize)
                .copied()
                .unwrap_or(f64::NAN);
            chk.f64(
                "dual_edge_weight_[row_out] (HEkkDual::updateFtranDSE, HEkk::updateDualSteepestEdgeWeights, devex)",
                w,
                a,
            );
        }
        chk.hash(
            "h_dual_edge_weight (DSE/Devex weight update)",
            d.h_dual_edge_weight,
            fnv_f64(&e.dual_edge_weight_),
        );
    }
    if let Some((_, p)) = &it.p {
        chk.int("info_.update_count", p.update_count, i.update_count);
        chk.f64(
            "info_.updated_primal_objective_value (HEkkPrimal::update / phase2UpdatePrimal)",
            p.updated_primal_objective_value,
            i.updated_primal_objective_value,
        );
        chk.f64(
            "HEkk::total_synthetic_tick_ (tick contract)",
            p.total_synthetic_tick,
            e.total_synthetic_tick_,
        );
    }
    if chk.ok() {
        return None;
    }
    if !it.between.is_empty() {
        chk.push(
            "context",
            String::new(),
            String::new(),
            format!(
                "records since iteration {}: {}",
                it.iteration - 1,
                it.between.join("; ")
            ),
        );
    }
    Some(chk)
}

/// Run to iteration k and check it; a panic inside the port becomes a
/// failing checker.
pub fn probe(
    fx: &Fixture,
    enters: &[EkkSolveEnter],
    iters: &[IterExp],
    sets: &[Vec<i32>],
    k: usize,
) -> Option<Checker> {
    let ctx = format!("{} k={k}", fx.file);
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let it = &iters[k - 1];
        let run = if it.solve_call == 0 {
            run_happ(enters, &ctx, it.iteration, None, false)
        } else {
            // Let the first solve finish as in the trace, stop the second.
            run_happ(
                enters,
                &ctx,
                full_run_limit(iters),
                Some(it.iteration),
                false,
            )
        };
        check_iteration(fx, iters, sets, k, &run.ekk)
    }));
    match r {
        Ok(res) => res,
        Err(p) => {
            let msg = p
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "<non-string panic>".into());
            let it = &iters[k - 1];
            let mut chk = Checker::new(format!(
                "{}:{} iteration={} (simplex_iteration_limit={})",
                fx.file, it.line, it.iteration, it.iteration
            ));
            chk.push("PANIC inside the port", String::new(), msg, String::new());
            Some(chk)
        }
    }
}

/// First iteration whose observable state differs. `linear` scans
/// k = 1, 2, ...; otherwise bisects (assumes a divergence persists, which
/// holds for the basis and tick fields once they differ).
pub fn first_diverging_iteration(fx: &Fixture, linear: bool) -> Option<(usize, Checker)> {
    let enter = enters_of(fx);
    let iters = iters_of(fx);
    let sets = replay_sets(fx, &iters);
    if iters.is_empty() {
        return None;
    }
    if linear {
        for k in 1..=iters.len() {
            if let Some(found) = probe(fx, &enter, &iters, &sets, k) {
                return Some((k, found));
            }
        }
        return None;
    }
    let n = iters.len();
    let last = probe(fx, &enter, &iters, &sets, n)?;
    let (mut lo, mut hi, mut bad) = (0usize, n, last);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        match probe(fx, &enter, &iters, &sets, mid) {
            Some(found) => {
                hi = mid;
                bad = found;
            }
            None => lo = mid,
        }
    }
    Some((hi, bad))
}

/// `simplex_iteration_limit` for a full run: above the fixture's count, so
/// a port that needs more iterations stops with kIterationLimit instead of
/// looping.
pub fn full_run_limit(iters: &[IterExp]) -> HighsInt {
    2 * iters.last().map_or(0, |i| i.iteration) + 100
}

/// Full HApp run (one or two `HEkk::solve` calls), then compare the
/// `ekk_solve_enter`, `ekk_solve_initialised` and `ekk_solve_exit` records
/// of every solve call.
pub fn exit_divergence(fx: &Fixture) -> Result<(), String> {
    let enters = enters_of(fx);
    let limit = full_run_limit(&iters_of(fx));
    let ctx = fx.file.clone();
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        run_happ(&enters, &ctx, limit, None, true).records
    }));
    let actual = match r {
        Ok(v) => v,
        Err(p) => {
            let msg = p
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "<non-string panic>".into());
            return Err(format!(
                "{}: PANIC inside the port during HEkk::solve: {msg}",
                fx.file
            ));
        }
    };
    compare_traces(
        fx,
        &actual,
        &TraceCompare::default().only(&[
            "ekk_solve_enter",
            "ekk_solve_initialised",
            "ekk_solve_exit",
        ]),
    )
    .map_err(|d| d.to_string())
}

/// The gate of a full-trace target: exit record, and on a mismatch the
/// first diverging iteration.
pub fn trace_gate(fx: &Fixture, linear: bool) -> Option<String> {
    let exit = exit_divergence(fx);
    let Err(exit_msg) = exit else { return None };
    let first = first_diverging_iteration(fx, linear);
    let mut s = String::new();
    match first {
        Some((_, chk)) => {
            s.push_str(&format!(
                "FIRST DIVERGING ITERATION ({}):\n  {}: {} mismatch(es)",
                if linear { "linear scan" } else { "bisection" },
                chk.ctx,
                chk.failures.len()
            ));
            for m in &chk.failures {
                s.push_str(&format!("\n    {}", m.to_string().replace('\n', "\n    ")));
            }
            s.push_str("\n\nexit record:\n");
        }
        None => s.push_str(
            "every iteration matches; the divergence is after the last iteration (final rebuild / cleanup / returnFromSolve):\n",
        ),
    }
    s.push_str(&exit_msg);
    Some(s)
}

/// The gate of a per-iteration target (`t_dual_rhs`, `t_dual_row`): every
/// iteration's observable state, scanned k = 1, 2, ... (`linear`) or
/// bisected. `component` names what the target owns; the message says how
/// many iterations were verified before the first divergence and which
/// HiGHS function the mismatching fields point at.
pub fn iteration_gate(fx: &Fixture, linear: bool, component: &str) -> Option<String> {
    let (k, chk) = first_diverging_iteration(fx, linear)?;
    let mut s = format!(
        "{component}: verified for iterations 1..{} only; FIRST DIVERGING ITERATION ({}):\n  {}: {} mismatch(es)",
        k - 1,
        if linear { "linear scan" } else { "bisection" },
        chk.ctx,
        chk.failures.len()
    );
    for m in &chk.failures {
        s.push_str(&format!("\n    {}", m.to_string().replace('\n', "\n    ")));
    }
    Some(s)
}

/// Run `gate` on every file, collecting failures.
pub fn run_trace_files(files: Vec<String>, what: &str, gate: impl Fn(&Fixture) -> Option<String>) {
    let mut fails = Failures::default();
    for file in files {
        let Some(fx) = load(&file) else { continue };
        match gate(&fx) {
            None => fails.cases += 1,
            Some(m) => fails.fail(m),
        }
    }
    fails.finish(what);
}
