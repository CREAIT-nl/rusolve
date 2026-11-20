//! LS integer assignment + LP repair of the continuous part (`mip_ls_int_lp_repair`).
//!
//! Idea: Rigal, "Résolution de problèmes à variables mixtes dans LocalSolver", EasyChair
//! Preprint 2601, 2020; Gardi, Benoist, Darlay, Estellon, Megel, *Mathematical
//! Programming Solver Based on Local Search*, Wiley 2014, ch. 3. Local search moves the
//! integer variables; the continuous part is then solved as an LP with the integers fixed.
//!
//! The fix / propagate / LP sequence and its constants follow HiGHS `tryRoundedPoint`
//! (HighsPrimalHeuristics.cpp:971-1059, ported in `primal_heuristics/part2.rs`;
//! MIT License, Copyright (c) 2026 HiGHS, see THIRD_PARTY_NOTICES.md). What is new here
//! is where it is applied (the LS end points of FJ, FPR and Local-MIP, which HiGHS submits
//! with `trySolution` only), the usefulness/duplicate/row filters and the budget.

use crate::prelude::*;

/// Minimum share of LP columns (not in `integer_cols`), as `1 / K_MIN_LP_COL_SHARE_INV`.
/// TODO(perf): first guess (10 %).
const K_MIN_LP_COL_SHARE_INV: usize = 10;
/// Calls per solve (each is one LP at most).
const K_MAX_CALLS: u32 = 8;
/// Total LP iterations per solve, in units of the per-call limit.
const K_TOTAL_ITERS_PER_CALL_LIMIT: i64 = 4;
/// Assignments remembered by the duplicate filter.
const K_MAX_SEEN: usize = 64;

/// Outcome of one `repair` call (counted in `LsIntLpRepair`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RepairOutcome {
    /// Model gate failed (flag off, sub-MIP, too few LP columns): nothing done.
    Gated,
    /// Per-solve call or iteration budget exhausted.
    Budget,
    /// Feasible point and no LP column carries cost: the LP cannot improve it.
    NoCost,
    /// Integer assignment already repaired in this solve.
    Duplicate,
    /// Point has the wrong length or an integer column far from integral.
    BadPoint,
    /// A row cannot be satisfied over the LP columns' bounds (single-row Farkas test).
    RowInfeasible,
    /// Fixing and propagating the integers proved infeasibility.
    PropInfeasible,
    /// The LP with fixed integers is infeasible.
    LpInfeasible,
    /// The LP stopped without a primal-feasible point (iteration limit, numerics).
    LpFailed,
    /// LP point found but the incumbent did not improve (not better, or rejected by
    /// `try_solution`'s explicit check).
    NotImproving,
    /// New incumbent.
    Improved,
}

/// Per-solve state and counters; `mipdata_.ls_int_lp_repair`, `None` with the flag off.
#[derive(Clone, Debug, Default)]
pub(crate) struct LsIntLpRepair {
    pub(crate) calls: u32,
    pub(crate) lp_solves: u32,
    pub(crate) lp_iterations: i64,
    /// O(nnz) work of the row prefilter, in matrix entries visited.
    pub(crate) prefilter_effort: u64,
    pub(crate) improved: u32,
    pub(crate) outcomes: Vec<(i32, RepairOutcome)>,
    seen: Vec<u64>,
}

/// Whether the model qualifies: flag on, not a sub-MIP, and at least 10 % LP columns.
pub(crate) fn gate(ms: &HighsMipSolver) -> bool {
    if !ms.options_mip_.mip_ls_int_lp_repair || ms.submip {
        return false;
    }
    let ncol = ms.num_col().max(0) as usize;
    let nint = ms.mipdata_.integer_cols.len();
    ncol > nint && (ncol - nint) * K_MIN_LP_COL_SHARE_INV >= ncol
}

/// FNV-1a over (column, rounded value) of the integer columns.
fn assignment_hash(integer_cols: &[HighsInt], point: &[f64]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &col in integer_cols {
        let v = point[col as usize].round() as i64;
        for b in (col as u64)
            .to_le_bytes()
            .iter()
            .chain(v.to_le_bytes().iter())
        {
            h ^= *b as u64;
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    }
    h
}

/// Single-row Farkas test with the integers fixed at `fixed` (indexed by column; only
/// `is_int` columns read) and the LP columns free in their global bounds. Returns false
/// when some row cannot be satisfied. `effort` counts matrix entries visited.
fn rows_satisfiable(ms: &HighsMipSolver, is_int: &[bool], fixed: &[f64], effort: &mut u64) -> bool {
    let feastol = ms.mipdata_.feastol;
    let gd = ms.mipdata_.workers[0].get_global_domain();
    let dom = &ms.mipdata_.domains[gd];
    for i in 0..ms.num_row() as usize {
        let start = ms.mipdata_.a_rstart_[i] as usize;
        let end = ms.mipdata_.a_rstart_[i + 1] as usize;
        let mut minact = 0.0f64;
        let mut maxact = 0.0f64;
        for k in start..end {
            let col = ms.mipdata_.a_rindex_[k] as usize;
            let a = ms.mipdata_.a_rvalue_[k];
            if is_int[col] {
                minact += a * fixed[col];
                maxact += a * fixed[col];
            } else {
                let (lo, up) = (dom.col_lower_[col], dom.col_upper_[col]);
                if a > 0.0 {
                    minact += a * lo;
                    maxact += a * up;
                } else {
                    minact += a * up;
                    maxact += a * lo;
                }
            }
        }
        *effort += (end - start) as u64;
        if minact > ms.row_upper_2(i as HighsInt) + feastol
            || maxact < ms.row_lower_2(i as HighsInt) - feastol
        {
            return false;
        }
    }
    true
}

/// Fix the integer columns of the LS point `point`, propagate, solve the LP over the
/// remaining columns and submit the result through `try_solution` with tag `source`.
/// `feasible` says whether `point` itself was already a feasible solution (polish) or an
/// infeasible LS end state (repair). Returns the outcome; only `Improved` changed the
/// incumbent.
pub(crate) fn repair(
    ms: &mut HighsMipSolver,
    point: &[f64],
    source: i32,
    feasible: bool,
) -> RepairOutcome {
    let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_LS_REPAIR);
    if !gate(ms) {
        return RepairOutcome::Gated;
    }
    let mut st = ms.mipdata_.ls_int_lp_repair.take().unwrap_or_default();
    let outcome = repair_body(ms, &mut st, point, feasible, source);
    st.outcomes.push((source, outcome));
    if outcome == RepairOutcome::Improved {
        st.improved += 1;
    }
    crate::highs_log_user!(
        &ms.options_mip_.log_options,
        HighsLogType::kInfo,
        "LS int-LP repair: source {} {} -> {:?} (calls {}, LP solves {}, LP iters {})\n",
        source,
        if feasible { "polish" } else { "repair" },
        outcome,
        st.calls,
        st.lp_solves,
        st.lp_iterations
    );
    ms.mipdata_.ls_int_lp_repair = Some(st);
    outcome
}

fn per_call_iteration_limit(ms: &HighsMipSolver) -> i64 {
    // tryRoundedPoint's limit (HighsPrimalHeuristics.cpp:1015-1016).
    cpp_max(10000i64, 2 * ms.mipdata_.firstrootlpiters)
}

fn repair_body(
    ms: &mut HighsMipSolver,
    st: &mut LsIntLpRepair,
    point: &[f64],
    feasible: bool,
    source: i32,
) -> RepairOutcome {
    let ncol = ms.num_col() as usize;
    if point.len() != ncol {
        return RepairOutcome::BadPoint;
    }
    let per_call = per_call_iteration_limit(ms);
    let remaining = K_TOTAL_ITERS_PER_CALL_LIMIT * per_call - st.lp_iterations;
    if st.calls >= K_MAX_CALLS || remaining <= 0 {
        return RepairOutcome::Budget;
    }
    let integer_cols = ms.mipdata_.integer_cols.clone();
    let mut is_int = vec![false; ncol];
    for &col in &integer_cols {
        is_int[col as usize] = true;
    }
    if feasible && !(0..ncol).any(|j| !is_int[j] && ms.col_cost_2(j as HighsInt) != 0.0) {
        return RepairOutcome::NoCost;
    }
    let h = assignment_hash(&integer_cols, point);
    if st.seen.contains(&h) {
        return RepairOutcome::Duplicate;
    }
    st.calls += 1;
    if st.seen.len() < K_MAX_SEEN {
        st.seen.push(h);
    }

    // Rounded integer values clamped to the global domain (tryRoundedPoint lines 978-992,
    // with the domain bounds of the global domain; a far-from-integral value means the
    // point is not an LS integer assignment).
    let feastol = ms.mipdata_.feastol;
    let mut fixed = vec![0.0f64; ncol];
    {
        let gd = ms.mipdata_.workers[0].get_global_domain();
        let dom = &ms.mipdata_.domains[gd];
        for &col in &integer_cols {
            let cu = col as usize;
            let mut rounded = 0.0;
            if fractionality(point[cu], Some(&mut rounded)) > feastol {
                return RepairOutcome::BadPoint;
            }
            fixed[cu] = cpp_max(dom.col_lower_[cu], cpp_min(dom.col_upper_[cu], rounded));
        }
    }
    if !feasible && !rows_satisfiable(ms, &is_int, &fixed, &mut st.prefilter_effort) {
        return RepairOutcome::RowInfeasible;
    }

    // tryRoundedPoint lines 974-1006: local copy of the global domain, fix, propagate.
    let gd = ms.mipdata_.workers[0].get_global_domain();
    let mut localdom = detach!(ms.mipdata_.domains[gd], |g| HighsDomain::copy_from(ms, g));
    let outcome = fix_and_solve(
        ms,
        st,
        &mut localdom,
        &integer_cols,
        &fixed,
        per_call.min(remaining),
        source,
    );
    localdom.release(ms);
    outcome
}

fn fix_and_solve(
    ms: &mut HighsMipSolver,
    st: &mut LsIntLpRepair,
    localdom: &mut HighsDomain,
    integer_cols: &[HighsInt],
    fixed: &[f64],
    iteration_limit: i64,
    source: i32,
) -> RepairOutcome {
    for &col in integer_cols {
        let cu = col as usize;
        let val = cpp_max(
            localdom.col_lower_[cu],
            cpp_min(localdom.col_upper_[cu], fixed[cu]),
        );
        localdom.fix_col(ms, col, val, HighsDomain_Reason::branching());
        if localdom.infeasible() {
            return RepairOutcome::PropInfeasible;
        }
        let _ = localdom.propagate(ms);
        if localdom.infeasible() {
            return RepairOutcome::PropInfeasible;
        }
    }

    // tryRoundedPoint lines 1008-1033: LP relaxation with the local bounds.
    let numintcols = integer_cols.len() as HighsInt;
    let mut lprelax = HighsLpRelaxation::new(ms);
    lprelax.set_mip_worker(0);
    lprelax.set_profiling();
    lprelax.load_model(ms);
    lprelax.set_iteration_limit(iteration_limit.min(HighsInt::MAX as i64) as HighsInt);
    let _ = lprelax.get_lp_solver().change_cols_bounds_interval(
        0,
        ms.num_col() - 1,
        &localdom.col_lower_,
        &localdom.col_upper_,
    );
    if ms.options_mip_.mip_root_presolve_only {
        let _ = lprelax
            .get_lp_solver()
            .set_option_value_str("presolve", kHighsOffString);
    }
    if !ms.options_mip_.mip_root_presolve_only && (5 * numintcols) / ms.num_col() >= 1 {
        let _ = lprelax
            .get_lp_solver()
            .set_option_value_str("presolve", kHighsOnString);
    } else if ms.mipdata_.firstrootbasis.valid {
        let _ = lprelax
            .get_lp_solver()
            .set_basis(&ms.mipdata_.firstrootbasis, "ls_int_lp_repair");
    }

    let status = lprelax.resolve_lp(ms, None);
    let iters = lprelax.get_num_lp_iterations();
    st.lp_solves += 1;
    st.lp_iterations += iters;
    ms.mipdata_.heuristic_lp_iterations += iters;
    ms.mipdata_.total_lp_iterations += iters;

    if status == HighsLpRelaxation_Status::kInfeasible {
        return RepairOutcome::LpInfeasible;
    }
    if !HighsLpRelaxation::unscaled_primal_feasible(status) {
        return RepairOutcome::LpFailed;
    }
    let mut lpsol: Vec<f64> = lprelax.get_lp_solver().get_solution().col_value.clone();
    for &col in integer_cols {
        lpsol[col as usize] = localdom.col_lower_[col as usize];
    }
    let ub_before = ms.mipdata_.upper_bound;
    let _ = ms.try_solution(&lpsol, source);
    if ms.mipdata_.upper_bound < ub_before {
        RepairOutcome::Improved
    } else {
        RepairOutcome::NotImproving
    }
}
