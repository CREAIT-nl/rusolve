//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP source: conflict_resolution.c (tightenCoefs 504-724, linearCombRows 975-1008,
//!   getQuotLargestSmallestCoef 1647-1674, weakenVarConflictRow 1511-1539,
//!   weakenConflictRow 1543-1596, weakenContinuousVarsConflictRow 1600-1642,
//!   computeScaleReason 2165-2184).

use crate::prelude::*;
use crate::util::highs_cdouble::HighsCDouble;
use crate::util::highs_pdqsort::pdqsort_by;
use crate::util::highs_utils::{cpp_max, cpp_min};

use super::{
    S2_EPSILON, S2_FEASTOL, S2ConflictRow, S2Context, S2DebugSol, remove_zero_vars,
    s2_debugsol_viol_ge, s2_is_eq, s2_is_ge, s2_is_gt, s2_is_infinity, s2_is_le, s2_is_lt,
    s2_is_negative, s2_is_positive, s2_is_zero, s2_lb_at, s2_trace_enabled, s2_ub_at,
};

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::tightenCoefs (lines 504-724)
///
/// Returns nchgcoefs.
pub fn tighten_coefs<C: S2Context>(
    ctx: &C,
    vals: &mut [f64],
    inds: &mut [HighsInt],
    nnz: &mut usize,
    lhs: &mut f64,
) -> HighsInt {
    let mut nchgcoefs: HighsInt = 0;
    let mut maxabsval: f64 = 0.0;
    let mut minacttmp = HighsCDouble::from(0.0);

    // compute activity of the row and get the absolute values of the coefficients
    for i in 0..*nnz {
        let idx = inds[i];
        let coefidx = idx as usize;

        debug_assert!(idx >= 0);
        debug_assert!(idx < ctx.num_col() as HighsInt);

        if vals[coefidx] > 0.0 {
            let lb = ctx.global_lb(idx);

            if s2_is_infinity(-lb) {
                return 0;
            }

            if ctx.is_integral(idx) {
                maxabsval = cpp_max(maxabsval, vals[coefidx]);
            }

            minacttmp += lb * vals[coefidx];
        } else {
            let ub = ctx.global_ub(idx);

            if s2_is_infinity(ub) {
                return 0;
            }

            if ctx.is_integral(idx) {
                maxabsval = cpp_max(maxabsval, -vals[coefidx]);
            }

            minacttmp += ub * vals[coefidx];
        }
    }

    let mut minact: f64 = f64::from(minacttmp);

    // row is redundant if minact is infinity
    if s2_is_infinity(minact) {
        return 0;
    }
    // no coefficients can be tightened
    if s2_is_infinity(-minact) {
        return 0;
    }

    // SCIP conflict_resolution.c:604 `assert(!SCIPsetIsGE(set, minact, *rowlhs))` —
    // "propagating constraint cannot be redundant". SCIP guarantees that invariant by
    // construction (its conflict row starts as the row of the constraint that detected
    // infeasibility), so under NDEBUG the assert is never reached; the port has no such
    // guarantee, because reaches this routine with HiGHS rows: the LP Farkas proof row of
    // `conflict_analysis_2` (the non-SCIP `S2_RESOLVE_LP_PROOFS` hook) and the HiGHS reason
    // rows. Those are valid for the MIP but can be redundant w.r.t. the variable bounds alone,
    // i.e. minact >= lhs. Write the row in bound-slack form: with b_j = lb for coef > 0 and ub
    // otherwise, t_j = coef_j (x_j - b_j) >= 0 and r = lhs - minact, so the row is
    // sum_j t_j >= r. Tightening variable i to newcoef replaces t_i by |newcoef| s_i with
    // s_i = |x_i - b_i| integral, and moves lhs by delta * b_i. For s_i >= 1 the guard already
    // makes the original row satisfied, so it carries no information there and validity needs
    // |newcoef| >= r; both branches set |newcoef| = |minact - lhs|, which equals r only for
    // r >= 0. With minact > lhs, newcoef flips sign, the code takes its "coefficient became
    // zero" exit with delta > 0, and pushes lhs += delta * bound: an arbitrary strengthening.
    // Measured on pg5_34: minact 6.5607e3 against lhs 1.1768e2 moved lhs to 6.5714e3, cutting
    // off the true optimum -14339.35345. Enforce the invariant instead of trusting it, with the
    // assert's own tolerance: `SCIPsetIsGE(minact, lhs)` is `minact - lhs >= -S2_EPSILON`, so
    // `s2_is_ge` is the negation of what the assert allows.
    //
    // The guard has to be non-strict, i.e. also reject `minact == lhs`, because the two
    // "coefficient became zero" exits below are reachable only when |newcoef| = |minact - lhs|
    // is zero or sign-wrong: under the assert they are dead code (the negative branch keeps the
    // coefficient exactly when newcoef = minact - lhs < -S2_EPSILON, the positive branch exactly
    // when newcoef = lhs - minact > S2_EPSILON). A drop moves lhs by -coef_i * bound_i but leaves
    // `minact` untouched, so the next iteration works with an effective r = lhs - minact that has
    // drifted from the row's real one, and the case split above no longer applies to it.
    // Measured on icir97_tension: minact == lhs == 1.8866667666953063e1 exactly, two negative
    // coefficients dropped, lhs 1.8866667666953063e1 -> 1.4309442332440668e2, which cut off the
    // optimum 6375 and let the solver close the gap at 6394.
    if s2_is_ge(minact, *lhs) {
        return 0;
    }

    // terminate, because coefficient tightening cannot be performed; also
    // excludes the case in which no integral variable is present
    // for lhs terminate if minact + maxabsval < rowlhs
    if s2_is_lt(minact + maxabsval, *lhs) {
        return 0;
    }

    // SCIPsortDownRealInt(absvals, rowinds, *rownnz): sort by absval descending,
    // permuting rowinds with it. Tie order (col asc) is our definition.
    pdqsort_by(&mut inds[..*nnz], |a, b| {
        let av = if ctx.is_integral(*a) {
            vals[*a as usize].abs()
        } else {
            0.0
        };
        let bv = if ctx.is_integral(*b) {
            vals[*b as usize].abs()
        } else {
            0.0
        };
        if av != bv { av > bv } else { a < b }
    });

    // loop over the integral variables and try to tighten the coefficients; see cons_linear for more details
    let mut i: usize = 0;
    while i < *nnz {
        let idx = inds[i];
        let coefidx = idx as usize;

        debug_assert!(idx >= 0);
        debug_assert!(idx < ctx.num_col() as HighsInt);

        // due to sorting, we can exit if we reached a continuous variable: all further integral variables have 0 coefficients anyway
        if !ctx.is_integral(idx) {
            break;
        }

        debug_assert!(ctx.is_integral(idx));

        if vals[coefidx] < 0.0 && s2_is_ge(minact - vals[coefidx], *lhs) {
            let newcoef = minact - (*lhs);
            let ub = ctx.global_ub(idx);

            debug_assert!(
                s2_is_ge(newcoef + S2_EPSILON, vals[coefidx]) || s2_is_ge(newcoef, vals[coefidx])
            );
            debug_assert!(!s2_is_positive(newcoef));

            if newcoef > vals[coefidx] {
                let mut delta = HighsCDouble::from(newcoef) + (-vals[coefidx]);
                delta = delta * ub;

                let newlhs = delta + *lhs;

                *lhs = f64::from(newlhs);

                nchgcoefs += 1;

                if s2_is_negative(newcoef) {
                    minacttmp += delta;
                    minact = f64::from(minacttmp);
                    vals[coefidx] = newcoef;
                } else {
                    *nnz -= 1;
                    vals[coefidx] = 0.0;
                    inds[i] = inds[*nnz];
                    // the C for-loop's ++i runs on continue: the swapped-in entry is skipped
                    i += 1;
                    continue;
                }
            }
        } else if vals[coefidx] > 0.0 && s2_is_ge(minact + vals[coefidx], *lhs) {
            let newcoef = (*lhs) - minact;
            let lb = ctx.global_lb(idx);

            debug_assert!(
                s2_is_le(newcoef, vals[coefidx] + S2_EPSILON) || s2_is_le(newcoef, vals[coefidx])
            );
            debug_assert!(!s2_is_negative(newcoef));

            if newcoef < vals[coefidx] {
                let mut delta = HighsCDouble::from(newcoef) + (-vals[coefidx]);
                delta = delta * lb;

                let newlhs = delta + *lhs;

                *lhs = f64::from(newlhs);

                nchgcoefs += 1;

                if s2_is_positive(newcoef) {
                    minacttmp += delta;
                    minact = f64::from(minacttmp);
                    vals[coefidx] = newcoef;
                } else {
                    *nnz -= 1;
                    vals[coefidx] = 0.0;
                    inds[i] = inds[*nnz];
                    // the C for-loop's ++i runs on continue: the swapped-in entry is skipped
                    i += 1;
                    continue;
                }
            }
        } else {
            // due to sorting we can stop completely if the precondition was not fulfilled for this variable
            break;
        }

        i += 1;
    }

    nchgcoefs
}

/// Diagnostic helper for `tighten_coefs_traced`: the row as a sparse term list.
fn s2_trace_terms(vals: &[f64], inds: &[HighsInt], n: usize) -> Vec<(HighsInt, f64)> {
    inds[..n]
        .iter()
        .map(|&col| (col, vals[col as usize]))
        .collect()
}

/// Diagnostic wrapper around `tighten_coefs` (`RUSOLVE_S2_TRACE=1`, inert otherwise): when the
/// row was satisfied by the debug solution before the call and violates it afterwards, dumps the
/// full inputs (coefficients, global bounds, integrality, minact, maxabsval) so the individual
/// tightening steps can be checked by hand.
pub fn tighten_coefs_traced<C: S2Context>(
    ctx: &C,
    vals: &mut [f64],
    inds: &mut [HighsInt],
    nnz: &mut usize,
    lhs: &mut f64,
    dbg: &S2DebugSol,
) -> HighsInt {
    if !s2_trace_enabled() {
        return tighten_coefs(ctx, vals, inds, nnz, lhs);
    }
    let before = s2_trace_terms(vals, inds, *nnz);
    let lhs0 = *lhs;
    let v0 = s2_debugsol_viol_ge(dbg, &before, lhs0);
    // recompute the two quantities the routine derives, for the dump only
    let mut minactq = HighsCDouble::from(0.0);
    let mut maxabsval = 0.0f64;
    let mut skip = false;
    for &(col, coef) in &before {
        if coef > 0.0 {
            let lb = ctx.global_lb(col);
            if s2_is_infinity(-lb) {
                skip = true;
            }
            if ctx.is_integral(col) {
                maxabsval = cpp_max(maxabsval, coef);
            }
            minactq += HighsCDouble::from(lb) * coef;
        } else {
            let ub = ctx.global_ub(col);
            if s2_is_infinity(ub) {
                skip = true;
            }
            if ctx.is_integral(col) {
                maxabsval = cpp_max(maxabsval, -coef);
            }
            minactq += HighsCDouble::from(ub) * coef;
        }
    }
    let minact0 = f64::from(minactq);

    let n = tighten_coefs(ctx, vals, inds, nnz, lhs);

    let after = s2_trace_terms(vals, inds, *nnz);
    let v1 = s2_debugsol_viol_ge(dbg, &after, *lhs);
    if v0 <= S2_FEASTOL && v1 > S2_FEASTOL {
        eprintln!(
            "s2-trace tighten: INVALIDATES debug solution: viol {v0:e} -> {v1:e}; \
             lhs {lhs0:e} -> {:e}; nnz {} -> {}; nchg={n}; minact={minact0:e}; \
             maxabsval={maxabsval:e}; early_return={skip}",
            *lhs,
            before.len(),
            after.len(),
        );
        // Every input column, with the global bound the routine multiplies changes by, its
        // integrality, and the debug solution's value: enough to replay the loop by hand.
        let aftercols: Vec<HighsInt> = after.iter().map(|&(col, _)| col).collect();
        for &(col, coef) in &before {
            eprintln!(
                "  col {col} coef {coef:e} gbound [{:e}, {:e}] int {} x={:e}{}",
                ctx.global_lb(col),
                ctx.global_ub(col),
                ctx.is_integral(col),
                dbg.x
                    .as_deref()
                    .and_then(|x| x.get(col as usize))
                    .copied()
                    .unwrap_or(f64::NAN),
                if aftercols.contains(&col) {
                    ""
                } else {
                    "  <-- dropped"
                },
            );
        }
    }
    n
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::linearCombRows (lines 975-1008)
pub fn linear_comb_rows(row1: &mut S2ConflictRow, row2: &S2ConflictRow, scale: f64) {
    debug_assert!(scale > 0.0);

    // add conflict and reason conflict rows;
    for i in 0..row2.nnz() {
        let idx = row2.inds[i];
        debug_assert!(idx >= 0);
        let idxu = idx as usize;
        if s2_is_zero(row1.vals[idxu]) {
            row1.inds.push(idx);
            row1.vals[idxu] = scale * row2.vals[idxu];
        } else {
            row1.vals[idxu] += scale * row2.vals[idxu];
        }
    }
    row1.lhs += scale * row2.lhs;

    // remove variables with zero coefficient.
    remove_zero_vars(row1);
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::getQuotLargestSmallestCoef (lines 1647-1674)
pub fn get_quot_largest_smallest_coef(inds: &[HighsInt], vals: &[f64], nnz: usize) -> f64 {
    debug_assert!(!vals.is_empty() || nnz == 0);

    if nnz == 0 {
        return 0.0;
    }

    let mut minval: f64 = kHighsInf;
    let mut maxval: f64 = -kHighsInf;

    for i in 0..nnz {
        let idx = inds[i] as usize;
        minval = cpp_min(minval, vals[idx]);
        maxval = cpp_max(maxval, vals[idx]);
    }
    (maxval / minval).abs()
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::weakenVarConflictRow (lines 1511-1539)
pub fn weaken_var_conflict_row<C: S2Context>(ctx: &C, row: &mut S2ConflictRow, pos: usize) {
    let idx = row.inds[pos] as usize;

    debug_assert!(pos < row.nnz());
    debug_assert!(!s2_is_zero(row.vals[idx]));

    // weaken with global upper bound
    if s2_is_gt(row.vals[idx], 0.0) {
        row.lhs -= row.vals[idx] * ctx.global_ub(idx as HighsInt);
    }
    // weaken with global lower bound
    else if s2_is_lt(row.vals[idx], 0.0) {
        row.lhs -= row.vals[idx] * ctx.global_lb(idx as HighsInt);
    }

    row.vals[idx] = 0.0;
    row.inds.swap_remove(pos);
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::weakenConflictRow (lines 1543-1596)
pub fn weaken_conflict_row<C: S2Context>(
    ctx: &C,
    row: &mut S2ConflictRow,
    pos: HighsInt,
    fixsides: Option<&[i8]>,
) {
    let mut i: isize = row.nnz() as isize - 1;
    while i >= 0 {
        let iu = i as usize;
        let idx = row.inds[iu];
        let idxu = idx as usize;

        if row.vals[idxu] > 0.0 {
            let ub = s2_ub_at(ctx, idx, pos, true);

            if s2_is_eq(ub, ctx.global_ub(idx))
                && (fixsides.is_none() || fixsides.unwrap()[idxu] == 0)
            {
                weaken_var_conflict_row(ctx, row, iu);
            }
        } else {
            let lb = s2_lb_at(ctx, idx, pos, true);

            if s2_is_eq(lb, ctx.global_lb(idx))
                && (fixsides.is_none() || fixsides.unwrap()[idxu] == 0)
            {
                weaken_var_conflict_row(ctx, row, iu);
            }
        }

        i -= 1;
    }
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::weakenContinuousVarsConflictRow (lines 1600-1642)
pub fn weaken_continuous_vars_conflict_row<C: S2Context>(
    ctx: &C,
    row: &mut S2ConflictRow,
    residx: HighsInt,
) {
    let mut i: isize = row.nnz() as isize - 1;
    while i >= 0 {
        let iu = i as usize;
        let idx = row.inds[iu];
        let idxu = idx as usize;

        // avoid weakening the variable we are resolving on
        if idx == residx {
            // continue
        } else if row.vals[idxu] > 0.0 && !ctx.is_integral(idx) {
            weaken_var_conflict_row(ctx, row, iu);
        } else if row.vals[idxu] < 0.0 && !ctx.is_integral(idx) {
            weaken_var_conflict_row(ctx, row, iu);
        }

        i -= 1;
    }
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::computeScaleReason (lines 2165-2184)
pub fn compute_scale_reason(
    conflict: &S2ConflictRow,
    reason: &S2ConflictRow,
    residx: HighsInt,
) -> f64 {
    let coefconf = conflict.vals[residx as usize];
    let coefreas = reason.vals[residx as usize];

    debug_assert!(!s2_is_zero(coefreas) && !s2_is_zero(coefconf));
    debug_assert!(coefconf * coefreas < 0.0);

    (coefconf / coefreas).abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::conflict_resolution::{S2TestChg, S2TestCtx};
    use crate::prelude::*;

    fn test_ctx(n: usize, integral: Vec<bool>) -> S2TestCtx {
        let num_integral = integral.iter().filter(|&&b| b).count();
        S2TestCtx {
            glb: vec![0.0; n],
            gub: vec![1.0; n],
            integral,
            num_integral,
            stack: Vec::new(),
        }
    }

    #[test]
    fn s2_rowops_tighten_redundant_row_untouched() {
        // A row that is redundant w.r.t. the variable bounds (minact > lhs) must not be
        // tightened: SCIP's precondition (conflict_resolution.c:604) does not hold.
        // -4x0 >= -10 with x0 integer in [0,2]: minact = -4*2 = -8 > lhs = -10.
        // Without the guard the negative branch computes newcoef = minact - lhs = +2 (sign
        // flip), delta = 2 - (-4) = 6, newlhs = -10 + 6*2 = 2, and drops the only variable,
        // returning the infeasible empty row 0 >= 2, which cuts off every solution.
        let mut ctx = test_ctx(1, vec![true]);
        ctx.gub[0] = 2.0;
        let mut vals = vec![-4.0];
        let mut inds = vec![0];
        let mut nnz: usize = 1;
        let mut lhs = -10.0;
        let nchg = tighten_coefs(&ctx, &mut vals, &mut inds, &mut nnz, &mut lhs);
        assert_eq!(nchg, 0);
        assert_eq!(nnz, 1);
        assert_eq!(vals[0], -4.0);
        assert_eq!(lhs, -10.0);
    }

    #[test]
    fn s2_rowops_tighten_pos() {
        // 3x0 + x1 >= 1, binary. Gives x0 + x1 >= 1, nchg 1, lhs 1.0.
        let ctx = test_ctx(2, vec![true, true]);
        let mut vals = vec![3.0, 1.0];
        let mut inds = vec![0, 1];
        let mut nnz: usize = 2;
        let mut lhs = 1.0;
        let nchg = tighten_coefs(&ctx, &mut vals, &mut inds, &mut nnz, &mut lhs);
        assert_eq!(nchg, 1);
        assert_eq!(vals[0], 1.0);
        assert_eq!(lhs, 1.0);
        assert_eq!(&inds[..nnz], &[0, 1]);
    }

    #[test]
    fn s2_rowops_tighten_neg() {
        // -4x0 + x1 >= -1, binary. Gives -3x0 + x1 >= 0, nchg 1.
        let ctx = test_ctx(2, vec![true, true]);
        let mut vals = vec![-4.0, 1.0];
        let mut inds = vec![0, 1];
        let mut nnz: usize = 2;
        let mut lhs = -1.0;
        let nchg = tighten_coefs(&ctx, &mut vals, &mut inds, &mut nnz, &mut lhs);
        assert_eq!(nchg, 1);
        assert_eq!(vals[0], -3.0);
        assert_eq!(lhs, 0.0);
        assert_eq!(vals[1], 1.0);
    }

    #[test]
    fn s2_rowops_tighten_cont() {
        // x0 + 2y >= 1, x0 binary, y continuous [0,1]. Unchanged, nchg 0.
        let ctx = test_ctx(2, vec![true, false]);
        let mut vals = vec![1.0, 2.0];
        let mut inds = vec![0, 1];
        let mut nnz: usize = 2;
        let mut lhs = 1.0;
        let nchg = tighten_coefs(&ctx, &mut vals, &mut inds, &mut nnz, &mut lhs);
        assert_eq!(nchg, 0);
        assert_eq!(vals[0], 1.0);
        assert_eq!(vals[1], 2.0);
        assert_eq!(lhs, 1.0);
    }

    #[test]
    fn s2_rowops_tighten_quirk() {
        // x0 + x1 >= 0, binary: minact == lhs, so the row is redundant w.r.t. the bounds and
        // SCIP's precondition (conflict_resolution.c:604) does not hold.
        //
        // This fixture used to pin the "coefficient became zero" exit this input reaches
        // (inds [1], vals[0] == 0, nchg 1, lhs 0.0). That exit is dead code in SCIP: the assert
        // makes |newcoef| = |minact - lhs| > S2_EPSILON hold, and a drop moves lhs while leaving
        // minact stale, so the next iteration would tighten against a drifted r. Keeping the
        // exit reachable is what cut the optimum off icir97_tension (6394 reported against the
        // true 6375), so the port now refuses the input instead, as the assert does.
        let ctx = test_ctx(2, vec![true, true]);
        let mut vals = vec![1.0, 1.0];
        let mut inds = vec![0, 1];
        let mut nnz: usize = 2;
        let mut lhs = 0.0;
        let nchg = tighten_coefs(&ctx, &mut vals, &mut inds, &mut nnz, &mut lhs);
        assert_eq!(nchg, 0);
        assert_eq!(nnz, 2);
        assert_eq!(&inds[..nnz], &[0, 1]);
        assert_eq!(vals[0], 1.0);
        assert_eq!(vals[1], 1.0);
        assert_eq!(lhs, 0.0);
    }

    #[test]
    fn s2_rowops_lincomb() {
        // (x0 + x1 >= 1) + 1.0 * (-x0 + x2 >= 0) gives inds [2,1], vals[1] = 1, vals[2] = 1, lhs 1.0.
        let mut row1 = S2ConflictRow {
            vals: vec![1.0, 1.0, 0.0],
            inds: vec![0, 1],
            lhs: 1.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        let row2 = S2ConflictRow {
            vals: vec![-1.0, 0.0, 1.0],
            inds: vec![0, 2],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        linear_comb_rows(&mut row1, &row2, 1.0);
        assert_eq!(row1.inds, vec![2, 1]);
        assert_eq!(row1.vals[1], 1.0);
        assert_eq!(row1.vals[2], 1.0);
        assert_eq!(row1.vals[0], 0.0);
        assert_eq!(row1.lhs, 1.0);
    }

    #[test]
    fn s2_rowops_coefquot() {
        let inds = vec![0, 1];
        let vals = vec![1.0, 4.0];
        assert_eq!(get_quot_largest_smallest_coef(&inds, &vals, 2), 4.0);
        let vals2 = vec![2.0, -3.0];
        assert_eq!(
            get_quot_largest_smallest_coef(&inds, &vals2, 2),
            0.6666666666666666
        );
        assert_eq!(get_quot_largest_smallest_coef(&inds, &vals, 0), 0.0);
    }

    #[test]
    fn s2_rowops_scale() {
        let conflict = S2ConflictRow {
            vals: vec![2.0],
            inds: vec![0],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        let reason = S2ConflictRow {
            vals: vec![-4.0],
            inds: vec![0],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        assert_eq!(compute_scale_reason(&conflict, &reason, 0), 0.5);
    }

    #[test]
    fn s2_rowops_weaken() {
        // x0 + x1 >= 1, binary global [0,1], change x1 ub 1 to 0 at pos 0.
        // weaken_conflict_row(pos 0, None): x0 is weakened (ub at pos 0 is still 1 == gub),
        // x1 is NOT weakened (ub at pos 0 is 0 != gub 1).
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0],
            gub: vec![1.0, 1.0],
            integral: vec![true, true],
            num_integral: 2,
            stack: vec![S2TestChg {
                col: 1,
                is_upper: true,
                newbound: 0.0,
                depth: 1,
                branching: false,
                reason: None,
            }],
        };
        let mut row = S2ConflictRow {
            vals: vec![1.0, 1.0],
            inds: vec![0, 1],
            lhs: 1.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        weaken_conflict_row(&ctx, &mut row, 0, None);
        assert_eq!(row.lhs, 0.0);
        assert_eq!(row.nnz(), 1);
        assert_eq!(&row.inds[..row.nnz()], &[1]);
    }
}
