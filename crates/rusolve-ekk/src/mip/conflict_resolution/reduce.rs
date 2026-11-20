//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP source: conflict_resolution.c (isBinaryConflictRow 755-778, ComplementedMirLhs 821-971).

use super::*;
use crate::prelude::*;

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::isBinaryConflictRow (lines 755-778)
pub fn is_binary_conflict_row<C: S2Context>(ctx: &C, row: &S2ConflictRow) -> bool {
    for i in 0..row.nnz() {
        let col = row.inds[i];
        if !s2_is_binary(ctx, col) {
            return false;
        }
    }
    true
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::ComplementedMirLhs (lines 821-971)
pub fn complemented_mir_lhs<C: S2Context>(
    ctx: &C,
    row: &mut S2ConflictRow,
    fixsides: &[i8],
    pos: HighsInt,
    residx: HighsInt,
    divisor: f64,
) {
    debug_assert!(s2_is_gt(divisor, 0.0));

    if !is_binary_conflict_row(ctx, row) {
        return;
    }

    // first handle x_k and initialize lhs deltas
    let mut varidx = residx;
    let resvarcoef = row.vals[varidx as usize];
    let mut deltaoldlhs;
    let mut deltanewlhs;
    if resvarcoef < 0.0 {
        deltaoldlhs = -row.vals[varidx as usize];
        deltanewlhs = -1.0;
        row.vals[varidx as usize] = -1.0;
    } else {
        deltaoldlhs = 0.0;
        deltanewlhs = 0.0;
        row.vals[varidx as usize] = 1.0;
    }

    // compute the delta for the left hand side after complementation in order to apply MIR
    // In a second loop set the new coefficients for the other variables and compute the lhs delta after complementation
    for i in 0..row.nnz() {
        varidx = row.inds[i];

        if varidx == residx {
            continue;
        }

        let coef = row.vals[varidx as usize];

        // complementation offset
        if (coef > 0.0
            && (s2_ub_at(ctx, varidx, pos, true) > 0.5 && fixsides[varidx as usize] != 1))
            || (coef < 0.0
                && (s2_lb_at(ctx, varidx, pos, true) > 0.5 || fixsides[varidx as usize] == -1))
        {
            deltaoldlhs += -coef;
        }
    }
    let oldlhs = row.lhs;
    let mut newlhs = (oldlhs + deltaoldlhs) / divisor;
    let fraclhs = newlhs - s2_floor(newlhs);

    if fraclhs < S2_MINFRAC || fraclhs > S2_MAXFRAC {
        row.vals[residx as usize] = resvarcoef;
        return;
    }

    // set the new coefficients for the other variables and compute the lhs deltas
    for i in 0..row.nnz() {
        varidx = row.inds[i];

        if varidx == residx {
            continue;
        }

        let coef = row.vals[varidx as usize] / divisor;
        let fraccoef = coef - s2_floor(coef);

        if (coef > 0.0
            && (s2_ub_at(ctx, varidx, pos, true) > 0.5 && fixsides[varidx as usize] != 1))
            || (coef < 0.0
                && (s2_lb_at(ctx, varidx, pos, true) > 0.5 || fixsides[varidx as usize] == -1))
        {
            let newcoef;
            if (1.0 - fraccoef) >= fraclhs {
                newcoef = -s2_ceil(-coef);

                row.vals[varidx as usize] = newcoef;
                deltanewlhs += newcoef;
            } else {
                newcoef = -s2_floor(-coef) - (1.0 - fraccoef) / fraclhs;

                row.vals[varidx as usize] = newcoef;
                deltanewlhs += newcoef;
            }
        } else {
            if fraccoef >= fraclhs {
                row.vals[varidx as usize] = s2_ceil(coef);
            } else {
                row.vals[varidx as usize] = s2_floor(coef) + fraccoef / fraclhs;
            }
        }
    }

    newlhs = s2_ceil(newlhs) + deltanewlhs;
    row.lhs = newlhs;

    // remove variables with zero coefficient.
    remove_zero_vars(row);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[test]
    fn s2_reduce_cmir_basic() {
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0, 0.0],
            gub: vec![1.0, 1.0, 1.0],
            integral: vec![true, true, true],
            num_integral: 3,
            stack: vec![
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                S2TestChg {
                    col: 2,
                    is_upper: true,
                    newbound: 1.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 1.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
            ],
        };
        let mut row = S2ConflictRow {
            vals: vec![3.0, 2.0, 2.0],
            inds: vec![0, 1, 2],
            lhs: 3.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        complemented_mir_lhs(&ctx, &mut row, &[0i8, 0, 0], 2, 0, 3.0);
        assert_eq!(row.vals[0], 1.0);
        assert_eq!(row.vals[1], 1.0);
        assert_eq!(row.vals[2], 0.0);
        assert_eq!(row.inds, vec![0, 1]);
        assert_eq!(row.lhs, 1.0);
    }

    #[test]
    fn s2_reduce_cmir_nochange() {
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0, 0.0],
            gub: vec![1.0, 1.0, 1.0],
            integral: vec![true, true, true],
            num_integral: 3,
            stack: vec![],
        };
        let mut row = S2ConflictRow {
            vals: vec![2.0, 2.0, 3.0],
            inds: vec![0, 1, 2],
            lhs: 3.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        complemented_mir_lhs(&ctx, &mut row, &[0i8, 0, 0], kHighsIInf, 0, 2.0);
        assert_eq!(row.vals, vec![2.0, 2.0, 3.0]);
        assert_eq!(row.inds, vec![0, 1, 2]);
        assert_eq!(row.lhs, 3.0);
    }

    #[test]
    fn s2_reduce_not_binary() {
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0],
            gub: vec![2.0, 1.0],
            integral: vec![true, true],
            num_integral: 2,
            stack: vec![],
        };
        let mut row = S2ConflictRow {
            vals: vec![3.0, 2.0],
            inds: vec![0, 1],
            lhs: 3.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        assert!(!is_binary_conflict_row(&ctx, &row));
        complemented_mir_lhs(&ctx, &mut row, &[0i8, 0], kHighsIInf, 0, 3.0);
        assert_eq!(row.vals, vec![3.0, 2.0]);
        assert_eq!(row.inds, vec![0, 1]);
        assert_eq!(row.lhs, 3.0);
    }
}
