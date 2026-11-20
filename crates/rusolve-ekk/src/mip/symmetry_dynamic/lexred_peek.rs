//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! peekStaticLexredIsFeasible and propagateStaticLexred (static LexRed
//! propagation with peeking) behind `mip_symmetry_dynamic`.

use super::*;
use crate::mip::symmetry_dynamic::lexred_data::var_order_get_index;
use crate::mip::symmetry_dynamic::lexred_prop::{
    always_lt_shifted_vars, can_gt_shifted_vars, get_var_bounds, propagate_lower_bound_var,
    propagate_upper_bound_sym_var, propagate_variable_pair,
};

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::peekStaticLexredIsFeasible (lines 1076-1200)
///
/// Checks whether the static lexred with a certain variable ordering is feasible in the
/// hypothetical scenario where variables with indices `fixi` and `fixj` are fixed to
/// `fixvaluei` / `fixvaluej` (peeking). Returns `peekfeasible`.
///
/// Signed permutations are not supported: `centeri = centerj = 0.0` always, and the
/// `issigned` / `isnegated` branches are dropped.
pub(crate) fn peek_static_lexred_is_feasible(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    lexdata: &S3LexData,
    varorder: Option<&[HighsInt]>,
    nselvars: HighsInt,
    fixi: HighsInt,
    fixj: HighsInt,
    fixrow: HighsInt,
    fixvaluei: f64,
    fixvaluej: f64,
    peeklbs: &mut [f64],
    peekubs: &mut [f64],
    peekbdset: &mut [bool],
) -> bool {
    let nvars = lexdata.vars.len() as HighsInt;
    // SCIP: `SCIP_Bool infeasible = FALSE; int nreductions;` declared once for the whole
    // loop (only written in non-peek mode, so never read here).
    let mut infeasible = false;
    let mut nreductions: HighsInt = 0;

    debug_assert!(nvars >= 0);
    debug_assert!(nselvars <= nvars);
    debug_assert!(nselvars > 0);
    debug_assert!(fixi >= 0);
    debug_assert!(fixi < nvars);
    debug_assert!(fixj < nvars);
    debug_assert!(fixi != fixj);
    debug_assert!(fixrow >= 0);
    debug_assert!(fixrow < nselvars);
    debug_assert!(fixi == var_order_get_index(varorder, fixrow));
    debug_assert!(fixj == lexdata.invperm[var_order_get_index(varorder, fixrow) as usize]);
    debug_assert!(fixi == lexdata.perm[fixj as usize]);

    let mut peekfeasible = true;

    // intialize peekbdset
    for i in 0..nvars {
        peekbdset[i as usize] = false;
    }

    peeklbs[fixi as usize] = fixvaluei;
    peeklbs[fixj as usize] = fixvaluej;
    peekubs[fixi as usize] = fixvaluei;
    peekubs[fixj as usize] = fixvaluej;
    peekbdset[fixi as usize] = true;
    peekbdset[fixj as usize] = true;

    let mut row = fixrow + 1;
    while row < nselvars {
        // get left and right column indices
        let i = var_order_get_index(varorder, row);
        let j = lexdata.invperm[i as usize];
        debug_assert!(i == lexdata.perm[j as usize]);

        // no fixed points
        debug_assert!(i != j);

        debug_assert!(0 <= i && i < nvars);
        debug_assert!(0 <= j && j < nvars);

        let vari = lexdata.vars[i as usize];
        let varj = lexdata.vars[j as usize];

        // Signed permutations not supported: centeri = centerj = 0.0 always.
        let centeri = 0.0;
        let centerj = 0.0;

        // propagate that vari >= varj

        // vari >= varj can never hold if the maximal value of vari is smaller than the
        // minimal value of varj
        if always_lt_shifted_vars(
            ms, localdom, vari, varj, centeri, centerj, false, true, i as usize, j as usize,
            peeklbs, peekubs, peekbdset,
        ) {
            peekfeasible = false;
            break;
        }

        propagate_lower_bound_var(
            ms,
            localdom,
            vari,
            varj,
            centeri,
            centerj,
            false,
            &mut infeasible,
            &mut nreductions,
            true,
            i as usize,
            j as usize,
            peeklbs,
            peekubs,
            peekbdset,
        );

        propagate_upper_bound_sym_var(
            ms,
            localdom,
            vari,
            varj,
            centeri,
            centerj,
            false,
            &mut infeasible,
            &mut nreductions,
            true,
            i as usize,
            j as usize,
            peeklbs,
            peekubs,
            peekbdset,
        );

        // if there exists a solution with vari > varj, reductions are feasible w.r.t. lexred
        if can_gt_shifted_vars(
            ms, localdom, vari, varj, centeri, centerj, false, true, i as usize, j as usize,
            peeklbs, peekubs, peekbdset,
        ) {
            break;
        }

        row += 1;
    }

    peekfeasible
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::propagateStaticLexred (lines 1204-1514)
///
/// Propagates static lexicographic reduction with the specified variable ordering.
///
/// Signed permutations are not supported: the `vari == varj` self-reflection branch
/// and the `isnegated` / `issigned` handling are dropped; `centeri = centerj = 0.0` always.
pub(crate) fn propagate_static_lexred(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    lexdata: &S3LexData,
    varorder: Option<&[HighsInt]>,
    nselvars: HighsInt,
    infeasible: &mut bool,
    nreductions: &mut HighsInt,
) {
    let nvars = lexdata.vars.len() as HighsInt;
    let eps = ms.mipdata_.epsilon;
    let feastol = ms.mipdata_.feastol;

    debug_assert!(nselvars >= 0);
    debug_assert!(!*infeasible);
    debug_assert!(*nreductions >= 0);

    // avoid trivial cases
    if nselvars <= 0 {
        return;
    }

    // iterate over the variable array entrywise
    //
    // We see this as two columns, with the left vector being the variable ordering,
    // and the right column the permuted variables of this var ordering.
    let mut row: HighsInt = 0;
    let mut i: HighsInt = -1;
    let mut j: HighsInt = -1;
    let mut vari: HighsInt = -1;
    let mut varj: HighsInt = -1;

    while row < nselvars {
        // left and right column indices
        i = var_order_get_index(varorder, row);
        j = lexdata.invperm[i as usize];
        debug_assert!(i == lexdata.perm[j as usize]);

        // no fixed points
        debug_assert!(i != j);

        debug_assert!(0 <= i && i < nvars);
        debug_assert!(0 <= j && j < nvars);

        vari = lexdata.vars[i as usize];
        varj = lexdata.vars[j as usize];

        // Signed permutations not supported: centeri = centerj = 0.0 always.
        let centeri = 0.0;
        let centerj = 0.0;

        propagate_variable_pair(
            ms,
            localdom,
            vari,
            varj,
            centeri,
            centerj,
            false,
            infeasible,
            nreductions,
            false,
            0,
            0,
            &mut [],
            &mut [],
            &mut [],
        );

        if *infeasible {
            return;
        }

        // terminate if there exists a solution being lexicographically strictly larger than its image
        if can_gt_shifted_vars(
            ms,
            localdom,
            vari,
            varj,
            centeri,
            centerj,
            false,
            false,
            0,
            0,
            &[],
            &[],
            &[],
        ) {
            break;
        }

        row += 1;
    }

    debug_assert!(vari != -1);
    debug_assert!(varj != -1);
    debug_assert!(0 <= i && i < nvars);
    debug_assert!(0 <= j && j < nvars);

    // if the variables in "row" are fixed to the same value, we might find further propagations
    if row < nselvars {
        let (lb1, ub1, lb2, ub2) = get_var_bounds(localdom, vari, varj, false, 0, 0, &[], &[], &[]);

        // Signed permutations not supported: the vari == varj self-reflection branch
        // is dropped; vari != varj always (no fixed points).

        // The previous loop is broken at row "row", which allows for choosing vari > varj.
        //
        // Now check if vari == varj is permitted, and if not, tighten the domain further.
        //
        // Case 1: vari is minimal (lbi).
        // Then, propagation of lbi = vari >= varj can yield two situations:
        //   Option 1: varj can take a value < lbi. Then no further reductions can be detected.
        //   Option 2: varj gets fixed to lbi. Then, we must check if feasibility is found, still.
        //     If it turns out infeasible, then we know vari cannot take value lbi, so we can increase the lower bound.

        // Signed permutations not supported: centeri = centerj = 0.0 always.
        let centeri = 0.0;
        let centerj = 0.0;

        // translate variable bounds to shifted variable domain; isnegated is always false
        let lbi = lb1 - centeri;
        let ubi = ub1 - centeri;
        let lbj = lb2 - centerj;
        let ubj = ub2 - centerj;

        // check whether peek is called
        if !sym_eq(lbi, lbj, eps, kHighsInf) && !sym_eq(ubi, ubj, eps, kHighsInf) {
            return;
        }

        let mut peeklbs = vec![0.0; nvars as usize];
        let mut peekubs = vec![0.0; nvars as usize];
        let mut peekbdset = vec![false; nvars as usize];

        if sym_eq(lbj, lbi, eps, kHighsInf) {
            // translate lbj back to original variable domain of variable j; isnegated false
            let fixvalj = lbj + centerj;

            // this is Option 2: varj gets fixed to lbi by propagation.
            let peekfeasible = peek_static_lexred_is_feasible(
                ms,
                localdom,
                lexdata,
                varorder,
                nselvars,
                i,
                j,
                row,
                lbi + centeri,
                fixvalj,
                &mut peeklbs,
                &mut peekubs,
                &mut peekbdset,
            );
            if !peekfeasible {
                // vari cannot take value lb1, so we increase the lower bound (do not use lbi as this is a shifted domain bound)
                if ms.is_col_integral(vari) {
                    // discrete variable type: increase lower bound by 1.
                    debug_assert!(is_integral(lb1, eps));
                    let (inf, success) =
                        s3_tighten_lb(ms, localdom, vari, lb1 + 1.0, true, feastol);
                    if success {
                        *nreductions += 1;
                    }
                    if inf {
                        *infeasible = true;
                        return;
                    }
                    let lb1 = lb1 + 1.0;
                    debug_assert!(sym_le(lb1, ub1, eps, kHighsInf));
                } else {
                    // continuous variable type: act as if we increase the variable by a very little bit.
                    // That is only possible if we're able to increase the variable bound by a bit.
                    if sym_eq(lbi, ubi, eps, kHighsInf) {
                        *infeasible = true;
                        return;
                    }
                }
            }
        }

        // Case 2: varj is maximal (ubj).
        // Then, propagation of vari >= varj = ubj can yield two situations:
        //   Option 1: vari can take a value > ubj. Then, no further reductions can be detected.
        //   Option 2: vari gets fixed to ubj. Then, we must check if feasibility is found, still.
        //     If it turns out infeasible, then we know varj cannot take value ubj, so we can decrease the upper bound.
        debug_assert!(sym_ge(ubi, ubj, eps, kHighsInf)); // this must be the case after reductions in the for-loop
        if sym_eq(ubi, ubj, eps, kHighsInf) {
            // translate ubj back to original variable domain of variable j; isnegated false
            let fixvalj = ubj + centerj;

            // this is Option 2: vari gets fixed to ubj by propagation.
            let peekfeasible = peek_static_lexred_is_feasible(
                ms,
                localdom,
                lexdata,
                varorder,
                nselvars,
                i,
                j,
                row,
                ubi + centeri,
                fixvalj,
                &mut peeklbs,
                &mut peekubs,
                &mut peekbdset,
            );
            if !peekfeasible {
                // varj cannot take value ub2, so we decrease the upper bound (do not use ubj as this is a shifted domain bound)
                if ms.is_col_integral(varj) {
                    // discrete variable type: decrease upper bound by 1. isnegated false
                    debug_assert!(is_integral(ub2, eps));
                    let (inf, success) =
                        s3_tighten_ub(ms, localdom, varj, ub2 - 1.0, true, feastol);
                    if success {
                        *nreductions += 1;
                    }
                    if inf {
                        *infeasible = true;
                        return;
                    }
                    let ubj = ubj - 1.0;
                    debug_assert!(sym_le(lbj, ubj, eps, kHighsInf));
                } else {
                    // continuous variable type: act as if we decrease the variable by a very little bit.
                    // that is only possible if we're able to decrease the variable bound by a bit.
                    if sym_eq(lbj, ubj, eps, kHighsInf) {
                        *infeasible = true;
                    }
                }
            }
        }
    }
}

/// `SCIPisIntegral` (EPSISINT), used only in debug asserts.
fn is_integral(v: f64, eps: f64) -> bool {
    (v - v.round()).abs() <= eps
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int_solver(num_col: HighsInt) -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = num_col;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0; num_col as usize];
        lp.col_lower_ = vec![0.0; num_col as usize];
        lp.col_upper_ = vec![10.0; num_col as usize];
        lp.integrality_ = vec![HighsVarType::kInteger; num_col as usize];
        lp.a_matrix_.num_col_ = num_col;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0; num_col as usize + 1];
        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0)
    }

    fn make_domain(lb0: f64, ub0: f64, lb1: f64, ub1: f64) -> HighsDomain {
        let mut localdom = HighsDomain::default();
        localdom.col_lower_ = vec![lb0, lb1];
        localdom.col_upper_ = vec![ub0, ub1];
        localdom.col_lower_pos_ = vec![-1; 2];
        localdom.col_upper_pos_ = vec![-1; 2];
        localdom.changedcolsflags_ = vec![0; 2];
        localdom
    }

    fn make_domain3(lb0: f64, ub0: f64, lb1: f64, ub1: f64, lb2: f64, ub2: f64) -> HighsDomain {
        let mut localdom = HighsDomain::default();
        localdom.col_lower_ = vec![lb0, lb1, lb2];
        localdom.col_upper_ = vec![ub0, ub1, ub2];
        localdom.col_lower_pos_ = vec![-1; 3];
        localdom.col_upper_pos_ = vec![-1; 3];
        localdom.changedcolsflags_ = vec![0; 3];
        localdom
    }

    fn swap_lexdata() -> S3LexData {
        S3LexData {
            isdynamic: false,
            vars: vec![0, 1],
            perm: vec![1, 0],
            invperm: vec![1, 0],
            varmap: std::collections::HashMap::new(),
            symvarmap: vec![],
            nsymvars: 0,
            hasdynamicperm: false,
        }
    }

    fn perm3_lexdata() -> S3LexData {
        S3LexData {
            isdynamic: false,
            vars: vec![0, 1, 2],
            perm: vec![1, 2, 0],
            invperm: vec![2, 0, 1],
            varmap: std::collections::HashMap::new(),
            symvarmap: vec![],
            nsymvars: 0,
            hasdynamicperm: false,
        }
    }

    /// 2 integer cols in [0,1], swap permutation, peek fixi=0,fixj=1,fixrow=0, values 0,0 -> true.
    #[test]
    fn s3_peek_feasible_equal_fixing() {
        let mut ms = int_solver(2);
        let mut localdom = make_domain(0.0, 1.0, 0.0, 1.0);
        let lexdata = swap_lexdata();
        let nvars = 2;
        let mut peeklbs = vec![0.0; nvars];
        let mut peekubs = vec![0.0; nvars];
        let mut peekbdset = vec![false; nvars];
        let feasible = peek_static_lexred_is_feasible(
            &mut ms,
            &mut localdom,
            &lexdata,
            None,
            2,
            0,
            1,
            0,
            0.0,
            0.0,
            &mut peeklbs,
            &mut peekubs,
            &mut peekbdset,
        );
        assert!(feasible);
    }

    /// 3 integer cols, perm=[1,2,0], invperm=[2,0,1], varorder None, col1 bounds [0,0],
    /// col0 and col2 bounds [0,1], peek fixi=0,fixj=2,fixrow=0 values 1,1 -> false.
    /// row 1 is vari=x1 (ub 0), varj=x0 (peek lb 1): alwaysLTshiftedVars -> infeasible.
    #[test]
    fn s3_peek_infeasible_row_after_fix() {
        let mut ms = int_solver(3);
        let mut localdom = make_domain3(0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let lexdata = perm3_lexdata();
        let nvars = 3;
        let mut peeklbs = vec![0.0; nvars];
        let mut peekubs = vec![0.0; nvars];
        let mut peekbdset = vec![false; nvars];
        let feasible = peek_static_lexred_is_feasible(
            &mut ms,
            &mut localdom,
            &lexdata,
            None,
            3,
            0,
            2,
            0,
            1.0,
            1.0,
            &mut peeklbs,
            &mut peekubs,
            &mut peekbdset,
        );
        assert!(!feasible);
    }

    /// 2 integer cols, swap, col0 in [0,1], col1 in [1,1], propagate_static_lexred ->
    /// col0 lb becomes 1, nreductions 1, infeasible false.
    #[test]
    fn s3_peek_static_propagation_tightens_lb() {
        let mut ms = int_solver(2);
        // Attach the clique table, implications and domain, as the clique-table tests do,
        // because a tightening of an integer column runs add_implications.
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.post_solve_stack.initialize_index_maps(0, 2);
        ms.mipdata_.cliquetable = crate::mip::Slot::new(HighsCliqueTable::new(2));
        let imp = crate::mip::implications::HighsImplications::new(&ms);
        *ms.mipdata_.implications = imp;
        ms.mipdata_.domains = vec![crate::mip::Slot::new(HighsDomain::new(&ms))];
        let mut localdom = make_domain(0.0, 1.0, 1.0, 1.0);
        let lexdata = swap_lexdata();
        let mut infeasible = false;
        let mut nreductions: HighsInt = 0;
        propagate_static_lexred(
            &mut ms,
            &mut localdom,
            &lexdata,
            None,
            2,
            &mut infeasible,
            &mut nreductions,
        );
        assert!(!infeasible);
        assert_eq!(nreductions, 1);
        assert_eq!(localdom.col_lower_[0], 1.0);
    }

    /// 2 integer cols, swap, col0 in [0,0], col1 in [1,1] -> infeasible true.
    #[test]
    fn s3_peek_static_detects_infeasible() {
        let mut ms = int_solver(2);
        let mut localdom = make_domain(0.0, 0.0, 1.0, 1.0);
        let lexdata = swap_lexdata();
        let mut infeasible = false;
        let mut nreductions: HighsInt = 0;
        propagate_static_lexred(
            &mut ms,
            &mut localdom,
            &lexdata,
            None,
            2,
            &mut infeasible,
            &mut nreductions,
        );
        assert!(infeasible);
    }
}
