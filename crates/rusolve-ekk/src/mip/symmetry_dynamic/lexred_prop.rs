//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! Ported from SCIP.

use super::*;

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::getVarBounds (lines 590-637)
pub(crate) fn get_var_bounds(
    localdom: &HighsDomain,
    var1: HighsInt,
    var2: HighsInt,
    peekmode: bool,
    varidx1: usize,
    varidx2: usize,
    peeklbs: &[f64],
    peekubs: &[f64],
    peekbdset: &[bool],
) -> (f64, f64, f64, f64) {
    let (lb1, ub1) = if peekmode && peekbdset[varidx1] {
        (peeklbs[varidx1], peekubs[varidx1])
    } else {
        (
            localdom.col_lower_[var1 as usize],
            localdom.col_upper_[var1 as usize],
        )
    };

    let (lb2, ub2) = if peekmode && peekbdset[varidx2] {
        (peeklbs[varidx2], peekubs[varidx2])
    } else {
        (
            localdom.col_lower_[var2 as usize],
            localdom.col_upper_[var2 as usize],
        )
    };

    (lb1, ub1, lb2, ub2)
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::alwaysLTshiftedVars (lines 644-683)
pub(crate) fn always_lt_shifted_vars(
    ms: &HighsMipSolver,
    localdom: &HighsDomain,
    var1: HighsInt,
    var2: HighsInt,
    shift1: f64,
    shift2: f64,
    isnegated: bool,
    peekmode: bool,
    varidx1: usize,
    varidx2: usize,
    peeklbs: &[f64],
    peekubs: &[f64],
    peekbdset: &[bool],
) -> bool {
    let eps = ms.mipdata_.epsilon;

    let (lb1, ub1, lb2, ub2) = get_var_bounds(
        localdom, var1, var2, peekmode, varidx1, varidx2, peeklbs, peekubs, peekbdset,
    );

    // for negated variables, check: var1 - shift1 < shift2 - var2
    if isnegated && scip_is_lt(ub1, shift1 + shift2 - ub2, eps) {
        return true;
    }

    // for non-negated variables, check: var1 - center1 < var2 - center2
    if (!isnegated) && scip_is_lt(ub1, shift1 - shift2 + lb2, eps) {
        return true;
    }

    false
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::canGTshiftedVars (lines 691-730)
pub(crate) fn can_gt_shifted_vars(
    ms: &HighsMipSolver,
    localdom: &HighsDomain,
    var1: HighsInt,
    var2: HighsInt,
    shift1: f64,
    shift2: f64,
    isnegated: bool,
    peekmode: bool,
    varidx1: usize,
    varidx2: usize,
    peeklbs: &[f64],
    peekubs: &[f64],
    peekbdset: &[bool],
) -> bool {
    let eps = ms.mipdata_.epsilon;

    let (lb1, ub1, lb2, ub2) = get_var_bounds(
        localdom, var1, var2, peekmode, varidx1, varidx2, peeklbs, peekubs, peekbdset,
    );

    // for negated variables, check: var1 - shift1 > shift2 - var2
    if isnegated && scip_is_gt(ub1, shift1 + shift2 - ub2, eps) {
        return true;
    }

    // for non-negated variables, check: var1 - center1 > var2 - center2
    if (!isnegated) && scip_is_gt(ub1, shift1 - shift2 + lb2, eps) {
        return true;
    }

    false
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::propagateLowerBoundVar (lines 735-817)
pub(crate) fn propagate_lower_bound_var(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    var1: HighsInt,
    var2: HighsInt,
    center1: f64,
    center2: f64,
    isnegated: bool,
    infeasible: &mut bool,
    nreductions: &mut HighsInt,
    peekmode: bool,
    varidx1: usize,
    varidx2: usize,
    peeklbs: &mut [f64],
    peekubs: &mut [f64],
    peekbdset: &mut [bool],
) {
    let eps = ms.mipdata_.epsilon;

    let (lb1, ub1, lb2, ub2) = get_var_bounds(
        localdom, var1, var2, peekmode, varidx1, varidx2, peeklbs, peekubs, peekbdset,
    );

    let mut tighten = false;
    let mut newbd = 0.0;

    // tighten domain of var1 to ensure that var1 - center1 >= isnegated * (var2 - center2)
    if isnegated {
        if scip_is_lt(lb1 - center1, center2 - ub2, eps) {
            tighten = true;
            newbd = center1 + center2 - ub2;
        }
    } else {
        if scip_is_lt(lb1 - center1, lb2 - center2, eps) {
            tighten = true;
            newbd = center1 + lb2 - center2;
        }
    }

    if tighten {
        // in peek mode, only store updated bounds
        if peekmode {
            peeklbs[varidx1] = newbd;
            peekubs[varidx1] = ub1;
            peekbdset[varidx1] = true;
        } else {
            let is_int = ms.is_col_integral(var1);
            let feastol = ms.mipdata_.feastol;
            let (inf, tight) = s3_tighten_lb(ms, localdom, var1, newbd, is_int, feastol);
            *infeasible = inf;
            if tight {
                *nreductions += 1;
            }
        }
    }
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::propagateUpperBoundSymVar (lines 822-927)
pub(crate) fn propagate_upper_bound_sym_var(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    var1: HighsInt,
    var2: HighsInt,
    center1: f64,
    center2: f64,
    isnegated: bool,
    infeasible: &mut bool,
    nreductions: &mut HighsInt,
    peekmode: bool,
    varidx1: usize,
    varidx2: usize,
    peeklbs: &mut [f64],
    peekubs: &mut [f64],
    peekbdset: &mut [bool],
) {
    let eps = ms.mipdata_.epsilon;

    let (lb1, ub1, lb2, ub2) = get_var_bounds(
        localdom, var1, var2, peekmode, varidx1, varidx2, peeklbs, peekubs, peekbdset,
    );

    let mut tighten = false;
    let mut newbd = 0.0;

    // tighten domain of var2 to ensure that var1 - center1 >= isnegated * (var2 - center2)
    if isnegated {
        if scip_is_lt(ub1 - center1, center2 - lb2, eps) {
            tighten = true;
            newbd = center1 + center2 - ub1;
        }
    } else {
        if scip_is_lt(ub1 - center1, ub2 - center2, eps) {
            tighten = true;
            newbd = center2 - center1 + ub1;
        }
    }

    if tighten {
        // in peek mode, only store updated bounds
        if peekmode {
            if isnegated {
                peeklbs[varidx2] = newbd;
                peekubs[varidx2] = ub2;
            } else {
                peeklbs[varidx2] = lb2;
                peekubs[varidx2] = newbd;
            }
            peekbdset[varidx2] = true;
        } else {
            let is_int = ms.is_col_integral(var2);
            let feastol = ms.mipdata_.feastol;
            let (inf, tight) = if isnegated {
                s3_tighten_lb(ms, localdom, var2, newbd, is_int, feastol)
            } else {
                s3_tighten_ub(ms, localdom, var2, newbd, is_int, feastol)
            };
            *infeasible = inf;
            if tight {
                *nreductions += 1;
            }
        }
    }
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::propagateVariablePair (lines 1003-1070)
pub(crate) fn propagate_variable_pair(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    var1: HighsInt,
    var2: HighsInt,
    center1: f64,
    center2: f64,
    isnegated: bool,
    infeasible: &mut bool,
    nreductions: &mut HighsInt,
    peekmode: bool,
    varidx1: usize,
    varidx2: usize,
    peeklbs: &mut [f64],
    peekubs: &mut [f64],
    peekbdset: &mut [bool],
) {
    // perform lexicographic comparison: var1 - center1 >= +/- (var2 - center2)
    if always_lt_shifted_vars(
        ms, localdom, var1, var2, center1, center2, isnegated, peekmode, varidx1, varidx2, peeklbs,
        peekubs, peekbdset,
    ) {
        // SCIP_DEBUG block (lines 1030-1046) dropped
        *infeasible = true;
        return;
    }

    // for signed permutations, a variable might be mapped to itself
    if var1 == var2 {
        // self-reflection omitted
        debug_assert!(!isnegated);
    } else {
        propagate_lower_bound_var(
            ms,
            localdom,
            var1,
            var2,
            center1,
            center2,
            isnegated,
            infeasible,
            nreductions,
            peekmode,
            varidx1,
            varidx2,
            peeklbs,
            peekubs,
            peekbdset,
        );
        if *infeasible {
            return;
        }

        propagate_upper_bound_sym_var(
            ms,
            localdom,
            var1,
            var2,
            center1,
            center2,
            isnegated,
            infeasible,
            nreductions,
            peekmode,
            varidx1,
            varidx2,
            peeklbs,
            peekubs,
            peekbdset,
        );
    }
}

fn scip_is_lt(a: f64, b: f64, eps: f64) -> bool {
    a - b < -eps
}

fn scip_is_gt(a: f64, b: f64, eps: f64) -> bool {
    a - b > eps
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_col_continuous_solver() -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0; 2];
        lp.col_lower_ = vec![0.0; 2];
        lp.col_upper_ = vec![10.0; 2];
        lp.integrality_ = vec![HighsVarType::kContinuous; 2];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0; 3];
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

    #[test]
    fn s3_lxp_get_var_bounds_reads_peek() {
        let localdom = make_domain(0.0, 10.0, 2.0, 5.0);
        let peeklbs = vec![3.0, 0.0];
        let peekubs = vec![7.0, 0.0];
        let peekbdset = vec![true, false];

        let (lb1, ub1, lb2, ub2) =
            get_var_bounds(&localdom, 0, 1, true, 0, 1, &peeklbs, &peekubs, &peekbdset);

        assert_eq!(lb1, 3.0);
        assert_eq!(ub1, 7.0);
        assert_eq!(lb2, 2.0);
        assert_eq!(ub2, 5.0);
    }

    #[test]
    fn s3_lxp_always_lt_detects_infeasible() {
        let mut ms = two_col_continuous_solver();
        let mut localdom = make_domain(0.0, 1.0, 5.0, 6.0);
        let mut infeasible = false;
        let mut nreductions: HighsInt = 0;
        let mut peeklbs = vec![0.0; 2];
        let mut peekubs = vec![0.0; 2];
        let mut peekbdset = vec![false; 2];

        propagate_variable_pair(
            &mut ms,
            &mut localdom,
            0,
            1,
            0.0,
            0.0,
            false,
            &mut infeasible,
            &mut nreductions,
            false,
            0,
            1,
            &mut peeklbs,
            &mut peekubs,
            &mut peekbdset,
        );

        assert!(infeasible);
        assert_eq!(localdom.col_lower_[0], 0.0);
        assert_eq!(localdom.col_upper_[0], 1.0);
        assert_eq!(localdom.col_lower_[1], 5.0);
        assert_eq!(localdom.col_upper_[1], 6.0);
    }

    #[test]
    fn s3_lxp_propagate_pair_tightens_lb() {
        let mut ms = two_col_continuous_solver();
        let mut localdom = make_domain(0.0, 10.0, 2.0, 5.0);
        let mut infeasible = false;
        let mut nreductions: HighsInt = 0;
        let mut peeklbs = vec![0.0; 2];
        let mut peekubs = vec![0.0; 2];
        let mut peekbdset = vec![false; 2];

        propagate_variable_pair(
            &mut ms,
            &mut localdom,
            0,
            1,
            0.0,
            0.0,
            false,
            &mut infeasible,
            &mut nreductions,
            false,
            0,
            1,
            &mut peeklbs,
            &mut peekubs,
            &mut peekbdset,
        );

        assert!(!infeasible);
        assert_eq!(localdom.col_lower_[0], 2.0);
        assert_eq!(localdom.col_upper_[0], 10.0);
        assert_eq!(localdom.col_lower_[1], 2.0);
        assert_eq!(localdom.col_upper_[1], 5.0);
        assert_eq!(nreductions, 1);
    }

    #[test]
    fn s3_lxp_propagate_pair_negated_tightens_lb() {
        let mut ms = two_col_continuous_solver();
        let mut localdom = make_domain(0.0, 2.0, 0.0, 4.0);
        let mut infeasible = false;
        let mut nreductions: HighsInt = 0;
        let mut peeklbs = vec![0.0; 2];
        let mut peekubs = vec![0.0; 2];
        let mut peekbdset = vec![false; 2];

        propagate_variable_pair(
            &mut ms,
            &mut localdom,
            0,
            1,
            0.0,
            3.0,
            true,
            &mut infeasible,
            &mut nreductions,
            false,
            0,
            1,
            &mut peeklbs,
            &mut peekubs,
            &mut peekbdset,
        );

        assert!(!infeasible);
        assert_eq!(localdom.col_lower_[0], 0.0);
        assert_eq!(localdom.col_upper_[0], 2.0);
        assert_eq!(localdom.col_lower_[1], 1.0);
        assert_eq!(localdom.col_upper_[1], 4.0);
        assert_eq!(nreductions, 1);
    }

    #[test]
    fn s3_lxp_propagate_pair_peek_stores_bounds() {
        let mut ms = two_col_continuous_solver();
        let mut localdom = make_domain(0.0, 10.0, 2.0, 5.0);
        let mut infeasible = false;
        let mut nreductions: HighsInt = 0;
        let mut peeklbs = vec![0.0; 2];
        let mut peekubs = vec![0.0; 2];
        let mut peekbdset = vec![false; 2];

        propagate_variable_pair(
            &mut ms,
            &mut localdom,
            0,
            1,
            0.0,
            0.0,
            false,
            &mut infeasible,
            &mut nreductions,
            true,
            0,
            1,
            &mut peeklbs,
            &mut peekubs,
            &mut peekbdset,
        );

        assert_eq!(peeklbs[0], 2.0);
        assert_eq!(peekubs[0], 10.0);
        assert!(peekbdset[0]);
        assert_eq!(localdom.col_lower_[0], 0.0);
        assert_eq!(localdom.col_upper_[0], 10.0);
        assert_eq!(nreductions, 0);
    }

    #[test]
    fn s3_lxp_propagate_pair_no_change() {
        let mut ms = two_col_continuous_solver();
        let mut localdom = make_domain(0.0, 10.0, 0.0, 10.0);
        let mut infeasible = false;
        let mut nreductions: HighsInt = 0;
        let mut peeklbs = vec![0.0; 2];
        let mut peekubs = vec![0.0; 2];
        let mut peekbdset = vec![false; 2];

        propagate_variable_pair(
            &mut ms,
            &mut localdom,
            0,
            1,
            0.0,
            0.0,
            false,
            &mut infeasible,
            &mut nreductions,
            false,
            0,
            1,
            &mut peeklbs,
            &mut peekubs,
            &mut peekbdset,
        );

        assert!(!infeasible);
        assert_eq!(localdom.col_lower_[0], 0.0);
        assert_eq!(localdom.col_upper_[0], 10.0);
        assert_eq!(localdom.col_lower_[1], 0.0);
        assert_eq!(localdom.col_upper_[1], 10.0);
        assert_eq!(nreductions, 0);
    }
}
