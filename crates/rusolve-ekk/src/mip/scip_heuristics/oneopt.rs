//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! Ported from SCIP.

use super::sort::s6_sort_real_ptr;
use super::*;
use crate::lp_data::hconst::HighsVarType;
use crate::mip::HighsMipSolver;
use crate::mip::mip_solver_data::kSolutionSourceHeuristic;

/// SCIP 11.0.0 (32dd425): heur_oneopt.c::calcShiftVal (lines 99-209)
///
/// Computes the maximal shift value for a variable such that all rows stay valid.
/// Returns the shift value (newsolval - solval), or 0.0 if no shift is possible.
pub(crate) fn s6_calc_shift_val(
    ms: &HighsMipSolver,
    col: HighsInt,
    solval: f64,
    activities: &[f64],
) -> f64 {
    let feastol = ms.mipdata_.feastol;
    let model = mip_model!(ms);
    let domains = &ms.mipdata_.domains[0];

    // C lines 120-122: get variable's global bounds and objective coefficient
    let lb = domains.col_lower_[col as usize];
    let ub = domains.col_upper_[col as usize];
    let obj = model.col_cost_[col as usize];
    let mut shiftdown = true;

    // C lines 126-134: determine shifting direction and maximal possible shifting w.r.t. corresponding bound
    let mut newsolval;
    if obj > 0.0 {
        newsolval = s6_feas_ceil(lb, feastol);
    } else if obj < 0.0 {
        newsolval = s6_feas_floor(ub, feastol);
        shiftdown = false;
    } else {
        newsolval = solval;
    }

    // C lines 137-139: newsolval might be bounding solval -> avoid numerical shifting
    if (shiftdown && s6_feas_ge(newsolval, solval, feastol))
        || (!shiftdown && s6_feas_le(newsolval, solval, feastol))
    {
        return 0.0;
    }

    // C lines 145-148: get data of LP column (column-wise from a_matrix_)
    let a_start = &model.a_matrix_.start_;
    let a_index = &model.a_matrix_.index_;
    let a_value = &model.a_matrix_.value_;
    let col_start = a_start[col as usize] as usize;
    let col_end = a_start[(col + 1) as usize] as usize;
    let ncolrows = col_end - col_start;

    debug_assert!(col_end >= col_start);

    // C lines 153-205: find maximal shift value, st. all rows stay valid
    for i in 0..ncolrows {
        let row = a_index[col_start + i];
        let colval = a_value[col_start + i];
        let rowpos = row as usize;

        // C line 163: only global rows need to be valid (all model rows are global)
        // C line 168: left = shiftdown == (colvals[i] > 0)
        let left = shiftdown == (colval > 0.0);
        // C line 169: side = left ? lhs : rhs
        let side = if left {
            model.row_lower_[rowpos]
        } else {
            model.row_upper_[rowpos]
        };

        // C line 172: only finite bounds need to be considered
        if !s6_is_inf(if left { -side } else { side }) {
            // C line 177: newsolvalrow = solval + (side - activities[rowpos]) / colvals[i]
            let newsolvalrow = solval + (side - activities[rowpos]) / colval;

            // C line 180: update shifting value
            if (shiftdown && newsolvalrow > newsolval) || (!shiftdown && newsolvalrow < newsolval) {
                // C line 184: activity = activities[rowpos] + colvals[i] * ((shiftdown ? floor(newsolvalrow) : ceil(newsolvalrow)) - solval)
                let activity = activities[rowpos]
                    + colval
                        * ((if shiftdown {
                            newsolvalrow.floor()
                        } else {
                            newsolvalrow.ceil()
                        }) - solval);

                // C lines 187-188: ensure that shifting preserves feasibility
                if shiftdown
                    == ((left && s6_feas_ge(activity, side, feastol))
                        || (!left && s6_feas_le(activity, side, feastol)))
                {
                    newsolval = newsolvalrow.floor();
                } else {
                    newsolval = newsolvalrow.ceil();
                }

                // C lines 199-201: newsolval might have reached solval -> avoid numerical shifting
                if (shiftdown && s6_feas_ge(newsolval, solval, feastol))
                    || (!shiftdown && s6_feas_le(newsolval, solval, feastol))
                {
                    return 0.0;
                }
            }
        }
    }

    // C line 208: we must not shift variables to infinity
    if s6_is_inf(if shiftdown { -newsolval } else { newsolval }) {
        0.0
    } else {
        newsolval - solval
    }
}

/// SCIP 11.0.0 (32dd425): heur_oneopt.c::updateRowActivities (lines 210-259)
///
/// Update row activities after a variable's solution value changed.
pub(crate) fn s6_update_row_activities(
    ms: &HighsMipSolver,
    activities: &mut [f64],
    col: HighsInt,
    shiftval: f64,
) {
    let model = mip_model!(ms);
    let a_start = &model.a_matrix_.start_;
    let a_index = &model.a_matrix_.index_;
    let a_value = &model.a_matrix_.value_;
    let col_start = a_start[col as usize] as usize;
    let col_end = a_start[(col + 1) as usize] as usize;
    let ncolrows = col_end - col_start;

    debug_assert!(col_end >= col_start);

    // C lines 237-256: enumerate all rows with nonzero entry in this column
    for i in 0..ncolrows {
        let row = a_index[col_start + i];
        let colval = a_value[col_start + i];
        let rowpos = row as usize;

        // C line 247: only regard global rows in the LP (all model rows are global)
        // C line 249: activities[rowpos] += shiftval * colvals[i]
        activities[rowpos] += shiftval * colval;

        // C lines 251-254: clamp to ±infinity
        if s6_is_inf(activities[rowpos]) {
            activities[rowpos] = S6_INFINITY;
        } else if s6_is_inf(-activities[rowpos]) {
            activities[rowpos] = -S6_INFINITY;
        }
    }
}

/// SCIP 11.0.0 (32dd425): heur_oneopt.c::heurExecOneopt (lines 476-917)
///
/// The 1-opt heuristic: shifts integer variables toward their objective-optimal
/// bound while maintaining row feasibility.
pub(crate) fn s6_oneopt_exec(ms: &mut HighsMipSolver) -> S6Result {
    // C lines 500-507: get heuristic's data, check if we already processed this solution
    let s6 = match ms.mipdata_.s6.as_mut() {
        Some(s) => s,
        None => return S6Result::DidNotRun,
    };

    // C line 508: *result = SCIP_DELAYED
    // C lines 511-513: we only want to process each solution once
    let bestsol = match s6.pool.best() {
        Some(sol) => sol,
        None => return S6Result::Delayed,
    };
    let bestsol_index = bestsol.index;
    let bestsol_x = bestsol.x.clone();

    if s6.oneopt_lastsolindex == bestsol_index {
        return S6Result::Delayed;
    }

    // ncalls is counted when the exec gets past step 1
    s6.oneopt_ncalls += 1;

    // C lines 519-523: get problem variables
    let model = mip_model!(ms);
    let nvars = model.num_col_ as usize;
    let integrality = &model.integrality_;

    // C line 522: nenfovars = nvars - ncont - nimplint (integer columns only)
    let int_cols: Vec<HighsInt> = (0..model.num_col_)
        .filter(|&i| integrality[i as usize] == HighsVarType::kInteger)
        .collect();
    let nenfovars = int_cols.len();

    // C lines 526-530: do not run if there are no discrete variables
    if nenfovars == 0 {
        return S6Result::DidNotRun;
    }

    // C line 566: we can only work on solutions valid in the transformed space (always true here)
    // C lines 594-596: we need an LP (Q9: the LP rows are the model rows); *result is still SCIP_DELAYED
    if model.num_row_ == 0 {
        return S6Result::Delayed;
    }

    // C line 598: *result = SCIP_DIDNOTFIND
    // C line 600: heurdata->lastsolindex = SCIPsolGetIndex(bestsol)
    s6.oneopt_lastsolindex = bestsol_index;

    // C line 601: create working solution (copy of bestsol)
    let mut x = bestsol_x;

    // C line 604: starting bound adjustment
    // C lines 606-634: change solution values due to possible global bound changes
    // nchgbound counter (heur_oneopt.c:606,623,631) is debug-only (SCIPdebugMsg), omitted.
    let feastol = ms.mipdata_.feastol;
    let domains = &ms.mipdata_.domains[0];
    for i in (0..nvars).rev() {
        let lb = domains.col_lower_[i];
        let ub = domains.col_upper_[i];
        let solval = x[i];
        // C line 619: old solution value is smaller than the actual lower bound
        if s6_feas_lt(solval, lb, feastol) {
            x[i] = lb;
        }
        // C line 627: old solution value is greater than the actual upper bound
        else if s6_feas_gt(solval, ub, feastol) {
            x[i] = ub;
        }
    }

    // C lines 638-663: initialize LP row activities
    let nlprows = model.num_row_ as usize;
    let mut activities = vec![0.0f64; nlprows];
    let mut valid = true;

    // Compute activities: iterate all columns, for each column iterate its rows
    let a_start = &model.a_matrix_.start_;
    let a_index = &model.a_matrix_.index_;
    let a_value = &model.a_matrix_.value_;
    for j in 0..model.num_col_ as usize {
        let col_start = a_start[j] as usize;
        let col_end = a_start[j + 1] as usize;
        for k in col_start..col_end {
            let row = a_index[k] as usize;
            activities[row] += a_value[k] * x[j];
        }
    }

    // lp.c:6730-6732 (SCIProwGetSolActivity): activity = MAX(activity, -inf); activity = MIN(activity, +inf)
    for act in activities.iter_mut() {
        *act = cpp_max(*act, -S6_INFINITY);
        *act = cpp_min(*act, S6_INFINITY);
    }

    // C lines 655-661: check feasibility of activities
    for i in 0..nlprows {
        if s6_feas_lt(activities[i], model.row_lower_[i], feastol)
            || s6_feas_gt(activities[i], model.row_upper_[i], feastol)
        {
            valid = false;
            break;
        }
    }

    // C lines 665-674: if not valid, return
    if !valid {
        return S6Result::DidNotFind;
    }

    // Compute ncont before the shift loop so we don't hold `model` across mutable calls
    let ncont = (model.num_col_ as usize) - nenfovars;

    // C lines 677-681: allocate buffer storage for possible shift candidates
    let mut shiftcandssize = 8usize;
    let mut shiftcands: Vec<HighsInt> = vec![0; shiftcandssize];
    let mut shiftvals: Vec<f64> = vec![0.0; shiftcandssize];
    let mut nsuccessfulshifts = 0;
    // niterations counter (heur_oneopt.c:681,687) is debug-only (SCIPdebugMsg), omitted.

    // C lines 682-792: do { ... } while (useloop && shifted)
    let mut shifted;
    loop {
        // C lines 685-688: initialize data
        shifted = false;
        let mut nshiftcands = 0usize;

        // C lines 694-721: enumerate all integer variables and find out which of them are shiftable
        for &col in &int_cols {
            let solval = x[col as usize];
            let shiftval = s6_calc_shift_val(ms, col, solval, &activities);

            // C line 706: insert the variable into the list of shifting candidates
            if !s6_feas_zero(shiftval, feastol) {
                if nshiftcands == shiftcandssize {
                    shiftcandssize *= 8;
                    shiftcands.resize(shiftcandssize, 0);
                    shiftvals.resize(shiftcandssize, 0.0);
                }
                shiftcands[nshiftcands] = col;
                shiftvals[nshiftcands] = shiftval;
                nshiftcands += 1;
            }
        }

        // C lines 724-790: if at least one variable can be shifted, shift variables sorted by their objective
        if nshiftcands > 0 {
            if nshiftcands == 1 {
                // C lines 731-743: the case that exactly one variable can be shifted
                let var = shiftcands[0];
                let solval = x[var as usize];
                let shiftval = shiftvals[0];
                x[var as usize] = solval + shiftval;
                s6_update_row_activities(ms, &mut activities, var, shiftval);
                nsuccessfulshifts += 1;
            } else {
                // C lines 744-788: multiple candidates, sort by objective
                let mut objcoeffs: Vec<f64> = vec![0.0; nshiftcands];
                let mut ptrs: Vec<usize> = (0..nshiftcands).collect();

                // C lines 753-762: compute objective coefficients (weighted by shiftval)
                let model_inner = mip_model!(ms);
                for i in 0..nshiftcands {
                    objcoeffs[i] = model_inner.col_cost_[shiftcands[i] as usize] * shiftvals[i];
                }

                // C line 765: sort arrays with respect to the first one
                s6_sort_real_ptr(&mut objcoeffs, &mut ptrs);

                // C lines 768-785: try to shift each variable
                for i in 0..nshiftcands {
                    let var = shiftcands[ptrs[i]];
                    let solval = x[var as usize];
                    let shiftval = s6_calc_shift_val(ms, var, solval, &activities);

                    // C line 778: update data structures for nonzero shift value
                    if !s6_feas_zero(shiftval, feastol) {
                        x[var as usize] = solval + shiftval;
                        s6_update_row_activities(ms, &mut activities, var, shiftval);
                        nsuccessfulshifts += 1;
                    }
                }
            }
            shifted = true;
        }

        // C line 792: while (heurdata->useloop && shifted)
        if !(S6_ONEOPT_USELOOP && shifted) {
            break;
        }
    }

    // C lines 794-900: if nsuccessfulshifts > 0, try to install the solution
    let mut found = false;
    if nsuccessfulshifts > 0 {
        // C line 799: check if pure IP or need LP
        if ncont == 0 {
            // C lines 801-813: pure IP, try solution directly
            found = ms.try_solution(&x, kSolutionSourceHeuristic);
        } else {
            // C lines 815-899: mixed problem, use try_rounded_point (stands in for SCIP dive LP)
            found = detach!(ms.mipdata_.heuristics, |h| {
                h.try_rounded_point(ms, 0, &x, kSolutionSourceHeuristic)
            });
        }
    }

    // C lines 905-906: useloop: lastsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip))
    if S6_ONEOPT_USELOOP {
        if let Some(s6) = ms.mipdata_.s6.as_mut() {
            if let Some(sol) = s6.pool.best() {
                s6.oneopt_lastsolindex = sol.index;
            }
        }
    }

    // C lines 806-889: *result = SCIP_FOUNDSOL on success
    if found {
        S6Result::FoundSol
    } else {
        S6Result::DidNotFind
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;
    use crate::lp_data::hconst::kHighsInf;
    use crate::lp_data::highs_callback::HighsCallback;
    use crate::lp_data::highs_lp::HighsLp;
    use crate::lp_data::highs_options::HighsOptions;
    use crate::mip::HighsMipSolver;
    use crate::mip::Slot;
    use crate::mip::mip_worker::HighsMipWorker;
    use crate::mip::scip_heuristics::S6Origin;
    use crate::mip::scip_heuristics::S6Sol;
    use crate::mip::scip_heuristics::S6State;

    fn make_ms_with_lp(
        num_col: HighsInt,
        num_row: HighsInt,
        col_lower: Vec<f64>,
        col_upper: Vec<f64>,
        row_lower: Vec<f64>,
        row_upper: Vec<f64>,
        col_cost: Vec<f64>,
        integrality: Vec<HighsVarType>,
        a_start: Vec<HighsInt>,
        a_index: Vec<HighsInt>,
        a_value: Vec<f64>,
    ) -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = num_col;
        lp.num_row_ = num_row;
        lp.col_cost_ = col_cost;
        lp.col_lower_ = col_lower;
        lp.col_upper_ = col_upper;
        lp.row_lower_ = row_lower;
        lp.row_upper_ = row_upper;
        lp.integrality_ = integrality;
        lp.a_matrix_.num_col_ = num_col;
        lp.a_matrix_.num_row_ = num_row;
        lp.a_matrix_.start_ = a_start;
        lp.a_matrix_.index_ = a_index;
        lp.a_matrix_.value_ = a_value;
        let options = HighsOptions {
            output_flag: false,
            ..Default::default()
        };
        let mut ms = HighsMipSolver::new(
            &mut HighsCallback::default(),
            &options,
            &lp,
            &HighsSolution::default(),
            false,
            0,
        );
        // Build the solver data (lps, domains[0], heuristics) as run_mip_presolve does (new_data, init),
        // without presolve, so that model_ == Orig and domains[0] holds the global bounds of `lp`.
        ms.new_data();
        ms.init();
        // Row-wise copy of the matrix (C++ ARstart_ / ARindex_ / ARvalue_), used by try_solution.
        let nrow = num_row as usize;
        let mut rowcount = vec![0usize; nrow];
        for &r in lp.a_matrix_.index_.iter() {
            rowcount[r as usize] += 1;
        }
        let mut ar_start: Vec<HighsInt> = vec![0; nrow + 1];
        for i in 0..nrow {
            ar_start[i + 1] = ar_start[i] + rowcount[i] as HighsInt;
        }
        let mut fill = ar_start.clone();
        let nnz = lp.a_matrix_.index_.len();
        let mut ar_index: Vec<HighsInt> = vec![0; nnz];
        let mut ar_value: Vec<f64> = vec![0.0; nnz];
        for j in 0..num_col as usize {
            for k in lp.a_matrix_.start_[j] as usize..lp.a_matrix_.start_[j + 1] as usize {
                let r = lp.a_matrix_.index_[k] as usize;
                let pos = fill[r] as usize;
                ar_index[pos] = j as HighsInt;
                ar_value[pos] = lp.a_matrix_.value_[k];
                fill[r] += 1;
            }
        }
        ms.mipdata_.a_rstart_ = ar_start;
        ms.mipdata_.a_rindex_ = ar_index;
        ms.mipdata_.a_rvalue_ = ar_value;
        ms
    }

    #[test]
    fn s6_oneopt_calc_shift_val_obj_positive() {
        // 1 int var, 1 row: 0 <= x <= 5, obj > 0 -> shift down toward lb
        // x in [0, 10], solval = 3
        // obj > 0: newsolval = ceil(lb) = ceil(0) = 0
        // shiftdown = true, newsolval(0) < solval(3) -> proceed
        // Row: left = shiftdown == (colval>0) = true == true = true
        // side = row_lower = 0, not inf
        // newsolvalrow = 3 + (0 - 3) / 1 = 0
        // shiftdown && 0 > 0? No. No update.
        // Return 0 - 3 = -3
        let ms = make_ms_with_lp(
            1,
            1,
            vec![0.0],
            vec![10.0],
            vec![0.0],
            vec![5.0],
            vec![1.0],
            vec![HighsVarType::kInteger],
            vec![0, 1],
            vec![0],
            vec![1.0],
        );
        let activities = vec![3.0];
        let shiftval = s6_calc_shift_val(&ms, 0, 3.0, &activities);
        assert!((shiftval - (-3.0)).abs() < 1e-9);
    }

    #[test]
    fn s6_oneopt_calc_shift_val_obj_negative() {
        // 1 int var, 1 row: x >= 2 (row_lower=2, row_upper=inf), obj < 0 -> shift up toward ub
        // x in [0, 10], solval = 5
        // obj < 0: newsolval = floor(ub) = floor(10) = 10, shiftdown=false
        // 10 > 5 so proceed. Row: left = false == (1>0) = false, side = row_upper = inf -> skip
        // No row limits. Return 10 - 5 = 5.
        let ms = make_ms_with_lp(
            1,
            1,
            vec![0.0],
            vec![10.0],
            vec![2.0],
            vec![kHighsInf],
            vec![-1.0],
            vec![HighsVarType::kInteger],
            vec![0, 1],
            vec![0],
            vec![1.0],
        );
        let activities = vec![5.0];
        let shiftval = s6_calc_shift_val(&ms, 0, 5.0, &activities);
        assert!((shiftval - 5.0).abs() < 1e-9);
    }

    #[test]
    fn s6_oneopt_calc_shift_val_obj_zero() {
        // obj == 0: newsolval = solval, then feas_ge/le check returns 0
        let ms = make_ms_with_lp(
            1,
            1,
            vec![0.0],
            vec![10.0],
            vec![0.0],
            vec![10.0],
            vec![0.0],
            vec![HighsVarType::kInteger],
            vec![0, 1],
            vec![0],
            vec![1.0],
        );
        let activities = vec![3.0];
        let shiftval = s6_calc_shift_val(&ms, 0, 3.0, &activities);
        assert!((shiftval - 0.0).abs() < 1e-9);
    }

    #[test]
    fn s6_oneopt_calc_shift_val_infinite_side_skipped() {
        // Row with both sides infinite: x free. Should skip the row.
        // obj > 0, shiftdown=true, newsolval=ceil(lb)=ceil(0)=0
        // Row: left=true, side=row_lower=-inf -> s6_is_inf(-(-inf))=s6_is_inf(inf)=true -> skip
        // No row constraint. Return 0 - 3 = -3.
        let ms = make_ms_with_lp(
            1,
            1,
            vec![0.0],
            vec![10.0],
            vec![-kHighsInf],
            vec![kHighsInf],
            vec![1.0],
            vec![HighsVarType::kInteger],
            vec![0, 1],
            vec![0],
            vec![1.0],
        );
        let activities = vec![3.0];
        let shiftval = s6_calc_shift_val(&ms, 0, 3.0, &activities);
        assert!((shiftval - (-3.0)).abs() < 1e-9);
    }

    #[test]
    fn s6_oneopt_calc_shift_val_feasibility_returns_zero() {
        // When the floor/ceil feasibility test fails and newsolval reaches solval, return 0.
        // x in [0,10], solval=3, obj>0 -> newsolval=ceil(0)=0
        // Row: row_lower=3, row_upper=inf, colval=1
        // left=true, side=3, not inf. newsolvalrow = 3 + (3-3)/1 = 3.
        // shiftdown && 3 > 0? Yes. activity = 3 + 1*(floor(3)-3) = 3+0=3.
        // left && feas_ge(3,3)? Yes. shiftdown==true==true. newsolval=floor(3)=3.
        // Then check: shiftdown && feas_ge(3,3)? Yes -> return 0.
        let ms = make_ms_with_lp(
            1,
            1,
            vec![0.0],
            vec![10.0],
            vec![3.0],
            vec![kHighsInf],
            vec![1.0],
            vec![HighsVarType::kInteger],
            vec![0, 1],
            vec![0],
            vec![1.0],
        );
        let activities = vec![3.0];
        let shiftval = s6_calc_shift_val(&ms, 0, 3.0, &activities);
        assert!((shiftval - 0.0).abs() < 1e-9);
    }

    #[test]
    fn s6_oneopt_update_row_activities_clamp() {
        // Test that activities are clamped at ±S6_INFINITY
        let ms = make_ms_with_lp(
            1,
            1,
            vec![0.0],
            vec![10.0],
            vec![0.0],
            vec![10.0],
            vec![1.0],
            vec![HighsVarType::kInteger],
            vec![0, 1],
            vec![0],
            vec![1.0],
        );
        let mut activities = vec![S6_INFINITY - 1.0];
        s6_update_row_activities(&ms, &mut activities, 0, 2.0);
        // activities[0] = (1e20 - 1) + 2*1 = 1e20 + 1, which is inf -> clamp to S6_INFINITY
        assert!((activities[0] - S6_INFINITY).abs() < 1e-9);

        let mut activities2 = vec![-(S6_INFINITY - 1.0)];
        s6_update_row_activities(&ms, &mut activities2, 0, -2.0);
        // activities2[0] = -(1e20-1) + (-2)*1 = -1e20 - 1, which is -inf -> clamp to -S6_INFINITY
        assert!((activities2[0] - (-S6_INFINITY)).abs() < 1e-9);
    }

    #[test]
    fn s6_oneopt_exec_early_return_lastsolindex() {
        // exec returns early when lastsolindex equals the best index
        let mut ms = make_ms_with_lp(
            2,
            1,
            vec![0.0, 0.0],
            vec![10.0, 10.0],
            vec![0.0],
            vec![10.0],
            vec![1.0, 1.0],
            vec![HighsVarType::kInteger, HighsVarType::kInteger],
            vec![0, 1, 2],
            vec![0, 0],
            vec![1.0, 1.0],
        );
        ms.mipdata_.s6 = Some(Box::new(S6State::default()));
        let s6 = ms.mipdata_.s6.as_mut().unwrap();
        s6.pool.sols.push(S6Sol {
            index: 42,
            obj: 5.0,
            x: vec![3.0, 2.0],
            origin: S6Origin::Highs(0),
            nodenum: 0,
        });
        s6.oneopt_lastsolindex = 42;
        let result = s6_oneopt_exec(&mut ms);
        assert_eq!(result, S6Result::Delayed);
        assert_eq!(ms.mipdata_.s6.as_ref().unwrap().oneopt_ncalls, 0);
    }

    #[test]
    fn s6_oneopt_exec_single_candidate() {
        // 2 int vars, 1 row: x1 + x2 <= 10, obj = [1, 0]
        // sol = [5, 5], obj[0]>0 -> shift x1 down. obj[1]==0 -> no shift for x2.
        // x1: lb=0, ub=10, solval=5, obj=1>0 -> newsolval=ceil(0)=0, shiftdown=true
        // Row: left=true, side=row_lower=0, not inf.
        // newsolvalrow = 5 + (0 - 10)/1 = -5. shiftdown && -5 > 0? No.
        // So no row update. Return 0 - 5 = -5.
        // x2: obj=0 -> newsolval=solval=5, feas_ge(5,5) -> return 0.
        // Only x1 is a candidate with shiftval=-5.
        // After shift: x1 = 5 + (-5) = 0. Activity = 0+5=5 <= 10. OK.
        // try_solution called with [0, 5].
        let mut ms = make_ms_with_lp(
            2,
            1,
            vec![0.0, 0.0],
            vec![10.0, 10.0],
            vec![0.0],
            vec![10.0],
            vec![1.0, 0.0],
            vec![HighsVarType::kInteger, HighsVarType::kInteger],
            vec![0, 1, 2],
            vec![0, 0],
            vec![1.0, 1.0],
        );
        ms.mipdata_.s6 = Some(Box::new(S6State::default()));
        let s6 = ms.mipdata_.s6.as_mut().unwrap();
        s6.pool.sols.push(S6Sol {
            index: 1,
            obj: 5.0,
            x: vec![5.0, 5.0],
            origin: S6Origin::Highs(0),
            nodenum: 0,
        });
        s6.oneopt_lastsolindex = -1;
        let result = s6_oneopt_exec(&mut ms);
        assert_eq!(result, S6Result::FoundSol);
    }

    #[test]
    fn s6_oneopt_exec_multi_candidate_sort() {
        // 2 int vars, 1 row: x1 + x2 <= 10, obj = [2, 1]
        // sol = [5, 5]. Both have obj > 0, both can shift down.
        // x1: shiftval = 0 - 5 = -5 (limited by lb=0)
        // x2: shiftval = 0 - 5 = -5 (limited by lb=0)
        // weightedobj: objcoeffs = [2*(-5), 1*(-5)] = [-10, -5]
        // sort ascending: [-10, -5] -> ptrs = [0, 1]
        // Process x1 first: recompute shiftval with current activities.
        // After x1 shift: x1=0, activity=0+5=5.
        // x2: solval=5, obj>0, newsolval=ceil(0)=0.
        // Row: left=true, side=0. newsolvalrow = 5 + (0-5)/1 = 0.
        // shiftdown && 0 > 0? No. Return 0-5=-5.
        // x2 shifts to 0. Activity = 0+0=0.
        // Both shifted. try_solution with [0, 0].
        let mut ms = make_ms_with_lp(
            2,
            1,
            vec![0.0, 0.0],
            vec![10.0, 10.0],
            vec![0.0],
            vec![10.0],
            vec![2.0, 1.0],
            vec![HighsVarType::kInteger, HighsVarType::kInteger],
            vec![0, 1, 2],
            vec![0, 0],
            vec![1.0, 1.0],
        );
        ms.mipdata_.s6 = Some(Box::new(S6State::default()));
        let s6 = ms.mipdata_.s6.as_mut().unwrap();
        s6.pool.sols.push(S6Sol {
            index: 1,
            obj: 15.0,
            x: vec![5.0, 5.0],
            origin: S6Origin::Highs(0),
            nodenum: 0,
        });
        s6.oneopt_lastsolindex = -1;
        let result = s6_oneopt_exec(&mut ms);
        assert_eq!(result, S6Result::FoundSol);
    }

    #[test]
    fn s6_oneopt_exec_mixed_uses_try_rounded_point() {
        // 1 int + 1 cont, 1 row: x_int + x_cont <= 10
        // obj = [1, 0]. sol = [5, 5].
        // x_int: obj>0, shift down. x_cont: not integer, not a candidate.
        // After shift: x_int=0. Mixed problem -> try_rounded_point.
        let mut ms = make_ms_with_lp(
            2,
            1,
            vec![0.0, 0.0],
            vec![10.0, 10.0],
            vec![0.0],
            vec![10.0],
            vec![1.0, 0.0],
            vec![HighsVarType::kInteger, HighsVarType::kContinuous],
            vec![0, 1, 2],
            vec![0, 0],
            vec![1.0, 1.0],
        );
        ms.mipdata_.s6 = Some(Box::new(S6State::default()));
        let s6 = ms.mipdata_.s6.as_mut().unwrap();
        s6.pool.sols.push(S6Sol {
            index: 1,
            obj: 5.0,
            x: vec![5.0, 5.0],
            origin: S6Origin::Highs(0),
            nodenum: 0,
        });
        s6.oneopt_lastsolindex = -1;
        // Master worker, as run() creates it (mip_solver.rs:343), needed by try_rounded_point.
        let master = HighsMipWorker::new(&mut ms, 0, 0, 0, 0, 0);
        ms.mipdata_.workers.push(Slot::new(master));
        let result = s6_oneopt_exec(&mut ms);
        assert_eq!(result, S6Result::FoundSol);
    }
}
