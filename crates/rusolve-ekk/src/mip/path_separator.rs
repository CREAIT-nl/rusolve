//! HiGHS source: mip/HighsPathSeparator.h, mip/HighsPathSeparator.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use std::collections::HashMap;

use crate::prelude::*;

/// HiGHS: HighsPathSeparator.h::HighsPathSeparator (lines 23-37)
#[derive(Clone, Default, Debug)]
pub struct HighsPathSeparator {
    /// C++ base class `HighsSeparator`.
    pub base: HighsSeparator,
    /// C++: `HighsRandom randgen`
    pub randgen: HighsRandom,
}

impl HighsPathSeparator {
    /// HiGHS: HighsPathSeparator.cpp::HighsPathSeparator::separateLpSolution (lines 26-550)
    /// C++: `void separateLpSolution(HighsLpRelaxation & lpRelaxation, HighsLpAggregator & lpAggregator, HighsTransformedLp & transLp, HighsCutPool & cutpool)`
    // PORT NOTE: `ms` is `lpRelaxation.getMipSolver()`; `cutpool` is an index into
    // `ms.mipdata_.cutpools` (see `HighsSeparatorBox` in separator.rs).
    // HiGHS quirk: `success` is overwritten by the first `generateCut` before any read (the path
    // loop always runs once), and the final `success |=` of the path-mixing block is never read in
    // C++ either, so `unused_assignments` is allowed for the whole body.
    #[allow(unused_assignments)]
    pub fn separate_lp_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &mut HighsLpRelaxation,
        lp_aggregator: &mut HighsLpAggregator,
        trans_lp: &mut HighsTransformedLp,
        cutpool: usize,
    ) {
        let feastol = ms.mipdata_.feastol;
        let epsilon = ms.mipdata_.epsilon;
        let lp = lp_relaxation.get_lp();
        let lp_solution = lp_relaxation.get_solution();

        let num_row = lp.num_row_ as usize;
        let num_col = lp.num_col_ as usize;

        let mut rowtype: Vec<PathSeparatorRowType> = vec![PathSeparatorRowType::kEq; num_row];
        for i in 0..lp.num_row_ {
            let iu = i as usize;
            if lp.row_lower_[iu] == lp.row_upper_[iu] {
                rowtype[iu] = PathSeparatorRowType::kEq;
                continue;
            }

            let mut lowerslack = kHighsInf;
            let mut upperslack = kHighsInf;

            if lp.row_lower_[iu] != -kHighsInf {
                lowerslack = lp_solution.row_value[iu] - lp.row_lower_[iu];
            }

            if lp.row_upper_[iu] != kHighsInf {
                upperslack = lp.row_upper_[iu] - lp_solution.row_value[iu];
            }

            rowtype[iu] = if lowerslack > feastol && upperslack > feastol {
                PathSeparatorRowType::kUnusuable
            } else if lowerslack < upperslack {
                PathSeparatorRowType::kGeq
            } else {
                PathSeparatorRowType::kLeq
            };
        }

        let mut num_continuous: Vec<HighsInt> = vec![0; num_row];

        let mut max_aggr_row_size: usize = 0;
        for &col in ms.mipdata_.continuous_cols.iter() {
            if trans_lp.bound_distance(col) == 0.0 {
                continue;
            }

            let cu = col as usize;
            max_aggr_row_size += (lp.a_matrix_.start_[cu + 1] - lp.a_matrix_.start_[cu]) as usize;
            for i in lp.a_matrix_.start_[cu]..lp.a_matrix_.start_[cu + 1] {
                num_continuous[lp.a_matrix_.index_[i as usize] as usize] += 1;
            }
        }

        let mut col_substitutions: Vec<(HighsInt, f64)> = vec![(-1, 0.0); num_col];

        // identify equality rows where only a single continuous variable with nonzero
        // transformed solution value is present. Mark those columns and remember the
        // rows so that we can always substitute such columns away using this equation
        // and block the equation from being used as a start row
        for i in 0..lp.num_row_ {
            let iu = i as usize;
            if rowtype[iu] != PathSeparatorRowType::kEq || num_continuous[iu] != 1 {
                continue;
            }

            let mut len: HighsInt = 0;
            let mut rowinds: &[HighsInt] = &[];
            let mut rowvals: &[f64] = &[];
            lp_relaxation.get_row(ms, i, &mut len, &mut rowinds, &mut rowvals);

            // find continuous variable
            let mut col: HighsInt = -1;
            let mut val = 0.0;
            for j in 0..len as usize {
                if ms.is_col_integral(rowinds[j]) {
                    continue;
                }
                if trans_lp.bound_distance(rowinds[j]) == 0.0 {
                    continue;
                }
                col = rowinds[j];
                val = rowvals[j];
                break;
            }

            debug_assert!(col != -1);
            debug_assert!(ms.is_col_continuous(col));
            debug_assert!(trans_lp.bound_distance(col) > 0.0);

            if col_substitutions[col as usize].0 != -1 {
                continue;
            }

            col_substitutions[col as usize] = (i, val);
            rowtype[iu] = PathSeparatorRowType::kUnusuable;
        }

        // for each continuous variable with nonzero transformed solution value
        // remember the <= and == rows where it is present with a positive coefficient
        // in its set of in-arc rows. Treat >= rows as <= rows with reversed sign
        // The reason to only store one set of rows for one sign of the coefficients
        // is that this directs the selection to be more diverse. Consider
        // aggregations of 2 rows where we allow both directions. When one of the rows
        // is used as start row we can always select the other one. When we only
        // project out variables with negative coefficients we give the aggregation
        // path an orientation and avoid this situation. For each aggregation of rows
        // the cut generation will try the reversed orientation in any case too.

        let mut in_arc_rows: Vec<(HighsInt, f64)> = Vec::with_capacity(max_aggr_row_size);
        let mut col_in_arcs: Vec<(HighsInt, HighsInt)> = vec![(0, 0); num_col];

        let mut out_arc_rows: Vec<(HighsInt, f64)> = Vec::with_capacity(max_aggr_row_size);
        let mut col_out_arcs: Vec<(HighsInt, HighsInt)> = vec![(0, 0); num_col];

        for &col in ms.mipdata_.continuous_cols.iter() {
            if trans_lp.bound_distance(col) == 0.0 {
                continue;
            }
            let cu = col as usize;
            if col_substitutions[cu].0 != -1 {
                continue;
            }

            col_in_arcs[cu].0 = in_arc_rows.len() as HighsInt;
            col_out_arcs[cu].0 = out_arc_rows.len() as HighsInt;
            for i in lp.a_matrix_.start_[cu]..lp.a_matrix_.start_[cu + 1] {
                let iu = i as usize;
                let row = lp.a_matrix_.index_[iu];
                let value = lp.a_matrix_.value_[iu];
                match rowtype[row as usize] {
                    PathSeparatorRowType::kUnusuable => continue,
                    PathSeparatorRowType::kLeq => {
                        if value < 0.0 {
                            in_arc_rows.push((row, value));
                        } else {
                            out_arc_rows.push((row, value));
                        }
                    }
                    PathSeparatorRowType::kGeq => {
                        if value > 0.0 {
                            in_arc_rows.push((row, value));
                        } else {
                            out_arc_rows.push((row, value));
                        }
                    }
                    PathSeparatorRowType::kEq => {
                        in_arc_rows.push((row, value));
                        out_arc_rows.push((row, value));
                    }
                }
            }

            col_in_arcs[cu].1 = in_arc_rows.len() as HighsInt;
            col_out_arcs[cu].1 = out_arc_rows.len() as HighsInt;
        }

        let mut cut_gen = HighsCutGeneration::new(ms, lp_relaxation, cutpool);
        let mut base_row_inds: Vec<HighsInt> = Vec::new();
        let mut base_row_vals: Vec<f64> = Vec::new();
        const MAX_PATH_LEN: usize = 6;
        let mut current_path: [HighsInt; MAX_PATH_LEN] = [0; MAX_PATH_LEN];
        let mut aggregated_path: Vec<(Vec<HighsInt>, Vec<f64>)> = Vec::new();
        let mut scales: [f64; 2] = [0.0; 2];
        let max_weight = 1.0 / feastol;
        let min_weight = feastol;

        for i in 0..lp.num_row_ {
            let iu = i as usize;
            match rowtype[iu] {
                PathSeparatorRowType::kUnusuable => continue,
                PathSeparatorRowType::kEq => {
                    if lp_solution.row_dual[iu] <= epsilon {
                        scales[0] = 1.0;
                        scales[1] = -1.0;
                    } else {
                        scales[0] = -1.0;
                        scales[1] = 1.0;
                    }
                }
                PathSeparatorRowType::kLeq => {
                    scales[0] = 1.0;
                    scales[1] = -1.0;
                }
                PathSeparatorRowType::kGeq => {
                    scales[0] = -1.0;
                    scales[1] = 1.0;
                }
            }

            let mut success = false;

            for s in 0..2usize {
                debug_assert!(lp_aggregator.is_empty());
                lp_aggregator.add_row(ms, lp_relaxation, i, scales[s]);

                current_path[0] = i;
                let mut curr_path_len: HighsInt = 1;

                let mut try_negated_scale = false;

                aggregated_path.clear();

                while curr_path_len != MAX_PATH_LEN as HighsInt {
                    lp_aggregator.get_current_aggregation(
                        ms,
                        lp_relaxation,
                        &mut base_row_inds,
                        &mut base_row_vals,
                        false,
                    );
                    let base_row_len = base_row_inds.len() as HighsInt;
                    let mut added_substitution_rows = false;

                    let mut best_out_arc_col: HighsInt = -1;
                    let mut out_arc_col_val = 0.0;
                    let mut out_arc_col_bound_dist = 0.0;

                    let mut best_in_arc_col: HighsInt = -1;
                    let mut in_arc_col_val = 0.0;
                    let mut in_arc_col_bound_dist = 0.0;

                    for j in 0..base_row_len as usize {
                        let col = base_row_inds[j];
                        if col >= lp.num_col_
                            || trans_lp.bound_distance(col) == 0.0
                            || lp_relaxation.is_col_integral(ms, col)
                        {
                            continue;
                        }

                        if col_substitutions[col as usize].0 != -1 {
                            added_substitution_rows = true;
                            let (sub_row, sub_val) = col_substitutions[col as usize];
                            lp_aggregator.add_row(
                                ms,
                                lp_relaxation,
                                sub_row,
                                -base_row_vals[j] / sub_val,
                            );
                            continue;
                        }

                        if added_substitution_rows {
                            continue;
                        }

                        if base_row_vals[j] < 0.0 {
                            if skip_col(
                                i,
                                &current_path[..curr_path_len as usize],
                                &mut try_negated_scale,
                                col,
                                &col_out_arcs,
                                &out_arc_rows,
                                &col_in_arcs,
                                &in_arc_rows,
                            ) {
                                continue;
                            }

                            if best_out_arc_col == -1
                                || trans_lp.bound_distance(col) > out_arc_col_bound_dist
                            {
                                best_out_arc_col = col;
                                out_arc_col_val = base_row_vals[j];
                                out_arc_col_bound_dist = trans_lp.bound_distance(col);
                            }
                        } else {
                            if skip_col(
                                i,
                                &current_path[..curr_path_len as usize],
                                &mut try_negated_scale,
                                col,
                                &col_in_arcs,
                                &in_arc_rows,
                                &col_out_arcs,
                                &out_arc_rows,
                            ) {
                                continue;
                            }

                            if best_in_arc_col == -1
                                || trans_lp.bound_distance(col) > in_arc_col_bound_dist
                            {
                                best_in_arc_col = col;
                                in_arc_col_val = base_row_vals[j];
                                in_arc_col_bound_dist = trans_lp.bound_distance(col);
                            }
                        }
                    }

                    if added_substitution_rows {
                        continue;
                    }

                    // generate cut
                    let mut rhs: f64 = 0.0;
                    success = cut_gen.generate_cut(
                        ms,
                        lp_relaxation,
                        trans_lp,
                        &mut base_row_inds,
                        &mut base_row_vals,
                        &mut rhs,
                        false,
                    );

                    lp_aggregator.get_current_aggregation(
                        ms,
                        lp_relaxation,
                        &mut base_row_inds,
                        &mut base_row_vals,
                        true,
                    );
                    if !aggregated_path.is_empty()
                        || best_out_arc_col != -1
                        || best_in_arc_col != -1
                    {
                        aggregated_path.push((base_row_inds.clone(), base_row_vals.clone()));
                    }

                    // generate reverse cut
                    rhs = 0.0;
                    success |= cut_gen.generate_cut(
                        ms,
                        lp_relaxation,
                        trans_lp,
                        &mut base_row_inds,
                        &mut base_row_vals,
                        &mut rhs,
                        false,
                    );

                    if success || (best_out_arc_col == -1 && best_in_arc_col == -1) {
                        break;
                    }

                    // we prefer to use an out edge if the bound distances are equal in
                    // feasibility tolerance otherwise we choose an inArc. This tie
                    // breaking is arbitrary, but we should direct the substitution to
                    // prefer one direction to increase diversity.
                    let found = if best_in_arc_col == -1
                        || (best_out_arc_col != -1
                            && out_arc_col_bound_dist >= in_arc_col_bound_dist - feastol)
                    {
                        match find_row(
                            best_out_arc_col,
                            out_arc_col_val,
                            &col_in_arcs,
                            &in_arc_rows,
                            &current_path[..curr_path_len as usize],
                            &mut self.randgen,
                            max_weight,
                            min_weight,
                        ) {
                            Some(found) => Some(found),
                            None => {
                                if best_in_arc_col == -1 {
                                    None
                                } else {
                                    find_row(
                                        best_in_arc_col,
                                        in_arc_col_val,
                                        &col_out_arcs,
                                        &out_arc_rows,
                                        &current_path[..curr_path_len as usize],
                                        &mut self.randgen,
                                        max_weight,
                                        min_weight,
                                    )
                                }
                            }
                        }
                    } else {
                        find_row(
                            best_in_arc_col,
                            in_arc_col_val,
                            &col_out_arcs,
                            &out_arc_rows,
                            &current_path[..curr_path_len as usize],
                            &mut self.randgen,
                            max_weight,
                            min_weight,
                        )
                    };
                    let Some((row, weight)) = found else {
                        break;
                    };

                    current_path[curr_path_len as usize] = row;
                    lp_aggregator.add_row(ms, lp_relaxation, row, weight);
                    curr_path_len += 1;
                }

                // if the path has length at least 2 try to separate a path mixing cut
                let mut path_len = aggregated_path.len() as HighsInt;
                if path_len > 1 {
                    // generate path mixing cut
                    // PORT NOTE: `HighsHashTable<HighsInt, HighsInt> indexPos` is a plain
                    // lookup table here (only `operator[]` and reads are used), so a
                    // `HashMap` gives the same results.
                    let mut index_pos: HashMap<HighsInt, HighsInt> = HashMap::new();

                    let mut inds: Vec<HighsInt> = Vec::with_capacity(num_col + num_row);
                    let mut solval: Vec<f64> = Vec::with_capacity(num_col + num_row);
                    let mut upper: Vec<f64> = Vec::with_capacity(num_col + num_row);
                    let mut is_integral: Vec<bool> = Vec::with_capacity(num_col + num_row);

                    let mut rhs: Vec<f64> = vec![0.0; path_len as usize];
                    let mut tmp_upper: Vec<f64> = Vec::new();
                    let mut tmp_solval: Vec<f64> = Vec::new();

                    let mut delta: f64 = 1.0;

                    // HiGHS #3179 (mip_path_same_bds): HighsPathSeparator.cpp lines 415-422 at 28339ce4b
                    let path_same_bds = ms.options_mip_.mip_path_same_bds;
                    if path_same_bds {
                        trans_lp.init_multi_row_transform();
                    }

                    for k in 0..path_len as usize {
                        let mut integral_positive = false;

                        // PORT NOTE: C++ `transLp.transform` takes the path entry's vectors by
                        // non-const reference and may rewrite them, and the loops below read
                        // the rewritten data, so the stored entries are passed directly.
                        let transformed = {
                            let (path_inds, path_vals) = &mut aggregated_path[k];
                            trans_lp.transform(
                                ms,
                                lp_relaxation,
                                path_vals,
                                &mut tmp_upper,
                                &mut tmp_solval,
                                path_inds,
                                &mut rhs[k],
                                &mut integral_positive,
                                false,
                                path_same_bds,
                            )
                        };
                        if !transformed {
                            path_len = k as HighsInt;
                            break;
                        }

                        if k == 0 {
                            if rhs[k] > kHighsTiny {
                                path_len = k as HighsInt;
                                break;
                            }
                            if rhs[k] >= -feastol {
                                rhs[k] = 0.0;
                            }
                        } else if rhs[k] >= rhs[k - 1] - feastol {
                            path_len = k as HighsInt;
                            break;
                        }

                        delta = cpp_max(rhs[k].abs(), delta);

                        let len = aggregated_path[k].0.len();
                        for j in 0..len {
                            let index = aggregated_path[k].0[j];
                            let pos = *index_pos.entry(index).or_insert(0);
                            if pos == 0 {
                                inds.push(index);
                                solval.push(tmp_solval[j]);
                                upper.push(tmp_upper[j]);
                                is_integral.push(lp_relaxation.is_col_integral(ms, index));
                                if *is_integral.last().unwrap() {
                                    delta = cpp_max(aggregated_path[k].1[j].abs(), delta);
                                }
                                index_pos.insert(index, inds.len() as HighsInt);
                            } else {
                                let pu = (pos - 1) as usize;
                                debug_assert!(inds[pu] == index);
                                debug_assert!(solval[pu] == tmp_solval[j]);
                                debug_assert!(upper[pu] == tmp_upper[j]);
                                if is_integral[pu] {
                                    delta = cpp_max(aggregated_path[k].1[j].abs(), delta);
                                }
                            }
                        }
                    }

                    if path_len > 1 {
                        delta = (delta + 1.0).log2().ceil().exp2();

                        let mut num_inds = inds.len() as HighsInt;

                        let mut cut_rhs = HighsCDouble::from(0.0);
                        let mut cut_vals: Vec<f64> = vec![0.0; num_inds as usize];
                        let mut max_frac: Vec<f64> = vec![0.0; num_inds as usize];
                        let mut down_sum: Vec<HighsCDouble> =
                            vec![HighsCDouble::from(0.0); num_inds as usize];
                        let mut f_sum: Vec<HighsCDouble> =
                            vec![HighsCDouble::from(0.0); num_inds as usize];

                        let mut f_last: f64 = 0.0;
                        let mut scale: f64 = -1.0 / delta;

                        for k in 0..path_len as usize {
                            let f = rhs[k] * scale;
                            let f_diff = HighsCDouble::from(f) - f_last;
                            let len = aggregated_path[k].0.len();
                            cut_rhs += f_diff;
                            for j in 0..len {
                                let pi = (index_pos[&aggregated_path[k].0[j]] - 1) as usize;

                                let gj = aggregated_path[k].1[j] * scale;

                                if !is_integral[pi] {
                                    cut_vals[pi] = cpp_max(cut_vals[pi], gj);
                                } else {
                                    let gjdown = gj.floor();
                                    let hj = gj - gjdown;
                                    max_frac[pi] = cpp_max(max_frac[pi], hj);
                                    down_sum[pi] += f_diff * gjdown;
                                    f_sum[pi] += f_diff;

                                    if f_sum[pi] < max_frac[pi] {
                                        cut_vals[pi] = f64::from(down_sum[pi] + f_sum[pi]);
                                    } else {
                                        cut_vals[pi] = f64::from(down_sum[pi] + max_frac[pi]);
                                    }
                                }
                            }

                            if k > 0 {
                                let mut viol = f64::from(cut_rhs);
                                for j in 0..num_inds as usize {
                                    viol -= solval[j] * cut_vals[j];
                                }

                                viol *= delta;

                                if viol > 10.0 * feastol {
                                    scale = -delta;
                                    let mut rhs_cut = f64::from(cut_rhs) * scale;
                                    for cv in cut_vals.iter_mut().take(num_inds as usize) {
                                        *cv *= scale;
                                    }

                                    let mut j = num_inds - 1;
                                    while j >= 0 {
                                        let ju = j as usize;
                                        if cut_vals[ju].abs() <= epsilon {
                                            num_inds -= 1;
                                            let nu = num_inds as usize;
                                            cut_vals.swap(ju, nu);
                                            inds.swap(ju, nu);
                                        }
                                        j -= 1;
                                    }

                                    cut_vals.truncate(num_inds as usize);
                                    inds.truncate(num_inds as usize);

                                    if trans_lp.untransform(
                                        ms,
                                        lp_relaxation,
                                        &mut cut_vals,
                                        &mut inds,
                                        &mut rhs_cut,
                                        false,
                                    ) {
                                        success |= cut_gen.finalize_and_add_cut(
                                            ms,
                                            lp_relaxation,
                                            trans_lp.globaldom_,
                                            &mut inds,
                                            &mut cut_vals,
                                            &mut rhs_cut,
                                        );
                                    }

                                    // printf("cut is violated for k = %d\n", k);
                                    break;
                                }
                            }
                            f_last = f;
                        }
                    }
                }

                lp_aggregator.clear();

                if !try_negated_scale {
                    break;
                }
            }
        }
    }

    /// HiGHS: HighsPathSeparator.h::HighsPathSeparator::HighsPathSeparator (lines 33-36)
    /// C++: `HighsPathSeparator(const HighsMipSolver & mipsolver)`
    pub fn new(ms: &HighsMipSolver) -> Self {
        let mut randgen = HighsRandom::default();
        randgen.initialise(ms.options_mip_.random_seed as HighsUInt);
        HighsPathSeparator {
            base: HighsSeparator::new(ms, kPathAggrSepaString),
            randgen,
        }
    }
}

/// C++ `isRowInCurrentPath` lambda (cpp 209-213); `current_path` is `currentPath[0..currPathLen]`.
fn is_row_in_current_path(current_path: &[HighsInt], row: HighsInt) -> bool {
    current_path.contains(&row)
}

/// C++ `checkWeight` lambda (cpp 219-222).
fn check_weight(w: f64, max_weight: f64, min_weight: f64) -> bool {
    let w = w.abs();
    w <= max_weight && w >= min_weight
}

/// C++ `skipCol` lambda (cpp 224-254). `current_path` is `currentPath[0..currPathLen]` and
/// `try_negated_scale` is the C++ `tryNegatedScale` captured by reference.
#[allow(clippy::too_many_arguments)]
fn skip_col(
    i: HighsInt,
    current_path: &[HighsInt],
    try_negated_scale: &mut bool,
    col: HighsInt,
    col_arcs: &[(HighsInt, HighsInt)],
    arc_rows: &[(HighsInt, f64)],
    other_col_arcs: &[(HighsInt, HighsInt)],
    other_arc_rows: &[(HighsInt, f64)],
) -> bool {
    let curr_path_len = current_path.len() as HighsInt;
    if curr_path_len == 1 && !*try_negated_scale {
        let (first, second) = col_arcs[col as usize];
        if second - first <= curr_path_len {
            for k in first..second {
                if arc_rows[k as usize].0 != i {
                    *try_negated_scale = true;
                    break;
                }
            }
        } else {
            *try_negated_scale = true;
        }
    }

    let (other_first, other_second) = other_col_arcs[col as usize];
    if other_first == other_second {
        return true;
    }
    if other_second - other_first <= curr_path_len {
        for k in other_first..other_second {
            if !is_row_in_current_path(current_path, other_arc_rows[k as usize].0) {
                return false;
            }
        }
        return true;
    }
    false
}

/// C++ `findRow` lambda (cpp 256-294). The C++ out-parameters `row` / `weight` are returned as
/// `Some((row, weight))` when it returns true; the random draw happens before any check, as in C++.
#[allow(clippy::too_many_arguments)]
fn find_row(
    best_arc_col: HighsInt,
    val: f64,
    col_arcs: &[(HighsInt, HighsInt)],
    arc_rows: &[(HighsInt, f64)],
    current_path: &[HighsInt],
    randgen: &mut HighsRandom,
    max_weight: f64,
    min_weight: f64,
) -> Option<(HighsInt, f64)> {
    let (first, second) = col_arcs[best_arc_col as usize];
    let arc_row = randgen.integer_range(first, second);
    let r = arc_rows[arc_row as usize].0;
    let w = -val / arc_rows[arc_row as usize].1;
    if !is_row_in_current_path(current_path, r) && check_weight(w, max_weight, min_weight) {
        return Some((r, w));
    }

    for next_row in (arc_row + 1)..second {
        let r = arc_rows[next_row as usize].0;
        let w = -val / arc_rows[next_row as usize].1;
        if !is_row_in_current_path(current_path, r) && check_weight(w, max_weight, min_weight) {
            return Some((r, w));
        }
    }

    for next_row in first..arc_row {
        let r = arc_rows[next_row as usize].0;
        let w = -val / arc_rows[next_row as usize].1;
        if !is_row_in_current_path(current_path, r) && check_weight(w, max_weight, min_weight) {
            return Some((r, w));
        }
    }

    None
}

/// HiGHS: HighsPathSeparator.cpp::RowType (lines 19-24; file-local `enum class RowType : int8_t`)
#[repr(i8)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum PathSeparatorRowType {
    kUnusuable = -2,
    kGeq = -1,
    #[default]
    kEq = 0,
    kLeq = 1,
}
