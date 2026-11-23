//! HiGHS source: lp_data/HighsLpUtils.cpp (considerScaling, scaleLp,
//! equilibrationScaleMatrix, maxValueScaleMatrix, unscaleSolution)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! LP scaling as run by `solveLpSimplex`. `analyseLp` (only reached when
//! `highs_analysis_level & kHighsAnalysisLevelModelData`) is out of scope
//! and is a no-op here.

#![allow(
    clippy::int_plus_one,
    clippy::needless_late_init,
    clippy::collapsible_match,
    clippy::assertions_on_constants
)]

use crate::prelude::*;

/// HiGHS: HighsLpUtils.cpp::considerScaling (lines 926-964, 39 lines)
pub fn consider_scaling(options: &HighsOptions, lp: &mut HighsLp) -> bool {
    // Indicate whether new scaling has been determined in the return value.
    let mut new_scaling = false;
    // Consider scaling the LP - either by finding new factors or by
    // applying any existing factors
    let allow_scaling =
        lp.num_col_ > 0 && options.simplex_scale_strategy != kSimplexScaleStrategyOff;
    if lp.scale_.has_scaling && !allow_scaling {
        // LP had scaling before, but now it is not permitted, clear any
        // scaling. Return true as scaling position has changed
        lp.clear_scale();
        return true;
    }
    let scaling_not_tried = lp.scale_.strategy == kSimplexScaleStrategyOff;
    let new_scaling_strategy = options.simplex_scale_strategy != lp.scale_.strategy
        && options.simplex_scale_strategy != kSimplexScaleStrategyChoose;
    let try_scaling = allow_scaling && (scaling_not_tried || new_scaling_strategy);
    if try_scaling {
        // Scaling will be tried, so ensure that any previous scaling is not applied
        lp.unapply_scale();
        let analyse_lp_data = (kHighsAnalysisLevelModelData & options.highs_analysis_level) != 0;
        // analyseLp: out of scope (no-op)
        scale_lp(options, lp, false);
        // If the LP is now scaled, then the scaling is new
        new_scaling = lp.is_scaled_;
    } else if lp.scale_.has_scaling {
        // Scaling factors are known, so ensure that they are applied
        lp.apply_scale();
    }
    // Ensure that either the LP has scale factors and is scaled, or
    // it doesn't have scale factors and isn't scaled
    debug_assert!(lp.scale_.has_scaling == lp.is_scaled_);
    new_scaling
}

/// HiGHS: HighsLpUtils.cpp::scaleLp (lines 966-1062, 97 lines)
/// C++ default arguments: force_scaling = false
#[allow(unused_assignments)]
pub fn scale_lp(options: &HighsOptions, lp: &mut HighsLp, force_scaling: bool) {
    lp.clear_scaling();
    let num_col = lp.num_col_;
    let num_row = lp.num_row_;
    // Scaling not well defined for models with no columns
    debug_assert!(num_col > 0);

    // Save the simplex_scale_strategy so that the option can be
    // modified for the course of this method
    let simplex_scale_strategy = options.simplex_scale_strategy;
    // Determine the actual strategy to use
    let mut use_scale_strategy = simplex_scale_strategy;
    if use_scale_strategy == kSimplexScaleStrategyChoose {
        // HiGHS is left to choose: currently use forced equilibration, but maybe do
        // something more intelligent
        use_scale_strategy = kSimplexScaleStrategyForcedEquilibration;
    }
    // Find out range of matrix values and skip matrix scaling if all
    // |values| are in [0.2, 5]
    let no_scaling_original_matrix_min_value = 0.2;
    let no_scaling_original_matrix_max_value = 5.0;
    let mut original_matrix_min_value = kHighsInf;
    let mut original_matrix_max_value = 0.0;
    lp.a_matrix_.range(
        &mut original_matrix_min_value,
        &mut original_matrix_max_value,
    );
    // Possibly force scaling, otherwise base the decision on the range
    // of values in the matrix, values that will be used later for
    // reporting
    let no_scaling = if force_scaling {
        false
    } else {
        (original_matrix_min_value >= no_scaling_original_matrix_min_value)
            && (original_matrix_max_value <= no_scaling_original_matrix_max_value)
    };
    let mut scaled_matrix = false;
    if no_scaling {
        // No matrix scaling, but possible cost scaling
        if options.highs_analysis_level != 0 {
            crate::highs_log_dev!(
                &options.log_options,
                HighsLogType::kInfo,
                "Scaling: Matrix has [min, max] values of [{}, {}] within [{}, {}] so no scaling performed\n",
                original_matrix_min_value,
                original_matrix_max_value,
                no_scaling_original_matrix_min_value,
                no_scaling_original_matrix_max_value
            );
        }
    } else {
        // Try scaling, so assign unit factors - partly because initial
        // factors may be assumed by the scaling method, but also because
        // scaling factors may not be computed for empty rows/columns
        lp.scale_.col.clear();
        lp.scale_.col.resize(num_col as usize, 1.0);
        lp.scale_.row.clear();
        lp.scale_.row.resize(num_row as usize, 1.0);
        let equilibration_scaling = use_scale_strategy == kSimplexScaleStrategyEquilibration
            || use_scale_strategy == kSimplexScaleStrategyForcedEquilibration;
        // Try scaling. Value of scaled_matrix indicates whether scaling
        // was considered valuable (and performed). If it's not valuable
        // then the matrix remains unscaled
        if equilibration_scaling {
            scaled_matrix = equilibration_scale_matrix(options, lp, use_scale_strategy);
        } else {
            scaled_matrix = max_value_scale_matrix(options, lp, use_scale_strategy);
        }
        if scaled_matrix {
            // Matrix is scaled, so scale the bounds and costs
            for i_col in 0..num_col as usize {
                lp.col_lower_[i_col] /= lp.scale_.col[i_col];
                lp.col_upper_[i_col] /= lp.scale_.col[i_col];
                lp.col_cost_[i_col] *= lp.scale_.col[i_col];
            }
            for i_row in 0..num_row as usize {
                lp.row_lower_[i_row] *= lp.scale_.row[i_row];
                lp.row_upper_[i_row] *= lp.scale_.row[i_row];
            }
            lp.scale_.has_scaling = true;
            lp.scale_.num_col = num_col;
            lp.scale_.num_row = num_row;
            lp.scale_.cost = 1.0;
            lp.is_scaled_ = true;
        } else {
            // Matrix is not scaled, so clear the scaling
            lp.clear_scaling();
        }
    }
    // Record the scaling strategy used
    lp.scale_.strategy = use_scale_strategy;
}

/// HiGHS: HighsLpUtils.cpp::equilibrationScaleMatrix (lines 1064-1385, 322 lines)
pub fn equilibration_scale_matrix(
    options: &HighsOptions,
    lp: &mut HighsLp,
    use_scale_strategy: HighsInt,
) -> bool {
    let num_col = lp.num_col_;
    let num_row = lp.num_row_;
    let nc = num_col as usize;
    let nr = num_row as usize;
    let col_cost = &lp.col_cost_;
    let scale = &mut lp.scale_;
    let col_scale = &mut scale.col;
    let row_scale = &mut scale.row;
    let a_start = &lp.a_matrix_.start_;
    let a_index = &lp.a_matrix_.index_;
    let a_value = &mut lp.a_matrix_.value_;

    let simplex_scale_strategy = use_scale_strategy;

    let mut original_matrix_min_value = kHighsInf;
    let mut original_matrix_max_value: f64 = 0.0;
    let an_x = a_start[nc];
    for k in 0..an_x as usize {
        let value = a_value[k].abs();
        original_matrix_min_value = cpp_min(original_matrix_min_value, value);
        original_matrix_max_value = cpp_max(original_matrix_max_value, value);
    }

    // Include cost in scaling if minimum nonzero cost is less than 0.1
    let mut min_nonzero_cost = kHighsInf;
    for i in 0..nc {
        if col_cost[i] != 0.0 {
            min_nonzero_cost = cpp_min(col_cost[i].abs(), min_nonzero_cost);
        }
    }
    let include_cost_in_scaling = min_nonzero_cost < 0.1;

    // Limits on scaling factors
    // Now that kHighsInf =
    // std::numeric_limits<double>::infinity(), this Qi-trick doesn't
    // work so, in recognition, use the old value of kHighsInf
    let finite_infinity = 1e200;
    let max_allow_scale = 2.0_f64.powf(options.allowed_matrix_scale_factor as f64);
    let min_allow_scale = 1.0 / max_allow_scale;

    let min_allow_col_scale = min_allow_scale;
    let max_allow_col_scale = max_allow_scale;
    let min_allow_row_scale = min_allow_scale;
    let max_allow_row_scale = max_allow_scale;

    // Search up to 6 times
    let mut row_min_value: Vec<f64> = vec![finite_infinity; nr];
    let mut row_max_value: Vec<f64> = vec![1.0 / finite_infinity; nr];
    for _search_count in 0..6 {
        // Find column scale, prepare row data
        for i_col in 0..nc {
            // For column scale (find)
            let mut col_min_value = finite_infinity;
            let mut col_max_value = 1.0 / finite_infinity;
            let abs_col_cost = col_cost[i_col].abs();
            if include_cost_in_scaling && abs_col_cost != 0.0 {
                col_min_value = cpp_min(col_min_value, abs_col_cost);
                col_max_value = cpp_max(col_max_value, abs_col_cost);
            }
            for k in a_start[i_col] as usize..a_start[i_col + 1] as usize {
                let value = a_value[k].abs() * row_scale[a_index[k] as usize];
                col_min_value = cpp_min(col_min_value, value);
                col_max_value = cpp_max(col_max_value, value);
            }
            let col_equilibration = 1.0 / (col_min_value * col_max_value).sqrt();
            // Ensure that column scale factor is not excessively large or small
            col_scale[i_col] = cpp_min(
                cpp_max(min_allow_col_scale, col_equilibration),
                max_allow_col_scale,
            );
            // For row scale (only collect)
            for k in a_start[i_col] as usize..a_start[i_col + 1] as usize {
                let i_row = a_index[k] as usize;
                let value = a_value[k].abs() * col_scale[i_col];
                row_min_value[i_row] = cpp_min(row_min_value[i_row], value);
                row_max_value[i_row] = cpp_max(row_max_value[i_row], value);
            }
        }
        // For row scale (find)
        for i_row in 0..nr {
            let row_equilibration = 1.0 / (row_min_value[i_row] * row_max_value[i_row]).sqrt();
            // Ensure that row scale factor is not excessively large or small
            row_scale[i_row] = cpp_min(
                cpp_max(min_allow_row_scale, row_equilibration),
                max_allow_row_scale,
            );
        }
        row_min_value.clear();
        row_min_value.resize(nr, finite_infinity);
        row_max_value.clear();
        row_max_value.resize(nr, 1.0 / finite_infinity);
    }
    // Make it numerically better
    // Also determine the max and min row and column scaling factors
    let mut min_col_scale = finite_infinity;
    let mut max_col_scale = 1.0 / finite_infinity;
    let mut min_row_scale = finite_infinity;
    let mut max_row_scale = 1.0 / finite_infinity;
    let log2 = 2.0_f64.ln();
    for i_col in 0..nc {
        col_scale[i_col] = 2.0_f64.powf((col_scale[i_col].ln() / log2 + 0.5).floor());
        min_col_scale = cpp_min(col_scale[i_col], min_col_scale);
        max_col_scale = cpp_max(col_scale[i_col], max_col_scale);
    }
    for i_row in 0..nr {
        row_scale[i_row] = 2.0_f64.powf((row_scale[i_row].ln() / log2 + 0.5).floor());
        min_row_scale = cpp_min(row_scale[i_row], min_row_scale);
        max_row_scale = cpp_max(row_scale[i_row], max_row_scale);
    }
    // Apply scaling to matrix and bounds
    let mut matrix_min_value = finite_infinity;
    let mut matrix_max_value: f64 = 0.0;
    let mut min_original_col_equilibration = finite_infinity;
    let mut sum_original_log_col_equilibration = 0.0;
    let mut max_original_col_equilibration: f64 = 0.0;
    let mut min_original_row_equilibration = finite_infinity;
    let mut sum_original_log_row_equilibration = 0.0;
    let mut max_original_row_equilibration: f64 = 0.0;
    let mut min_col_equilibration = finite_infinity;
    let mut sum_log_col_equilibration = 0.0;
    let mut max_col_equilibration: f64 = 0.0;
    let mut min_row_equilibration = finite_infinity;
    let mut sum_log_row_equilibration = 0.0;
    let mut max_row_equilibration: f64 = 0.0;
    let mut original_row_min_value: Vec<f64> = vec![finite_infinity; nr];
    let mut original_row_max_value: Vec<f64> = vec![1.0 / finite_infinity; nr];
    row_min_value.clear();
    row_min_value.resize(nr, finite_infinity);
    row_max_value.clear();
    row_max_value.resize(nr, 1.0 / finite_infinity);
    for i_col in 0..nc {
        let mut original_col_min_value = finite_infinity;
        let mut original_col_max_value = 1.0 / finite_infinity;
        let mut col_min_value = finite_infinity;
        let mut col_max_value = 1.0 / finite_infinity;
        for k in a_start[i_col] as usize..a_start[i_col + 1] as usize {
            let i_row = a_index[k] as usize;
            let original_value = a_value[k].abs();
            original_col_min_value = cpp_min(original_value, original_col_min_value);
            original_col_max_value = cpp_max(original_value, original_col_max_value);
            original_row_min_value[i_row] = cpp_min(original_row_min_value[i_row], original_value);
            original_row_max_value[i_row] = cpp_max(original_row_max_value[i_row], original_value);
            a_value[k] *= col_scale[i_col] * row_scale[i_row];
            let value = a_value[k].abs();
            col_min_value = cpp_min(value, col_min_value);
            col_max_value = cpp_max(value, col_max_value);
            row_min_value[i_row] = cpp_min(row_min_value[i_row], value);
            row_max_value[i_row] = cpp_max(row_max_value[i_row], value);
        }
        matrix_min_value = cpp_min(matrix_min_value, col_min_value);
        matrix_max_value = cpp_max(matrix_max_value, col_max_value);

        let original_col_equilibration =
            1.0 / (original_col_min_value * original_col_max_value).sqrt();
        min_original_col_equilibration =
            cpp_min(original_col_equilibration, min_original_col_equilibration);
        sum_original_log_col_equilibration += original_col_equilibration.ln();
        max_original_col_equilibration =
            cpp_max(original_col_equilibration, max_original_col_equilibration);
        let col_equilibration = 1.0 / (col_min_value * col_max_value).sqrt();
        min_col_equilibration = cpp_min(col_equilibration, min_col_equilibration);
        sum_log_col_equilibration += col_equilibration.ln();
        max_col_equilibration = cpp_max(col_equilibration, max_col_equilibration);
    }

    for i_row in 0..nr {
        let original_row_equilibration =
            1.0 / (original_row_min_value[i_row] * original_row_max_value[i_row]).sqrt();
        min_original_row_equilibration =
            cpp_min(original_row_equilibration, min_original_row_equilibration);
        sum_original_log_row_equilibration += original_row_equilibration.ln();
        max_original_row_equilibration =
            cpp_max(original_row_equilibration, max_original_row_equilibration);
        let row_equilibration = 1.0 / (row_min_value[i_row] * row_max_value[i_row]).sqrt();
        min_row_equilibration = cpp_min(row_equilibration, min_row_equilibration);
        sum_log_row_equilibration += row_equilibration.ln();
        max_row_equilibration = cpp_max(row_equilibration, max_row_equilibration);
    }
    let geomean_original_col_equilibration =
        (sum_original_log_col_equilibration / num_col as f64).exp();
    let geomean_original_row_equilibration =
        (sum_original_log_row_equilibration / num_row as f64).exp();
    let geomean_col_equilibration = (sum_log_col_equilibration / num_col as f64).exp();
    let geomean_row_equilibration = (sum_log_row_equilibration / num_row as f64).exp();
    if options.log_dev_level != 0 {
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "Scaling: Original equilibration: min/mean/max {:11.4}/{:11.4}/{:11.4} (cols); min/mean/max {:11.4}/{:11.4}/{:11.4} (rows)\n",
            min_original_col_equilibration,
            geomean_original_col_equilibration,
            max_original_col_equilibration,
            min_original_row_equilibration,
            geomean_original_row_equilibration,
            max_original_row_equilibration
        );
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "Scaling: Final    equilibration: min/mean/max {:11.4}/{:11.4}/{:11.4} (cols); min/mean/max {:11.4}/{:11.4}/{:11.4} (rows)\n",
            min_col_equilibration,
            geomean_col_equilibration,
            max_col_equilibration,
            min_row_equilibration,
            geomean_row_equilibration,
            max_row_equilibration
        );
    }

    // Compute the mean equilibration improvement
    let geomean_original_col = cpp_max(
        geomean_original_col_equilibration,
        1.0 / geomean_original_col_equilibration,
    );
    let geomean_original_row = cpp_max(
        geomean_original_row_equilibration,
        1.0 / geomean_original_row_equilibration,
    );
    let geomean_col = cpp_max(geomean_col_equilibration, 1.0 / geomean_col_equilibration);
    let geomean_row = cpp_max(geomean_row_equilibration, 1.0 / geomean_row_equilibration);
    let mean_equilibration_improvement =
        ((geomean_original_col * geomean_original_row) / (geomean_col * geomean_row)).sqrt();
    // Compute the extreme equilibration improvement
    let original_col_ratio = max_original_col_equilibration / min_original_col_equilibration;
    let original_row_ratio = max_original_row_equilibration / min_original_row_equilibration;
    let col_ratio = max_col_equilibration / min_col_equilibration;
    let row_ratio = max_row_equilibration / min_row_equilibration;
    let extreme_equilibration_improvement =
        (original_col_ratio + original_row_ratio) / (col_ratio + row_ratio);
    // Compute the max/min matrix value improvement
    let matrix_value_ratio = matrix_max_value / matrix_min_value;
    let original_matrix_value_ratio = original_matrix_max_value / original_matrix_min_value;
    let matrix_value_ratio_improvement = original_matrix_value_ratio / matrix_value_ratio;
    if options.log_dev_level != 0 {
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "Scaling: Improves    mean equilibration by a factor {:.4}\n",
            mean_equilibration_improvement
        );
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "Scaling: Improves extreme equilibration by a factor {:.4}\n",
            extreme_equilibration_improvement
        );
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "Scaling: Improves max/min matrix values by a factor {:.4}\n",
            matrix_value_ratio_improvement
        );
    }
    let possibly_abandon_scaling =
        simplex_scale_strategy != kSimplexScaleStrategyForcedEquilibration;
    let improvement_factor = extreme_equilibration_improvement
        * mean_equilibration_improvement
        * matrix_value_ratio_improvement;

    let improvement_factor_required = 1.0;
    let poor_improvement = improvement_factor < improvement_factor_required;

    // Possibly abandon scaling if it's not improved equilibration significantly
    if possibly_abandon_scaling && poor_improvement {
        // Unscale the matrix
        for i_col in 0..nc {
            for k in a_start[i_col] as usize..a_start[i_col + 1] as usize {
                let i_row = a_index[k] as usize;
                a_value[k] /= col_scale[i_col] * row_scale[i_row];
            }
        }
        if options.log_dev_level != 0 {
            crate::highs_log_dev!(
                &options.log_options,
                HighsLogType::kInfo,
                "Scaling: Improvement factor {:.4} < {:.4} required, so no scaling applied\n",
                improvement_factor,
                improvement_factor_required
            );
        }
        return false;
    } else if options.log_dev_level != 0 {
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "Scaling: Factors are in [{:.4}, {:.4}] for columns and in [{:.4}, {:.4}] for rows\n",
            min_col_scale,
            max_col_scale,
            min_row_scale,
            max_row_scale
        );
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "Scaling: Improvement factor is {:.4} >= {:.4} so scale LP\n",
            improvement_factor,
            improvement_factor_required
        );
    }
    true
}

/// HiGHS: HighsLpUtils.cpp::maxValueScaleMatrix (lines 1387-1531, 145 lines)
pub fn max_value_scale_matrix(
    options: &HighsOptions,
    lp: &mut HighsLp,
    use_scale_strategy: HighsInt,
) -> bool {
    let num_col = lp.num_col_;
    let num_row = lp.num_row_;
    let nc = num_col as usize;
    let nr = num_row as usize;
    let scale = &mut lp.scale_;
    let col_scale = &mut scale.col;
    let row_scale = &mut scale.row;
    let a_start = &lp.a_matrix_.start_;
    let a_index = &lp.a_matrix_.index_;
    let a_value = &mut lp.a_matrix_.value_;

    debug_assert!(options.simplex_scale_strategy == kSimplexScaleStrategyMaxValue);
    debug_assert!(kSimplexScaleStrategyMaxValue015 == kSimplexScaleStrategyMaxValue);
    debug_assert!(kSimplexScaleStrategyMaxValue0157 == kSimplexScaleStrategyMaxValue);

    let log2 = 2.0_f64.ln();
    let max_allow_scale = 2.0_f64.powf(options.allowed_matrix_scale_factor as f64);
    let min_allow_scale = 1.0 / max_allow_scale;

    let min_allow_col_scale = min_allow_scale;
    let max_allow_col_scale = max_allow_scale;
    let min_allow_row_scale = min_allow_scale;
    let max_allow_row_scale = max_allow_scale;

    let mut min_row_scale = kHighsInf;
    let mut max_row_scale: f64 = 0.0;
    let mut original_matrix_min_value = kHighsInf;
    let mut original_matrix_max_value: f64 = 0.0;
    // Determine the row scaling. Also determine the max/min row scaling
    // factors, and max/min original matrix values
    let mut row_max_value: Vec<f64> = vec![0.0; nr];
    for i_col in 0..nc {
        for k in a_start[i_col] as usize..a_start[i_col + 1] as usize {
            let i_row = a_index[k] as usize;
            let value = a_value[k].abs();
            row_max_value[i_row] = cpp_max(row_max_value[i_row], value);
            original_matrix_min_value = cpp_min(original_matrix_min_value, value);
            original_matrix_max_value = cpp_max(original_matrix_max_value, value);
        }
    }
    for i_row in 0..nr {
        if row_max_value[i_row] != 0.0 {
            let mut row_scale_value = 1.0 / row_max_value[i_row];
            // Convert the row scale factor to the nearest power of two, and
            // ensure that it is not excessively large or small
            row_scale_value = 2.0_f64.powf((row_scale_value.ln() / log2 + 0.5).floor());
            row_scale_value = cpp_min(
                cpp_max(min_allow_row_scale, row_scale_value),
                max_allow_row_scale,
            );
            min_row_scale = cpp_min(row_scale_value, min_row_scale);
            max_row_scale = cpp_max(row_scale_value, max_row_scale);
            row_scale[i_row] = row_scale_value;
        }
    }
    // Determine the column scaling, whilst applying the row scaling
    // Also determine the max/min column scaling factors, and max/min
    // matrix values
    let mut min_col_scale = kHighsInf;
    let mut max_col_scale: f64 = 0.0;
    let mut matrix_min_value = kHighsInf;
    let mut matrix_max_value: f64 = 0.0;
    for i_col in 0..nc {
        let mut col_max_value: f64 = 0.0;
        for k in a_start[i_col] as usize..a_start[i_col + 1] as usize {
            let i_row = a_index[k] as usize;
            a_value[k] *= row_scale[i_row];
            let value = a_value[k].abs();
            col_max_value = cpp_max(col_max_value, value);
        }
        if col_max_value != 0.0 {
            let mut col_scale_value = 1.0 / col_max_value;
            // Convert the col scale factor to the nearest power of two, and
            // ensure that it is not excessively large or small
            col_scale_value = 2.0_f64.powf((col_scale_value.ln() / log2 + 0.5).floor());
            col_scale_value = cpp_min(
                cpp_max(min_allow_col_scale, col_scale_value),
                max_allow_col_scale,
            );
            min_col_scale = cpp_min(col_scale_value, min_col_scale);
            max_col_scale = cpp_max(col_scale_value, max_col_scale);
            col_scale[i_col] = col_scale_value;
            for k in a_start[i_col] as usize..a_start[i_col + 1] as usize {
                a_value[k] *= col_scale[i_col];
                let value = a_value[k].abs();
                matrix_min_value = cpp_min(matrix_min_value, value);
                matrix_max_value = cpp_max(matrix_max_value, value);
            }
        }
    }
    let matrix_value_ratio = matrix_max_value / matrix_min_value;
    let original_matrix_value_ratio = original_matrix_max_value / original_matrix_min_value;
    let matrix_value_ratio_improvement = original_matrix_value_ratio / matrix_value_ratio;

    let improvement_factor = matrix_value_ratio_improvement;

    let improvement_factor_required = 1.0;
    let poor_improvement = improvement_factor <= improvement_factor_required;

    if poor_improvement {
        // Unscale the matrix
        for i_col in 0..nc {
            for k in a_start[i_col] as usize..a_start[i_col + 1] as usize {
                let i_row = a_index[k] as usize;
                a_value[k] /= col_scale[i_col] * row_scale[i_row];
            }
        }
        if options.log_dev_level != 0 {
            crate::highs_log_dev!(
                &options.log_options,
                HighsLogType::kInfo,
                "Scaling: Improvement factor {:.4} < {:.4} required, so no scaling applied\n",
                improvement_factor,
                improvement_factor_required
            );
        }
        false
    } else {
        if options.log_dev_level != 0 {
            crate::highs_log_dev!(
                &options.log_options,
                HighsLogType::kInfo,
                "Scaling: Factors are in [{:.4}, {:.4}] for columns and in [{:.4}, {:.4}] for rows\n",
                min_col_scale,
                max_col_scale,
                min_row_scale,
                max_row_scale
            );
        }
        true
    }
}

/// HiGHS: HighsLpUtils.cpp::unscaleSolution (lines 1571-1586, 16 lines)
pub fn unscale_solution(solution: &mut HighsSolution, scale: &HighsScale) {
    debug_assert!(scale.has_scaling);
    debug_assert!(solution.col_value.len() == scale.num_col as usize);
    debug_assert!(solution.col_dual.len() == scale.num_col as usize);
    debug_assert!(solution.row_value.len() == scale.num_row as usize);
    debug_assert!(solution.row_dual.len() == scale.num_row as usize);

    for i_col in 0..scale.num_col as usize {
        solution.col_value[i_col] *= scale.col[i_col];
        solution.col_dual[i_col] /= scale.col[i_col] / scale.cost;
    }
    for i_row in 0..scale.num_row as usize {
        solution.row_value[i_row] /= scale.row[i_row];
        solution.row_dual[i_row] *= scale.row[i_row] * scale.cost;
    }
}
