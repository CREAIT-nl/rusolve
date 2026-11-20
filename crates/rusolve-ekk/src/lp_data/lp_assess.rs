//! HiGHS source: lp_data/HighsLpUtils.cpp (assessLp, lpDimensionsOk, assessCosts, assessBounds)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! LP validation and normalisation, as run by `Highs::passModel`: costs with
//! |c| >= infinite_cost become +-Inf, bounds beyond +-infinite_bound become
//! +-Inf, and the matrix drops |a| <= small_matrix_value.

#![allow(
    clippy::int_plus_one,
    clippy::needless_late_init,
    clippy::collapsible_match,
    clippy::assertions_on_constants
)]

use crate::lp_data::hstruct::HighsScale;
use crate::prelude::*;
use crate::util::highs_matrix_utils::assess_matrix_dimensions;

/// HiGHS: HighsLpUtils.cpp::assessLp (lines 33-105, 73 lines)
pub fn assess_lp(lp: &mut HighsLp, options: &HighsOptions) -> HighsStatus {
    let mut return_status = HighsStatus::kOk;
    let mut call_status = if lp_dimensions_ok("assessLp", lp, &options.log_options) {
        HighsStatus::kOk
    } else {
        HighsStatus::kError
    };
    return_status = interpret_call_status(
        &options.log_options,
        call_status,
        return_status,
        "assessLpDimensions",
    );
    if return_status == HighsStatus::kError {
        return return_status;
    }

    if lp.num_col_ != 0 {
        // Assess the LP column costs
        let index_collection = HighsIndexCollection {
            dimension_: lp.num_col_,
            is_interval_: true,
            from_: 0,
            to_: lp.num_col_ - 1,
            ..HighsIndexCollection::default()
        };
        call_status = assess_costs(
            options,
            0,
            &index_collection,
            &mut lp.col_cost_,
            &mut lp.has_infinite_cost_,
            options.infinite_cost,
        );
        return_status = interpret_call_status(
            &options.log_options,
            call_status,
            return_status,
            "assessCosts",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        // Assess the LP column bounds
        let is_mip = lp.is_mip();
        call_status = assess_bounds(
            options,
            "Col",
            0,
            &index_collection,
            &mut lp.col_lower_,
            &mut lp.col_upper_,
            options.infinite_bound,
            if is_mip {
                Some(&lp.integrality_[..])
            } else {
                None
            },
        );
        return_status = interpret_call_status(
            &options.log_options,
            call_status,
            return_status,
            "assessBounds",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
    }
    if lp.num_row_ != 0 {
        // Assess the LP row bounds
        let index_collection = HighsIndexCollection {
            dimension_: lp.num_row_,
            is_interval_: true,
            from_: 0,
            to_: lp.num_row_ - 1,
            ..HighsIndexCollection::default()
        };
        call_status = assess_bounds(
            options,
            "Row",
            0,
            &index_collection,
            &mut lp.row_lower_,
            &mut lp.row_upper_,
            options.infinite_bound,
            None,
        );
        return_status = interpret_call_status(
            &options.log_options,
            call_status,
            return_status,
            "assessBounds",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
    }
    // If the LP has no columns the matrix must be empty and there is
    // nothing left to test
    if lp.num_col_ == 0 {
        debug_assert!(lp.a_matrix_.num_nz() == 0);
        return HighsStatus::kOk;
    }
    // From here, any LP has lp.num_col_ > 0 and lp.a_matrix_.start_[lp.num_col_]
    // exists (as the number of nonzeros)
    debug_assert!(lp.num_col_ > 0);

    // Assess the LP matrix - even if there are no rows!
    call_status = lp.a_matrix_.assess(
        &options.log_options,
        "LP",
        options.small_matrix_value,
        options.large_matrix_value,
        options.h16_assess_small,
    );
    return_status = interpret_call_status(
        &options.log_options,
        call_status,
        return_status,
        "assessMatrix",
    );
    if return_status == HighsStatus::kError {
        return return_status;
    }
    // If entries have been removed from the matrix, resize the index
    // and value vectors to prevent bug in presolve
    let lp_num_nz = lp.a_matrix_.num_nz();
    if (lp.a_matrix_.index_.len() as HighsInt) > lp_num_nz {
        lp.a_matrix_.index_.resize(lp_num_nz as usize, 0);
    }
    if (lp.a_matrix_.value_.len() as HighsInt) > lp_num_nz {
        lp.a_matrix_.value_.resize(lp_num_nz as usize, 0.0);
    }
    if return_status != HighsStatus::kOk {
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "assessLp returns HighsStatus = {}\n",
            highs_status_to_string(return_status)
        );
    }
    return_status
}

/// HiGHS: HighsLpUtils.cpp::lpDimensionsOk (lines 107-268, 162 lines)
pub fn lp_dimensions_ok(message: &str, lp: &HighsLp, log_options: &HighsLogOptions) -> bool {
    let mut ok = true;
    let num_col = lp.num_col_;
    let num_row = lp.num_row_;
    if num_col < 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on num_col = {} >= 0\n",
            message,
            num_col
        );
    }
    ok = num_col >= 0 && ok;
    if num_row < 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on num_row = {} >= 0\n",
            message,
            num_row
        );
    }
    ok = num_row >= 0 && ok;
    if !ok {
        return ok;
    }

    let col_cost_size = lp.col_cost_.len() as HighsInt;
    let col_lower_size = lp.col_lower_.len() as HighsInt;
    let col_upper_size = lp.col_upper_.len() as HighsInt;
    let legal_col_cost_size = col_cost_size >= num_col;
    let legal_col_lower_size = col_lower_size >= num_col;
    let legal_col_upper_size = col_upper_size >= num_col;
    if !legal_col_cost_size {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on col_cost.size() = {} < {} = num_col\n",
            message,
            col_cost_size,
            num_col
        );
    }
    ok = legal_col_cost_size && ok;
    if !legal_col_lower_size {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on col_lower.size() = {} < {} = num_col\n",
            message,
            col_lower_size,
            num_col
        );
    }
    ok = legal_col_lower_size && ok;
    if !legal_col_upper_size {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on col_upper.size() = {} < {} = num_col\n",
            message,
            col_upper_size,
            num_col
        );
    }
    ok = legal_col_upper_size && ok;

    let legal_format = lp.a_matrix_.format_ == MatrixFormat::kColwise
        || lp.a_matrix_.format_ == MatrixFormat::kRowwise;
    if !legal_format {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on a_matrix_.format\n",
            message
        );
    }
    ok = legal_format && ok;
    let num_vec = if lp.a_matrix_.is_colwise() {
        num_col
    } else {
        num_row
    };
    let partitioned = false;
    let a_matrix_p_end: Vec<HighsInt> = Vec::new();
    let legal_matrix_dimensions = assess_matrix_dimensions(
        log_options,
        num_vec,
        partitioned,
        &lp.a_matrix_.start_,
        &a_matrix_p_end,
        &lp.a_matrix_.index_,
        &lp.a_matrix_.value_,
    ) == HighsStatus::kOk;
    if !legal_matrix_dimensions {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on a_matrix dimensions\n",
            message
        );
    }
    ok = legal_matrix_dimensions && ok;

    let row_lower_size = lp.row_lower_.len() as HighsInt;
    let row_upper_size = lp.row_upper_.len() as HighsInt;
    let legal_row_lower_size = row_lower_size >= num_row;
    let legal_row_upper_size = row_upper_size >= num_row;
    if !legal_row_lower_size {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on row_lower.size() = {} < {} = num_row\n",
            message,
            row_lower_size,
            num_row
        );
    }
    ok = legal_row_lower_size && ok;
    if !legal_row_upper_size {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on row_upper.size() = {} < {} = num_row\n",
            message,
            row_upper_size,
            num_row
        );
    }
    ok = legal_row_upper_size && ok;

    let legal_a_matrix_num_col = lp.a_matrix_.num_col_ == num_col;
    let legal_a_matrix_num_row = lp.a_matrix_.num_row_ == num_row;
    if !legal_a_matrix_num_col {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on a_matrix.num_col_ = {} != {} = num_col\n",
            message,
            lp.a_matrix_.num_col_,
            num_col
        );
    }
    ok = legal_a_matrix_num_col && ok;
    if !legal_a_matrix_num_row {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on a_matrix.num_row_ = {} != {} = num_row\n",
            message,
            lp.a_matrix_.num_row_,
            num_row
        );
    }
    ok = legal_a_matrix_num_row && ok;

    let scale: &HighsScale = &lp.scale_;
    let scale_strategy = scale.strategy;
    let legal_scale_strategy = scale_strategy >= 0;
    if !legal_scale_strategy {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on scale_.scale_strategy\n",
            message
        );
    }
    ok = legal_scale_strategy && ok;
    let scale_row_size = scale.row.len() as HighsInt;
    let scale_col_size = scale.col.len() as HighsInt;
    let legal_scale_num_col;
    let legal_scale_num_row;
    let legal_scale_row_size;
    let legal_scale_col_size;
    if scale.has_scaling {
        legal_scale_num_col = scale.num_col == num_col;
        legal_scale_num_row = scale.num_row == num_row;
        legal_scale_row_size = scale_row_size >= num_row;
        legal_scale_col_size = scale_col_size >= num_col;
    } else {
        legal_scale_num_col = scale.num_col == 0;
        legal_scale_num_row = scale.num_row == 0;
        legal_scale_row_size = scale_row_size == 0;
        legal_scale_col_size = scale_col_size == 0;
    }
    if !legal_scale_num_col {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on scale_.num_col = {} != {}\n",
            message,
            scale.num_col,
            if scale.has_scaling { num_col } else { 0 }
        );
    }
    ok = legal_scale_num_col && ok;
    if !legal_scale_num_row {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on scale_.num_row = {} != {}\n",
            message,
            scale.num_row,
            if scale.has_scaling { num_row } else { 0 }
        );
    }
    ok = legal_scale_num_row && ok;
    if !legal_scale_col_size {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on scale_.col.size() = {} {} {}\n",
            message,
            scale_col_size,
            if scale.has_scaling { ">=" } else { "==" },
            if scale.has_scaling { num_col } else { 0 }
        );
    }
    ok = legal_scale_col_size && ok;
    if !legal_scale_row_size {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails on scale_.row.size() = {} {} {}\n",
            message,
            scale_row_size,
            if scale.has_scaling { ">=" } else { "==" },
            if scale.has_scaling { num_row } else { 0 }
        );
    }
    ok = legal_scale_row_size && ok;
    if !ok {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "LP dimension validation ({}) fails\n",
            message
        );
    }

    ok
}

/// HiGHS: HighsLpUtils.cpp::assessCosts (lines 270-336, 67 lines)
pub fn assess_costs(
    options: &HighsOptions,
    ml_col_os: HighsInt,
    index_collection: &HighsIndexCollection,
    cost: &mut [f64],
    has_infinite_cost: &mut bool,
    infinite_cost: f64,
) -> HighsStatus {
    let _ = ml_col_os;
    let return_status = HighsStatus::kOk;
    debug_assert!(ok(index_collection));
    let mut from_k: HighsInt = 0;
    let mut to_k: HighsInt = 0;
    limits(index_collection, &mut from_k, &mut to_k);
    if from_k > to_k {
        return return_status;
    }

    let mut local_col: HighsInt;
    let mut usr_col: HighsInt = -1;
    let mut num_infinite_cost: HighsInt = 0;
    for k in from_k..to_k + 1 {
        if index_collection.is_interval_ || index_collection.is_mask_ {
            local_col = k;
        } else {
            local_col = index_collection.set_[k as usize];
        }
        if index_collection.is_interval_ {
            usr_col += 1;
        } else {
            usr_col = k;
        }
        if index_collection.is_mask_ && index_collection.mask_[local_col as usize] == 0 {
            continue;
        }
        let uc = usr_col as usize;
        if cost[uc] >= infinite_cost {
            num_infinite_cost += 1;
            cost[uc] = kHighsInf;
        } else if cost[uc] <= -infinite_cost {
            num_infinite_cost += 1;
            cost[uc] = -kHighsInf;
        }
    }
    if num_infinite_cost > 0 {
        *has_infinite_cost = true;
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "{} |cost| values greater than or equal to {:12} are treated as Infinity\n",
            num_infinite_cost,
            infinite_cost
        );
    }
    return_status
}

/// HiGHS: HighsLpUtils.cpp::assessBounds (lines 338-470, 133 lines)
/// C++ default arguments: integrality = nullptr
#[allow(clippy::too_many_arguments)]
pub fn assess_bounds(
    options: &HighsOptions,
    type_: &str,
    ml_ix_os: HighsInt,
    index_collection: &HighsIndexCollection,
    lower: &mut [f64],
    upper: &mut [f64],
    infinite_bound: f64,
    integrality: Option<&[HighsVarType]>,
) -> HighsStatus {
    debug_assert!(ok(index_collection));
    let mut from_k: HighsInt = 0;
    let mut to_k: HighsInt = 0;
    limits(index_collection, &mut from_k, &mut to_k);
    if from_k > to_k {
        return HighsStatus::kOk;
    }

    let mut error_found = false;
    let mut warning_found = false;
    let mut num_infinite_lower_bound: HighsInt = 0;
    let mut num_infinite_upper_bound: HighsInt = 0;
    let mut local_ix: HighsInt;
    let mut ml_ix: HighsInt;
    let mut usr_ix: HighsInt = -1;
    for k in from_k..to_k + 1 {
        if index_collection.is_interval_ || index_collection.is_mask_ {
            local_ix = k;
        } else {
            local_ix = index_collection.set_[k as usize];
        }
        if index_collection.is_interval_ {
            usr_ix += 1;
        } else {
            usr_ix = k;
        }
        ml_ix = ml_ix_os + local_ix;
        if index_collection.is_mask_ && index_collection.mask_[local_ix as usize] == 0 {
            continue;
        }
        let ux = usr_ix as usize;

        if !highs_is_infinity(-lower[ux]) {
            // Check whether a finite lower bound will be treated as -Infinity
            let infinite_lower_bound = lower[ux] <= -infinite_bound;
            if infinite_lower_bound {
                lower[ux] = -kHighsInf;
                num_infinite_lower_bound += 1;
            }
        }
        if !highs_is_infinity(upper[ux]) {
            // Check whether a finite upper bound will be treated as Infinity
            let infinite_upper_bound = upper[ux] >= infinite_bound;
            if infinite_upper_bound {
                upper[ux] = kHighsInf;
                num_infinite_upper_bound += 1;
            }
        }
        // Check that the lower bound does not exceed the upper bound
        let mut legal_lower_upper_bound = lower[ux] <= upper[ux];
        if let Some(integrality) = integrality {
            // Legal for semi-variables to have inconsistent bounds
            if integrality[ux] == HighsVarType::kSemiContinuous
                || integrality[ux] == HighsVarType::kSemiInteger
            {
                legal_lower_upper_bound = true;
            }
        }
        if !legal_lower_upper_bound {
            // Leave inconsistent bounds to be used to deduce infeasibility
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kWarning,
                "{:3}  {:12} has inconsistent bounds [{:12}, {:12}]\n",
                type_,
                ml_ix,
                lower[ux],
                upper[ux]
            );
            warning_found = true;
        }
        // Check that the lower bound is not as much as +Infinity
        let legal_lower_bound = lower[ux] < infinite_bound;
        if !legal_lower_bound {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "{:3}  {:12} has lower bound of {:12} >= {:12}\n",
                type_,
                ml_ix,
                lower[ux],
                infinite_bound
            );
            error_found = true;
        }
        // Check that the upper bound is not as little as -Infinity
        let legal_upper_bound = upper[ux] > -infinite_bound;
        if !legal_upper_bound {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "{:3}  {:12} has upper bound of {:12} <= {:12}\n",
                type_,
                ml_ix,
                upper[ux],
                -infinite_bound
            );
            error_found = true;
        }
    }
    if num_infinite_lower_bound != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "{:3}s:{:12} lower bounds    less than or equal to {:12} are treated as -Infinity\n",
            type_,
            num_infinite_lower_bound,
            -infinite_bound
        );
    }
    if num_infinite_upper_bound != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "{:3}s:{:12} upper bounds greater than or equal to {:12} are treated as +Infinity\n",
            type_,
            num_infinite_upper_bound,
            infinite_bound
        );
    }

    if error_found {
        HighsStatus::kError
    } else if warning_found {
        HighsStatus::kWarning
    } else {
        HighsStatus::kOk
    }
}
