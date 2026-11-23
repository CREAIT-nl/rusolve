//! HiGHS source: simplex/HSimplex.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (ci) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Free simplex helpers plus `SimplexBasis::clear/setup`.

use crate::prelude::*;

impl SimplexBasis {
    /// HiGHS: HSimplex.cpp::SimplexBasis::clear (lines 26-34, 9 lines)
    pub fn clear(&mut self) {
        self.hash = 0;
        self.basic_index_.clear();
        self.nonbasic_flag_.clear();
        self.nonbasic_move_.clear();
        self.debug_id = -1;
        self.debug_update_count = -1;
        self.debug_origin_name = "None".to_string();
    }

    /// HiGHS: HSimplex.cpp::SimplexBasis::setup (lines 36-44, 9 lines)
    pub fn setup(&mut self, num_col: HighsInt, num_row: HighsInt) {
        self.hash = 0;
        self.basic_index_.resize(num_row as usize, 0);
        self.nonbasic_flag_.resize((num_col + num_row) as usize, 0);
        self.nonbasic_move_.resize((num_col + num_row) as usize, 0);
        self.debug_id = -1;
        self.debug_update_count = -1;
        self.debug_origin_name = "None".to_string();
    }
}

/// HiGHS: HSimplex.cpp::appendNonbasicColsToBasis (lines 46-66, 21 lines)
/// Overload taking `HighsBasis&`.
pub fn append_nonbasic_cols_to_basis(
    lp: &mut HighsLp,
    highs_basis: &mut HighsBasis,
    xnum_new_col: HighsInt,
) {
    debug_assert!(highs_basis.valid);
    if !highs_basis.valid {
        print!("\n!!Appending columns to invalid basis!!\n\n");
    }
    // Add nonbasic structurals
    if xnum_new_col == 0 {
        return;
    }
    let new_num_col = lp.num_col_ + xnum_new_col;
    highs_basis
        .col_status
        .resize(new_num_col as usize, HighsBasisStatus::default());
    // Make any new columns nonbasic
    for i_col in lp.num_col_..new_num_col {
        let ci = i_col as usize;
        if !highs_is_infinity(-lp.col_lower_[ci]) {
            highs_basis.col_status[ci] = HighsBasisStatus::kLower;
        } else if !highs_is_infinity(lp.col_upper_[ci]) {
            highs_basis.col_status[ci] = HighsBasisStatus::kUpper;
        } else {
            highs_basis.col_status[ci] = HighsBasisStatus::kZero;
        }
    }
}

/// HiGHS: HSimplex.cpp::appendNonbasicColsToBasis (lines 68-121, 54 lines)
/// Overload taking `SimplexBasis&`.
pub fn append_nonbasic_cols_to_basis_simplex(
    lp: &mut HighsLp,
    basis: &mut SimplexBasis,
    xnum_new_col: HighsInt,
) {
    // Add nonbasic structurals
    if xnum_new_col == 0 {
        return;
    }
    let new_num_col = lp.num_col_ + xnum_new_col;
    let new_num_tot = new_num_col + lp.num_row_;
    basis.nonbasic_flag_.resize(new_num_tot as usize, 0);
    basis.nonbasic_move_.resize(new_num_tot as usize, 0);
    // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
    // necessary
    for i_row in (0..lp.num_row_).rev() {
        let r = i_row as usize;
        let i_col = basis.basic_index_[r];
        if i_col >= lp.num_col_ {
            // This basic variable is a row, so shift its index
            basis.basic_index_[r] += xnum_new_col;
        }
        basis.nonbasic_flag_[(new_num_col + i_row) as usize] =
            basis.nonbasic_flag_[(lp.num_col_ + i_row) as usize];
        basis.nonbasic_move_[(new_num_col + i_row) as usize] =
            basis.nonbasic_move_[(lp.num_col_ + i_row) as usize];
    }
    // Make any new columns nonbasic
    for i_col in lp.num_col_..new_num_col {
        let ci = i_col as usize;
        basis.nonbasic_flag_[ci] = kNonbasicFlagTrue;
        let lower = lp.col_lower_[ci];
        let upper = lp.col_upper_[ci];
        let move_ = if lower == upper {
            // Fixed
            kNonbasicMoveZe
        } else if !highs_is_infinity(-lower) {
            // Finite lower bound so boxed or lower
            if !highs_is_infinity(upper) {
                // Finite upper bound so boxed
                if lower.abs() < upper.abs() {
                    kNonbasicMoveUp
                } else {
                    kNonbasicMoveDn
                }
            } else {
                // Lower (since upper bound is infinite)
                kNonbasicMoveUp
            }
        } else if !highs_is_infinity(upper) {
            // Upper
            kNonbasicMoveDn
        } else {
            // FREE
            kNonbasicMoveZe
        };
        debug_assert!(move_ != kIllegalMoveValue);
        basis.nonbasic_move_[ci] = move_;
    }
}

/// HiGHS: HSimplex.cpp::appendBasicRowsToBasis (lines 123-137, 15 lines)
/// Overload taking `HighsBasis&`.
pub fn append_basic_rows_to_basis(
    lp: &mut HighsLp,
    highs_basis: &mut HighsBasis,
    xnum_new_row: HighsInt,
) {
    debug_assert!(highs_basis.valid);
    if !highs_basis.valid {
        print!("\n!!Appending rows to invalid basis!!\n\n");
    }
    // Add basic logicals
    if xnum_new_row == 0 {
        return;
    }
    let new_num_row = lp.num_row_ + xnum_new_row;
    highs_basis
        .row_status
        .resize(new_num_row as usize, HighsBasisStatus::default());
    // Make the new rows basic
    for i_row in lp.num_row_..new_num_row {
        highs_basis.row_status[i_row as usize] = HighsBasisStatus::kBasic;
    }
}

/// HiGHS: HSimplex.cpp::appendBasicRowsToBasis (lines 139-155, 17 lines)
/// Overload taking `SimplexBasis&`.
pub fn append_basic_rows_to_basis_simplex(
    lp: &mut HighsLp,
    basis: &mut SimplexBasis,
    xnum_new_row: HighsInt,
) {
    // Add basic logicals
    if xnum_new_row == 0 {
        return;
    }

    let new_num_row = lp.num_row_ + xnum_new_row;
    let new_num_tot = lp.num_col_ + new_num_row;
    basis.nonbasic_flag_.resize(new_num_tot as usize, 0);
    basis.nonbasic_move_.resize(new_num_tot as usize, 0);
    basis.basic_index_.resize(new_num_row as usize, 0);
    // Make the new rows basic
    for i_row in lp.num_row_..new_num_row {
        basis.nonbasic_flag_[(lp.num_col_ + i_row) as usize] = kNonbasicFlagFalse;
        basis.nonbasic_move_[(lp.num_col_ + i_row) as usize] = 0;
        basis.basic_index_[i_row as usize] = lp.num_col_ + i_row;
    }
}

/// HiGHS: HSimplex.cpp::getUnscaledInfeasibilities (lines 157-254, 98 lines)
pub fn get_unscaled_infeasibilities(
    options: &HighsOptions,
    scale: &HighsScale,
    basis: &SimplexBasis,
    info: &HighsSimplexInfo,
    highs_info: &mut HighsInfo,
) {
    let primal_feasibility_tolerance = options.primal_feasibility_tolerance;
    let dual_feasibility_tolerance = options.dual_feasibility_tolerance;

    // Zero the counts of unscaled primal and dual infeasibilities
    highs_info.num_primal_infeasibilities = 0;
    highs_info.max_primal_infeasibility = 0.0;
    highs_info.sum_primal_infeasibilities = 0.0;
    highs_info.num_dual_infeasibilities = 0;
    highs_info.max_dual_infeasibility = 0.0;
    highs_info.sum_dual_infeasibilities = 0.0;

    debug_assert!((scale.col.len() as HighsInt) == scale.num_col);
    debug_assert!((scale.row.len() as HighsInt) == scale.num_row);
    for i_var in 0..(scale.num_col + scale.num_row) {
        let vi = i_var as usize;
        // Look at the dual infeasibilities of nonbasic variables
        if basis.nonbasic_flag_[vi] == kNonbasicFlagFalse {
            continue;
        }
        // No dual infeasibility for fixed rows and columns
        if info.work_lower_[vi] == info.work_upper_[vi] {
            continue;
        }
        let col = i_var < scale.num_col;
        let scale_mu = if col {
            let i_col = i_var;
            debug_assert!((scale.col.len() as HighsInt) > i_col);
            1.0 / (scale.col[i_col as usize] / scale.cost)
        } else {
            let i_row = i_var - scale.num_col;
            debug_assert!((scale.row.len() as HighsInt) > i_row);
            scale.row[i_row as usize] * scale.cost
        };
        let dual = info.work_dual_[vi];
        let lower = info.work_lower_[vi];
        let upper = info.work_upper_[vi];
        let unscaled_dual = dual * scale_mu;

        let dual_infeasibility = if highs_is_infinity(-lower) && highs_is_infinity(upper) {
            // Free: any nonzero dual value is infeasible
            unscaled_dual.abs()
        } else {
            // Not fixed: any dual infeasibility is given by value signed by
            // nonbasicMove. This assumes that nonbasicMove=0 for fixed
            // variables
            (-(basis.nonbasic_move_[vi] as HighsInt)) as f64 * unscaled_dual
        };
        if dual_infeasibility > 0.0 {
            if dual_infeasibility >= dual_feasibility_tolerance {
                highs_info.num_dual_infeasibilities += 1;
            }
            highs_info.max_dual_infeasibility =
                cpp_max(dual_infeasibility, highs_info.max_dual_infeasibility);
            highs_info.sum_dual_infeasibilities += dual_infeasibility;
        }
    }
    // Look at the primal infeasibilities of basic variables
    for ix in 0..scale.num_row {
        let i_var = basis.basic_index_[ix as usize];
        let col = i_var < scale.num_col;
        let scale_mu = if col {
            let i_col = i_var;
            scale.col[i_col as usize]
        } else {
            let i_row = i_var - scale.num_col;
            1.0 / scale.row[i_row as usize]
        };
        let unscaled_lower = info.base_lower_[ix as usize] * scale_mu;
        let unscaled_value = info.base_value_[ix as usize] * scale_mu;
        let unscaled_upper = info.base_upper_[ix as usize] * scale_mu;
        // @primal_infeasibility calculation
        let mut primal_infeasibility = 0.0;
        if unscaled_value < unscaled_lower - primal_feasibility_tolerance {
            primal_infeasibility = unscaled_lower - unscaled_value;
        } else if unscaled_value > unscaled_upper + primal_feasibility_tolerance {
            primal_infeasibility = unscaled_value - unscaled_upper;
        }
        if primal_infeasibility > 0.0 {
            highs_info.num_primal_infeasibilities += 1;
            highs_info.max_primal_infeasibility =
                cpp_max(primal_infeasibility, highs_info.max_primal_infeasibility);
            highs_info.sum_primal_infeasibilities += primal_infeasibility;
        }
    }
    set_solution_status(highs_info);
}

/// HiGHS: HSimplex.cpp::setSolutionStatus (lines 256-271, 16 lines)
pub fn set_solution_status(highs_info: &mut HighsInfo) {
    if highs_info.num_primal_infeasibilities < 0 {
        highs_info.primal_solution_status = kSolutionStatusNone;
    } else if highs_info.num_primal_infeasibilities == 0 {
        highs_info.primal_solution_status = kSolutionStatusFeasible;
    } else {
        highs_info.primal_solution_status = kSolutionStatusInfeasible;
    }
    if highs_info.num_dual_infeasibilities < 0 {
        highs_info.dual_solution_status = kSolutionStatusNone;
    } else if highs_info.num_dual_infeasibilities == 0 {
        highs_info.dual_solution_status = kSolutionStatusFeasible;
    } else {
        highs_info.dual_solution_status = kSolutionStatusInfeasible;
    }
}
// SCALING

/// HiGHS: HSimplex.cpp::scaleSimplexCost (lines 274-313, 40 lines)
pub fn scale_simplex_cost(options: &HighsOptions, lp: &mut HighsLp, cost_scale: &mut f64) {
    // Scale the costs by no less than minAlwCostScale
    let max_allowed_cost_scale = 2.0_f64.powf(options.allowed_cost_scale_factor as f64);
    let mut max_nonzero_cost: f64 = 0.0;
    for i_col in 0..lp.num_col_ {
        let ci = i_col as usize;
        if lp.col_cost_[ci] != 0.0 {
            max_nonzero_cost = cpp_max(lp.col_cost_[ci].abs(), max_nonzero_cost);
        }
    }
    // Scaling the costs up effectively increases the dual tolerance to
    // which the problem is solved - so, if the max cost is small the
    // scaling factor pushes it up by a power of 2 so it's close to 1
    // Scaling the costs down effectively decreases the dual tolerance
    // to which the problem is solved - so this can't be done too much
    *cost_scale = 1.0;
    let ln2 = 2.0_f64.ln();
    // Scale if the max cost is positive and outside the range [1/16, 16]
    if (max_nonzero_cost > 0.0) && ((max_nonzero_cost < (1.0 / 16.0)) || (max_nonzero_cost > 16.0))
    {
        *cost_scale = max_nonzero_cost;
        *cost_scale = 2.0_f64.powf((cost_scale.ln() / ln2 + 0.5).floor());
        *cost_scale = cpp_min(*cost_scale, max_allowed_cost_scale);
    }
    if *cost_scale == 1.0 {
        if highs_log_user_enabled(&options.log_options) {
            highs_log_user(
                &options.log_options,
                HighsLogType::kInfo,
                &format!(
                    "LP cost vector not scaled down: max cost is {}\n",
                    max_nonzero_cost
                ),
            );
        }
        return;
    }
    // Scale the costs (and record of max_nonzero_cost) by cost_scale, being at
    // most max_allowed_cost_scale
    for i_col in 0..lp.num_col_ {
        lp.col_cost_[i_col as usize] /= *cost_scale;
    }
    max_nonzero_cost /= *cost_scale;
    if highs_log_user_enabled(&options.log_options) {
        highs_log_user(
            &options.log_options,
            HighsLogType::kInfo,
            &format!(
                "LP cost vector scaled down by {}: max cost is {}\n",
                *cost_scale, max_nonzero_cost
            ),
        );
    }
}

/// HiGHS: HSimplex.cpp::unscaleSimplexCost (lines 315-318, 4 lines)
pub fn unscale_simplex_cost(lp: &mut HighsLp, cost_scale: f64) {
    for i_col in 0..lp.num_col_ {
        lp.col_cost_[i_col as usize] *= cost_scale;
    }
}

/// HiGHS: HSimplex.cpp::isBasisRightSize (lines 320-330, 11 lines)
pub fn is_basis_right_size(lp: &HighsLp, basis: &SimplexBasis) -> bool {
    let mut right_size = true;
    right_size =
        (basis.nonbasic_flag_.len() as HighsInt) == lp.num_col_ + lp.num_row_ && right_size;
    right_size =
        (basis.nonbasic_move_.len() as HighsInt) == lp.num_col_ + lp.num_row_ && right_size;
    right_size = (basis.basic_index_.len() as HighsInt) == lp.num_row_ && right_size;
    right_size
}
