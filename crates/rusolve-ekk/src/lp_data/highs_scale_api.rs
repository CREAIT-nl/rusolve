//! HiGHS source: lp_data/Highs.cpp (scaleCol, scaleRow),
//! lp_data/HighsInterface.cpp (scaleColInterface, scaleRowInterface),
//! lp_data/HighsLpUtils.cpp (applyScalingToLpCol, applyScalingToLpRow)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! User scaling of one column / row of the incumbent model.
//! The interfaces also scale the solution, dual values and basis status
//! (a negative scale flips nonbasic bounds), and keep the simplex LP valid
//! only when it is unscaled (`ekk_instance_.lp_` is updated in place).
//! `HighsSparseMatrix::scaleCol` / `scaleRow` are util::sparse_matrix.

use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::scaleCol (lines 3581-3589)
    pub fn scale_col(&mut self, col: HighsInt, scale_value: f64) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;
        self.clear_derived_model_properties();
        let call_status = self.scale_col_interface(col, scale_value);
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "scaleCol",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::scaleRow (lines 3591-3599)
    pub fn scale_row(&mut self, row: HighsInt, scale_value: f64) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;
        self.clear_derived_model_properties();
        let call_status = self.scale_row_interface(row, scale_value);
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "scaleRow",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: HighsInterface.cpp::Highs::scaleColInterface (lines 1095-1138)
    pub(crate) fn scale_col_interface(&mut self, col: HighsInt, scale_value: f64) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;

        // Ensure that the LP is column-wise
        self.model_lp_.ensure_colwise();
        if col < 0 {
            return HighsStatus::kError;
        }
        if col >= self.model_lp_.num_col_ {
            return HighsStatus::kError;
        }
        if scale_value == 0.0 {
            return HighsStatus::kError;
        }

        let call_status = apply_scaling_to_lp_col(&mut self.model_lp_, col, scale_value);
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "applyScalingToLpCol",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }

        if scale_value < 0.0 && self.basis_.valid {
            // Negative, so flip any nonbasic status
            if self.basis_.col_status[col as usize] == HighsBasisStatus::kLower {
                self.basis_.col_status[col as usize] = HighsBasisStatus::kUpper;
            } else if self.basis_.col_status[col as usize] == HighsBasisStatus::kUpper {
                self.basis_.col_status[col as usize] = HighsBasisStatus::kLower;
            }
        }
        if self.ekk_instance_.status_.initialised_for_solve {
            if scale_value < 0.0 && self.ekk_instance_.status_.has_basis {
                // Negative, so flip any nonbasic status
                if self.ekk_instance_.basis_.nonbasic_move_[col as usize] == kNonbasicMoveUp {
                    self.ekk_instance_.basis_.nonbasic_move_[col as usize] = kNonbasicMoveDn;
                } else if self.ekk_instance_.basis_.nonbasic_move_[col as usize] == kNonbasicMoveDn
                {
                    self.ekk_instance_.basis_.nonbasic_move_[col as usize] = kNonbasicMoveUp;
                }
            }
        }
        // Deduce the consequences of a scaled column
        self.invalidate_model_status_solution_and_info();

        // Determine any implications for simplex data
        self.ekk_instance_.update_status(LpAction::kScaledCol);
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::scaleRowInterface (lines 1140-1185)
    pub(crate) fn scale_row_interface(&mut self, row: HighsInt, scale_value: f64) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;

        // Ensure that the LP is column-wise
        self.model_lp_.ensure_colwise();

        if row < 0 {
            return HighsStatus::kError;
        }
        if row >= self.model_lp_.num_row_ {
            return HighsStatus::kError;
        }
        if scale_value == 0.0 {
            return HighsStatus::kError;
        }

        let call_status = apply_scaling_to_lp_row(&mut self.model_lp_, row, scale_value);
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "applyScalingToLpRow",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }

        if scale_value < 0.0 && self.basis_.valid {
            // Negative, so flip any nonbasic status
            if self.basis_.row_status[row as usize] == HighsBasisStatus::kLower {
                self.basis_.row_status[row as usize] = HighsBasisStatus::kUpper;
            } else if self.basis_.row_status[row as usize] == HighsBasisStatus::kUpper {
                self.basis_.row_status[row as usize] = HighsBasisStatus::kLower;
            }
        }
        if self.ekk_instance_.status_.initialised_for_solve {
            if scale_value < 0.0 && self.ekk_instance_.status_.has_basis {
                // Negative, so flip any nonbasic status
                let var = self.model_lp_.num_col_ + row;
                if self.ekk_instance_.basis_.nonbasic_move_[var as usize] == kNonbasicMoveUp {
                    self.ekk_instance_.basis_.nonbasic_move_[var as usize] = kNonbasicMoveDn;
                } else if self.ekk_instance_.basis_.nonbasic_move_[var as usize] == kNonbasicMoveDn
                {
                    self.ekk_instance_.basis_.nonbasic_move_[var as usize] = kNonbasicMoveUp;
                }
            }
        }
        // Deduce the consequences of a scaled row
        self.invalidate_model_status_solution_and_info();

        // Determine any implications for simplex data
        self.ekk_instance_.update_status(LpAction::kScaledRow);
        HighsStatus::kOk
    }
}

/// HiGHS: HighsLpUtils.cpp::applyScalingToLpCol (lines 1534-1551)
pub fn apply_scaling_to_lp_col(lp: &mut HighsLp, col: HighsInt, col_scale: f64) -> HighsStatus {
    if col < 0 {
        return HighsStatus::kError;
    }
    if col >= lp.num_col_ {
        return HighsStatus::kError;
    }
    if col_scale == 0.0 {
        return HighsStatus::kError;
    }

    lp.a_matrix_.scale_col(col, col_scale);
    lp.col_cost_[col as usize] *= col_scale;
    if col_scale > 0.0 {
        lp.col_lower_[col as usize] /= col_scale;
        lp.col_upper_[col as usize] /= col_scale;
    } else {
        let new_upper = lp.col_lower_[col as usize] / col_scale;
        lp.col_lower_[col as usize] = lp.col_upper_[col as usize] / col_scale;
        lp.col_upper_[col as usize] = new_upper;
    }
    HighsStatus::kOk
}

/// HiGHS: HighsLpUtils.cpp::applyScalingToLpRow (lines 1553-1569)
pub fn apply_scaling_to_lp_row(lp: &mut HighsLp, row: HighsInt, row_scale: f64) -> HighsStatus {
    if row < 0 {
        return HighsStatus::kError;
    }
    if row >= lp.num_row_ {
        return HighsStatus::kError;
    }
    if row_scale == 0.0 {
        return HighsStatus::kError;
    }

    lp.a_matrix_.scale_row(row, row_scale);
    if row_scale > 0.0 {
        lp.row_lower_[row as usize] *= row_scale;
        lp.row_upper_[row as usize] *= row_scale;
    } else {
        let new_upper = lp.row_lower_[row as usize] * row_scale;
        lp.row_lower_[row as usize] = lp.row_upper_[row as usize] * row_scale;
        lp.row_upper_[row as usize] = new_upper;
    }
    HighsStatus::kOk
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// x0 in [1, 4], x1 in [-inf, 2]; row 2 x0 + x1 in [0, 6].
    fn lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, -1.0];
        lp.col_lower_ = vec![1.0, -kHighsInf];
        lp.col_upper_ = vec![4.0, 2.0];
        lp.row_lower_ = vec![0.0];
        lp.row_upper_ = vec![6.0];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![2.0, 1.0];
        lp
    }

    #[test]
    fn negative_column_scale_swaps_bounds() {
        let mut lp = lp();
        assert_eq!(apply_scaling_to_lp_col(&mut lp, 0, -2.0), HighsStatus::kOk);
        assert_eq!(lp.col_cost_[0], -2.0);
        assert_eq!((lp.col_lower_[0], lp.col_upper_[0]), (-2.0, -0.5));
        assert_eq!(lp.a_matrix_.value_[0], -4.0);
        assert_eq!(
            apply_scaling_to_lp_col(&mut lp, 0, 0.0),
            HighsStatus::kError
        );
        assert_eq!(
            apply_scaling_to_lp_col(&mut lp, 2, 1.0),
            HighsStatus::kError
        );
    }

    #[test]
    fn negative_row_scale_swaps_bounds() {
        let mut lp = lp();
        assert_eq!(apply_scaling_to_lp_row(&mut lp, 0, -0.5), HighsStatus::kOk);
        assert_eq!((lp.row_lower_[0], lp.row_upper_[0]), (-3.0, -0.0));
        assert_eq!(lp.a_matrix_.value_, vec![-1.0, -0.5]);
    }

    #[test]
    fn highs_scale_col_and_row_check_indices() {
        let mut highs = Highs::default();
        highs.model_lp_ = lp();
        assert_eq!(highs.scale_col(-1, 2.0), HighsStatus::kError);
        assert_eq!(highs.scale_row(1, 2.0), HighsStatus::kError);
        assert_eq!(highs.scale_col(1, 2.0), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.col_upper_[1], 1.0);
        assert_eq!(highs.scale_row(0, 2.0), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.row_upper_[0], 12.0);
    }
}
