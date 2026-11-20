//! HiGHS source: lp_data/Highs.cpp (setSolution, both overloads; setBasis()),
//! lp_data/HighsInterface.cpp (basisForSolution), lp_data/HighsLpUtils.cpp
//! (calculateColDualsQuad), lp_data/HighsSolution.cpp (HighsSolution::hasUndefined)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! A user-supplied (warm start) solution. `set_solution` installs
//! full primal and/or dual values (row values / column duals recomputed in
//! quad precision); `set_solution_sparse` builds a partial primal solution
//! with `kHighsUndefined` for unset columns, which a MIP completes in
//! `complete_solution_from_discrete_assignment`. For an LP with a
//! basis-using solver, `calledOptimizeModel` turns a valid solution without a
//! valid basis into one with `basis_for_solution` (wired in highs_run.rs).
//!
//! `calculateRowValuesQuad` is already ported
//! (`presolve_support::calculate_row_values_quad`).

use crate::highs_log_dev;
use crate::highs_log_user;
use crate::lp_data::highs_solution::{assess_user_col_values, assign_user_col_values};
use crate::lp_data::presolve_support::calculate_row_values_quad;
use crate::prelude::*;
use crate::util::highs_cdouble::HighsCDouble;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::setSolution (lines 2492-2548), the overload
    /// `setSolution(const HighsSolution& solution)`.
    pub fn set_solution(&mut self, solution: &HighsSolution) -> HighsStatus {
        // HiGHS 1.16 #3556 (h16_user_solution_assess): Highs.cpp lines 2708-2775 (1.16.0)
        if self.options_.h16_user_solution_assess {
            return self.set_solution_h16(solution);
        }
        let mut return_status = HighsStatus::kOk;
        // Determine whether a new solution will be defined. If so,
        // the old solution and any basis are cleared
        let new_primal_solution = self.model_lp_.num_col_ > 0
            && solution.col_value.len() >= self.model_lp_.num_col_ as usize;
        let new_dual_solution = self.model_lp_.num_row_ > 0
            && solution.row_dual.len() >= self.model_lp_.num_row_ as usize;
        let new_solution = new_primal_solution || new_dual_solution;

        if new_solution {
            self.invalidate_solver_data();
        } else {
            // Solution is rejected, so give a logging message and error
            // return
            highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "setSolution: User solution is rejected due to mismatch between \
                 size of col_value and row_dual vectors ({}, {}) and number \
                 of columns and rows in the model ({}, {})\n",
                solution.col_value.len() as i32,
                solution.row_dual.len() as i32,
                self.model_lp_.num_col_,
                self.model_lp_.num_row_
            );
            return_status = HighsStatus::kError;
        }

        if new_primal_solution {
            self.solution_.col_value = solution.col_value.clone();
            if self.model_lp_.num_row_ > 0 {
                // Worth computing the row values
                self.solution_
                    .row_value
                    .resize(self.model_lp_.num_row_ as usize, 0.0);
                // Matrix must be column-wise
                self.model_lp_.a_matrix_.ensure_colwise();
                let mut sol = std::mem::take(&mut self.solution_);
                let call_status = calculate_row_values_quad(&self.model_lp_, &mut sol, -1);
                self.solution_ = sol;
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "calculateRowValuesQuad",
                );
                if return_status == HighsStatus::kError {
                    return return_status;
                }
            }
            self.solution_.value_valid = true;
        }
        if new_dual_solution {
            self.solution_.row_dual = solution.row_dual.clone();
            if self.model_lp_.num_col_ > 0 {
                // Worth computing the column duals
                self.solution_
                    .col_dual
                    .resize(self.model_lp_.num_col_ as usize, 0.0);
                // Matrix must be column-wise
                self.model_lp_.a_matrix_.ensure_colwise();
                let mut sol = std::mem::take(&mut self.solution_);
                let call_status = calculate_col_duals_quad(&self.model_lp_, &mut sol);
                self.solution_ = sol;
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "calculateColDuals",
                );
                if return_status == HighsStatus::kError {
                    return return_status;
                }
            }
            self.solution_.dual_valid = true;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::setSolution (lines 2576-2625), the overload
    /// `setSolution(const HighsInt num_entries, const HighsInt* index, const
    /// double* value)`. `std::vector<bool> is_set` -> `Vec<bool>`.
    pub fn set_solution_sparse(
        &mut self,
        num_entries: HighsInt,
        index: &[HighsInt],
        value: &[f64],
    ) -> HighsStatus {
        // HiGHS 1.16 #3556 (h16_user_solution_assess): Highs.cpp lines 2804-2824 (1.16.0)
        if self.options_.h16_user_solution_assess {
            return self.set_solution_sparse_h16(num_entries, index, value);
        }
        let mut return_status = HighsStatus::kOk;
        if self.model_lp_.num_col_ == 0 {
            return return_status;
        }
        // Warn about duplicates in index
        let mut num_duplicates: HighsInt = 0;
        let mut is_set: Vec<bool>;
        is_set = vec![false; self.model_lp_.num_col_ as usize];
        for i_x in 0..num_entries {
            let i_col = index[i_x as usize];
            if i_col < 0 || i_col >= self.model_lp_.num_col_ {
                highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "setSolution: User solution index {} has value {} out of \
                     range [0, {})\n",
                    i_x,
                    i_col,
                    self.model_lp_.num_col_
                );
                return HighsStatus::kError;
            } else if value[i_x as usize]
                < self.model_lp_.col_lower_[i_col as usize]
                    - self.options_.primal_feasibility_tolerance
                || self.model_lp_.col_upper_[i_col as usize]
                    + self.options_.primal_feasibility_tolerance
                    < value[i_x as usize]
            {
                highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "setSolution: User solution value {} of {} is infeasible \
                     for bounds [{}, {}]\n",
                    i_x,
                    crate::util::c_printf::c_fmt_g(value[i_x as usize], 6),
                    crate::util::c_printf::c_fmt_g(self.model_lp_.col_lower_[i_col as usize], 6),
                    crate::util::c_printf::c_fmt_g(self.model_lp_.col_upper_[i_col as usize], 6)
                );
                return HighsStatus::kError;
            }
            if is_set[i_col as usize] {
                num_duplicates += 1;
            }
            is_set[i_col as usize] = true;
        }
        if num_duplicates > 0 {
            highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "setSolution: User set of indices has {} duplicate{}: last \
                 value used\n",
                num_duplicates,
                if num_duplicates == 1 { "" } else { "s" }
            );
            return_status = HighsStatus::kWarning;
        }

        // Clear the solution, indicate the values not determined by the
        // user, and insert the values determined by the user
        let mut new_solution = HighsSolution::default();
        new_solution.col_value = vec![kHighsUndefined; self.model_lp_.num_col_ as usize];
        for i_x in 0..num_entries {
            let i_col = index[i_x as usize];
            new_solution.col_value[i_col as usize] = value[i_x as usize];
        }
        let call_status = self.set_solution(&new_solution);
        interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "setSolution",
        )
    }

    /// HiGHS 1.16 (#3556): Highs.cpp::Highs::setSolution (lines 2708-2775), the
    /// overload `setSolution(const HighsSolution& solution)`. Called only when
    /// `h16_user_solution_assess` is on. `HighsLp& lp` is `self.model_lp_`.
    fn set_solution_h16(&mut self, solution: &HighsSolution) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;
        // Determine whether a new solution will be defined. If so,
        // the old solution and any basis are cleared
        let new_primal_solution = self.model_lp_.num_col_ > 0
            && solution.col_value.len() >= self.model_lp_.num_col_ as usize;
        let new_dual_solution = self.model_lp_.num_row_ > 0
            && solution.row_dual.len() >= self.model_lp_.num_row_ as usize;
        let new_solution = new_primal_solution || new_dual_solution;

        if new_solution {
            self.invalidate_solver_data();
        } else {
            // Solution is rejected, so give a logging message and error
            // return
            highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "setSolution: User solution is rejected due to mismatch between \
                 size of col_value and row_dual vectors ({}, {}) and number \
                 of columns and rows in the model ({}, {})\n",
                solution.col_value.len() as i32,
                solution.row_dual.len() as i32,
                self.model_lp_.num_col_,
                self.model_lp_.num_row_
            );
            return HighsStatus::kError;
        }

        let mut modifiable_solution = solution.clone();
        if new_primal_solution {
            // There should be no error return from assessUserColValues, since
            // this corresponds to index illegalities, but use it to warn about
            // infeasibilities in value
            return_status = assess_user_col_values(
                &self.model_lp_,
                &self.options_,
                self.model_lp_.num_col_,
                None,
                &solution.col_value,
                true,
            );
            debug_assert!(return_status != HighsStatus::kError);
            if self.model_lp_.num_row_ > 0 {
                // Worth computing the row values
                modifiable_solution
                    .row_value
                    .resize(self.model_lp_.num_row_ as usize, 0.0);
                // Matrix must be column-wise
                self.model_lp_.a_matrix_.ensure_colwise();
                let call_status =
                    calculate_row_values_quad(&self.model_lp_, &mut modifiable_solution, -1);
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "calculateRowValuesQuad",
                );
                if return_status == HighsStatus::kError {
                    return return_status;
                }
            }
            self.solution_.col_value = modifiable_solution.col_value.clone();
            self.solution_.row_value = modifiable_solution.row_value.clone();
            self.solution_.value_valid = true;
        }
        if new_dual_solution {
            if self.model_lp_.num_col_ > 0 {
                // Worth computing the column duals
                modifiable_solution
                    .col_dual
                    .resize(self.model_lp_.num_col_ as usize, 0.0);
                // Matrix must be column-wise
                self.model_lp_.a_matrix_.ensure_colwise();
                let call_status =
                    calculate_col_duals_quad(&self.model_lp_, &mut modifiable_solution);
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "calculateColDuals",
                );
                if return_status == HighsStatus::kError {
                    return return_status;
                }
            }
            self.solution_.col_dual = modifiable_solution.col_dual.clone();
            self.solution_.row_dual = modifiable_solution.row_dual.clone();
            self.solution_.dual_valid = true;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS 1.16 (#3556): Highs.cpp::Highs::setSolution (lines 2804-2824), the
    /// overload `setSolution(const HighsInt num_entries, const HighsInt* index,
    /// const double* value)`. Called only when `h16_user_solution_assess` is
    /// on. The Rust API always passes an index slice (C++ `index != nullptr`).
    fn set_solution_sparse_h16(
        &mut self,
        num_entries: HighsInt,
        index: &[HighsInt],
        value: &[f64],
    ) -> HighsStatus {
        if self.model_lp_.num_col_ == 0 {
            return HighsStatus::kOk;
        }
        // Check for illegalities, and warn about duplicates in index, but
        // don't check for feasibility of values, as this will be done when
        // setSolution is called below.
        let perform_feasibility_check = false;
        let return_status = assess_user_col_values(
            &self.model_lp_,
            &self.options_,
            num_entries,
            Some(index),
            value,
            perform_feasibility_check,
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }

        // Clear the solution, indicate the values not determined by the
        // user, and insert the values determined by the user
        let mut new_solution = HighsSolution::default();
        assign_user_col_values(
            self.model_lp_.num_col_,
            &mut new_solution.col_value,
            num_entries,
            Some(index),
            value,
        );
        let call_status = self.set_solution(&new_solution);
        interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "setSolution",
        )
    }

    /// HiGHS: Highs.cpp::Highs::setBasis (lines 2807-2817), the no-argument
    /// overload `setBasis()`: `invalidateBasis()` then `newHighsBasis()`.
    pub fn set_basis_invalid(&mut self) -> HighsStatus {
        // Invalidate the basis for HiGHS
        //
        // Don't set to logical basis since that causes presolve to be
        // skipped
        self.invalidate_basis();
        // Follow implications of a new HiGHS basis
        self.new_highs_basis();
        // Can't use returnFromHighs since...
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::basisForSolution (lines 352-393).
    /// The local `HighsBasis basis` is passed to `set_basis(&basis, "")`.
    pub fn basis_for_solution(&mut self) -> HighsStatus {
        debug_assert!(!self.model_lp_.is_mip() || self.options_.solve_relaxation);
        debug_assert!(self.solution_.value_valid);
        self.invalidate_basis();
        let mut num_basic: HighsInt = 0;
        let mut basis = HighsBasis::default();
        for i_col in 0..self.model_lp_.num_col_ {
            if (self.model_lp_.col_lower_[i_col as usize]
                - self.solution_.col_value[i_col as usize])
                .abs()
                <= self.options_.primal_feasibility_tolerance
            {
                basis.col_status.push(HighsBasisStatus::kLower);
            } else if (self.model_lp_.col_upper_[i_col as usize]
                - self.solution_.col_value[i_col as usize])
                .abs()
                <= self.options_.primal_feasibility_tolerance
            {
                basis.col_status.push(HighsBasisStatus::kUpper);
            } else {
                num_basic += 1;
                basis.col_status.push(HighsBasisStatus::kBasic);
            }
        }
        let num_basic_col = num_basic;
        for i_row in 0..self.model_lp_.num_row_ {
            if (self.model_lp_.row_lower_[i_row as usize]
                - self.solution_.row_value[i_row as usize])
                .abs()
                <= self.options_.primal_feasibility_tolerance
            {
                basis.row_status.push(HighsBasisStatus::kLower);
            } else if (self.model_lp_.row_upper_[i_row as usize]
                - self.solution_.row_value[i_row as usize])
                .abs()
                <= self.options_.primal_feasibility_tolerance
            {
                basis.row_status.push(HighsBasisStatus::kUpper);
            } else {
                num_basic += 1;
                basis.row_status.push(HighsBasisStatus::kBasic);
            }
        }
        let num_basic_row = num_basic - num_basic_col;
        debug_assert!(basis.col_status.len() as i32 == self.model_lp_.num_col_);
        debug_assert!(basis.row_status.len() as i32 == self.model_lp_.num_row_);
        highs_log_dev!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "LP has {} rows and solution yields {} possible basic variables \
             ({} / {}; {} / {})\n",
            self.model_lp_.num_row_,
            num_basic,
            num_basic_col,
            self.model_lp_.num_col_,
            num_basic_row,
            self.model_lp_.num_row_
        );
        self.set_basis(&basis, "")
    }
}

/// HiGHS: HighsLpUtils.cpp::calculateColDualsQuad (lines 2980-3009).
/// `HighsCDouble` accumulators (util::hcdouble).
pub fn calculate_col_duals_quad(lp: &HighsLp, solution: &mut HighsSolution) -> HighsStatus {
    let correct_size = solution.row_dual.len() as i32 == lp.num_row_;
    let is_colwise = lp.a_matrix_.is_colwise();
    let data_error = !correct_size || !is_colwise;
    debug_assert!(!data_error);
    if data_error {
        return HighsStatus::kError;
    }

    let mut col_dual_quad: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); lp.num_col_ as usize];

    for col in 0..lp.num_col_ {
        let mut i = lp.a_matrix_.start_[col as usize];
        while i < lp.a_matrix_.start_[(col + 1) as usize] {
            let row = lp.a_matrix_.index_[i as usize];
            debug_assert!(row >= 0);
            debug_assert!(row < lp.num_row_);
            col_dual_quad[col as usize] += HighsCDouble::from(solution.row_dual[row as usize])
                * lp.a_matrix_.value_[i as usize];
            i += 1;
        }
        col_dual_quad[col as usize] += lp.col_cost_[col as usize];
    }

    // assign quad values to double vector
    solution.col_dual.resize(lp.num_col_ as usize, 0.0);
    for i in 0..lp.num_col_ as usize {
        solution.col_dual[i] = f64::from(col_dual_quad[i]);
    }

    HighsStatus::kOk
}

impl HighsSolution {
    /// HiGHS: HighsSolution.cpp::HighsSolution::hasUndefined (lines 2148-2152)
    pub fn has_undefined(&self) -> bool {
        for value in &self.col_value {
            if *value == kHighsUndefined {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// min x0 + x1, x0 + 2 x1 in [1, 8], x0 in [0, 4], x1 in [0, 3].
    fn lp_highs() -> Highs {
        let mut highs = Highs::default();
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![4.0, 3.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![8.0];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 2.0];
        highs.model_lp_ = lp;
        highs
    }

    #[test]
    fn full_primal_solution_computes_row_values() {
        let mut highs = lp_highs();
        let mut solution = HighsSolution::default();
        solution.col_value = vec![1.0, 1.5];
        assert_eq!(highs.set_solution(&solution), HighsStatus::kOk);
        assert!(highs.solution_.value_valid);
        assert!(!highs.solution_.dual_valid);
        assert_eq!(highs.solution_.row_value, vec![4.0]);
    }

    #[test]
    fn dual_solution_computes_col_duals() {
        let mut highs = lp_highs();
        let mut solution = HighsSolution::default();
        solution.row_dual = vec![0.5];
        assert_eq!(highs.set_solution(&solution), HighsStatus::kOk);
        assert!(highs.solution_.dual_valid);
        // col_dual = c + A^T y (calculateColDualsQuad adds col_cost_ after the row-dual sum)
        assert_eq!(highs.solution_.col_dual, vec![1.5, 2.0]);
    }

    #[test]
    fn mismatched_sizes_are_rejected() {
        let mut highs = lp_highs();
        let mut solution = HighsSolution::default();
        solution.col_value = vec![1.0];
        assert_eq!(highs.set_solution(&solution), HighsStatus::kError);
        assert!(!highs.solution_.value_valid);
    }

    #[test]
    fn sparse_solution_marks_unset_columns_undefined() {
        let mut highs = lp_highs();
        assert_eq!(highs.set_solution_sparse(1, &[1], &[2.0]), HighsStatus::kOk);
        assert_eq!(highs.solution_.col_value[0], kHighsUndefined);
        assert!(highs.solution_.has_undefined());
        // Duplicates: last value used, with a warning
        assert_eq!(
            highs.set_solution_sparse(2, &[0, 0], &[1.0, 2.0]),
            HighsStatus::kWarning
        );
        assert_eq!(highs.solution_.col_value[0], 2.0);
        // Infeasible for the bounds
        assert_eq!(
            highs.set_solution_sparse(1, &[1], &[9.0]),
            HighsStatus::kError
        );
    }

    #[test]
    fn basis_for_solution_classifies_bounds() {
        let mut highs = lp_highs();
        let mut solution = HighsSolution::default();
        solution.col_value = vec![0.0, 3.0];
        assert_eq!(highs.set_solution(&solution), HighsStatus::kOk);
        assert_eq!(highs.basis_for_solution(), HighsStatus::kOk);
        assert!(highs.basis_.valid);
        assert_eq!(
            highs.basis_.col_status,
            vec![HighsBasisStatus::kLower, HighsBasisStatus::kUpper]
        );
        // Row value 6 is strictly inside [1, 8]
        assert_eq!(highs.basis_.row_status, vec![HighsBasisStatus::kBasic]);
    }

    // ---- HiGHS 1.16 (#3556, h16_user_solution_assess) ----

    fn lp_highs_on() -> Highs {
        let mut highs = lp_highs();
        highs.options_.h16_user_solution_assess = true;
        highs
    }

    /// 1.16 Highs.cpp 2814-2817 + 2741-2742: the sparse overload no longer
    /// rejects a bound-infeasible value; the dense overload it calls warns
    /// (assessUserColValues, HighsSolution.cpp 2340-2356) and stores it.
    #[test]
    fn h16_usersol_on_sparse_bound_infeasible_value_is_stored_with_warning() {
        let mut highs = lp_highs_on();
        assert_eq!(
            highs.set_solution_sparse(1, &[1], &[9.0]),
            HighsStatus::kWarning
        );
        assert!(highs.solution_.value_valid);
        assert_eq!(highs.solution_.col_value[0], kHighsUndefined);
        assert_eq!(highs.solution_.col_value[1], 9.0);
    }

    /// 1.15.1 Highs.cpp 2600-2611: the same value is rejected.
    #[test]
    fn h16_usersol_off_sparse_bound_infeasible_value_is_error() {
        let mut highs = lp_highs();
        assert_eq!(
            highs.set_solution_sparse(1, &[1], &[9.0]),
            HighsStatus::kError
        );
        assert!(!highs.solution_.value_valid);
    }

    /// Dense setSolution with x0 = 5 > 4: 1.16 kWarning (assessUserColValues),
    /// values and row values still set; 1.15.1 kOk.
    #[test]
    fn h16_usersol_dense_infeasible_value_warns_only_when_on() {
        let mut solution = HighsSolution::default();
        solution.col_value = vec![5.0, 1.0];
        let mut highs = lp_highs_on();
        assert_eq!(highs.set_solution(&solution), HighsStatus::kWarning);
        assert!(highs.solution_.value_valid);
        assert_eq!(highs.solution_.col_value, vec![5.0, 1.0]);
        assert_eq!(highs.solution_.row_value, vec![7.0]);
        let mut highs = lp_highs();
        assert_eq!(highs.set_solution(&solution), HighsStatus::kOk);
        assert_eq!(highs.solution_.row_value, vec![7.0]);
    }

    /// Dense primal and dual on path: same values as the 1.15.1 path.
    #[test]
    fn h16_usersol_on_dense_primal_and_dual_match_off() {
        let mut solution = HighsSolution::default();
        solution.col_value = vec![1.0, 1.5];
        solution.row_dual = vec![0.5];
        let mut on = lp_highs_on();
        let mut off = lp_highs();
        assert_eq!(on.set_solution(&solution), HighsStatus::kOk);
        assert_eq!(off.set_solution(&solution), HighsStatus::kOk);
        assert_eq!(on.solution_.row_value, vec![4.0]);
        assert_eq!(on.solution_.col_dual, vec![1.5, 2.0]);
        assert_eq!(on.solution_.row_value, off.solution_.row_value);
        assert_eq!(on.solution_.col_dual, off.solution_.col_dual);
        assert!(on.solution_.value_valid && on.solution_.dual_valid);
        // Size mismatch: early kError (Highs.cpp 2724-2734)
        let mut bad = HighsSolution::default();
        bad.col_value = vec![1.0];
        assert_eq!(on.set_solution(&bad), HighsStatus::kError);
    }

    /// TestCheckSolution.cpp check-set-mip-solution test7 (3556.diff 271-292),
    /// with index slices for the C++ non-null pointers.
    #[test]
    fn h16_usersol_on_illegal_empty_and_duplicate_entries() {
        let mut highs = lp_highs_on();
        let num_col = highs.model_lp_.num_col_;
        assert_eq!(
            highs.set_solution_sparse(-1, &[num_col], &[]),
            HighsStatus::kError
        );
        // num_entries = 0 is now just a warning. The 1.16 setSolution
        // (Highs.cpp 2814-2823) returns early only on kError, so it goes on
        // to set an all-undefined solution
        assert_eq!(
            highs.set_solution_sparse(0, &[num_col], &[]),
            HighsStatus::kWarning
        );
        assert!(highs.solution_.value_valid);
        assert_eq!(
            highs.solution_.col_value,
            vec![kHighsUndefined, kHighsUndefined]
        );
        // index out of range
        assert_eq!(
            highs.set_solution_sparse(1, &[num_col], &[0.0]),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_solution_sparse(1, &[-1], &[0.0]),
            HighsStatus::kError
        );
        // Check that duplicate values are spotted
        assert_eq!(
            highs.set_solution_sparse(3, &[0, 1, 0], &[0.0, 1.0, 2.0]),
            HighsStatus::kWarning
        );
        assert_eq!(highs.solution_.col_value, vec![2.0, 1.0]);
    }

    /// assessUserColValues tests integrality only when assess_as_mip
    /// (HighsSolution.cpp 2321, 2346-2348; HighsLpUtils.cpp 2622).
    #[test]
    fn h16_usersol_on_integrality_assessed_only_as_mip() {
        let mut highs = lp_highs_on();
        highs.model_lp_.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kContinuous];
        let mut solution = HighsSolution::default();
        solution.col_value = vec![1.5, 1.0];
        assert_eq!(highs.set_solution(&solution), HighsStatus::kWarning);
        highs.options_.solve_relaxation = true;
        assert_eq!(highs.set_solution(&solution), HighsStatus::kOk);
        // assessUserColValues directly: > max_warning infeasible entries
        let options = highs.options_.clone();
        let mut lp = HighsLp::default();
        lp.num_col_ = 12;
        lp.col_lower_ = vec![0.0; 12];
        lp.col_upper_ = vec![1.0; 12];
        let value = vec![2.0; 12];
        assert_eq!(
            crate::lp_data::highs_solution::assess_user_col_values(
                &lp, &options, 12, None, &value, true
            ),
            HighsStatus::kWarning
        );
        assert_eq!(
            crate::lp_data::highs_solution::assess_user_col_values(
                &lp, &options, 11, None, &value, true
            ),
            HighsStatus::kError
        );
    }
}
