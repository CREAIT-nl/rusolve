//! HiGHS source: lp_data/HighsLpUtils.cpp (withoutSemiVariables,
//! activeModifiedUpperBounds)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The MIP side of semi-variables. `Highs::callSolveMip`
//! (highs_run.rs, already wired) solves `without_semi_variables(model_.lp_,
//! solution_, primal_feasibility_tolerance)` when the model has
//! semi-variables: each one gets a binary indicator column and two rows
//! (`x - l*y >= 0`, `x - u*y <= 0`), and any user solution is extended to
//! match. After the solve, `active_modified_upper_bounds` rejects a solution
//! that sits at an upper bound `assess_semi_variables` tightened.

use crate::prelude::*;

/// HiGHS: HighsLpUtils.cpp::withoutSemiVariables (lines 3202-3365).
/// `HighsLp lp = lp_` -> `let mut lp = lp_.clone()`. The column and row
/// names use `format!("semi_binary_{n}")`, `semi_lb_{n}`, `semi_ub_{n}`.
/// HiGHS quirk: the third parameter is named `mip_feasibility_tolerance` but
/// `callSolveMip` passes `options_.primal_feasibility_tolerance`.
pub fn without_semi_variables(
    lp_: &HighsLp,
    solution: &mut HighsSolution,
    mip_feasibility_tolerance: f64,
) -> HighsLp {
    let mut lp = lp_.clone();
    let num_col = lp.num_col_;
    let num_row = lp.num_row_;
    let mut num_semi_variables: HighsInt = 0;
    for i_col in 0..num_col {
        if lp.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
            || lp.integrality_[i_col as usize] == HighsVarType::kSemiInteger
        {
            num_semi_variables += 1;
        }
    }
    debug_assert!(num_semi_variables != 0);
    // Insert spaces for index/value of new coefficients for
    // semi-variables
    // HiGHS: vector<HighsInt>& start = lp.a_matrix_.start_;
    // HiGHS: vector<HighsInt>& index = lp.a_matrix_.index_;
    // HiGHS: vector<double>& value = lp.a_matrix_.value_;
    let num_nz = lp.a_matrix_.start_[num_col as usize];
    let new_num_nz = num_nz + 2 * num_semi_variables;
    let mut new_el = new_num_nz;
    lp.a_matrix_.index_.resize(new_num_nz as usize, 0);
    lp.a_matrix_.value_.resize(new_num_nz as usize, 0.0);
    for i_col in (0..num_col).rev() {
        let from_el = lp.a_matrix_.start_[(i_col + 1) as usize] - 1;
        lp.a_matrix_.start_[(i_col + 1) as usize] = new_el;
        if lp.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
            || lp.integrality_[i_col as usize] == HighsVarType::kSemiInteger
        {
            new_el -= 2;
        }
        let mut i_el = from_el;
        while i_el >= lp.a_matrix_.start_[i_col as usize] {
            new_el -= 1;
            lp.a_matrix_.index_[new_el as usize] = lp.a_matrix_.index_[i_el as usize];
            lp.a_matrix_.value_[new_el as usize] = lp.a_matrix_.value_[i_el as usize];
            i_el -= 1;
        }
    }
    debug_assert!(new_el == 0);
    // Insert the new coefficients for semi-variables
    let mut row_num = num_row;
    for i_col in 0..num_col {
        if lp.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
            || lp.integrality_[i_col as usize] == HighsVarType::kSemiInteger
        {
            let i_el = lp.a_matrix_.start_[(i_col + 1) as usize] - 2;
            lp.a_matrix_.index_[i_el as usize] = row_num;
            row_num += 1;
            lp.a_matrix_.value_[i_el as usize] = 1.0;
            let i_el = i_el + 1;
            lp.a_matrix_.index_[i_el as usize] = row_num;
            row_num += 1;
            lp.a_matrix_.value_[i_el as usize] = 1.0;
        }
    }
    let num_nz = lp.a_matrix_.start_[num_col as usize];
    let new_num_nz = num_nz + 2 * num_semi_variables;
    let mut row_num = num_row;
    let mut semi_col_num: HighsInt = 0;
    let mut semi_row_num: HighsInt = 0;
    // Insert the new variables and their coefficients
    // HiGHS: std::stringstream ss; (replaced by format! strings)
    let have_col_names = lp.col_names_.len() == lp.num_col_ as usize;
    let have_row_names = lp.row_names_.len() == lp.num_row_ as usize;
    let have_solution = solution.value_valid;
    if have_solution {
        // Create zeroed row values for the new rows
        debug_assert!(solution.row_value.len() as HighsInt == lp_.num_row_);
        for _i_col in 0..2 * num_semi_variables {
            solution.row_value.push(0.0);
        }
        debug_assert!(solution.col_value.len() as HighsInt == lp_.num_col_);
        debug_assert!(
            solution.row_value.len() as HighsInt == lp_.num_row_ + 2 * num_semi_variables
        );
    }
    for i_col in 0..num_col {
        if lp.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
            || lp.integrality_[i_col as usize] == HighsVarType::kSemiInteger
        {
            let semi_lower_bound = lp.col_lower_[i_col as usize];
            let semi_upper_bound = lp.col_upper_[i_col as usize];
            // Add a binary variable with zero cost
            lp.col_cost_.push(0.0);
            lp.col_lower_.push(0.0);
            lp.col_upper_.push(1.0);
            // Complete x - l*y >= 0
            lp.row_lower_.push(0.0);
            lp.row_upper_.push(kHighsInf);
            if have_col_names {
                // Create a column name
                lp.col_names_.push(format!("semi_binary_{}", semi_col_num));
                semi_col_num += 1;
            }
            if have_row_names {
                // Create a row name
                lp.row_names_.push(format!("semi_lb_{}", semi_row_num));
            }
            lp.a_matrix_.index_.push(row_num);
            row_num += 1;
            lp.a_matrix_.value_.push(-semi_lower_bound);
            // Accommodate any primal solution
            if have_solution {
                // Record the previous solution value so any change can be
                // determined
                let prev_primal = solution.col_value[i_col as usize];
                if solution.col_value[i_col as usize] <= mip_feasibility_tolerance {
                    // Currently at or below zero, so binary is 0
                    solution.col_value[i_col as usize] = 0.0;
                    solution.col_value.push(0.0);
                } else {
                    // Otherwise, solution is at least lower bound, and binary
                    // is 1
                    solution.col_value[i_col as usize] =
                        cpp_max(semi_lower_bound, solution.col_value[i_col as usize]);
                    solution.col_value.push(1.0);
                }
                let dl_primal = solution.col_value[i_col as usize] - prev_primal;
                if dl_primal != 0.0 {
                    // Change in primal value, so update row values. NB start
                    // has been extended to incorporate the values in this
                    // column for the new rows. Their solution values are zero,
                    // but will be set later
                    let mut i_el = lp.a_matrix_.start_[i_col as usize];
                    while i_el < lp.a_matrix_.start_[(i_col + 1) as usize] {
                        solution.row_value[lp.a_matrix_.index_[i_el as usize] as usize] +=
                            dl_primal * lp.a_matrix_.value_[i_el as usize];
                        i_el += 1;
                    }
                }
                let new_col = lp.col_cost_.len() as HighsInt - 1;
                let binary_value = solution.col_value[new_col as usize];
                solution.row_value[row_num as usize - 1] =
                    solution.col_value[i_col as usize] - semi_lower_bound * binary_value;
                solution.row_value[row_num as usize] =
                    solution.col_value[i_col as usize] - semi_upper_bound * binary_value;
            }
            // Complete x - u*y <= 0
            lp.row_lower_.push(-kHighsInf);
            lp.row_upper_.push(0.0);
            if have_row_names {
                // Create a row name
                lp.row_names_.push(format!("semi_ub_{}", semi_row_num));
                semi_row_num += 1;
            }
            lp.a_matrix_.index_.push(row_num);
            row_num += 1;
            lp.a_matrix_.value_.push(-semi_upper_bound);
            // Add the next start
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
            lp.integrality_.push(HighsVarType::kInteger);
            if lp.integrality_[i_col as usize] == HighsVarType::kSemiContinuous {
                lp.integrality_[i_col as usize] = HighsVarType::kContinuous;
            } else if lp.integrality_[i_col as usize] == HighsVarType::kSemiInteger {
                lp.integrality_[i_col as usize] = HighsVarType::kInteger;
            }
            // Change the lower bound on the semi-variable to zero. Cannot
            // do this earlier, as its original value is used in constraint
            // 0 <= x-l*y
            lp.col_lower_[i_col as usize] = 0.0;
        }
    }
    let num_col = num_col + num_semi_variables;
    lp.num_col_ += num_semi_variables;
    lp.num_row_ += 2 * num_semi_variables;
    debug_assert!(lp.a_matrix_.index_.len() as HighsInt == new_num_nz);
    // Ensure that the matrix dimensions are consistent with the LP
    // dimensions
    lp.a_matrix_.num_col_ = lp.num_col_;
    lp.a_matrix_.num_row_ = lp.num_row_;
    // Clear any modifications inherited from lp_
    lp.mods_.clear();
    lp
}

/// HiGHS: HighsLpUtils.cpp::activeModifiedUpperBounds (lines 687-722).
pub fn active_modified_upper_bounds(
    options: &HighsOptions,
    lp: &HighsLp,
    col_value: &[f64],
) -> bool {
    // HiGHS: const std::vector<HighsInt>& tightened_semi_variable_upper_bound_index =
    //     lp.mods_.save_tightened_semi_variable_upper_bound_index;
    let tightened_semi_variable_upper_bound_index =
        &lp.mods_.save_tightened_semi_variable_upper_bound_index;
    let num_tightened_upper = tightened_semi_variable_upper_bound_index.len() as HighsInt;
    let mut num_active_modified_upper: HighsInt = 0;
    let mut min_semi_variable_margin = kHighsInf;
    for k in 0..num_tightened_upper {
        let value = col_value[tightened_semi_variable_upper_bound_index[k as usize] as usize];
        let upper = lp.col_upper_[tightened_semi_variable_upper_bound_index[k as usize] as usize];
        let semi_variable_margin = upper - value;
        if value > upper - options.primal_feasibility_tolerance {
            min_semi_variable_margin = 0.0;
            num_active_modified_upper += 1;
        } else {
            min_semi_variable_margin = cpp_min(semi_variable_margin, min_semi_variable_margin);
        }
    }
    if num_active_modified_upper != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "{} semi-variables are active at modified upper bounds\n",
            num_active_modified_upper
        );
    } else if num_tightened_upper != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "No semi-variables are active at modified upper bounds: a large minimum margin ({}) suggests optimality, but there is no guarantee\n",
            crate::util::c_printf::c_fmt_g(min_semi_variable_margin, 6)
        );
    }
    num_active_modified_upper != 0
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// One row `x0 + x1 <= 10`; x0 semi-continuous in {0} U [2, 6], x1 continuous.
    fn semi_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![2.0, 0.0];
        lp.col_upper_ = vec![6.0, 4.0];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![10.0];
        lp.integrality_ = vec![HighsVarType::kSemiContinuous, HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        lp
    }

    #[test]
    fn adds_binary_column_and_two_rows_per_semi_variable() {
        let mut solution = HighsSolution::default();
        let lp = without_semi_variables(&semi_lp(), &mut solution, 1e-7);
        assert_eq!((lp.num_col_, lp.num_row_), (3, 3));
        assert_eq!(lp.integrality_[0], HighsVarType::kContinuous);
        assert_eq!(lp.integrality_[2], HighsVarType::kInteger);
        assert_eq!(lp.col_lower_[0], 0.0);
        assert_eq!((lp.col_lower_[2], lp.col_upper_[2]), (0.0, 1.0));
        assert_eq!(lp.a_matrix_.start_, vec![0, 3, 4, 6]);
        assert_eq!(lp.a_matrix_.index_, vec![0, 1, 2, 0, 1, 2]);
        assert_eq!(lp.a_matrix_.value_, vec![1.0, 1.0, 1.0, 1.0, -2.0, -6.0]);
        assert_eq!(lp.row_lower_, vec![-kHighsInf, 0.0, -kHighsInf]);
        assert_eq!(lp.row_upper_, vec![10.0, kHighsInf, 0.0]);
        assert!(!solution.value_valid);
    }

    #[test]
    fn extends_a_user_solution() {
        let mut solution = HighsSolution::default();
        solution.col_value = vec![1.0, 3.0];
        solution.row_value = vec![4.0];
        solution.value_valid = true;
        let _lp = without_semi_variables(&semi_lp(), &mut solution, 1e-7);
        // 1.0 > tolerance: raised to the semi lower bound 2, indicator 1
        assert_eq!(solution.col_value, vec![2.0, 3.0, 1.0]);
        assert_eq!(solution.row_value, vec![5.0, 0.0, -4.0]);
    }

    #[test]
    fn no_tightened_upper_bounds_is_inactive() {
        let lp = semi_lp();
        assert!(!active_modified_upper_bounds(
            &HighsOptions::default(),
            &lp,
            &[0.0, 0.0]
        ));
    }
}
