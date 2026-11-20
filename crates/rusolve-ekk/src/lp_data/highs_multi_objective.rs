//! HiGHS source: lp_data/Highs.cpp (optimizeHighs, passLinearObjectives,
//! addLinearObjective, clearLinearObjectives, clearSolverDualData,
//! invalidateSolverDualData), lp_data/HighsInterface.cpp (validLinearObjective,
//! hasRepeatedLinearObjectivePriorities, comparison,
//! returnFromLexicographicOptimization, multiobjectiveSolve), Highs.h
//! (getNumLinearObjectives, getLinearObjective)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Multiple linear objectives, solved by blending (one weighted
//! objective, `blend_multi_objectives = true`) or lexicographically (by
//! decreasing priority, each optimum fixed by an added row before the next).
//! `Highs::run` calls `optimize_highs`, which is in place: with no linear
//! objectives (the default) it is `optimize_model`, so `run` is unchanged.
//!

use crate::prelude::*;

/// HiGHS: HighsInterface.cpp::comparison (lines 3913-3916), the `static`
/// comparator for `std::sort` in `multiobjectiveSolve`: priority descending.
/// HiGHS quirk: `std::sort` is not stable; the priorities are distinct
/// (`hasRepeatedLinearObjectivePriorities` is checked first), so the order
/// is fully determined.
fn comparison(x1: (HighsInt, HighsInt), x2: (HighsInt, HighsInt)) -> bool {
    x1.0 > x2.0
}

impl Highs {
    /// HiGHS: Highs.cpp::Highs::optimizeHighs (lines 1137-1143)
    pub fn optimize_highs(&mut self) -> HighsStatus {
        // Level 1 of Highs::run()
        //
        // Move the "mods" to here
        if !self.multi_linear_objective_.is_empty() {
            self.multiobjective_solve()
        } else {
            self.optimize_model()
        }
    }

    /// HiGHS: Highs.cpp::Highs::passLinearObjectives (lines 713-723).
    /// `const HighsLinearObjective* linear_objective` -> a slice of
    /// `num_linear_objective` entries.
    pub fn pass_linear_objectives(
        &mut self,
        num_linear_objective: HighsInt,
        linear_objective: &[HighsLinearObjective],
    ) -> HighsStatus {
        if num_linear_objective < 0 {
            return HighsStatus::kOk;
        }
        self.multi_linear_objective_.clear();
        for i_obj in 0..num_linear_objective {
            if self.add_linear_objective(&linear_objective[i_obj as usize], i_obj)
                != HighsStatus::kOk
            {
                return HighsStatus::kError;
            }
        }
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::addLinearObjective (lines 725-736).
    /// C++ default argument: `iObj = -1`. `model_.isQp()` is
    /// `model_hessian_.dim_ != 0` (HighsModel.h line 27).
    pub fn add_linear_objective(
        &mut self,
        linear_objective: &HighsLinearObjective,
        i_obj: HighsInt,
    ) -> HighsStatus {
        if self.model_hessian_.dim_ != 0 {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot define additional linear objective for QP\n"
            );
            return HighsStatus::kError;
        }
        if !self.valid_linear_objective(linear_objective, i_obj) {
            return HighsStatus::kError;
        }
        self.multi_linear_objective_.push(linear_objective.clone());
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::clearLinearObjectives (lines 738-741)
    pub fn clear_linear_objectives(&mut self) -> HighsStatus {
        self.multi_linear_objective_.clear();
        HighsStatus::kOk
    }

    /// HiGHS: Highs.h::Highs::getNumLinearObjectives (lines 181-183)
    pub fn get_num_linear_objectives(&self) -> HighsInt {
        self.multi_linear_objective_.len() as HighsInt
    }

    /// HiGHS: Highs.h::Highs::getLinearObjective (lines 188-192)
    pub fn get_linear_objective(&self, idx: HighsInt) -> &HighsLinearObjective {
        debug_assert!(idx >= 0 && idx < self.multi_linear_objective_.len() as HighsInt);
        &self.multi_linear_objective_[idx as usize]
    }

    /// HiGHS: HighsInterface.cpp::Highs::validLinearObjective (lines 3867-3889).
    /// `iObj >= 0 ? std::to_string(iObj) : ""` -> `format!("{i_obj}")` or "".
    pub fn valid_linear_objective(
        &self,
        linear_objective: &HighsLinearObjective,
        i_obj: HighsInt,
    ) -> bool {
        let linear_objective_coefficients_size = linear_objective.coefficients.len() as HighsInt;
        if linear_objective_coefficients_size != self.model_lp_.num_col_ {
            let obj_str = if i_obj >= 0 {
                format!("{}", i_obj)
            } else {
                String::new()
            };
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Coefficient vector for linear objective {} has size {} != {} = lp.num_col_\n",
                obj_str,
                linear_objective_coefficients_size,
                self.model_lp_.num_col_
            );
            return false;
        }
        if !self.options_.blend_multi_objectives
            && self.has_repeated_linear_objective_priorities(Some(linear_objective))
        {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Repeated priorities for lexicographic optimization is illegal\n"
            );
            return false;
        }
        true
    }

    /// HiGHS: HighsInterface.cpp::Highs::hasRepeatedLinearObjectivePriorities
    /// (lines 3891-3911). C++ default argument: `linear_objective = nullptr`.
    pub fn has_repeated_linear_objective_priorities(
        &self,
        linear_objective: Option<&HighsLinearObjective>,
    ) -> bool {
        // Look for repeated values in the linear objective priorities, also
        // comparing linear_objective if it's not a null pointer. Cost is
        // O(n^2), but who will have more than O(1) linear objectives!
        let num_linear_objective = self.multi_linear_objective_.len() as HighsInt;
        if num_linear_objective <= 0 || (num_linear_objective <= 1 && linear_objective.is_none()) {
            return false;
        }
        for i_obj0 in 0..num_linear_objective {
            let priority0 = self.multi_linear_objective_[i_obj0 as usize].priority;
            for i_obj1 in (i_obj0 + 1)..num_linear_objective {
                let priority1 = self.multi_linear_objective_[i_obj1 as usize].priority;
                if priority1 == priority0 {
                    return true;
                }
            }
            if let Some(lo) = linear_objective {
                if lo.priority == priority0 {
                    return true;
                }
            }
        }
        false
    }

    /// HiGHS: HighsInterface.cpp::Highs::returnFromLexicographicOptimization
    /// (lines 3918-3938). `deleteRows(from, to)` -> `delete_rows_interval`.
    pub fn return_from_lexicographic_optimization(
        &mut self,
        return_status: HighsStatus,
        original_lp_num_row: HighsInt,
    ) -> HighsStatus {
        // Save model_status_ and info_ since they are cleared by calling
        // deleteRows
        let model_status = self.model_status_;
        let info = self.info_.clone();
        let num_linear_objective = self.multi_linear_objective_.len() as HighsInt;
        if num_linear_objective > 1 {
            self.delete_rows_interval(original_lp_num_row, self.model_lp_.num_row_ - 1);
            // Recover model_status_ and info_, and then account for lack of basis or
            // dual solution
            self.model_status_ = model_status;
            self.info_ = info;
            self.info_.objective_function_value = 0.0;
            self.info_.basis_validity = kBasisValidityInvalid;
            self.info_.invalidate_dual_kkt();
            self.solution_.value_valid = true;
            self.model_lp_
                .col_cost_
                .resize(self.model_lp_.num_col_ as usize, 0.0);
            for i in 0..self.model_lp_.num_col_ as usize {
                self.model_lp_.col_cost_[i] = 0.0;
            }
        }
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::multiobjectiveSolve (lines 3940-4208).
    /// `std::stringstream` reports -> `String` with `format!`; `std::sort(...,
    /// comparison)` -> `sort_by` with [`comparison`] mapped to an `Ordering`.
    pub fn multiobjective_solve(&mut self) -> HighsStatus {
        let coeff_logging_size_limit: HighsInt = 10;
        let num_linear_objective = self.multi_linear_objective_.len() as HighsInt;

        debug_assert!(num_linear_objective > 0);

        for i_obj in 0..num_linear_objective {
            let multi_linear_objective = &self.multi_linear_objective_[i_obj as usize];
            if multi_linear_objective.coefficients.len() != self.model_lp_.num_col_ as usize {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Multiple linear objective coefficient vector {} has size incompatible with model\n",
                    i_obj
                );
                return HighsStatus::kError;
            }
        }

        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "Solving with {} multiple linear objectives, {}\n",
            num_linear_objective,
            if self.options_.blend_multi_objectives {
                "blending objectives by weight"
            } else {
                "using lexicographic optimization by priority"
            }
        );
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "Ix      weight      offset     abs_tol     rel_tol    priority{}\n",
            if self.model_lp_.num_col_ < coeff_logging_size_limit {
                "   coefficients"
            } else {
                ""
            }
        );
        for i_obj in 0..num_linear_objective {
            let linear_objective = &self.multi_linear_objective_[i_obj as usize];
            let mut multi_objective_log = String::new();
            multi_objective_log.push_str(&format!(
                "{:>2} {} {} {} {} {:>11}  ",
                i_obj,
                crate::util::c_printf::c_fmt_g_width(linear_objective.weight, 6, 11),
                crate::util::c_printf::c_fmt_g_width(linear_objective.offset, 6, 11),
                crate::util::c_printf::c_fmt_g_width(linear_objective.abs_tolerance, 6, 11),
                crate::util::c_printf::c_fmt_g_width(linear_objective.rel_tolerance, 6, 11),
                linear_objective.priority
            ));
            if self.model_lp_.num_col_ < coeff_logging_size_limit {
                for i_col in 0..self.model_lp_.num_col_ {
                    multi_objective_log.push_str(&format!(
                        "{} c_{{{}}} = {}",
                        if i_col == 0 { "" } else { "," },
                        i_col,
                        crate::util::c_printf::c_fmt_g(
                            linear_objective.coefficients[i_col as usize],
                            6
                        )
                    ));
                }
            }
            multi_objective_log.push('\n');
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "{}",
                multi_objective_log
            );
        }
        // Solving with a different objective, but don't call
        // this->clearSolver() since this loses the current solution - that
        // may have been provided by the user (#2419). Just clear the dual
        // data.
        //
        self.clear_solver_dual_data();
        if self.options_.blend_multi_objectives {
            // Objectives are blended by weight and minimized
            self.model_lp_.offset_ = 0.0;
            self.model_lp_
                .col_cost_
                .resize(self.model_lp_.num_col_ as usize, 0.0);
            for i in 0..self.model_lp_.num_col_ as usize {
                self.model_lp_.col_cost_[i] = 0.0;
            }
            for i_obj in 0..num_linear_objective {
                let multi_linear_objective = &self.multi_linear_objective_[i_obj as usize];
                self.model_lp_.offset_ +=
                    multi_linear_objective.weight * multi_linear_objective.offset;
                for i_col in 0..self.model_lp_.num_col_ {
                    self.model_lp_.col_cost_[i_col as usize] += multi_linear_objective.weight
                        * multi_linear_objective.coefficients[i_col as usize];
                }
            }
            self.model_lp_.sense_ = ObjSense::kMinimize;

            let mut multi_objective_log = String::new();
            multi_objective_log.push_str("Solving with blended objective");
            if self.model_lp_.num_col_ < coeff_logging_size_limit {
                multi_objective_log.push_str(&format!(
                    ": {} {}",
                    if self.model_lp_.sense_ == ObjSense::kMinimize {
                        "min"
                    } else {
                        "max"
                    },
                    crate::util::c_printf::c_fmt_g(self.model_lp_.offset_, 6)
                ));
                for i_col in 0..self.model_lp_.num_col_ {
                    multi_objective_log.push_str(&format!(
                        " + ({}) x[{}]",
                        crate::util::c_printf::c_fmt_g(self.model_lp_.col_cost_[i_col as usize], 6),
                        i_col
                    ));
                }
            }
            multi_objective_log.push('\n');
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "{}",
                multi_objective_log
            );
            return self.optimize_model();
        }

        // Objectives are applied lexicographically
        if self.model_hessian_.dim_ != 0 && num_linear_objective > 1 {
            // Lexicographic optimization with a single linear objective is
            // trivially standard optimization, so is OK
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot perform non-trivial lexicographic optimization for QP\n"
            );
            return HighsStatus::kError;
        }
        // Check whether there are repeated linear objective priorities
        if self.has_repeated_linear_objective_priorities(None) {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Repeated priorities for lexicographic optimization is illegal\n"
            );
            return HighsStatus::kError;
        }
        let mut priority_objective: Vec<(HighsInt, HighsInt)> = Vec::new();

        for i_obj in 0..num_linear_objective {
            priority_objective.push((self.multi_linear_objective_[i_obj as usize].priority, i_obj));
        }
        // PORT NOTE: std::sort, keys distinct
        priority_objective.sort_by(|a, b| {
            if comparison(*a, *b) {
                std::cmp::Ordering::Less
            } else if comparison(*b, *a) {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });
        // Clear LP objective
        self.model_lp_.offset_ = 0.0;
        self.model_lp_
            .col_cost_
            .resize(self.model_lp_.num_col_ as usize, 0.0);
        for i in 0..self.model_lp_.num_col_ as usize {
            self.model_lp_.col_cost_[i] = 0.0;
        }
        let original_lp_num_row = self.model_lp_.num_row_;
        let mut index: Vec<HighsInt> = vec![0; self.model_lp_.num_col_ as usize];
        let mut value: Vec<f64> = vec![0.0; self.model_lp_.num_col_ as usize];
        // Use the solution of one MIP to provide an integer feasible
        // solution of the next
        let mut solution = HighsSolution::default();
        for i_ix in 0..num_linear_objective {
            let priority = priority_objective[i_ix as usize].0;
            let i_obj = priority_objective[i_ix as usize].1;
            // Use this objective
            let linear_objective = self.multi_linear_objective_[i_obj as usize].clone();
            self.model_lp_.offset_ = linear_objective.offset;
            self.model_lp_.col_cost_ = linear_objective.coefficients.clone();
            self.model_lp_.sense_ = if linear_objective.weight > 0.0 {
                ObjSense::kMinimize
            } else {
                ObjSense::kMaximize
            };
            if self.model_lp_.is_mip() && solution.value_valid {
                let set_solution_status = self.set_solution(&solution);
                if set_solution_status == HighsStatus::kError {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Failure to use one MIP to provide an integer feasible solution of the next\n"
                    );
                    return self.return_from_lexicographic_optimization(
                        HighsStatus::kError,
                        original_lp_num_row,
                    );
                }
                let mut valid = false;
                let mut integral = false;
                let mut feasible = false;
                let _assess_primal_solution =
                    self.assess_primal_solution(&mut valid, &mut integral, &mut feasible);
                if !valid || !integral || !feasible {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kWarning,
                        "Failure to use one MIP to provide an integer feasible solution of the next: status is valid = {}, integral = {}, feasible = {}\n",
                        highs_bool_to_string(valid, 2),
                        highs_bool_to_string(integral, 2),
                        highs_bool_to_string(feasible, 2)
                    );
                }
            }
            let mut multi_objective_log = String::new();
            multi_objective_log.push_str(&format!("Solving with objective {}", i_obj));
            if self.model_lp_.num_col_ < coeff_logging_size_limit {
                multi_objective_log.push_str(&format!(
                    ": {} {}",
                    if self.model_lp_.sense_ == ObjSense::kMinimize {
                        "min"
                    } else {
                        "max"
                    },
                    crate::util::c_printf::c_fmt_g(self.model_lp_.offset_, 6)
                ));
                for i_col in 0..self.model_lp_.num_col_ {
                    multi_objective_log.push_str(&format!(
                        " + ({}) x[{}]",
                        crate::util::c_printf::c_fmt_g(self.model_lp_.col_cost_[i_col as usize], 6),
                        i_col
                    ));
                }
            }
            multi_objective_log.push('\n');
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "{}",
                multi_objective_log
            );
            let optimize_model_status = self.optimize_model();
            if optimize_model_status == HighsStatus::kError {
                return self.return_from_lexicographic_optimization(
                    HighsStatus::kError,
                    original_lp_num_row,
                );
            }
            if self.model_status_ != HighsModelStatus::kOptimal {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kWarning,
                    "After priority {} solve, model status is {}\n",
                    priority,
                    crate::lp_data::model_utils::util_model_status_to_string(self.model_status_)
                );
                return self.return_from_lexicographic_optimization(
                    HighsStatus::kWarning,
                    original_lp_num_row,
                );
            }
            if i_ix == num_linear_objective - 1 {
                break;
            }
            if self.model_lp_.is_mip() {
                // Save the solution to provide an integer feasible solution of
                // the next MIP
                solution.col_value = self.solution_.col_value.clone();
                solution.value_valid = true;
            }
            // Add the constraint
            let mut nnz: HighsInt = 0;
            for i_col in 0..self.model_lp_.num_col_ {
                if self.model_lp_.col_cost_[i_col as usize] != 0.0 {
                    index[nnz as usize] = i_col;
                    value[nnz as usize] = self.model_lp_.col_cost_[i_col as usize];
                    nnz += 1;
                }
            }
            let objective = self.info_.objective_function_value;
            let mut lower_bound = -kHighsInf;
            let mut upper_bound = kHighsInf;
            if self.model_lp_.sense_ == ObjSense::kMinimize {
                // Minimizing, so set a greater upper bound than the objective
                if linear_objective.abs_tolerance >= 0.0 {
                    upper_bound = objective + linear_objective.abs_tolerance;
                }
                if linear_objective.rel_tolerance >= 0.0 {
                    if objective >= 0.0 {
                        // Guarantees objective of at least (1+t).f^*
                        //
                        // so ((1+t).f^*-f^*)/f^* = t
                        upper_bound = cpp_min(
                            objective * (1.0 + linear_objective.rel_tolerance),
                            upper_bound,
                        );
                    } else if objective < 0.0 {
                        // Guarantees objective of at least (1-t).f^*
                        //
                        // so ((1-t).f^*-f^*)/f^* = -t
                        upper_bound = cpp_min(
                            objective * (1.0 - linear_objective.rel_tolerance),
                            upper_bound,
                        );
                    }
                }
                upper_bound -= self.model_lp_.offset_;
            } else {
                // Maximizing, so set a lesser lower bound than the objective
                if linear_objective.abs_tolerance >= 0.0 {
                    lower_bound = objective - linear_objective.abs_tolerance;
                }
                if linear_objective.rel_tolerance >= 0.0 {
                    if objective >= 0.0 {
                        // Guarantees objective of at most (1-t).f^*
                        //
                        // so ((1-t).f^*-f^*)/f^* = -t
                        lower_bound = cpp_max(
                            objective * (1.0 - linear_objective.rel_tolerance),
                            lower_bound,
                        );
                    } else if objective < 0.0 {
                        // Guarantees objective of at least (1+t).f^*
                        //
                        // so ((1+t).f^*-f^*)/f^* = t
                        lower_bound = cpp_max(
                            objective * (1.0 + linear_objective.rel_tolerance),
                            lower_bound,
                        );
                    }
                }
                lower_bound -= self.model_lp_.offset_;
            }
            if lower_bound == -kHighsInf && upper_bound == kHighsInf {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kWarning,
                    "After priority {} solve, no objective constraint due to absolute tolerance being {} < 0, and relative tolerance being {} < 0\n",
                    priority,
                    crate::util::c_printf::c_fmt_g(linear_objective.abs_tolerance, 6),
                    crate::util::c_printf::c_fmt_g(linear_objective.rel_tolerance, 6)
                );
            }
            let mut multi_objective_log = String::new();
            multi_objective_log.push_str(&format!("Add constraint for objective {}: ", i_obj));
            if nnz < coeff_logging_size_limit {
                multi_objective_log.push_str(&format!(
                    "{} <= ",
                    crate::util::c_printf::c_fmt_g(lower_bound, 6)
                ));
                for i_el in 0..nnz {
                    multi_objective_log.push_str(&format!(
                        "{}({}) x[{}]",
                        if i_el > 0 { " + " } else { "" },
                        crate::util::c_printf::c_fmt_g(value[i_el as usize], 6),
                        index[i_el as usize]
                    ));
                }
                multi_objective_log.push_str(&format!(
                    " <= {}\n",
                    crate::util::c_printf::c_fmt_g(upper_bound, 6)
                ));
            } else {
                multi_objective_log.push_str(&format!(
                    "Bounds [{}, {}]\n",
                    crate::util::c_printf::c_fmt_g(lower_bound, 6),
                    crate::util::c_printf::c_fmt_g(upper_bound, 6)
                ));
            }
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "{}",
                multi_objective_log
            );
            let add_row_status = self.add_row(
                lower_bound,
                upper_bound,
                nnz,
                &index[..nnz as usize],
                &value[..nnz as usize],
            );
            debug_assert!(add_row_status == HighsStatus::kOk);
        }
        self.return_from_lexicographic_optimization(HighsStatus::kOk, original_lp_num_row)
    }

    /// HiGHS: Highs.cpp::Highs::clearSolverDualData (lines 78-83)
    pub fn clear_solver_dual_data(&mut self) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;
        self.clear_derived_model_properties();
        self.invalidate_solver_dual_data();
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::invalidateSolverDualData (lines 3925-3930)
    pub(crate) fn invalidate_solver_dual_data(&mut self) {
        self.invalidate_model_status();
        self.invalidate_ranging();
        self.info_.invalidate();
        self.run_data_.invalidate();
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    fn objective(coefficients: Vec<f64>, priority: HighsInt) -> HighsLinearObjective {
        let mut o = HighsLinearObjective::default();
        o.coefficients = coefficients;
        o.priority = priority;
        o.weight = 1.0;
        o
    }

    fn two_col_highs() -> Highs {
        let mut highs = Highs::default();
        highs.model_lp_.num_col_ = 2;
        highs.model_lp_.col_cost_ = vec![0.0; 2];
        highs.model_lp_.col_lower_ = vec![0.0; 2];
        highs.model_lp_.col_upper_ = vec![1.0; 2];
        highs.model_lp_.a_matrix_.start_ = vec![0, 0, 0];
        highs
    }

    #[test]
    fn comparison_orders_by_decreasing_priority() {
        assert!(comparison((3, 0), (1, 1)));
        assert!(!comparison((1, 0), (3, 1)));
    }

    #[test]
    fn optimize_highs_without_objectives_is_optimize_model() {
        let mut highs = Highs::default();
        assert_eq!(highs.get_num_linear_objectives(), 0);
        // An empty model: optimize_model sets kModelEmpty
        assert_eq!(highs.optimize_highs(), HighsStatus::kOk);
        assert_eq!(highs.model_status_, HighsModelStatus::kModelEmpty);
    }

    #[test]
    fn add_checks_size_and_priorities() {
        let mut highs = two_col_highs();
        assert_eq!(
            highs.add_linear_objective(&objective(vec![1.0], 0), -1),
            HighsStatus::kError
        );
        assert_eq!(
            highs.add_linear_objective(&objective(vec![1.0, 0.0], 2), -1),
            HighsStatus::kOk
        );
        assert!(!highs.has_repeated_linear_objective_priorities(None));
        assert!(
            highs.has_repeated_linear_objective_priorities(Some(&objective(vec![0.0, 1.0], 2)))
        );
        // Blending allows repeated priorities; lexicographic does not
        assert_eq!(
            highs.add_linear_objective(&objective(vec![0.0, 1.0], 2), -1),
            HighsStatus::kOk
        );
        highs.options_.blend_multi_objectives = false;
        assert_eq!(
            highs.add_linear_objective(&objective(vec![1.0, 1.0], 2), -1),
            HighsStatus::kError
        );
        assert_eq!(highs.get_num_linear_objectives(), 2);
        assert_eq!(highs.clear_linear_objectives(), HighsStatus::kOk);
        assert_eq!(highs.get_num_linear_objectives(), 0);
    }

    #[test]
    fn pass_replaces_and_negative_count_is_a_no_op() {
        let mut highs = two_col_highs();
        let objs = vec![objective(vec![1.0, 0.0], 1), objective(vec![0.0, 1.0], 0)];
        assert_eq!(highs.pass_linear_objectives(2, &objs), HighsStatus::kOk);
        assert_eq!(highs.get_num_linear_objectives(), 2);
        assert_eq!(highs.pass_linear_objectives(-1, &[]), HighsStatus::kOk);
        assert_eq!(highs.get_num_linear_objectives(), 2);
        assert_eq!(
            highs.pass_linear_objectives(1, &objs[1..]),
            HighsStatus::kOk
        );
        assert_eq!(highs.get_linear_objective(0).priority, 0);
    }
}
