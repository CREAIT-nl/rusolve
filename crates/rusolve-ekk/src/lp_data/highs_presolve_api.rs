//! HiGHS source: lp_data/Highs.cpp (presolve, postsolve(solution),
//! crossover, the MIP branches of runPresolve and callRunPostsolve),
//! lp_data/HighsSolution.cpp (computeObjectiveValue), Highs.h
//! (getPresolvedLp, getModelPresolveStatus)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The public presolve / postsolve / crossover entry points.
//! `run` never reaches them: `runPresolve` and `callRunPostsolve` are
//! already ported (highs_presolve.rs, highs_postsolve.rs) for the LP path,
//! and their MIP branches call `run_presolve_mip` / `call_run_postsolve_mip`
//! here, which `run` cannot reach (a MIP goes to `callSolveMip` before any
//! `runPresolve`, and `postsolve` is never called by `run`).
//!
//! `run_presolve_mip` builds a `HighsMipSolver` on `model_lp_` exactly as
//! `call_solve_mip` does (highs_run.rs: callback moved in and back) and runs
//! `run_mip_presolve(options_.presolve_reduction_limit)`.

use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::presolve (lines 987-1086).
    pub fn presolve(&mut self) -> HighsStatus {
        // Highs.cpp:988 — const HighsLogOptions& log_options = options_.log_options;
        // (accessed via self.options_.log_options in the log macros)

        // Highs.cpp:989-994
        if self.model_lp_.needs_mods(self.options_.infinite_cost) {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Model contains infinite costs or semi-variables, so cannot be presolved independently\n"
            );
            return HighsStatus::kError;
        }

        // Highs.cpp:995-1000
        if self.model_hessian_.dim_ != 0 {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Model is a QP, for which no presolve techniques are implemented\n"
            );
            return HighsStatus::kError;
        }

        // Highs.cpp:1001
        // Every switch arm below assigns return_status (Highs.cpp:1001 init is dead).
        let return_status: HighsStatus;

        // Highs.cpp:1002 — reportModelStats() is print-only; dropped.
        // self.report_model_stats();

        // Highs.cpp:1003
        self.clear_presolve();

        // Highs.cpp:1004-1024
        if self.model_lp_.num_col_ == 0 && self.model_lp_.num_row_ == 0 {
            self.model_presolve_status_ = HighsPresolveStatus::kNotReduced;
        } else {
            // Highs.cpp:1007
            let force_presolve = true;

            // Highs.cpp:1011-1012
            let status = self.initialize_multi_threading();
            if status != HighsStatus::kOk {
                return status;
            }

            // Highs.cpp:1013-1015 — profiling is not ported; dropped.
            // Highs.cpp:1023 — clearProfiling() dropped.

            // Highs.cpp:1021
            let force_lp_presolve = self.options_.solve_relaxation;

            // Highs.cpp:1022
            self.model_presolve_status_ = self.run_presolve(force_lp_presolve, force_presolve);
        }

        // Highs.cpp:1026
        let mut using_reduced_lp = false;

        // Highs.cpp:1027-1028
        crate::lp_data::presolve_support::report_presolve_reductions(
            &self.options_.log_options,
            self.model_presolve_status_,
            &self.model_lp_,
            self.presolve_.get_reduced_problem(),
        );

        // Highs.cpp:1029-1077
        match self.model_presolve_status_ {
            HighsPresolveStatus::kNotPresolved => {
                // Highs.cpp:1032
                debug_assert!(self.model_presolve_status_ != HighsPresolveStatus::kNotPresolved);
                return_status = HighsStatus::kError;
            }
            HighsPresolveStatus::kNotReduced
            | HighsPresolveStatus::kInfeasible
            | HighsPresolveStatus::kReduced
            | HighsPresolveStatus::kReducedToEmpty
            | HighsPresolveStatus::kUnboundedOrInfeasible => {
                // Highs.cpp:1042-1057
                if self.model_presolve_status_ == HighsPresolveStatus::kInfeasible {
                    self.set_highs_model_status_and_clear_solution_and_basis(
                        HighsModelStatus::kInfeasible,
                    );
                } else if self.model_presolve_status_ == HighsPresolveStatus::kNotReduced {
                    self.presolved_model_lp_ = self.model_lp_.clone();
                } else if self.model_presolve_status_ == HighsPresolveStatus::kReduced
                    || self.model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty
                {
                    using_reduced_lp = true;
                }
                return_status = HighsStatus::kOk;
            }
            HighsPresolveStatus::kTimeout => {
                using_reduced_lp = true;
                return_status = HighsStatus::kWarning;
            }
            _ => {
                // Highs.cpp:1070
                debug_assert!(self.model_presolve_status_ == HighsPresolveStatus::kOutOfMemory);
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Presolve fails due to memory allocation error\n"
                );
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kPresolveError,
                );
                return_status = HighsStatus::kError;
            }
        }

        // Highs.cpp:1078-1081
        if using_reduced_lp {
            self.presolved_model_lp_ = self.presolve_.get_reduced_problem().clone();
            self.presolved_model_lp_.set_matrix_dimensions();
        }

        // Highs.cpp:1083-1084
        let status_str = self.presolve_status_to_string(self.model_presolve_status_);
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "Presolve status: {}\n",
            status_str
        );

        // Highs.cpp:1085
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::runPresolve (lines 3793-3811), the
    /// `model_.isMip() && !force_lp_presolve` branch, called from
    /// `run_presolve` (highs_presolve.rs). `presolve_.data_.reduced_lp_ =
    /// solver.getPresolvedModel()` and `presolve_.data_.postSolveStack =
    /// solver.getPostsolveStack()` are clones out of the solver.
    pub(crate) fn run_presolve_mip(&mut self) -> HighsPresolveStatus {
        // Highs.cpp:3798 — HighsMipSolver solver(callback_, options_, original_lp, solution_);
        // Rule (f): callback_ is moved into the solver and restored after.
        let mut solver = HighsMipSolver::new(
            &mut self.callback_,
            &self.options_,
            &self.model_lp_,
            &self.solution_,
            false,
            0,
        );

        // Highs.cpp:3799 — solver.setProfiling(this->profiling_); profiling is a
        // no-op stub (as in call_solve_mip).
        solver.set_profiling();

        // Highs.cpp:3800-3802: start the MIP solver's timer so that timeout in
        // presolve can be identified
        solver.timer_.start();

        // Highs.cpp:3803-3804
        solver.run_mip_presolve(self.options_.presolve_reduction_limit);

        // Rule (f): restore callback_ first (as call_solve_mip does); C++ shares
        // it by reference and nothing below touches it.
        self.callback_ = std::mem::take(&mut solver.callback_);

        // Highs.cpp:3805
        let presolve_return_status = solver.get_presolve_status();

        // Highs.cpp:3807
        self.presolve_.data_.reduced_lp_ = solver.get_presolved_model().clone();

        // Highs.cpp:3808
        self.presolve_.data_.post_solve_stack = solver.get_postsolve_stack();

        // Highs.cpp:3809
        self.presolve_.presolve_status_ = presolve_return_status;

        // Highs.cpp:3810 — commented-out line dropped.

        presolve_return_status
    }

    /// HiGHS: Highs.cpp::Highs::callRunPostsolve (lines 4483-4527), the
    /// `isMip() && !basis.valid` branch, called from `call_run_postsolve`
    /// (highs_postsolve.rs). `fractionality` is util::highs_utils.
    pub(crate) fn call_run_postsolve_mip(&mut self) {
        // Highs.cpp:4489-4491
        self.presolve_.data_.recovered_solution_.dual_valid = false;
        self.presolve_.data_.recovered_solution_.col_dual.clear();
        self.presolve_.data_.recovered_solution_.row_dual.clear();

        // Highs.cpp:4493
        self.presolve_.data_.recovered_basis_.valid = false;

        // Highs.cpp:4495
        let postsolve_status = self.run_postsolve();

        // Highs.cpp:4497-4521
        if postsolve_status == HighsPostsolveStatus::kSolutionRecovered {
            // Highs.cpp:4498
            self.solution_ = self.presolve_.data_.recovered_solution_.clone();

            // Highs.cpp:4499
            self.model_status_ = HighsModelStatus::kUnknown;

            // Highs.cpp:4500
            self.info_.invalidate();

            // Highs.cpp:4501-4503
            // HighsLp& lp = this->model_.lp_;
            // this->info_.objective_function_value = computeObjectiveValue(lp, this->solution_);
            self.info_.objective_function_value =
                compute_objective_value(&self.model_lp_, &self.solution_);

            // Highs.cpp:4504-4505
            let is_qp = self.model_hessian_.dim_ != 0;
            debug_assert!(!is_qp);

            // Highs.cpp:4506-4509
            let get_residuals = true;
            crate::lp_data::presolve_support::get_kkt_failures(
                &self.options_,
                is_qp,
                &self.model_lp_,
                &self.model_lp_.col_cost_,
                &self.solution_,
                &mut self.info_,
                get_residuals,
            );

            // Highs.cpp:4510-4518
            let max_integrality_violation = &mut self.info_.max_integrality_violation;
            *max_integrality_violation = 0.0;
            for i_col in 0..self.model_lp_.num_col_ {
                if self.model_lp_.integrality_[i_col as usize] == HighsVarType::kInteger {
                    *max_integrality_violation = crate::util::highs_utils::cpp_max(
                        crate::util::highs_utils::fractionality(
                            self.solution_.col_value[i_col as usize],
                            None,
                        ),
                        *max_integrality_violation,
                    );
                }
            }

            // Highs.cpp:4519-4521
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "Postsolve performed for MIP, but model status cannot be known\n"
            );
        } else {
            // Highs.cpp:4523-4526
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Postsolve return status is {}\n",
                postsolve_status as i32
            );
            self.set_highs_model_status_and_clear_solution_and_basis(
                HighsModelStatus::kPostsolveError,
            );
        }
    }

    /// HiGHS: Highs.cpp::Highs::postsolve (lines 3601-3604), the overload
    /// `postsolve(const HighsSolution& solution)`: an empty basis.
    pub fn postsolve_solution(&mut self, solution: &HighsSolution) -> HighsStatus {
        // Highs.cpp:3602-3603
        let basis = HighsBasis::default();
        self.postsolve(solution, &basis)
    }

    /// HiGHS: Highs.cpp::Highs::crossover (lines 5064-5089).
    /// `callCrossover` -> `presolve::icrash_x::call_crossover`;
    /// `getLpKktFailures` -> `presolve_support::get_lp_kkt_failures`.
    pub fn crossover(&mut self, user_solution: &HighsSolution) -> HighsStatus {
        // Highs.cpp:5065: every branch below assigns return_status (init is dead).
        let return_status: HighsStatus;

        // Highs.cpp:5066-5067
        // HighsLogOptions& log_options = options_.log_options;
        // HighsLp& lp = model_.lp_;

        // Highs.cpp:5068-5076
        if self.model_lp_.is_mip() {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot apply crossover to solve MIP\n"
            );
            return_status = HighsStatus::kError;
        } else if self.model_hessian_.dim_ != 0 {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot apply crossover to solve QP\n"
            );
            return_status = HighsStatus::kError;
        } else {
            // Highs.cpp:5077 — clearSolver() returns HighsStatus; C++ discards it.
            self.clear_solver();

            // Highs.cpp:5078
            self.solution_ = user_solution.clone();

            // Highs.cpp:5080-5081
            return_status = crate::presolve::icrash_x::call_crossover(
                &self.options_,
                &self.model_lp_,
                &mut self.basis_,
                &mut self.solution_,
                &mut self.model_status_,
                &mut self.info_,
                &mut self.callback_,
            );

            // Highs.cpp:5082
            if return_status == HighsStatus::kError {
                return return_status;
            }

            // Highs.cpp:5084-5085
            self.info_.objective_function_value = self
                .model_lp_
                .objective_value(&self.solution_.col_value, self.options_.h16_obj_cdouble);

            // Highs.cpp:5086
            crate::lp_data::presolve_support::get_lp_kkt_failures(
                &self.options_,
                &self.model_lp_,
                &self.solution_,
                &mut self.info_,
                false,
            );
        }

        // Highs.cpp:5088
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.h::Highs::getPresolvedLp (line 504)
    pub fn get_presolved_lp(&self) -> &HighsLp {
        &self.presolved_model_lp_
    }

    /// HiGHS: Highs.h::Highs::getModelPresolveStatus (lines 592-594)
    pub fn get_model_presolve_status(&self) -> HighsPresolveStatus {
        self.model_presolve_status_
    }
}

/// HiGHS: HighsSolution.cpp::computeObjectiveValue (lines 1413-1419).
/// Plain (not quad) accumulation, offset added last.
pub fn compute_objective_value(lp: &HighsLp, solution: &HighsSolution) -> f64 {
    // HighsSolution.cpp:1414
    let mut objective_value = 0.0;

    // HighsSolution.cpp:1415-1416
    for i_col in 0..lp.num_col_ {
        objective_value += lp.col_cost_[i_col as usize] * solution.col_value[i_col as usize];
    }

    // HighsSolution.cpp:1417
    objective_value += lp.offset_;

    // HighsSolution.cpp:1418
    objective_value
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn objective_value_adds_offset_last() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.col_cost_ = vec![2.0, -1.0];
        lp.offset_ = 0.5;
        let mut solution = HighsSolution::default();
        solution.col_value = vec![1.5, 4.0];
        assert_eq!(compute_objective_value(&lp, &solution), -0.5);
    }

    #[test]
    fn presolve_rejects_models_needing_mods() {
        let mut highs = Highs::default();
        highs.model_lp_.num_col_ = 1;
        highs.model_lp_.col_cost_ = vec![1.0];
        highs.model_lp_.col_lower_ = vec![1.0];
        highs.model_lp_.col_upper_ = vec![2.0];
        highs.model_lp_.integrality_ = vec![HighsVarType::kSemiContinuous];
        highs.model_lp_.a_matrix_.start_ = vec![0, 0];
        assert_eq!(highs.presolve(), HighsStatus::kError);
    }

    #[test]
    fn crossover_rejects_a_mip() {
        let mut highs = Highs::default();
        highs.model_lp_.num_col_ = 1;
        highs.model_lp_.col_cost_ = vec![1.0];
        highs.model_lp_.col_lower_ = vec![0.0];
        highs.model_lp_.col_upper_ = vec![2.0];
        highs.model_lp_.integrality_ = vec![HighsVarType::kInteger];
        highs.model_lp_.a_matrix_.num_col_ = 1;
        highs.model_lp_.a_matrix_.start_ = vec![0, 0];
        assert_eq!(
            highs.crossover(&HighsSolution::default()),
            HighsStatus::kError
        );
    }

    #[test]
    fn presolve_getters_read_the_fields() {
        let highs = Highs::default();
        assert_eq!(
            highs.get_model_presolve_status(),
            highs.model_presolve_status_
        );
        assert_eq!(highs.get_presolved_lp().num_col_, 0);
    }
}
