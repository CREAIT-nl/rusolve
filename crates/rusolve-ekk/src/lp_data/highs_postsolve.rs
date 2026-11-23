//! HiGHS source: lp_data/Highs.cpp (postsolve, callRunPostsolve)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The public `Highs::postsolve(solution, basis)`: postsolve a solution of
//! the presolved LP supplied by the caller. Not on the `Highs::run` path.
//! The translation rules of `highs_presolve.rs` (module doc) apply.

#![allow(unused_variables, unused_mut, unused_imports)]

use crate::lp_data::highs_presolve::HighsSolveTarget;
use crate::lp_data::highs_run::Highs;
use crate::lp_data::highs_solution::{
    is_col_dual_solution_right_size, is_row_dual_solution_right_size,
};
use crate::lp_data::highs_solve::highs_status_from_highs_model_status;
use crate::lp_data::presolve_support::*;
use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::postsolve(solution, basis) (lines 3606-3627, 22 lines)
    ///
    /// The one-argument overload (lines 3601-3604) is
    /// `postsolve(solution, &HighsBasis::default())`. Profiling is dropped.
    pub fn postsolve(&mut self, solution: &HighsSolution, basis: &HighsBasis) -> HighsStatus {
        let can_run_postsolve = self.model_presolve_status_ == HighsPresolveStatus::kNotPresolved
            || self.model_presolve_status_ == HighsPresolveStatus::kNotReduced
            || self.model_presolve_status_ == HighsPresolveStatus::kReduced
            || self.model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty
            || self.model_presolve_status_ == HighsPresolveStatus::kTimeout
            || self.model_presolve_status_ == HighsPresolveStatus::kOutOfMemory;
        if !can_run_postsolve {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "Cannot run postsolve with presolve status: {}\n",
                self.presolve_status_to_string(self.model_presolve_status_)
            );
            return HighsStatus::kWarning;
        }
        // Profiling is dropped (module doc).
        let return_status = self.call_run_postsolve(solution, basis);
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::callRunPostsolve (lines 4451-4673, 223 lines)
    ///
    /// Only reached from `postsolve`, never from `run`. The `isMip() &&
    /// !basis.valid` branch is `call_run_postsolve_mip`. `getKktFailures` is an
    /// out-of-scope stub (`get_kkt_failures` panics): port the call as
    /// written. `callSolveLp(incumbent_lp, msg)` ->
    /// `self.call_solve_lp_target(HighsSolveTarget::Incumbent, msg)`;
    /// `callLpKktCheck(model_.lp_)` -> `self.call_lp_kkt_check()`.
    /// `return returnFromOptimizeModel(s, undo_mods)` -> `return self.return_from_optimize_model(s, undo_mods)`.
    pub fn call_run_postsolve(
        &mut self,
        solution: &HighsSolution,
        basis: &HighsBasis,
    ) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;
        let mut call_status: HighsStatus;
        let presolved_lp = &self.presolve_.data_.reduced_lp_;

        // Must at least have a primal column solution of the right size
        if solution.col_value.len() as HighsInt != presolved_lp.num_col_ {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Primal solution provided to postsolve is of size {} rather than {}\n",
                solution.col_value.len() as HighsInt,
                presolved_lp.num_col_
            );
            return HighsStatus::kError;
        }
        // Check any basis that is supplied
        let basis_supplied =
            !basis.col_status.is_empty() || !basis.row_status.is_empty() || basis.valid;
        if basis_supplied {
            if !is_basis_consistent(presolved_lp, basis) {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Basis provided to postsolve is incorrect size or inconsistent\n"
                );
                return HighsStatus::kError;
            }
        }
        // Copy in the solution provided
        self.presolve_.data_.recovered_solution_ = solution.clone();
        // Ignore any row values
        self.presolve_.data_.recovered_solution_.row_value =
            vec![0.0; presolved_lp.num_row_ as usize];
        self.presolve_.data_.recovered_solution_.value_valid = true;

        if self.model_lp_.is_mip() && !basis.valid {
            // Postsolving a MIP without a valid basis (Highs.cpp 4483-4527)
            self.call_run_postsolve_mip();
        } else {
            // Postsolving an LP, or a MIP after solving the relaxation
            // (identified by passing a valid basis).
            //
            // If there are dual values, make sure that both vectors are the
            // right size
            let dual_supplied = self.presolve_.data_.recovered_solution_.col_dual.len() > 0
                || self.presolve_.data_.recovered_solution_.row_dual.len() > 0
                || self.presolve_.data_.recovered_solution_.dual_valid;
            if dual_supplied {
                if !is_row_dual_solution_right_size(
                    presolved_lp,
                    &self.presolve_.data_.recovered_solution_,
                ) {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Row dual solution provided to postsolve is of size {} rather than {}\n",
                        solution.row_dual.len() as HighsInt,
                        presolved_lp.num_row_
                    );
                    return HighsStatus::kError;
                }
                if !is_col_dual_solution_right_size(
                    presolved_lp,
                    &self.presolve_.data_.recovered_solution_,
                ) {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Column dual solution provided to postsolve is of size {} rather than {}\n",
                        solution.col_dual.len() as HighsInt,
                        presolved_lp.num_col_
                    );
                    return HighsStatus::kError;
                }
                self.presolve_.data_.recovered_solution_.dual_valid = true;
            } else {
                self.presolve_.data_.recovered_solution_.dual_valid = false;
            }
            // Copy in the basis provided. It's already been checked for
            // consistency, so the basis is valid iff it was supplied
            self.presolve_.data_.recovered_basis_ = basis.clone();
            self.presolve_.data_.recovered_basis_.valid = basis_supplied;

            let postsolve_status = self.run_postsolve();

            if postsolve_status == HighsPostsolveStatus::kSolutionRecovered {
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kVerbose,
                    "Postsolve finished\n"
                );
                // Set solution and its status
                self.solution_.clear();
                self.solution_ = self.presolve_.data_.recovered_solution_.clone();
                debug_assert!(self.solution_.value_valid);
                if !self.solution_.dual_valid {
                    self.solution_.col_dual = vec![0.0; self.model_lp_.num_col_ as usize];
                    self.solution_.row_dual = vec![0.0; self.model_lp_.num_row_ as usize];
                }
                self.basis_ = self.presolve_.data_.recovered_basis_.clone();
                // Validity of the solution and basis is inherited (the C++
                // assignments of valid/useful/status are commented out there).
                self.basis_.debug_origin_name.push_str(": after postsolve");
                if self.basis_.valid {
                    // Save the options to allow the best simplex strategy to be used
                    let save_options = self.options_.clone();
                    self.options_.simplex_strategy = kSimplexStrategyChoose;
                    // Ensure that the parallel solver isn't used
                    self.options_.simplex_min_concurrency = 1;
                    self.options_.simplex_max_concurrency = 1;
                    // The basis returned from postsolve is just basic/nonbasic
                    // and EKK expects a refined basis, so set it up now
                    refine_basis(&self.model_lp_, &self.solution_, &mut self.basis_);
                    if self.options_.lp_postsolve_snap_basis {
                        snap_postsolve_basis(&self.model_lp_, &self.solution_, &mut self.basis_);
                    }
                    let fixed_tiny = if self.options_.lp_postsolve_fix_tiny {
                        let tol = self.options_.primal_feasibility_tolerance;
                        fix_tiny_range_cols(&mut self.model_lp_, &mut self.basis_, tol)
                    } else {
                        Vec::new()
                    };
                    // Scrap the EKK data from solving the presolved LP
                    self.ekk_instance_.invalidate();
                    self.ekk_instance_.lp_name_ = "Postsolve LP".to_string();
                    call_status = self.call_solve_lp_target(
                        HighsSolveTarget::Incumbent,
                        "Solving the original LP from the solution after postsolve",
                    );
                    return_status = interpret_call_status(
                        &self.options_.log_options,
                        call_status,
                        return_status,
                        "callSolveLp",
                    );
                    restore_tiny_range_cols(&mut self.model_lp_, &fixed_tiny);
                    // Recover the options
                    self.options_ = save_options;
                    // model_.isQp() is false: this port has no Hessian.
                    let is_qp = false;
                    debug_assert!(!is_qp);
                    let get_residuals = true;
                    get_kkt_failures(
                        &self.options_,
                        is_qp,
                        &self.model_lp_,
                        &self.model_lp_.col_cost_,
                        &self.solution_,
                        &mut self.info_,
                        get_residuals,
                    );
                    if return_status == HighsStatus::kError {
                        // Set undo_mods = false, since passing models requiring
                        // modification to Highs::presolve is illegal
                        let undo_mods = false;
                        return self.return_from_optimize_model(return_status, undo_mods);
                    }
                } else {
                    self.basis_.clear();
                    // callLpKktCheck sets LPs with model status kUnknown to
                    // kOptimal if the tolerances are met (see module doc).
                    self.model_status_ = if self.model_lp_.is_mip() {
                        HighsModelStatus::kNotset
                    } else {
                        HighsModelStatus::kUnknown
                    };
                    self.call_lp_kkt_check();
                    self.info_.valid = true;
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "\nPure postsolve yields primal {} basis: model status is {}\n",
                        if self.solution_.dual_valid {
                            "and dual solution, but no"
                        } else {
                            "but no dual solution or"
                        },
                        util_model_status_to_string(self.model_status_)
                    );
                }
            } else {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Postsolve return status is {}\n",
                    postsolve_status as i32
                );
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kPostsolveError,
                );
                // Set undo_mods = false, since passing models requiring
                // modification to Highs::presolve is illegal
                let undo_mods = false;
                return self.return_from_optimize_model(HighsStatus::kError, undo_mods);
            }
        }
        call_status = highs_status_from_highs_model_status(self.model_status_);
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "highsStatusFromHighsModelStatus",
        );
        return_status
    }
}
