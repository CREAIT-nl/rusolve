//! HiGHS source: lp_data/Highs.cpp (the presolve branch of calledOptimizeModel,
//! runPresolve, runPostsolve, clearPresolve,
//! presolveStatusToString, the `solveLp` lambda of calledOptimizeModel)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The LP presolve path of `Highs::run`. `run` (highs_run.rs) keeps the
//! C++ gate verbatim:
//!
//! ```text
//! if ((unconstrained_lp || has_basis || without_presolve) && solver_will_use_basis)
//!     <solve the incumbent LP: unchanged, bit-exact path>
//! else
//!     called_optimize_model_presolve(...)   // Highs.cpp 1554-1983
//! if (!no_incumbent_lp_solution_or_basis) { callLpKktCheck(model_.lp_); info_.valid = true; }
//! ```
//!
//! # Translation rules for this file
//!
//! * `incumbent_lp` (`model_.lp_`) is `self.model_lp_`; `reduced_lp`
//!   (`presolve_.getReducedProblem()`) is `self.presolve_.data_.reduced_lp_`.
//!   Because `callSolveLp(HighsLp& lp, ...)` borrows `lp` next to the other
//!   `Highs` members, the LP is named by [`HighsSolveTarget`]:
//!   `solveLp(incumbent_lp, msg, t)` -> `self.solve_lp_lambda(HighsSolveTarget::Incumbent, msg, &mut t)`,
//!   `solveLp(reduced_lp, msg, t)` -> `self.solve_lp_lambda(HighsSolveTarget::Reduced, msg, &mut t)`.
//!   `callSolveLp(incumbent_lp, msg)` -> `self.call_solve_lp_target(HighsSolveTarget::Incumbent, msg)`.
//! * The `solveLp` lambda assigns the enclosing `call_status`; here it
//!   returns it: `call_status = self.solve_lp_lambda(...)`.
//! * `return returnFromOptimizeModel(s, undo_mods)` inside
//!   `called_optimize_model_presolve` -> `return Err(self.return_from_optimize_model(s, undo_mods))`;
//!   a bare `return HighsStatus::kError` (after `cleanBounds`) -> `return Err(HighsStatus::kError)`.
//!   Falling out of the `else` branch -> `Ok(no_incumbent_lp_solution_or_basis)`.
//!   `undo_mods` is the enclosing `calledOptimizeModel` local, passed in.
//! * `presolve_.options_ = &options_` -> `self.presolve_.options_ = self.options_.clone()`
//!   (see presolve_component.rs). `presolve_.init(original_lp, timer_)` ->
//!   `self.presolve_.init(&self.model_lp_, &self.timer_, false)`.
//! * `model_.isEmpty()` -> `self.model_lp_.num_col_ == 0 && self.model_lp_.num_row_ == 0`;
//!   `model_.isMip()` -> `self.model_lp_.is_mip()` (always false from `run`, which
//!   routes a MIP to `call_solve_mip`; the MIP branch of `runPresolve` is
//!   `run_presolve_mip` and the `isMip` branch of `callRunPostsolve` is
//!   `call_run_postsolve_mip`, reached from the public
//!   `presolve` / `postsolve`).
//! * Dropped, each a no-op in this port: `timer_.start/stop/read(timer_.*_clock)`
//!   (sub-clocks; use 0.0 for the `this_*_time` reads), the `run_data_` time fields
//!   (highs_run_data.rs),
//!   `possibly_use_log_dev_level_2`, `kAllowDeveloperAssert` (false in the
//!   reference build: take the `assessSmallValues` branch), the
//!   `consider_pdlp_cleanup = false` block and its `presolved_lp_info` /
//!   `presolved_lp_pdlp_iteration_count` copies, profiling in `postsolve`,
//!   the final timing `highsLogDev` report (Highs.cpp 1994-2053: logging only).
//! * `callLpKktCheck(presolve_.getReducedProblem(), "Before postsolve")` ->
//!   `self.call_lp_kkt_check_target(HighsSolveTarget::Reduced)`.
//! * `HighsOptions save_options = options_; ... options_ = save_options;`
//!   -> `let save_options = self.options_.clone(); ... self.options_ = save_options;`.

#![allow(
    unused_variables,
    unused_mut,
    unused_imports,
    clippy::needless_late_init
)]

use crate::lp_data::highs_run::Highs;
use crate::lp_data::highs_solve::{highs_status_from_highs_model_status, solve_lp};
use crate::lp_data::presolve_support::*;
use crate::prelude::*;

/// Which LP a `callSolveLp(HighsLp& lp, ...)` / `solveLp(lp, ...)` call
/// solves (see the module doc).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum HighsSolveTarget {
    /// `model_.lp_` (`incumbent_lp`)
    Incumbent,
    /// `presolve_.data_.reduced_lp_` (`reduced_lp`)
    Reduced,
}

impl Highs {
    /// HiGHS: Highs.cpp::Highs::clearPresolve (lines 3903-3907, 5 lines)
    ///
    /// `presolved_model_.clear()` clears [`Highs::presolved_model_lp_`] (the port
    /// keeps no presolved Hessian).
    pub fn clear_presolve(&mut self) {
        self.model_presolve_status_ = HighsPresolveStatus::kNotPresolved;
        self.presolved_model_lp_.clear();
        self.presolve_.clear();
    }

    /// HiGHS: Highs.cpp::Highs::presolveStatusToString (lines 3687-3711, 25 lines)
    pub fn presolve_status_to_string(&self, presolve_status: HighsPresolveStatus) -> String {
        match presolve_status {
            HighsPresolveStatus::kNotPresolved => "Not presolved",
            HighsPresolveStatus::kNotReduced => "Not reduced",
            HighsPresolveStatus::kInfeasible => "Infeasible",
            HighsPresolveStatus::kUnboundedOrInfeasible => "Unbounded or infeasible",
            HighsPresolveStatus::kReduced => "Reduced",
            HighsPresolveStatus::kReducedToEmpty => "Reduced to empty",
            HighsPresolveStatus::kTimeout => "Timeout",
            HighsPresolveStatus::kOutOfMemory => "Memory allocation error",
            _ => {
                if cfg!(debug_assertions) {
                    panic!("presolveStatusToString: assert(1 == 0)");
                }
                "Unrecognised presolve status"
            }
        }
        .to_string()
    }

    /// HiGHS: Highs.cpp::Highs::callSolveLp(HighsLp& lp, message) (lines 4122-4137, 16 lines)
    /// for an explicit target LP. `Highs::call_solve_lp` (highs_run.rs) is the
    /// `Incumbent` case and stays the one the presolve-off path calls.
    pub fn call_solve_lp_target(&mut self, target: HighsSolveTarget, message: &str) -> HighsStatus {
        // Check that the model is column-wise
        debug_assert!(self.model_lp_.a_matrix_.is_colwise());
        let lp = match target {
            HighsSolveTarget::Incumbent => &mut self.model_lp_,
            HighsSolveTarget::Reduced => &mut self.presolve_.data_.reduced_lp_,
        };
        let mut solver_object = HighsLpSolverObject::new(
            lp,
            &mut self.basis_,
            &mut self.solution_,
            &mut self.info_,
            &mut self.callback_,
            &mut self.options_,
            &mut self.timer_,
        );
        // Solve the LP
        let return_status = solve_lp(&mut self.ekk_instance_, &mut solver_object, message);
        // Extract the model status
        self.model_status_ = solver_object.model_status_;
        return_status
    }

    /// HiGHS: Highs.cpp::Highs::calledOptimizeModel, the `solveLp` lambda
    /// (lines 1492-1507, 16 lines). Returns the `call_status` the lambda
    /// assigns. The solve sub-clock and the log-dev-level juggling are
    /// no-ops in this port, so `time` is left at 0.
    pub fn solve_lp_lambda(
        &mut self,
        target: HighsSolveTarget,
        lp_solve_description: &str,
        time: &mut f64,
    ) -> HighsStatus {
        *time = 0.0;
        self.call_solve_lp_target(target, lp_solve_description)
    }

    /// HiGHS: HighsInterface.cpp::Highs::callLpKktCheck(lp, message) (lines
    /// 2948-2951) for an explicit target LP -> HighsSolution.cpp::lpKktCheck
    /// (presolve_support.rs `lp_kkt_check`).
    pub fn call_lp_kkt_check_target(&mut self, target: HighsSolveTarget, message: &str) {
        let lp = match target {
            HighsSolveTarget::Incumbent => &self.model_lp_,
            HighsSolveTarget::Reduced => &self.presolve_.data_.reduced_lp_,
        };
        lp_kkt_check(
            &mut self.model_status_,
            &mut self.info_,
            lp,
            &self.solution_,
            &self.basis_,
            &self.options_,
            message,
        );
    }

    /// HiGHS: Highs.cpp::Highs::calledOptimizeModel, the presolve `else`
    /// branch (lines 1554-1983, 430 lines): presolve, solve the reduced LP,
    /// postsolve, and the clean-up re-solve of the incumbent LP from the
    /// postsolved basis.
    ///
    /// `return_status` is the enclosing local on entry; `lp_no_solution_basis`
    /// is the local of lines 1410-1411. Returns
    /// `Ok(no_incumbent_lp_solution_or_basis)` when control falls out of the
    /// branch (run() then does the shared KKT/`info_.valid` tail) or
    /// `Err(status)` for every early `return` (see the module doc).
    ///
    /// Determinism: the clean-up solve sets `options_.solver = kSimplexString`,
    /// `simplex_strategy = kSimplexStrategyChoose`, min/max concurrency 1 and
    /// (if > 0) `factor_pivot_threshold` from
    /// `ekk_instance_.info_.factor_pivot_threshold` of the presolved solve;
    /// `refine_basis` then `ekk_instance_.invalidate()` before it. The
    /// kReduced case calls `ekk_instance_.clear()` before solving the
    /// reduced LP and sets `options_.objective_bound = kHighsInf` around it.
    pub fn called_optimize_model_presolve(
        &mut self,
        return_status: HighsStatus,
        lp_no_solution_basis: bool,
        undo_mods: bool,
    ) -> Result<bool, HighsStatus> {
        let mut return_status = return_status;
        // Enclosing local of calledOptimizeModel, false on entry (run() sets it so).
        let mut no_incumbent_lp_solution_or_basis = false;
        // Enclosing local assigned by the solveLp lambda. Every read below
        // follows an assignment in the same arm.
        let mut call_status: HighsStatus;
        let mut this_solve_original_lp_time: f64 = 0.0;
        let mut this_solve_presolved_lp_time: f64 = 0.0;

        // If using IPM to solve the reduced LP, but not crossover, set
        // lp_presolve_requires_basis_postsolve so that presolve can use
        // rules for which postsolve does not generate a basis.
        let lp_presolve_requires_basis_postsolve =
            self.options_.lp_presolve_requires_basis_postsolve;
        if lp_no_solution_basis {
            self.options_.lp_presolve_requires_basis_postsolve = false;
        }
        // Presolve, forced to LP presolve
        let force_lp_presolve = true;
        let this_presolve_time: f64 = 0.0;
        self.model_presolve_status_ = self.run_presolve(force_lp_presolve, false);
        // Recover any modified options
        self.options_.lp_presolve_requires_basis_postsolve = lp_presolve_requires_basis_postsolve;

        // Set an illegal local pivot threshold value that's updated after
        // solving the presolved LP - if simplex is used
        let mut factor_pivot_threshold: f64 = -1.0;

        // Log the presolve reductions
        report_presolve_reductions(
            &self.options_.log_options,
            self.model_presolve_status_,
            &self.model_lp_,
            self.presolve_.get_reduced_problem(),
        );
        self.run_data_.presolved_model_num_col = self.presolve_.get_reduced_problem().num_col_;
        self.run_data_.presolved_model_num_row = self.presolve_.get_reduced_problem().num_row_;
        self.run_data_.presolved_model_num_nz =
            self.presolve_.get_reduced_problem().a_matrix_.num_nz();
        match self.model_presolve_status_ {
            HighsPresolveStatus::kNotPresolved => {
                self.ekk_instance_.lp_name_ = "Original LP".to_string();
                call_status = self.solve_lp_lambda(
                    HighsSolveTarget::Incumbent,
                    "Not presolved: solving the LP",
                    &mut this_solve_original_lp_time,
                );
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "callSolveLp",
                );
                if return_status == HighsStatus::kError {
                    return Err(self.return_from_optimize_model(return_status, undo_mods));
                }
            }
            HighsPresolveStatus::kNotReduced => {
                self.ekk_instance_.lp_name_ = "Unreduced LP".to_string();
                // Log the presolve reductions
                call_status = self.solve_lp_lambda(
                    HighsSolveTarget::Incumbent,
                    "Problem not reduced by presolve: solving the LP",
                    &mut this_solve_original_lp_time,
                );
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "callSolveLp",
                );
                if return_status == HighsStatus::kError {
                    return Err(self.return_from_optimize_model(return_status, undo_mods));
                }
            }
            HighsPresolveStatus::kReduced => {
                self.presolve_.data_.reduced_lp_.origin_name_ = "Reduced LP".to_string();
                self.presolve_.data_.reduced_lp_.set_matrix_dimensions();
                // kAllowDeveloperAssert is false in the reference build, so
                // take the assessSmallValues branch.
                self.presolve_
                    .data_
                    .reduced_lp_
                    .a_matrix_
                    .assess_small_values(
                        &self.options_.log_options,
                        self.options_.small_matrix_value,
                    );
                call_status = clean_bounds(&self.options_, &mut self.presolve_.data_.reduced_lp_);
                // Ignore any warning from clean bounds since the original LP
                // is still solved after presolve
                if interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "cleanBounds",
                ) == HighsStatus::kError
                {
                    return Err(HighsStatus::kError);
                }
                // Solving the presolved LP with strictly reduced dimensions
                // so ensure that the Ekk instance is cleared
                self.ekk_instance_.clear();
                self.ekk_instance_.lp_name_ = "Presolved LP".to_string();
                // Don't try dual cut-off when solving the presolved LP, as the
                // objective values aren't correct
                let save_objective_bound = self.options_.objective_bound;
                self.options_.objective_bound = kHighsInf;
                call_status = self.solve_lp_lambda(
                    HighsSolveTarget::Reduced,
                    "Solving the presolved LP",
                    &mut this_solve_presolved_lp_time,
                );
                if self.ekk_instance_.status_.initialised_for_solve {
                    // Record the pivot threshold resulting from solving the presolved LP
                    // with simplex
                    factor_pivot_threshold = self.ekk_instance_.info_.factor_pivot_threshold;
                }
                // Restore the dual objective cut-off
                self.options_.objective_bound = save_objective_bound;
                return_status = interpret_call_status(
                    &self.options_.log_options,
                    call_status,
                    return_status,
                    "callSolveLp",
                );
                if return_status == HighsStatus::kError {
                    return Err(self.return_from_optimize_model(return_status, undo_mods));
                }
                no_incumbent_lp_solution_or_basis = matches!(
                    self.model_status_,
                    HighsModelStatus::kInfeasible
                        | HighsModelStatus::kUnbounded
                        | HighsModelStatus::kUnboundedOrInfeasible
                        | HighsModelStatus::kTimeLimit
                        | HighsModelStatus::kIterationLimit
                        | HighsModelStatus::kInterrupt
                );
                if no_incumbent_lp_solution_or_basis {
                    // Postsolve won't be performed, so clear the HEkk data
                    // corresponding to the (strictly reduced) presolved LP here
                    self.ekk_instance_.clear();
                    self.set_highs_model_status_and_clear_solution_and_basis(self.model_status_);
                }
            }
            HighsPresolveStatus::kReducedToEmpty => {
                // Create a trivial optimal solution for postsolve to use
                self.solution_.clear();
                self.basis_.clear();
                self.basis_.debug_origin_name = "Presolve to empty".to_string();
                self.basis_.valid = true;
                self.basis_.alien = false;
                self.basis_.useful = true;
                self.basis_.was_alien = false;
                self.solution_.value_valid = true;
                self.solution_.dual_valid = true;
                // Optimality will not be spotted if there's no basis
                // postsolve
                self.model_status_ = HighsModelStatus::kOptimal;
            }
            HighsPresolveStatus::kInfeasible => {
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kInfeasible,
                );
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Problem status detected on presolve: {}\n",
                    util_model_status_to_string(self.model_status_)
                );
                return Err(self.return_from_optimize_model(return_status, undo_mods));
            }
            HighsPresolveStatus::kUnboundedOrInfeasible => {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Problem status detected on presolve: {}\n",
                    util_model_status_to_string(HighsModelStatus::kUnboundedOrInfeasible)
                );
                if self.options_.allow_unbounded_or_infeasible {
                    self.set_highs_model_status_and_clear_solution_and_basis(
                        HighsModelStatus::kUnboundedOrInfeasible,
                    );
                    return Err(self.return_from_optimize_model(return_status, undo_mods));
                }
                // Presolve has returned kUnboundedOrInfeasible, but HiGHS
                // can't return this. Use primal simplex solver on the original
                // LP
                // HiGHS 1.16 #3421 (h16_unresolved_unbounded_or_infeasible): Highs.cpp line 1902 (1.16.0)
                let save_return_status = return_status;
                let save_options = self.options_.clone();
                self.options_.solver = kSimplexString.to_string();
                self.options_.simplex_strategy = kSimplexStrategyPrimal;
                // The solveLp result is not read here (the C++ return_status
                // is not updated after this solve).
                self.solve_lp_lambda(
                    HighsSolveTarget::Incumbent,
                    "Solving the original LP with primal simplex to determine infeasible or unbounded",
                    &mut this_solve_original_lp_time,
                );
                // Recover the options
                self.options_ = save_options;
                // HiGHS 1.16 #3421 (h16_unresolved_unbounded_or_infeasible): Highs.cpp lines 1912-1936 (1.16.0)
                if self.options_.h16_unresolved_unbounded_or_infeasible {
                    // Primal simplex doesn't necessarily resolve the
                    // kUnboundedOrInfeasible case, as it may time out or reach
                    // the iteration count (#3383), or even fail
                    debug_assert!(
                        self.model_status_ == HighsModelStatus::kInfeasible
                            || self.model_status_ == HighsModelStatus::kUnbounded
                            || self.model_status_ == HighsModelStatus::kTimeLimit
                            || self.model_status_ == HighsModelStatus::kIterationLimit
                            || self.model_status_ == HighsModelStatus::kInterrupt
                    );
                    let resolved = self.model_status_ == HighsModelStatus::kInfeasible
                        || self.model_status_ == HighsModelStatus::kUnbounded;
                    if resolved {
                        self.info_.valid = true;
                    } else {
                        // Unresolved, so return with result from presolve
                        crate::highs_log_user!(
                            &self.options_.log_options,
                            HighsLogType::kInfo,
                            "Unable to resolve whether LP is unbounded or infeasible: primal simplex yields model status of {}\n",
                            util_model_status_to_string(self.model_status_)
                        );
                        self.model_status_ = HighsModelStatus::kUnboundedOrInfeasible;
                        return_status = save_return_status;
                        self.set_highs_model_status_and_clear_solution_and_basis(
                            HighsModelStatus::kUnboundedOrInfeasible,
                        );
                    }
                    return Err(self.return_from_optimize_model(return_status, undo_mods));
                }
                if return_status == HighsStatus::kError {
                    return Err(self.return_from_optimize_model(return_status, undo_mods));
                }
                self.info_.valid = true;
                debug_assert!(
                    self.model_status_ == HighsModelStatus::kInfeasible
                        || self.model_status_ == HighsModelStatus::kUnbounded
                );
                return Err(self.return_from_optimize_model(return_status, undo_mods));
            }
            HighsPresolveStatus::kTimeout => {
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kTimeLimit,
                );
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kWarning,
                    "Presolve reached timeout\n"
                );
                return Err(self.return_from_optimize_model(HighsStatus::kWarning, undo_mods));
            }
            HighsPresolveStatus::kOutOfMemory => {
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kMemoryLimit,
                );
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Presolve fails due to memory allocation error\n"
                );
                return Err(self.return_from_optimize_model(HighsStatus::kError, undo_mods));
            }
            _ => {
                debug_assert!(self.model_presolve_status_ == HighsPresolveStatus::kNullError);
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kPresolveError,
                );
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Presolve returned status {}\n",
                    self.model_presolve_status_ as i32
                );
                return Err(self.return_from_optimize_model(HighsStatus::kError, undo_mods));
            }
        }
        // End of presolve
        //
        // Cases of infeasibility/unboundedness timeout and memory errors
        // all handled, so just the successes remain
        debug_assert!(matches!(
            self.model_presolve_status_,
            HighsPresolveStatus::kNotPresolved
                | HighsPresolveStatus::kNotReduced
                | HighsPresolveStatus::kReduced
                | HighsPresolveStatus::kReducedToEmpty
        ));

        // Postsolve. Does nothing if there were no reductions during presolve.

        // If presolve has been run assuming that there's no basis
        // postsolve - allowing sparsify to be used in presolve - so
        // invalidate any basis
        if lp_no_solution_basis {
            self.invalidate_basis();
        }
        let have_optimal_reduced_solution = self.model_presolve_status_
            == HighsPresolveStatus::kReducedToEmpty
            || (self.model_presolve_status_ == HighsPresolveStatus::kReduced
                && self.model_status_ == HighsModelStatus::kOptimal);
        let have_unknown_reduced_solution = self.model_presolve_status_
            == HighsPresolveStatus::kReduced
            && self.model_status_ == HighsModelStatus::kUnknown;
        if have_optimal_reduced_solution || have_unknown_reduced_solution {
            debug_assert!(
                self.model_status_ == HighsModelStatus::kOptimal
                    || self.model_status_ == HighsModelStatus::kUnknown
                    || self.model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty
            );
            debug_assert!(
                self.model_presolve_status_ == HighsPresolveStatus::kReduced
                    || self.model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty
            );
            if have_unknown_reduced_solution {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kWarning,
                    "Running postsolve on non-optimal solution of reduced LP\n\n"
                );
            }
            // If presolve is nontrivial, extract the optimal solution
            // and (any) basis for the presolved problem in order to generate
            // the solution and basis for postsolve to use to generate a
            // solution and basis that is, hopefully, optimal.
            self.presolve_.data_.recovered_solution_ = self.solution_.clone();
            self.presolve_.data_.recovered_basis_ = self.basis_.clone();

            if self.model_presolve_status_ == HighsPresolveStatus::kReduced {
                self.call_lp_kkt_check_target(HighsSolveTarget::Reduced, "Before postsolve");
            }

            // Timing is not ported: the postsolve clock reads are logging only.
            let this_postsolve_time: f64 = 0.0;
            let postsolve_status = self.run_postsolve();
            self.presolve_.info_.presolve_time = this_presolve_time;
            self.presolve_.info_.postsolve_time = this_postsolve_time;

            if postsolve_status == HighsPostsolveStatus::kSolutionRecovered {
                // Indicate that nontrivial postsolve has been performed
                if self.model_presolve_status_ == HighsPresolveStatus::kReduced
                    || self.model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty
                {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "Performed postsolve\n"
                    );
                }
                // Set solution and its status
                self.solution_.clear();
                self.solution_ = self.presolve_.data_.recovered_solution_.clone();
                self.solution_.value_valid = true;
                if !self.basis_.valid {
                    // Have a primal-dual solution, but no basis, since PDLP or
                    // IPX without crossover were used. In the case of IPX, this
                    // was because either run_crossover was "off" or "choose"
                    // and IPX determined optimality
                    self.solution_.dual_valid = true;
                    self.invalidate_basis();
                } else {
                    //
                    // Hot-start the simplex solver for the incumbent LP
                    //
                    self.solution_.dual_valid = true;
                    // Set basis and its status
                    self.basis_.valid = true;
                    self.basis_.useful = true;
                    self.basis_.col_status =
                        self.presolve_.data_.recovered_basis_.col_status.clone();
                    self.basis_.row_status =
                        self.presolve_.data_.recovered_basis_.row_status.clone();
                    self.basis_.debug_origin_name += ": after postsolve";
                    // Basic primal activities are wrong after postsolve, so
                    // possibly skip KKT check
                    let perform_kkt_check = true;
                    if perform_kkt_check {
                        // Possibly force debug to perform KKT check on what's
                        // returned from postsolve
                        let force_debug = false;
                        let save_highs_debug_level = self.options_.highs_debug_level;
                        if force_debug {
                            self.options_.highs_debug_level = kHighsDebugLevelCostly;
                        }
                        if debug_highs_solution(
                            "After returning from postsolve",
                            &self.options_,
                            &self.model_lp_,
                            &self.solution_,
                            &self.basis_,
                        ) == HighsDebugStatus::kLogicalError
                        {
                            return Err(
                                self.return_from_optimize_model(HighsStatus::kError, undo_mods)
                            );
                        }
                        self.options_.highs_debug_level = save_highs_debug_level;
                    }
                    // Save the options to allow the best simplex strategy to
                    // be used. full_logging is false in the reference build,
                    // so the log_dev_level override is not applied.
                    let save_options = self.options_.clone();
                    // Force the use of simplex to clean up
                    self.options_.solver = kSimplexString.to_string();
                    self.options_.simplex_strategy = kSimplexStrategyChoose;
                    // Ensure that the parallel solver isn't used
                    self.options_.simplex_min_concurrency = 1;
                    self.options_.simplex_max_concurrency = 1;
                    // Use any pivot threshold resulting from solving the presolved LP
                    if factor_pivot_threshold > 0.0 {
                        self.options_.factor_pivot_threshold = factor_pivot_threshold;
                    }
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
                    // Set up the iteration count so that adding the value after
                    // callSolveLp gives the difference
                    let mut postsolve_iteration_count = -self.info_.simplex_iteration_count;
                    call_status = self.solve_lp_lambda(
                        HighsSolveTarget::Incumbent,
                        "Solving the original LP from the solution after postsolve",
                        &mut this_solve_original_lp_time,
                    );
                    restore_tiny_range_cols(&mut self.model_lp_, &fixed_tiny);
                    // Determine the iteration count
                    postsolve_iteration_count += self.info_.simplex_iteration_count;
                    //
                    return_status = HighsStatus::kOk;
                    return_status = interpret_call_status(
                        &self.options_.log_options,
                        call_status,
                        return_status,
                        "callSolveLp",
                    );
                    // Recover the options
                    self.options_ = save_options;
                    if return_status == HighsStatus::kError {
                        return Err(self.return_from_optimize_model(return_status, undo_mods));
                    }
                    self.run_data_.num_simplex_iterations_after_postsolve =
                        postsolve_iteration_count;
                    if postsolve_iteration_count > 0 {
                        crate::highs_log_user!(
                            &self.options_.log_options,
                            HighsLogType::kInfo,
                            "Required {} simplex iterations after postsolve\n",
                            postsolve_iteration_count
                        );
                    }
                }
            } else {
                // HiGHS 1.16 #3350 (h16_postsolve_basis_check): Highs.cpp lines 2123-2125 (1.16.0)
                if self.options_.h16_postsolve_basis_check {
                    // Postsolve has failed, so any simplex data for the presolved
                    // problem is useless
                    self.ekk_instance_.clear();
                }
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Postsolve return status is {}\n",
                    postsolve_status as i32
                );
                self.set_highs_model_status_and_clear_solution_and_basis(
                    HighsModelStatus::kPostsolveError,
                );
                return Err(self.return_from_optimize_model(HighsStatus::kError, undo_mods));
            }
        }
        // The PDLP clean-up block is not ported: consider_pdlp_cleanup is
        // false in the reference build. The shared KKT/info_.valid tail is
        // done by `run`, which reads no_incumbent_lp_solution_or_basis.
        Ok(no_incumbent_lp_solution_or_basis)
    }

    /// HiGHS: Highs.cpp::Highs::runPresolve (lines 3750-3868, 119 lines)
    ///
    /// C++ default arguments: `force_presolve = false`.
    ///
    /// The `model_.isMip() && !force_lp_presolve` branch (HighsMipSolver
    /// presolve) is `run_presolve_mip`. `presolve_log_ = presolve_.getPresolveLog()`
    /// is a clone.
    pub fn run_presolve(
        &mut self,
        force_lp_presolve: bool,
        force_presolve: bool,
    ) -> HighsPresolveStatus {
        self.presolve_.clear();
        // Exit if presolve is set to off (unless presolve is forced)
        if self.options_.presolve == kHighsOffString && !force_presolve {
            return HighsPresolveStatus::kNotPresolved;
        }

        if self.model_lp_.num_col_ == 0 && self.model_lp_.num_row_ == 0 {
            // Empty models shouldn't reach here, but this status would cause
            // no harm if one did
            if cfg!(debug_assertions) {
                panic!("runPresolve: assert(1 == 0)");
            }
            return HighsPresolveStatus::kNotReduced;
        }

        // Ensure that the LP is column-wise
        self.model_lp_.ensure_colwise();

        if self.model_lp_.num_col_ == 0 && self.model_lp_.num_row_ == 0 {
            return HighsPresolveStatus::kNullError;
        }

        // Ensure that the timer is running
        if !self.timer_.running() {
            self.timer_.start();
        }
        let start_presolve = self.timer_.read();

        // Set time limit.
        if self.options_.time_limit > 0.0 && self.options_.time_limit < kHighsInf {
            let left = self.options_.time_limit - start_presolve;
            if left <= 0.0 {
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Time limit reached while reading in matrix\n"
                );
                return HighsPresolveStatus::kTimeout;
            }

            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kVerbose,
                "Time limit set: reading matrix took {}, presolve time left: {}\n",
                start_presolve,
                left
            );
        }

        // Presolve.
        let presolve_return_status: HighsPresolveStatus;
        if self.model_lp_.is_mip() && !force_lp_presolve {
            // Use presolve for MIP (Highs.cpp 3793-3811)
            presolve_return_status = self.run_presolve_mip();
        } else {
            // Use presolve for LP
            self.presolve_.init(&self.model_lp_, &self.timer_, false);
            self.presolve_.options_ = self.options_.clone();
            if self.options_.time_limit > 0.0 && self.options_.time_limit < kHighsInf {
                let current = self.timer_.read();
                let time_init = current - start_presolve;
                let left = self.presolve_.options_.time_limit - time_init;
                if left <= 0.0 {
                    crate::highs_log_dev!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Time limit reached while copying matrix into presolve.\n"
                    );
                    return HighsPresolveStatus::kTimeout;
                }
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kVerbose,
                    "Time limit set: copying matrix took {}, presolve time left: {}\n",
                    time_init,
                    left
                );
            }

            presolve_return_status = self.presolve_.run();
        }

        let presolve_status_string = self.presolve_status_to_string(presolve_return_status);
        crate::highs_log_dev!(
            &self.options_.log_options,
            HighsLogType::kVerbose,
            "presolve_.run() returns status: {}\n",
            presolve_status_string
        );

        // Update reduction counts.
        debug_assert!(presolve_return_status == self.presolve_.presolve_status_);
        self.presolve_log_ = self.presolve_.get_presolve_log().clone();
        match self.presolve_.presolve_status_ {
            HighsPresolveStatus::kReduced => {
                self.presolve_.info_.n_cols_removed =
                    self.model_lp_.num_col_ - self.presolve_.data_.reduced_lp_.num_col_;
                self.presolve_.info_.n_rows_removed =
                    self.model_lp_.num_row_ - self.presolve_.data_.reduced_lp_.num_row_;
                self.presolve_.info_.n_nnz_removed = self.model_lp_.a_matrix_.num_nz()
                    - self.presolve_.data_.reduced_lp_.a_matrix_.num_nz();
                // Clear any scaling information inherited by the reduced LP
                self.presolve_.data_.reduced_lp_.clear_scale();
                debug_assert!(crate::lp_data::lp_assess::lp_dimensions_ok(
                    "RunPresolve: reduced_lp",
                    &self.presolve_.data_.reduced_lp_,
                    &self.options_.log_options
                ));
            }
            HighsPresolveStatus::kReducedToEmpty => {
                self.presolve_.info_.n_cols_removed = self.model_lp_.num_col_;
                self.presolve_.info_.n_rows_removed = self.model_lp_.num_row_;
                self.presolve_.info_.n_nnz_removed = self.model_lp_.a_matrix_.num_nz();
            }
            _ => {}
        }
        // Presolve creates integrality vector for an LP, so clear it
        if !self.model_lp_.is_mip() {
            self.presolve_.data_.reduced_lp_.integrality_.clear();
        }

        presolve_return_status
    }

    /// HiGHS: Highs.cpp::Highs::runPostsolve (lines 3870-3895, 26 lines)
    /// HiGHS 1.16 (#3248): Highs.cpp::runPostsolve (lines 4105-4147)
    ///
    /// HiGHS 1.16 (#3350): Highs.cpp::runPostsolve (lines 4105-4147); the basis
    /// consistency check after undo() is behind `h16_postsolve_basis_check`.
    ///
    /// `postSolveStack.undo(options_, recovered_solution_, recovered_basis_)`
    /// borrows three fields of `self.presolve_.data_` and `self.options_`
    /// disjointly: `self.presolve_.data_.post_solve_stack.undo(&self.options_,
    /// &mut self.presolve_.data_.recovered_solution_, &mut self.presolve_.data_.recovered_basis_)`.
    /// `calculateRowValuesQuad(model_.lp_, recovered_solution_)` ->
    /// `calculate_row_values_quad(&self.model_lp_, &mut self.presolve_.data_.recovered_solution_, -1)`.
    pub fn run_postsolve(&mut self) -> HighsPostsolveStatus {
        // assert(presolve_.has_run_) is commented out in the C++ source.
        let have_primal_solution = self.presolve_.data_.recovered_solution_.value_valid;
        // Need at least a primal solution
        if !have_primal_solution {
            return HighsPostsolveStatus::kNoPrimalSolutionError;
        }
        let have_dual_solution = self.presolve_.data_.recovered_solution_.dual_valid;
        // HiGHS 1.16 #3350 (h16_postsolve_basis_check): Highs.cpp lines 4115, 4119-4128 (1.16.0)
        let h16_postsolve_basis_check = self.options_.h16_postsolve_basis_check;
        let have_basis = if h16_postsolve_basis_check {
            self.presolve_.data_.recovered_basis_.valid
        } else {
            false
        };
        // undo() takes the C++ default arguments report_col = -1, thread_safe = false.
        self.presolve_.data_.post_solve_stack.undo(
            &self.options_,
            &mut self.presolve_.data_.recovered_solution_,
            &mut self.presolve_.data_.recovered_basis_,
            -1,
            false,
        );
        // HiGHS 1.16 #3350 (h16_postsolve_basis_check): Highs.cpp lines 4119-4128 (1.16.0)
        if h16_postsolve_basis_check
            && have_basis
            && !is_basis_consistent(&self.model_lp_, &self.presolve_.data_.recovered_basis_)
        {
            // Recovered basis was assumed to be consistent, but #3323
            // exposes that in exceptional circumstances it may not be,
            // so check here and return solver error if it is
            // inconsistent
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Highs::runPostsolve: Error in basis after postsolve\n"
            );
            return HighsPostsolveStatus::kBasisError;
        }
        // Compute the row activities
        debug_assert!(self.model_lp_.a_matrix_.is_colwise());
        calculate_row_values_quad(
            &self.model_lp_,
            &mut self.presolve_.data_.recovered_solution_,
            -1,
        );

        // HiGHS 1.16 #3248 (h16_max_duals): Highs.cpp lines 4133-4140 (1.16.0)
        if have_dual_solution && self.model_lp_.sense_ == ObjSense::kMaximize {
            if self.options_.h16_max_duals {
                // Negate the dual values since the incumbent LP is a maximization
                for i_col in 0..self.model_lp_.num_col_ {
                    self.presolve_.data_.recovered_solution_.col_dual[i_col as usize] *= -1.0;
                }
                for i_row in 0..self.model_lp_.num_row_ {
                    self.presolve_.data_.recovered_solution_.row_dual[i_row as usize] *= -1.0;
                }
            } else {
                self.presolve_.negate_reduced_lp_col_duals();
            }
        }

        // Ensure that the postsolve status is used to set
        // presolve_.postsolve_status_, as well as being returned
        let postsolve_status = HighsPostsolveStatus::kSolutionRecovered;
        self.presolve_.postsolve_status_ = postsolve_status;
        postsolve_status
    }
}

#[cfg(test)]
mod h16_max_duals_tests {
    //! HiGHS 1.16 #3248 tests: the blending LP of `TestPresolve.cpp` (max-lp-dual-postsolve),
    //! solved with IPM and no crossover, presolve on, for a maximization model.
    use super::*;

    /// The blending LP of 3248.diff: max 8x + 10y, 0.3x + 0.5y <= 120,
    /// 0.7x + 0.5y <= 210, x, y >= 0. Optimum x = 225, y = 105, objective 2850;
    /// both rows binding, so the row duals (max convention) are 15 and 5.
    fn blending_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![8.0, 10.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![120.0, 210.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![0.3, 0.7, 0.5, 0.5];
        lp.sense_ = ObjSense::kMaximize;
        lp.offset_ = 0.0;
        lp
    }

    /// Highs with the blending LP, IPM without crossover, presolve on or off,
    /// and `h16_max_duals` set to `h16_max_duals`. Not yet run.
    fn blending_highs(presolve_on: bool, h16_max_duals: bool) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_max_duals", h16_max_duals),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", "ipm"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("run_crossover", "off"),
            HighsStatus::kOk
        );
        let presolve = if presolve_on { "on" } else { "off" };
        assert_eq!(
            highs.set_option_value_str("presolve", presolve),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(blending_lp()), HighsStatus::kOk);
        highs
    }

    /// Solve the blending LP and check it is optimal.
    fn solve_blending(presolve_on: bool, h16_max_duals: bool) -> Highs {
        let mut highs = blending_highs(presolve_on, h16_max_duals);
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        highs
    }

    #[test]
    fn h16_max_duals_on_matches_presolve_off() {
        // Presolve off never reaches run_postsolve, so it is the reference.
        let reference = solve_blending(false, false);
        let on = solve_blending(true, true);
        let ref_sol = reference.get_solution();
        let on_sol = on.get_solution();
        assert_eq!(on_sol.col_dual.len(), 2);
        assert_eq!(on_sol.row_dual.len(), 2);
        for i in 0..2 {
            assert!(
                (on_sol.col_dual[i] - ref_sol.col_dual[i]).abs() <= 1e-6,
                "col_dual[{i}]: on {} vs off {}",
                on_sol.col_dual[i],
                ref_sol.col_dual[i]
            );
            assert!(
                (on_sol.row_dual[i] - ref_sol.row_dual[i]).abs() <= 1e-6,
                "row_dual[{i}]: on {} vs off {}",
                on_sol.row_dual[i],
                ref_sol.row_dual[i]
            );
        }
        // The hand-derived duals of the blending LP (IPM tolerance).
        assert!((on_sol.row_dual[0] - 15.0).abs() <= 1e-6);
        assert!((on_sol.row_dual[1] - 5.0).abs() <= 1e-6);
    }

    #[test]
    fn h16_max_duals_on_redundant_row_keeps_duals() {
        // 3248.diff 5-70: solve, add a redundant row that presolve removes so that
        // postsolve runs, re-solve; the original duals must be unchanged.
        let mut highs = solve_blending(true, true);
        let first = highs.get_solution().clone();
        let optimal_objective = highs.get_info().objective_function_value;
        assert_eq!(
            highs.add_row(-kHighsInf, 1000.0, 2, &[0, 1], &[1.0, 1.0]),
            HighsStatus::kOk
        );
        highs.clear_solver();
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        let objective = highs.get_info().objective_function_value;
        assert!((optimal_objective - objective).abs() <= 1e-6 * optimal_objective.abs().max(1.0));
        let second = highs.get_solution();
        assert_eq!(second.row_dual.len(), 3);
        for i in 0..2 {
            assert!((first.col_dual[i] - second.col_dual[i]).abs() <= 1e-6);
            assert!((first.row_dual[i] - second.row_dual[i]).abs() <= 1e-6);
        }
        assert!(second.row_dual[2].abs() <= 1e-6);
    }

    #[test]
    fn h16_max_duals_off_is_1_15_1() {
        // Flag off keeps the 1.15.1 path: negateReducedLpColDuals (columns only) at
        // Highs.cpp:3886-3887. Presolve-on without the redundant row still solves
        // to the same duals as the reference.
        let first = solve_blending(true, false);
        assert!((first.get_solution().row_dual[0] - 15.0).abs() <= 1e-6);
        assert!((first.get_solution().row_dual[1] - 5.0).abs() <= 1e-6);
        // With the redundant row, the 1.15.1 postsolve does not reach an optimal
        // status (the flag-on test above does); the status is kept as the 1.15.1
        // result.
        let mut highs = first;
        assert_eq!(
            highs.add_row(-kHighsInf, 1000.0, 2, &[0, 1], &[1.0, 1.0]),
            HighsStatus::kOk
        );
        highs.clear_solver();
        assert_eq!(highs.run(), HighsStatus::kWarning);
    }
}

#[cfg(test)]
mod h16_pstat_tests {
    use super::*;

    /// Two columns, two rows, column-wise matrix with a unit diagonal.
    fn h16_pstat_tiny_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 2,
            col_cost_: vec![1.0, 1.0],
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![1.0, 1.0],
            row_lower_: vec![0.0, 0.0],
            row_upper_: vec![1.0, 1.0],
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        lp
    }

    /// Runs postsolve on an empty postsolve stack with the given recovered basis
    /// status for each of the two columns and two rows.
    fn h16_pstat_run_postsolve(
        flag: bool,
        col_status: [HighsBasisStatus; 2],
        row_status: [HighsBasisStatus; 2],
    ) -> HighsPostsolveStatus {
        let mut highs = Highs::default();
        highs.options_.h16_postsolve_basis_check = flag;
        highs.model_lp_ = h16_pstat_tiny_lp();
        // The original dimensions that presolve would have recorded.
        highs
            .presolve_
            .data_
            .post_solve_stack
            .initialize_index_maps(2, 2);
        let sol = &mut highs.presolve_.data_.recovered_solution_;
        sol.value_valid = true;
        sol.dual_valid = false;
        sol.col_value = vec![0.0, 0.0];
        sol.row_value = vec![0.0, 0.0];
        let basis = &mut highs.presolve_.data_.recovered_basis_;
        basis.valid = true;
        basis.col_status = col_status.to_vec();
        basis.row_status = row_status.to_vec();
        highs.run_postsolve()
    }

    /// HiGHS 1.16 #3350 (h16_postsolve_basis_check): an inconsistent recovered basis (four
    /// basic variables for two rows) gives kBasisError on the on path. Off keeps the
    /// 1.15.1 result kSolutionRecovered.
    #[test]
    fn h16_pstat_postsolve_basis_check_on_off() {
        let all_basic = [HighsBasisStatus::kBasic; 2];
        assert_eq!(
            h16_pstat_run_postsolve(true, all_basic, all_basic),
            HighsPostsolveStatus::kBasisError
        );
        assert_eq!(
            h16_pstat_run_postsolve(false, all_basic, all_basic),
            HighsPostsolveStatus::kSolutionRecovered
        );
    }

    /// HiGHS 1.16 #3350: a consistent basis (two basic variables for two rows) passes the
    /// check on the on path, and the off path agrees.
    #[test]
    fn h16_pstat_postsolve_basis_consistent_on() {
        let col_status = [HighsBasisStatus::kBasic, HighsBasisStatus::kNonbasic];
        let row_status = [HighsBasisStatus::kBasic, HighsBasisStatus::kNonbasic];
        assert_eq!(
            h16_pstat_run_postsolve(true, col_status, row_status),
            HighsPostsolveStatus::kSolutionRecovered
        );
        assert_eq!(
            h16_pstat_run_postsolve(false, col_status, row_status),
            HighsPostsolveStatus::kSolutionRecovered
        );
    }

    /// Feasible and unbounded LP (column 0 has cost -1 and no upper bound, and row 0
    /// is x - y >= 0) on which presolve returns kUnboundedOrInfeasible. The primal
    /// simplex solve on the original LP stops at the zero iteration limit.
    fn h16_pstat_unbounded_lp() -> HighsLp {
        let inf = kHighsInf;
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 1,
            col_cost_: vec![-1.0, 0.0],
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![inf, 1.0],
            row_lower_: vec![0.0],
            row_upper_: vec![inf],
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, -1.0];
        lp
    }

    /// HiGHS 1.16 #3421 (h16_unresolved_unbounded_or_infeasible): with presolve reporting
    /// kUnboundedOrInfeasible and no simplex iterations, the on path returns the
    /// presolve result without a panic. Off keeps the 1.15.1 result: the primal
    /// simplex model status kIterationLimit with info_.valid set. Release only: on
    /// both paths a C++ assert (debug_assert) fires, as it does in a debug HiGHS.
    #[test]
    #[cfg_attr(
        debug_assertions,
        ignore = "C++ asserts on allow_unbounded_or_infeasible fire in debug builds"
    )]
    fn h16_pstat_unresolved_unbounded_or_infeasible_on_off() {
        for flag in [true, false] {
            let mut highs = Highs::default();
            highs.options_.output_flag = false;
            highs.options_.presolve = kHighsOnString.to_string();
            highs.options_.allow_unbounded_or_infeasible = false;
            highs.options_.simplex_iteration_limit = 0;
            highs.options_.h16_unresolved_unbounded_or_infeasible = flag;
            assert_eq!(highs.pass_model(h16_pstat_unbounded_lp()), HighsStatus::kOk);
            let status = highs.run();
            if flag {
                assert_eq!(status, HighsStatus::kError);
                assert_eq!(
                    highs.get_model_status(),
                    HighsModelStatus::kUnboundedOrInfeasible
                );
            } else {
                assert_eq!(status, HighsStatus::kWarning);
                assert_eq!(highs.get_model_status(), HighsModelStatus::kIterationLimit);
                assert!(highs.info_.valid);
            }
        }
    }
}
