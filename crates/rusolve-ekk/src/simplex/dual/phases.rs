//! HiGHS source: simplex/HEkkDual.cpp (lines 573-1183)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Phase 1 / phase 2 loops, rebuild, cleanup.
//! See `HEkkDual` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkDual {
    /// HiGHS: HEkkDual.cpp::HEkkDual::solvePhase1 (lines 573-815, 243 lines)
    /// HiGHS 1.16 (#3353): HEkkDual.cpp::solvePhase1 (lines 612-614)
    pub fn solve_phase1(&mut self, ekk: &mut HEkk) {
        // When starting a new phase the (updated) dual objective function
        // value isn't known.
        ekk.status_.has_primal_objective_value = false;
        ekk.status_.has_dual_objective_value = false;
        // Set rebuild_reason so that it's assigned when first tested
        self.rebuild_reason = kRebuildReasonNo;
        debug_assert!(self.solve_phase == kSolvePhase1);
        debug_assert!(!ekk.solve_bailout_);
        if ekk.bailout() {
            return;
        }
        // Report the phase start
        crate::highs_log_dev!(
            &ekk.options_.log_options,
            HighsLogType::kDetailed,
            "dual-phase-1-start\n"
        );
        // Switch to dual phase 1 bounds
        ekk.initialise_bound(SimplexAlgorithm::kDual, self.solve_phase, false);
        ekk.initialise_nonbasic_value_and_move();
        // HiGHS 1.16 #3353 (h16_dual_free_list): HEkkDual.cpp lines 612-614 (1.16.0)
        if ekk.options_.h16_dual_free_list {
            // Free variables are boxed in phase 1, so the list of free variables
            // created at the start of any earlier phase 2 must not be used (#2300)
            self.dual_row.free_list.clear();
        }

        // If there's no backtracking basis, save the initial basis in case of
        // backtracking
        if !ekk.info_.valid_backtracking_basis_ {
            ekk.put_backtracking_basis();
        }

        // Main solving structure
        ekk.analysis_.simplex_timer_start(IterateClock, 0);
        loop {
            ekk.analysis_
                .simplex_timer_start(IterateDualRebuildClock, 0);
            self.rebuild(ekk);
            ekk.analysis_.simplex_timer_stop(IterateDualRebuildClock, 0);
            if self.solve_phase == kSolvePhaseError {
                ekk.model_status_ = HighsModelStatus::kSolveError;
                return;
            }
            if self.solve_phase == kSolvePhaseUnknown {
                // If backtracking, may change phase, so drop out
                ekk.analysis_.simplex_timer_stop(IterateClock, 0);
                return;
            }
            if ekk.bailout() {
                break;
            }
            loop {
                if self.debug_dual_simplex(ekk, "Before iteration", false)
                    == HighsDebugStatus::kLogicalError
                {
                    self.solve_phase = kSolvePhaseError;
                    return;
                }
                match ekk.info_.simplex_strategy {
                    kSimplexStrategyDualTasks => self.iterate_tasks(ekk),
                    kSimplexStrategyDualMulti => self.iterate_multi(ekk),
                    _ => self.iterate(ekk),
                }
                if ekk.bailout() {
                    break;
                }
                debug_assert!(self.solve_phase != kSolvePhaseTabooBasis);
                if self.rebuild_reason != 0 {
                    break;
                }
            }
            if ekk.solve_bailout_ {
                break;
            }
            // If the data are fresh from rebuild(), possibly break out of the
            // outer loop to see what's occurred
            let finished =
                ekk.status_.has_fresh_rebuild && !ekk.rebuild_refactor(self.rebuild_reason);
            if finished && ekk.taboo_bad_basis_change() {
                // A bad basis change has had to be made taboo without any other
                // basis changes having been performed from a fresh rebuild.
                self.solve_phase = kSolvePhaseTabooBasis;
                return;
            }
            if finished {
                break;
            }
        }
        ekk.analysis_.simplex_timer_stop(IterateClock, 0);
        // Possibly return due to bailing out, having now stopped IterateClock
        if ekk.solve_bailout_ {
            return;
        }

        // If bailing out, should have done so already
        debug_assert!(!ekk.solve_bailout_);
        // Assess outcome of dual phase 1
        if self.row_out == kNoRowChosen {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "dual-phase-1-optimal\n"
            );
            // Optimal in phase 1
            if ekk.info_.dual_objective_value == 0.0 {
                // Zero phase 1 objective so go to phase 2
                self.solve_phase = kSolvePhase2;
            } else {
                // A nonzero dual objective value at an optimal solution of
                // phase 1 means that there may be dual infeasibilities.
                // If assessPhase1Optimality() is happy that there are no dual
                // infeasibilities, it will set solve_phase = kSolvePhase2;
                self.assess_phase1_optimality(ekk);
            }
        } else if self.rebuild_reason == kRebuildReasonChooseColumnFail
            || self.rebuild_reason == kRebuildReasonExcessivePrimalValue
        {
            // chooseColumn has failed or excessive primal values have been created
            self.solve_phase = kSolvePhaseError;
            // Solve error is opaque to users, so put in some logging
            if self.rebuild_reason == kRebuildReasonChooseColumnFail {
                crate::highs_log_user!(
                    &ekk.options_.log_options,
                    HighsLogType::kError,
                    "Dual simplex ratio test failed due to excessive dual values: consider scaling down the LP objective coefficients\n"
                );
            } else {
                crate::highs_log_user!(
                    &ekk.options_.log_options,
                    HighsLogType::kError,
                    "Dual simplex detected excessive primal values: consider scaling down the LP bounds\n"
                );
            }
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "dual-phase-1-not-solved\n"
            );
            ekk.model_status_ = HighsModelStatus::kSolveError;
        } else if self.variable_in == -1 {
            // We got dual phase 1 unbounded - strange
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "dual-phase-1-unbounded\n"
            );
            if ekk.info_.costs_perturbed {
                // Clean up perturbation
                self.cleanup(ekk);
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "Cleaning up cost perturbation when unbounded in phase 1\n"
                );
                if self.dual_infeas_count == 0 {
                    // No dual infeasibilities and (since unbounded) at least zero
                    // phase 1 objective so go to phase 2
                    self.solve_phase = kSolvePhase2;
                }
            } else {
                // Report strange issues
                self.solve_phase = kSolvePhaseError;
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kInfo,
                    "dual-phase-1-not-solved\n"
                );
                ekk.model_status_ = HighsModelStatus::kSolveError;
            }
        }

        // Debug here since not all simplex data will be correct until after
        // rebuild() when switching to Phase 2
        let no_debug = ekk.info_.num_dual_infeasibilities > 0
            && ekk.model_status_ == HighsModelStatus::kNotset;
        if !no_debug
            && self.debug_dual_simplex(ekk, "End of solvePhase1", false)
                == HighsDebugStatus::kLogicalError
        {
            self.solve_phase = kSolvePhaseError;
            return;
        }

        // Can also get here with solve_phase = kSolvePhaseError
        let solve_phase_ok = self.solve_phase == kSolvePhase1
            || self.solve_phase == kSolvePhase2
            || self.solve_phase == kSolvePhaseExit
            || self.solve_phase == kSolvePhaseError;
        if !solve_phase_ok {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "HEkkDual::solvePhase1 solve_phase == {} (solve call {}; iter {})\n",
                self.solve_phase,
                ekk.debug_solve_call_num_,
                ekk.iteration_count_
            );
        }
        debug_assert!(solve_phase_ok);
        if self.solve_phase == kSolvePhase2
            || self.solve_phase == kSolvePhaseExit
            || self.solve_phase == kSolvePhaseError
        {
            // Moving to phase 2 or exiting, so make sure that the simplex
            // bounds and nonbasic value/move correspond to the LP
            ekk.initialise_bound(SimplexAlgorithm::kDual, kSolvePhase2, false);
            ekk.initialise_nonbasic_value_and_move();
            if self.solve_phase == kSolvePhase2 {
                // Moving to phase 2 so possibly reinstate cost perturbation
                if ekk.dual_simplex_phase1_cleanup_level_
                    < ekk.options_.max_dual_simplex_phase1_cleanup_level
                {
                    // Allow cost perturbation for now, but may have to prevent it
                    // to avoid cleanup-perturbation loops
                    ekk.info_.allow_cost_shifting = true;
                    ekk.info_.allow_cost_perturbation = true;
                }
                // Comment if cost perturbation is not permitted
                if !ekk.info_.allow_cost_perturbation {
                    crate::highs_log_dev!(
                        &ekk.options_.log_options,
                        HighsLogType::kWarning,
                        "Moving to phase 2, but not allowing cost perturbation\n"
                    );
                }
            }
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::solvePhase2 (lines 817-1014, 198 lines)
    pub fn solve_phase2(&mut self, ekk: &mut HEkk) {
        // Have to set multi_chooseAgain = 1 since it may come in set to 0
        // from the end of Phase 1
        self.multi_choose_again = 1;
        // When starting a new phase the (updated) dual objective function
        // value isn't known.
        ekk.status_.has_primal_objective_value = false;
        ekk.status_.has_dual_objective_value = false;
        // Set rebuild_reason so that it's assigned when first tested
        self.rebuild_reason = kRebuildReasonNo;
        // Set solve_phase = kSolvePhase2 and ekk_instance_.solve_bailout_ = false so
        // they are set if solvePhase2() is called directly
        self.solve_phase = kSolvePhase2;
        ekk.solve_bailout_ = false;
        if ekk.bailout() {
            return;
        }
        // Report the phase start
        crate::highs_log_dev!(
            &ekk.options_.log_options,
            HighsLogType::kDetailed,
            "dual-phase-2-start\n"
        );
        // Collect free variables
        self.dual_row.create_freelist(ekk);

        // If there's no backtracking basis, save the initial basis in case
        // of backtracking
        if !ekk.info_.valid_backtracking_basis_ {
            ekk.put_backtracking_basis();
        }

        // Main solving structure
        ekk.analysis_.simplex_timer_start(IterateClock, 0);
        loop {
            // Outer loop of solvePhase2()
            // Rebuild all values, reinverting B if updates have been performed
            ekk.analysis_
                .simplex_timer_start(IterateDualRebuildClock, 0);
            self.rebuild(ekk);
            ekk.analysis_.simplex_timer_stop(IterateDualRebuildClock, 0);
            if self.solve_phase == kSolvePhaseError {
                ekk.model_status_ = HighsModelStatus::kSolveError;
                return;
            }
            if self.solve_phase == kSolvePhaseUnknown {
                // If backtracking, may change phase, so drop out
                ekk.analysis_.simplex_timer_stop(IterateClock, 0);
                return;
            }
            if ekk.bailout() {
                break;
            }
            if self.bailout_on_dual_objective(ekk) {
                break;
            }
            if self.dual_infeas_count > 0 {
                break;
            }
            loop {
                // Inner loop of solvePhase2()
                if self.debug_dual_simplex(ekk, "Before iteration", false)
                    == HighsDebugStatus::kLogicalError
                {
                    self.solve_phase = kSolvePhaseError;
                    return;
                }
                match ekk.info_.simplex_strategy {
                    kSimplexStrategyDualTasks => self.iterate_tasks(ekk),
                    kSimplexStrategyDualMulti => self.iterate_multi(ekk),
                    _ => self.iterate(ekk),
                }
                if ekk.bailout() {
                    break;
                }
                if self.bailout_on_dual_objective(ekk) {
                    break;
                }
                debug_assert!(self.solve_phase != kSolvePhaseTabooBasis);

                // If possibly dual unbounded, assess whether this implies
                // primal infeasibility.
                if self.rebuild_reason == kRebuildReasonPossiblyDualUnbounded {
                    self.assess_possibly_dual_unbounded(ekk);
                }

                if self.rebuild_reason != 0 {
                    break;
                }
            }
            if ekk.solve_bailout_ {
                break;
            }
            // If the data are fresh from rebuild(), possibly break out of the
            // outer loop to see what's occurred
            let finished =
                ekk.status_.has_fresh_rebuild && !ekk.rebuild_refactor(self.rebuild_reason);
            if finished && ekk.taboo_bad_basis_change() {
                // A bad basis change has had to be made taboo without any other
                // basis changes having been performed from a fresh rebuild.
                self.solve_phase = kSolvePhaseTabooBasis;
                return;
            }
            if finished {
                break;
            }
        }
        ekk.analysis_.simplex_timer_stop(IterateClock, 0);
        // Possibly return due to bailing out, having now stopped IterateClock
        if ekk.solve_bailout_ {
            return;
        }

        // If bailing out, should have done so already
        debug_assert!(!ekk.solve_bailout_);
        // Assess outcome of dual phase 2
        if self.dual_infeas_count > 0 {
            // There are dual infeasibilities so possibly switch to Phase 1 and
            // return.
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "dual-phase-2-found-free\n"
            );
            self.solve_phase = kSolvePhase1;
        } else if self.row_out == kNoRowChosen {
            // There is no candidate in CHUZR, even after rebuild so probably optimal
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "dual-phase-2-optimal\n"
            );
            // Remove any cost perturbations and see if basis is still dual feasible
            self.cleanup(ekk);
            if self.dual_infeas_count > 0 {
                // There are dual infeasibilities, so consider performing primal
                // simplex iterations to get dual feasibility
                self.solve_phase = kSolvePhaseOptimalCleanup;
            } else {
                // There are no dual infeasibilities so optimal!
                self.solve_phase = kSolvePhaseOptimal;
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kDetailed,
                    "problem-optimal\n"
                );
                ekk.model_status_ = HighsModelStatus::kOptimal;
            }
        } else if self.rebuild_reason == kRebuildReasonChooseColumnFail
            || self.rebuild_reason == kRebuildReasonExcessivePrimalValue
        {
            // chooseColumn has failed or excessive primal values have been created
            self.solve_phase = kSolvePhaseError;
            // Solve error is opaque to users, so put in some logging
            if self.rebuild_reason == kRebuildReasonChooseColumnFail {
                crate::highs_log_user!(
                    &ekk.options_.log_options,
                    HighsLogType::kError,
                    "Dual simplex ratio test failed due to excessive dual values: consider scaling down the LP objective coefficients\n"
                );
            } else {
                crate::highs_log_user!(
                    &ekk.options_.log_options,
                    HighsLogType::kError,
                    "Dual simplex detected excessive primal values: consider scaling down the LP bounds\n"
                );
            }
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "dual-phase-2-not-solved\n"
            );
            ekk.model_status_ = HighsModelStatus::kSolveError;
        } else {
            // Can only be that primal infeasibility has been detected
            debug_assert!(ekk.model_status_ == HighsModelStatus::kInfeasible);
            debug_assert!(self.solve_phase == kSolvePhaseExit);
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "problem-primal-infeasible\n"
            );
        }
        // Possibly debug unless before primal simplex clean-up (in which
        // case there will be dual infeasibilities).
        if self.solve_phase != kSolvePhaseOptimalCleanup
            && self.debug_dual_simplex(ekk, "End of solvePhase2", false)
                == HighsDebugStatus::kLogicalError
        {
            self.solve_phase = kSolvePhaseError;
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::rebuild (lines 1016-1128, 113 lines)
    pub fn rebuild(&mut self, ekk: &mut HEkk) {
        // Clear taboo flag from any bad basis changes
        ekk.clear_bad_basis_change_taboo_flag();

        // Decide whether refactorization should be performed
        let refactor_basis_matrix = ekk.rebuild_refactor(self.rebuild_reason);

        // Take a local copy of the rebuild reason and then reset the global value
        let local_rebuild_reason = self.rebuild_reason;
        self.rebuild_reason = kRebuildReasonNo;
        if refactor_basis_matrix {
            // Get a nonsingular inverse if possible.
            if !ekk.get_nonsingular_inverse(self.solve_phase) {
                self.solve_phase = kSolvePhaseError;
                return;
            }
            // Record the synthetic clock for INVERT, and zero it for UPDATE
            ekk.reset_synthetic_clock();
        }

        let alt_debug_level: HighsInt = -1;
        ekk.debug_nla_check_invert("HEkkDual::rebuild", alt_debug_level);

        if !ekk.status_.has_ar_matrix {
            // Don't have the row-wise matrix, so reinitialise it
            //
            // Should only happen when backtracking
            debug_assert!(ekk.info_.backtracking_);
            ekk.initialise_partitioned_rowwise_matrix();
            debug_assert!(
                ekk.ar_matrix_
                    .debug_partition_ok(&ekk.basis_.nonbasic_flag_)
            );
        }
        // Record whether the update objective value should be tested.
        let check_updated_objective_value = ekk.status_.has_dual_objective_value;
        let mut previous_dual_objective_value = -kHighsInf;
        if check_updated_objective_value {
            previous_dual_objective_value = ekk.info_.updated_dual_objective_value;
        }
        // Recompute dual solution
        ekk.compute_dual();

        if ekk.info_.backtracking_ {
            // If backtracking, may change phase, so drop out
            self.solve_phase = kSolvePhaseUnknown;
            return;
        }
        ekk.analysis_.simplex_timer_start(CorrectDualClock, 0);
        let mut free_infeasibility_count = self.dual_infeas_count;
        self.correct_dual_infeasibilities(ekk, &mut free_infeasibility_count);
        self.dual_infeas_count = free_infeasibility_count;
        ekk.analysis_.simplex_timer_stop(CorrectDualClock, 0);

        // Recompute primal solution
        ekk.compute_primal();

        // Collect primal infeasible as a list
        ekk.analysis_.simplex_timer_start(CollectPrIfsClock, 0);
        self.dual_rhs.create_array_of_primal_infeasibilities(ekk);
        let col_aq_density = ekk.info_.col_aq_density;
        self.dual_rhs.create_infeas_list(ekk, col_aq_density);
        ekk.analysis_.simplex_timer_stop(CollectPrIfsClock, 0);

        // Dual objective section
        ekk.compute_dual_objective_value(self.solve_phase);

        if check_updated_objective_value {
            // Apply the objective value correction due to computing duals
            // from scratch.
            let dual_objective_value_correction =
                ekk.info_.dual_objective_value - previous_dual_objective_value;
            ekk.info_.updated_dual_objective_value += dual_objective_value_correction;
        }
        // Now that there's a new dual_objective_value, reset the updated value
        ekk.info_.updated_dual_objective_value = ekk.info_.dual_objective_value;

        if !ekk.info_.run_quiet {
            ekk.compute_infeasibilities_for_reporting(SimplexAlgorithm::kDual, self.solve_phase);
            self.report_rebuild(ekk, local_rebuild_reason);
        }

        // Record the synthetic clock for INVERT, and zero it for UPDATE
        ekk.reset_synthetic_clock();

        // Dual simplex doesn't maintain the number of primal
        // infeasibilities, so set it to an illegal value now
        ekk.invalidate_primal_infeasibility_record();
        // Although dual simplex should always be dual feasible,
        // infeasibilities are only corrected in rebuild
        ekk.invalidate_dual_infeasibility_record();

        // Data are fresh from rebuild
        ekk.status_.has_fresh_rebuild = true;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::cleanup (lines 1130-1183, 54 lines)
    pub fn cleanup(&mut self, ekk: &mut HEkk) {
        if self.solve_phase == kSolvePhase1 {
            // Take action to prevent infinite loop.
            ekk.dual_simplex_phase1_cleanup_level_ += 1;
            let excessive_cleanup_calls = ekk.dual_simplex_phase1_cleanup_level_
                > ekk.options_.max_dual_simplex_phase1_cleanup_level;
            if excessive_cleanup_calls {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kError,
                    "Dual simplex cleanup level has exceeded limit of {}\n",
                    ekk.options_.max_dual_simplex_phase1_cleanup_level
                );
                debug_assert!(!excessive_cleanup_calls);
            }
        }
        crate::highs_log_dev!(
            &ekk.options_.log_options,
            HighsLogType::kDetailed,
            "dual-cleanup-shift\n"
        );
        // Remove perturbation and don't permit further perturbation
        ekk.initialise_cost(SimplexAlgorithm::kDual, kSolvePhaseUnknown, false);
        ekk.info_.allow_cost_perturbation = false;
        // No solve_phase term in initialiseBound is surely an omission -
        // when cleanup called in phase 1
        ekk.initialise_bound(SimplexAlgorithm::kDual, self.solve_phase, false);
        // Possibly take a copy of the original duals before recomputing them
        let mut _original_work_dual: Vec<f64> = Vec::new();
        if ekk.options_.highs_debug_level > kHighsDebugLevelCheap {
            _original_work_dual = ekk.info_.work_dual_.clone();
        }
        // Compute the dual values
        ekk.compute_dual();
        // Compute the dual infeasibilities
        ekk.compute_simplex_dual_infeasible();
        self.dual_infeas_count = ekk.info_.num_dual_infeasibilities;

        // Compute the dual objective value
        ekk.compute_dual_objective_value(self.solve_phase);
        // Now that there's a new dual_objective_value, reset the updated value
        ekk.info_.updated_dual_objective_value = ekk.info_.dual_objective_value;

        if !ekk.info_.run_quiet {
            // Report the primal infeasibilities
            ekk.compute_simplex_primal_infeasible();
            // In phase 1, report the simplex LP dual infeasibilities
            // In phase 2, report the simplex dual infeasibilities (known)
            if self.solve_phase == kSolvePhase1 {
                ekk.compute_simplex_lp_dual_infeasible();
            }
            self.report_rebuild(ekk, kRebuildReasonCleanup);
        }
    }
}

#[cfg(test)]
mod h16_dual_free_list_tests {
    use super::*;
    use crate::prelude::*;

    /// HiGHS 1.16 (#3353): check/instances/issue-2300.mps (diff lines 29-258), built by
    /// hand from its ROWS / COLUMNS / RHS / BOUNDS sections (50 rows, 35 columns, MAX).
    fn issue_2300_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 35,
            num_row_: 50,
            sense_: ObjSense::kMaximize,
            ..Default::default()
        };
        lp.col_cost_ = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0,
        ];
        lp.col_lower_ = vec![
            0.0, -kHighsInf, 0.0, -kHighsInf, 0.0, -kHighsInf, 0.0, -kHighsInf, 0.0, 0.0,
            -kHighsInf, -kHighsInf, -kHighsInf, 0.0, -kHighsInf, -kHighsInf, 0.0, 0.0, 0.0, 0.0,
            0.0, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf, 0.0, -kHighsInf, 0.0, -kHighsInf,
            -kHighsInf, 0.0, -kHighsInf, -kHighsInf, 0.0, 0.0,
        ];
        lp.col_upper_ = vec![
            kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
            kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
            kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
            kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
            kHighsInf, kHighsInf, kHighsInf,
        ];
        lp.row_lower_ = vec![
            -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf, 0.0, 0.0, -kHighsInf, -kHighsInf,
            -kHighsInf, 0.0, -kHighsInf, 0.0, 0.0, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf,
            -kHighsInf, -kHighsInf, 0.0, -kHighsInf, -kHighsInf, -kHighsInf, 0.0, -kHighsInf,
            -kHighsInf, 0.0, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf,
            -kHighsInf, -kHighsInf, -kHighsInf, 0.0, -kHighsInf, 0.0, 0.0, 0.0, -kHighsInf,
            -kHighsInf, -kHighsInf, -kHighsInf, 0.0, -kHighsInf, -kHighsInf, -kHighsInf,
            -kHighsInf, 0.0,
        ];
        lp.row_upper_ = vec![
            0.0, 0.0, 0.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        lp.a_matrix_.num_col_ = 35;
        lp.a_matrix_.num_row_ = 50;
        lp.a_matrix_.start_ = vec![
            0, 3, 7, 10, 17, 20, 24, 27, 34, 37, 40, 45, 49, 55, 58, 62, 67, 70, 73, 76, 79, 82,
            87, 104, 109, 113, 116, 120, 123, 128, 129, 132, 137, 139, 142, 145,
        ];
        lp.a_matrix_.index_ = vec![
            3, 10, 22, 10, 11, 22, 33, 3, 14, 48, 8, 12, 19, 25, 33, 35, 44, 3, 15, 20, 4, 26, 29,
            34, 3, 13, 47, 1, 5, 16, 23, 33, 35, 44, 3, 18, 24, 3, 7, 30, 5, 6, 33, 36, 49, 0, 33,
            40, 44, 23, 32, 33, 39, 44, 45, 3, 8, 25, 9, 17, 33, 41, 5, 18, 24, 26, 33, 3, 29, 34,
            3, 32, 45, 1, 3, 16, 3, 28, 46, 3, 31, 42, 12, 21, 27, 33, 37, 4, 5, 9, 11, 12, 13, 19,
            23, 26, 33, 35, 37, 38, 39, 44, 47, 49, 2, 26, 33, 37, 43, 28, 33, 39, 46, 3, 6, 36,
            14, 33, 48, 49, 2, 3, 43, 4, 15, 20, 33, 49, 33, 0, 3, 40, 7, 9, 19, 30, 33, 31, 42, 3,
            17, 41, 3, 21, 27,
        ];
        lp.a_matrix_.value_ = vec![
            1.0,
            -1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            1.0,
            2.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -2.0,
            -1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            2.0,
            1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            -2.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            3.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            -1.0,
            1.0,
            1.0,
            1.0,
            -2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            -2.928694748e-05,
            -1.999853566,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            5.0,
            1.0,
            1.0,
            -1.0,
            1.0,
            -1.0,
            1.0,
            1.0,
            -1.0,
            -1.0,
            2.0,
            1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
            1.0,
            -1.0,
            -1.0,
        ];
        lp
    }

    /// Solves the LP with the option set to `flag`, default options (presolve on).
    /// Returns (run status, model status, objective).
    fn solve_issue_2300(flag: bool) -> (HighsStatus, HighsModelStatus, f64) {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("log_to_console", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_dual_free_list", flag),
            HighsStatus::kOk
        );
        // the 1.15.1 oracle solves with simplex; lp_solver_auto would hand this LP to HiPO
        assert_eq!(
            highs.set_option_value_bool("lp_solver_auto", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(issue_2300_lp()), HighsStatus::kOk);
        let run_status = highs.run();
        (
            run_status,
            highs.get_model_status(),
            highs.get_info().objective_function_value,
        )
    }

    /// Flag on: the phase-1 clear reaches the new code and the LP is Optimal with the
    /// PR's objective (C++ check: 4.0735475989e-04 within 1e-12).
    #[test]
    fn h16_dual_free_list_on_optimal() {
        let (run_status, status, objective) = solve_issue_2300(true);
        assert_eq!(run_status, HighsStatus::kOk);
        assert_eq!(status, HighsModelStatus::kOptimal);
        assert!((objective - 4.0735475989e-04).abs() < 1e-12, "{objective}");
    }

    /// Flag off: the 1.15.1 result is unchanged. The 1.15.1 oracle (highspy 1.15.1,
    /// presolve on) gives kSolveError with objective 0 on this LP.
    #[test]
    fn h16_dual_free_list_off_unchanged() {
        let (run_status, status, objective) = solve_issue_2300(false);
        assert_eq!(run_status, HighsStatus::kError);
        assert_eq!(status, HighsModelStatus::kSolveError);
        assert_eq!(objective, 0.0);
    }
}
