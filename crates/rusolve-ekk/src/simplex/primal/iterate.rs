//! HiGHS source: simplex/HEkkPrimal.cpp (lines 686-1136)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Rebuild, iterate, CHUZC (`chuzc`/`chooseColumn`), entering-variable checks.
//! See `HEkkPrimal` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkPrimal {
    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::rebuild (lines 686-806, 121 lines)
    pub fn rebuild(&mut self, ekk: &mut HEkk) {
        // Clear taboo flag from any bad basis changes
        ekk.clear_bad_basis_change_taboo_flag();

        // Record whether the update objective value should be tested.
        let check_updated_objective_value = ekk.status_.has_primal_objective_value;
        let mut previous_primal_objective_value = -kHighsInf;
        if check_updated_objective_value {
            previous_primal_objective_value = ekk.info_.updated_primal_objective_value;
        }

        // Decide whether refactorization should be performed
        let refactor_basis_matrix = ekk.rebuild_refactor(self.rebuild_reason);

        // Take a local copy of the rebuild reason and then reset the global value
        let local_rebuild_reason = self.rebuild_reason;
        self.rebuild_reason = kRebuildReasonNo;
        if refactor_basis_matrix {
            if !ekk.get_nonsingular_inverse(self.solve_phase) {
                self.solve_phase = kSolvePhaseError;
                return;
            }
            // Record the synthetic clock for INVERT, and zero it for UPDATE
            ekk.reset_synthetic_clock();
        }
        if !ekk.status_.has_ar_matrix {
            // Should only happen when backtracking
            debug_assert!(ekk.info_.backtracking_);
            ekk.initialise_partitioned_rowwise_matrix();
            debug_assert!(
                ekk.ar_matrix_
                    .debug_partition_ok(&ekk.basis_.nonbasic_flag_)
            );
        }

        if ekk.info_.backtracking_ {
            // If backtracking, may change phase, so drop out
            self.solve_phase = kSolvePhaseUnknown;
            return;
        }

        ekk.compute_primal();
        if self.solve_phase == kSolvePhase2 {
            let correct_primal_ok = self.correct_primal(ekk, false);
            if kAllowDeveloperAssert {
                debug_assert!(correct_primal_ok);
            }
        }
        self.get_basic_primal_infeasibility(ekk);
        if ekk.info_.num_primal_infeasibilities > 0 {
            // Primal infeasibilities so should be in phase 1
            if self.solve_phase == kSolvePhase2 {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "HEkkPrimal::rebuild switching back to phase 1 from phase 2\n"
                );
                self.solve_phase = kSolvePhase1;
            }
            self.phase1_compute_dual(ekk);
        } else {
            // No primal infeasibilities so in phase 2. Reset costs if was
            // previously in phase 1
            if self.solve_phase == kSolvePhase1 {
                ekk.initialise_cost(SimplexAlgorithm::kPrimal, self.solve_phase, false);
                self.solve_phase = kSolvePhase2;
            }
            ekk.compute_dual();
        }
        ekk.compute_simplex_dual_infeasible();
        ekk.compute_primal_objective_value();
        if check_updated_objective_value {
            // Apply the objective value correction due to computing primal
            // values from scratch.
            let primal_objective_value_correction =
                ekk.info_.primal_objective_value - previous_primal_objective_value;
            ekk.info_.updated_primal_objective_value += primal_objective_value_correction;
        }
        // Now that there's a new dual_objective_value, reset the updated value
        ekk.info_.updated_primal_objective_value = ekk.info_.primal_objective_value;

        self.report_rebuild(ekk, local_rebuild_reason);

        // Record the synthetic clock for INVERT, and zero it for UPDATE
        ekk.reset_synthetic_clock();

        // Determine whether to use hyper-sparse CHUZC
        if self.solve_phase == kSolvePhase1 {
            self.use_hyper_chuzc = false;
        } else {
            self.use_hyper_chuzc = false; // true;
        }
        self.hyper_choose_column_clear(ekk);

        self.num_flip_since_rebuild = 0;
        // Data are fresh from rebuild
        ekk.status_.has_fresh_rebuild = true;
        debug_assert!(self.solve_phase == kSolvePhase1 || self.solve_phase == kSolvePhase2);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::iterate (lines 808-931, 124 lines)
    pub fn iterate(&mut self, ekk: &mut HEkk) {
        let from_check_iter: HighsInt = 15;
        let to_check_iter: HighsInt = from_check_iter + 10;
        if ekk.debug_solve_report_ {
            ekk.debug_iteration_report_ =
                ekk.iteration_count_ >= from_check_iter && ekk.iteration_count_ <= to_check_iter;
            if ekk.debug_iteration_report_ {
                println!(
                    "HEkkPrimal::iterate Debug iteration {}",
                    ekk.iteration_count_
                );
            }
        }

        if self.debug_primal_simplex(ekk, "Before iteration", false)
            == HighsDebugStatus::kLogicalError
        {
            self.solve_phase = kSolvePhaseError;
            return;
        }
        // Initialise row_out so that aborting iteration before CHUZR due to
        // numerical test of chosen reduced cost can be spotted
        self.row_out = kNoRowSought;
        // Perform CHUZC
        self.chuzc(ekk);
        if self.variable_in == -1 {
            self.rebuild_reason = kRebuildReasonPossiblyOptimal;
            return;
        }

        // Perform FTRAN - and dual value cross-check to decide whether to use the
        // variable
        if !self.use_variable_in(ekk) {
            if self.rebuild_reason != 0 {
                debug_assert!(self.rebuild_reason == kRebuildReasonPossiblySingularBasis);
            }
            return;
        }
        debug_assert!(self.rebuild_reason == 0);

        // Perform CHUZR
        if self.solve_phase == kSolvePhase1 {
            self.phase1_choose_row(ekk);
            debug_assert!(self.row_out != kNoRowSought);
            if self.row_out == kNoRowChosen {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kError,
                    "Primal phase 1 choose row failed\n"
                );
                self.solve_phase = kSolvePhaseError;
                return;
            }
        } else {
            self.choose_row(ekk);
        }
        debug_assert!(self.rebuild_reason == 0);

        // Consider whether to perform a bound swap
        debug_assert!(self.solve_phase == kSolvePhase2 || self.row_out >= 0);
        self.consider_bound_swap(ekk);
        if self.rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded {
            return;
        }
        debug_assert!(self.rebuild_reason == 0);

        if self.row_out >= 0 {
            // Perform unit BTRAN and PRICE to get pivotal row - and do a
            // numerical check.
            self.assess_pivot(ekk);
            if self.rebuild_reason != 0 {
                debug_assert!(self.rebuild_reason == kRebuildReasonPossiblySingularBasis);
                return;
            }
        }

        if self.is_bad_basis_change(ekk) {
            return;
        }

        // Any pivoting is numerically acceptable, so perform update.
        self.update(ekk);
        // Force rebuild if there are no infeasibilities in phase 1
        if ekk.info_.num_primal_infeasibilities == 0 && self.solve_phase == kSolvePhase1 {
            self.rebuild_reason = kRebuildReasonPossiblyPhase1Feasible;
        }

        let ok_rebuild_reason = self.rebuild_reason == kRebuildReasonNo
            || self.rebuild_reason == kRebuildReasonPossiblyPhase1Feasible
            || self.rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex
            || self.rebuild_reason == kRebuildReasonSyntheticClockSaysInvert
            || self.rebuild_reason == kRebuildReasonUpdateLimitReached;
        if !ok_rebuild_reason {
            println!(
                "HEkkPrimal::rebuild Solve {}; Iter {}: rebuild_reason = {}",
                ekk.debug_solve_call_num_, ekk.iteration_count_, self.rebuild_reason
            );
        }
        debug_assert!(ok_rebuild_reason);
        debug_assert!(self.solve_phase == kSolvePhase1 || self.solve_phase == kSolvePhase2);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::chuzc (lines 933-975, 43 lines)
    pub fn chuzc(&mut self, ekk: &mut HEkk) {
        if self.done_next_chuzc {
            debug_assert!(self.use_hyper_chuzc);
        }
        // workDual aliases ekk.info_.work_dual_; moved out only for the taboo calls
        let mut work_dual = std::mem::take(&mut ekk.info_.work_dual_);
        ekk.apply_taboo_variable_in(&mut work_dual, 0.0);
        ekk.info_.work_dual_ = work_dual;
        if self.use_hyper_chuzc {
            // Perform hyper-sparse CHUZC and then check result using full CHUZC
            if !self.done_next_chuzc {
                self.choose_column(ekk, true);
            }
            let check_hyper_chuzc = true;
            if check_hyper_chuzc {
                let hyper_sparse_variable_in = self.variable_in;
                self.choose_column(ekk, false);
                let mut hyper_sparse_measure: f64 = 0.0;
                if hyper_sparse_variable_in >= 0 {
                    let wd = ekk.info_.work_dual_[hyper_sparse_variable_in as usize];
                    let squared_dual_infeasibility = wd * wd;
                    hyper_sparse_measure = squared_dual_infeasibility
                        / self.edge_weight_[hyper_sparse_variable_in as usize];
                }
                let mut measure: f64 = 0.0;
                if self.variable_in >= 0 {
                    let wd = ekk.info_.work_dual_[self.variable_in as usize];
                    let squared_dual_infeasibility = wd * wd;
                    measure =
                        squared_dual_infeasibility / self.edge_weight_[self.variable_in as usize];
                }
                let abs_measure_error = (hyper_sparse_measure - measure).abs();
                let measure_error = abs_measure_error > 1e-12;

                // todo this fails on some rarer occasions, e.g. on glass4
                debug_assert!(!measure_error);
                self.variable_in = hyper_sparse_variable_in;
            }
        } else {
            self.choose_column(ekk, false);
        }
        let mut work_dual = std::mem::take(&mut ekk.info_.work_dual_);
        ekk.unapply_taboo_variable_in(&mut work_dual);
        ekk.info_.work_dual_ = work_dual;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::chooseColumn (lines 977-1072, 96 lines)
    /// C++ default arguments: hyper_sparse = false
    pub fn choose_column(&mut self, ekk: &mut HEkk, hyper_sparse: bool) {
        debug_assert!(!hyper_sparse || !self.done_next_chuzc);
        let mut best_measure: f64 = 0.0;
        self.variable_in = -1;

        let local_use_hyper_chuzc = hyper_sparse;
        // Consider nonbasic free columns first
        let num_nonbasic_free_col: HighsInt = *self.nonbasic_free_col_set.count();
        if local_use_hyper_chuzc {
            if !self.initialise_hyper_chuzc {
                self.hyper_choose_column(ekk);
            }
            if self.initialise_hyper_chuzc {
                ekk.analysis_
                    .simplex_timer_start(ChuzcHyperInitialiselClock, 0);
                self.num_hyper_chuzc_candidates = 0;
                if num_nonbasic_free_col != 0 {
                    for ix in 0..num_nonbasic_free_col {
                        let i_col = self.nonbasic_free_col_set.entry()[ix as usize];
                        let dual_infeasibility = ekk.info_.work_dual_[i_col as usize].abs();
                        if dual_infeasibility > self.dual_feasibility_tolerance {
                            let measure = dual_infeasibility * dual_infeasibility
                                / self.edge_weight_[i_col as usize];
                            add_to_decreasing_heap(
                                &mut self.num_hyper_chuzc_candidates,
                                self.max_num_hyper_chuzc_candidates,
                                &mut self.hyper_chuzc_measure,
                                &mut self.hyper_chuzc_candidate,
                                measure,
                                i_col,
                            );
                        }
                    }
                }
                // Now look at other columns
                for i_col in 0..self.num_tot {
                    let dual_infeasibility = (-(ekk.basis_.nonbasic_move_[i_col as usize]
                        as HighsInt)) as f64
                        * ekk.info_.work_dual_[i_col as usize];
                    if dual_infeasibility > self.dual_feasibility_tolerance {
                        let measure = dual_infeasibility * dual_infeasibility
                            / self.edge_weight_[i_col as usize];
                        add_to_decreasing_heap(
                            &mut self.num_hyper_chuzc_candidates,
                            self.max_num_hyper_chuzc_candidates,
                            &mut self.hyper_chuzc_measure,
                            &mut self.hyper_chuzc_candidate,
                            measure,
                            i_col,
                        );
                    }
                }
                // Sort the heap
                sort_decreasing_heap(
                    self.num_hyper_chuzc_candidates,
                    &mut self.hyper_chuzc_measure,
                    &mut self.hyper_chuzc_candidate,
                );
                self.initialise_hyper_chuzc = false;
                ekk.analysis_
                    .simplex_timer_stop(ChuzcHyperInitialiselClock, 0);
                // Choose the first entry - if there is one
                if self.num_hyper_chuzc_candidates != 0 {
                    self.variable_in = self.hyper_chuzc_candidate[1];
                    best_measure = self.hyper_chuzc_measure[1];
                    self.max_hyper_chuzc_non_candidate_measure =
                        self.hyper_chuzc_measure[self.num_hyper_chuzc_candidates as usize];
                    if self.report_hyper_chuzc {
                        println!(
                            "Full CHUZC: Max         measure is {:9.4e} for column {:4}, and max non-candidate measure of  {:9.4e}",
                            best_measure,
                            self.variable_in,
                            self.max_hyper_chuzc_non_candidate_measure
                        );
                    }
                }
            }
        } else {
            ekk.analysis_.simplex_timer_start(ChuzcPrimalClock, 0);
            // Choose any attractive nonbasic free column
            if num_nonbasic_free_col != 0 {
                for ix in 0..num_nonbasic_free_col {
                    let i_col = self.nonbasic_free_col_set.entry()[ix as usize];
                    let dual_infeasibility = ekk.info_.work_dual_[i_col as usize].abs();
                    if dual_infeasibility > self.dual_feasibility_tolerance
                        && dual_infeasibility * dual_infeasibility
                            > best_measure * self.edge_weight_[i_col as usize]
                    {
                        self.variable_in = i_col;
                        best_measure = dual_infeasibility * dual_infeasibility
                            / self.edge_weight_[i_col as usize];
                    }
                }
            }
            // Now look at other columns
            for i_col in 0..self.num_tot {
                let dual_infeasibility = (-(ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt))
                    as f64
                    * ekk.info_.work_dual_[i_col as usize];
                if dual_infeasibility > self.dual_feasibility_tolerance
                    && dual_infeasibility * dual_infeasibility
                        > best_measure * self.edge_weight_[i_col as usize]
                {
                    self.variable_in = i_col;
                    best_measure =
                        dual_infeasibility * dual_infeasibility / self.edge_weight_[i_col as usize];
                }
            }
            ekk.analysis_.simplex_timer_stop(ChuzcPrimalClock, 0);
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::useVariableIn (lines 1074-1136, 63 lines)
    pub fn use_variable_in(&mut self, ekk: &mut HEkk) -> bool {
        // rebuild_reason = kRebuildReasonPossiblySingularBasis is set if
        // numerical trouble is detected
        let updated_theta_dual = ekk.info_.work_dual_[self.variable_in as usize];
        // Determine the move direction - can't use nonbasicMove_[variable_in]
        // due to free columns
        self.move_in = if updated_theta_dual > 0.0 { -1 } else { 1 };
        // Unless the variable is free, nonbasicMove[variable_in] should be the same
        // as move_in
        if ekk.basis_.nonbasic_move_[self.variable_in as usize] != 0 {
            debug_assert!(
                ekk.basis_.nonbasic_move_[self.variable_in as usize] as HighsInt == self.move_in
            );
        }
        //
        // FTRAN
        //
        // Compute pivot column
        ekk.pivot_column_ftran(self.variable_in, &mut self.col_aq);
        // Compute the dual for the pivot column and compare it with the
        // updated value
        let computed_theta_dual =
            ekk.compute_dual_for_tableau_column(self.variable_in, &self.col_aq);
        ekk.debug_updated_dual(updated_theta_dual, computed_theta_dual);
        // Feed in the computed dual value.
        ekk.info_.work_dual_[self.variable_in as usize] = computed_theta_dual;
        // Reassign theta_dual to be the computed value
        self.theta_dual = ekk.info_.work_dual_[self.variable_in as usize];
        // Determine whether theta_dual is too small or has changed sign
        let theta_dual_small = self.theta_dual.abs() <= self.dual_feasibility_tolerance;
        let theta_dual_sign_error = updated_theta_dual * computed_theta_dual <= 0.0;

        // If theta_dual is small, then it's no longer a dual infeasibility,
        // so reduce the number of dual infeasibilities.
        if theta_dual_small {
            ekk.info_.num_dual_infeasibilities -= 1;
        }
        if theta_dual_small || theta_dual_sign_error {
            // The computed dual is small or has a sign error, so don't use it
            let mut theta_dual_size = "";
            if theta_dual_small {
                theta_dual_size = "; too small";
            }
            let mut theta_dual_sign = "";
            if theta_dual_sign_error {
                theta_dual_sign = "; sign error";
            }
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "Chosen entering variable {} (Iter = {}; Update = {}) has computed (updated) dual of {:10.4e} ({:10.4e}) so don't use it{}{}\n",
                self.variable_in,
                ekk.iteration_count_,
                ekk.info_.update_count,
                computed_theta_dual,
                updated_theta_dual,
                theta_dual_size,
                theta_dual_sign
            );
            // If a significant computed dual has sign error, consider reinverting
            if !theta_dual_small && ekk.info_.update_count > 0 {
                self.rebuild_reason = kRebuildReasonPossiblySingularBasis;
            }
            self.hyper_choose_column_clear(ekk);
            return false;
        }
        true
    }
}
