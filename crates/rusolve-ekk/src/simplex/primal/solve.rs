//! HiGHS source: simplex/HEkkPrimal.cpp (lines 20-684)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Primal simplex driver, phase 1 / phase 2 loops, cleanup.
//! See `HEkkPrimal` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkPrimal {
    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::solve (lines 20-327, 308 lines)
    /// C++ default arguments: pass_force_phase2 = false
    pub fn solve(&mut self, ekk: &mut HEkk, pass_force_phase2: bool) -> HighsStatus {
        // Initialise control data for a particular solve
        self.initialise_solve(ekk);
        // Assumes that the LP has a positive number of rows
        if ekk.is_unconstrained_lp() {
            return ekk.return_from_solve(HighsStatus::kError);
        }

        if !ekk.status_.has_invert {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kError,
                "HEkkPrimal::solve called without INVERT\n"
            );
            debug_assert!(ekk.status_.has_fresh_invert);
            return ekk.return_from_solve(HighsStatus::kError);
        }

        if self.debug_primal_simplex(ekk, "Initialise", true) == HighsDebugStatus::kLogicalError {
            return ekk.return_from_solve(HighsStatus::kError);
        }

        // Get the nonbasic free column set
        self.get_nonbasic_free_column_set(ekk);

        let primal_feasible_with_unperturbed_bounds = ekk.info_.num_primal_infeasibilities == 0;
        let force_phase2 = pass_force_phase2
            || ekk.info_.max_primal_infeasibility * ekk.info_.max_primal_infeasibility
                < ekk.options_.primal_feasibility_tolerance;
        // Determine whether the solution is near-optimal. Values 1000 and
        // 1e-3 (ensuring sum<1) are unimportant, as the sum of dual
        // infeasibilities for near-optimal solutions is typically many
        // orders of magnitude smaller than 1, and the sum of dual
        // infeasibilities will be very much larger for non-trivial LPs
        // that are primal feasible for a logical or crash basis.
        //
        // Consider there to be no primal infeasibilities if there are none,
        // or if phase 2 is forced, in which case any primal infeasibilities
        // will be shifted
        let no_simplex_primal_infeasibilities =
            primal_feasible_with_unperturbed_bounds || force_phase2;
        let near_optimal = ekk.info_.num_dual_infeasibilities < 1000
            && ekk.info_.max_dual_infeasibility < 1e-3
            && no_simplex_primal_infeasibilities;
        // For reporting, save primal infeasibility data for the LP without
        // bound perturbations
        let unperturbed_num_infeasibilities = ekk.info_.num_primal_infeasibilities;
        let unperturbed_max_infeasibility = ekk.info_.max_primal_infeasibility;
        let unperturbed_sum_infeasibilities = ekk.info_.sum_primal_infeasibilities;
        if near_optimal {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "Primal feasible and num / max / sum dual infeasibilities of {} / {} / {}, so near-optimal\n",
                ekk.info_.num_dual_infeasibilities,
                ekk.info_.max_dual_infeasibility,
                ekk.info_.sum_dual_infeasibilities
            );
        }

        // Perturb bounds according to whether the solution is near-optimal
        let perturb_bounds = !near_optimal;
        if !perturb_bounds {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "Near-optimal, so don't use bound perturbation\n"
            );
        }
        if perturb_bounds && ekk.info_.primal_simplex_bound_perturbation_multiplier != 0.0 {
            ekk.initialise_bound(
                SimplexAlgorithm::kPrimal,
                kSolvePhaseUnknown,
                perturb_bounds,
            );
            ekk.initialise_nonbasic_value_and_move();
            ekk.compute_primal();
            ekk.compute_simplex_primal_infeasible();
        }

        // Check whether the time/iteration limit has been reached. First
        // point at which a non-error return can occur
        if ekk.bailout() {
            return ekk.return_from_solve(HighsStatus::kWarning);
        }

        // Now to do some iterations!
        let mut num_primal_infeasibility = ekk.info_.num_primal_infeasibilities;
        self.solve_phase = if num_primal_infeasibility > 0 {
            kSolvePhase1
        } else {
            kSolvePhase2
        };
        if force_phase2 {
            // Dual infeasibilities without cost perturbation involved
            // fixed variables or were (at most) small, so can easily be
            // removed by flips for and fixed variables shifts for the rest
            self.solve_phase = kSolvePhase2;
            if !pass_force_phase2 {
                let local_report = false; // true;
                if !primal_feasible_with_unperturbed_bounds && local_report {
                    print!(
                        "Solve {}: Forcing phase 2 since near primal feasible with unperturbed costs\n\
                         num / max / sum primal infeasibilities\n\
                         {} / {:11.4} / {:11.4} (  perturbed bounds)\n\
                         {} / {:11.4} / {:11.4} (unperturbed bounds)\n",
                        ekk.debug_solve_call_num_,
                        ekk.info_.num_primal_infeasibilities,
                        ekk.info_.max_primal_infeasibility,
                        ekk.info_.sum_primal_infeasibilities,
                        unperturbed_num_infeasibilities,
                        unperturbed_max_infeasibility,
                        unperturbed_sum_infeasibilities
                    );
                }
            }
        }
        if ekk.debug_ok_for_solve(SimplexAlgorithm::kPrimal, self.solve_phase)
            == HighsDebugStatus::kLogicalError
        {
            return ekk.return_from_solve(HighsStatus::kError);
        }
        // Resize the copy of scattered edge weights for backtracking
        ekk.info_
            .backtracking_basis_edge_weight_
            .resize(self.num_tot as usize, 0.0);

        // The major solving loop
        //
        // Possibly write out the column header for iteration reports, and
        // initialise records for primal correction reporting
        self.local_report_iter(ekk, true);
        self.correct_primal(ekk, true);
        while self.solve_phase != 0 {
            let it0 = ekk.iteration_count_;
            // When starting a new phase the (updated) primal objective function
            // value isn't known. Indicate this so that when the value
            // computed from scratch in rebuild() isn't checked against the
            // updated value
            ekk.status_.has_primal_objective_value = false;
            if self.solve_phase == kSolvePhaseUnknown {
                // Determine the number of primal infeasibilities, and hence the solve
                // phase
                ekk.compute_simplex_primal_infeasible();
                num_primal_infeasibility = ekk.info_.num_primal_infeasibilities;
                self.solve_phase = if num_primal_infeasibility > 0 {
                    kSolvePhase1
                } else {
                    kSolvePhase2
                };
                if ekk.info_.backtracking_ {
                    // Backtracking
                    ekk.initialise_cost(SimplexAlgorithm::kPrimal, self.solve_phase, false);
                    ekk.initialise_nonbasic_value_and_move();
                    // Can now forget that we might have been backtracking
                    ekk.info_.backtracking_ = false;
                }
            }
            debug_assert!(self.solve_phase == kSolvePhase1 || self.solve_phase == kSolvePhase2);
            if self.solve_phase == kSolvePhase1 {
                //
                // Phase 1
                //
                // solve_phase = kSolvePhase1 if the iteration or time limit has
                // been reached
                //
                // solve_phase = kSolvePhase2 if there are no primal infeasibilities
                //
                // solve_phase = kSolvePhaseUnknown if backtracking
                //
                // solve_phase = kSolvePhaseExit if primal infeasibility is
                // detected, in which case model_status_ =
                // HighsModelStatus::kInfeasible is set
                //
                // solve_phase = kSolvePhaseTabooBasis is set if only basis change is
                // taboo
                //
                // solve_phase = kSolvePhaseError is set if an error occurs
                self.solve_phase1(ekk);
                debug_assert!(
                    self.solve_phase == kSolvePhase1
                        || self.solve_phase == kSolvePhase2
                        || self.solve_phase == kSolvePhaseUnknown
                        || self.solve_phase == kSolvePhaseExit
                        || self.solve_phase == kSolvePhaseTabooBasis
                        || self.solve_phase == kSolvePhaseError
                );
                ekk.info_.primal_phase1_iteration_count += ekk.iteration_count_ - it0;
            } else if self.solve_phase == kSolvePhase2 {
                //
                // Phase 2
                //
                // solve_phase = kSolvePhaseOptimal if there are no dual
                // infeasibilities
                //
                // solve_phase = kSolvePhase1 if there are primal
                // infeasibilities
                //
                // solve_phase = kSolvePhase2 if the iteration or time limit has
                // been reached
                //
                // solve_phase = kSolvePhaseOptimalCleanup if, after removing bound
                // shifts, there are primal infeasibilities to clean up
                //
                // solve_phase = kSolvePhaseUnknown if backtracking
                //
                // solve_phase = kSolvePhaseExit if primal unboundedness is
                // detected, in which case model_status_ =
                // HighsModelStatus::kUnbounded is set
                //
                // solve_phase = kSolvePhaseTabooBasis is set if only basis change is
                // taboo
                //
                // solve_phase = kSolvePhaseError is set if an error occurs
                self.solve_phase2(ekk);
                debug_assert!(
                    self.solve_phase == kSolvePhaseOptimal
                        || self.solve_phase == kSolvePhase1
                        || self.solve_phase == kSolvePhase2
                        || self.solve_phase == kSolvePhaseOptimalCleanup
                        || self.solve_phase == kSolvePhaseUnknown
                        || self.solve_phase == kSolvePhaseExit
                        || self.solve_phase == kSolvePhaseTabooBasis
                        || self.solve_phase == kSolvePhaseError
                );
                debug_assert!(
                    self.solve_phase != kSolvePhaseExit
                        || ekk.model_status_ == HighsModelStatus::kUnbounded
                );
                ekk.info_.primal_phase2_iteration_count += ekk.iteration_count_ - it0;
            } else {
                // Should only be kSolvePhase1 or kSolvePhase2
                ekk.model_status_ = HighsModelStatus::kSolveError;
                return ekk.return_from_solve(HighsStatus::kError);
            }
            // Return if bailing out from solve
            if ekk.solve_bailout_ {
                return ekk.return_from_solve(HighsStatus::kWarning);
            }
            // Can have all possible cases of solve_phase
            debug_assert!(self.solve_phase >= kSolvePhaseMin && self.solve_phase <= kSolvePhaseMax);
            // Look for scenarios when the major solving loop ends
            if self.solve_phase == kSolvePhaseTabooBasis {
                // Only basis change is taboo so return HighsStatus::kWarning
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kInfo,
                    "HEkkPrimal::solve Only basis change is taboo\n"
                );
                ekk.model_status_ = HighsModelStatus::kUnknown;
                return ekk.return_from_solve(HighsStatus::kWarning);
            }
            if self.solve_phase == kSolvePhaseError {
                // Solver error so return HighsStatus::kError
                ekk.model_status_ = HighsModelStatus::kSolveError;
                return ekk.return_from_solve(HighsStatus::kError);
            }
            if self.solve_phase == kSolvePhaseExit {
                // LP identified as not having an optimal solution
                debug_assert!(
                    ekk.model_status_ == HighsModelStatus::kInfeasible
                        || ekk.model_status_ == HighsModelStatus::kUnbounded
                );
                // If infeasible, save the primal phase 1 dual values before
                // they are overwritten with the duals for the original
                // objective
                if ekk.model_status_ == HighsModelStatus::kInfeasible {
                    ekk.primal_phase1_dual_ = ekk.info_.work_dual_.clone();
                }
                break;
            }
            if self.solve_phase == kSolvePhaseOptimalCleanup {
                // Primal infeasibilities after phase 2. Dual feasible with
                // primal infeasibilities so use dual simplex to clean up
                break;
            }
            // If solve_phase == kSolvePhaseOptimal == 0 then major solving
            // loop ends naturally since solve_phase is false
        }
        // If bailing out, should have returned already
        debug_assert!(!ekk.solve_bailout_);
        // Should only have these cases
        debug_assert!(
            self.solve_phase == kSolvePhaseExit
                || self.solve_phase == kSolvePhaseOptimal
                || self.solve_phase == kSolvePhaseOptimalCleanup
        );
        if self.solve_phase == kSolvePhaseOptimal {
            ekk.model_status_ = HighsModelStatus::kOptimal;
        }

        if self.solve_phase == kSolvePhaseOptimalCleanup {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "HEkkPrimal:: Using dual simplex to try to clean up num / max / sum = {} / {} / {} primal infeasibilities\n",
                ekk.info_.num_primal_infeasibilities,
                ekk.info_.max_primal_infeasibility,
                ekk.info_.sum_primal_infeasibilities
            );
            ekk.compute_primal_objective_value();
            // Use dual to clean up. This almost always yields optimality,
            // and shouldn't yield infeasibility - since the current point
            // is dual feasible - but can yield
            // unboundedness. Time/iteration limit return is, of course,
            // possible, as are solver error
            let mut return_status = HighsStatus::kOk;
            // Switch off any bound perturbation
            let save_dual_simplex_cost_perturbation_multiplier =
                ekk.info_.dual_simplex_cost_perturbation_multiplier;
            ekk.info_.dual_simplex_cost_perturbation_multiplier = 0.0;
            let simplex_strategy = ekk.info_.simplex_strategy;
            ekk.info_.simplex_strategy = kSimplexStrategyDualPlain;
            let mut dual_solver = HEkkDual::new(ekk);
            let call_status = dual_solver.solve(ekk, true);
            // Restore any bound perturbation
            ekk.info_.dual_simplex_cost_perturbation_multiplier =
                save_dual_simplex_cost_perturbation_multiplier;
            ekk.info_.simplex_strategy = simplex_strategy;
            debug_assert!(ekk.called_return_from_solve_);
            return_status = interpret_call_status(
                &ekk.options_.log_options,
                call_status,
                return_status,
                "HEkkDual::solve",
            );
            // Reset called_return_from_solve_ to be false, since it's
            // called for this solve
            ekk.called_return_from_solve_ = false;
            if return_status != HighsStatus::kOk {
                return ekk.return_from_solve(return_status);
            }
            if ekk.model_status_ == HighsModelStatus::kOptimal
                && ekk.info_.num_primal_infeasibilities + ekk.info_.num_dual_infeasibilities != 0
            {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "HEkkPrimal:: Dual simplex clean up yields  optimality, but with {} (max {}) primal infeasibilities and {} (max {}) dual infeasibilities\n",
                    ekk.info_.num_primal_infeasibilities,
                    ekk.info_.max_primal_infeasibility,
                    ekk.info_.num_dual_infeasibilities,
                    ekk.info_.max_dual_infeasibility
                );
            }
        }
        if ekk.debug_ok_for_solve(SimplexAlgorithm::kPrimal, self.solve_phase)
            == HighsDebugStatus::kLogicalError
        {
            return ekk.return_from_solve(HighsStatus::kError);
        }
        ekk.return_from_solve(HighsStatus::kOk)
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::initialiseInstance (lines 329-375, 47 lines)
    pub fn initialise_instance(&mut self, ekk: &mut HEkk) {
        // Called in constructor for HEkkPrimal class
        // analysis = &ekk.analysis_: pointer alias removed (pattern a); use ekk.analysis_ directly.

        self.num_col = ekk.lp_.num_col_;
        self.num_row = ekk.lp_.num_row_;
        self.num_tot = self.num_col + self.num_row;

        // Setup local vectors
        self.col_aq.setup(self.num_row);
        self.row_ep.setup(self.num_row);
        self.row_ap.setup(self.num_col);
        self.col_basic_feasibility_change.setup(self.num_row);
        self.row_basic_feasibility_change.setup(self.num_col);
        self.col_steepest_edge.setup(self.num_row);

        self.ph1_sorter_r.reserve(self.num_row as usize);
        self.ph1_sorter_t.reserve(self.num_row as usize);

        self.num_free_col = 0;
        for i_col in 0..self.num_tot {
            if ekk.info_.work_lower_[i_col as usize] == -kHighsInf
                && ekk.info_.work_upper_[i_col as usize] == kHighsInf
            {
                // Free column
                self.num_free_col += 1;
            }
        }
        // Set up the HSet instances, possibly using the internal error reporting and
        // debug option
        let debug = ekk.options_.highs_debug_level > kHighsDebugLevelCheap;
        if self.num_free_col != 0 {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "HEkkPrimal:: LP has {} free columns\n",
                self.num_free_col
            );
            self.nonbasic_free_col_set.setup(
                self.num_free_col,
                self.num_tot,
                ekk.options_.output_flag,
                ekk.options_.log_options.log_stream,
                debug,
                true,
            );
        }
        // Set up the hyper-sparse CHUZC data
        self.hyper_chuzc_candidate
            .resize((1 + self.max_num_hyper_chuzc_candidates) as usize, 0);
        self.hyper_chuzc_measure
            .resize((1 + self.max_num_hyper_chuzc_candidates) as usize, 0.0);
        self.hyper_chuzc_candidate_set.setup(
            self.max_num_hyper_chuzc_candidates,
            self.num_tot,
            ekk.options_.output_flag,
            ekk.options_.log_options.log_stream,
            debug,
            true,
        );
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::initialiseSolve (lines 377-421, 45 lines)
    pub fn initialise_solve(&mut self, ekk: &mut HEkk) {
        // Copy values of simplex solver options to dual simplex options
        self.primal_feasibility_tolerance = ekk.options_.primal_feasibility_tolerance;
        self.dual_feasibility_tolerance = ekk.options_.dual_feasibility_tolerance;
        self.objective_target = ekk.options_.objective_target;

        ekk.status_.has_primal_objective_value = false;
        ekk.status_.has_dual_objective_value = false;

        ekk.model_status_ = HighsModelStatus::kNotset;
        ekk.solve_bailout_ = false;
        ekk.called_return_from_solve_ = false;
        ekk.exit_algorithm_ = SimplexAlgorithm::kPrimal;

        self.rebuild_reason = kRebuildReasonNo;
        if !ekk.status_.has_dual_steepest_edge_weights {
            // No dual weights to maintain, so ensure that the vectors are
            // assigned since they are used around factorization and when
            // setting up the backtracking information. ToDo Eliminate this
            // opacity
            ekk.dual_edge_weight_ = vec![1.0; self.num_row as usize];
            ekk.scattered_dual_edge_weight_
                .resize(self.num_tot as usize, 0.0);
        }
        let edge_weight_strategy = ekk.options_.simplex_primal_edge_weight_strategy;
        if edge_weight_strategy == kSimplexEdgeWeightStrategyChoose
            || edge_weight_strategy == kSimplexEdgeWeightStrategyDevex
        {
            // By default, use Devex
            self.edge_weight_mode = EdgeWeightMode::kDevex;
        } else if edge_weight_strategy == kSimplexEdgeWeightStrategyDantzig {
            self.edge_weight_mode = EdgeWeightMode::kDantzig;
        } else {
            debug_assert!(edge_weight_strategy == kSimplexEdgeWeightStrategySteepestEdge);
            self.edge_weight_mode = EdgeWeightMode::kSteepestEdge;
        }
        if self.edge_weight_mode == EdgeWeightMode::kDantzig {
            self.edge_weight_ = vec![1.0; self.num_tot as usize];
        } else if self.edge_weight_mode == EdgeWeightMode::kDevex {
            self.initialise_devex_framework(ekk);
        } else if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            self.compute_primal_steepest_edge_weights(ekk);
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::solvePhase1 (lines 423-513, 91 lines)
    pub fn solve_phase1(&mut self, ekk: &mut HEkk) {
        // When starting a new phase the (updated) primal objective function
        // value isn't known. Indicate this so that when the value
        // computed from scratch in build() isn't checked against the
        // updated value
        ekk.status_.has_primal_objective_value = false;
        ekk.status_.has_dual_objective_value = false;
        // Possibly bail out immediately if iteration limit is current value
        if ekk.bailout() {
            return;
        }
        crate::highs_log_dev!(
            &ekk.options_.log_options,
            HighsLogType::kDetailed,
            "primal-phase1-start\n"
        );
        // If there's no backtracking basis, save the initial basis in case of
        // backtracking
        if !ekk.info_.valid_backtracking_basis_ {
            ekk.put_backtracking_basis();
        }

        // Main solving structure
        loop {
            //
            // Rebuild
            //
            // solve_phase = kSolvePhaseError is set if the basis matrix is singular
            self.rebuild(ekk);
            if self.solve_phase == kSolvePhaseError {
                return;
            }
            if self.solve_phase == kSolvePhaseUnknown {
                return;
            }
            if ekk.bailout() {
                return;
            }
            debug_assert!(self.solve_phase == kSolvePhase1 || self.solve_phase == kSolvePhase2);
            //
            // solve_phase = kSolvePhase2 is set if no primal infeasibilities
            // are found in rebuild(), in which case return for phase 2
            if self.solve_phase == kSolvePhase2 {
                break;
            }

            loop {
                self.iterate(ekk);
                if ekk.bailout() {
                    return;
                }
                if self.solve_phase == kSolvePhaseError {
                    return;
                }
                debug_assert!(self.solve_phase == kSolvePhase1);
                if self.rebuild_reason != 0 {
                    break;
                }
            }
            // If the data are fresh from rebuild() and no flips have
            // occurred, possibly break out of the outer loop to see what's
            // occurred
            let finished = ekk.status_.has_fresh_rebuild
                && self.num_flip_since_rebuild == 0
                && !ekk.rebuild_refactor(self.rebuild_reason);
            if finished && ekk.taboo_bad_basis_change() {
                // A bad basis change has had to be made taboo without any other
                // basis changes or flips having been performed from a fresh
                // rebuild. In other words, the only basis change that could be
                // made is not permitted, so no definitive statement about the
                // LP can be made.
                self.solve_phase = kSolvePhaseTabooBasis;
                return;
            }
            if finished {
                break;
            }
        }
        // If bailing out, should have returned already
        debug_assert!(!ekk.solve_bailout_);
        // Will only have accurate simplex info if moving to phase 2 - but
        // should check primal feasibility and residual information if LP
        // is primal infeasible
        if self.debug_primal_simplex(ekk, "End of solvePhase1", false)
            == HighsDebugStatus::kLogicalError
        {
            self.solve_phase = kSolvePhaseError;
            return;
        }
        if self.solve_phase == kSolvePhase1 {
            // Determine whether primal infeasibility has been identified
            if self.variable_in < 0 {
                // Optimal in phase 1, so should have primal infeasibilities
                debug_assert!(ekk.info_.num_primal_infeasibilities > 0);
                if ekk.info_.bounds_shifted || ekk.info_.bounds_perturbed {
                    // Remove any bound shifts or perturbations and return to
                    // phase 1
                    self.cleanup(ekk);
                } else {
                    ekk.model_status_ = HighsModelStatus::kInfeasible;
                    self.solve_phase = kSolvePhaseExit;
                }
            }
        }
        if self.solve_phase == kSolvePhase2 {
            // Moving to phase 2 so comment if bound perturbation is not permitted
            //
            // It may have been prevented to avoid cleanup-perturbation loops
            if !ekk.info_.allow_bound_perturbation {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "Moving to phase 2, but not allowing bound perturbation\n"
                );
            }
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::solvePhase2 (lines 515-645, 131 lines)
    pub fn solve_phase2(&mut self, ekk: &mut HEkk) {
        // When starting a new phase the (updated) primal objective function
        // value isn't known. Indicate this so that when the value
        // computed from scratch in build() isn't checked against the
        // updated value
        ekk.status_.has_primal_objective_value = false;
        ekk.status_.has_dual_objective_value = false;
        // Possibly bail out immediately if iteration limit is current value
        if ekk.bailout() {
            return;
        }
        crate::highs_log_dev!(
            &ekk.options_.log_options,
            HighsLogType::kDetailed,
            "primal-phase2-start\n"
        );
        self.phase2_update_primal(ekk, true);

        // If there's no backtracking basis Save the initial basis in case of
        // backtracking
        if !ekk.info_.valid_backtracking_basis_ {
            ekk.put_backtracking_basis();
        }

        // Main solving structure
        loop {
            //
            // Rebuild
            //
            // solve_phase = kSolvePhaseError is set if the basis matrix is singular
            self.rebuild(ekk);
            if self.solve_phase == kSolvePhaseError {
                return;
            }
            if self.solve_phase == kSolvePhaseUnknown {
                return;
            }
            if ekk.bailout() {
                return;
            }
            debug_assert!(self.solve_phase == kSolvePhase1 || self.solve_phase == kSolvePhase2);
            //
            // solve_phase = kSolvePhase1 is set if primal infeasibilities
            // are found in rebuild(), in which case return for phase 1
            if self.solve_phase == kSolvePhase1 {
                break;
            }

            loop {
                self.iterate(ekk);
                if ekk.bailout() {
                    return;
                }
                if self.solve_phase == kSolvePhaseError {
                    return;
                }
                debug_assert!(self.solve_phase == kSolvePhase2);
                if self.rebuild_reason != 0 {
                    break;
                }
            }
            // If the data are fresh from rebuild() and no flips have
            // occurred, possibly break out of the outer loop to see what's
            // occurred
            let finished = ekk.status_.has_fresh_rebuild
                && self.num_flip_since_rebuild == 0
                && !ekk.rebuild_refactor(self.rebuild_reason);
            if finished && ekk.taboo_bad_basis_change() {
                // A bad basis change has had to be made taboo without any other
                // basis changes or flips having been performed from a fresh
                // rebuild. In other words, the only basis change that could be
                // made is not permitted, so no definitive statement about the
                // LP can be made.
                self.solve_phase = kSolvePhaseTabooBasis;
                return;
            }
            if finished {
                break;
            }
        }
        // If bailing out, should have returned already
        debug_assert!(!ekk.solve_bailout_);
        if self.debug_primal_simplex(ekk, "End of solvePhase2", false)
            == HighsDebugStatus::kLogicalError
        {
            self.solve_phase = kSolvePhaseError;
            return;
        }
        if self.solve_phase == kSolvePhase1 {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "primal-return-phase1\n"
            );
        } else if self.variable_in == -1 {
            // There is no candidate in CHUZC, even after rebuild so probably optimal
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "primal-phase-2-optimal\n"
            );
            // Remove any bound perturbations and see if basis is still primal feasible
            self.cleanup(ekk);
            if ekk.info_.num_primal_infeasibilities > 0 {
                // There are primal infeasibilities, so consider performing dual
                // simplex iterations to get primal feasibility
                self.solve_phase = kSolvePhaseOptimalCleanup;
            } else {
                // There are no primal infeasibilities so optimal!
                self.solve_phase = kSolvePhaseOptimal;
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kDetailed,
                    "problem-optimal\n"
                );
                ekk.model_status_ = HighsModelStatus::kOptimal;
                ekk.compute_dual_objective_value(2); // Why?
            }
        } else if self.row_out == kNoRowSought {
            // CHUZR has not been performed - because the chosen reduced cost
            // was unattractive when computed from scratch and no rebuild was
            // required. This is very rare and should be handled otherwise
            //
            println!(
                "HEkkPrimal::solvePhase2 row_out = {} solve {}",
                self.row_out, ekk.debug_solve_call_num_
            );
            debug_assert!(self.row_out != kNoRowSought);
        } else {
            // No candidate in CHUZR
            if self.row_out >= 0 {
                println!(
                    "HEkkPrimal::solvePhase2 row_out = {} solve {}",
                    self.row_out, ekk.debug_solve_call_num_
                );
            }
            // Ensure that CHUZR was performed and found no row
            debug_assert!(self.row_out == kNoRowChosen);

            // There is no candidate in CHUZR, so probably primal unbounded
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "primal-phase-2-unbounded\n"
            );
            if ekk.info_.bounds_shifted || ekk.info_.bounds_perturbed {
                // If the bounds have been shifted or perturbed, clean up and
                // return
                self.cleanup(ekk);
                // If there are primal infeasibilities, go back to phase 1
                if ekk.info_.num_primal_infeasibilities > 0 {
                    self.solve_phase = kSolvePhase1;
                }
            } else {
                // The bounds have not been perturbed, so primal unbounded
                self.solve_phase = kSolvePhaseExit;
                // Primal unbounded, so save primal ray
                self.save_primal_ray(ekk);
                // Model status should be unset
                debug_assert!(ekk.model_status_ == HighsModelStatus::kNotset);
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kInfo,
                    "problem-primal-unbounded\n"
                );
                ekk.model_status_ = HighsModelStatus::kUnbounded;
            }
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::cleanup (lines 647-684, 38 lines)
    pub fn cleanup(&mut self, ekk: &mut HEkk) {
        if !ekk.info_.bounds_shifted && !ekk.info_.bounds_perturbed {
            return;
        }
        crate::highs_log_dev!(
            &ekk.options_.log_options,
            HighsLogType::kDetailed,
            "primal-cleanup-shift\n"
        );
        // Remove perturbation and don't permit further perturbation
        ekk.initialise_bound(SimplexAlgorithm::kPrimal, self.solve_phase, false);
        ekk.initialise_nonbasic_value_and_move();
        ekk.info_.allow_bound_perturbation = false;
        // Possibly take a copy of the original duals before recomputing them
        // (HiGHS: debugging copy of baseValue_ is commented out in the C++ source.)
        // Compute the primal values
        ekk.compute_primal();
        // Possibly analyse the change in duals
        // (HiGHS: debugCleanup is commented out in the C++ source.)
        // Compute the primal infeasibilities
        ekk.compute_simplex_primal_infeasible();

        // Compute the primal objective value
        ekk.compute_primal_objective_value();
        // Now that there's a new primal_objective_value, reset the updated
        // value
        ekk.info_.updated_primal_objective_value = ekk.info_.primal_objective_value;

        // Report the dual infeasibilities
        ekk.compute_simplex_dual_infeasible();
        // In phase 1, report the simplex LP dual infeasibilities
        // In phase 2, report the simplex dual infeasibilities (known)
        self.report_rebuild(ekk, kRebuildReasonCleanup);
    }
}
