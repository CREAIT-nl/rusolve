//! HiGHS source: simplex/HEkkDual.cpp (lines 2341-3003)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Dual infeasibility correction (flips/shifts), phase-1 optimality assessment, objective
//! bound checks, unboundedness and infeasibility proofs.
//! See `HEkkDual` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkDual {
    /// HiGHS: HEkkDual.cpp::HEkkDual::computeDualInfeasibilitiesWithFixedVariableFlips (lines 2341-2386, 46 lines)
    pub fn compute_dual_infeasibilities_with_fixed_variable_flips(&mut self, ekk: &mut HEkk) {
        let mut num_dual_infeasibility: HighsInt = 0;
        let mut max_dual_infeasibility: f64 = 0.0;
        let mut sum_dual_infeasibility: f64 = 0.0;
        let num_tot = ekk.lp_.num_col_ + ekk.lp_.num_row_;
        for i_var in 0..num_tot {
            let iv = i_var as usize;
            if ekk.basis_.nonbasic_flag_[iv] == 0 {
                continue;
            }
            // Nonbasic column
            let lower = ekk.info_.work_lower_[iv];
            let upper = ekk.info_.work_upper_[iv];
            let dual = ekk.info_.work_dual_[iv];

            let dual_infeasibility: f64 = if lower == -kHighsInf && upper == kHighsInf {
                dual.abs()
            } else {
                (-(ekk.basis_.nonbasic_move_[iv] as HighsInt)) as f64 * dual
            };
            if dual_infeasibility > 0.0 {
                if dual_infeasibility >= ekk.options_.dual_feasibility_tolerance {
                    num_dual_infeasibility += 1;
                }
                max_dual_infeasibility = cpp_max(dual_infeasibility, max_dual_infeasibility);
                sum_dual_infeasibility += dual_infeasibility;
            }
        }
        ekk.info_.num_dual_infeasibilities = num_dual_infeasibility;
        ekk.info_.max_dual_infeasibility = max_dual_infeasibility;
        ekk.info_.sum_dual_infeasibilities = sum_dual_infeasibility;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::correctDualInfeasibilities (lines 2388-2539, 152 lines)
    /// MUST-PORT (determinism, random stream): `random.fraction()` (on `ekk.random_`) per cost shift, in loop order; same stream as CHUZR.
    pub fn correct_dual_infeasibilities(
        &mut self,
        ekk: &mut HEkk,
        free_infeasibility_count: &mut HighsInt,
    ) {
        *free_infeasibility_count = 0;
        let dual_feasibility_tolerance = ekk.options_.dual_feasibility_tolerance;
        let mut flip_dual_objective_value_change: f64 = 0.0;
        let mut shift_dual_objective_value_change: f64 = 0.0;
        let mut num_flip: HighsInt = 0;
        let mut num_shift: HighsInt = 0;
        let mut sum_flip: f64 = 0.0;
        let mut sum_shift: f64 = 0.0;
        let mut max_flip: f64 = 0.0;
        let mut max_shift: f64 = 0.0;
        let mut min_dual_infeasibility_for_flip: f64 = kHighsInf;
        let mut max_dual_infeasibility_for_flip: f64 = 0.0;
        let mut num_dual_infeasibilities_for_flip: HighsInt = 0;
        let mut sum_dual_infeasibilities_for_flip: f64 = 0.0;
        let mut num_dual_infeasibilities_for_shift: HighsInt = 0;
        let mut max_dual_infeasibility_for_shift: f64 = 0.0;
        let mut sum_dual_infeasibilities_for_shift: f64 = 0.0;
        let num_tot = ekk.lp_.num_col_ + ekk.lp_.num_row_;
        for i_var in 0..num_tot {
            let iv = i_var as usize;
            if ekk.basis_.nonbasic_flag_[iv] == 0 {
                continue;
            }
            // Nonbasic column
            let lower = ekk.info_.work_lower_[iv];
            let upper = ekk.info_.work_upper_[iv];
            let current_dual = ekk.info_.work_dual_[iv];
            let mv: HighsInt = ekk.basis_.nonbasic_move_[iv] as HighsInt;
            let fixed = lower == upper;
            let boxed = lower > -kHighsInf && upper < kHighsInf;
            let free = lower == -kHighsInf && upper == kHighsInf;
            if free {
                let free_dual_infeasibility = current_dual.abs();
                if free_dual_infeasibility >= dual_feasibility_tolerance {
                    *free_infeasibility_count += 1;
                }
                continue;
            }
            let dual_infeasibility: f64 = (-mv) as f64 * current_dual;
            if dual_infeasibility < dual_feasibility_tolerance {
                continue;
            }
            // There is a dual infeasibility to remove
            //
            // force_phase2 is set true to prevent flipping of non-fixed
            // (boxed) variables when correcting infeasibilities in the first
            // set of duals computed after cost perturbation
            if fixed || (boxed && !self.force_phase2) {
                // Flip for fixed variables and boxed variables when not forcing phase 2
                ekk.flip_bound(i_var);
                let flip = upper - lower;
                let mut local_dual_objective_change = mv as f64 * flip * current_dual;
                local_dual_objective_change *= ekk.cost_scale_;
                flip_dual_objective_value_change += local_dual_objective_change;
                num_flip += 1;
                max_flip = cpp_max(flip.abs(), max_flip);
                sum_flip += flip.abs();
                // Flipping fixed variables is trivial, so only track the
                // infeasibilities involved when flipping boxed variables
                if !fixed {
                    min_dual_infeasibility_for_flip =
                        cpp_min(dual_infeasibility, min_dual_infeasibility_for_flip);
                    if dual_infeasibility >= dual_feasibility_tolerance {
                        num_dual_infeasibilities_for_flip += 1;
                    }
                    sum_dual_infeasibilities_for_flip += dual_infeasibility;
                    max_dual_infeasibility_for_flip =
                        cpp_max(dual_infeasibility, max_dual_infeasibility_for_flip);
                }
                continue;
            }
            // Either boxed but not fixed, of one-sided, so shift
            //
            // Cost shifting must always be possible
            debug_assert!(ekk.info_.allow_cost_shifting);
            // Other variable = shift
            if dual_infeasibility >= dual_feasibility_tolerance {
                num_dual_infeasibilities_for_shift += 1;
            }
            sum_dual_infeasibilities_for_shift += dual_infeasibility;
            max_dual_infeasibility_for_shift =
                cpp_max(dual_infeasibility, max_dual_infeasibility_for_shift);
            ekk.info_.costs_shifted = true;
            let shift: f64;
            if mv == kNonbasicMoveUp as HighsInt {
                // MUST-PORT: one random draw per shift
                let new_dual = (1.0 + ekk.random_.fraction()) * dual_feasibility_tolerance;
                shift = new_dual - current_dual;
                ekk.info_.work_dual_[iv] = new_dual;
                ekk.info_.work_cost_[iv] = ekk.info_.work_cost_[iv] + shift;
            } else {
                // MUST-PORT: one random draw per shift
                let new_dual = -(1.0 + ekk.random_.fraction()) * dual_feasibility_tolerance;
                shift = new_dual - current_dual;
                ekk.info_.work_dual_[iv] = new_dual;
                ekk.info_.work_cost_[iv] = ekk.info_.work_cost_[iv] + shift;
            }
            let mut local_dual_objective_change = shift * ekk.info_.work_value_[iv];
            local_dual_objective_change *= ekk.cost_scale_;
            shift_dual_objective_value_change += local_dual_objective_change;
            num_shift += 1;
            max_shift = cpp_max(shift.abs(), max_shift);
            sum_shift += shift.abs();
            let direction = if mv == kNonbasicMoveUp as HighsInt {
                "  up"
            } else {
                "down"
            };
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kVerbose,
                "Move {}: cost shift = {}; objective change = {}\n",
                direction,
                shift,
                local_dual_objective_change
            );
        }
        ekk.analysis_.num_correct_dual_primal_flip += num_flip;
        ekk.analysis_.max_correct_dual_primal_flip =
            cpp_max(max_flip, ekk.analysis_.max_correct_dual_primal_flip);
        ekk.analysis_
            .min_correct_dual_primal_flip_dual_infeasibility = cpp_min(
            min_dual_infeasibility_for_flip,
            ekk.analysis_
                .min_correct_dual_primal_flip_dual_infeasibility,
        );
        if num_flip != 0 && self.force_phase2 {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "Performed num / max / sum = {} / {} / {} flip(s) for num / min / max / sum dual infeasibility of {} / {} / {} / {}; objective change = {}\n",
                num_flip,
                max_flip,
                sum_flip,
                num_dual_infeasibilities_for_flip,
                min_dual_infeasibility_for_flip,
                max_dual_infeasibility_for_flip,
                sum_dual_infeasibilities_for_flip,
                flip_dual_objective_value_change
            );
        }
        ekk.analysis_.num_correct_dual_cost_shift += num_shift;
        ekk.analysis_.max_correct_dual_cost_shift =
            cpp_max(max_shift, ekk.analysis_.max_correct_dual_cost_shift);
        ekk.analysis_.max_correct_dual_cost_shift_dual_infeasibility = cpp_max(
            max_dual_infeasibility_for_shift,
            ekk.analysis_.max_correct_dual_cost_shift_dual_infeasibility,
        );
        if num_shift != 0 {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "Performed num / max / sum = {} / {} / {} shift(s) for num / max / sum dual infeasibility of {} / {} / {}; objective change = {}\n",
                num_shift,
                max_shift,
                sum_shift,
                num_dual_infeasibilities_for_shift,
                max_dual_infeasibility_for_shift,
                sum_dual_infeasibilities_for_shift,
                shift_dual_objective_value_change
            );
        }
        self.force_phase2 = false;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::proofOfPrimalInfeasibility (lines 2541-2543, 3 lines)
    pub fn proof_of_primal_infeasibility(&mut self, ekk: &mut HEkk) -> bool {
        let mut row_ep = std::mem::take(&mut self.row_ep);
        let result = ekk.proof_of_primal_infeasibility_ep(&mut row_ep, self.move_out, self.row_out);
        self.row_ep = row_ep;
        result
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::saveDualRay (lines 2545-2551, 7 lines)
    pub fn save_dual_ray(&mut self, ekk: &mut HEkk) {
        debug_assert!(self.row_out >= 0);
        debug_assert!(self.move_out != kNoRaySign);
        ekk.dual_ray_record_.clear();
        ekk.dual_ray_record_.index = self.row_out;
        ekk.dual_ray_record_.sign = self.move_out;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::assessPhase1Optimality (lines 2553-2601, 49 lines)
    pub fn assess_phase1_optimality(&mut self, ekk: &mut HEkk) {
        // Should only be called when optimal in phase 1 (row_out == kNoRowChosen)
        // with nonzero dual activity, and after a fresh rebuild - so
        // "final" decisions can be made.
        debug_assert!(self.solve_phase == kSolvePhase1);
        debug_assert!(self.row_out == kNoRowChosen);
        debug_assert!(ekk.info_.dual_objective_value != 0.0);
        debug_assert!(ekk.status_.has_fresh_rebuild);

        crate::highs_log_dev!(
            &ekk.options_.log_options,
            HighsLogType::kInfo,
            "Optimal in phase 1 but not jumping to phase 2 since dual objective is {:10.4e}: Costs perturbed = {}\n",
            ekk.info_.dual_objective_value,
            ekk.info_.costs_perturbed
        );
        if ekk.info_.costs_perturbed {
            // Clean up perturbation
            self.cleanup(ekk);
            self.assess_phase1_optimality_unperturbed(ekk);
        } else {
            debug_assert!(self.dual_infeas_count == 0);
            debug_assert!(ekk.info_.dual_objective_value != 0.0);
            self.assess_phase1_optimality_unperturbed(ekk);
        }
        if self.dual_infeas_count > 0 {
            // Must still be solve_phase = kSolvePhase1 since dual
            // infeasibilities with respect to phase 1 bounds mean that primal
            // values must change, so primal feasibility is unknown
            debug_assert!(self.solve_phase == kSolvePhase1);
        } else {
            // Optimal in dual phase 1, so either dual feasible wrt Phase 2
            // bounds and going to phase 2, or identified dual infeasibility and exiting
            debug_assert!(
                self.solve_phase == kSolvePhase2
                    || (self.solve_phase == kSolvePhaseExit
                        && ekk.model_status_ == HighsModelStatus::kUnboundedOrInfeasible)
            );
            if self.solve_phase == kSolvePhase2 {
                // Reset the duals, if necessary shifting costs of free variables
                // so that their duals are zero
                self.exit_phase1_reset_duals(ekk);
            }
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::assessPhase1OptimalityUnperturbed (lines 2603-2654, 52 lines)
    pub fn assess_phase1_optimality_unperturbed(&mut self, ekk: &mut HEkk) {
        debug_assert!(!ekk.info_.costs_perturbed);
        if self.dual_infeas_count == 0 {
            // No dual infeasibilities with respect to phase 1 bounds.
            if ekk.info_.dual_objective_value == 0.0 {
                // No dual infeasibilities with respect to phase 2 bounds so
                // go to phase 2
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kInfo,
                    "LP is dual feasible wrt Phase 2 bounds after removing cost perturbations so go to phase 2\n"
                );
                self.solve_phase = kSolvePhase2;
            } else {
                // Nonzero dual objective value: could be insignificant dual
                // infeasibilities
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kInfo,
                    "LP is dual feasible wrt Phase 1 bounds after removing cost perturbations: dual objective is {:10.4e}\n",
                    ekk.info_.dual_objective_value
                );
                ekk.compute_simplex_lp_dual_infeasible();
                let num_lp_dual_infeasibilities =
                    ekk.analysis_.num_dual_phase_1_lp_dual_infeasibility;
                if num_lp_dual_infeasibilities == 0 {
                    crate::highs_log_dev!(
                        &ekk.options_.log_options,
                        HighsLogType::kInfo,
                        "LP is dual feasible wrt Phase 2 bounds after removing cost perturbations so go to phase 2\n"
                    );
                    self.solve_phase = kSolvePhase2;
                } else {
                    // LP is dual infeasible if the dual objective is sufficiently
                    // negative, so no conclusions on the primal LP can be deduced
                    // - could be primal unbounded or primal infeasible.
                    //
                    // Indicate the conclusion of dual infeasibility by setting
                    // the scaled model status
                    self.report_on_possible_lp_dual_infeasibility(ekk);
                    ekk.model_status_ = HighsModelStatus::kUnboundedOrInfeasible;
                    self.solve_phase = kSolvePhaseExit;
                }
            }
        } else {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "LP has {} dual feasibilities wrt Phase 1 bounds after removing cost perturbations so return to phase 1\n",
                self.dual_infeas_count
            );
            debug_assert!(self.solve_phase == kSolvePhase1);
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::exitPhase1ResetDuals (lines 2656-2709, 54 lines)
    pub fn exit_phase1_reset_duals(&mut self, ekk: &mut HEkk) {
        if ekk.info_.costs_perturbed {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "Costs are already perturbed in exitPhase1ResetDuals\n"
            );
        } else {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "Re-perturbing costs when optimal in phase 1\n"
            );
            ekk.initialise_cost(SimplexAlgorithm::kDual, kSolvePhase2, true);
            ekk.compute_dual();
        }

        let num_tot = ekk.lp_.num_col_ + ekk.lp_.num_row_;
        let mut num_shift: HighsInt = 0;
        let mut sum_shift: f64 = 0.0;
        for i_var in 0..num_tot {
            let iv = i_var as usize;
            if ekk.basis_.nonbasic_flag_[iv] != 0 {
                let lp_lower: f64;
                let lp_upper: f64;
                if i_var < ekk.lp_.num_col_ {
                    lp_lower = ekk.lp_.col_lower_[iv];
                    lp_upper = ekk.lp_.col_upper_[iv];
                } else {
                    let i_row = i_var - ekk.lp_.num_col_;
                    lp_lower = ekk.lp_.row_lower_[i_row as usize];
                    lp_upper = ekk.lp_.row_upper_[i_row as usize];
                }
                if lp_lower <= -kHighsInf && lp_upper >= kHighsInf {
                    let shift = -ekk.info_.work_dual_[iv];
                    ekk.info_.work_dual_[iv] = 0.0;
                    ekk.info_.work_cost_[iv] = ekk.info_.work_cost_[iv] + shift;
                    num_shift += 1;
                    sum_shift += shift.abs();
                    crate::highs_log_dev!(
                        &ekk.options_.log_options,
                        HighsLogType::kVerbose,
                        "Variable {} is free: shift cost to zero dual of {}\n",
                        i_var,
                        shift
                    );
                }
            }
        }
        if num_shift != 0 {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "Performed {} cost shift(s) for free variables to zero dual values: total = {}\n",
                num_shift,
                sum_shift
            );
            ekk.info_.costs_shifted = true;
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::reportOnPossibleLpDualInfeasibility (lines 2711-2732, 22 lines)
    pub fn report_on_possible_lp_dual_infeasibility(&mut self, ekk: &mut HEkk) {
        debug_assert!(self.solve_phase == kSolvePhase1);
        debug_assert!(self.row_out == kNoRowChosen);
        debug_assert!(!ekk.info_.costs_perturbed);
        let lp_dual_status: &str = if ekk.analysis_.num_dual_phase_1_lp_dual_infeasibility != 0 {
            "infeasible"
        } else {
            "feasible"
        };
        crate::highs_log_dev!(
            &ekk.options_.log_options,
            HighsLogType::kInfo,
            "LP is dual {} with dual phase 1 objective {:10.4e} and num / max / sum dual infeasibilities = {} / {:9.4e} / {:9.4e}\n",
            lp_dual_status,
            ekk.info_.dual_objective_value,
            ekk.analysis_.num_dual_phase_1_lp_dual_infeasibility,
            ekk.analysis_.max_dual_phase_1_lp_dual_infeasibility,
            ekk.analysis_.sum_dual_phase_1_lp_dual_infeasibility
        );
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::dualInfoOk (lines 2734-2761, 28 lines)
    pub fn dual_info_ok(&self, ekk: &HEkk, lp: &HighsLp) -> bool {
        let lp_num_col = lp.num_col_;
        let lp_num_row = lp.num_row_;
        let mut dimensions_ok: bool;
        dimensions_ok = lp_num_col == self.solver_num_col && lp_num_row == self.solver_num_row;
        debug_assert!(dimensions_ok);
        if !dimensions_ok {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kError,
                "LP-Solver dimension incompatibility ({}, {}) != ({}, {})\n",
                lp_num_col,
                self.solver_num_col,
                lp_num_row,
                self.solver_num_row
            );
            return false;
        }
        // simplex_nla->lp_ is the HEkk's LP (NlaRefs::lp_ is built from ekk.lp_)
        dimensions_ok = lp_num_col == ekk.lp_.num_col_ && lp_num_row == ekk.lp_.num_row_;
        debug_assert!(dimensions_ok);
        if !dimensions_ok {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kError,
                "LP-Factor dimension incompatibility ({}, {}) != ({}, {})\n",
                lp_num_col,
                ekk.lp_.num_col_,
                lp_num_row,
                ekk.lp_.num_row_
            );
            return false;
        }
        true
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::bailoutOnDualObjective (lines 2763-2777, 15 lines)
    pub fn bailout_on_dual_objective(&mut self, ekk: &mut HEkk) -> bool {
        if ekk.solve_bailout_ {
            // Bailout has already been decided: check that it's for one of these
            // reasons
            debug_assert!(
                ekk.model_status_ == HighsModelStatus::kTimeLimit
                    || ekk.model_status_ == HighsModelStatus::kIterationLimit
                    || ekk.model_status_ == HighsModelStatus::kObjectiveBound
            );
        } else if ekk.lp_.sense_ == ObjSense::kMinimize && self.solve_phase == kSolvePhase2 {
            if ekk.info_.updated_dual_objective_value > ekk.options_.objective_bound {
                ekk.solve_bailout_ = self.reached_exact_objective_bound(ekk);
            }
        }
        ekk.solve_bailout_
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::reachedExactObjectiveBound (lines 2779-2849, 71 lines)
    pub fn reached_exact_objective_bound(&mut self, ekk: &mut HEkk) -> bool {
        let mut reached_exact_objective_bound = false;
        let use_row_ap_density: f64 = cpp_min(cpp_max(ekk.info_.row_ap_density, 0.01), 1.0);
        let check_frequency: HighsInt = (1.0 / use_row_ap_density) as HighsInt;
        debug_assert!(check_frequency > 0);

        let check_exact_dual_objective_value = ekk.info_.update_count % check_frequency == 0;

        if check_exact_dual_objective_value {
            let objective_bound = ekk.options_.objective_bound;
            let perturbed_dual_objective_value = ekk.info_.updated_dual_objective_value;
            let perturbed_value_residual = perturbed_dual_objective_value - objective_bound;
            let mut dual_col = HVector::default();
            let mut dual_row = HVector::default();
            let exact_dual_objective_value =
                self.compute_exact_dual_objective_value(ekk, &mut dual_col, &mut dual_row);
            let exact_value_residual = exact_dual_objective_value - objective_bound;
            let action: &str;
            if exact_dual_objective_value > objective_bound {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kDetailed,
                    "HEkkDual::solvePhase2: {:12} = Objective > ObjectiveUB = {:12}\n",
                    ekk.info_.updated_dual_objective_value,
                    objective_bound
                );
                action = "Have DualUB bailout";
                if ekk.info_.costs_perturbed || ekk.info_.costs_shifted {
                    // Remove cost perturbation/shifting
                    ekk.initialise_cost(SimplexAlgorithm::kDual, kSolvePhase2, false);
                }

                // Set the duals as computed in the computeExactDualObjective call
                for i in 0..self.solver_num_col {
                    ekk.info_.work_dual_[i as usize] =
                        ekk.info_.work_cost_[i as usize] - dual_row.array[i as usize];
                }
                for i in self.solver_num_col..self.solver_num_tot {
                    ekk.info_.work_dual_[i as usize] =
                        -dual_col.array[(i - self.solver_num_col) as usize];
                }

                // Since the computeExactDualObjectiveValue() call succeeded, if there are
                // any dual infeasibilities they can be removed by a bound flip
                self.force_phase2 = false;
                let mut dual_infeas_count = self.dual_infeas_count;
                self.correct_dual_infeasibilities(ekk, &mut dual_infeas_count);
                self.dual_infeas_count = dual_infeas_count;

                // no shifts should have occurred
                debug_assert!(!ekk.info_.costs_shifted);
                reached_exact_objective_bound = true;
                ekk.model_status_ = HighsModelStatus::kObjectiveBound;
            } else {
                action = "No   DualUB bailout";
            }
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "{} on iteration {}: Density {:11.4e}; Frequency {}: Residual(Perturbed = {}; Exact = {})\n",
                action,
                ekk.iteration_count_,
                use_row_ap_density,
                check_frequency,
                perturbed_value_residual,
                exact_value_residual
            );
        }
        reached_exact_objective_bound
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::computeExactDualObjectiveValue (lines 2851-2961, 111 lines)
    pub fn compute_exact_dual_objective_value(
        &mut self,
        ekk: &mut HEkk,
        dual_col: &mut HVector,
        dual_row: &mut HVector,
    ) -> f64 {
        // Create a local buffer for the pi vector
        dual_col.setup(ekk.lp_.num_row_);
        dual_col.clear();
        for i_row in 0..ekk.lp_.num_row_ {
            let i_var = ekk.basis_.basic_index_[i_row as usize];
            if i_var < ekk.lp_.num_col_ {
                let value = ekk.lp_.col_cost_[i_var as usize];
                if value != 0.0 {
                    dual_col.array[i_row as usize] = value;
                    dual_col.index[dual_col.count as usize] = i_row;
                    dual_col.count += 1;
                }
            }
        }
        // Create a local buffer for the dual vector
        let num_tot = ekk.lp_.num_col_ + ekk.lp_.num_row_;
        dual_row.setup(ekk.lp_.num_col_);
        dual_row.clear();
        if dual_col.count != 0 {
            let quad_precision = false;
            let expected_density: f64 = 1.0;
            {
                let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
                ekk.simplex_nla_
                    .btran(&refs, dual_col, expected_density, None);
            }
            ekk.lp_
                .a_matrix_
                .price_by_column(quad_precision, dual_row, dual_col, kDebugReportOff);
        }
        // Compute dual infeasibilities
        ekk.compute_simplex_dual_infeasible();
        if ekk.info_.num_dual_infeasibilities > 0 {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "When computing exact dual objective, the unperturbed costs yield num / max / sum dual infeasibilities = {} / {} / {}\n",
                ekk.info_.num_dual_infeasibilities,
                ekk.info_.max_dual_infeasibility,
                ekk.info_.sum_dual_infeasibilities
            );
        }
        let mut dual_objective = HighsCDouble::from(ekk.lp_.offset_);
        let mut norm_dual: f64 = 0.0;
        let mut norm_delta_dual: f64 = 0.0;
        for i_col in 0..ekk.lp_.num_col_ {
            let ic = i_col as usize;
            if ekk.basis_.nonbasic_flag_[ic] == 0 {
                continue;
            }
            let exact_dual = ekk.lp_.col_cost_[ic] - dual_row.array[ic];
            // The active value must be decided based on the exact dual. For a nonbasic
            // column the bound that must be used may flip due to cost perturbation
            // flipping the sign of its dual and for a basic variable we may need to add
            // to the dual objective using one of the bounds when its dual is not zero.

            let active_value: f64 = if exact_dual > ekk.options_.small_matrix_value {
                ekk.lp_.col_lower_[ic]
            } else if exact_dual < -ekk.options_.small_matrix_value {
                ekk.lp_.col_upper_[ic]
            } else {
                ekk.info_.work_value_[ic]
            };

            // when the active value is infinite the dual objective lower bound is
            // -infinity
            if highs_is_infinity(active_value.abs()) {
                return -kHighsInf;
            }

            let residual = (exact_dual - ekk.info_.work_dual_[ic]).abs();
            norm_dual += exact_dual.abs();
            norm_delta_dual += residual;
            if residual > 1e10 {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "Col {:4}: ExactDual = {:11.4e}; WorkDual = {:11.4e}; Residual = {:11.4e}\n",
                    i_col,
                    exact_dual,
                    ekk.info_.work_dual_[ic],
                    residual
                );
            }
            dual_objective += active_value * exact_dual;
        }

        for i_var in ekk.lp_.num_col_..num_tot {
            let iv = i_var as usize;
            if ekk.basis_.nonbasic_flag_[iv] == 0 {
                continue;
            }
            let i_row = i_var - ekk.lp_.num_col_;
            let ir = i_row as usize;
            let exact_dual = dual_col.array[ir];
            // Similarly to the column case above the active value must be decided based
            // on the exact dual.

            let active_value: f64 = if exact_dual > ekk.options_.small_matrix_value {
                ekk.lp_.row_lower_[ir]
            } else if exact_dual < -ekk.options_.small_matrix_value {
                ekk.lp_.row_upper_[ir]
            } else {
                -ekk.info_.work_value_[iv]
            };

            // when the active value is infinite the dual objective lower bound is
            // -infinity
            if highs_is_infinity(active_value.abs()) {
                return -kHighsInf;
            }

            let residual = (exact_dual + ekk.info_.work_dual_[iv]).abs();
            norm_dual += exact_dual.abs();
            norm_delta_dual += residual;
            if residual > 1e10 {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "Row {:4}: ExactDual = {:11.4e}; WorkDual = {:11.4e}; Residual = {:11.4e}\n",
                    i_row,
                    exact_dual,
                    ekk.info_.work_dual_[iv],
                    residual
                );
            }
            dual_objective += active_value * exact_dual;
        }
        let relative_delta = norm_delta_dual / cpp_max(norm_dual, 1.0);
        if relative_delta > 1e-3 {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kWarning,
                "||exact dual vector|| = {}; ||delta dual vector|| = {}: ratio = {}\n",
                norm_dual,
                norm_delta_dual,
                relative_delta
            );
        }
        f64::from(dual_objective)
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::debugDualSimplex (lines 2963-2970, 8 lines)
    /// C++ default arguments: initialise = false
    /// Debug check: out of scope; stub returns HighsDebugStatus::kNotChecked (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_dual_simplex(
        &mut self,
        ekk: &mut HEkk,
        message: &str,
        initialise: bool,
    ) -> HighsDebugStatus {
        HighsDebugStatus::kNotChecked
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::isBadBasisChange (lines 2972-2975, 4 lines)
    pub fn is_bad_basis_change(&mut self, ekk: &mut HEkk) -> bool {
        ekk.is_bad_basis_change(
            SimplexAlgorithm::kDual,
            self.variable_in,
            self.row_out,
            self.rebuild_reason,
        )
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::assessPossiblyDualUnbounded (lines 2977-3003, 27 lines)
    pub fn assess_possibly_dual_unbounded(&mut self, ekk: &mut HEkk) {
        debug_assert!(self.rebuild_reason == kRebuildReasonPossiblyDualUnbounded);
        if self.solve_phase != kSolvePhase2 {
            return;
        }
        if !ekk.status_.has_fresh_rebuild {
            return;
        }
        // Appears to be dual unbounded in phase 2 after fresh
        // rebuild. Normally this implies primal infeasibility, but only
        // allow this to be claimed if the proof of primal infeasibility
        // is true.
        let proof_of_infeasibility = self.proof_of_primal_infeasibility(ekk);
        if proof_of_infeasibility {
            // There is a proof of primal infeasibility
            self.solve_phase = kSolvePhaseExit;
            // Save dual ray information
            self.save_dual_ray(ekk);
            // Model status should be unset?
            debug_assert!(ekk.model_status_ == HighsModelStatus::kNotset);
            ekk.model_status_ = HighsModelStatus::kInfeasible;
        } else {
            // No proof of primal infeasibility, so assume dual unbounded
            // claim is spurious. Make row_out taboo, and prevent rebuild
            ekk.add_bad_basis_change(
                self.row_out,
                self.variable_out,
                self.variable_in,
                BadBasisChangeReason::kFailedInfeasibilityProof,
                true,
            );
            self.rebuild_reason = kRebuildReasonNo;
        }
    }
}
