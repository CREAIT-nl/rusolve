//! HiGHS source: simplex/HEkkPrimal.cpp (lines 1552-2138)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Hyper-sparse CHUZC, dual and primal updates (phase 1 and 2), primal correction.
//! See `HEkkPrimal` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkPrimal {
    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::hyperChooseColumn (lines 1552-1623, 72 lines)
    pub fn hyper_choose_column(&mut self, ekk: &mut HEkk) {
        if !self.use_hyper_chuzc {
            return;
        }
        if self.initialise_hyper_chuzc {
            return;
        }
        ekk.analysis_.simplex_timer_start(ChuzcHyperClock, 0);
        if self.report_hyper_chuzc {
            print!(
                "H-S  CHUZC: Max changed measure is {:9.4e} for column {:4}",
                self.max_changed_measure_value, self.max_changed_measure_column
            );
        }
        let mut best_measure = self.max_changed_measure_value;
        self.variable_in = -1;
        if self.max_changed_measure_column >= 0 {
            // Use max_changed_measure_column if it is well defined and has
            // nonzero dual. It may have been zeroed because it is taboo
            if ekk.info_.work_dual_[self.max_changed_measure_column as usize] != 0.0 {
                self.variable_in = self.max_changed_measure_column;
            }
        }
        let consider_nonbasic_free_column = *self.nonbasic_free_col_set.count() != 0;
        if self.num_hyper_chuzc_candidates != 0 {
            for i_entry in 1..=self.num_hyper_chuzc_candidates {
                let i_col = self.hyper_chuzc_candidate[i_entry as usize];
                if ekk.basis_.nonbasic_flag_[i_col as usize] == kNonbasicFlagFalse {
                    debug_assert!(ekk.basis_.nonbasic_move_[i_col as usize] == 0);
                    continue;
                }
                // Assess any dual infeasibility
                let mut dual_infeasibility = (-(ekk.basis_.nonbasic_move_[i_col as usize]
                    as HighsInt)) as f64
                    * ekk.info_.work_dual_[i_col as usize];
                if consider_nonbasic_free_column {
                    if self.nonbasic_free_col_set.r#in(i_col) {
                        dual_infeasibility = ekk.info_.work_dual_[i_col as usize].abs();
                    }
                }
                if dual_infeasibility > self.dual_feasibility_tolerance {
                    if dual_infeasibility * dual_infeasibility
                        > best_measure * self.edge_weight_[i_col as usize]
                    {
                        best_measure = dual_infeasibility * dual_infeasibility
                            / self.edge_weight_[i_col as usize];
                        self.variable_in = i_col;
                    }
                }
            }
        }
        if self.variable_in != self.max_changed_measure_column {
            if self.report_hyper_chuzc {
                print!(
                    ", and after HS CHUZC set it is now {:9.4e} for column {:4}",
                    best_measure, self.variable_in
                );
            }
            self.max_hyper_chuzc_non_candidate_measure = cpp_max(
                self.max_changed_measure_value,
                self.max_hyper_chuzc_non_candidate_measure,
            );
        }
        if best_measure >= self.max_hyper_chuzc_non_candidate_measure {
            // Candidate is at least as good as any unknown column, so accept it
            self.done_next_chuzc = true;
            if self.report_hyper_chuzc {
                println!(
                    ", and no       has  measure >  {:9.4e}",
                    self.max_hyper_chuzc_non_candidate_measure
                );
            }
        } else {
            // Candidate isn't as good as best unknown column, so do a full CHUZC
            // Shouldn't claim to have done the next CHUZC
            debug_assert!(!self.done_next_chuzc);
            self.done_next_chuzc = false;
            self.initialise_hyper_chuzc = true;
            if self.report_hyper_chuzc {
                println!(
                    ", but some may have measure >= {:9.4e}",
                    self.max_hyper_chuzc_non_candidate_measure
                );
            }
        }
        ekk.analysis_.simplex_timer_stop(ChuzcHyperClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::hyperChooseColumnStart (lines 1625-1629, 5 lines)
    pub fn hyper_choose_column_start(&mut self, ekk: &mut HEkk) {
        self.max_changed_measure_value = 0.0;
        self.max_changed_measure_column = -1;
        self.done_next_chuzc = false;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::hyperChooseColumnClear (lines 1631-1635, 5 lines)
    pub fn hyper_choose_column_clear(&mut self, ekk: &mut HEkk) {
        self.initialise_hyper_chuzc = self.use_hyper_chuzc;
        self.max_hyper_chuzc_non_candidate_measure = -1.0;
        self.done_next_chuzc = false;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::hyperChooseColumnChangedInfeasibility (lines 1637-1651, 15 lines)
    pub fn hyper_choose_column_changed_infeasibility(
        &mut self,
        ekk: &mut HEkk,
        infeasibility: f64,
        i_col: HighsInt,
    ) {
        if infeasibility * infeasibility
            > self.max_changed_measure_value * self.edge_weight_[i_col as usize]
        {
            self.max_hyper_chuzc_non_candidate_measure = cpp_max(
                self.max_changed_measure_value,
                self.max_hyper_chuzc_non_candidate_measure,
            );
            self.max_changed_measure_value =
                infeasibility * infeasibility / self.edge_weight_[i_col as usize];
            self.max_changed_measure_column = i_col;
        } else if infeasibility * infeasibility
            > self.max_hyper_chuzc_non_candidate_measure * self.edge_weight_[i_col as usize]
        {
            self.max_hyper_chuzc_non_candidate_measure =
                infeasibility * infeasibility / self.edge_weight_[i_col as usize];
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::hyperChooseColumnBasicFeasibilityChange (lines 1653-1693, 41 lines)
    pub fn hyper_choose_column_basic_feasibility_change(&mut self, ekk: &mut HEkk) {
        if !self.use_hyper_chuzc {
            return;
        }
        ekk.analysis_
            .simplex_timer_start(ChuzcHyperBasicFeasibilityChangeClock, 0);
        let mut to_entry: HighsInt = 0;
        let use_row_indices = ekk.simplex_nla_.sparse_loop_style(
            self.row_basic_feasibility_change.count,
            self.num_col,
            &mut to_entry,
        );
        for i_entry in 0..to_entry {
            let i_col = if use_row_indices {
                self.row_basic_feasibility_change.index[i_entry as usize]
            } else {
                i_entry
            };
            let dual_infeasibility = (-(ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt))
                as f64
                * ekk.info_.work_dual_[i_col as usize];
            if dual_infeasibility > self.dual_feasibility_tolerance {
                self.hyper_choose_column_changed_infeasibility(ekk, dual_infeasibility, i_col);
            }
        }
        let use_col_indices = ekk.simplex_nla_.sparse_loop_style(
            self.col_basic_feasibility_change.count,
            self.num_row,
            &mut to_entry,
        );
        for i_entry in 0..to_entry {
            let i_row = if use_col_indices {
                self.col_basic_feasibility_change.index[i_entry as usize]
            } else {
                i_entry
            };
            let i_col = self.num_col + i_row;
            let dual_infeasibility = (-(ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt))
                as f64
                * ekk.info_.work_dual_[i_col as usize];
            if dual_infeasibility > self.dual_feasibility_tolerance {
                self.hyper_choose_column_changed_infeasibility(ekk, dual_infeasibility, i_col);
            }
        }
        // Any nonbasic free columns will be handled explicitly in
        // hyperChooseColumnDualChange, so only look at them here if not
        // flipping
        let num_nonbasic_free_col = *self.nonbasic_free_col_set.count();
        if self.row_out < 0 && num_nonbasic_free_col != 0 {
            for i_entry in 0..num_nonbasic_free_col {
                let i_col = self.nonbasic_free_col_set.entry()[i_entry as usize];
                let dual_infeasibility = ekk.info_.work_dual_[i_col as usize].abs();
                if dual_infeasibility > self.dual_feasibility_tolerance {
                    self.hyper_choose_column_changed_infeasibility(ekk, dual_infeasibility, i_col);
                }
            }
        }
        ekk.analysis_
            .simplex_timer_stop(ChuzcHyperBasicFeasibilityChangeClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::hyperChooseColumnDualChange (lines 1695-1759, 65 lines)
    pub fn hyper_choose_column_dual_change(&mut self, ekk: &mut HEkk) {
        if !self.use_hyper_chuzc {
            return;
        }
        ekk.analysis_.simplex_timer_start(ChuzcHyperDualClock, 0);
        let mut to_entry: HighsInt = 0;
        // Look at changes in the columns and assess any dual infeasibility
        let use_row_indices =
            ekk.simplex_nla_
                .sparse_loop_style(self.row_ap.count, self.num_col, &mut to_entry);
        for i_entry in 0..to_entry {
            let i_col = if use_row_indices {
                self.row_ap.index[i_entry as usize]
            } else {
                i_entry
            };
            let dual_infeasibility = (-(ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt))
                as f64
                * ekk.info_.work_dual_[i_col as usize];
            if i_col == self.check_column && ekk.iteration_count_ >= self.check_iter {
                let measure =
                    dual_infeasibility * dual_infeasibility / self.edge_weight_[i_col as usize];
                if self.report_hyper_chuzc {
                    println!(
                        "Changing column {}: measure = {} ",
                        self.check_column, measure
                    );
                }
            }
            if dual_infeasibility > self.dual_feasibility_tolerance {
                self.hyper_choose_column_changed_infeasibility(ekk, dual_infeasibility, i_col);
            }
        }
        // Look at changes in the rows and assess any dual infeasibility
        let use_col_indices =
            ekk.simplex_nla_
                .sparse_loop_style(self.row_ep.count, self.num_row, &mut to_entry);
        for i_entry in 0..to_entry {
            let i_row = if use_col_indices {
                self.row_ep.index[i_entry as usize]
            } else {
                i_entry
            };
            let i_col = i_row + self.num_col;
            let dual_infeasibility = (-(ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt))
                as f64
                * ekk.info_.work_dual_[i_col as usize];
            if i_col == self.check_column && ekk.iteration_count_ >= self.check_iter {
                let measure =
                    dual_infeasibility * dual_infeasibility / self.edge_weight_[i_col as usize];
                if self.report_hyper_chuzc {
                    println!(
                        "Changing column {}: measure = {} ",
                        self.check_column, measure
                    );
                }
            }
            if dual_infeasibility > self.dual_feasibility_tolerance {
                self.hyper_choose_column_changed_infeasibility(ekk, dual_infeasibility, i_col);
            }
        }
        // Look for measure changes in any nonbasic free columns and assess
        // any dual infeasibility
        let num_nonbasic_free_col = *self.nonbasic_free_col_set.count();
        if num_nonbasic_free_col != 0 {
            for i_entry in 0..num_nonbasic_free_col {
                let i_col = self.nonbasic_free_col_set.entry()[i_entry as usize];
                let dual_infeasibility = ekk.info_.work_dual_[i_col as usize].abs();
                if dual_infeasibility > self.dual_feasibility_tolerance {
                    self.hyper_choose_column_changed_infeasibility(ekk, dual_infeasibility, i_col);
                }
            }
        }
        // Assess any dual infeasibility for the leaving column - should be dual
        // feasible!
        let i_col = self.variable_out;
        let dual_infeasibility = (-(ekk.basis_.nonbasic_move_[i_col as usize] as HighsInt)) as f64
            * ekk.info_.work_dual_[i_col as usize];
        if dual_infeasibility > self.dual_feasibility_tolerance {
            println!(
                "Dual infeasibility {} for leaving column!",
                dual_infeasibility
            );
            debug_assert!(dual_infeasibility <= self.dual_feasibility_tolerance);
            self.hyper_choose_column_changed_infeasibility(ekk, dual_infeasibility, i_col);
        }
        ekk.analysis_.simplex_timer_stop(ChuzcHyperDualClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::updateDual (lines 1761-1787, 27 lines)
    pub fn update_dual(&mut self, ekk: &mut HEkk) {
        ekk.analysis_.simplex_timer_start(UpdateDualClock, 0);
        debug_assert!(self.alpha_col != 0.0);
        debug_assert!(self.row_out >= 0);
        // Update the duals
        self.theta_dual = ekk.info_.work_dual_[self.variable_in as usize] / self.alpha_col;
        for i_el in 0..self.row_ap.count {
            let i_col = self.row_ap.index[i_el as usize];
            ekk.info_.work_dual_[i_col as usize] -=
                self.theta_dual * self.row_ap.array[i_col as usize];
        }
        for i_el in 0..self.row_ep.count {
            let i_row = self.row_ep.index[i_el as usize];
            let i_col = i_row + self.num_col;
            ekk.info_.work_dual_[i_col as usize] -=
                self.theta_dual * self.row_ep.array[i_row as usize];
        }
        // Dual for the pivot
        ekk.info_.work_dual_[self.variable_in as usize] = 0.0;
        ekk.info_.work_dual_[self.variable_out as usize] = -self.theta_dual;

        ekk.invalidate_dual_infeasibility_record();
        // After dual update in primal simplex the dual objective value is not known
        ekk.status_.has_dual_objective_value = false;
        ekk.analysis_.simplex_timer_stop(UpdateDualClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::phase1ComputeDual (lines 1789-1842, 54 lines)
    /// MUST-PORT (determinism, random stream): reads `num_tot_random_value_`.
    pub fn phase1_compute_dual(&mut self, ekk: &mut HEkk) {
        let mut buffer = HVector::default();
        buffer.setup(self.num_row);
        buffer.clear();
        buffer.count = 0;
        // Accumulate costs for checking
        ekk.info_.work_cost_ = vec![0.0; self.num_tot as usize];
        // Zero the dual values
        ekk.info_.work_dual_ = vec![0.0; self.num_tot as usize];
        // Determine the base value for cost perturbation
        let base = ekk.info_.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
        for i_row in 0..self.num_row {
            let value = ekk.info_.base_value_[i_row as usize];
            let lower = ekk.info_.base_lower_[i_row as usize];
            let upper = ekk.info_.base_upper_[i_row as usize];
            let mut bound_violated: HighsInt = 0;
            if value < lower - self.primal_feasibility_tolerance {
                bound_violated = -1;
            } else if value > upper + self.primal_feasibility_tolerance {
                bound_violated = 1;
            }
            if bound_violated == 0 {
                continue;
            }
            let mut cost = bound_violated as f64;
            if base != 0.0 {
                cost *= 1.0 + base * ekk.info_.num_tot_random_value_[i_row as usize];
            }
            buffer.array[i_row as usize] = cost;
            buffer.index[buffer.count as usize] = i_row;
            buffer.count += 1;
        }
        if buffer.count <= 0 {
            // Strange, should be a non-trivial RHS
            debug_assert!(buffer.count > 0);
            return;
        }
        for i_row in 0..self.num_row {
            let basic = ekk.basis_.basic_index_[i_row as usize];
            ekk.info_.work_cost_[basic as usize] = buffer.array[i_row as usize];
        }
        //
        // Full BTRAN
        //
        ekk.full_btran(&mut buffer);
        //
        // Full PRICE
        //
        let mut buffer_long = HVector::default();
        buffer_long.setup(self.num_col);
        ekk.full_price(&buffer, &mut buffer_long);

        for i_col in 0..self.num_col {
            ekk.info_.work_dual_[i_col as usize] = (-(ekk.basis_.nonbasic_flag_[i_col as usize]
                as HighsInt)) as f64
                * buffer_long.array[i_col as usize];
        }
        let mut i_row = 0;
        let mut i_col = self.num_col;
        while i_row < self.num_row {
            ekk.info_.work_dual_[i_col as usize] = (-(ekk.basis_.nonbasic_flag_[i_col as usize]
                as HighsInt)) as f64
                * buffer.array[i_row as usize];
            i_row += 1;
            i_col += 1;
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::phase1UpdatePrimal (lines 1844-1891, 48 lines)
    /// MUST-PORT (determinism, random stream): reads `num_tot_random_value_`.
    pub fn phase1_update_primal(&mut self, ekk: &mut HEkk) {
        ekk.analysis_.simplex_timer_start(UpdatePrimalClock, 0);
        self.col_basic_feasibility_change.clear();
        //
        // Update basic primal values, identifying all the feasibility
        // changes giving a value to col_basic_feasibility_change so that the duals
        // can be updated.
        //
        // Determine the base value for cost perturbation
        let base = ekk.info_.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
        for i_el in 0..self.col_aq.count {
            let i_row = self.col_aq.index[i_el as usize];
            ekk.info_.base_value_[i_row as usize] -=
                self.theta_primal * self.col_aq.array[i_row as usize];
            let i_col = ekk.basis_.basic_index_[i_row as usize];
            let was_cost = ekk.info_.work_cost_[i_col as usize];
            let value = ekk.info_.base_value_[i_row as usize];
            let lower = ekk.info_.base_lower_[i_row as usize];
            let upper = ekk.info_.base_upper_[i_row as usize];
            let mut bound_violated: HighsInt = 0;
            if value < lower - self.primal_feasibility_tolerance {
                bound_violated = -1;
            } else if value > upper + self.primal_feasibility_tolerance {
                bound_violated = 1;
            }
            let mut cost = bound_violated as f64;
            if base != 0.0 {
                cost *= 1.0 + base * ekk.info_.num_tot_random_value_[i_row as usize];
            }
            ekk.info_.work_cost_[i_col as usize] = cost;
            if was_cost != 0.0 {
                if cost == 0.0 {
                    ekk.info_.num_primal_infeasibilities -= 1;
                }
            } else {
                if cost != 0.0 {
                    ekk.info_.num_primal_infeasibilities += 1;
                }
            }
            let delta_cost = cost - was_cost;
            if delta_cost != 0.0 {
                self.col_basic_feasibility_change.array[i_row as usize] = delta_cost;
                let cnt = self.col_basic_feasibility_change.count;
                self.col_basic_feasibility_change.index[cnt as usize] = i_row;
                self.col_basic_feasibility_change.count += 1;
                if i_col >= self.num_col {
                    ekk.info_.work_dual_[i_col as usize] += delta_cost;
                }
            }
        }
        // Don't set baseValue[row_out] yet so that dual update due to
        // feasibility changes is done correctly
        ekk.invalidate_primal_max_sum_infeasibility_record();
        ekk.analysis_.simplex_timer_stop(UpdatePrimalClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::considerInfeasibleValueIn (lines 1893-1951, 59 lines)
    /// MUST-PORT (determinism, random stream): reads `num_tot_random_value_`.
    pub fn consider_infeasible_value_in(&mut self, ekk: &mut HEkk) {
        debug_assert!(self.row_out >= 0);
        // Determine the base value for cost perturbation
        let base = ekk.info_.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
        let lower = ekk.info_.work_lower_[self.variable_in as usize];
        let upper = ekk.info_.work_upper_[self.variable_in as usize];
        let mut bound_violated: HighsInt = 0;
        if self.value_in < lower - self.primal_feasibility_tolerance {
            bound_violated = -1;
        } else if self.value_in > upper + self.primal_feasibility_tolerance {
            bound_violated = 1;
        }
        if bound_violated == 0 {
            return;
        }
        // The primal value of the entering variable is not feasible
        if self.solve_phase == kSolvePhase1 {
            ekk.info_.num_primal_infeasibilities += 1;
            let mut cost = bound_violated as f64;
            if base != 0.0 {
                cost *= 1.0 + base * ekk.info_.num_tot_random_value_[self.row_out as usize];
            }
            ekk.info_.work_cost_[self.variable_in as usize] = cost;
            ekk.info_.work_dual_[self.variable_in as usize] += cost;
        } else if self.primal_correction_strategy == kSimplexPrimalCorrectionStrategyNone {
            // @primal_infeasibility calculation
            let primal_infeasibility = if bound_violated < 0 {
                lower - self.value_in
            } else {
                self.value_in - upper
            };
            ekk.info_.num_primal_infeasibilities += 1;
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kWarning,
                "Entering variable has primal infeasibility of {} for [{}, {}, {}]\n",
                primal_infeasibility,
                lower,
                self.value_in,
                upper
            );
            self.rebuild_reason = kRebuildReasonPrimalInfeasibleInPrimalSimplex;
        } else {
            let mut bound_shift: f64 = 0.0;
            let vi = self.variable_in as usize;
            if bound_violated > 0 {
                // Perturb the upper bound to accommodate the infeasibility
                let mut bound = ekk.info_.work_upper_[vi];
                let random_value = ekk.info_.num_tot_random_value_[vi];
                self.shift_bound(
                    ekk,
                    false,
                    self.variable_in,
                    self.value_in,
                    random_value,
                    &mut bound,
                    &mut bound_shift,
                );
                ekk.info_.work_upper_[vi] = bound;
                ekk.info_.work_upper_shift_[vi] += bound_shift;
            } else {
                // Perturb the lower bound to accommodate the infeasibility
                let mut bound = ekk.info_.work_lower_[vi];
                let random_value = ekk.info_.num_tot_random_value_[vi];
                self.shift_bound(
                    ekk,
                    true,
                    self.variable_in,
                    self.value_in,
                    random_value,
                    &mut bound,
                    &mut bound_shift,
                );
                ekk.info_.work_lower_[vi] = bound;
                ekk.info_.work_lower_shift_[vi] += bound_shift;
            }
            // Surely better to record this
            //
            //    info.bounds_shifted = true;
            ekk.info_.bounds_perturbed = true;
        }
        ekk.invalidate_primal_max_sum_infeasibility_record();
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::phase2UpdatePrimal (lines 1953-2060, 108 lines)
    /// MUST-PORT (determinism, random stream): reads `num_tot_random_value_`.
    /// C++ default arguments: initialise = false
    pub fn phase2_update_primal(&mut self, ekk: &mut HEkk, initialise: bool) {
        if initialise {
            self.max_max_local_primal_infeasibility_ = 0.0;
            self.max_max_ignored_violation_ = 0.0;
            return;
        }
        ekk.analysis_.simplex_timer_start(UpdatePrimalClock, 0);
        let mut primal_infeasible = false;
        let mut max_local_primal_infeasibility: f64 = 0.0;
        let mut max_ignored_violation: f64 = 0.0;
        // If shifts are only identified in rebuild() the bounds can be
        // ignored. If they aren't ignored, then violations lead to either
        // identification of infeasibilities (and return to Phase 1) or
        // shifting of bounds to accommodate them.
        //
        // primal_correction_strategy is defined (const) as
        // kSimplexPrimalCorrectionStrategyAlways
        let ignore_bounds =
            self.primal_correction_strategy == kSimplexPrimalCorrectionStrategyInRebuild;
        debug_assert!(self.primal_correction_strategy == kSimplexPrimalCorrectionStrategyAlways);
        let mut to_entry: HighsInt = 0;
        let use_col_indices =
            ekk.simplex_nla_
                .sparse_loop_style(self.col_aq.count, self.num_row, &mut to_entry);
        for i_entry in 0..to_entry {
            let i_row = if use_col_indices {
                self.col_aq.index[i_entry as usize]
            } else {
                i_entry
            };
            ekk.info_.base_value_[i_row as usize] -=
                self.theta_primal * self.col_aq.array[i_row as usize];
            // Determine whether a bound is violated and take action
            let lower = ekk.info_.base_lower_[i_row as usize];
            let upper = ekk.info_.base_upper_[i_row as usize];
            let value = ekk.info_.base_value_[i_row as usize];
            let mut bound_violated: HighsInt = 0;
            if value < lower - self.primal_feasibility_tolerance {
                bound_violated = -1;
            } else if value > upper + self.primal_feasibility_tolerance {
                bound_violated = 1;
            }
            if bound_violated == 0 {
                continue;
            }
            // A bound is violated
            if self.primal_correction_strategy == kSimplexPrimalCorrectionStrategyNone {
                // Not used
                if cfg!(debug_assertions) {
                    panic!("assert(111 == 222)");
                }
                // @primal_infeasibility calculation
                let primal_infeasibility = if bound_violated < 0 {
                    lower - value
                } else {
                    value - upper
                };
                max_local_primal_infeasibility =
                    cpp_max(primal_infeasibility, max_local_primal_infeasibility);
                if primal_infeasibility > self.primal_feasibility_tolerance {
                    ekk.info_.num_primal_infeasibilities += 1;
                    primal_infeasible = true;
                }
            } else if ignore_bounds {
                // Not used
                if cfg!(debug_assertions) {
                    panic!("assert(111 == 333)");
                }
                let ignored_violation = if bound_violated < 0 {
                    lower - value
                } else {
                    value - upper
                };
                max_ignored_violation = cpp_max(ignored_violation, max_ignored_violation);
            } else {
                let i_col = ekk.basis_.basic_index_[i_row as usize];
                let ic = i_col as usize;
                let mut bound_shift: f64 = 0.0;
                if bound_violated > 0 {
                    // Perturb the upper bound to accommodate the infeasibility
                    let mut bound = ekk.info_.work_upper_[ic];
                    let base_value = ekk.info_.base_value_[i_row as usize];
                    let random_value = ekk.info_.num_tot_random_value_[ic];
                    self.shift_bound(
                        ekk,
                        false,
                        i_col,
                        base_value,
                        random_value,
                        &mut bound,
                        &mut bound_shift,
                    );
                    ekk.info_.work_upper_[ic] = bound;
                    ekk.info_.base_upper_[i_row as usize] = ekk.info_.work_upper_[ic];
                    ekk.info_.work_upper_shift_[ic] += bound_shift;
                } else {
                    // Perturb the lower bound to accommodate the infeasibility
                    let mut bound = ekk.info_.work_lower_[ic];
                    let base_value = ekk.info_.base_value_[i_row as usize];
                    let random_value = ekk.info_.num_tot_random_value_[ic];
                    self.shift_bound(
                        ekk,
                        true,
                        i_col,
                        base_value,
                        random_value,
                        &mut bound,
                        &mut bound_shift,
                    );
                    ekk.info_.work_lower_[ic] = bound;
                    ekk.info_.base_lower_[i_row as usize] = ekk.info_.work_lower_[ic];
                    ekk.info_.work_lower_shift_[ic] += bound_shift;
                }
                ekk.info_.bounds_shifted = true;
                debug_assert!(bound_shift > 0.0);
            }
        }
        if primal_infeasible {
            self.rebuild_reason = kRebuildReasonPrimalInfeasibleInPrimalSimplex;
            if max_local_primal_infeasibility > self.max_max_local_primal_infeasibility_ * 2.0 {
                self.max_max_local_primal_infeasibility_ = max_local_primal_infeasibility;
                println!(
                    "phase2UpdatePrimal: max_local_primal_infeasibility = {}",
                    max_local_primal_infeasibility
                );
            }
            ekk.invalidate_primal_max_sum_infeasibility_record();
        }
        if max_ignored_violation > self.max_max_ignored_violation_ * 2.0 {
            self.max_max_ignored_violation_ = max_ignored_violation;
            println!(
                "phase2UpdatePrimal: max_ignored_violation = {}",
                max_ignored_violation
            );
        }
        ekk.info_.updated_primal_objective_value +=
            ekk.info_.work_dual_[self.variable_in as usize] * self.theta_primal;

        ekk.analysis_.simplex_timer_stop(UpdatePrimalClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::correctPrimal (lines 2062-2138, 77 lines)
    /// MUST-PORT (determinism, random stream): reads `num_tot_random_value_`.
    /// C++ default arguments: initialise = false
    pub fn correct_primal(&mut self, ekk: &mut HEkk, initialise: bool) -> bool {
        if self.primal_correction_strategy == kSimplexPrimalCorrectionStrategyNone {
            return true;
        }
        if initialise {
            self.max_max_primal_correction_ = 0.0;
            return true;
        }
        debug_assert!(self.solve_phase == kSolvePhase2);
        let mut num_primal_correction: HighsInt = 0;
        let mut max_primal_correction: f64 = 0.0;
        let mut sum_primal_correction: f64 = 0.0;
        let mut num_primal_correction_skipped: HighsInt = 0;
        for i_row in 0..self.num_row {
            let lower = ekk.info_.base_lower_[i_row as usize];
            let upper = ekk.info_.base_upper_[i_row as usize];
            let value = ekk.info_.base_value_[i_row as usize];
            let mut bound_violated: HighsInt = 0;
            if value < lower - self.primal_feasibility_tolerance {
                bound_violated = -1;
            } else if value > upper + self.primal_feasibility_tolerance {
                bound_violated = 1;
            }
            if bound_violated != 0 {
                if ekk.info_.allow_bound_perturbation {
                    let i_col = ekk.basis_.basic_index_[i_row as usize];
                    let ic = i_col as usize;
                    let mut bound_shift: f64 = 0.0;
                    if bound_violated > 0 {
                        // Perturb the upper bound to accommodate the infeasibility
                        let mut bound = ekk.info_.work_upper_[ic];
                        let base_value = ekk.info_.base_value_[i_row as usize];
                        let random_value = ekk.info_.num_tot_random_value_[ic];
                        self.shift_bound(
                            ekk,
                            false,
                            i_col,
                            base_value,
                            random_value,
                            &mut bound,
                            &mut bound_shift,
                        );
                        ekk.info_.work_upper_[ic] = bound;
                        ekk.info_.base_upper_[i_row as usize] = ekk.info_.work_upper_[ic];
                        ekk.info_.work_upper_shift_[ic] += bound_shift;
                    } else {
                        // Perturb the lower bound to accommodate the infeasibility
                        let mut bound = ekk.info_.work_lower_[ic];
                        let base_value = ekk.info_.base_value_[i_row as usize];
                        let random_value = ekk.info_.num_tot_random_value_[ic];
                        self.shift_bound(
                            ekk,
                            true,
                            i_col,
                            base_value,
                            random_value,
                            &mut bound,
                            &mut bound_shift,
                        );
                        ekk.info_.work_lower_[ic] = bound;
                        ekk.info_.base_lower_[i_row as usize] = ekk.info_.work_lower_[ic];
                        ekk.info_.work_lower_shift_[ic] += bound_shift;
                    }
                    debug_assert!(bound_shift > 0.0);
                    num_primal_correction += 1;
                    max_primal_correction = cpp_max(bound_shift, max_primal_correction);
                    sum_primal_correction += bound_shift;
                    // Surely better to record this
                    //
                    // info.bounds_shifted = true;
                    ekk.info_.bounds_perturbed = true;
                } else {
                    // Bound perturbation is not permitted
                    num_primal_correction_skipped += 1;
                }
            }
        }
        if num_primal_correction_skipped != 0 {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kError,
                "correctPrimal: Missed {} bound shifts\n",
                num_primal_correction_skipped
            );
            if kAllowDeveloperAssert {
                debug_assert!(num_primal_correction_skipped == 0);
            }
            return false;
        }
        if max_primal_correction > 2.0 * self.max_max_primal_correction_ {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "phase2CorrectPrimal: num / max / sum primal corrections = {} / {} / {}\n",
                num_primal_correction,
                max_primal_correction,
                sum_primal_correction
            );
            self.max_max_primal_correction_ = max_primal_correction;
        }
        true
    }
}
