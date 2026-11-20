//! HiGHS source: simplex/HEkkDual.cpp (lines 1942-2339)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! FTRANs (column, BFRT, DSE), dual and primal updates, cost shifting, pivots, Devex,
//! edge-weight strategy.
//! See `HEkkDual` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkDual {
    /// HiGHS: HEkkDual.cpp::HEkkDual::updateFtran (lines 1942-1968, 27 lines)
    pub fn update_ftran(&mut self, ekk: &mut HEkk) {
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }
        ekk.analysis_.simplex_timer_start(FtranClock, 0);
        // Clear the pivotal column and indicate that its values should be packed
        self.col_aq.clear();
        self.col_aq.pack_flag = true;
        // Get the constraint matrix column by combining just one column
        // with unit multiplier
        ekk.lp_
            .a_matrix_
            .collect_aj(&mut self.col_aq, self.variable_in, 1.0);
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.operation_record_before(
                kSimplexNlaFtran,
                &self.col_aq,
                ekk.info_.col_aq_density,
            );
        }
        // Perform FTRAN
        {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            ekk.simplex_nla_
                .ftran(&refs, &mut self.col_aq, ekk.info_.col_aq_density, None);
        }
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_
                .operation_record_after(kSimplexNlaFtran, &self.col_aq);
        }
        let local_col_aq_density = (self.col_aq.count as f64) * self.inv_solver_num_row;
        let mut col_aq_density = ekk.info_.col_aq_density;
        ekk.update_operation_result_density(local_col_aq_density, &mut col_aq_density);
        ekk.info_.col_aq_density = col_aq_density;
        // Save the pivot value computed column-wise - used for numerical checking
        self.alpha_col = self.col_aq.array[self.row_out as usize];
        ekk.analysis_.simplex_timer_stop(FtranClock, 0);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::updateFtranBFRT (lines 1970-2007, 38 lines)
    pub fn update_ftran_bfrt(&mut self, ekk: &mut HEkk) {
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }

        // Only time updateFtranBFRT if dualRow.workCount > 0;
        let time_update_ftran_bfrt = self.dual_row.work_count > 0;

        if time_update_ftran_bfrt {
            ekk.analysis_.simplex_timer_start(FtranBfrtClock, 0);
        }

        self.dual_row.update_flip(ekk, &mut self.col_bfrt);

        if self.col_bfrt.count != 0 {
            if ekk.analysis_.analyse_simplex_summary_data {
                ekk.analysis_.operation_record_before(
                    kSimplexNlaFtranBfrt,
                    &self.col_bfrt,
                    ekk.info_.col_bfrt_density,
                );
            }
            {
                let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
                ekk.simplex_nla_
                    .ftran(&refs, &mut self.col_bfrt, ekk.info_.col_bfrt_density, None);
            }
            if ekk.analysis_.analyse_simplex_summary_data {
                ekk.analysis_
                    .operation_record_after(kSimplexNlaFtranBfrt, &self.col_bfrt);
            }
        }
        if time_update_ftran_bfrt {
            ekk.analysis_.simplex_timer_stop(FtranBfrtClock, 0);
        }
        let local_col_bfrt_density = (self.col_bfrt.count as f64) * self.inv_solver_num_row;
        let mut col_bfrt_density = ekk.info_.col_bfrt_density;
        ekk.update_operation_result_density(local_col_bfrt_density, &mut col_bfrt_density);
        ekk.info_.col_bfrt_density = col_bfrt_density;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::updateFtranDSE (lines 2009-2048, 40 lines)
    pub fn update_ftran_dse(&mut self, ekk: &mut HEkk, dse_vector: &mut HVector) {
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }
        ekk.analysis_.simplex_timer_start(FtranDseClock, 0);
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.operation_record_before(
                kSimplexNlaFtranDse,
                dse_vector,
                ekk.info_.row_dse_density,
            );
        }
        {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            // Apply R{-1}
            ekk.simplex_nla_
                .unapply_basis_matrix_row_scale(&refs, dse_vector);

            // Perform FTRAN DSE
            ekk.simplex_nla_
                .ftran_in_scaled_space(dse_vector, ekk.info_.row_dse_density, None);
        }
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_
                .operation_record_after(kSimplexNlaFtranDse, dse_vector);
        }
        ekk.analysis_.simplex_timer_stop(FtranDseClock, 0);
        let local_row_dse_density = (dse_vector.count as f64) * self.inv_solver_num_row;
        let mut row_dse_density = ekk.info_.row_dse_density;
        ekk.update_operation_result_density(local_row_dse_density, &mut row_dse_density);
        ekk.info_.row_dse_density = row_dse_density;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::updateVerify (lines 2050-2063, 14 lines)
    pub fn update_verify(&mut self, ekk: &mut HEkk) {
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }

        // Use the two pivot values to identify numerical trouble
        if ekk.reinvert_on_numerical_trouble(
            "HEkkDual::updateVerify",
            &mut self.numerical_trouble,
            self.alpha_col,
            self.alpha_row,
            kNumericalTroubleTolerance,
        ) {
            self.rebuild_reason = kRebuildReasonPossiblySingularBasis;
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::updateDual (lines 2065-2127, 63 lines)
    pub fn update_dual(&mut self, ekk: &mut HEkk) {
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }

        // Update - dual (shift and back)
        if self.theta_dual == 0.0 {
            // Little to do if theta_dual is zero
            let amount = -ekk.info_.work_dual_[self.variable_in as usize];
            self.shift_cost(ekk, self.variable_in, amount);
        } else {
            // Update the whole vector of dual values
            self.dual_row.update_dual(ekk, self.theta_dual);
            if ekk.info_.simplex_strategy != kSimplexStrategyDualPlain && self.slice_price != 0 {
                // Update the slice-by-slice copy of dual variables
                for i in 0..self.slice_num {
                    self.slice_dual_row[i as usize].update_dual(ekk, self.theta_dual);
                }
            }
        }
        // Identify the changes in the dual objective
        let mut dual_objective_value_change: f64;
        let variable_in_delta_dual = ekk.info_.work_dual_[self.variable_in as usize];
        let variable_in_value = ekk.info_.work_value_[self.variable_in as usize];
        let variable_in_nonbasic_flag: HighsInt =
            ekk.basis_.nonbasic_flag_[self.variable_in as usize] as HighsInt;
        dual_objective_value_change =
            (variable_in_nonbasic_flag as f64) * (-variable_in_value * variable_in_delta_dual);
        dual_objective_value_change *= ekk.cost_scale_;
        ekk.info_.updated_dual_objective_value += dual_objective_value_change;
        // Surely variable_out_nonbasicFlag is always 0 since it's basic - so there's
        // no dual objective change
        let variable_out_nonbasic_flag: HighsInt =
            ekk.basis_.nonbasic_flag_[self.variable_out as usize] as HighsInt;
        debug_assert!(variable_out_nonbasic_flag == 0);
        if variable_out_nonbasic_flag != 0 {
            let variable_out_delta_dual =
                ekk.info_.work_dual_[self.variable_out as usize] - self.theta_dual;
            let variable_out_value = ekk.info_.work_value_[self.variable_out as usize];
            dual_objective_value_change = (variable_out_nonbasic_flag as f64)
                * (-variable_out_value * variable_out_delta_dual);
            dual_objective_value_change *= ekk.cost_scale_;
            ekk.info_.updated_dual_objective_value += dual_objective_value_change;
        }
        ekk.info_.work_dual_[self.variable_in as usize] = 0.0;
        ekk.info_.work_dual_[self.variable_out as usize] = -self.theta_dual;

        self.shift_back(ekk, self.variable_out);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::updatePrimal (lines 2129-2196, 68 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion): accumulates `total_synthetic_tick_` from the vectors' `synthetic_tick`.
    pub fn update_primal(&mut self, ekk: &mut HEkk, dse_vector: &mut HVector) {
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }
        if self.edge_weight_mode == EdgeWeightMode::kDevex {
            let updated_edge_weight = ekk.dual_edge_weight_[self.row_out as usize];
            ekk.dual_edge_weight_[self.row_out as usize] = self.computed_edge_weight;
            self.new_devex_framework = self.new_devex_framework(ekk, updated_edge_weight);
        }
        // DSE_Vector is either col_DSE = B^{-1}B^{-T}e_p (if using dual
        // steepest edge weights) or row_ep = B^{-T}e_p.
        //
        // Update - primal and weight
        self.dual_rhs.update_primal(ekk, &self.col_bfrt, 1.0);
        self.dual_rhs.update_infeas_list(ekk, &self.col_bfrt);
        let x_out = ekk.info_.base_value_[self.row_out as usize];
        let l_out = ekk.info_.base_lower_[self.row_out as usize];
        let u_out = ekk.info_.base_upper_[self.row_out as usize];
        self.theta_primal = (x_out
            - (if self.delta_primal < 0.0 {
                l_out
            } else {
                u_out
            }))
            / self.alpha_col;
        let ok_update_primal = self
            .dual_rhs
            .update_primal(ekk, &self.col_aq, self.theta_primal);
        if !ok_update_primal {
            self.rebuild_reason = kRebuildReasonExcessivePrimalValue;
            return;
        }
        ekk.update_bad_basis_change(&self.col_aq, self.theta_primal);
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            let pivot_in_scaled_space = {
                let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
                ekk.simplex_nla_.pivot_in_scaled_space(
                    &refs,
                    &self.col_aq,
                    self.variable_in,
                    self.row_out,
                )
            };
            if ekk.simplex_in_scaled_space_ {
                debug_assert!(pivot_in_scaled_space == self.alpha_col);
            }
            let new_pivotal_edge_weight = ekk.dual_edge_weight_[self.row_out as usize]
                / (pivot_in_scaled_space * pivot_in_scaled_space);
            let kai = -2.0 / pivot_in_scaled_space;
            ekk.update_dual_steepest_edge_weights(
                self.row_out,
                self.variable_in,
                &self.col_aq,
                new_pivotal_edge_weight,
                kai,
                &dse_vector.array,
            );
            ekk.dual_edge_weight_[self.row_out as usize] = new_pivotal_edge_weight;
        } else if self.edge_weight_mode == EdgeWeightMode::kDevex {
            // Pivotal row is for the current basis: weights are required for
            // the next basis so have to divide the current (exact) weight by
            // the pivotal value
            let mut new_pivotal_edge_weight =
                ekk.dual_edge_weight_[self.row_out as usize] / (self.alpha_col * self.alpha_col);
            new_pivotal_edge_weight = cpp_max(1.0, new_pivotal_edge_weight);
            // Update rest of weights
            ekk.update_dual_devex_weights(&self.col_aq, new_pivotal_edge_weight);
            ekk.dual_edge_weight_[self.row_out as usize] = new_pivotal_edge_weight;
            self.num_devex_iterations += 1;
        }
        self.dual_rhs.update_infeas_list(ekk, &self.col_aq);

        // Whether or not dual steepest edge weights are being used, have to
        // add in DSE_Vector->synthetic_tick_ since this contains the
        // contribution from forming row_ep = B^{-T}e_p.
        ekk.total_synthetic_tick_ += self.col_aq.synthetic_tick;
        ekk.total_synthetic_tick_ += dse_vector.synthetic_tick;
        // reporting-only work counter, never read by the search.
        crate::work::add_simplex_ticks(self.col_aq.synthetic_tick + dse_vector.synthetic_tick);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::shiftCost (lines 2199-2212, 14 lines)
    pub fn shift_cost(&mut self, ekk: &mut HEkk, i_col: HighsInt, amount: f64) {
        ekk.info_.costs_shifted = true;
        debug_assert!(ekk.info_.work_shift_[i_col as usize] == 0.0);
        if amount == 0.0 {
            return;
        }
        let use_amount = amount;
        ekk.info_.work_shift_[i_col as usize] = use_amount;
        // Analysis
        let shift = use_amount.abs();
        ekk.analysis_.net_num_single_cost_shift += 1;
        ekk.analysis_.num_single_cost_shift += 1;
        ekk.analysis_.sum_single_cost_shift += shift;
        ekk.analysis_.max_single_cost_shift = cpp_max(shift, ekk.analysis_.max_single_cost_shift);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::shiftBack (lines 2215-2222, 8 lines)
    pub fn shift_back(&mut self, ekk: &mut HEkk, i_col: HighsInt) {
        if ekk.info_.work_shift_[i_col as usize] == 0.0 {
            return;
        }
        ekk.info_.work_dual_[i_col as usize] -= ekk.info_.work_shift_[i_col as usize];
        ekk.info_.work_shift_[i_col as usize] = 0.0;
        // Analysis
        ekk.analysis_.net_num_single_cost_shift -= 1;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::updatePivots (lines 2224-2267, 44 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn update_pivots(&mut self, ekk: &mut HEkk) {
        // If reinversion is needed then skip this method
        if self.rebuild_reason != 0 {
            return;
        }
        // Transform the vectors used in updateFactor if the simplex NLA involves
        // scaling
        ekk.transform_for_update(
            &mut self.col_aq,
            &mut self.row_ep,
            self.variable_in,
            &mut self.row_out,
        );
        //
        // Update the sets of indices of basic and nonbasic variables
        //
        ekk.update_pivots(self.variable_in, self.row_out, self.move_out);
        //
        ekk.iteration_count_ += 1;
        //
        // Update the invertible representation of the basis matrix
        ekk.update_factor(
            &mut self.col_aq,
            &mut self.row_ep,
            &mut self.row_out,
            &mut self.rebuild_reason,
        );
        //
        // Update the row-wise representation of the nonbasic columns
        ekk.update_matrix(self.variable_in, self.variable_out);
        //
        // Delete Freelist entry for variable_in
        self.dual_row.delete_freelist(ekk, self.variable_in);
        //
        // Update the primal value for the row where the basis change has
        // occurred, and set the corresponding primal infeasibility value in
        // dualRHS.work_infeasibility
        let value = ekk.info_.work_value_[self.variable_in as usize] + self.theta_primal;
        self.dual_rhs.update_pivots(ekk, self.row_out, value);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::initialiseDevexFramework (lines 2269-2298, 30 lines)
    pub fn initialise_devex_framework(&mut self, ekk: &mut HEkk) {
        // Initialise the Devex framework: reference set is all basic
        // variables
        ekk.analysis_.simplex_timer_start(DevexIzClock, 0);
        // Resize in case this is the first call
        ekk.info_
            .devex_index_
            .resize(self.solver_num_tot as usize, 0);
        // The devex reference set is initialised to be the current set of
        // basic variables: devex_index = 1-nonbasicFlag^2
        for vr_n in 0..self.solver_num_tot {
            let flag = ekk.basis_.nonbasic_flag_[vr_n as usize] as HighsInt;
            ekk.info_.devex_index_[vr_n as usize] = 1 - flag * flag;
        }
        // Set all initial weights to 1, zero the count of iterations with
        // this Devex framework, increment the number of Devex frameworks
        // and indicate that there's no need for a new Devex framework
        ekk.dual_edge_weight_.clear();
        ekk.dual_edge_weight_
            .resize(self.solver_num_row as usize, 1.0);
        self.num_devex_iterations = 0;
        self.new_devex_framework = false;
        self.minor_new_devex_framework = false;
        ekk.analysis_.simplex_timer_stop(DevexIzClock, 0);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::interpretDualEdgeWeightStrategy (lines 2300-2324, 25 lines)
    pub fn interpret_dual_edge_weight_strategy(
        &mut self,
        ekk: &mut HEkk,
        dual_edge_weight_strategy: HighsInt,
    ) {
        if dual_edge_weight_strategy == kSimplexEdgeWeightStrategyChoose {
            self.edge_weight_mode = EdgeWeightMode::kSteepestEdge;
            self.allow_dual_steepest_edge_to_devex_switch = true;
        } else if dual_edge_weight_strategy == kSimplexEdgeWeightStrategyDantzig {
            self.edge_weight_mode = EdgeWeightMode::kDantzig;
        } else if dual_edge_weight_strategy == kSimplexEdgeWeightStrategyDevex {
            self.edge_weight_mode = EdgeWeightMode::kDevex;
        } else if dual_edge_weight_strategy == kSimplexEdgeWeightStrategySteepestEdge {
            self.edge_weight_mode = EdgeWeightMode::kSteepestEdge;
            self.allow_dual_steepest_edge_to_devex_switch = false;
        } else {
            if cfg!(debug_assertions) {
                panic!("assert(1 == 0)");
            }
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "HEkkDual::interpretDualEdgeWeightStrategy: unrecognised dual_edge_weight_strategy = {} - using dual steepest edge with possible switch to Devex\n",
                dual_edge_weight_strategy
            );
            self.edge_weight_mode = EdgeWeightMode::kSteepestEdge;
            self.allow_dual_steepest_edge_to_devex_switch = true;
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::possiblyUseLiDualSteepestEdge (lines 2326-2339, 14 lines)
    pub fn possibly_use_li_dual_steepest_edge(&mut self, ekk: &mut HEkk) {
        // Decide whether to use LiDSE by not storing squared primal infeasibilities
        ekk.info_.store_squared_primal_infeasibility = true;
        if ekk.options_.less_infeasible_DSE_check
            && is_less_infeasible_dse_candidate(&ekk.options_.log_options, &ekk.lp_)
            && ekk.options_.less_infeasible_DSE_choose_row
        {
            // LP is a candidate for LiDSE: use LiDSE
            ekk.info_.store_squared_primal_infeasibility = false;
        }
    }
}
