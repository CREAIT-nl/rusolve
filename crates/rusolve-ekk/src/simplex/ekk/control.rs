//! HiGHS source: simplex/HEkkControl.cpp, simplex/HEkkInterface.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Edge-weight control (DSE -> Devex switch) and vector appends.

use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkkControl.cpp::HEkk::initialiseControl (lines 14-42, 29 lines)
    pub fn initialise_control(&mut self) {
        // Copy tolerances from options
        self.info_.allow_dual_steepest_edge_to_devex_switch =
            self.options_.simplex_dual_edge_weight_strategy == kSimplexEdgeWeightStrategyChoose;
        self.info_.dual_steepest_edge_weight_log_error_threshold =
            self.options_.dual_steepest_edge_weight_log_error_threshold;
        self.info_.control_iteration_count0 = self.iteration_count_;
        // Initialise the densities
        self.info_.col_aq_density = 0.0;
        self.info_.row_ep_density = 0.0;
        self.info_.row_ap_density = 0.0;
        self.info_.row_dse_density = 0.0;
        self.info_.col_steepest_edge_density = 0.0;
        self.info_.col_basic_feasibility_change_density = 0.0;
        self.info_.row_basic_feasibility_change_density = 0.0;
        self.info_.col_bfrt_density = 0.0;
        self.info_.primal_col_density = 0.0;
        // Set the row_dual_density to 1 since it's assumed all costs are at
        // least perturbed from zero, if not initially nonzero
        self.info_.dual_col_density = 1.0;
        // Initialise the data used to determine the switch from DSE to
        // Devex
        self.info_.costly_dse_frequency = 0.0;
        self.info_.num_costly_dse_iteration = 0;
        self.info_.costly_dse_measure = 0.0;
        self.info_.average_log_low_dse_weight_error = 0.0;
        self.info_.average_log_high_dse_weight_error = 0.0;
    }

    /// HiGHS: HEkkControl.cpp::HEkk::updateOperationResultDensity (lines 44-48, 5 lines)
    pub fn update_operation_result_density(&self, local_density: f64, density: &mut f64) {
        *density = (1.0 - kRunningAverageMultiplier) * *density
            + kRunningAverageMultiplier * local_density;
    }

    /// HiGHS: HEkkControl.cpp::HEkk::assessDSEWeightError (lines 50-75, 26 lines)
    pub fn assess_dse_weight_error(&mut self, computed_edge_weight: f64, updated_edge_weight: f64) {
        // Compute the (relative) dual steepest edge weight error for
        // analysis and debugging
        self.edge_weight_error_ =
            (updated_edge_weight - computed_edge_weight).abs() / cpp_max(1.0, computed_edge_weight);
        if self.edge_weight_error_ > self.options_.dual_steepest_edge_weight_error_tolerance {
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Dual steepest edge weight error is {}\n",
                self.edge_weight_error_
            );
        }
        // Compute the relative deviation in the updated weight compared
        // with the computed weight
        let weight_relative_deviation: f64;
        if updated_edge_weight < computed_edge_weight {
            // Updated weight is low
            weight_relative_deviation = computed_edge_weight / updated_edge_weight;
            self.info_.average_log_low_dse_weight_error = 0.99
                * self.info_.average_log_low_dse_weight_error
                + 0.01 * weight_relative_deviation.ln();
        } else {
            // Updated weight is correct or high
            weight_relative_deviation = updated_edge_weight / computed_edge_weight;
            self.info_.average_log_high_dse_weight_error = 0.99
                * self.info_.average_log_high_dse_weight_error
                + 0.01 * weight_relative_deviation.ln();
        }
    }

    /// HiGHS: HEkkControl.cpp::HEkk::switchToDevex (lines 77-146, 70 lines)
    pub fn switch_to_devex(&mut self) -> bool {
        // Parameters controlling switch from DSE to Devex on cost
        let kCostlyDseMeasureLimit: f64 = 1000.0;
        let kCostlyDseMinimumDensity: f64 = 0.01;
        let kCostlyDseFractionNumTotalIterationBeforeSwitch: f64 = 0.1;
        let kCostlyDseFractionNumCostlyDseIterationBeforeSwitch: f64 = 0.05;
        let mut switch_to_devex = false;
        // Firstly consider switching on the basis of NLA cost
        let costly_dse_measure_denominator: f64 = cpp_max(
            cpp_max(self.info_.row_ep_density, self.info_.col_aq_density),
            self.info_.row_ap_density,
        );
        if costly_dse_measure_denominator > 0.0 {
            self.info_.costly_dse_measure =
                self.info_.row_dse_density / costly_dse_measure_denominator;
            self.info_.costly_dse_measure =
                self.info_.costly_dse_measure * self.info_.costly_dse_measure;
        } else {
            self.info_.costly_dse_measure = 0.0;
        }
        let costly_dse_iteration = self.info_.costly_dse_measure > kCostlyDseMeasureLimit
            && self.info_.row_dse_density > kCostlyDseMinimumDensity;
        self.info_.costly_dse_frequency =
            (1.0 - kRunningAverageMultiplier) * self.info_.costly_dse_frequency;
        if costly_dse_iteration {
            self.info_.num_costly_dse_iteration += 1;
            self.info_.costly_dse_frequency += kRunningAverageMultiplier * 1.0;
            // What if non-dual iterations have been performed: need to think about this
            let local_iteration_count: HighsInt =
                self.iteration_count_ - self.info_.control_iteration_count0;
            let local_num_tot: HighsInt = self.lp_.num_col_ + self.lp_.num_row_;
            // Switch to Devex if at least 5% of the (at least) 0.1NumTot iterations
            // have been costly
            switch_to_devex = self.info_.allow_dual_steepest_edge_to_devex_switch
                && (self.info_.num_costly_dse_iteration as f64
                    > (local_iteration_count as f64)
                        * kCostlyDseFractionNumCostlyDseIterationBeforeSwitch)
                && ((local_iteration_count as f64)
                    > kCostlyDseFractionNumTotalIterationBeforeSwitch * (local_num_tot as f64));

            if switch_to_devex {
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Switch from DSE to Devex after {} costly DSE iterations of {} with densities C_Aq = {:11.4}; R_Ep = {:11.4}; R_Ap = {:11.4}; DSE = {:11.4}\n",
                    self.info_.num_costly_dse_iteration,
                    local_iteration_count,
                    self.info_.col_aq_density,
                    self.info_.row_ep_density,
                    self.info_.row_ap_density,
                    self.info_.row_dse_density
                );
            }
        }
        if !switch_to_devex {
            // Secondly consider switching on the basis of weight accuracy
            let local_measure = self.info_.average_log_low_dse_weight_error
                + self.info_.average_log_high_dse_weight_error;
            let local_threshold = self.info_.dual_steepest_edge_weight_log_error_threshold;
            switch_to_devex = self.info_.allow_dual_steepest_edge_to_devex_switch
                && local_measure > local_threshold;
            if switch_to_devex {
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Switch from DSE to Devex with log error measure of {} > {} = threshold\n",
                    local_measure,
                    local_threshold
                );
            }
        }
        switch_to_devex
    }
}

impl HEkk {
    /// HiGHS: HEkkInterface.cpp::HEkk::appendColsToVectors (lines 15-20, 6 lines)
    pub fn append_cols_to_vectors(
        &mut self,
        num_new_col: HighsInt,
        col_cost: &[f64],
        col_lower: &[f64],
        col_upper: &[f64],
    ) {
        append_cols_to_lp_vectors(&mut self.lp_, num_new_col, col_cost, col_lower, col_upper);
    }

    /// HiGHS: HEkkInterface.cpp::HEkk::appendRowsToVectors (lines 22-26, 5 lines)
    pub fn append_rows_to_vectors(
        &mut self,
        num_new_row: HighsInt,
        row_lower: &[f64],
        row_upper: &[f64],
    ) {
        append_rows_to_lp_vectors(&mut self.lp_, num_new_row, row_lower, row_upper);
    }
}
