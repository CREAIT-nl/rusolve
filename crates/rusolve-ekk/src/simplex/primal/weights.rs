//! HiGHS source: simplex/HEkkPrimal.cpp (lines 2140-2573)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Basic-feasibility-change updates, Devex and primal steepest-edge weights.
//! See `HEkkPrimal` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkPrimal {
    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::basicFeasibilityChangeUpdateDual (lines 2140-2178, 39 lines)
    pub fn basic_feasibility_change_update_dual(&mut self, ekk: &mut HEkk) {
        ekk.analysis_
            .simplex_timer_start(UpdateDualBasicFeasibilityChangeClock, 0);
        // See HEkkPrimal.cpp for the comment on basic logicals/structurals.
        self.basic_feasibility_change_btran(ekk);
        self.basic_feasibility_change_price(ekk);
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
            ekk.info_.work_dual_[i_col as usize] -=
                self.row_basic_feasibility_change.array[i_col as usize];
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
            ekk.info_.work_dual_[i_col as usize] -=
                self.col_basic_feasibility_change.array[i_row as usize];
        }
        ekk.invalidate_dual_infeasibility_record();
        ekk.analysis_
            .simplex_timer_stop(UpdateDualBasicFeasibilityChangeClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::basicFeasibilityChangeBtran (lines 2180-2205, 26 lines)
    pub fn basic_feasibility_change_btran(&mut self, ekk: &mut HEkk) {
        ekk.analysis_
            .simplex_timer_start(BtranBasicFeasibilityChangeClock, 0);
        let solver_num_row = ekk.lp_.num_row_;
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.operation_record_before(
                kSimplexNlaBtranBasicFeasibilityChange,
                &self.col_basic_feasibility_change,
                ekk.info_.col_basic_feasibility_change_density,
            );
        }
        {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            ekk.simplex_nla_.btran(
                &refs,
                &mut self.col_basic_feasibility_change,
                ekk.info_.col_basic_feasibility_change_density,
                None,
            );
        }
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.operation_record_after(
                kSimplexNlaBtranBasicFeasibilityChange,
                &self.col_basic_feasibility_change,
            );
        }
        let local_col_basic_feasibility_change_density =
            self.col_basic_feasibility_change.count as f64 / solver_num_row as f64;
        let mut density = ekk.info_.col_basic_feasibility_change_density;
        ekk.update_operation_result_density(
            local_col_basic_feasibility_change_density,
            &mut density,
        );
        ekk.info_.col_basic_feasibility_change_density = density;
        ekk.analysis_
            .simplex_timer_stop(BtranBasicFeasibilityChangeClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::basicFeasibilityChangePrice (lines 2207-2276, 70 lines)
    pub fn basic_feasibility_change_price(&mut self, ekk: &mut HEkk) {
        ekk.analysis_
            .simplex_timer_start(PriceBasicFeasibilityChangeClock, 0);
        let local_density =
            1.0 * self.col_basic_feasibility_change.count as f64 / self.num_row as f64;
        let mut use_col_price = false;
        let mut use_row_price_w_switch = false;
        ekk.choose_price_technique(
            ekk.info_.price_strategy,
            local_density,
            &mut use_col_price,
            &mut use_row_price_w_switch,
        );
        if ekk.analysis_.analyse_simplex_summary_data {
            if use_col_price {
                let expected_density = 1.0;
                ekk.analysis_.operation_record_before(
                    kSimplexNlaPriceBasicFeasibilityChange,
                    &self.col_basic_feasibility_change,
                    expected_density,
                );
                ekk.analysis_.num_col_price += 1;
            } else if use_row_price_w_switch {
                ekk.analysis_.operation_record_before(
                    kSimplexNlaPriceBasicFeasibilityChange,
                    &self.col_basic_feasibility_change,
                    ekk.info_.col_basic_feasibility_change_density,
                );
                ekk.analysis_.num_row_price_with_switch += 1;
            } else {
                ekk.analysis_.operation_record_before(
                    kSimplexNlaPriceBasicFeasibilityChange,
                    &self.col_basic_feasibility_change,
                    ekk.info_.col_basic_feasibility_change_density,
                );
                ekk.analysis_.num_row_price += 1;
            }
        }
        self.row_basic_feasibility_change.clear();
        let quad_precision = false;
        if use_col_price {
            // Perform column-wise PRICE
            ekk.lp_.a_matrix_.price_by_column(
                quad_precision,
                &mut self.row_basic_feasibility_change,
                &self.col_basic_feasibility_change,
                kDebugReportOff,
            );
        } else if use_row_price_w_switch {
            // Perform hyper-sparse row-wise PRICE, but switch if the density of
            // row_basic_feasibility_change becomes extreme
            let switch_density = kHyperPriceDensity;
            ekk.ar_matrix_.price_by_row_with_switch(
                quad_precision,
                &mut self.row_basic_feasibility_change,
                &self.col_basic_feasibility_change,
                ekk.info_.row_basic_feasibility_change_density,
                0,
                switch_density,
                kDebugReportOff,
            );
        } else {
            // Perform hyper-sparse row-wise PRICE
            ekk.ar_matrix_.price_by_row(
                quad_precision,
                &mut self.row_basic_feasibility_change,
                &self.col_basic_feasibility_change,
                kDebugReportOff,
            );
        }
        if use_col_price {
            // Column-wise PRICE computes components corresponding to basic
            // variables, so zero these by exploiting the fact that, for basic
            // variables, nonbasicFlag[*]=0
            let nonbasic_flag = &ekk.basis_.nonbasic_flag_;
            for i_col in 0..self.num_col {
                self.row_basic_feasibility_change.array[i_col as usize] *=
                    nonbasic_flag[i_col as usize] as f64;
            }
        }
        // Update the record of average row_basic_feasibility_change density
        let local_row_basic_feasibility_change_density =
            self.row_basic_feasibility_change.count as f64 / self.num_col as f64;
        let mut density = ekk.info_.row_basic_feasibility_change_density;
        ekk.update_operation_result_density(
            local_row_basic_feasibility_change_density,
            &mut density,
        );
        ekk.info_.row_basic_feasibility_change_density = density;
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.operation_record_after(
                kSimplexNlaPriceBasicFeasibilityChange,
                &self.row_basic_feasibility_change,
            );
        }
        ekk.analysis_
            .simplex_timer_stop(PriceBasicFeasibilityChangeClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::initialiseDevexFramework (lines 2278-2289, 12 lines)
    pub fn initialise_devex_framework(&mut self, ekk: &mut HEkk) {
        self.edge_weight_.clear();
        self.edge_weight_.resize(self.num_tot as usize, 1.0);
        self.devex_index_.clear();
        self.devex_index_.resize(self.num_tot as usize, 0);
        for i_col in 0..self.num_tot {
            let nonbasic_flag = ekk.basis_.nonbasic_flag_[i_col as usize] as HighsInt;
            self.devex_index_[i_col as usize] = nonbasic_flag * nonbasic_flag;
        }
        self.num_devex_iterations_ = 0;
        self.num_bad_devex_weight_ = 0;
        if self.report_hyper_chuzc {
            println!("initialiseDevexFramework");
        }
        self.hyper_choose_column_clear(ekk);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::updateDevex (lines 2291-2338, 48 lines)
    pub fn update_devex(&mut self, ekk: &mut HEkk) {
        ekk.analysis_.simplex_timer_start(DevexUpdateWeightClock, 0);
        // Compute the pivot weight from the reference set
        let mut d_pivot_weight = 0.0;
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
            let i_col = ekk.basis_.basic_index_[i_row as usize];
            let d_alpha =
                self.devex_index_[i_col as usize] as f64 * self.col_aq.array[i_row as usize];
            d_pivot_weight += d_alpha * d_alpha;
        }
        d_pivot_weight += self.devex_index_[self.variable_in as usize] as f64 * 1.0;

        // Check if the saved weight is too large
        if self.edge_weight_[self.variable_in as usize] > kBadDevexWeightFactor * d_pivot_weight {
            self.num_bad_devex_weight_ += 1;
        }

        // Update the devex weight for all
        let d_pivot = self.col_aq.array[self.row_out as usize];
        d_pivot_weight /= d_pivot * d_pivot;

        for i_el in 0..self.row_ap.count {
            let i_col = self.row_ap.index[i_el as usize];
            let alpha = self.row_ap.array[i_col as usize];
            let mut devex = d_pivot_weight * alpha * alpha;
            devex += self.devex_index_[i_col as usize] as f64 * 1.0;
            if self.edge_weight_[i_col as usize] < devex {
                self.edge_weight_[i_col as usize] = devex;
            }
        }
        for i_el in 0..self.row_ep.count {
            let i_row = self.row_ep.index[i_el as usize];
            let i_col = i_row + self.num_col;
            let alpha = self.row_ep.array[i_row as usize];
            let mut devex = d_pivot_weight * alpha * alpha;
            devex += self.devex_index_[i_col as usize] as f64 * 1.0;
            if self.edge_weight_[i_col as usize] < devex {
                self.edge_weight_[i_col as usize] = devex;
            }
        }
        // Update devex weight for the pivots
        self.edge_weight_[self.variable_out as usize] = cpp_max(1.0, d_pivot_weight);
        self.edge_weight_[self.variable_in as usize] = 1.0;
        self.num_devex_iterations_ += 1;
        ekk.analysis_.simplex_timer_stop(DevexUpdateWeightClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::computePrimalSteepestEdgeWeights (lines 2340-2368, 29 lines)
    pub fn compute_primal_steepest_edge_weights(&mut self, ekk: &mut HEkk) {
        let report_var: HighsInt = -16;
        self.edge_weight_.resize(self.num_tot as usize, 0.0);
        if ekk.logical_basis() {
            let a_matrix = &ekk.lp_.a_matrix_;
            for i_col in 0..self.num_col {
                self.edge_weight_[i_col as usize] = 1.0;
                for i_el in a_matrix.start_[i_col as usize]..a_matrix.start_[(i_col + 1) as usize] {
                    self.edge_weight_[i_col as usize] +=
                        a_matrix.value_[i_el as usize] * a_matrix.value_[i_el as usize];
                }
            }
        } else {
            let mut local_col_aq = HVector::default();
            local_col_aq.setup(self.num_row);
            for i_var in 0..self.num_tot {
                if ekk.basis_.nonbasic_flag_[i_var as usize] != 0 {
                    self.edge_weight_[i_var as usize] =
                        self.compute_primal_steepest_edge_weight(ekk, i_var, &mut local_col_aq);
                    if i_var == report_var {
                        println!("Tableau column {}\nRow       Value", report_var);
                        for i_row in 0..self.num_row {
                            if local_col_aq.array[i_row as usize] != 0.0 {
                                println!(
                                    "{:3}  {:10.7}",
                                    i_row, local_col_aq.array[i_row as usize]
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::computePrimalSteepestEdgeWeight (lines 2370-2383, 14 lines)
    pub fn compute_primal_steepest_edge_weight(
        &mut self,
        ekk: &mut HEkk,
        i_var: HighsInt,
        local_col_aq: &mut HVector,
    ) -> f64 {
        local_col_aq.clear();
        ekk.lp_.a_matrix_.collect_aj(local_col_aq, i_var, 1.0);
        local_col_aq.pack_flag = false;
        {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            ekk.simplex_nla_
                .ftran(&refs, local_col_aq, ekk.info_.col_aq_density, None);
        }
        let local_col_aq_density = (1.0 * local_col_aq.count as f64) / ekk.lp_.num_row_ as f64;
        let mut density = ekk.info_.col_aq_density;
        ekk.update_operation_result_density(local_col_aq_density, &mut density);
        ekk.info_.col_aq_density = density;
        1.0 + local_col_aq.norm2()
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::updatePrimalSteepestEdgeWeights (lines 2385-2470, 86 lines)
    pub fn update_primal_steepest_edge_weights(&mut self, ekk: &mut HEkk) {
        // Compute, for all j: lambda_j = hat{a}_{pj} / hat{a}_{pq}; see HEkkPrimal.cpp.
        self.col_steepest_edge.copy(&self.col_aq);
        let mut col_steepest_edge = std::mem::take(&mut self.col_steepest_edge);
        self.update_btran_pse(ekk, &mut col_steepest_edge);
        let col_aq_squared_2norm = self.col_aq.norm2();
        let report_col_aq = false;
        if report_col_aq {
            println!(
                "updatePrimalSteepestEdgeWeights: in = {}; out = {}; ||col_aq||^2 = {}",
                self.variable_in, self.variable_out, col_aq_squared_2norm
            );
            println!("Pivotal column {}\nRow       Value", self.variable_in);
            for i_row in 0..self.num_row {
                if self.col_aq.array[i_row as usize] != 0.0 {
                    println!("{:3}  {:10.7}", i_row, self.col_aq.array[i_row as usize]);
                }
            }
        }
        debug_assert!(ekk.basis_.nonbasic_flag_[self.variable_in as usize] != 0);
        let a_matrix = &ekk.lp_.a_matrix_;
        let mut i_var: HighsInt;
        let mut pivotal_row_value: f64;
        for i_x in 0..(self.row_ap.count + self.row_ep.count) {
            if i_x < self.row_ap.count {
                i_var = self.row_ap.index[i_x as usize];
                pivotal_row_value = self.row_ap.array[i_var as usize];
            } else {
                let i_row = self.row_ep.index[(i_x - self.row_ap.count) as usize];
                i_var = self.num_col + i_row;
                pivotal_row_value = self.row_ep.array[i_row as usize];
            }
            if i_var == self.variable_in {
                continue;
            }
            if ekk.basis_.nonbasic_flag_[i_var as usize] == 0 {
                continue;
            }
            let lambda = pivotal_row_value / self.alpha_col;
            let mut mu_aj = 0.0;
            if i_var < self.num_col {
                for i_el in a_matrix.start_[i_var as usize]..a_matrix.start_[(i_var + 1) as usize] {
                    mu_aj += col_steepest_edge.array[a_matrix.index_[i_el as usize] as usize]
                        * a_matrix.value_[i_el as usize];
                }
            } else {
                mu_aj = col_steepest_edge.array[(i_var - self.num_col) as usize];
            }
            let min_weight = 1.0 + lambda * lambda;
            self.edge_weight_[i_var as usize] +=
                lambda * lambda * col_aq_squared_2norm - 2.0 * lambda * mu_aj;
            self.edge_weight_[i_var as usize] += lambda * lambda;
            if self.edge_weight_[i_var as usize] < min_weight {
                self.edge_weight_[i_var as usize] = min_weight;
            }
        }
        // New weight of the leaving variable: (col_aq_squared_2norm + 1) / pivot^2
        self.edge_weight_[self.variable_out as usize] =
            (1.0 + col_aq_squared_2norm) / (self.alpha_col * self.alpha_col);
        self.edge_weight_[self.variable_in as usize] = 0.0;
        self.col_steepest_edge = col_steepest_edge;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::updateDualSteepestEdgeWeights (lines 2472-2495, 24 lines)
    pub fn update_dual_steepest_edge_weights(&mut self, ekk: &mut HEkk) {
        self.col_steepest_edge.copy(&self.row_ep);
        let mut col_steepest_edge = std::mem::take(&mut self.col_steepest_edge);
        self.update_ftran_dse(ekk, &mut col_steepest_edge);
        // Compute the weight from row_ep and over-write the updated weight
        let row_out = self.row_out;
        if ekk.simplex_in_scaled_space_ {
            ekk.dual_edge_weight_[row_out as usize] = self.row_ep.norm2();
        } else {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            ekk.dual_edge_weight_[row_out as usize] = ekk
                .simplex_nla_
                .row_ep2_norm_in_scaled_space(&refs, row_out, &self.row_ep);
        }
        let pivot_in_scaled_space = {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            ekk.simplex_nla_
                .pivot_in_scaled_space(&refs, &self.col_aq, self.variable_in, row_out)
        };
        if ekk.simplex_in_scaled_space_ {
            debug_assert!(pivot_in_scaled_space == self.alpha_col);
        }
        let new_pivotal_edge_weight = ekk.dual_edge_weight_[row_out as usize]
            / (pivot_in_scaled_space * pivot_in_scaled_space);
        let kai = -2.0 / pivot_in_scaled_space;
        ekk.update_dual_steepest_edge_weights(
            row_out,
            self.variable_in,
            &self.col_aq,
            new_pivotal_edge_weight,
            kai,
            &col_steepest_edge.array,
        );
        ekk.dual_edge_weight_[row_out as usize] = new_pivotal_edge_weight;
        self.col_steepest_edge = col_steepest_edge;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::updateFtranDSE (lines 2497-2517, 21 lines)
    pub fn update_ftran_dse(&mut self, ekk: &mut HEkk, col_steepest_edge: &mut HVector) {
        // For comments on scaling actions, see HEkkDual::updateFtranDSE
        ekk.analysis_.simplex_timer_start(FtranDseClock, 0);
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.operation_record_before(
                kSimplexNlaFtranDse,
                col_steepest_edge,
                ekk.info_.row_dse_density,
            );
        }
        {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            // Apply R{-1}
            ekk.simplex_nla_
                .unapply_basis_matrix_row_scale(&refs, col_steepest_edge);
            // Perform FTRAN DSE
            ekk.simplex_nla_.ftran_in_scaled_space(
                col_steepest_edge,
                ekk.info_.row_dse_density,
                None,
            );
        }
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_
                .operation_record_after(kSimplexNlaFtranDse, col_steepest_edge);
        }
        ekk.analysis_.simplex_timer_stop(FtranDseClock, 0);
        let local_row_dse_density = (1.0 * col_steepest_edge.count as f64) / self.num_row as f64;
        let mut density = ekk.info_.row_dse_density;
        ekk.update_operation_result_density(local_row_dse_density, &mut density);
        ekk.info_.row_dse_density = density;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::updateBtranPSE (lines 2519-2537, 19 lines)
    pub fn update_btran_pse(&mut self, ekk: &mut HEkk, col_steepest_edge: &mut HVector) {
        ekk.analysis_.simplex_timer_start(BtranPseClock, 0);
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.operation_record_before(
                kSimplexNlaBtranPse,
                col_steepest_edge,
                ekk.info_.col_steepest_edge_density,
            );
        }
        // Perform BTRAN PSE
        {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            ekk.simplex_nla_.btran(
                &refs,
                col_steepest_edge,
                ekk.info_.col_steepest_edge_density,
                None,
            );
        }
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_
                .operation_record_after(kSimplexNlaBtranPse, col_steepest_edge);
        }
        ekk.analysis_.simplex_timer_stop(BtranPseClock, 0);
        let local_col_steepest_edge_density =
            (1.0 * col_steepest_edge.count as f64) / self.num_row as f64;
        let mut density = ekk.info_.col_steepest_edge_density;
        ekk.update_operation_result_density(local_col_steepest_edge_density, &mut density);
        ekk.info_.col_steepest_edge_density = density;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::updateVerify (lines 2539-2573, 35 lines)
    pub fn update_verify(&mut self, ekk: &mut HEkk) {
        // updateVerify for primal
        let numerical_trouble_tolerance = 1e-7;
        self.numerical_trouble = 0.0;
        let abs_alpha_from_col = self.alpha_col.abs();
        let alpha_row_source: &str = if self.variable_in < self.num_col {
            self.alpha_row = self.row_ap.array[self.variable_in as usize];
            "Col"
        } else {
            self.alpha_row = self.row_ep.array[(self.variable_in - self.num_col) as usize];
            "Row"
        };
        let abs_alpha_from_row = self.alpha_row.abs();
        let abs_alpha_diff = (abs_alpha_from_col - abs_alpha_from_row).abs();
        let min_abs_alpha = cpp_min(abs_alpha_from_col, abs_alpha_from_row);
        self.numerical_trouble = abs_alpha_diff / min_abs_alpha;
        if self.numerical_trouble > numerical_trouble_tolerance {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "Numerical check: Iter {:4}: alpha_col = {:12}, (From {:3} alpha_row = {:12}), aDiff = {:12}: measure = {:12}\n",
                ekk.iteration_count_,
                self.alpha_col,
                alpha_row_source,
                self.alpha_row,
                abs_alpha_diff,
                self.numerical_trouble
            );
        }
        if kAllowDeveloperAssert {
            debug_assert!(self.numerical_trouble < 1e-3);
        }
        // Reinvert if the relative difference is large enough, and updates have been
        // performed
        if self.numerical_trouble > 1e-7 && ekk.info_.update_count > 0 {
            self.rebuild_reason = kRebuildReasonPossiblySingularBasis;
        }
    }
}
