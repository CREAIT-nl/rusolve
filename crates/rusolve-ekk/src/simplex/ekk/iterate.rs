//! HiGHS source: simplex/HEkk.cpp (lines 2774-3442)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Per-iteration kernels: FTRAN/BTRAN/PRICE, primal/dual computation, factor and
//! pivot updates, bad-basis-change detection, infeasibility computation.

use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkk.cpp::HEkk::pivotColumnFtran (lines 2774-2790, 17 lines)
    pub fn pivot_column_ftran(&mut self, i_col: HighsInt, col_aq: &mut HVector) {
        self.analysis_.simplex_timer_start(FtranClock, 0);
        col_aq.clear();
        col_aq.pack_flag = true;
        self.lp_.a_matrix_.collect_aj(col_aq, i_col, 1.0);
        if self.analysis_.analyse_simplex_summary_data {
            self.analysis_.operation_record_before(
                kSimplexNlaFtran,
                col_aq,
                self.info_.col_aq_density,
            );
        }
        {
            let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
            self.simplex_nla_
                .ftran(&refs, col_aq, self.info_.col_aq_density, None);
        }
        if self.analysis_.analyse_simplex_summary_data {
            self.analysis_
                .operation_record_after(kSimplexNlaFtran, col_aq);
        }
        let num_row = self.lp_.num_row_;
        let local_col_aq_density = col_aq.count as f64 / num_row as f64;
        let mut density = self.info_.col_aq_density;
        self.update_operation_result_density(local_col_aq_density, &mut density);
        self.info_.col_aq_density = density;
        self.analysis_.simplex_timer_stop(FtranClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::unitBtran (lines 2792-2810, 19 lines)
    pub fn unit_btran(&mut self, i_row: HighsInt, row_ep: &mut HVector) {
        self.analysis_.simplex_timer_start(BtranClock, 0);
        row_ep.clear();
        row_ep.count = 1;
        row_ep.index[0] = i_row;
        row_ep.array[i_row as usize] = 1.0;
        row_ep.pack_flag = true;
        if self.analysis_.analyse_simplex_summary_data {
            self.analysis_.operation_record_before(
                kSimplexNlaBtranEp,
                row_ep,
                self.info_.row_ep_density,
            );
        }
        {
            let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
            self.simplex_nla_
                .btran(&refs, row_ep, self.info_.row_ep_density, None);
        }
        if self.analysis_.analyse_simplex_summary_data {
            self.analysis_
                .operation_record_after(kSimplexNlaBtranEp, row_ep);
        }
        let num_row = self.lp_.num_row_;
        let local_row_ep_density = row_ep.count as f64 / num_row as f64;
        let mut density = self.info_.row_ep_density;
        self.update_operation_result_density(local_row_ep_density, &mut density);
        self.info_.row_ep_density = density;
        self.analysis_.simplex_timer_stop(BtranClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::fullBtran (lines 2812-2828, 17 lines)
    pub fn full_btran(&mut self, buffer: &mut HVector) {
        // Performs BTRAN on the buffer supplied. Make sure that
        // buffer.count is large (>lp_.num_row_ to be sure) rather
        // than 0 if the indices of the RHS (and true value of buffer.count)
        // isn't known.
        self.analysis_.simplex_timer_start(BtranFullClock, 0);
        if self.analysis_.analyse_simplex_summary_data {
            self.analysis_.operation_record_before(
                kSimplexNlaBtranFull,
                buffer,
                self.info_.dual_col_density,
            );
        }
        {
            let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
            self.simplex_nla_
                .btran(&refs, buffer, self.info_.dual_col_density, None);
        }
        if self.analysis_.analyse_simplex_summary_data {
            self.analysis_
                .operation_record_after(kSimplexNlaBtranFull, buffer);
        }
        let local_dual_col_density = buffer.count as f64 / self.lp_.num_row_ as f64;
        let mut density = self.info_.dual_col_density;
        self.update_operation_result_density(local_dual_col_density, &mut density);
        self.info_.dual_col_density = density;
        self.analysis_.simplex_timer_stop(BtranFullClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::choosePriceTechnique (lines 2830-2843, 14 lines)
    pub fn choose_price_technique(
        &self,
        price_strategy: HighsInt,
        row_ep_density: f64,
        use_col_price: &mut bool,
        use_row_price_w_switch: &mut bool,
    ) {
        // By default switch to column PRICE when pi_p has at least this
        // density
        let density_for_column_price_switch = 0.75;
        *use_col_price = (price_strategy == kSimplexPriceStrategyCol)
            || (price_strategy == kSimplexPriceStrategyRowSwitchColSwitch
                && row_ep_density > density_for_column_price_switch);
        *use_row_price_w_switch = price_strategy == kSimplexPriceStrategyRowSwitch
            || price_strategy == kSimplexPriceStrategyRowSwitchColSwitch;
    }

    /// HiGHS: HEkk.cpp::HEkk::tableauRowPrice (lines 2845-2900, 56 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    pub fn tableau_row_price(
        &mut self,
        quad_precision: bool,
        row_ep: &HVector,
        row_ap: &mut HVector,
        debug_report: HighsInt,
    ) {
        self.analysis_.simplex_timer_start(PriceClock, 0);
        let solver_num_row = self.lp_.num_row_;
        let solver_num_col = self.lp_.num_col_;
        let local_density = 1.0 * row_ep.count as f64 / solver_num_row as f64;
        let mut use_col_price = false;
        let mut use_row_price_w_switch = false;
        self.choose_price_technique(
            self.info_.price_strategy,
            local_density,
            &mut use_col_price,
            &mut use_row_price_w_switch,
        );
        if self.analysis_.analyse_simplex_summary_data {
            if use_col_price {
                let expected_density = 1.0;
                self.analysis_.operation_record_before(
                    kSimplexNlaPriceAp,
                    row_ep,
                    expected_density,
                );
                self.analysis_.num_col_price += 1;
            } else if use_row_price_w_switch {
                self.analysis_.operation_record_before(
                    kSimplexNlaPriceAp,
                    row_ep,
                    self.info_.row_ep_density,
                );
                self.analysis_.num_row_price_with_switch += 1;
            } else {
                self.analysis_.operation_record_before(
                    kSimplexNlaPriceAp,
                    row_ep,
                    self.info_.row_ep_density,
                );
                self.analysis_.num_row_price += 1;
            }
        }
        row_ap.clear();
        if use_col_price {
            // Perform column-wise PRICE
            self.lp_
                .a_matrix_
                .price_by_column(quad_precision, row_ap, row_ep, debug_report);
        } else if use_row_price_w_switch {
            // Perform hyper-sparse row-wise PRICE, but switch if the density of row_ap
            // becomes extreme
            let switch_density = kHyperPriceDensity;
            self.ar_matrix_.price_by_row_with_switch(
                quad_precision,
                row_ap,
                row_ep,
                self.info_.row_ap_density,
                0,
                switch_density,
                debug_report,
            );
        } else {
            // Perform hyper-sparse row-wise PRICE
            self.ar_matrix_
                .price_by_row(quad_precision, row_ap, row_ep, debug_report);
        }
        if use_col_price {
            // Column-wise PRICE computes components corresponding to basic
            // variables, so zero these by exploiting the fact that, for basic
            // variables, nonbasicFlag[*]=0
            let nonbasic_flag = &self.basis_.nonbasic_flag_;
            for i_col in 0..solver_num_col as usize {
                row_ap.array[i_col] *= nonbasic_flag[i_col] as f64;
            }
        }
        // Update the record of average row_ap density
        let local_row_ap_density = row_ap.count as f64 / solver_num_col as f64;
        let mut density = self.info_.row_ap_density;
        self.update_operation_result_density(local_row_ap_density, &mut density);
        self.info_.row_ap_density = density;
        if self.analysis_.analyse_simplex_summary_data {
            self.analysis_
                .operation_record_after(kSimplexNlaPriceAp, row_ap);
        }
        self.analysis_.simplex_timer_stop(PriceClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::fullPrice (lines 2902-2915, 14 lines)
    pub fn full_price(&mut self, full_col: &HVector, full_row: &mut HVector) {
        self.analysis_.simplex_timer_start(PriceFullClock, 0);
        full_row.clear();
        if self.analysis_.analyse_simplex_summary_data {
            let expected_density = 1.0;
            self.analysis_.operation_record_before(
                kSimplexNlaPriceFull,
                full_col,
                expected_density,
            );
        }
        let quad_precision = false;
        self.lp_
            .a_matrix_
            .price_by_column(quad_precision, full_row, full_col, kDebugReportOff);
        if self.analysis_.analyse_simplex_summary_data {
            self.analysis_
                .operation_record_after(kSimplexNlaPriceFull, full_row);
        }
        self.analysis_.simplex_timer_stop(PriceFullClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::computePrimal (lines 2917-2952, 36 lines)
    pub fn compute_primal(&mut self) {
        self.analysis_.simplex_timer_start(ComputePrimalClock, 0);
        let num_row = self.lp_.num_row_;
        let num_col = self.lp_.num_col_;
        // Setup a local buffer for the values of basic variables
        // PORT NOTE (perf): recycled storage, same values as a fresh vector (HVectorPool).
        let mut primal_col = self.hvector_pool_.take();
        primal_col.setup(num_row);
        primal_col.clear();
        for i in 0..num_col + num_row {
            if self.basis_.nonbasic_flag_[i as usize] != 0
                && self.info_.work_value_[i as usize] != 0.0
            {
                self.lp_.a_matrix_.collect_aj(
                    &mut primal_col,
                    i,
                    self.info_.work_value_[i as usize],
                );
            }
        }
        // It's possible that the buffer has no nonzeros, so performing
        // FTRAN is unnecessary. Not much of a saving, but the zero density
        // looks odd in the analysis!
        if primal_col.count != 0 {
            {
                let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
                self.simplex_nla_.ftran(
                    &refs,
                    &mut primal_col,
                    self.info_.primal_col_density,
                    None,
                );
            }
            let local_primal_col_density = primal_col.count as f64 / num_row as f64;
            let mut density = self.info_.primal_col_density;
            self.update_operation_result_density(local_primal_col_density, &mut density);
            self.info_.primal_col_density = density;
        }
        for i in 0..num_row as usize {
            let i_col = self.basis_.basic_index_[i] as usize;
            self.info_.base_value_[i] = -primal_col.array[i];
            self.info_.base_lower_[i] = self.info_.work_lower_[i_col];
            self.info_.base_upper_[i] = self.info_.work_upper_[i_col];
        }
        self.hvector_pool_.give(primal_col);
        // Indicate that the primal infeasibility information isn't known
        self.info_.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.info_.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;

        self.analysis_.simplex_timer_stop(ComputePrimalClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::computeDual (lines 2954-3000, 47 lines)
    pub fn compute_dual(&mut self) {
        self.analysis_.simplex_timer_start(ComputeDualClock, 0);
        // Create a local buffer for the pi vector
        // PORT NOTE (perf): recycled storage, same values as a fresh vector (HVectorPool).
        let mut dual_col = self.hvector_pool_.take();
        dual_col.setup(self.lp_.num_row_);
        dual_col.clear();
        for i_row in 0..self.lp_.num_row_ {
            let b = self.basis_.basic_index_[i_row as usize] as usize;
            let value = self.info_.work_cost_[b] + self.info_.work_shift_[b];
            if value != 0.0 {
                dual_col.index[dual_col.count as usize] = i_row;
                dual_col.count += 1;
                dual_col.array[i_row as usize] = value;
            }
        }
        // If debugging, save the current duals
        let debug_compute_dual = false;
        if debug_compute_dual {
            self.debug_compute_dual(true);
            self.debug_simplex_dual_infeasible("(old duals)", true);
        }
        // Copy the costs in case the basic costs are all zero. Zipped subslices: same
        // element-wise operations as the C++ index loops, without a bounds check per element.
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        let num_tot_usize = num_tot as usize;
        {
            let info = &mut self.info_;
            let cost = &info.work_cost_[..num_tot_usize];
            let shift = &info.work_shift_[..num_tot_usize];
            for ((dual, cost_i), shift_i) in info.work_dual_[..num_tot_usize]
                .iter_mut()
                .zip(cost)
                .zip(shift)
            {
                *dual = *cost_i + *shift_i;
            }
        }

        if dual_col.count != 0 {
            self.full_btran(&mut dual_col);
            // Create a local buffer for the values of reduced costs
            let mut dual_row = self.hvector_pool_.take();
            dual_row.setup(self.lp_.num_col_);
            self.full_price(&dual_col, &mut dual_row);
            let num_col = self.lp_.num_col_ as usize;
            for (dual, &r) in self.info_.work_dual_[..num_col]
                .iter_mut()
                .zip(&dual_row.array[..num_col])
            {
                *dual -= r;
            }
            let num_row = self.lp_.num_row_ as usize;
            for (dual, &r) in self.info_.work_dual_[num_col..num_tot_usize]
                .iter_mut()
                .zip(&dual_col.array[..num_row])
            {
                *dual -= r;
            }
            if debug_compute_dual {
                self.debug_compute_dual(false);
                self.debug_simplex_dual_infeasible("(new duals)", true);
            }
            self.hvector_pool_.give(dual_row);
        }
        self.hvector_pool_.give(dual_col);
        // Indicate that the dual infeasibility information isn't known
        self.info_.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.info_.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;

        self.analysis_.simplex_timer_stop(ComputeDualClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::computeDualForTableauColumn (lines 3002-3013, 12 lines)
    pub fn compute_dual_for_tableau_column(
        &self,
        i_var: HighsInt,
        tableau_column: &HVector,
    ) -> f64 {
        let work_cost = &self.info_.work_cost_;
        let basic_index = &self.basis_.basic_index_;

        let mut dual = self.info_.work_cost_[i_var as usize];
        for i in 0..tableau_column.count as usize {
            let i_row = tableau_column.index[i] as usize;
            dual -= tableau_column.array[i_row] * work_cost[basic_index[i_row] as usize];
        }
        dual
    }

    /// HiGHS: HEkk.cpp::HEkk::reinvertOnNumericalTrouble (lines 3015-3059, 45 lines)
    pub fn reinvert_on_numerical_trouble(
        &mut self,
        method_name: &str,
        numerical_trouble_measure: &mut f64,
        alpha_from_col: f64,
        alpha_from_row: f64,
        numerical_trouble_tolerance: f64,
    ) -> bool {
        let abs_alpha_from_col = alpha_from_col.abs();
        let abs_alpha_from_row = alpha_from_row.abs();
        let min_abs_alpha = cpp_min(abs_alpha_from_col, abs_alpha_from_row);
        let abs_alpha_diff = (abs_alpha_from_col - abs_alpha_from_row).abs();
        *numerical_trouble_measure = abs_alpha_diff / min_abs_alpha;
        let update_count = self.info_.update_count;
        // Reinvert if the relative difference is large enough, and updates have been
        // performed
        let numerical_trouble = *numerical_trouble_measure > numerical_trouble_tolerance;
        let reinvert = numerical_trouble && update_count > 0;
        self.debug_report_reinvert_on_numerical_trouble(
            method_name,
            *numerical_trouble_measure,
            alpha_from_col,
            alpha_from_row,
            numerical_trouble_tolerance,
            reinvert,
        );
        if reinvert {
            // Consider increasing the Markowitz multiplier
            let current_pivot_threshold = self.info_.factor_pivot_threshold;
            let mut new_pivot_threshold = 0.0;
            if current_pivot_threshold < kDefaultPivotThreshold {
                // Threshold is below default value, so increase it
                new_pivot_threshold = cpp_min(
                    current_pivot_threshold * kPivotThresholdChangeFactor,
                    kDefaultPivotThreshold,
                );
            } else if current_pivot_threshold < kMaxPivotThreshold {
                // Threshold is below max value, so increase it if few updates have been
                // performed
                if update_count < 10 {
                    new_pivot_threshold = cpp_min(
                        current_pivot_threshold * kPivotThresholdChangeFactor,
                        kMaxPivotThreshold,
                    );
                }
            }
            if new_pivot_threshold != 0.0 {
                highs_log_user(
                    &self.options_.log_options,
                    HighsLogType::kWarning,
                    &format!(
                        "   Increasing Markowitz threshold to {}\n",
                        new_pivot_threshold
                    ),
                );
                self.info_.factor_pivot_threshold = new_pivot_threshold;
                self.simplex_nla_.set_pivot_threshold(new_pivot_threshold);
            }
        }
        reinvert
    }

    /// HiGHS: HEkk.cpp::HEkk::transformForUpdate (lines 3064-3067, 4 lines)
    pub fn transform_for_update(
        &mut self,
        column: &mut HVector,
        row_ep: &mut HVector,
        variable_in: HighsInt,
        row_out: &mut HighsInt,
    ) {
        let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
        self.simplex_nla_
            .transform_for_update(&refs, column, row_ep, variable_in, *row_out);
    }

    /// HiGHS: HEkk.cpp::HEkk::flipBound (lines 3069-3073, 5 lines)
    pub fn flip_bound(&mut self, i_col: HighsInt) {
        self.basis_.nonbasic_move_[i_col as usize] = -self.basis_.nonbasic_move_[i_col as usize];
        let r#move = self.basis_.nonbasic_move_[i_col as usize];
        self.info_.work_value_[i_col as usize] = if r#move == 1 {
            self.info_.work_lower_[i_col as usize]
        } else {
            self.info_.work_upper_[i_col as usize]
        };
    }

    /// HiGHS: HEkk.cpp::HEkk::updateFactor (lines 3075-3103, 29 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion): `total_synthetic_tick_ >= build_synthetic_tick_` sets kRebuildReasonSyntheticClockSaysInvert.
    /// `simplex_nla_.update(column, row_ep, iRow, hint)` -> `self.simplex_nla_.update(&mut self.basis_.basic_index_, ...)`.
    pub fn update_factor(
        &mut self,
        column: &mut HVector,
        row_ep: &mut HVector,
        i_row: &mut HighsInt,
        hint: &mut HighsInt,
    ) {
        self.analysis_.simplex_timer_start(UpdateFactorClock, 0);
        self.simplex_nla_
            .update(&mut self.basis_.basic_index_, column, row_ep, i_row, hint);
        // Now have a representation of B^{-1}, but it is not fresh
        self.status_.has_invert = true;
        if self.info_.update_count >= self.info_.update_limit {
            *hint = kRebuildReasonUpdateLimitReached;
        }

        // Determine whether to reinvert based on the synthetic clock
        let reinvert_synthetic_clock = self.total_synthetic_tick_ >= self.build_synthetic_tick_;
        let performed_min_updates =
            self.info_.update_count >= kSyntheticTickReinversionMinUpdateCount;
        if reinvert_synthetic_clock && performed_min_updates {
            *hint = kRebuildReasonSyntheticClockSaysInvert;
        }
        self.analysis_.simplex_timer_stop(UpdateFactorClock, 0);
        // Use the next level down for the debug level, since the cost of
        // checking the INVERT every iteration is an order more expensive
        // than checking after factorization.
        let alt_debug_level = self.options_.highs_debug_level - 1;
        let debug_status = self.debug_nla_check_invert("HEkk::updateFactor", alt_debug_level);
        if debug_status == HighsDebugStatus::kError {
            *hint = kRebuildReasonPossiblySingularBasis;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::updateFactor (lines 3075-3103), called by
    /// `HEkkDual::majorUpdateFactor` (HEkkDualMulti.cpp 906-933) with
    /// `multi_finish[0].col_aq` / `row_ep` chained through `HVector::next` and the `iRows`
    /// array. Same body as [`HEkk::update_factor`] with
    /// `simplex_nla_.update_multi`.
    ///
    /// # Safety
    /// As [`HFactor::update_cft`].
    pub unsafe fn update_factor_multi(
        &mut self,
        column: *mut HVector,
        row_ep: *mut HVector,
        i_row: &mut [HighsInt],
        hint: &mut HighsInt,
    ) {
        self.analysis_.simplex_timer_start(UpdateFactorClock, 0);
        // SAFETY: caller guarantees column/row_ep are valid HVector pointers (the
        // multi_finish chain) and i_row is a valid slice of HighsInt for the chain length.
        // Mirrors C++ HEkkDualMulti.cpp 906-933 passing raw pointers to updateFactor.
        unsafe {
            self.simplex_nla_.update_multi(
                &mut self.basis_.basic_index_,
                column,
                row_ep,
                i_row,
                hint,
            );
        }
        // Now have a representation of B^{-1}, but it is not fresh
        self.status_.has_invert = true;
        if self.info_.update_count >= self.info_.update_limit {
            *hint = kRebuildReasonUpdateLimitReached;
        }

        // Determine whether to reinvert based on the synthetic clock
        let reinvert_synthetic_clock = self.total_synthetic_tick_ >= self.build_synthetic_tick_;
        let performed_min_updates =
            self.info_.update_count >= kSyntheticTickReinversionMinUpdateCount;
        if reinvert_synthetic_clock && performed_min_updates {
            *hint = kRebuildReasonSyntheticClockSaysInvert;
        }
        self.analysis_.simplex_timer_stop(UpdateFactorClock, 0);
        // Use the next level down for the debug level, since the cost of
        // checking the INVERT every iteration is an order more expensive
        // than checking after factorization.
        let alt_debug_level = self.options_.highs_debug_level - 1;
        let debug_status = self.debug_nla_check_invert("HEkk::updateFactor", alt_debug_level);
        if debug_status == HighsDebugStatus::kError {
            *hint = kRebuildReasonPossiblySingularBasis;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::updatePivots (lines 3105-3150, 46 lines)
    pub fn update_pivots(&mut self, variable_in: HighsInt, row_out: HighsInt, move_out: HighsInt) {
        self.analysis_.simplex_timer_start(UpdatePivotsClock, 0);
        let variable_out = self.basis_.basic_index_[row_out as usize];

        // update hash value of basis
        sparse_inverse_combine(&mut self.basis_.hash, variable_out);
        sparse_combine(&mut self.basis_.hash, variable_in);
        self.visited_basis_.insert(self.basis_.hash);

        let vin = variable_in as usize;
        let vout = variable_out as usize;
        // Incoming variable
        self.basis_.basic_index_[row_out as usize] = variable_in;
        self.basis_.nonbasic_flag_[vin] = 0;
        self.basis_.nonbasic_move_[vin] = 0;
        self.info_.base_lower_[row_out as usize] = self.info_.work_lower_[vin];
        self.info_.base_upper_[row_out as usize] = self.info_.work_upper_[vin];

        // Outgoing variable
        self.basis_.nonbasic_flag_[vout] = 1;
        if self.info_.work_lower_[vout] == self.info_.work_upper_[vout] {
            self.info_.work_value_[vout] = self.info_.work_lower_[vout];
            self.basis_.nonbasic_move_[vout] = 0;
        } else if move_out == -1 {
            self.info_.work_value_[vout] = self.info_.work_lower_[vout];
            self.basis_.nonbasic_move_[vout] = 1;
        } else {
            self.info_.work_value_[vout] = self.info_.work_upper_[vout];
            self.basis_.nonbasic_move_[vout] = -1;
        }
        // Update the dual objective value
        let nw_value = self.info_.work_value_[vout];
        let vr_dual = self.info_.work_dual_[vout];
        let dl_dual_objective_value = nw_value * vr_dual;
        self.info_.updated_dual_objective_value += dl_dual_objective_value;
        self.info_.update_count += 1;
        // Update the number of basic logicals
        if variable_out < self.lp_.num_col_ {
            self.info_.num_basic_logicals += 1;
        }
        if variable_in < self.lp_.num_col_ {
            self.info_.num_basic_logicals -= 1;
        }
        // No longer have a representation of B^{-1}, and certainly not
        // fresh!
        self.status_.has_invert = false;
        self.status_.has_fresh_invert = false;
        // Data are no longer fresh from rebuild
        self.status_.has_fresh_rebuild = false;
        self.analysis_.simplex_timer_stop(UpdatePivotsClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::isBadBasisChange (lines 3152-3203, 52 lines)
    pub fn is_bad_basis_change(
        &mut self,
        algorithm: SimplexAlgorithm,
        variable_in: HighsInt,
        row_out: HighsInt,
        rebuild_reason: HighsInt,
    ) -> bool {
        if rebuild_reason != 0 {
            return false;
        }
        if variable_in == -1 || row_out == -1 {
            return false;
        }
        let mut currhash: u64 = self.basis_.hash;
        let variable_out = self.basis_.basic_index_[row_out as usize];

        sparse_inverse_combine(&mut currhash, variable_out);
        sparse_combine(&mut currhash, variable_in);

        let mut cycling_detected = false;
        let posible_cycling = self.visited_basis_.contains(&currhash);
        if posible_cycling {
            if self.iteration_count_ == self.previous_iteration_cycling_detected + 1 {
                // Cycling detected on successive iterations suggests infinite cycling
                cycling_detected = true;
            } else {
                self.previous_iteration_cycling_detected = self.iteration_count_;
            }
        }
        if cycling_detected {
            if algorithm == SimplexAlgorithm::kDual {
                self.analysis_.num_dual_cycling_detections += 1;
            } else {
                self.analysis_.num_primal_cycling_detections += 1;
            }
            highs_log_dev(
                &self.options_.log_options,
                HighsLogType::kWarning,
                &format!(
                    " basis change ({} out; {} in) is bad\n",
                    variable_out, variable_in
                ),
            );
            self.add_bad_basis_change(
                row_out,
                variable_out,
                variable_in,
                BadBasisChangeReason::kCycling,
                true,
            );
            return true;
        } else {
            // Look to see whether this basis change is in the list of bad
            // ones
            for change in self.bad_basis_change_.iter_mut() {
                if change.variable_out == variable_out
                    && change.variable_in == variable_in
                    && change.row_out == row_out
                {
                    change.taboo = true;
                    return true;
                }
            }
        }

        false
    }

    /// HiGHS: HEkk.cpp::HEkk::updateMatrix (lines 3205-3211, 7 lines)
    pub fn update_matrix(&mut self, variable_in: HighsInt, variable_out: HighsInt) {
        self.analysis_.simplex_timer_start(UpdateMatrixClock, 0);
        self.ar_matrix_
            .update(variable_in, variable_out, &self.lp_.a_matrix_);
        self.analysis_.simplex_timer_stop(UpdateMatrixClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::computeInfeasibilitiesForReporting (lines 3213-3229, 17 lines)
    /// C++ default arguments: solve_phase = kSolvePhase2
    pub fn compute_infeasibilities_for_reporting(
        &mut self,
        algorithm: SimplexAlgorithm,
        solve_phase: HighsInt,
    ) {
        if algorithm == SimplexAlgorithm::kPrimal {
            // Report the primal and dual infeasibilities
            self.compute_simplex_infeasible();
        } else {
            // Report the primal infeasibilities
            self.compute_simplex_primal_infeasible();
            if solve_phase == kSolvePhase1 {
                // In phase 1, report the simplex LP dual infeasibilities
                self.compute_simplex_lp_dual_infeasible();
            } else {
                // In phase 2, report the simplex dual infeasibilities
                self.compute_simplex_dual_infeasible();
            }
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::computeSimplexInfeasible (lines 3231-3234, 4 lines)
    pub fn compute_simplex_infeasible(&mut self) {
        self.compute_simplex_primal_infeasible();
        self.compute_simplex_dual_infeasible();
    }

    /// HiGHS: HEkk.cpp::HEkk::computeSimplexPrimalInfeasible (lines 3236-3294, 59 lines)
    pub fn compute_simplex_primal_infeasible(&mut self) {
        // Computes num/max/sum of primal infeasibilities according to the
        // simplex bounds. This is used to determine optimality in dual
        // phase 1 and dual phase 2, albeit using different bounds in
        // workLower/Upper.
        self.analysis_.simplex_timer_start(ComputePrIfsClock, 0);
        let scaled_primal_feasibility_tolerance = self.options_.primal_feasibility_tolerance;
        let mut num_primal_infeasibility: HighsInt = 0;
        let mut max_primal_infeasibility: f64 = 0.0;
        let mut sum_primal_infeasibility: f64 = 0.0;

        for i in 0..(self.lp_.num_col_ + self.lp_.num_row_) as usize {
            if self.basis_.nonbasic_flag_[i] != 0 {
                // Nonbasic column
                let value = self.info_.work_value_[i];
                let lower = self.info_.work_lower_[i];
                let upper = self.info_.work_upper_[i];
                // @primal_infeasibility calculation
                let mut primal_infeasibility = 0.0;
                if value < lower - scaled_primal_feasibility_tolerance {
                    primal_infeasibility = lower - value;
                } else if value > upper + scaled_primal_feasibility_tolerance {
                    primal_infeasibility = value - upper;
                }
                if primal_infeasibility > 0.0 {
                    if primal_infeasibility > scaled_primal_feasibility_tolerance {
                        num_primal_infeasibility += 1;
                    }
                    max_primal_infeasibility =
                        cpp_max(primal_infeasibility, max_primal_infeasibility);
                    sum_primal_infeasibility += primal_infeasibility;
                }
            }
        }
        for i in 0..self.lp_.num_row_ as usize {
            // Basic variable
            let value = self.info_.base_value_[i];
            let lower = self.info_.base_lower_[i];
            let upper = self.info_.base_upper_[i];
            // @primal_infeasibility calculation
            let mut primal_infeasibility = 0.0;
            if value < lower - scaled_primal_feasibility_tolerance {
                primal_infeasibility = lower - value;
            } else if value > upper + scaled_primal_feasibility_tolerance {
                primal_infeasibility = value - upper;
            }
            if primal_infeasibility > 0.0 {
                if primal_infeasibility > scaled_primal_feasibility_tolerance {
                    num_primal_infeasibility += 1;
                }
                max_primal_infeasibility = cpp_max(primal_infeasibility, max_primal_infeasibility);
                sum_primal_infeasibility += primal_infeasibility;
            }
        }
        self.info_.num_primal_infeasibilities = num_primal_infeasibility;
        self.info_.max_primal_infeasibility = max_primal_infeasibility;
        self.info_.sum_primal_infeasibilities = sum_primal_infeasibility;
        self.analysis_.simplex_timer_stop(ComputePrIfsClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::computeSimplexDualInfeasible (lines 3296-3336, 41 lines)
    pub fn compute_simplex_dual_infeasible(&mut self) {
        self.analysis_.simplex_timer_start(ComputeDuIfsClock, 0);
        // Computes num/max/sum of dual infeasibilities in phase 1 and phase
        // 2 according to nonbasicMove. The bounds are only used to identify
        // free variables. Fixed variables are assumed to have
        // nonbasicMove=0 so that no dual infeasibility is counted for them.
        let scaled_dual_feasibility_tolerance = self.options_.dual_feasibility_tolerance;
        let mut num_dual_infeasibility: HighsInt = 0;
        let mut max_dual_infeasibility: f64 = 0.0;
        let mut sum_dual_infeasibility: f64 = 0.0;

        for i_col in 0..(self.lp_.num_col_ + self.lp_.num_row_) as usize {
            if self.basis_.nonbasic_flag_[i_col] == 0 {
                continue;
            }
            // Nonbasic column
            let dual = self.info_.work_dual_[i_col];
            let lower = self.info_.work_lower_[i_col];
            let upper = self.info_.work_upper_[i_col];
            let dual_infeasibility = if highs_is_infinity(-lower) && highs_is_infinity(upper) {
                // Free: any nonzero dual value is infeasible
                dual.abs()
            } else {
                // Not free: any dual infeasibility is given by the dual value
                // signed by nonbasicMove
                (-(self.basis_.nonbasic_move_[i_col] as HighsInt)) as f64 * dual
            };
            if dual_infeasibility > 0.0 {
                if dual_infeasibility >= scaled_dual_feasibility_tolerance {
                    num_dual_infeasibility += 1;
                }
                max_dual_infeasibility = cpp_max(dual_infeasibility, max_dual_infeasibility);
                sum_dual_infeasibility += dual_infeasibility;
            }
        }
        self.info_.num_dual_infeasibilities = num_dual_infeasibility;
        self.info_.max_dual_infeasibility = max_dual_infeasibility;
        self.info_.sum_dual_infeasibilities = sum_dual_infeasibility;
        self.analysis_.simplex_timer_stop(ComputeDuIfsClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::computeSimplexLpDualInfeasible (lines 3338-3422, 85 lines)
    pub fn compute_simplex_lp_dual_infeasible(&mut self) {
        // Compute num/max/sum of dual infeasibilities according to the
        // bounds of the simplex LP. Assumes that boxed variables have
        // primal variable at the bound corresponding to the sign of the
        // dual so should only be used in dual phase 1 - where it's only
        // used for reporting after rebuilds.
        let scaled_dual_feasibility_tolerance = self.options_.dual_feasibility_tolerance;
        let mut num_dual_infeasibility: HighsInt = 0;
        let mut max_dual_infeasibility: f64 = 0.0;
        let mut sum_dual_infeasibility: f64 = 0.0;

        for i_col in 0..self.lp_.num_col_ {
            let i_var = i_col;
            if self.basis_.nonbasic_flag_[i_var as usize] == 0 {
                continue;
            }
            // Nonbasic column
            let dual = self.info_.work_dual_[i_var as usize];
            let lower = self.lp_.col_lower_[i_col as usize];
            let upper = self.lp_.col_upper_[i_col as usize];
            let dual_infeasibility = if highs_is_infinity(upper) {
                if highs_is_infinity(-lower) {
                    // Free: any nonzero dual value is infeasible
                    dual.abs()
                } else {
                    // Only lower bounded: a negative dual is infeasible
                    -dual
                }
            } else if highs_is_infinity(-lower) {
                // Only upper bounded: a positive dual is infeasible
                dual
            } else {
                // Boxed or fixed: any dual value is feasible
                0.0
            };
            if dual_infeasibility > 0.0 {
                if dual_infeasibility >= scaled_dual_feasibility_tolerance {
                    num_dual_infeasibility += 1;
                }
                max_dual_infeasibility = cpp_max(dual_infeasibility, max_dual_infeasibility);
                sum_dual_infeasibility += dual_infeasibility;
            }
        }
        for i_row in 0..self.lp_.num_row_ {
            let i_var = self.lp_.num_col_ + i_row;
            if self.basis_.nonbasic_flag_[i_var as usize] == 0 {
                continue;
            }
            // Nonbasic row
            let dual = -self.info_.work_dual_[i_var as usize];
            let lower = self.lp_.row_lower_[i_row as usize];
            let upper = self.lp_.row_upper_[i_row as usize];
            let dual_infeasibility = if highs_is_infinity(upper) {
                if highs_is_infinity(-lower) {
                    // Free: any nonzero dual value is infeasible
                    dual.abs()
                } else {
                    // Only lower bounded: a negative dual is infeasible
                    -dual
                }
            } else if highs_is_infinity(-lower) {
                // Only upper bounded: a positive dual is infeasible
                dual
            } else {
                // Boxed or fixed: any dual value is feasible
                0.0
            };
            if dual_infeasibility > 0.0 {
                if dual_infeasibility >= scaled_dual_feasibility_tolerance {
                    num_dual_infeasibility += 1;
                }
                max_dual_infeasibility = cpp_max(dual_infeasibility, max_dual_infeasibility);
                sum_dual_infeasibility += dual_infeasibility;
            }
        }
        self.analysis_.num_dual_phase_1_lp_dual_infeasibility = num_dual_infeasibility;
        self.analysis_.max_dual_phase_1_lp_dual_infeasibility = max_dual_infeasibility;
        self.analysis_.sum_dual_phase_1_lp_dual_infeasibility = sum_dual_infeasibility;
    }

    /// HiGHS: HEkk.cpp::HEkk::invalidatePrimalMaxSumInfeasibilityRecord (lines 3424-3427, 4 lines)
    pub fn invalidate_primal_max_sum_infeasibility_record(&mut self) {
        self.info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.info_.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
    }

    /// HiGHS: HEkk.cpp::HEkk::invalidatePrimalInfeasibilityRecord (lines 3429-3432, 4 lines)
    pub fn invalidate_primal_infeasibility_record(&mut self) {
        self.info_.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.invalidate_primal_max_sum_infeasibility_record();
    }

    /// HiGHS: HEkk.cpp::HEkk::invalidateDualMaxSumInfeasibilityRecord (lines 3434-3437, 4 lines)
    pub fn invalidate_dual_max_sum_infeasibility_record(&mut self) {
        self.info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
        self.info_.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
    }

    /// HiGHS: HEkk.cpp::HEkk::invalidateDualInfeasibilityRecord (lines 3439-3442, 4 lines)
    pub fn invalidate_dual_infeasibility_record(&mut self) {
        self.info_.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
        self.invalidate_dual_max_sum_infeasibility_record();
    }
}
