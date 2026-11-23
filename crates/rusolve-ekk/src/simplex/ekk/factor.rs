//! HiGHS source: simplex/HEkk.cpp (lines 1695-2325)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Strategy/threads choice, (backtracking) INVERT, objective values, DSE/Devex weights.

use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkk.cpp::HEkk::chooseSimplexStrategyThreads (lines 1695-1785, 91 lines)
    // PORT NOTE: the C++ member reads no other member, so this takes no `self`; the caller
    // can then pass `&self.options_, &mut self.info_` without copying either.
    pub fn choose_simplex_strategy_threads(options: &HighsOptions, info: &mut HighsSimplexInfo) {
        debug_assert!(info.num_dual_infeasibilities > 0 || info.num_primal_infeasibilities > 0);
        // Set the internal simplex strategy and number of threads for dual
        // simplex
        // By default, use the HighsOptions strategy. If this is
        // kSimplexStrategyChoose, then the strategy used will depend on
        // whether the current basis is primal feasible.
        info.simplex_strategy = options.simplex_strategy;
        if info.simplex_strategy == kSimplexStrategyChoose {
            // HiGHS is left to choose the simplex strategy
            if info.num_primal_infeasibilities > 0 {
                // Not primal feasible, so use dual simplex
                info.simplex_strategy = kSimplexStrategyDual;
            } else {
                // Primal feasible. so use primal simplex
                info.simplex_strategy = kSimplexStrategyPrimal;
            }
        }
        // Set min/max_threads to correspond to serial code. They will be
        // set to other values if parallel options are used.
        info.min_concurrency = 1;
        info.max_concurrency = 1;
        // Record the min/max minimum concurrency in the options
        let simplex_min_concurrency = options.simplex_min_concurrency;
        let simplex_max_concurrency = options.simplex_max_concurrency;
        // C++ HEkk.cpp 1724: HighsInt max_threads = highs::parallel::num_threads();
        let max_threads: HighsInt = crate::parallel::highs_parallel::num_threads();

        // Don't allow the parallel option to switch from
        // kSimplexStrategyDual to kSimplexStrategyDualMulti so the MIP
        // solver only ever uses serial simplex

        if info.simplex_strategy == kSimplexStrategyDualTasks {
            info.min_concurrency = cpp_max(kDualTasksMinConcurrency, simplex_min_concurrency);
            info.max_concurrency = cpp_max(info.min_concurrency, simplex_max_concurrency);
        } else if info.simplex_strategy == kSimplexStrategyDualMulti {
            info.min_concurrency = cpp_max(kDualMultiMinConcurrency, simplex_min_concurrency);
            info.max_concurrency = cpp_max(info.min_concurrency, simplex_max_concurrency);
        }

        // Set the concurrency to be used to be the maximum number
        info.num_concurrency = info.max_concurrency;
        // Give a warning if the concurrency to be used is less than the
        // minimum concurrency allowed
        if info.num_concurrency < simplex_min_concurrency {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kWarning,
                "Using concurrency of {} for parallel strategy rather than minimum number ({}) specified in options\n",
                info.num_concurrency,
                simplex_min_concurrency
            );
        }
        // Give a warning if the concurrency to be used is more than the
        // maximum concurrency allowed
        if info.num_concurrency > simplex_max_concurrency {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kWarning,
                "Using concurrency of {} for parallel strategy rather than maximum number ({}) specified in options\n",
                info.num_concurrency,
                simplex_max_concurrency
            );
        }
        // Give a warning if the concurrency to be used is less than the
        // number of threads available
        if info.num_concurrency > max_threads {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kWarning,
                "Number of threads available = {} < {} = Simplex concurrency to be used: Parallel performance may be less than anticipated\n",
                max_threads,
                info.num_concurrency
            );
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::getNonsingularInverse (lines 1787-1880, 94 lines)
    /// C++ default arguments: solve_phase = 0
    pub fn get_nonsingular_inverse(&mut self, solve_phase: HighsInt) -> bool {
        debug_assert!(self.status_.has_basis);
        // Take a copy of basicIndex from before INVERT to be used as the
        // saved ordering of basic variables - so reinvert will run
        // identically.
        let basic_index_before_compute_factor: Vec<HighsInt> = self.basis_.basic_index_.clone();
        // Save the number of updates performed in case it has to be used to determine
        // a limit
        let simplex_update_count = self.info_.update_count;
        // Dual simplex edge weights are identified with rows, so must be
        // permuted according to INVERT. Scatter the edge weights so that,
        // after INVERT, they can be gathered according to the new
        // permutation of basicIndex
        self.analysis_.simplex_timer_start(PermWtClock, 0);
        for i in 0..self.lp_.num_row_ {
            self.scattered_dual_edge_weight_[self.basis_.basic_index_[i as usize] as usize] =
                self.dual_edge_weight_[i as usize];
        }
        self.analysis_.simplex_timer_stop(PermWtClock, 0);

        // Call computeFactor to perform INVERT
        let mut rank_deficiency = self.compute_factor();
        if rank_deficiency != 0 {
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "HEkk::getNonsingularInverse Rank_deficiency: solve {} (Iteration {})\n",
                self.debug_solve_call_num_,
                self.iteration_count_
            );
        }

        let artificial_rank_deficiency = false;
        if artificial_rank_deficiency {
            if self.info_.phase1_backtracking_test_done == 0 && solve_phase == kSolvePhase1 {
                // Claim rank deficiency to test backtracking
                rank_deficiency = 1;
                self.info_.phase1_backtracking_test_done = 1;
            } else if self.info_.phase2_backtracking_test_done == 0 && solve_phase == kSolvePhase2 {
                // Claim rank deficiency to test backtracking
                rank_deficiency = 1;
                self.info_.phase2_backtracking_test_done = 1;
            }
        }
        if rank_deficiency != 0 {
            // Rank deficient basis, so backtrack to last full rank basis
            //
            // Get the last nonsingular basis - so long as there is one
            let deficient_hash: u64 = self.basis_.hash;
            if !self.get_backtracking_basis() {
                return false;
            }
            // Record that backtracking is taking place
            self.info_.backtracking_ = true;
            self.visited_basis_.clear();
            self.visited_basis_.insert(self.basis_.hash);
            self.visited_basis_.insert(deficient_hash);
            self.update_status(LpAction::kBacktracking);
            let backtrack_rank_deficiency = self.compute_factor();
            // This basis has previously been inverted successfully, so it shouldn't be
            // singular
            if backtrack_rank_deficiency != 0 {
                return false;
            }
            // simplex update limit will be half of the number of updates
            // performed, so make sure that at least one update was performed
            if simplex_update_count <= 1 {
                return false;
            }
            let use_simplex_update_limit = self.info_.update_limit;
            let new_simplex_update_limit = simplex_update_count / 2;
            self.info_.update_limit = new_simplex_update_limit;
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "Rank deficiency of {} after {} simplex updates, so backtracking: max updates reduced from {} to {}\n",
                rank_deficiency,
                simplex_update_count,
                use_simplex_update_limit,
                new_simplex_update_limit
            );
        } else {
            // Current basis is full rank so save it
            self.put_backtracking_basis_from(&basic_index_before_compute_factor);
            // Indicate that backtracking is not taking place
            self.info_.backtracking_ = false;
            // Reset the update limit in case this is the first successful
            // inversion after backtracking
            self.info_.update_limit = self.options_.simplex_update_limit;
        }
        // Gather the edge weights according to the permutation of
        // basicIndex after INVERT
        self.analysis_.simplex_timer_start(PermWtClock, 0);
        for i in 0..self.lp_.num_row_ {
            self.dual_edge_weight_[i as usize] =
                self.scattered_dual_edge_weight_[self.basis_.basic_index_[i as usize] as usize];
        }
        self.analysis_.simplex_timer_stop(PermWtClock, 0);
        true
    }

    /// HiGHS: HEkk.cpp::HEkk::getBacktrackingBasis (lines 1882-1895, 14 lines)
    pub fn get_backtracking_basis(&mut self) -> bool {
        if !self.info_.valid_backtracking_basis_ {
            return false;
        }
        self.basis_ = self.info_.backtracking_basis_.clone();
        self.info_.costs_shifted = self.info_.backtracking_basis_costs_shifted_ != 0;
        self.info_.costs_perturbed = self.info_.backtracking_basis_costs_perturbed_ != 0;
        self.info_.bounds_shifted = self.info_.backtracking_basis_bounds_shifted_ != 0;
        self.info_.bounds_perturbed = self.info_.backtracking_basis_bounds_perturbed_ != 0;
        self.info_.work_shift_ = self.info_.backtracking_basis_work_shift_.clone();
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        for i_var in 0..num_tot {
            self.scattered_dual_edge_weight_[i_var as usize] =
                self.info_.backtracking_basis_edge_weight_[i_var as usize];
        }
        true
    }

    /// HiGHS: HEkk.cpp::HEkk::putBacktrackingBasis (lines 1897-1904, 8 lines)
    pub fn put_backtracking_basis(&mut self) {
        self.analysis_.simplex_timer_start(PermWtClock, 0);
        for i in 0..self.lp_.num_row_ {
            self.scattered_dual_edge_weight_[self.basis_.basic_index_[i as usize] as usize] =
                self.dual_edge_weight_[i as usize];
        }
        self.analysis_.simplex_timer_stop(PermWtClock, 0);
        let basic_index = self.basis_.basic_index_.clone();
        self.put_backtracking_basis_from(&basic_index);
    }

    /// HiGHS: HEkk.cpp::HEkk::putBacktrackingBasis (lines 1906-1920, 15 lines)
    pub fn put_backtracking_basis_from(&mut self, basic_index_before_compute_factor: &[HighsInt]) {
        self.info_.valid_backtracking_basis_ = true;
        self.info_.backtracking_basis_ = self.basis_.clone();
        self.info_.backtracking_basis_.basic_index_ = basic_index_before_compute_factor.to_vec();
        self.info_.backtracking_basis_costs_shifted_ = self.info_.costs_shifted as HighsInt;
        self.info_.backtracking_basis_costs_perturbed_ = self.info_.costs_perturbed as HighsInt;
        self.info_.backtracking_basis_bounds_shifted_ = self.info_.bounds_shifted as HighsInt;
        self.info_.backtracking_basis_bounds_perturbed_ = self.info_.bounds_perturbed as HighsInt;
        self.info_.backtracking_basis_work_shift_ = self.info_.work_shift_.clone();
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        for i_var in 0..num_tot {
            self.info_.backtracking_basis_edge_weight_[i_var as usize] =
                self.scattered_dual_edge_weight_[i_var as usize];
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::computePrimalObjectiveValue (lines 1922-1944, 23 lines)
    pub fn compute_primal_objective_value(&mut self) {
        self.analysis_.simplex_timer_start(ComputePrObjClock, 0);
        self.info_.primal_objective_value = 0.0;
        for i_row in 0..self.lp_.num_row_ {
            let i_var = self.basis_.basic_index_[i_row as usize];
            if i_var < self.lp_.num_col_ {
                self.info_.primal_objective_value +=
                    self.info_.base_value_[i_row as usize] * self.lp_.col_cost_[i_var as usize];
            }
        }
        for i_col in 0..self.lp_.num_col_ {
            if self.basis_.nonbasic_flag_[i_col as usize] != 0 {
                self.info_.primal_objective_value +=
                    self.info_.work_value_[i_col as usize] * self.lp_.col_cost_[i_col as usize];
            }
        }
        self.info_.primal_objective_value *= self.cost_scale_;
        // Objective value calculation is done using primal values and
        // original costs so offset is vanilla
        self.info_.primal_objective_value += self.lp_.offset_;
        // Now have primal objective value
        self.status_.has_primal_objective_value = true;
        self.analysis_.simplex_timer_stop(ComputePrObjClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::computeDualObjectiveValue (lines 1946-1967, 22 lines)
    /// C++ default arguments: phase = 2
    pub fn compute_dual_objective_value(&mut self, phase: HighsInt) {
        self.analysis_.simplex_timer_start(ComputeDuObjClock, 0);
        self.info_.dual_objective_value = 0.0;
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        for i_col in 0..num_tot {
            if self.basis_.nonbasic_flag_[i_col as usize] != 0 {
                self.info_.dual_objective_value +=
                    self.info_.work_value_[i_col as usize] * self.info_.work_dual_[i_col as usize];
            }
        }
        self.info_.dual_objective_value *= self.cost_scale_;
        if phase != 1 {
            // In phase 1 the dual objective has no objective
            // shift. Otherwise, if minimizing the shift is added. If
            // maximizing, workCost (and hence workDual) are negated, so the
            // shift is subtracted. Hence the shift is added according to the
            // sign implied by sense_
            self.info_.dual_objective_value +=
                ((self.lp_.sense_ as HighsInt) as f64) * self.lp_.offset_;
        }
        // Now have dual objective value
        self.status_.has_dual_objective_value = true;
        self.analysis_.simplex_timer_stop(ComputeDuObjClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::rebuildRefactor (lines 1969-2010, 42 lines)
    pub fn rebuild_refactor(&mut self, rebuild_reason: HighsInt) -> bool {
        // If no updates have been performed, then don't refactor!
        if self.info_.update_count == 0 {
            return false;
        }
        // Otherwise, refactor by default
        let mut refactor = true;
        let mut solution_error: f64 = 0.0;
        if self.options_.no_unnecessary_rebuild_refactor {
            // Consider whether not to refactor in rebuild
            //
            // Must have an INVERT just to consider this!
            debug_assert!(self.status_.has_invert);
            if rebuild_reason == kRebuildReasonNo
                || rebuild_reason == kRebuildReasonPossiblyOptimal
                || rebuild_reason == kRebuildReasonPossiblyPhase1Feasible
                || rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded
                || rebuild_reason == kRebuildReasonPossiblyDualUnbounded
                || rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex
            {
                // By default, don't refactor!
                refactor = false;
                // Possibly revise the decision based on accuracy when solving a
                // test system
                let error_tolerance = self.options_.rebuild_refactor_solution_error_tolerance;
                if error_tolerance > 0.0 {
                    solution_error = self.factor_solve_error();
                    refactor = solution_error > error_tolerance;
                }
            }
        }
        let _ = solution_error;
        // report_refactorization is the constant false in HiGHS; the printf block is omitted.
        refactor
    }

    /// HiGHS: HEkk.cpp::HEkk::computeFactor (lines 2012-2054, 43 lines)
    /// `simplex_nla_.invert()` -> `self.simplex_nla_.invert(&mut self.basis_.basic_index_)`.
    pub fn compute_factor(&mut self) -> HighsInt {
        debug_assert!(self.status_.has_nla);
        if self.status_.has_fresh_invert {
            return 0;
        }
        // Clear any bad basis changes
        self.clear_bad_basis_change(BadBasisChangeReason::kAll);
        debug_assert!(self.lp_factor_row_compatible());
        // Perform INVERT
        self.analysis_.simplex_timer_start(InvertClock, 0);
        let rank_deficiency = self.simplex_nla_.invert(&mut self.basis_.basic_index_);
        self.analysis_.simplex_timer_stop(InvertClock, 0);
        //
        // Set up hot start information
        self.hot_start_.refactor_info = self.simplex_nla_.factor_.refactor_info_.clone();
        self.hot_start_.nonbasic_move = self.basis_.nonbasic_move_.clone();
        self.hot_start_.valid = true;

        if self.analysis_.analyse_factor_data {
            self.analysis_
                .update_invert_form_data(&self.simplex_nla_.factor_);
        }

        let mut alt_debug_level: HighsInt = -1;
        if rank_deficiency != 0 {
            alt_debug_level = kHighsDebugLevelCostly;
        }
        self.debug_nla_check_invert("HEkk::computeFactor - original", alt_debug_level);

        if rank_deficiency != 0 {
            // Have an invertible representation, but of B with column(s)
            // replacements due to singularity. So no (fresh) representation of
            // B^{-1}
            self.status_.has_invert = false;
            self.status_.has_fresh_invert = false;
        } else {
            // Now have a representation of B^{-1}, and it is fresh!
            self.status_.has_invert = true;
            self.status_.has_fresh_invert = true;
        }
        // Set the update count to zero since the corrected invertible
        // representation may be used for an initial basis. In any case the
        // number of updates shouldn't be positive
        self.info_.update_count = 0;

        self.simplex_stats_.num_invert += 1;
        rank_deficiency
    }

    /// HiGHS: HEkk.cpp::HEkk::computeDualSteepestEdgeWeights (lines 2056-2077, 22 lines)
    /// C++ default arguments: initial = false
    pub fn compute_dual_steepest_edge_weights(&mut self, initial: bool) {
        if self.analysis_.analyse_simplex_time {
            self.analysis_.simplex_timer_start(SimplexIzDseWtClock, 0);
            self.analysis_.simplex_timer_start(DseIzClock, 0);
        }
        let num_row = self.lp_.num_row_;
        let mut row_ep = HVector::default();
        row_ep.setup(num_row);
        debug_assert!(self.dual_edge_weight_.len() as HighsInt >= num_row);
        for i_row in 0..num_row {
            self.dual_edge_weight_[i_row as usize] =
                self.compute_dual_steepest_edge_weight(i_row, &mut row_ep);
        }
        if self.analysis_.analyse_simplex_time {
            self.analysis_.simplex_timer_stop(SimplexIzDseWtClock, 0);
            self.analysis_.simplex_timer_stop(DseIzClock, 0);
            if initial {
                let iz_dse_wt_tt = self.analysis_.simplex_timer_read(SimplexIzDseWtClock, 0);
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kDetailed,
                    "Computed {} initial DSE weights in {}s\n",
                    num_row,
                    iz_dse_wt_tt
                );
            }
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::computeDualSteepestEdgeWeight (lines 2079-2091, 13 lines)
    pub fn compute_dual_steepest_edge_weight(
        &mut self,
        i_row: HighsInt,
        row_ep: &mut HVector,
    ) -> f64 {
        row_ep.clear();
        row_ep.count = 1;
        row_ep.index[0] = i_row;
        row_ep.array[i_row as usize] = 1.0;
        row_ep.pack_flag = false;
        self.simplex_nla_
            .btran_in_scaled_space(row_ep, self.info_.row_ep_density, None);
        let local_row_ep_density = (1.0 * row_ep.count as f64) / self.lp_.num_row_ as f64;
        let mut row_ep_density = self.info_.row_ep_density;
        self.update_operation_result_density(local_row_ep_density, &mut row_ep_density);
        self.info_.row_ep_density = row_ep_density;
        row_ep.norm2()
    }

    /// HiGHS: HEkk.cpp::HEkk::updateDualSteepestEdgeWeights (lines 2094-2238, 145 lines)
    pub fn update_dual_steepest_edge_weights(
        &mut self,
        row_out: HighsInt,
        variable_in: HighsInt,
        column: &HVector,
        new_pivotal_edge_weight: f64,
        kai: f64,
        dual_steepest_edge_array: &[f64],
    ) {
        self.analysis_.simplex_timer_start(DseUpdateWeightClock, 0);

        let num_row = self.lp_.num_row_;
        let column_count = column.count;
        let variable_index = &column.index;
        let column_array = &column.array;

        let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
        let col_aq_scale = self.simplex_nla_.variable_scale_factor(&refs, variable_in);
        let col_ap_scale = self.simplex_nla_.basic_col_scale_factor(&refs, row_out);
        let inv_col_ap_scale = 1.0 / col_ap_scale;

        let dse_check = false;
        let mut alt_dual_steepest_edge_column = HVector::default();
        let mut alt_pivotal_column = HVector::default();
        if dse_check {
            // Compute the DSE column otherwise to check
            alt_dual_steepest_edge_column.setup(num_row);
            alt_dual_steepest_edge_column.clear();
            alt_dual_steepest_edge_column.count = 1;
            alt_dual_steepest_edge_column.index[0] = row_out;
            alt_dual_steepest_edge_column.array[row_out as usize] = 1.0;
            alt_dual_steepest_edge_column.pack_flag = false;
            self.simplex_nla_.btran_in_scaled_space(
                &mut alt_dual_steepest_edge_column,
                self.info_.row_ep_density,
                None,
            );
            self.simplex_nla_.ftran_in_scaled_space(
                &mut alt_dual_steepest_edge_column,
                self.info_.row_dse_density,
                None,
            );
            // Compute the pivotal column in the scaled space otherwise to check
            //
            // Need \bar{B}^{-1}(R.aq.cq) = \bar{B}^{-1}R.(cq.aq)
            //
            alt_pivotal_column.setup(num_row);
            alt_pivotal_column.clear();
            //
            // Determine cq, and apply it in forming RHS
            //
            self.lp_
                .a_matrix_
                .collect_aj(&mut alt_pivotal_column, variable_in, col_aq_scale);
            self.simplex_nla_
                .apply_basis_matrix_row_scale(&refs, &mut alt_pivotal_column);
            self.simplex_nla_.ftran_in_scaled_space(
                &mut alt_pivotal_column,
                self.info_.col_aq_density,
                None,
            );
            let mut max_dse_column_error: f64 = 0.0;
            let mut sum_dse_column_error: f64 = 0.0;
            let mut num_dse_column_error: HighsInt = 0;
            let dse_column_value_tolerance = 1e-2;
            let dse_column_error_tolerance = 1e-4;
            let mut dse_array_count: HighsInt = 0;
            for i_row in 0..num_row {
                let dual_steepest_edge_array_value =
                    dual_steepest_edge_array[i_row as usize] * inv_col_ap_scale;
                if dual_steepest_edge_array_value != 0.0 {
                    dse_array_count += 1;
                }
                if dual_steepest_edge_array_value.abs() > dse_column_value_tolerance
                    || alt_dual_steepest_edge_column.array[i_row as usize].abs()
                        > dse_column_value_tolerance
                {
                    let dse_column_error = (alt_dual_steepest_edge_column.array[i_row as usize]
                        - dual_steepest_edge_array_value)
                        .abs()
                        / cpp_max(1.0, dual_steepest_edge_array_value.abs());
                    sum_dse_column_error += dse_column_error;
                    if dse_column_error > dse_column_error_tolerance {
                        max_dse_column_error = cpp_max(dse_column_error, max_dse_column_error);
                        num_dse_column_error += 1;
                    }
                }
            }
            if max_dse_column_error > dse_column_error_tolerance {
                println!(
                    "HEkk::updateDualSteepestEdgeWeights: Iter {:2} has num / max / sum = {} / {} / {} DSE column errors exceeding = {}",
                    self.iteration_count_,
                    num_dse_column_error,
                    max_dse_column_error,
                    sum_dse_column_error,
                    dse_column_error_tolerance
                );
                println!(
                    "DSE column count alt = {}; og = {})",
                    alt_dual_steepest_edge_column.count, dse_array_count
                );
                for i_row in 0..num_row {
                    let dual_steepest_edge_array_value =
                        dual_steepest_edge_array[i_row as usize] * inv_col_ap_scale;
                    if alt_dual_steepest_edge_column.array[i_row as usize] != 0.0
                        && dual_steepest_edge_array_value != 0.0
                    {
                        let dse_column_error = (alt_dual_steepest_edge_column.array
                            [i_row as usize]
                            - dual_steepest_edge_array_value)
                            .abs()
                            / cpp_max(1.0, dual_steepest_edge_array_value.abs());
                        if dse_column_error > 1e-10 {
                            println!(
                                "Row {:4}: DSE column (alt = {:11.4e}; og = {:11.4e}) difference {:10.4e}",
                                i_row,
                                alt_dual_steepest_edge_column.array[i_row as usize],
                                dual_steepest_edge_array_value,
                                dse_column_error
                            );
                        }
                    }
                }
                debug_assert!(max_dse_column_error < dse_column_error_tolerance);
            }
        }

        if (self.dual_edge_weight_.len() as HighsInt) < num_row {
            println!(
                "HEkk::updateDualSteepestEdgeWeights solve {}: dual_edge_weight_.size() = {} < {}",
                self.debug_solve_call_num_,
                self.dual_edge_weight_.len(),
                num_row
            );
        }
        debug_assert!(self.dual_edge_weight_.len() as HighsInt >= num_row);
        let mut to_entry: HighsInt = 0;
        let use_row_indices =
            self.simplex_nla_
                .sparse_loop_style(column_count, num_row, &mut to_entry);
        let convert_to_scaled_space = !self.simplex_in_scaled_space_;
        for i_entry in 0..to_entry {
            let i_row = if use_row_indices {
                variable_index[i_entry as usize]
            } else {
                i_entry
            };
            let mut aa_i_row = column_array[i_row as usize];
            if aa_i_row == 0.0 {
                continue;
            }
            let mut dual_steepest_edge_array_value = dual_steepest_edge_array[i_row as usize];
            if convert_to_scaled_space {
                let basic_col_scale = self.simplex_nla_.basic_col_scale_factor(&refs, i_row);
                aa_i_row /= basic_col_scale;
                aa_i_row *= col_aq_scale;
                dual_steepest_edge_array_value *= inv_col_ap_scale;
            }
            if dse_check {
                let pivotal_column_error =
                    (aa_i_row - alt_pivotal_column.array[i_row as usize]).abs();
                if pivotal_column_error > 1e-4 {
                    println!(
                        "HEkk::updateDualSteepestEdgeWeights Row {:2} of pivotal column has error {:10.4e}",
                        i_row, pivotal_column_error
                    );
                }
                debug_assert!(pivotal_column_error < 1e-4);
            }
            // One bounds check for the weight entry, not three. Same two statements as C++.
            let weight = &mut self.dual_edge_weight_[i_row as usize];
            *weight += aa_i_row
                * (new_pivotal_edge_weight * aa_i_row + kai * dual_steepest_edge_array_value);
            *weight = cpp_max(kMinDualSteepestEdgeWeight, *weight);
        }
        self.analysis_.simplex_timer_stop(DseUpdateWeightClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::updateDualDevexWeights (lines 2241-2269, 29 lines)
    pub fn update_dual_devex_weights(&mut self, column: &HVector, new_pivotal_edge_weight: f64) {
        self.analysis_
            .simplex_timer_start(DevexUpdateWeightClock, 0);

        let num_row = self.lp_.num_row_;
        let column_count = column.count;
        let variable_index = &column.index;
        let column_array = &column.array;

        if (self.dual_edge_weight_.len() as HighsInt) < num_row {
            println!(
                "HEkk::updateDualDevexWeights solve {}: dual_edge_weight_.size() = {} < {}",
                self.debug_solve_call_num_,
                self.dual_edge_weight_.len(),
                num_row
            );
        }
        debug_assert!(self.dual_edge_weight_.len() as HighsInt >= num_row);
        let mut to_entry: HighsInt = 0;
        let use_row_indices =
            self.simplex_nla_
                .sparse_loop_style(column_count, num_row, &mut to_entry);
        for i_entry in 0..to_entry {
            let i_row = if use_row_indices {
                variable_index[i_entry as usize]
            } else {
                i_entry
            };
            let aa_i_row = column_array[i_row as usize];
            self.dual_edge_weight_[i_row as usize] = cpp_max(
                self.dual_edge_weight_[i_row as usize],
                new_pivotal_edge_weight * aa_i_row * aa_i_row,
            );
        }
        self.analysis_.simplex_timer_stop(DevexUpdateWeightClock, 0);
    }

    /// HiGHS: HEkk.cpp::HEkk::resetSyntheticClock (lines 2271-2274, 4 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion): copies `simplex_nla_.build_synthetic_tick_` and zeroes `total_synthetic_tick_`.
    pub fn reset_synthetic_clock(&mut self) {
        self.build_synthetic_tick_ = self.simplex_nla_.build_synthetic_tick_;
        // the INVERT just finished cost `build_synthetic_tick_`; reporting-only.
        crate::work::add_simplex_ticks(self.build_synthetic_tick_);
        self.total_synthetic_tick_ = 0.0;
    }

    /// HiGHS: HEkk.cpp::HEkk::initialisePartitionedRowwiseMatrix (lines 2276-2284, 9 lines)
    pub fn initialise_partitioned_rowwise_matrix(&mut self) {
        if self.status_.has_ar_matrix {
            return;
        }
        self.analysis_.simplex_timer_start(matrixSetupClock, 0);
        self.ar_matrix_
            .create_rowwise_partitioned(&self.lp_.a_matrix_, Some(&self.basis_.nonbasic_flag_));
        debug_assert!(
            self.ar_matrix_
                .debug_partition_ok(&self.basis_.nonbasic_flag_)
        );
        self.analysis_.simplex_timer_stop(matrixSetupClock, 0);
        self.status_.has_ar_matrix = true;
    }

    /// HiGHS: HEkk.cpp::HEkk::lpFactorRowCompatible (lines 2286-2288, 3 lines)
    pub fn lp_factor_row_compatible(&self) -> bool {
        self.lp_factor_row_compatible_num(self.lp_.num_row_)
    }

    /// HiGHS: HEkk.cpp::HEkk::lpFactorRowCompatible (lines 2290-2302, 13 lines)
    pub fn lp_factor_row_compatible_num(&self, expected_num_row: HighsInt) -> bool {
        // Check for LP-HFactor row compatibility
        let consistent_num_row = self.simplex_nla_.factor_.num_row == expected_num_row;
        if !consistent_num_row {
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kError,
                "HEkk::initialiseSimplexLpBasisAndFactor: LP({:6}, {:6}) has factor_num_row = {}\n",
                self.lp_.num_col_,
                expected_num_row,
                self.simplex_nla_.factor_.num_row
            );
        }
        consistent_num_row
    }

    /// HiGHS: HEkk.cpp::HEkk::simplexStrategyToString (lines 2304-2320, 17 lines)
    pub fn simplex_strategy_to_string(&self, simplex_strategy: HighsInt) -> String {
        debug_assert!(
            kSimplexStrategyMin <= simplex_strategy && simplex_strategy <= kSimplexStrategyMax
        );
        if simplex_strategy == kSimplexStrategyChoose {
            return "choose simplex solver".to_string();
        }
        if simplex_strategy == kSimplexStrategyDual {
            return "dual simplex solver".to_string();
        }
        if simplex_strategy == kSimplexStrategyDualPlain {
            return "serial dual simplex solver".to_string();
        }
        if simplex_strategy == kSimplexStrategyDualTasks {
            return "parallel dual simplex solver - SIP".to_string();
        }
        if simplex_strategy == kSimplexStrategyDualMulti {
            return "parallel dual simplex solver - PAMI".to_string();
        }
        if simplex_strategy == kSimplexStrategyPrimal {
            return "primal simplex solver".to_string();
        }
        "Unknown".to_string()
    }

    /// HiGHS: HEkk.cpp::HEkk::zeroBasicDuals (lines 2322-2325, 4 lines)
    pub fn zero_basic_duals(&mut self) {
        for i_row in 0..self.lp_.num_row_ {
            self.info_.work_dual_[self.basis_.basic_index_[i_row as usize] as usize] = 0.0;
        }
    }
}
