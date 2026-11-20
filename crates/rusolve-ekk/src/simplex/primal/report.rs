//! HiGHS source: simplex/HEkkPrimal.cpp (lines 2575-2984)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Iteration analysis/reporting, free-column set, perturbed equations, bound shifting,
//! primal ray, debug stubs.
//! See `HEkkPrimal` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkPrimal {
    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::iterationAnalysisData (lines 2575-2621, 47 lines)
    pub fn iteration_analysis_data(&mut self, ekk: &mut HEkk) {
        // Possibly compute the infeasibility data
        if ekk.analysis_.analyse_simplex_runtime_data {
            ekk.compute_infeasibilities_for_reporting(SimplexAlgorithm::kPrimal, kSolvePhase2);
        }
        let info = &ekk.info_;
        ekk.analysis_.simplex_strategy = kSimplexStrategyPrimal;
        ekk.analysis_.edge_weight_mode = self.edge_weight_mode;
        ekk.analysis_.solve_phase = self.solve_phase;
        ekk.analysis_.simplex_iteration_count = ekk.iteration_count_;
        ekk.analysis_.devex_iteration_count = self.num_devex_iterations_;
        ekk.analysis_.pivotal_row_index = self.row_out;
        ekk.analysis_.leaving_variable = self.variable_out;
        ekk.analysis_.entering_variable = self.variable_in;
        ekk.analysis_.rebuild_reason = self.rebuild_reason;
        ekk.analysis_.reduced_rhs_value = 0.0;
        ekk.analysis_.reduced_cost_value = 0.0;
        ekk.analysis_.edge_weight = 0.0;
        ekk.analysis_.primal_delta = 0.0;
        ekk.analysis_.primal_step = self.theta_primal;
        ekk.analysis_.dual_step = self.theta_dual;
        ekk.analysis_.pivot_value_from_column = self.alpha_col;
        ekk.analysis_.pivot_value_from_row = self.alpha_row;
        ekk.analysis_.numerical_trouble = self.numerical_trouble;
        ekk.analysis_.edge_weight_error = ekk.edge_weight_error_;
        ekk.analysis_.objective_value = info.updated_primal_objective_value;
        ekk.analysis_.num_primal_infeasibility = info.num_primal_infeasibilities;
        ekk.analysis_.num_dual_infeasibility = info.num_dual_infeasibilities;
        ekk.analysis_.sum_primal_infeasibility = info.sum_primal_infeasibilities;
        ekk.analysis_.sum_dual_infeasibility = info.sum_dual_infeasibilities;
        if self.edge_weight_mode == EdgeWeightMode::kDevex && self.num_devex_iterations_ == 0 {
            ekk.analysis_.num_devex_framework += 1;
        }
        ekk.analysis_.col_aq_density = info.col_aq_density;
        ekk.analysis_.row_ep_density = info.row_ep_density;
        ekk.analysis_.row_ap_density = info.row_ap_density;
        ekk.analysis_.row_DSE_density = info.row_dse_density;
        ekk.analysis_.col_steepest_edge_density = info.col_steepest_edge_density;
        ekk.analysis_.col_basic_feasibility_change_density =
            info.col_basic_feasibility_change_density;
        ekk.analysis_.row_basic_feasibility_change_density =
            info.row_basic_feasibility_change_density;
        ekk.analysis_.col_BFRT_density = info.col_bfrt_density;
        ekk.analysis_.primal_col_density = info.primal_col_density;
        ekk.analysis_.dual_col_density = info.dual_col_density;
        ekk.analysis_.num_costly_DSE_iteration = info.num_costly_dse_iteration;
        ekk.analysis_.costly_DSE_measure = info.costly_dse_measure;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::iterationAnalysis (lines 2623-2627, 5 lines)
    pub fn iteration_analysis(&mut self, ekk: &mut HEkk) {
        self.iteration_analysis_data(ekk);
        ekk.analysis_.iteration_report();
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.iteration_record();
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::localReportIterHeader (lines 2629-2631, 3 lines)
    pub fn local_report_iter_header(&mut self, ekk: &mut HEkk) {
        println!(" Iter ColIn Row_Out ColOut");
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::localReportIter (lines 2633-2685, 53 lines)
    /// C++ default arguments: header = false
    pub fn local_report_iter(&mut self, ekk: &mut HEkk, header: bool) {
        if !self.report_hyper_chuzc {
            return;
        }
        let iteration_count = ekk.iteration_count_;
        if header {
            self.local_report_iter_header(ekk);
            self.last_header_iteration_count_ = iteration_count;
        } else {
            if ekk.iteration_count_ > self.last_header_iteration_count_ + 10 {
                self.local_report_iter_header(ekk);
                self.last_header_iteration_count_ = iteration_count;
            }
            if self.row_out >= 0 {
                print!(
                    "{:5} {:5}  {:5}  {:5}",
                    iteration_count, self.variable_in, self.row_out, self.variable_out
                );
            } else {
                print!("{:5} {:5} Bound flip   ", iteration_count, self.variable_in);
            }
            if self.check_column >= 0 && iteration_count >= self.check_iter {
                let flag = ekk.basis_.nonbasic_flag_[self.check_column as usize] as HighsInt;
                let move_ = ekk.basis_.nonbasic_move_[self.check_column as usize] as HighsInt;
                let info = &ekk.info_;
                let lower = info.work_lower_[self.check_column as usize];
                let upper = info.work_upper_[self.check_column as usize];
                let value: f64 = if flag == kNonbasicFlagTrue as HighsInt {
                    info.work_value_[self.check_column as usize]
                } else {
                    let mut i_row: HighsInt = 0;
                    while i_row < self.num_row {
                        if ekk.basis_.basic_index_[i_row as usize] == self.check_column {
                            break;
                        }
                        i_row += 1;
                    }
                    debug_assert!(i_row < self.num_row);
                    info.base_value_[i_row as usize]
                };
                print!(
                    ": Var {:2} ({:1}, {:2}) [{:9.4}, {:9.4}, {:9.4}]",
                    self.check_column, flag, move_, lower, value, upper
                );
                if flag == kNonbasicFlagTrue as HighsInt {
                    let dual = info.work_dual_[self.check_column as usize];
                    let weight = self.edge_weight_[self.check_column as usize];
                    let mut infeasibility = (-move_) as f64 * dual;
                    if lower == -kHighsInf && upper == kHighsInf {
                        infeasibility = dual.abs();
                    }
                    if infeasibility < self.dual_feasibility_tolerance {
                        infeasibility = 0.0;
                    }
                    let measure = infeasibility * infeasibility / weight;
                    print!(
                        " Du = {:9.4}; Wt = {:9.4}; Ms = {:9.4}",
                        dual, weight, measure
                    );
                }
            }
            println!();
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::reportRebuild (lines 2687-2695, 9 lines)
    /// C++ default arguments: reason_for_rebuild = -1
    pub fn report_rebuild(&mut self, ekk: &mut HEkk, reason_for_rebuild: HighsInt) {
        ekk.analysis_.simplex_timer_start(ReportRebuildClock, 0);
        self.iteration_analysis_data(ekk);
        ekk.analysis_.rebuild_reason = reason_for_rebuild;
        let rebuild_reason_string = ekk.rebuild_reason(reason_for_rebuild);
        ekk.analysis_.rebuild_reason_string = rebuild_reason_string;
        if ekk.options_.output_flag {
            ekk.analysis_.invert_report();
        }
        ekk.analysis_.simplex_timer_stop(ReportRebuildClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::getNonbasicFreeColumnSet (lines 2697-2710, 14 lines)
    pub fn get_nonbasic_free_column_set(&mut self, ekk: &mut HEkk) {
        if self.num_free_col == 0 {
            return;
        }
        debug_assert!(self.num_free_col > 0);
        let info = &ekk.info_;
        let basis = &ekk.basis_;
        self.nonbasic_free_col_set.clear();
        for i_col in 0..self.num_tot {
            let nonbasic_free = basis.nonbasic_flag_[i_col as usize] == kNonbasicFlagTrue
                && info.work_lower_[i_col as usize] <= -kHighsInf
                && info.work_upper_[i_col as usize] >= kHighsInf;
            if nonbasic_free {
                self.nonbasic_free_col_set.add(i_col);
            }
        }
        //  nonbasic_free_col_set.print();
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::removeNonbasicFreeColumn (lines 2712-2726, 15 lines)
    pub fn remove_nonbasic_free_column(&mut self, ekk: &mut HEkk) {
        let remove_nonbasic_free_column = ekk.basis_.nonbasic_move_[self.variable_in as usize] == 0;
        if remove_nonbasic_free_column {
            let removed_nonbasic_free_column = self.nonbasic_free_col_set.remove(self.variable_in);
            if !removed_nonbasic_free_column {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kError,
                    "HEkkPrimal::phase1update failed to remove nonbasic free column {}\n",
                    self.variable_in
                );
                debug_assert!(removed_nonbasic_free_column);
            }
        }
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::adjustPerturbedEquationOut (lines 2728-2757, 30 lines)
    pub fn adjust_perturbed_equation_out(&mut self, ekk: &mut HEkk) {
        if !ekk.info_.bounds_perturbed {
            return;
        }
        let lp = &ekk.lp_;
        let info = &mut ekk.info_;
        let (lp_lower, lp_upper) = if self.variable_out < self.num_col {
            (
                lp.col_lower_[self.variable_out as usize],
                lp.col_upper_[self.variable_out as usize],
            )
        } else {
            (
                -lp.row_upper_[(self.variable_out - self.num_col) as usize],
                -lp.row_lower_[(self.variable_out - self.num_col) as usize],
            )
        };
        if lp_lower < lp_upper {
            return;
        }
        // Leaving variable is fixed
        //  double save_theta_primal = theta_primal;
        let true_fixed_value = lp_lower;
        // Modify theta_primal so that variable leaves at true fixed value
        self.theta_primal =
            (info.base_value_[self.row_out as usize] - true_fixed_value) / self.alpha_col;
        /*
          printf("For equation %4" HIGHSINT_FORMAT " to be nonbasic at RHS %10.4g
          requires theta_primal to change by %10.4g from %10.4g to %10.4g\n",
          variable_out, true_fixed_value, theta_primal-save_theta_primal,
          save_theta_primal, theta_primal);
        */
        info.work_lower_[self.variable_out as usize] = true_fixed_value;
        info.work_upper_[self.variable_out as usize] = true_fixed_value;
        info.work_range_[self.variable_out as usize] = 0.0;
        self.value_in = info.work_value_[self.variable_in as usize] + self.theta_primal;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::getBasicPrimalInfeasibility (lines 2759-2808, 50 lines)
    pub fn get_basic_primal_infeasibility(&mut self, ekk: &mut HEkk) {
        // Gets the num/max/sum of basic primal infeasibilities,
        ekk.analysis_.simplex_timer_start(ComputePrIfsClock, 0);
        let primal_feasibility_tolerance = ekk.options_.primal_feasibility_tolerance;
        let info = &mut ekk.info_;
        let updated_num_primal_infeasibility = info.num_primal_infeasibilities;
        info.num_primal_infeasibilities = 0;
        info.max_primal_infeasibility = 0.0;
        info.sum_primal_infeasibilities = 0.0;

        for i_row in 0..self.num_row {
            let value = info.base_value_[i_row as usize];
            let lower = info.base_lower_[i_row as usize];
            let upper = info.base_upper_[i_row as usize];
            // @primal_infeasibility calculation
            let mut primal_infeasibility: f64 = 0.0;
            if value < lower - primal_feasibility_tolerance {
                primal_infeasibility = lower - value;
            } else if value > upper + primal_feasibility_tolerance {
                primal_infeasibility = value - upper;
            }
            if primal_infeasibility > 0.0 {
                if primal_infeasibility > primal_feasibility_tolerance {
                    info.num_primal_infeasibilities += 1;
                }
                info.max_primal_infeasibility =
                    cpp_max(primal_infeasibility, info.max_primal_infeasibility);
                info.sum_primal_infeasibilities += primal_infeasibility;
            }
        }
        if updated_num_primal_infeasibility >= 0 {
            // The number of primal infeasibilities should be correct
            let num_primal_infeasibility_ok =
                info.num_primal_infeasibilities == updated_num_primal_infeasibility;
            // if (!num_primal_infeasibility_ok)
            //   printf("In iteration %" HIGHSINT_FORMAT
            //          ": num_primal_infeasibility = %" HIGHSINT_FORMAT
            //          " != %" HIGHSINT_FORMAT
            //          " = "
            //          "updated_num_primal_infeasibility\n",
            //          ekk_instance_.iteration_count_, num_primal_infeasibility,
            //          updated_num_primal_infeasibility);

            debug_assert!(num_primal_infeasibility_ok);
        }
        ekk.analysis_.simplex_timer_stop(ComputePrIfsClock, 0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::shiftBound (lines 2810-2865, 56 lines)
    pub fn shift_bound(
        &mut self,
        ekk: &mut HEkk,
        lower: bool,
        i_var: HighsInt,
        value: f64,
        random_value: f64,
        bound: &mut f64,
        shift: &mut f64,
    ) {
        // If infeasibility is very large, then adding feasibility may not
        // yield a new value (see #1144) so new_infeasibility < 0 is false,
        // tripping the old assert
        //
        // Ambros proposed adding feasibility
        // *(infeasibility/scale_threshold) when
        // infeasibility/scale_threshold > 1, but this could lead to a large
        // value for new_infeasibility. Sounds better to accept degeneracy
        // in this edge case.
        // MUST-PORT (determinism): random_value is supplied by the caller
        // (numTotRandomValue_); shiftBound itself draws nothing from the stream.
        let feasibility = (1.0 + random_value) * self.primal_feasibility_tolerance;
        let old_bound = *bound;
        let type_name: &str = if lower { "lower" } else { "upper" };
        let (infeasibility, new_infeasibility) = if lower {
            // Bound to shift is lower
            debug_assert!(value < *bound - self.primal_feasibility_tolerance);
            let infeasibility = *bound - value;
            debug_assert!(infeasibility > 0.0);
            // Determine the amount by which value will be feasible - so that
            // (ideally) it's not degenerate
            *shift = infeasibility + feasibility;
            *bound -= *shift;
            let new_infeasibility = *bound - value;
            (infeasibility, new_infeasibility)
        } else {
            // Bound to shift is upper
            debug_assert!(value > *bound + self.primal_feasibility_tolerance);
            let infeasibility = value - *bound;
            debug_assert!(infeasibility > 0.0);
            // Determine the amount by which value will be feasible - so that
            // (ideally) it's not degenerate
            *shift = infeasibility + feasibility;
            *bound += *shift;
            let new_infeasibility = value - *bound;
            (infeasibility, new_infeasibility)
        };
        if new_infeasibility > 0.0 {
            // new_infeasibility should be non-positive, and negative unless
            // bound is excessively large, whereas feasibility is positive
            // C++ passes 12 arguments (the last is `error`) for 11 conversions;
            // printf consumes them in order, so the 11 printed values below
            // reproduce the C++ output, and `error` is never printed. The C++
            // fflush(stdout) is omitted (stdout flushing is not a solver effect).
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kInfo,
                "HEkkPrimal::shiftBound Value({:4}) = {:10.4} exceeds {}: \
                 random_value = {}; value = {}; \
                 feasibility = {}; infeasibility = {}; shift = {}; bound = {}; \
                 new_infeasibility = {} with error {}\n",
                i_var,
                value,
                type_name,
                old_bound,
                random_value,
                value,
                feasibility,
                infeasibility,
                *shift,
                *bound,
                new_infeasibility
            );
        }
        debug_assert!(new_infeasibility <= 0.0);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::savePrimalRay (lines 2867-2873, 7 lines)
    pub fn save_primal_ray(&mut self, ekk: &mut HEkk) {
        debug_assert!(self.variable_in >= 0);
        debug_assert!(self.move_in != kNoRaySign);
        ekk.primal_ray_record_.clear();
        ekk.primal_ray_record_.index = self.variable_in;
        ekk.primal_ray_record_.sign = -self.move_in;
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::debugPrimalSimplex (lines 2875-2885, 11 lines)
    /// C++ default arguments: initialise = false
    /// Debug check: out of scope; stub returns HighsDebugStatus::kNotChecked (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_primal_simplex(
        &mut self,
        ekk: &mut HEkk,
        message: &str,
        initialise: bool,
    ) -> HighsDebugStatus {
        HighsDebugStatus::kNotChecked
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::debugPrimalSteepestEdgeWeights (lines 2887-2902, 16 lines)
    /// Debug check: out of scope; stub returns HighsDebugStatus::kNotChecked (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_primal_steepest_edge_weights(
        &mut self,
        ekk: &mut HEkk,
        message: &str,
    ) -> HighsDebugStatus {
        HighsDebugStatus::kNotChecked
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::debugPrimalSteepestEdgeWeights (lines 2904-2979, 76 lines)
    /// C++ default arguments: alt_debug_level = -1
    /// Debug check: out of scope; stub returns HighsDebugStatus::kNotChecked (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_primal_steepest_edge_weights_level(
        &mut self,
        ekk: &mut HEkk,
        alt_debug_level: HighsInt,
    ) -> HighsDebugStatus {
        HighsDebugStatus::kNotChecked
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::isBadBasisChange (lines 2981-2984, 4 lines)
    pub fn is_bad_basis_change(&mut self, ekk: &mut HEkk) -> bool {
        ekk.is_bad_basis_change(
            SimplexAlgorithm::kPrimal,
            self.variable_in,
            self.row_out,
            self.rebuild_reason,
        )
    }
}
