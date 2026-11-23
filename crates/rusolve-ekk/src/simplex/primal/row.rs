//! HiGHS source: simplex/HEkkPrimal.cpp (lines 1138-1550)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! CHUZR: phase-1 and phase-2 ratio tests, bound swaps, pivot assessment, update.
//! See `HEkkPrimal` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkPrimal {
    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::phase1ChooseRow (lines 1138-1266, 129 lines)
    pub fn phase1_choose_row(&mut self, ekk: &mut HEkk) {
        // analysis->simplexTimerStart(Chuzr1Clock);
        // Collect phase 1 theta lists
        //

        let d_pivot_tol: f64 = if ekk.info_.update_count < 10 {
            1e-9
        } else if ekk.info_.update_count < 20 {
            1e-8
        } else {
            1e-7
        };
        self.ph1_sorter_r.clear();
        self.ph1_sorter_t.clear();
        for i in 0..self.col_aq.count {
            let i_row: HighsInt = self.col_aq.index[i as usize];
            let ru = i_row as usize;
            let d_alpha: f64 = self.col_aq.array[ru] * (self.move_in as f64);
            let base_lower = &ekk.info_.base_lower_;
            let base_upper = &ekk.info_.base_upper_;
            let base_value = &ekk.info_.base_value_;

            // When the basic variable x[i] decrease
            if d_alpha > d_pivot_tol {
                // Whether it can become feasible by going below its upper bound
                if base_value[ru] > base_upper[ru] + self.primal_feasibility_tolerance {
                    let d_feas_theta: f64 =
                        (base_value[ru] - base_upper[ru] - self.primal_feasibility_tolerance)
                            / d_alpha;
                    self.ph1_sorter_r.push((d_feas_theta, i_row));
                    self.ph1_sorter_t.push((d_feas_theta, i_row));
                }
                // Whether it can become infeasible (again) by going below its
                // lower bound
                if base_value[ru] > base_lower[ru] - self.primal_feasibility_tolerance
                    && base_lower[ru] > -kHighsInf
                {
                    let d_relax_theta: f64 = (base_value[ru] - base_lower[ru]
                        + self.primal_feasibility_tolerance)
                        / d_alpha;
                    let d_tight_theta: f64 = (base_value[ru] - base_lower[ru]) / d_alpha;
                    self.ph1_sorter_r
                        .push((d_relax_theta, i_row - self.num_row));
                    self.ph1_sorter_t
                        .push((d_tight_theta, i_row - self.num_row));
                }
            }

            // When the basic variable x[i] increase
            if d_alpha < -d_pivot_tol {
                // Whether it can become feasible by going above its lower bound
                if base_value[ru] < base_lower[ru] - self.primal_feasibility_tolerance {
                    let d_feas_theta: f64 = (base_value[ru] - base_lower[ru]
                        + self.primal_feasibility_tolerance)
                        / d_alpha;
                    self.ph1_sorter_r.push((d_feas_theta, i_row - self.num_row));
                    self.ph1_sorter_t.push((d_feas_theta, i_row - self.num_row));
                }
                // Whether it can become infeasible (again) by going above its
                // upper bound
                if base_value[ru] < base_upper[ru] + self.primal_feasibility_tolerance
                    && base_upper[ru] < kHighsInf
                {
                    let d_relax_theta: f64 =
                        (base_value[ru] - base_upper[ru] - self.primal_feasibility_tolerance)
                            / d_alpha;
                    let d_tight_theta: f64 = (base_value[ru] - base_upper[ru]) / d_alpha;
                    self.ph1_sorter_r.push((d_relax_theta, i_row));
                    self.ph1_sorter_t.push((d_tight_theta, i_row));
                }
            }
        }

        // analysis->simplexTimerStop(Chuzr1Clock);
        // When there are no candidates at all, we can leave it here
        if self.ph1_sorter_r.is_empty() {
            self.row_out = kNoRowChosen;
            self.variable_out = -1;
            return;
        }

        // Now sort the relaxed theta to find the final break point. TODO:
        // Consider partial sort. Or heapify [O(n)] and then pop k points
        // [kO(log(n))].

        // analysis->simplexTimerStart(Chuzr2Clock);
        pdqsort(&mut self.ph1_sorter_r[..]);
        let mut d_max_theta: f64 = self.ph1_sorter_r[0].0;
        let mut d_gradient: f64 = self.theta_dual.abs();
        for candidate in self.ph1_sorter_r.iter() {
            let d_my_theta: f64 = candidate.0;
            let index: HighsInt = candidate.1;
            let i_row: HighsInt = if index >= 0 {
                index
            } else {
                index + self.num_row
            };
            d_gradient -= self.col_aq.array[i_row as usize].abs();
            // Stop when the gradient start to decrease
            if d_gradient <= 0.0 {
                break;
            }
            d_max_theta = d_my_theta;
        }

        // Find out the biggest possible alpha for pivot
        pdqsort(&mut self.ph1_sorter_t[..]);
        let mut d_max_alpha: f64 = 0.0;
        let mut i_last: usize = self.ph1_sorter_t.len();
        for i in 0..self.ph1_sorter_t.len() {
            let d_my_theta: f64 = self.ph1_sorter_t[i].0;
            let index: HighsInt = self.ph1_sorter_t[i].1;
            let i_row: HighsInt = if index >= 0 {
                index
            } else {
                index + self.num_row
            };
            let d_abs_alpha: f64 = self.col_aq.array[i_row as usize].abs();
            // Stop when the theta is too large
            if d_my_theta > d_max_theta {
                i_last = i;
                break;
            }
            // Update the maximal possible alpha
            if d_max_alpha < d_abs_alpha {
                d_max_alpha = d_abs_alpha;
            }
        }

        // Finally choose a pivot with good enough alpha, working backwards
        self.row_out = kNoRowChosen;
        self.variable_out = -1;
        self.move_out = 0;
        let mut i = i_last;
        while i > 0 {
            let index: HighsInt = self.ph1_sorter_t[i - 1].1;
            let i_row: HighsInt = if index >= 0 {
                index
            } else {
                index + self.num_row
            };
            let d_abs_alpha: f64 = self.col_aq.array[i_row as usize].abs();
            if d_abs_alpha > d_max_alpha * 0.1 {
                self.row_out = i_row;
                self.move_out = if index >= 0 { 1 } else { -1 };
                break;
            }
            i -= 1;
        }
        // analysis->simplexTimerStop(Chuzr2Clock);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::chooseRow (lines 1268-1325, 58 lines)
    pub fn choose_row(&mut self, ekk: &mut HEkk) {
        let info = &ekk.info_;
        let base_lower = &info.base_lower_;
        let base_upper = &info.base_upper_;
        let base_value = &info.base_value_;
        // analysis->simplexTimerStart(Chuzr1Clock);
        // Initialize
        self.row_out = kNoRowChosen;

        // Choose row pass 1
        let alpha_tol: f64 = if info.update_count < 10 {
            1e-9
        } else if info.update_count < 20 {
            1e-8
        } else {
            1e-7
        };

        let mut relax_theta: f64 = 1e100;
        let mut relax_space: f64;
        for i in 0..self.col_aq.count {
            let i_row: HighsInt = self.col_aq.index[i as usize];
            let ru = i_row as usize;
            let alpha: f64 = self.col_aq.array[ru] * (self.move_in as f64);
            if alpha > alpha_tol {
                relax_space = base_value[ru] - base_lower[ru] + self.primal_feasibility_tolerance;
                if relax_space < relax_theta * alpha {
                    relax_theta = relax_space / alpha;
                }
            } else if alpha < -alpha_tol {
                relax_space = base_value[ru] - base_upper[ru] - self.primal_feasibility_tolerance;
                if relax_space > relax_theta * alpha {
                    relax_theta = relax_space / alpha;
                }
            }
        }
        // analysis->simplexTimerStop(Chuzr1Clock);

        // analysis->simplexTimerStart(Chuzr2Clock);
        let mut best_alpha: f64 = 0.0;
        for i in 0..self.col_aq.count {
            let i_row: HighsInt = self.col_aq.index[i as usize];
            let ru = i_row as usize;
            let alpha: f64 = self.col_aq.array[ru] * (self.move_in as f64);
            if alpha > alpha_tol {
                // Positive pivotal column entry
                let tight_space: f64 = base_value[ru] - base_lower[ru];
                if tight_space < relax_theta * alpha && best_alpha < alpha {
                    best_alpha = alpha;
                    self.row_out = i_row;
                }
            } else if alpha < -alpha_tol {
                // Negative pivotal column entry
                let tight_space: f64 = base_value[ru] - base_upper[ru];
                if tight_space > relax_theta * alpha && best_alpha < -alpha {
                    best_alpha = -alpha;
                    self.row_out = i_row;
                }
            }
        }
        // analysis->simplexTimerStop(Chuzr2Clock);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::considerBoundSwap (lines 1327-1394, 68 lines)
    pub fn consider_bound_swap(&mut self, ekk: &mut HEkk) {
        let info = &ekk.info_;
        let work_lower = &info.work_lower_;
        let work_upper = &info.work_upper_;
        let base_lower = &info.base_lower_;
        let base_upper = &info.base_upper_;
        let work_value = &info.work_value_;
        let base_value = &info.base_value_;

        // Compute the primal theta and see if we should have done a bound
        // flip instead
        if self.row_out == kNoRowChosen {
            debug_assert!(self.solve_phase == kSolvePhase2);
            // No binding ratio in CHUZR, so flip or unbounded
            self.theta_primal = (self.move_in as f64) * kHighsInf;
            self.move_out = 0;
        } else {
            debug_assert!(self.row_out >= 0);
            let ro = self.row_out as usize;
            // Determine the step to the leaving bound
            //
            self.alpha_col = self.col_aq.array[ro];
            // In Phase 1, move_out depends on whether the leaving variable is
            // becoming feasible - moves up to lower (down to upper) - or
            // remaining feasible - moves down to lower (up to upper) - so
            // can't be set so easily as in phase 2
            if self.solve_phase == kSolvePhase2 {
                self.move_out = if self.alpha_col * (self.move_in as f64) > 0.0 {
                    -1
                } else {
                    1
                };
            }
            self.theta_primal = 0.0;
            if self.move_out == 1 {
                self.theta_primal = (base_value[ro] - base_upper[ro]) / self.alpha_col;
            } else {
                self.theta_primal = (base_value[ro] - base_lower[ro]) / self.alpha_col;
            }
            debug_assert!(self.theta_primal > -kHighsInf && self.theta_primal < kHighsInf);
        }

        // Look to see if there is a bound flip
        let mut flipped = false;
        let lower_in: f64 = work_lower[self.variable_in as usize];
        let upper_in: f64 = work_upper[self.variable_in as usize];
        self.value_in = work_value[self.variable_in as usize] + self.theta_primal;
        if self.move_in > 0 {
            if self.value_in > upper_in + self.primal_feasibility_tolerance {
                flipped = true;
                self.row_out = kNoRowChosen;
                self.value_in = upper_in;
                self.theta_primal = upper_in - lower_in;
            }
        } else if self.value_in < lower_in - self.primal_feasibility_tolerance {
            flipped = true;
            self.row_out = kNoRowChosen;
            self.value_in = lower_in;
            self.theta_primal = lower_in - upper_in;
        }
        let pivot_or_flipped = self.row_out >= 0 || flipped;
        if self.solve_phase == kSolvePhase2 {
            // Check for possible unboundedness
            if !pivot_or_flipped {
                self.rebuild_reason = kRebuildReasonPossiblyPrimalUnbounded;
                return;
            }
        }
        // Check for possible error
        debug_assert!(pivot_or_flipped);
        debug_assert!(flipped == (self.row_out == kNoRowChosen));
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::assessPivot (lines 1396-1420, 25 lines)
    pub fn assess_pivot(&mut self, ekk: &mut HEkk) {
        debug_assert!(self.row_out >= 0);
        // Record the pivot entry
        self.alpha_col = self.col_aq.array[self.row_out as usize];
        self.variable_out = ekk.basis_.basic_index_[self.row_out as usize];

        // Compute the tableau row
        //
        // BTRAN
        //
        // Compute unit BTran for tableau row and FT update
        ekk.unit_btran(self.row_out, &mut self.row_ep);
        //
        // PRICE
        //
        let quad_precision = false;
        ekk.tableau_row_price(
            quad_precision,
            &self.row_ep,
            &mut self.row_ap,
            kDebugReportOff,
        );

        // Checks row-wise pivot against column-wise pivot for
        // numerical trouble
        //
        // rebuild_reason = kRebuildReasonPossiblySingularBasis is set if
        // numerical trouble is detected
        self.update_verify(ekk);
    }

    /// HiGHS: HEkkPrimal.cpp::HEkkPrimal::update (lines 1422-1550, 129 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion): accumulates `total_synthetic_tick_`.
    pub fn update(&mut self, ekk: &mut HEkk) {
        // Perform update operations that are independent of phase
        debug_assert!(self.rebuild_reason == 0);
        let flipped = self.row_out < 0;
        if flipped {
            self.variable_out = self.variable_in;
            self.alpha_col = 0.0;
            self.numerical_trouble = 0.0;
            ekk.info_.work_value_[self.variable_in as usize] = self.value_in;
            debug_assert!(
                ekk.basis_.nonbasic_move_[self.variable_in as usize] as HighsInt == self.move_in
            );
            ekk.basis_.nonbasic_move_[self.variable_in as usize] = (-self.move_in) as i8;
        } else {
            // Adjust perturbation if leaving equation
            self.adjust_perturbed_equation_out(ekk);
        }

        // Start hyper-sparse CHUZC, that takes place through phase1Update()
        self.hyper_choose_column_start(ekk);

        if self.solve_phase == kSolvePhase1 {
            // Update primal values
            self.phase1_update_primal(ekk);

            // Update the duals with respect to feasibility changes
            self.basic_feasibility_change_update_dual(ekk);

            // For hyper-sparse CHUZC, analyse the duals that have just changed
            self.hyper_choose_column_basic_feasibility_change(ekk);
        } else {
            // Update primal values, and identify any infeasibilities
            //
            // rebuild_reason =
            // kRebuildReasonPrimalInfeasibleInPrimalSimplex is set if a
            // primal infeasibility is found
            self.phase2_update_primal(ekk, false);
        }

        debug_assert!(
            self.rebuild_reason == kRebuildReasonNo
                || self.rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex
        );

        if flipped {
            ekk.info_.primal_bound_swap += 1;
            ekk.invalidate_dual_infeasibility_record();
            self.iteration_analysis(ekk);
            self.local_report_iter(ekk, false);
            self.num_flip_since_rebuild += 1;
            // Update the synthetic clock for UPDATE
            ekk.total_synthetic_tick_ += self.col_aq.synthetic_tick;
            // reporting-only work counter, never read by the search.
            crate::work::add_simplex_ticks(self.col_aq.synthetic_tick);
            return;
        }

        debug_assert!(self.row_out >= 0);

        // Now set the value of the entering variable
        ekk.info_.base_value_[self.row_out as usize] = self.value_in;
        // Consider whether the entering value is feasible and, if not, take
        // action
        //
        // rebuild_reason =
        // kRebuildReasonPrimalInfeasibleInPrimalSimplex is set in
        // phase 2 if a primal infeasibility is found
        self.consider_infeasible_value_in(ekk);

        // Update the dual values
        self.theta_dual = ekk.info_.work_dual_[self.variable_in as usize];
        self.update_dual(ekk);

        // Update any non-unit primal edge weights
        if self.edge_weight_mode == EdgeWeightMode::kDevex {
            self.update_devex(ekk);
        } else if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            self.debug_primal_steepest_edge_weights(ekk, "before update");
            self.update_primal_steepest_edge_weights(ekk);
        }

        // If entering column was nonbasic free, remove it from the set
        self.remove_nonbasic_free_column(ekk);

        // For hyper-sparse CHUZC, analyse the duals and weights that have
        // just changed
        self.hyper_choose_column_dual_change(ekk);

        if ekk.status_.has_dual_steepest_edge_weights {
            ekk.dev_debug_dual_steepest_edge_weights("before update");
            self.update_dual_steepest_edge_weights(ekk);
        }
        // Perform pivoting
        //
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
        ekk.update_pivots(self.variable_in, self.row_out, self.move_out);
        //
        // Update the invertible representation of the basis matrix
        ekk.update_factor(
            &mut self.col_aq,
            &mut self.row_ep,
            &mut self.row_out,
            &mut self.rebuild_reason,
        );

        if ekk.status_.has_dual_steepest_edge_weights {
            ekk.dev_debug_dual_steepest_edge_weights("after  update");
        }
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            self.debug_primal_steepest_edge_weights(ekk, "after update");
        }
        //
        // Update the row-wise representation of the nonbasic columns
        ekk.update_matrix(self.variable_in, self.variable_out);
        if ekk.info_.update_count >= ekk.info_.update_limit {
            self.rebuild_reason = kRebuildReasonUpdateLimitReached;
        }

        // Update the iteration count
        ekk.iteration_count_ += 1;

        // Reset the devex when there are too many errors
        if self.edge_weight_mode == EdgeWeightMode::kDevex
            && self.num_bad_devex_weight_ > kAllowedNumBadDevexWeight
        {
            self.initialise_devex_framework(ekk);
        }

        // Report on the iteration
        self.iteration_analysis(ekk);
        self.local_report_iter(ekk, false);

        // Update the synthetic clock for UPDATE
        ekk.total_synthetic_tick_ += self.col_aq.synthetic_tick;
        ekk.total_synthetic_tick_ += self.row_ep.synthetic_tick;
        // reporting-only work counter, never read by the search.
        crate::work::add_simplex_ticks(self.col_aq.synthetic_tick + self.row_ep.synthetic_tick);

        // Perform hyper-sparse CHUZC
        self.hyper_choose_column(ekk);
    }
}
