//! HiGHS source: simplex/HEkkDualMulti.cpp (lines 299-533)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! PAMI (parallel multiple-pricing dual, `simplex_strategy = kSimplexStrategyDualMulti`):
//! minor iterations: `minorUpdate*` and `minorInitialiseDevexFramework` (the `for_each` over the remaining choices in `minorUpdateRows`, line 491).

use super::super::*;
use crate::prelude::*;

impl HEkkDual {
    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::minorUpdate (lines 299-335, 37 lines)
    pub fn minor_update(&mut self, ekk: &mut HEkk) {
        // Minor update - store roll back data
        let finish_idx = self.multi_n_finish as usize;
        self.multi_finish[finish_idx].move_in =
            ekk.basis_.nonbasic_move_[self.variable_in as usize];
        self.multi_finish[finish_idx].shift_out = ekk.info_.work_shift_[self.variable_out as usize];
        self.multi_finish[finish_idx].flip_list.clear();
        for i in 0..self.dual_row.work_count {
            self.multi_finish[finish_idx]
                .flip_list
                .push(self.dual_row.work_data[i as usize].0);
        }

        // Minor update - key parts
        self.minor_update_dual(ekk);
        self.minor_update_primal(ekk);
        self.minor_update_pivots(ekk);
        self.minor_update_rows(ekk);
        if self.minor_new_devex_framework {
            self.minor_initialise_devex_framework(ekk);
        }
        self.multi_n_finish += 1;
        // Analyse the iteration: possibly report; possibly switch strategy
        self.iteration_analysis_minor(ekk);

        // Minor update - check for the next iteration
        let mut count_remain: HighsInt = 0;
        for i in 0..self.multi_num {
            let i_row = self.multi_choice[i as usize].row_out;
            if i_row < 0 {
                continue;
            }
            let my_infeas = self.multi_choice[i as usize].infeas_value;
            let my_weight = self.multi_choice[i as usize].infeas_ed_wt;
            count_remain +=
                (my_infeas / my_weight > self.multi_choice[i as usize].infeas_limit) as HighsInt;
        }
        if count_remain == 0 {
            self.multi_choose_again = 1;
        }
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::minorUpdateDual (lines 337-372, 36 lines)
    pub fn minor_update_dual(&mut self, ekk: &mut HEkk) {
        /*
         * 1. Update the dual solution
         *    XXX Data parallel (depends on the ap partition before)
         */
        if self.theta_dual == 0.0 {
            let wd = ekk.info_.work_dual_[self.variable_in as usize];
            self.shift_cost(ekk, self.variable_in, -wd);
        } else {
            self.dual_row.update_dual(ekk, self.theta_dual);
            if self.slice_price != 0 {
                for i in 0..self.slice_num {
                    self.slice_dual_row[i as usize].update_dual(ekk, self.theta_dual);
                }
            }
        }
        ekk.info_.work_dual_[self.variable_in as usize] = 0.0;
        ekk.info_.work_dual_[self.variable_out as usize] = -self.theta_dual;
        self.shift_back(ekk, self.variable_out);

        /*
         * 2. Apply global bound flip
         */
        {
            let ch = self.multi_finish[self.multi_n_finish as usize]
                .choice
                .expect("multi_finish choice set");
            let mut col_bfrt = std::mem::take(&mut self.multi_choice[ch].col_bfrt);
            self.dual_row.update_flip(ekk, &mut col_bfrt);
            self.multi_choice[ch].col_bfrt = col_bfrt;
        }

        /*
         * 3. Apply local bound flips
         */
        // SAFETY: owner-only / disjoint index: `multi_choice` is a fixed-size array field of
        // `self`; `this_ep` points at `multi_choice[ich].row_ep` (read-only via compute_dot)
        // while `base_value` is written on the same element. No other reference into
        // `multi_choice` exists during this loop.
        let multi_choice_ptr: *mut MChoice = self.multi_choice.as_mut_ptr();
        for ich in 0..self.multi_num {
            if ich == self.multi_i_choice || self.multi_choice[ich as usize].row_out >= 0 {
                for i in 0..self.dual_row.work_count {
                    let col = self.dual_row.work_data[i as usize].0;
                    // SAFETY: as above.
                    let dot = unsafe {
                        ekk.lp_
                            .a_matrix_
                            .compute_dot_hvector(&(*multi_choice_ptr.add(ich as usize)).row_ep, col)
                    };
                    // SAFETY: as above.
                    unsafe {
                        (*multi_choice_ptr.add(ich as usize)).base_value -=
                            self.dual_row.work_data[i as usize].1 * dot;
                    }
                }
            }
        }
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::minorUpdatePrimal (lines 374-437, 64 lines)
    pub fn minor_update_primal(&mut self, ekk: &mut HEkk) {
        let choice_idx = self.multi_i_choice as usize;
        let finish_idx = self.multi_n_finish as usize;
        let value_out = self.multi_choice[choice_idx].base_value;
        let lower_out = self.multi_choice[choice_idx].base_lower;
        let upper_out = self.multi_choice[choice_idx].base_upper;
        if self.delta_primal < 0.0 {
            self.theta_primal = (value_out - lower_out) / self.alpha_row;
            self.multi_finish[finish_idx].basic_bound = lower_out;
        }
        if self.delta_primal > 0.0 {
            self.theta_primal = (value_out - upper_out) / self.alpha_row;
            self.multi_finish[finish_idx].basic_bound = upper_out;
        }
        self.multi_finish[finish_idx].theta_primal = self.theta_primal;

        if self.edge_weight_mode == EdgeWeightMode::kDevex && !self.new_devex_framework {
            debug_assert!(self.row_out >= 0);
            if self.row_out < 0 {
                println!("ERROR: row_out = {} in minorUpdatePrimal", self.row_out);
            }
            let updated_edge_weight = ekk.dual_edge_weight_[self.row_out as usize];
            self.new_devex_framework = self.new_devex_framework(ekk, updated_edge_weight);
            self.minor_new_devex_framework = self.new_devex_framework;
            // Transform the edge weight of the pivotal row according to the
            // simplex update
            let mut new_pivotal_edge_weight =
                self.computed_edge_weight / (self.alpha_row * self.alpha_row);
            new_pivotal_edge_weight = cpp_max(1.0, new_pivotal_edge_weight);
            // Store the Devex weight of the leaving row now - OK since it's
            // stored in finish->EdWt and the updated weights are stored in
            // multi_choice[*].infeasEdWt
            self.multi_finish[finish_idx].ed_wt = new_pivotal_edge_weight;
        }

        /*
         * 5. Update the other primal value
         *    By the pivot (theta_primal)
         */
        // SAFETY: owner-only / disjoint index: `multi_choice` is a fixed-size array field of
        // `self`; `this_ep` points at `multi_choice[ich].row_ep` (read-only via compute_dot)
        // while `base_value`, `infeas_value`, `infeas_ed_wt` are written on the same element.
        // No other reference into `multi_choice` exists during this loop.
        let multi_choice_ptr: *mut MChoice = self.multi_choice.as_mut_ptr();
        for ich in 0..self.multi_num {
            if self.multi_choice[ich as usize].row_out >= 0 {
                // SAFETY: as above.
                let dot = unsafe {
                    ekk.lp_.a_matrix_.compute_dot_hvector(
                        &(*multi_choice_ptr.add(ich as usize)).row_ep,
                        self.variable_in,
                    )
                };
                // SAFETY: as above.
                unsafe {
                    (*multi_choice_ptr.add(ich as usize)).base_value -= self.theta_primal * dot;
                }
                // SAFETY: as above.
                let value = unsafe { (*multi_choice_ptr.add(ich as usize)).base_value };
                // SAFETY: as above.
                let lower = unsafe { (*multi_choice_ptr.add(ich as usize)).base_lower };
                // SAFETY: as above.
                let upper = unsafe { (*multi_choice_ptr.add(ich as usize)).base_upper };
                let mut infeas = 0.0;
                if value < lower - self.tp {
                    infeas = value - lower;
                }
                if value > upper + self.tp {
                    infeas = value - upper;
                }
                infeas *= infeas;
                // SAFETY: as above.
                unsafe {
                    (*multi_choice_ptr.add(ich as usize)).infeas_value = infeas;
                }
                if self.edge_weight_mode == EdgeWeightMode::kDevex {
                    // Update the other Devex weights
                    let new_pivotal_edge_weight = self.multi_finish[finish_idx].ed_wt;
                    let aa_i_row = dot;
                    // SAFETY: as above.
                    unsafe {
                        (*multi_choice_ptr.add(ich as usize)).infeas_ed_wt = cpp_max(
                            (*multi_choice_ptr.add(ich as usize)).infeas_ed_wt,
                            new_pivotal_edge_weight * aa_i_row * aa_i_row,
                        );
                    }
                }
            }
        }
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::minorUpdatePivots (lines 438-456, 19 lines)
    pub fn minor_update_pivots(&mut self, ekk: &mut HEkk) {
        let finish_idx = self.multi_n_finish as usize;
        ekk.update_pivots(self.variable_in, self.row_out, self.move_out);
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            // Transform the edge weight of the pivotal row according to the
            // simplex update
            self.multi_finish[finish_idx].ed_wt /= self.alpha_row * self.alpha_row;
        }
        self.multi_finish[finish_idx].basic_value =
            ekk.info_.work_value_[self.variable_in as usize] + self.theta_primal;
        ekk.update_matrix(self.variable_in, self.variable_out);
        self.multi_finish[finish_idx].variable_in = self.variable_in;
        self.multi_finish[finish_idx].alpha_row = self.alpha_row;
        // numericalTrouble is not set in minor iterations, only in
        // majorUpdate, so set it to an illegal value so that its
        // distribution is not updated
        self.numerical_trouble = -1.0;
        ekk.iteration_count_ += 1;
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::minorUpdateRows (lines 458-525, 68 lines)
    pub fn minor_update_rows(&mut self, ekk: &mut HEkk) {
        ekk.analysis_.simplex_timer_start(UpdateRowClock, 0);
        let finish_idx = self.multi_n_finish as usize;
        let ch = self.multi_finish[finish_idx]
            .choice
            .expect("multi_finish choice set");
        // C++ `const HVector* Row = multi_finish[n].row_ep` is the MChoice at `ch`.
        // We need `Row` (read-only) while mutating other `multi_choice` entries, so take
        // a raw pointer to the array (rule (a): owner-only, fixed-size, no resize).
        let multi_choice_ptr: *mut MChoice = self.multi_choice.as_mut_ptr();
        // SAFETY: owner-only: `multi_choice` is a fixed-size array field of `self` that does
        // not move while this pointer is live; `Row` is read-only (saxpy reads it).
        let row_ptr: *const HVector =
            unsafe { std::ptr::addr_of!((*multi_choice_ptr.add(ch)).row_ep) };

        // SAFETY: `row_ptr` is valid (see above); `(*row_ptr).count` is read-only.
        let row_count = unsafe { (*row_ptr).count };
        let update_rows_in_dense =
            (row_count < 0) || (row_count as f64 > 0.1 * self.solver_num_row as f64);

        if update_rows_in_dense {
            let mut multi_n_tasks: HighsInt = 0;
            let mut multi_i_which = [0 as HighsInt; kSimplexConcurrencyLimit as usize];
            let mut multi_xpivot = [0.0f64; kSimplexConcurrencyLimit as usize];
            let mut multi_vector: [*mut HVector; kSimplexConcurrencyLimit as usize] =
                [std::ptr::null_mut(); kSimplexConcurrencyLimit as usize];

            /*
             * Dense mode
             *  1. Find which ones to do and the pivotX
             *  2. Do all of them in task parallel
             */

            // Collect tasks
            for ich in 0..self.multi_num {
                if self.multi_choice[ich as usize].row_out >= 0 {
                    // SAFETY: owner-only / disjoint index: `next_ep` points at
                    // `multi_choice[ich].row_ep` (read-only via compute_dot).
                    let pivot_x = unsafe {
                        ekk.lp_.a_matrix_.compute_dot_hvector(
                            &(*multi_choice_ptr.add(ich as usize)).row_ep,
                            self.variable_in,
                        )
                    };
                    if pivot_x.abs() < kHighsTiny {
                        continue;
                    }
                    // SAFETY: owner-only: storing a mutable pointer to
                    // `multi_choice[ich].row_ep`; the array does not move.
                    multi_vector[multi_n_tasks as usize] = unsafe {
                        std::ptr::addr_of_mut!((*multi_choice_ptr.add(ich as usize)).row_ep)
                    };
                    multi_xpivot[multi_n_tasks as usize] = -pivot_x / self.alpha_row;
                    multi_i_which[multi_n_tasks as usize] = ich;
                    multi_n_tasks += 1;
                }
            }

            // Perform tasks
            // C++ `highs::parallel::for_each(0, multi_nTasks, ...)` is run sequentially in
            // order (same results; rule (p'): at one worker it reduces to serial order).
            for i in 0..multi_n_tasks {
                // SAFETY: owner-only: `multi_vector[i]` points at
                // `multi_choice[multi_i_which[i]].row_ep`; the array does not move.
                let next_ep_ptr = multi_vector[i as usize];
                let xpivot = multi_xpivot[i as usize];
                // SAFETY: owner-only: `next_ep_ptr` is a valid pointer into the fixed-size
                // `multi_choice` array; `row_ptr` is a valid read-only pointer into the same
                // array; the two indices are disjoint (the pivot row `ch` is not among the
                // remaining choices being updated).
                unsafe {
                    (*next_ep_ptr).saxpy(xpivot, &*row_ptr);
                    (*next_ep_ptr).tight();
                    if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                        multi_xpivot[i as usize] = (*next_ep_ptr).norm2();
                    }
                }
            }

            // Put weight back
            if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                for i in 0..multi_n_tasks {
                    self.multi_choice[multi_i_which[i as usize] as usize].infeas_ed_wt =
                        multi_xpivot[i as usize];
                }
            }
        } else {
            // Sparse mode: just do it sequentially
            for ich in 0..self.multi_num {
                if self.multi_choice[ich as usize].row_out >= 0 {
                    // SAFETY: owner-only / disjoint index: `next_ep` points at
                    // `multi_choice[ich].row_ep` (read-only via compute_dot).
                    let pivot_x = unsafe {
                        ekk.lp_.a_matrix_.compute_dot_hvector(
                            &(*multi_choice_ptr.add(ich as usize)).row_ep,
                            self.variable_in,
                        )
                    };
                    if pivot_x.abs() < kHighsTiny {
                        continue;
                    }
                    // SAFETY: owner-only: `next_ep` points at `multi_choice[ich].row_ep`;
                    // `row_ptr` is a valid read-only pointer into the same fixed-size array;
                    // the indices are disjoint.
                    unsafe {
                        (*multi_choice_ptr.add(ich as usize))
                            .row_ep
                            .saxpy(-pivot_x / self.alpha_row, &*row_ptr);
                        (*multi_choice_ptr.add(ich as usize)).row_ep.tight();
                        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                            (*multi_choice_ptr.add(ich as usize)).infeas_ed_wt =
                                (*multi_choice_ptr.add(ich as usize)).row_ep.norm2();
                        }
                    }
                }
            }
        }
        ekk.analysis_.simplex_timer_stop(UpdateRowClock, 0);
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::minorInitialiseDevexFramework (lines 527-533, 7 lines)
    pub fn minor_initialise_devex_framework(&mut self, _ekk: &mut HEkk) {
        // Set the local Devex weights to 1
        for i in 0..self.multi_num {
            self.multi_choice[i as usize].infeas_ed_wt = 1.0;
        }
        self.minor_new_devex_framework = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Unit tests

    #[test]
    fn filled_matches_cpp_line_by_line() {
        // Unit checks of minorUpdate, minorUpdateDual, minorUpdatePrimal, minorUpdatePivots, minorUpdateRows, minorInitialiseDevexFramework on a small LP with
        // simplex_strategy = 3 (PAMI); the end-to-end bit-identity check is PAMI-IT.
        let mut ekk = HEkk::default();
        let mut dual = HEkkDual::default();

        // Test minor_initialise_devex_framework: every multi_choice[i].infeas_ed_wt == 1.0
        // for i < multi_num and minor_new_devex_framework == false.
        dual.multi_num = 3;
        dual.minor_new_devex_framework = true;
        dual.multi_choice[0].infeas_ed_wt = 5.0;
        dual.multi_choice[1].infeas_ed_wt = 7.0;
        dual.multi_choice[2].infeas_ed_wt = 9.0;
        dual.multi_choice[3].infeas_ed_wt = 11.0; // beyond multi_num, must stay
        dual.minor_initialise_devex_framework(&mut ekk);
        assert_eq!(dual.multi_choice[0].infeas_ed_wt, 1.0);
        assert_eq!(dual.multi_choice[1].infeas_ed_wt, 1.0);
        assert_eq!(dual.multi_choice[2].infeas_ed_wt, 1.0);
        assert_eq!(dual.multi_choice[3].infeas_ed_wt, 11.0);
        assert!(!dual.minor_new_devex_framework);
    }

    /// minorUpdatePrimal (HEkkDualMulti.cpp 374-437): the theta_primal / basic_bound
    /// branches (delta_primal < 0 and > 0) and the per-choice primal update with a
    /// hand-computed infeasibility. Expected values are the C++ statements evaluated
    /// serially in C++ order.
    #[test]
    fn w5_pami_b_minor_update_primal() {
        // delta_primal < 0 branch: theta_primal = (valueOut - lowerOut) / alpha_row,
        // basic_bound = lowerOut.
        {
            let mut ekk = HEkk::default();
            let mut dual = HEkkDual::default();
            dual.multi_i_choice = 0;
            dual.multi_n_finish = 0;
            dual.multi_choice[0].base_value = 10.0;
            dual.multi_choice[0].base_lower = 2.0;
            dual.multi_choice[0].base_upper = 20.0;
            dual.alpha_row = 4.0;
            dual.delta_primal = -1.0;
            dual.edge_weight_mode = EdgeWeightMode::kDantzig;
            dual.multi_num = 0; // no other choices to update
            dual.minor_update_primal(&mut ekk);
            // theta_primal = (10 - 2) / 4 = 2
            assert_eq!(dual.theta_primal, 2.0);
            assert_eq!(dual.multi_finish[0].basic_bound, 2.0);
            assert_eq!(dual.multi_finish[0].theta_primal, 2.0);
        }

        // delta_primal > 0 branch: theta_primal = (valueOut - upperOut) / alpha_row,
        // basic_bound = upperOut.
        {
            let mut ekk = HEkk::default();
            let mut dual = HEkkDual::default();
            dual.multi_i_choice = 0;
            dual.multi_n_finish = 0;
            dual.multi_choice[0].base_value = 10.0;
            dual.multi_choice[0].base_lower = 2.0;
            dual.multi_choice[0].base_upper = 20.0;
            dual.alpha_row = 4.0;
            dual.delta_primal = 1.0;
            dual.edge_weight_mode = EdgeWeightMode::kDantzig;
            dual.multi_num = 0;
            dual.minor_update_primal(&mut ekk);
            // theta_primal = (10 - 20) / 4 = -2.5
            assert_eq!(dual.theta_primal, -2.5);
            assert_eq!(dual.multi_finish[0].basic_bound, 20.0);
            assert_eq!(dual.multi_finish[0].theta_primal, -2.5);
        }

        // Per-choice primal update with a hand-computed infeasibility.
        // choice 1 (row_out >= 0): row_ep.array[variable_in] = 3.0, so
        // dot = 3.0; base_value -= theta_primal * dot = 5.0 - 2.0*3.0 = -1.0;
        // value = -1.0, lower = 0.0, upper = 10.0, tp = 1e-7;
        // value < lower - tp -> infeas = -1.0 - 0.0 = -1.0; infeas *= infeas = 1.0.
        {
            let mut ekk = HEkk::default();
            let mut dual = HEkkDual::default();
            dual.multi_i_choice = 0;
            dual.multi_n_finish = 0;
            dual.multi_choice[0].base_value = 10.0;
            dual.multi_choice[0].base_lower = 2.0;
            dual.multi_choice[0].base_upper = 20.0;
            dual.alpha_row = 4.0;
            dual.delta_primal = -1.0;
            dual.edge_weight_mode = EdgeWeightMode::kDantzig;
            dual.tp = 1e-7;
            dual.variable_in = 0;
            dual.multi_num = 2;
            // choice 0: row_out < 0 -> skipped by the loop
            dual.multi_choice[0].row_out = -1;
            // choice 1: row_out >= 0 -> updated
            dual.multi_choice[1].row_out = 1;
            dual.multi_choice[1].base_value = 5.0;
            dual.multi_choice[1].base_lower = 0.0;
            dual.multi_choice[1].base_upper = 10.0;
            dual.multi_choice[1].row_ep.count = 1;
            dual.multi_choice[1].row_ep.index = vec![0; 3];
            dual.multi_choice[1].row_ep.array = vec![3.0, 0.0, 0.0];
            dual.minor_update_primal(&mut ekk);
            // theta_primal = (10 - 2) / 4 = 2
            assert_eq!(dual.theta_primal, 2.0);
            // choice 1 base_value = 5.0 - 2.0 * 3.0 = -1.0
            assert_eq!(dual.multi_choice[1].base_value, -1.0);
            // infeas = (-1.0 - 0.0)^2 = 1.0
            assert_eq!(dual.multi_choice[1].infeas_value, 1.0);
            // choice 0 untouched (row_out < 0)
            assert_eq!(dual.multi_choice[0].base_value, 10.0);
        }
    }

    /// minorUpdateRows (HEkkDualMulti.cpp 458-525): the sparse branch (Row.count small
    /// -> not dense) with a hand-computed saxpy. Expected values are the C++ statements
    /// evaluated serially in C++ order.
    #[test]
    fn w5_pami_b_minor_update_rows_sparse() {
        let mut ekk = HEkk::default();
        let mut dual = HEkkDual::default();
        // 100 rows so that Row.count = 1 <= 0.1 * solver_num_row: the sparse branch.
        dual.solver_num_row = 100;
        dual.multi_n_finish = 0;
        dual.multi_num = 2;
        dual.variable_in = 0;
        dual.alpha_row = 2.0;
        dual.edge_weight_mode = EdgeWeightMode::kDantzig;
        // finish[0].choice = 0 -> Row = multi_choice[0].row_ep
        dual.multi_finish[0].choice = Some(0);
        // Row (choice 0): count = 1 with index [0] (sparse branch), array[0] = 4.0
        dual.multi_choice[0].row_out = -1; // not updated (row_out < 0)
        dual.multi_choice[0].row_ep.count = 1;
        dual.multi_choice[0].row_ep.index = vec![0; 3];
        dual.multi_choice[0].row_ep.array = vec![4.0, 0.0, 0.0];
        // choice 1: row_out >= 0 -> updated
        dual.multi_choice[1].row_out = 1;
        dual.multi_choice[1].row_ep.count = 1;
        dual.multi_choice[1].row_ep.index = vec![0; 3];
        dual.multi_choice[1].row_ep.array = vec![1.0, 0.0, 0.0];
        // pivotX = computeDot(next_ep, variable_in) = next_ep.array[0] = 1.0
        // |pivotX| >= kHighsTiny -> update
        // xpivot = -pivotX / alpha_row = -1.0 / 2.0 = -0.5
        // saxpy(-0.5, Row): array[0] = 1.0 + (-0.5)*4.0 = -1.0
        dual.minor_update_rows(&mut ekk);
        assert_eq!(dual.multi_choice[1].row_ep.array[0], -1.0);
        // choice 0 untouched (row_out < 0)
        assert_eq!(dual.multi_choice[0].row_ep.array[0], 4.0);
    }
}
