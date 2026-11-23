//! HiGHS source: simplex/HEkkDualMulti.cpp (lines 26-297)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! PAMI (parallel multiple-pricing dual, `simplex_strategy = kSimplexStrategyDualMulti`):
//! CHUZR: `iterateMulti`, `majorChooseRow`, `majorChooseRowBtran` (the parallel BTRAN `for_each`, line 201) and `minorChooseRow`.

use super::super::*;
use crate::prelude::*;

// SAFETY: The raw pointers wrapped here are only dereferenced inside a
// `for_each` body where each task touches a disjoint index (rule a/s).
#[derive(Clone, Copy)]
struct SyncPtr<T>(*mut T);
unsafe impl<T> Sync for SyncPtr<T> {}

impl<T> SyncPtr<T> {
    /// Returns the wrapped raw pointer; a method call makes the closure capture the whole `SyncPtr`.
    #[inline]
    fn get(&self) -> *mut T {
        self.0
    }
}

impl HEkkDual {
    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::iterateMulti (lines 26-65, 40 lines)
    pub fn iterate_multi(&mut self, ekk: &mut HEkk) {
        self.slice_price = 1;

        // Report candidate
        self.major_choose_row(ekk);
        self.minor_choose_row(ekk);
        if self.row_out == kNoRowChosen {
            self.rebuild_reason = kRebuildReasonPossiblyOptimal;
            return;
        }

        // Assign the slice_row_ep, skip if possible
        let finish_choice = self.multi_finish[self.multi_n_finish as usize]
            .choice
            .unwrap();
        if (1.0 * self.multi_choice[finish_choice].row_ep.count as f64 / self.solver_num_row as f64)
            < 0.01
        {
            self.slice_price = 0;
        }

        if self.slice_price != 0 {
            // #pragma omp parallel
            // #pragma omp single
            let mut row_ep = std::mem::take(&mut self.multi_choice[finish_choice].row_ep);
            self.choose_column_slice(ekk, &mut row_ep);
            self.multi_choice[finish_choice].row_ep = row_ep;
        } else {
            let mut row_ep = std::mem::take(&mut self.multi_choice[finish_choice].row_ep);
            self.choose_column(ekk, &mut row_ep);
            self.multi_choice[finish_choice].row_ep = row_ep;
        }
        // If we failed.
        if self.rebuild_reason != 0 {
            if self.multi_n_finish != 0 {
                self.major_update(ekk);
            } else {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "PAMI skipping majorUpdate() due to multi_nFinish = {}; rebuild_reason = {}\n",
                    self.multi_n_finish,
                    self.rebuild_reason
                );
            }
            return;
        }

        self.minor_update(ekk);
        self.major_update(ekk);
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorChooseRow (lines 67-171, 105 lines)
    pub fn major_choose_row(&mut self, ekk: &mut HEkk) {
        // 0. Initial check to see if we need to do it again
        if ekk.info_.update_count == 0 {
            self.multi_choose_again = 1;
        }
        if self.multi_choose_again == 0 {
            return;
        }
        self.multi_choose_again = 0;
        self.multi_iteration += 1;
        // Major loop:
        //     repeat 1-5, until we found a good sets of choices
        let mut choice_index = vec![0 as HighsInt; self.multi_num as usize];
        loop {
            // 1. Multiple CHUZR
            let mut initial_count: HighsInt = 0;

            // Call the hyper-graph method, but partSwitch=0 so just uses
            // choose_multi_global
            self.dual_rhs.choose_multi_hyper_graph_auto(
                ekk,
                &mut choice_index,
                &mut initial_count,
                self.multi_num,
            );
            // dualRHS.chooseMultiGlobal(choiceIndex.data(), &initialCount, multi_num);
            if initial_count == 0 && self.dual_rhs.work_cutoff == 0.0 {
                // OPTIMAL
                return;
            }

            // 2. Shrink the size by cutoff
            let mut choice_count: HighsInt = 0;
            for i in 0..initial_count {
                let i_row = choice_index[i as usize];
                if self.dual_rhs.work_infeasibility[i_row as usize]
                    / ekk.dual_edge_weight_[i_row as usize]
                    >= self.dual_rhs.work_cutoff
                {
                    choice_index[choice_count as usize] = i_row;
                    choice_count += 1;
                }
            }

            if initial_count == 0 || choice_count <= initial_count / 3 {
                // Need to do the list again
                self.dual_rhs
                    .create_infeas_list(ekk, ekk.info_.col_aq_density);
                continue;
            }

            // 3. Store the choiceIndex to buffer
            for ich in 0..self.multi_num {
                self.multi_choice[ich as usize].row_out = kNoRowChosen;
            }
            for ich in 0..choice_count {
                self.multi_choice[ich as usize].row_out = choice_index[ich as usize];
            }

            // 4. Parallel BTRAN and compute weight
            self.major_choose_row_btran(ekk);

            // 5. Update row densities
            for ich in 0..self.multi_num {
                if self.multi_choice[ich as usize].row_out >= 0 {
                    let local_row_ep_density = self.multi_choice[ich as usize].row_ep.count as f64
                        / self.solver_num_row as f64;
                    let mut density = ekk.info_.row_ep_density;
                    ekk.update_operation_result_density(local_row_ep_density, &mut density);
                    ekk.info_.row_ep_density = density;
                }
            }

            if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                // 6. Check updated and computed weight - just for dual steepest edge
                let mut count_wrong_ed_wt: HighsInt = 0;
                for i in 0..self.multi_num {
                    let i_row = self.multi_choice[i as usize].row_out;
                    if i_row < 0 {
                        continue;
                    }
                    let updated_edge_weight = ekk.dual_edge_weight_[i_row as usize];
                    ekk.dual_edge_weight_[i_row as usize] =
                        self.multi_choice[i as usize].infeas_ed_wt;
                    self.computed_edge_weight = ekk.dual_edge_weight_[i_row as usize];
                    //      if (updated_edge_weight < 0.25 * computed_edge_weight) {
                    if !self.accept_dual_steepest_edge_weight(ekk, updated_edge_weight) {
                        self.multi_choice[i as usize].row_out = kNoRowChosen;
                        count_wrong_ed_wt += 1;
                    }
                }
                if count_wrong_ed_wt <= choice_count / 3 {
                    break;
                }
            } else {
                // No checking required if not using dual steepest edge so break
                break;
            }
        }

        // 6. Take other info associated with choices
        self.multi_chosen = 0;
        let k_pami_cutoff: f64 = 0.95;
        for i in 0..self.multi_num {
            let i_row = self.multi_choice[i as usize].row_out;
            if i_row < 0 {
                continue;
            }
            self.multi_chosen += 1;
            // Other info
            self.multi_choice[i as usize].base_value = ekk.info_.base_value_[i_row as usize];
            self.multi_choice[i as usize].base_lower = ekk.info_.base_lower_[i_row as usize];
            self.multi_choice[i as usize].base_upper = ekk.info_.base_upper_[i_row as usize];
            self.multi_choice[i as usize].infeas_value =
                self.dual_rhs.work_infeasibility[i_row as usize];
            self.multi_choice[i as usize].infeas_ed_wt = ekk.dual_edge_weight_[i_row as usize];
            self.multi_choice[i as usize].infeas_limit = self.dual_rhs.work_infeasibility
                [i_row as usize]
                / ekk.dual_edge_weight_[i_row as usize];
            self.multi_choice[i as usize].infeas_limit *= k_pami_cutoff;
        }

        // 6. Finish count
        self.multi_n_finish = 0;
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorChooseRowBtran (lines 173-246, 74 lines)
    pub fn major_choose_row_btran(&mut self, ekk: &mut HEkk) {
        ekk.analysis_.simplex_timer_start(BtranClock, 0);

        // 4.1. Prepare BTRAN buffer
        let mut multi_ntasks: HighsInt = 0;
        let mut multi_i_row = [0 as HighsInt; kSimplexConcurrencyLimit as usize];
        let mut multi_i_which = [0 as HighsInt; kSimplexConcurrencyLimit as usize];
        let mut multi_ed_wt = [0.0f64; kSimplexConcurrencyLimit as usize];
        let mut multi_vector: [*mut HVector; kSimplexConcurrencyLimit as usize] =
            [std::ptr::null_mut(); kSimplexConcurrencyLimit as usize];
        for ich in 0..self.multi_num {
            if self.multi_choice[ich as usize].row_out >= 0 {
                multi_i_row[multi_ntasks as usize] = self.multi_choice[ich as usize].row_out;
                multi_vector[multi_ntasks as usize] =
                    &mut self.multi_choice[ich as usize].row_ep as *mut HVector;
                multi_i_which[multi_ntasks as usize] = ich;
                multi_ntasks += 1;
            }
        }

        if ekk.analysis_.analyse_simplex_summary_data {
            for i in 0..multi_ntasks {
                ekk.analysis_.operation_record_before_count(
                    kSimplexNlaBtranEp,
                    1,
                    ekk.info_.row_ep_density,
                );
            }
        }
        // 4.2 Perform BTRAN
        // #pragma omp parallel for schedule(static, 1)
        // printf("start %d tasks for btran\n", multi_ntasks);
        // std::vector<HighsInt> tmp(multi_ntasks);

        // SAFETY: multi_vector[i] points to multi_choice[ich].row_ep for distinct ich
        // (disjoint index, rule a/s). multi_ed_wt[i] is written at disjoint index i.
        // The for_each body runs tasks on distinct indices, so no two tasks touch the
        // same multi_choice entry or the same multi_ed_wt slot.
        let multi_vector_ptr = SyncPtr(multi_vector.as_mut_ptr());
        let multi_ed_wt_ptr = SyncPtr(multi_ed_wt.as_mut_ptr());
        let multi_i_row_ptr = SyncPtr(multi_i_row.as_mut_ptr());
        let edge_weight_mode = self.edge_weight_mode;
        let row_ep_density = ekk.info_.row_ep_density;
        let dual_edge_weight_ptr = SyncPtr(ekk.dual_edge_weight_.as_mut_ptr());
        let simplex_nla_ptr = SyncPtr(&mut ekk.simplex_nla_ as *mut HSimplexNla);
        let lp_ptr = SyncPtr(&mut ekk.lp_ as *mut HighsLp);
        // `NlaRefs` borrows the whole `basis_.basic_index_`, as every other `btran` caller does.
        let basic_index_ptr = SyncPtr(ekk.basis_.basic_index_.as_mut_ptr());
        let basic_index_len = ekk.basis_.basic_index_.len();
        let options_ptr = SyncPtr(&mut ekk.options_ as *mut HighsOptions);
        let analysis_ptr = SyncPtr(&mut ekk.analysis_ as *mut HighsSimplexAnalysis);

        crate::parallel::highs_parallel::for_each(
            0,
            multi_ntasks,
            &|start: HighsInt, end: HighsInt| {
                for i in start..end {
                    let i_row = unsafe { *multi_i_row_ptr.get().offset(i as isize) };
                    let work_ep = unsafe { &mut **multi_vector_ptr.get().offset(i as isize) };
                    work_ep.clear();
                    work_ep.count = 1;
                    work_ep.index[0] = i_row;
                    work_ep.array[i_row as usize] = 1.0;
                    work_ep.pack_flag = true;
                    let factor_timer_clock_pointer = unsafe {
                        (*analysis_ptr.get()).get_thread_factor_timer_clock_ptr(
                            crate::parallel::highs_parallel::thread_num(),
                        )
                    };
                    // SAFETY: refs borrows lp, basic_index, options which are not
                    // mutated during the for_each (disjoint from multi_choice writes).
                    let refs = unsafe {
                        NlaRefs::new(
                            &*lp_ptr.get(),
                            std::slice::from_raw_parts(basic_index_ptr.get(), basic_index_len),
                            &*options_ptr.get(),
                        )
                    };
                    unsafe {
                        (*simplex_nla_ptr.get()).btran(
                            &refs,
                            work_ep,
                            row_ep_density,
                            factor_timer_clock_pointer,
                        )
                    };
                    if edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                        // For Dual steepest edge we know the exact weight as the 2-norm of
                        // work_ep
                        unsafe {
                            *multi_ed_wt_ptr.get().offset(i as isize) = work_ep.norm2();
                        }
                    } else {
                        // For Devex (and Dantzig) we take the updated edge weight
                        unsafe {
                            *multi_ed_wt_ptr.get().offset(i as isize) =
                                *dual_edge_weight_ptr.get().offset(i_row as isize);
                        }
                    }
                }
            },
            1,
        );

        // printf("multi_ntask task schedule:");
        // for (HighsInt i = 0; i < multi_ntasks; ++i) {
        //  printf(" %d->w%d", i, tmp[i]);
        //}
        // printf("\n");

        if ekk.analysis_.analyse_simplex_summary_data {
            for i in 0..multi_ntasks {
                // SAFETY: multi_vector[i] is a valid pointer into multi_choice[].row_ep
                // (disjoint index, set up before the for_each).
                let count = unsafe { (*multi_vector[i as usize]).count };
                ekk.analysis_
                    .operation_record_after_count(kSimplexNlaBtranEp, count);
            }
        }
        // 4.3 Put back edge weights: the edge weights for the chosen rows
        // are stored in multi_choice[*].infeasEdWt
        for i in 0..multi_ntasks {
            self.multi_choice[multi_i_which[i as usize] as usize].infeas_ed_wt =
                multi_ed_wt[i as usize];
        }

        ekk.analysis_.simplex_timer_stop(BtranClock, 0);
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::minorChooseRow (lines 248-297, 50 lines)
    pub fn minor_choose_row(&mut self, ekk: &mut HEkk) {
        // 1. Find which to go out
        //        Because we had other checking code
        //        We know current best is OK to be used
        self.multi_i_choice = -1;
        let mut best_merit: f64 = 0.0;
        for ich in 0..self.multi_num {
            let i_row = self.multi_choice[ich as usize].row_out;
            if i_row < 0 {
                continue;
            }
            let infeas_value = self.multi_choice[ich as usize].infeas_value;
            let infeas_ed_wt = self.multi_choice[ich as usize].infeas_ed_wt;
            let infeas_merit = infeas_value / infeas_ed_wt;
            if best_merit < infeas_merit {
                best_merit = infeas_merit;
                self.multi_i_choice = ich;
            }
        }

        // 2. Obtain other info for current sub-optimization choice
        self.row_out = kNoRowChosen;
        if self.multi_i_choice != -1 {
            let work_choice_idx = self.multi_i_choice as usize;

            // Assign useful variables
            self.row_out = self.multi_choice[work_choice_idx].row_out;
            self.variable_out = ekk.basis_.basic_index_[self.row_out as usize];
            let value_out = self.multi_choice[work_choice_idx].base_value;
            let lower_out = self.multi_choice[work_choice_idx].base_lower;
            let upper_out = self.multi_choice[work_choice_idx].base_upper;
            self.delta_primal = value_out
                - (if value_out < lower_out {
                    lower_out
                } else {
                    upper_out
                });
            self.move_out = if self.delta_primal < 0.0 { -1 } else { 1 };

            // Assign buffers
            let finish_idx = self.multi_n_finish as usize;
            self.multi_finish[finish_idx].row_out = self.row_out;
            self.multi_finish[finish_idx].variable_out = self.variable_out;
            self.multi_finish[finish_idx].choice = Some(work_choice_idx);
            // Save the edge weight - over-written later when using Devex
            self.multi_finish[finish_idx].ed_wt = self.multi_choice[work_choice_idx].infeas_ed_wt;

            // Disable current row
            self.multi_choice[work_choice_idx].row_out = kNoRowChosen;
        }
    }
}

#[cfg(test)]
mod tests {
    // Unit tests

    use crate::prelude::*;
    use crate::simplex::dual::HEkkDual;
    use crate::simplex::simplex_const::kNoRowChosen;

    /// Hand-derived from HEkkDualMulti.cpp minorChooseRow (lines 248-297): the merit
    /// `infeasValue / infeasEdWt` picks the best candidate (strict `<`, first wins ties), the
    /// chosen slot is disabled (`row_out = kNoRowChosen`) and the MFinish record is saved.
    #[test]
    fn filled_matches_cpp_line_by_line() {
        let mut ekk = HEkk::default();
        ekk.basis_.basic_index_ = vec![0; 8];
        ekk.basis_.basic_index_[3] = 13;
        ekk.basis_.basic_index_[5] = 42;
        let mut dual = HEkkDual::default();
        dual.multi_num = 3;
        dual.multi_n_finish = 0;
        // merits: slot 0 = 2.0 / 1.0, slot 1 = 6.0 / 2.0 = 3.0 (best), slot 2 = -1 (disabled)
        dual.multi_choice[0].row_out = 3;
        dual.multi_choice[0].infeas_value = 2.0;
        dual.multi_choice[0].infeas_ed_wt = 1.0;
        dual.multi_choice[0].base_value = 0.5;
        dual.multi_choice[0].base_lower = 1.0;
        dual.multi_choice[0].base_upper = 4.0;
        dual.multi_choice[1].row_out = 5;
        dual.multi_choice[1].infeas_value = 6.0;
        dual.multi_choice[1].infeas_ed_wt = 2.0;
        dual.multi_choice[1].base_value = 7.0;
        dual.multi_choice[1].base_lower = 1.0;
        dual.multi_choice[1].base_upper = 4.0;
        dual.multi_choice[2].row_out = kNoRowChosen;
        dual.multi_choice[2].infeas_value = 100.0;
        dual.multi_choice[2].infeas_ed_wt = 1.0;

        dual.minor_choose_row(&mut ekk);

        assert_eq!(dual.multi_i_choice, 1);
        assert_eq!(dual.row_out, 5);
        assert_eq!(dual.variable_out, 42);
        // valueOut (7.0) is not below lowerOut, so deltaPrimal = 7.0 - upperOut (4.0) = 3.0
        assert_eq!(dual.delta_primal, 3.0);
        assert_eq!(dual.move_out, 1);
        assert_eq!(dual.multi_finish[0].row_out, 5);
        assert_eq!(dual.multi_finish[0].variable_out, 42);
        assert_eq!(dual.multi_finish[0].choice, Some(1));
        assert_eq!(dual.multi_finish[0].ed_wt, 2.0);
        assert_eq!(dual.multi_choice[1].row_out, kNoRowChosen);
        assert_eq!(dual.multi_choice[0].row_out, 3);
    }

    /// C++ minorChooseRow lines 256-266 and 281-282: the strict `bestMerit < infeasMerit`
    /// keeps the first of two equal merits, and `valueOut < lowerOut` gives
    /// `deltaPrimal = valueOut - lowerOut < 0`, so `move_out = -1`. Writes go to
    /// `multi_finish[multi_nFinish]`.
    #[test]
    fn w5_pami_a_minor_choose_row_tie_and_below_lower() {
        let mut ekk = HEkk::default();
        ekk.basis_.basic_index_ = vec![0; 8];
        ekk.basis_.basic_index_[2] = 7;
        ekk.basis_.basic_index_[6] = 9;
        let mut dual = HEkkDual::default();
        dual.multi_num = 3;
        dual.multi_n_finish = 1;
        dual.multi_choice[0].row_out = kNoRowChosen;
        dual.multi_choice[1].row_out = 2;
        dual.multi_choice[1].infeas_value = 4.0;
        dual.multi_choice[1].infeas_ed_wt = 2.0;
        dual.multi_choice[1].base_value = -3.0;
        dual.multi_choice[1].base_lower = -1.0;
        dual.multi_choice[1].base_upper = 5.0;
        dual.multi_choice[2].row_out = 6;
        dual.multi_choice[2].infeas_value = 2.0;
        dual.multi_choice[2].infeas_ed_wt = 1.0;

        dual.minor_choose_row(&mut ekk);

        assert_eq!(dual.multi_i_choice, 1);
        assert_eq!(dual.row_out, 2);
        assert_eq!(dual.variable_out, 7);
        assert_eq!(dual.delta_primal, -2.0);
        assert_eq!(dual.move_out, -1);
        assert_eq!(dual.multi_finish[1].row_out, 2);
        assert_eq!(dual.multi_finish[1].variable_out, 7);
        assert_eq!(dual.multi_finish[1].choice, Some(1));
        assert_eq!(dual.multi_finish[1].ed_wt, 2.0);
        assert_eq!(dual.multi_finish[0].choice, None);
        assert_eq!(dual.multi_choice[1].row_out, kNoRowChosen);
        assert_eq!(dual.multi_choice[2].row_out, 6);
    }

    /// C++ majorChooseRow lines 71-72: with `update_count != 0` and `multi_chooseAgain == 0`
    /// it returns before touching `multi_iteration`, the choices or `multi_nFinish`.
    #[test]
    fn w5_pami_a_major_choose_row_skips_when_not_choosing_again() {
        let mut ekk = HEkk::default();
        ekk.info_.update_count = 3;
        let mut dual = HEkkDual::default();
        dual.multi_num = 2;
        dual.multi_choose_again = 0;
        dual.multi_iteration = 5;
        dual.multi_n_finish = 1;
        dual.multi_choice[0].row_out = 4;
        dual.major_choose_row(&mut ekk);
        assert_eq!(dual.multi_choose_again, 0);
        assert_eq!(dual.multi_iteration, 5);
        assert_eq!(dual.multi_n_finish, 1);
        assert_eq!(dual.multi_choice[0].row_out, 4);
    }

    /// C++ minorChooseRow with no candidate: multiICHOICE = -1 and rowOut = kNoRowChosen.
    #[test]
    fn w5_pami_a_minor_choose_row_no_candidate() {
        let mut ekk = HEkk::default();
        let mut dual = HEkkDual::default();
        dual.multi_num = 2;
        dual.multi_choice[0].row_out = kNoRowChosen;
        dual.multi_choice[1].row_out = kNoRowChosen;
        dual.minor_choose_row(&mut ekk);
        assert_eq!(dual.multi_i_choice, -1);
        assert_eq!(dual.row_out, kNoRowChosen);
        assert_eq!(dual.multi_n_finish, 0);
    }
}
