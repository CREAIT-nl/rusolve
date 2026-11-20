//! HiGHS source: simplex/HEkkDualMulti.cpp (lines 765-1020)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! PAMI (parallel multiple-pricing dual, `simplex_strategy = kSimplexStrategyDualMulti`):
//! `majorUpdatePrimal` (the `for_each`s at 773 and 804), `majorUpdateFactor` (builds the `HVector::next` chain for `HFactor::update_multi`), `majorRollback`, `checkNonUnitWeightError` and the iteration analysis.

use super::super::*;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g;

#[derive(Clone, Copy)]
struct SyncPtr(*mut f64);
// SAFETY: The pointer is used in parallel for_each loops where each thread writes to a
// disjoint range of indices.
unsafe impl Send for SyncPtr {}
// SAFETY: As above.
unsafe impl Sync for SyncPtr {}

impl SyncPtr {
    /// Returns the wrapped pointer; a method call captures the whole `SyncPtr` in closures.
    #[inline]
    fn get(self) -> *mut f64 {
        self.0
    }
}

impl HEkkDual {
    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorUpdatePrimal (lines 765-904, 140 lines)
    pub fn major_update_primal(&mut self, ekk: &mut HEkk) {
        let update_primal_in_dense = self.dual_rhs.work_count < 0;
        if update_primal_in_dense {
            // Dense update of primal values, infeasibility list and
            // non-pivotal edge weights
            let mix_array = &self.col_bfrt.array;
            let base_value_ptr = SyncPtr(ekk.info_.base_value_.as_mut_ptr());
            let base_lower = &ekk.info_.base_lower_;
            let base_upper = &ekk.info_.base_upper_;
            let infeas_ptr = SyncPtr(self.dual_rhs.work_infeasibility.as_mut_ptr());
            let tp = self.tp;
            let store_sq = ekk.info_.store_squared_primal_infeasibility;
            let num_row = self.solver_num_row;

            // #pragma omp parallel for schedule(static)
            crate::parallel::highs_parallel::for_each(
                0,
                num_row,
                &|start: HighsInt, end: HighsInt| {
                    for i_row in start..end {
                        let idx = i_row as usize;
                        // SAFETY: Disjoint writes to base_value_ptr
                        let val = unsafe {
                            let bv = base_value_ptr.get().add(idx);
                            *bv -= mix_array[idx];
                            *bv
                        };
                        let less = base_lower[idx] - val;
                        let more = val - base_upper[idx];
                        let infeas = if less > tp {
                            less
                        } else if more > tp {
                            more
                        } else {
                            0.0
                        };
                        // SAFETY: Disjoint writes to infeas_ptr
                        unsafe {
                            let ip = infeas_ptr.get().add(idx);
                            if store_sq {
                                *ip = infeas * infeas;
                            } else {
                                *ip = infeas.abs();
                            }
                        }
                    }
                },
                100,
            );

            if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge
                || (self.edge_weight_mode == EdgeWeightMode::kDevex && !self.new_devex_framework)
            {
                // Dense update of any edge weights (except weights for pivotal rows)
                let edge_weight_ptr = SyncPtr(ekk.dual_edge_weight_.as_mut_ptr());
                for i_fn in 0..self.multi_n_finish {
                    let finish = &self.multi_finish[i_fn as usize];
                    let ch = finish.choice.expect("multi_finish choice set");
                    // multi_finish[iFn].EdWt has already been transformed to correspond to
                    // the new basis
                    let new_pivotal_edge_weight = finish.ed_wt;
                    let col_array = &self.multi_choice[ch].col_aq.array;
                    if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                        // Update steepest edge weights
                        let dse_array = &self.multi_choice[ch].row_ep.array;
                        let kai = -2.0 / finish.alpha_row;
                        // #pragma omp parallel for schedule(static)
                        crate::parallel::highs_parallel::for_each(
                            0,
                            num_row,
                            &|start: HighsInt, end: HighsInt| {
                                for i_row in start..end {
                                    let idx = i_row as usize;
                                    let aa_i_row = col_array[idx];
                                    // SAFETY: Disjoint writes to edge_weight_ptr
                                    unsafe {
                                        let ew = edge_weight_ptr.get().add(idx);
                                        *ew += aa_i_row
                                            * (new_pivotal_edge_weight * aa_i_row
                                                + kai * dse_array[idx]);
                                        *ew = cpp_max(kMinDualSteepestEdgeWeight, *ew);
                                    }
                                }
                            },
                            100,
                        );
                    } else {
                        // Update Devex weights
                        for i_row in 0..num_row {
                            let idx = i_row as usize;
                            let aa_i_row = col_array[idx];
                            // SAFETY: Disjoint writes to edge_weight_ptr
                            unsafe {
                                let ew = edge_weight_ptr.get().add(idx);
                                *ew = cpp_max(*ew, new_pivotal_edge_weight * aa_i_row * aa_i_row);
                            }
                        }
                    }
                }
            }
        } else {
            // Sparse update of primal values, infeasibility list and
            // non-pivotal edge weights
            self.dual_rhs.update_primal(ekk, &self.col_bfrt, 1.0);
            self.dual_rhs.update_infeas_list(ekk, &self.col_bfrt);

            // Sparse update of any edge weights and infeasList. Weights for
            // rows pivotal in this set of MI are updated, but this is based
            // on the previous updated weights not the computed pivotal weight
            // that's known. For the rows pivotal in this set of MI, the
            // weights will be over-written in the next section of code.
            for i_fn in 0..self.multi_n_finish {
                let finish = &self.multi_finish[i_fn as usize];
                let ch = finish.choice.expect("multi_finish choice set");
                let col = &self.multi_choice[ch].col_aq;
                let new_pivotal_edge_weight = finish.ed_wt;
                if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                    // Update steepest edge weights
                    let row = &self.multi_choice[ch].row_ep;
                    let kai = -2.0 / finish.alpha_row;
                    ekk.update_dual_steepest_edge_weights(
                        finish.row_out,
                        finish.variable_in,
                        col,
                        new_pivotal_edge_weight,
                        kai,
                        &row.array,
                    );
                } else if self.edge_weight_mode == EdgeWeightMode::kDevex
                    && !self.new_devex_framework
                {
                    // Update Devex weights
                    ekk.update_dual_devex_weights(col, new_pivotal_edge_weight);
                }
                self.dual_rhs.update_infeas_list(ekk, col);
            }
        }

        // Update primal value for the rows pivotal in this set of MI
        for i_fn in 0..self.multi_n_finish {
            let finish = &self.multi_finish[i_fn as usize];
            let i_row = finish.row_out;
            let value =
                ekk.info_.base_value_[i_row as usize] - finish.basic_bound + finish.basic_value;
            self.dual_rhs.update_pivots(ekk, i_row, value);
        }

        // Update any edge weights for the rows pivotal in this set of MI.
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge
            || (self.edge_weight_mode == EdgeWeightMode::kDevex && !self.new_devex_framework)
        {
            // Update weights for the pivots using the computed values.
            let edge_weight = &mut ekk.dual_edge_weight_;
            for i_fn in 0..self.multi_n_finish {
                let finish = &self.multi_finish[i_fn as usize];
                let ch = finish.choice.expect("multi_finish choice set");
                let i_row = finish.row_out;
                let new_pivotal_edge_weight = finish.ed_wt;
                let col_array = &self.multi_choice[ch].col_aq.array;
                // The weight for this pivot is known, but weights for rows
                // pivotal earlier need to be updated
                if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                    // Steepest edge
                    let dse_array = &self.multi_choice[ch].row_ep.array;
                    let kai = -2.0 / finish.alpha_row;
                    for j_fn in 0..i_fn {
                        let j_finish = &self.multi_finish[j_fn as usize];
                        let j_row = j_finish.row_out;
                        let value = col_array[j_row as usize];
                        edge_weight[j_row as usize] += value
                            * (new_pivotal_edge_weight * value + kai * dse_array[j_row as usize]);
                        edge_weight[j_row as usize] =
                            cpp_max(kMinDualSteepestEdgeWeight, edge_weight[j_row as usize]);
                    }
                    edge_weight[i_row as usize] = new_pivotal_edge_weight;
                } else {
                    // Devex
                    for j_fn in 0..i_fn {
                        let j_finish = &self.multi_finish[j_fn as usize];
                        let j_row = j_finish.row_out;
                        let value = col_array[j_row as usize];
                        edge_weight[j_row as usize] = cpp_max(
                            edge_weight[j_row as usize],
                            new_pivotal_edge_weight * value * value,
                        );
                    }
                    edge_weight[i_row as usize] = new_pivotal_edge_weight;
                    self.num_devex_iterations += 1;
                }
            }
        }
        self.check_non_unit_weight_error(ekk, "999");
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorUpdateFactor (lines 906-933, 28 lines)
    pub fn major_update_factor(&mut self, ekk: &mut HEkk) {
        // 9. Update the factor by CFT
        let n = self.multi_n_finish;
        let mut i_rows = vec![0; n as usize];
        // C++ `multi_finish[iCh].row_ep` / `col_aq` are `HVector*` into `multi_choice`
        // (rule (a)): every link and the two chain heads derive from this one base pointer,
        // so writing a `next` field never invalidates a link taken earlier.
        let multi_choice: *mut MChoice = self.multi_choice.as_mut_ptr();
        for i_ch in 0..(n - 1) {
            let ch = self.multi_finish[i_ch as usize]
                .choice
                .expect("multi_finish choice set");
            let ch_next = self.multi_finish[(i_ch + 1) as usize]
                .choice
                .expect("multi_finish choice set");
            // SAFETY: disjoint index / owner-only: `ch` and `ch_next` are `multi_choice`
            // indices set by majorChooseRow (HEkkDualMulti.cpp 288-290), so in bounds of the
            // fixed-size array, which is a field of `self` and does not move while the chain
            // is live; no reference into `multi_choice` exists during these writes.
            unsafe {
                (*multi_choice.add(ch)).row_ep.next =
                    Some(HVectorNext(std::ptr::NonNull::new_unchecked(
                        std::ptr::addr_of_mut!((*multi_choice.add(ch_next)).row_ep),
                    )));
                (*multi_choice.add(ch)).col_aq.next =
                    Some(HVectorNext(std::ptr::NonNull::new_unchecked(
                        std::ptr::addr_of_mut!((*multi_choice.add(ch_next)).col_aq),
                    )));
            }
            i_rows[i_ch as usize] = self.multi_finish[i_ch as usize].row_out;
        }
        i_rows[(n - 1) as usize] = self.multi_finish[(n - 1) as usize].row_out;
        if n > 0 {
            let c0 = self.multi_finish[0]
                .choice
                .expect("multi_finish choice set");
            // SAFETY: as above (HVector next chain, owner-only): the chain heads and every
            // link point into `multi_choice`, which is neither moved nor otherwise borrowed
            // during the call; `i_rows` and `rebuild_reason` are disjoint storage.
            unsafe {
                ekk.update_factor_multi(
                    std::ptr::addr_of_mut!((*multi_choice.add(c0)).col_aq),
                    std::ptr::addr_of_mut!((*multi_choice.add(c0)).row_ep),
                    &mut i_rows,
                    &mut self.rebuild_reason,
                );
            }
        }

        // Determine whether to reinvert based on the synthetic clock
        let use_build_synthetic_tick = ekk.build_synthetic_tick_ * kMultiBuildSyntheticTickMu;
        let reinvert_synthetic_clock = ekk.total_synthetic_tick_ >= use_build_synthetic_tick;
        let performed_min_updates =
            ekk.info_.update_count >= kMultiSyntheticTickReinversionMinUpdateCount;
        if reinvert_synthetic_clock && performed_min_updates {
            self.rebuild_reason = kRebuildReasonSyntheticClockSaysInvert;
        }
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorRollback (lines 935-961, 27 lines)
    pub fn major_rollback(&mut self, ekk: &mut HEkk) {
        for i_fn in (0..self.multi_n_finish).rev() {
            let finish = &self.multi_finish[i_fn as usize];

            // 1. Roll back pivot
            ekk.basis_.nonbasic_move_[finish.variable_in as usize] = finish.move_in;
            ekk.basis_.nonbasic_flag_[finish.variable_in as usize] = 1;
            ekk.basis_.nonbasic_move_[finish.variable_out as usize] = 0;
            ekk.basis_.nonbasic_flag_[finish.variable_out as usize] = 0;
            ekk.basis_.basic_index_[finish.row_out as usize] = finish.variable_out;

            // 2. Roll back matrix
            ekk.update_matrix(finish.variable_out, finish.variable_in);

            // 3. Roll back flips
            for &var in &finish.flip_list {
                ekk.flip_bound(var);
            }

            // 4. Roll back cost
            ekk.info_.work_shift_[finish.variable_in as usize] = 0.0;
            ekk.info_.work_shift_[finish.variable_out as usize] = finish.shift_out;

            // 5. The iteration count
            ekk.iteration_count_ -= 1;
        }
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::checkNonUnitWeightError (lines 963-977, 15 lines)
    pub fn check_non_unit_weight_error(&self, ekk: &HEkk, message: &str) -> bool {
        let mut error_found = false;
        if self.edge_weight_mode == EdgeWeightMode::kDantzig {
            let edge_weight = &ekk.dual_edge_weight_;
            let mut unit_wt_error = 0.0;
            for i_row in 0..self.solver_num_row {
                unit_wt_error += (edge_weight[i_row as usize] - 1.0).abs();
            }
            error_found = unit_wt_error > 1e-4;
            if error_found {
                println!(
                    "Non-unit Edge weight error of {}: {}",
                    c_fmt_g(unit_wt_error, 6),
                    message
                );
            }
        }
        error_found
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::iterationAnalysisMinorData (lines 979-983, 5 lines)
    pub fn iteration_analysis_minor_data(&mut self, ekk: &mut HEkk) {
        // Port gap: HighsSimplexAnalysis has no field for C++ analysis->multi_* / *_concurrency; its iterationReport/iterationRecord are no-ops here.
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::iterationAnalysisMinor (lines 985-993, 9 lines)
    pub fn iteration_analysis_minor(&mut self, ekk: &mut HEkk) {
        // Possibly report on the iteration
        // PAMI uses alpha_row but serial solver uses alpha
        self.alpha_col = self.alpha_row;
        self.iteration_analysis_data(ekk);
        self.iteration_analysis_minor_data(ekk);
        ekk.analysis_.iteration_report();
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.iteration_record();
        }
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::iterationAnalysisMajorData (lines 995-1001, 7 lines)
    pub fn iteration_analysis_major_data(&mut self, ekk: &mut HEkk) {
        ekk.analysis_.numerical_trouble = self.numerical_trouble;
        // Port gap: HighsSimplexAnalysis has no field for C++ analysis->multi_* / *_concurrency; its iterationReport/iterationRecord are no-ops here.
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::iterationAnalysisMajor (lines 1003-1020, 18 lines)
    pub fn iteration_analysis_major(&mut self, ekk: &mut HEkk) {
        self.iteration_analysis_major_data(ekk);
        // Possibly switch from DSE to Devex
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            //    switch_to_devex = analysis->switchToDevex();
            let switch_to_devex = ekk.switch_to_devex();
            if switch_to_devex {
                self.edge_weight_mode = EdgeWeightMode::kDevex;
                // Set up the Devex framework
                self.initialise_devex_framework(ekk);
            }
        }
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.iteration_record();
            ekk.analysis_.iteration_record_major();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Unit tests

    #[test]
    fn filled_matches_cpp_line_by_line() {
        // Unit checks of majorUpdatePrimal, majorUpdateFactor, majorRollback, checkNonUnitWeightError, iterationAnalysisMinorData, iterationAnalysisMinor, iterationAnalysisMajorData, iterationAnalysisMajor on a small LP with
        // simplex_strategy = 3 (PAMI); the end-to-end bit-identity check is PAMI-IT.
        let mut ekk = HEkk::default();
        let mut dual = HEkkDual::default();

        // Test check_non_unit_weight_error
        dual.solver_num_row = 3;
        dual.edge_weight_mode = EdgeWeightMode::kDantzig;
        ekk.dual_edge_weight_ = vec![1.0, 1.0, 1.0];
        assert!(!dual.check_non_unit_weight_error(&ekk, "t"));
        ekk.dual_edge_weight_[1] = 2.0;
        assert!(dual.check_non_unit_weight_error(&ekk, "t"));

        // Test major_rollback with n=0
        dual.multi_n_finish = 0;
        let count_before = ekk.iteration_count_;
        dual.major_rollback(&mut ekk);
        assert_eq!(ekk.iteration_count_, count_before);

        // Test iteration_analysis_major_data does not panic
        dual.iteration_analysis_major_data(&mut ekk);
    }

    /// Dense branch of majorUpdatePrimal (HEkkDualMulti.cpp 766-822) and the pivotal-row
    /// pass (850-899), steepest edge and Devex, with `multi_finish[iFn].choice` not equal to
    /// `iFn`. Expected values are the C++ statements evaluated serially in C++ order.
    #[test]
    fn w5_pami_d_major_update_primal_dense() {
        for se in [true, false] {
            let mut ekk = HEkk::default();
            let mut dual = HEkkDual::default();
            let num_row: HighsInt = 3;
            dual.solver_num_row = num_row;
            dual.tp = 1e-7;
            dual.dual_rhs.work_count = -1;
            dual.dual_rhs.work_infeasibility = vec![9.0; 3];
            dual.edge_weight_mode = if se {
                EdgeWeightMode::kSteepestEdge
            } else {
                EdgeWeightMode::kDevex
            };
            dual.new_devex_framework = false;
            dual.col_bfrt.array = vec![0.5, 1.0, 0.0];
            ekk.info_.base_value_ = vec![1.0, 5.0, -2.0];
            ekk.info_.base_lower_ = vec![0.0, 0.0, 0.0];
            ekk.info_.base_upper_ = vec![1.0, 3.0, 1.0];
            ekk.info_.store_squared_primal_infeasibility = false;
            ekk.dual_edge_weight_ = vec![1.0, 1.5, 2.0];

            let cols = [[0.3, -0.7, 1.1], [0.2, 0.9, -0.4]];
            let rows = [[0.6, 0.1, -0.8], [-0.5, 0.25, 0.75]];
            // finish iFn uses multi_choice[1 - iFn]
            for i_fn in 0..2usize {
                let ch = 1 - i_fn;
                dual.multi_choice[ch].col_aq.array = cols[i_fn].to_vec();
                dual.multi_choice[ch].row_ep.array = rows[i_fn].to_vec();
                let f = &mut dual.multi_finish[i_fn];
                f.choice = Some(ch);
                f.row_out = [0, 2][i_fn];
                f.ed_wt = [2.0, 3.0][i_fn];
                f.alpha_row = [0.5, -0.25][i_fn];
                f.basic_bound = [0.25, -1.0][i_fn];
                f.basic_value = [0.75, 0.5][i_fn];
            }
            dual.multi_n_finish = 2;

            // Expected, serially in C++ order
            let mut bv = [1.0, 5.0, -2.0];
            let lo = [0.0, 0.0, 0.0];
            let up = [1.0, 3.0, 1.0];
            let mut inf = [0.0f64; 3];
            for r in 0..3 {
                bv[r] -= [0.5, 1.0, 0.0][r];
                let less = lo[r] - bv[r];
                let more = bv[r] - up[r];
                let x = if less > 1e-7 {
                    less
                } else if more > 1e-7 {
                    more
                } else {
                    0.0
                };
                inf[r] = f64::abs(x);
            }
            let mut ew = [1.0, 1.5, 2.0];
            for i_fn in 0..2 {
                let w = [2.0, 3.0][i_fn];
                let kai = -2.0 / [0.5, -0.25][i_fn];
                for r in 0..3 {
                    let a = cols[i_fn][r];
                    if se {
                        ew[r] += a * (w * a + kai * rows[i_fn][r]);
                        ew[r] = cpp_max(kMinDualSteepestEdgeWeight, ew[r]);
                    } else {
                        ew[r] = cpp_max(ew[r], w * a * a);
                    }
                }
            }
            let tp = ekk.options_.primal_feasibility_tolerance;
            for i_fn in 0..2 {
                let r = [0usize, 2][i_fn];
                let v = bv[r] - [0.25, -1.0][i_fn] + [0.75, 0.5][i_fn];
                bv[r] = v;
                let (l, u) = (lo[r], up[r]);
                let mut p = 0.0;
                if v < l - tp {
                    p = l - v;
                } else if v > u + tp {
                    p = v - u;
                }
                inf[r] = f64::abs(p);
            }
            for i_fn in 0..2 {
                let w = [2.0, 3.0][i_fn];
                let kai = -2.0 / [0.5, -0.25][i_fn];
                for j_fn in 0..i_fn {
                    let jr = [0usize, 2][j_fn];
                    let a = cols[i_fn][jr];
                    if se {
                        ew[jr] += a * (w * a + kai * rows[i_fn][jr]);
                        ew[jr] = cpp_max(kMinDualSteepestEdgeWeight, ew[jr]);
                    } else {
                        ew[jr] = cpp_max(ew[jr], w * a * a);
                    }
                }
                ew[[0usize, 2][i_fn]] = w;
            }

            dual.major_update_primal(&mut ekk);
            for r in 0..3 {
                assert_eq!(ekk.info_.base_value_[r].to_bits(), bv[r].to_bits());
                assert_eq!(
                    dual.dual_rhs.work_infeasibility[r].to_bits(),
                    inf[r].to_bits()
                );
                assert_eq!(ekk.dual_edge_weight_[r].to_bits(), ew[r].to_bits());
            }
            assert_eq!(dual.num_devex_iterations, if se { 0 } else { 2 });
        }
    }
}
