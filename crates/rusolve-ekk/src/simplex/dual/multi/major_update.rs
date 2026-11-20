//! HiGHS source: simplex/HEkkDualMulti.cpp (lines 535-763)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! PAMI (parallel multiple-pricing dual, `simplex_strategy = kSimplexStrategyDualMulti`):
//! `majorUpdate` and the major FTRANs (`majorUpdateFtranPrepare`, `majorUpdateFtranParallel` with its `for_each`s at 650, `majorUpdateFtranFinal` with 713 and 726).

use super::super::*;
use crate::prelude::*;

/// Raw `f64` pointer wrapper for the `for_each` bodies of `majorUpdateFtranFinal`. Each
/// task writes a disjoint index range of one dense vector (C++ `double*` capture).
#[derive(Clone, Copy)]
struct F64Ptr(*mut f64);
// SAFETY: The pointer is used in parallel for_each loops where each task writes to a
// disjoint range of indices of one vector, and reads only the other, unmodified vector.
unsafe impl Send for F64Ptr {}
// SAFETY: As above.
unsafe impl Sync for F64Ptr {}

impl F64Ptr {
    /// Returns the wrapped pointer; a method call captures the whole `F64Ptr` in closures.
    #[inline]
    fn get(self) -> *mut f64 {
        self.0
    }
}

/// Raw `HVector` pointer wrapper for the `multi_vector[]` list of `majorUpdateFtranParallel`
/// (C++ `HVector_ptr`). Each task owns exactly one vector of the list.
#[derive(Clone, Copy)]
struct HVectorRawPtr(*mut HVector);
// SAFETY: Every entry of the list is a distinct vector, and each task of the `for_each`
// takes its own index, so no two tasks touch the same `HVector`.
unsafe impl Send for HVectorRawPtr {}
// SAFETY: As above.
unsafe impl Sync for HVectorRawPtr {}

impl HVectorRawPtr {
    /// Returns the wrapped pointer; a method call captures the whole wrapper in closures.
    #[inline]
    fn get(self) -> *mut HVector {
        self.0
    }
}

impl HEkkDual {
    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorUpdate (lines 535-572, 38 lines)
    pub fn major_update(&mut self, ekk: &mut HEkk) {
        // 0. See if it's ready to perform a major update
        if self.rebuild_reason != 0 {
            self.multi_choose_again = 1;
        }
        if self.multi_choose_again == 0 {
            return;
        }

        // Major update - FTRANs
        self.major_update_ftran_prepare(ekk);
        self.major_update_ftran_parallel(ekk);
        self.major_update_ftran_final(ekk);

        // Major update - check for roll back
        for i_fn in 0..self.multi_n_finish {
            let i_finish = i_fn as usize;
            let ch = self.multi_finish[i_finish]
                .choice
                .expect("multi_finish choice set");
            let i_row_out = self.multi_finish[i_finish].row_out;
            let alpha_row = self.multi_finish[i_finish].alpha_row;
            let alpha_col = self.multi_choice[ch].col_aq.array[i_row_out as usize];

            // Use the two pivot values to identify numerical trouble
            if ekk.reinvert_on_numerical_trouble(
                "HEkkDual::majorUpdate",
                &mut self.numerical_trouble,
                alpha_col,
                alpha_row,
                kMultiNumericalTroubleTolerance,
            ) {
                self.rebuild_reason = kRebuildReasonPossiblySingularBasis;
                self.major_rollback(ekk);
                return;
            }
        }

        // Major update - primal and factor
        self.major_update_primal(ekk);
        self.major_update_factor(ekk);
        if self.new_devex_framework {
            self.initialise_devex_framework(ekk);
        }
        self.iteration_analysis_major(ekk);
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorUpdateFtranPrepare (lines 574-608, 35 lines)
    pub fn major_update_ftran_prepare(&mut self, ekk: &mut HEkk) {
        // Prepare FTRAN BFRT buffer
        self.col_bfrt.clear();
        for i_fn in 0..self.multi_n_finish {
            let ch = self.multi_finish[i_fn as usize]
                .choice
                .expect("multi_finish choice set");
            let variable_in = self.multi_finish[i_fn as usize].variable_in;
            let theta_primal = self.multi_finish[i_fn as usize].theta_primal;
            let mc: *mut MChoice = self.multi_choice.as_mut_ptr();
            // SAFETY: `ch` and every `jch` below are `multi_choice` indices set by
            // majorChooseRow (HEkkDualMulti.cpp 288-290); `jch != ch` because each minor
            // iteration owns its own choice. No other reference into `multi_choice` is live.
            let vec: &mut HVector = unsafe { &mut (*mc.add(ch)).col_bfrt };
            ekk.lp_.a_matrix_.collect_aj(vec, variable_in, theta_primal);

            // Update this buffer by previous Row_ep
            for j_fn in (0..i_fn).rev() {
                let j_finish = &self.multi_finish[j_fn as usize];
                let j_ch = j_finish.choice.expect("multi_finish choice set");
                let j_variable_in = j_finish.variable_in;
                let j_variable_out = j_finish.variable_out;
                let j_alpha_row = j_finish.alpha_row;
                // SAFETY: as above; `j_ch != ch`, so `row_ep` of `j_ch` and `col_bfrt` of
                // `ch` are disjoint fields of distinct choices.
                let j_row_ep_array: &Vec<f64> = unsafe { &(*mc.add(j_ch)).row_ep.array };
                let mut pivot_x: f64 = 0.0;
                for k in 0..(vec.count.max(0) as usize) {
                    let i_row = vec.index[k] as usize;
                    pivot_x += vec.array[i_row] * j_row_ep_array[i_row];
                }
                if pivot_x.abs() > kHighsTiny {
                    pivot_x /= j_alpha_row;
                    ekk.lp_.a_matrix_.collect_aj(vec, j_variable_in, -pivot_x);
                    ekk.lp_.a_matrix_.collect_aj(vec, j_variable_out, pivot_x);
                }
            }
            self.col_bfrt.saxpy(1.0, vec);
        }

        // Prepare regular FTRAN buffer
        for i_fn in 0..self.multi_n_finish {
            let ch = self.multi_finish[i_fn as usize]
                .choice
                .expect("multi_finish choice set");
            let variable_in = self.multi_finish[i_fn as usize].variable_in;
            let i_column = &mut self.multi_choice[ch].col_aq;
            i_column.clear();
            i_column.pack_flag = true;
            ekk.lp_.a_matrix_.collect_aj(i_column, variable_in, 1.0);
        }
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorUpdateFtranParallel (lines 610-691, 82 lines)
    pub fn major_update_ftran_parallel(&mut self, ekk: &mut HEkk) {
        ekk.analysis_.simplex_timer_start(FtranMixParClock, 0);

        // Prepare buffers
        let mut multi_ntasks: HighsInt = 0;
        let mut multi_density = [0.0f64; (kSimplexConcurrencyLimit * 2 + 1) as usize];
        let mut multi_vector: [HVectorRawPtr; (kSimplexConcurrencyLimit * 2 + 1) as usize] =
            [HVectorRawPtr(std::ptr::null_mut()); (kSimplexConcurrencyLimit * 2 + 1) as usize];
        // BFRT first
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_.operation_record_before_count(
                kSimplexNlaFtranBfrt,
                self.col_bfrt.count,
                ekk.info_.col_aq_density,
            );
        }
        multi_density[multi_ntasks as usize] = ekk.info_.col_aq_density;
        multi_vector[multi_ntasks as usize] = HVectorRawPtr(&mut self.col_bfrt as *mut HVector);
        multi_ntasks += 1;
        if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
            // Then DSE
            for i_fn in 0..self.multi_n_finish {
                let ch = self.multi_finish[i_fn as usize]
                    .choice
                    .expect("multi_finish choice set");
                if ekk.analysis_.analyse_simplex_summary_data {
                    ekk.analysis_.operation_record_before_count(
                        kSimplexNlaFtranDse,
                        self.multi_choice[ch].row_ep.count,
                        ekk.info_.row_dse_density,
                    );
                }
                multi_density[multi_ntasks as usize] = ekk.info_.row_dse_density;
                multi_vector[multi_ntasks as usize] =
                    HVectorRawPtr(&mut self.multi_choice[ch].row_ep as *mut HVector);
                multi_ntasks += 1;
            }
        }
        // Then Column
        for i_fn in 0..self.multi_n_finish {
            let ch = self.multi_finish[i_fn as usize]
                .choice
                .expect("multi_finish choice set");
            if ekk.analysis_.analyse_simplex_summary_data {
                ekk.analysis_.operation_record_before_count(
                    kSimplexNlaFtran,
                    self.multi_choice[ch].col_aq.count,
                    ekk.info_.col_aq_density,
                );
            }
            multi_density[multi_ntasks as usize] = ekk.info_.col_aq_density;
            multi_vector[multi_ntasks as usize] =
                HVectorRawPtr(&mut self.multi_choice[ch].col_aq as *mut HVector);
            multi_ntasks += 1;
        }

        // Perform FTRAN
        // #pragma omp parallel for schedule(dynamic, 1)
        // printf("majorUpdateFtranParallel: starting %d tasks\n", multi_ntasks);
        {
            let refs = NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_);
            let simplex_nla = &ekk.simplex_nla_;
            let analysis = &ekk.analysis_;
            crate::parallel::highs_parallel::for_each(
                0,
                multi_ntasks,
                &|start: HighsInt, end: HighsInt| {
                    for i in start..end {
                        let rhs: *mut HVector = multi_vector[i as usize].get();
                        let density = multi_density[i as usize];
                        // C++ `analysis->getThreadFactorTimerClockPointer()`; the HighsSimplexAnalysis
                        // port keeps no per-thread clock, so thread 0 is used.
                        let factor_timer_clock_pointer =
                            analysis.get_thread_factor_timer_clock_ptr(0);
                        // SAFETY: each `i` names one vector of `multi_vector`; the vectors are
                        // distinct, so tasks never alias.
                        simplex_nla.ftran(
                            &refs,
                            unsafe { &mut *rhs },
                            density,
                            factor_timer_clock_pointer,
                        );
                    }
                },
                1,
            );
        }

        // Update ticks
        for i_fn in 0..self.multi_n_finish {
            let ch = self.multi_finish[i_fn as usize]
                .choice
                .expect("multi_finish choice set");
            ekk.total_synthetic_tick_ += self.multi_choice[ch].col_aq.synthetic_tick;
            ekk.total_synthetic_tick_ += self.multi_choice[ch].row_ep.synthetic_tick;
            // reporting-only work counter, never read by the search.
            crate::work::add_simplex_ticks(
                self.multi_choice[ch].col_aq.synthetic_tick
                    + self.multi_choice[ch].row_ep.synthetic_tick,
            );
        }

        // Update rates
        if ekk.analysis_.analyse_simplex_summary_data {
            ekk.analysis_
                .operation_record_after_count(kSimplexNlaFtranBfrt, self.col_bfrt.count);
        }
        for i_fn in 0..self.multi_n_finish {
            let ch = self.multi_finish[i_fn as usize]
                .choice
                .expect("multi_finish choice set");
            let local_col_aq_density =
                self.multi_choice[ch].col_aq.count as f64 / self.solver_num_row as f64;
            let mut col_aq_density = ekk.info_.col_aq_density;
            ekk.update_operation_result_density(local_col_aq_density, &mut col_aq_density);
            ekk.info_.col_aq_density = col_aq_density;
            if ekk.analysis_.analyse_simplex_summary_data {
                ekk.analysis_.operation_record_after_count(
                    kSimplexNlaFtran,
                    self.multi_choice[ch].col_aq.count,
                );
            }
            if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                let local_row_dse_density =
                    self.multi_choice[ch].row_ep.count as f64 / self.solver_num_row as f64;
                let mut row_dse_density = ekk.info_.row_dse_density;
                ekk.update_operation_result_density(local_row_dse_density, &mut row_dse_density);
                ekk.info_.row_dse_density = row_dse_density;
                if ekk.analysis_.analyse_simplex_summary_data {
                    ekk.analysis_.operation_record_after_count(
                        kSimplexNlaFtranDse,
                        self.multi_choice[ch].row_ep.count,
                    );
                }
            }
        }
        ekk.analysis_.simplex_timer_stop(FtranMixParClock, 0);
    }

    /// HiGHS: HEkkDualMulti.cpp::HEkkDual::majorUpdateFtranFinal (lines 693-763, 71 lines)
    pub fn major_update_ftran_final(&mut self, ekk: &mut HEkk) {
        ekk.analysis_.simplex_timer_start(FtranMixFinalClock, 0);
        let update_ftran_in_dense = self.dual_rhs.work_count < 0;
        if update_ftran_in_dense {
            for i_fn in 0..self.multi_n_finish {
                let ch = self.multi_finish[i_fn as usize]
                    .choice
                    .expect("multi_finish choice set");
                let mc: *mut MChoice = self.multi_choice.as_mut_ptr();
                // SAFETY: `ch` is a `multi_choice` index; no reference into `multi_choice`
                // is live while the raw pointers below are used.
                unsafe {
                    (*mc.add(ch)).col_aq.count = -1;
                    (*mc.add(ch)).row_ep.count = -1;
                }
                let my_col = F64Ptr(unsafe { (*mc.add(ch)).col_aq.array.as_mut_ptr() });
                let my_row = F64Ptr(unsafe { (*mc.add(ch)).row_ep.array.as_mut_ptr() });
                for j_fn in 0..i_fn {
                    let j_finish = &self.multi_finish[j_fn as usize];
                    let j_ch = j_finish.choice.expect("multi_finish choice set");
                    let pivot_row = j_finish.row_out;
                    let pivot_alpha = j_finish.alpha_row;
                    // SAFETY: `j_ch != ch` (distinct minor iterations), so this reads a vector
                    // that the writes below never touch.
                    let pivot_array = F64Ptr(unsafe { (*mc.add(j_ch)).col_aq.array.as_mut_ptr() });
                    // SAFETY: `pivot_row` is a row index of the vectors.
                    let pivot_x1 = unsafe { *my_col.get().add(pivot_row as usize) };
                    let pivot_x2 = unsafe { *my_row.get().add(pivot_row as usize) };

                    // The FTRAN regular buffer
                    if pivot_x1.abs() > kHighsTiny {
                        let pivot = pivot_x1 / pivot_alpha;
                        // #pragma omp parallel for
                        crate::parallel::highs_parallel::for_each(
                            0,
                            self.solver_num_row,
                            &|start: HighsInt, end: HighsInt| {
                                for i in start..end {
                                    // SAFETY: each task writes a disjoint index range of `myCol`.
                                    unsafe {
                                        *my_col.get().add(i as usize) -=
                                            pivot * *pivot_array.get().add(i as usize);
                                    }
                                }
                            },
                            100,
                        );
                        // SAFETY: `pivot_row` is in bounds; the for_each has completed.
                        unsafe {
                            *my_col.get().add(pivot_row as usize) = pivot;
                        }
                    }
                    // The FTRAN-DSE buffer
                    if pivot_x2.abs() > kHighsTiny {
                        let pivot = pivot_x2 / pivot_alpha;
                        // #pragma omp parallel for
                        crate::parallel::highs_parallel::for_each(
                            0,
                            self.solver_num_row,
                            &|start: HighsInt, end: HighsInt| {
                                for i in start..end {
                                    // SAFETY: each task writes a disjoint index range of `myRow`.
                                    unsafe {
                                        *my_row.get().add(i as usize) -=
                                            pivot * *pivot_array.get().add(i as usize);
                                    }
                                }
                            },
                            100,
                        );
                        // SAFETY: `pivot_row` is in bounds; the for_each has completed.
                        unsafe {
                            *my_row.get().add(pivot_row as usize) = pivot;
                        }
                    }
                }
            }
        } else {
            for i_fn in 0..self.multi_n_finish {
                let ch = self.multi_finish[i_fn as usize]
                    .choice
                    .expect("multi_finish choice set");
                let mc: *mut MChoice = self.multi_choice.as_mut_ptr();
                for j_fn in 0..i_fn {
                    let j_finish = &self.multi_finish[j_fn as usize];
                    let j_ch = j_finish.choice.expect("multi_finish choice set");
                    let pivot_row = j_finish.row_out;
                    let j_alpha_row = j_finish.alpha_row;
                    // SAFETY: `ch` and `j_ch` are distinct `multi_choice` indices; the
                    // references below are to distinct vectors and are not held across iterations.
                    let col: &mut HVector = unsafe { &mut (*mc.add(ch)).col_aq };
                    let row: &mut HVector = unsafe { &mut (*mc.add(ch)).row_ep };
                    let j_col: &HVector = unsafe { &(*mc.add(j_ch)).col_aq };

                    // The FTRAN regular buffer
                    let mut pivot_x1 = col.array[pivot_row as usize];
                    if pivot_x1.abs() > kHighsTiny {
                        pivot_x1 /= j_alpha_row;
                        col.saxpy(-pivot_x1, j_col);
                        col.array[pivot_row as usize] = pivot_x1;
                    }
                    // The FTRAN-DSE buffer
                    let mut pivot_x2 = row.array[pivot_row as usize];
                    if pivot_x2.abs() > kHighsTiny {
                        pivot_x2 /= j_alpha_row;
                        row.saxpy(-pivot_x2, j_col);
                        row.array[pivot_row as usize] = pivot_x2;
                    }
                }
            }
        }
        ekk.analysis_.simplex_timer_stop(FtranMixFinalClock, 0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Unit tests
    //
    // Expected values are the C++ statements evaluated serially in C++ order on plain
    // arrays (`RefVec`: dense array plus index list, as `HVector`). majorUpdateFtranParallel
    // needs an INVERT for `simplex_nla_.ftran` and is covered end to end by PAMI-IT.

    /// Dense array plus index list, mirroring the `HVector` fields the C++ touches.
    #[derive(Clone)]
    struct RefVec {
        a: [f64; 3],
        idx: Vec<usize>,
    }

    impl RefVec {
        fn new(entries: &[(usize, f64)]) -> Self {
            let mut v = RefVec {
                a: [0.0; 3],
                idx: Vec::new(),
            };
            for &(r, x) in entries {
                v.a[r] = x;
                v.idx.push(r);
            }
            v
        }
    }

    /// C++ HighsSparseMatrix::collectAj (HighsSparseMatrix.cpp 1618-1637) on a 3x2 matrix.
    fn collect_ref(v: &mut RefVec, use_col: usize, multiplier: f64) {
        let cols: [&[(usize, f64)]; 2] = [&[(0, 1.0), (2, -2.0)], &[(1, 3.0), (2, 0.5)]];
        let entries: Vec<(usize, f64)> = if use_col < 2 {
            cols[use_col].to_vec()
        } else {
            vec![(use_col - 2, 1.0)]
        };
        for (r, val) in entries {
            let value0 = v.a[r];
            let value1 = if use_col < 2 {
                value0 + multiplier * val
            } else {
                value0 + multiplier
            };
            if value0 == 0.0 {
                v.idx.push(r);
            }
            v.a[r] = if value1.abs() < kHighsTiny {
                kHighsZero
            } else {
                value1
            };
        }
    }

    /// C++ HVectorBase::saxpy (HVectorBase.cpp 176-200).
    fn saxpy_ref(x: &mut RefVec, pivot_x: f64, p: &RefVec) {
        for &r in &p.idx {
            let x0 = x.a[r];
            let x1 = x0 + pivot_x * p.a[r];
            if x0 == 0.0 {
                x.idx.push(r);
            }
            x.a[r] = if x1.abs() < kHighsTiny {
                kHighsZero
            } else {
                x1
            };
        }
    }

    fn to_hvector(v: &RefVec) -> HVector {
        let mut h = HVector::default();
        h.setup(3);
        for &r in &v.idx {
            h.index[h.count as usize] = r as HighsInt;
            h.count += 1;
            h.array[r] = v.a[r];
        }
        h
    }

    fn assert_matches(h: &HVector, v: &RefVec) {
        assert_eq!(h.count as usize, v.idx.len());
        for (k, &r) in v.idx.iter().enumerate() {
            assert_eq!(h.index[k] as usize, r);
        }
        for r in 0..3 {
            assert_eq!(h.array[r].to_bits(), v.a[r].to_bits(), "row {r}");
        }
    }

    fn small_lp() -> HEkk {
        let mut ekk = HEkk::default();
        let a = &mut ekk.lp_.a_matrix_;
        a.num_col_ = 2;
        a.num_row_ = 3;
        a.start_ = vec![0, 2, 4];
        a.index_ = vec![0, 2, 1, 2];
        a.value_ = vec![1.0, -2.0, 3.0, 0.5];
        ekk
    }

    #[test]
    fn filled_matches_cpp_line_by_line() {
        // Unit checks of majorUpdate, majorUpdateFtranPrepare, majorUpdateFtranParallel, majorUpdateFtranFinal on a small LP with
        // simplex_strategy = 3 (PAMI); the end-to-end bit-identity check is PAMI-IT.

        // majorUpdate 539-540: no rebuild and no choose-again returns before the FTRANs
        // (multi_finish[0].choice is None, so reaching them would panic).
        let mut ekk = small_lp();
        let mut dual = HEkkDual::default();
        dual.rebuild_reason = 0;
        dual.multi_choose_again = 0;
        dual.multi_n_finish = 1;
        dual.major_update(&mut ekk);
        assert_eq!(dual.multi_choose_again, 0);
        assert_eq!(dual.rebuild_reason, 0);
    }

    /// majorUpdateFtranPrepare (HEkkDualMulti.cpp 574-608), with `multi_finish[iFn].choice`
    /// not equal to `iFn` and a `col_BFRT` that already holds an entry (the C++ does not
    /// clear the per-choice buffer here).
    #[test]
    fn w5_pami_c_major_update_ftran_prepare() {
        let mut ekk = small_lp();
        let mut dual = HEkkDual::default();
        dual.solver_num_row = 3;
        // Stale global BFRT buffer, cleared by line 576
        dual.col_bfrt = to_hvector(&RefVec::new(&[(0, 9.0)]));
        let bfrt1 = RefVec::new(&[(1, 0.25)]);
        let bfrt0 = RefVec::new(&[]);
        dual.multi_choice[1].col_bfrt = to_hvector(&bfrt1);
        dual.multi_choice[0].col_bfrt = to_hvector(&bfrt0);
        let row_ep1 = [0.2, 0.7, -0.3];
        dual.multi_choice[1].row_ep = to_hvector(&RefVec::new(&[(0, 0.2), (1, 0.7), (2, -0.3)]));
        // Stale column buffers, cleared by line 603
        dual.multi_choice[1].col_aq = to_hvector(&RefVec::new(&[(2, 4.0)]));
        dual.multi_choice[0].col_aq = to_hvector(&RefVec::new(&[(0, -1.0), (1, 2.0)]));
        {
            let f = &mut dual.multi_finish[0];
            f.choice = Some(1);
            f.variable_in = 0;
            f.variable_out = 3;
            f.theta_primal = 0.5;
            f.alpha_row = 0.8;
        }
        {
            let f = &mut dual.multi_finish[1];
            f.choice = Some(0);
            f.variable_in = 1;
            f.variable_out = 4;
            f.theta_primal = -1.5;
            f.alpha_row = -0.4;
        }
        dual.multi_n_finish = 2;

        // Expected, serially in C++ order
        let mut bfrt = RefVec::new(&[]);
        let mut v1 = bfrt1.clone();
        collect_ref(&mut v1, 0, 0.5);
        saxpy_ref(&mut bfrt, 1.0, &v1);
        let mut v0 = bfrt0.clone();
        collect_ref(&mut v0, 1, -1.5);
        let mut pivot_x: f64 = 0.0;
        for &r in &v0.idx {
            pivot_x += v0.a[r] * row_ep1[r];
        }
        assert!(pivot_x.abs() > kHighsTiny);
        pivot_x /= 0.8;
        collect_ref(&mut v0, 0, -pivot_x);
        collect_ref(&mut v0, 3, pivot_x);
        saxpy_ref(&mut bfrt, 1.0, &v0);
        let mut aq1 = RefVec::new(&[]);
        collect_ref(&mut aq1, 0, 1.0);
        let mut aq0 = RefVec::new(&[]);
        collect_ref(&mut aq0, 1, 1.0);

        dual.major_update_ftran_prepare(&mut ekk);
        assert_matches(&dual.multi_choice[1].col_bfrt, &v1);
        assert_matches(&dual.multi_choice[0].col_bfrt, &v0);
        assert_matches(&dual.col_bfrt, &bfrt);
        assert_matches(&dual.multi_choice[1].col_aq, &aq1);
        assert_matches(&dual.multi_choice[0].col_aq, &aq0);
        assert!(dual.multi_choice[1].col_aq.pack_flag);
        assert!(dual.multi_choice[0].col_aq.pack_flag);
    }

    /// Three finishes for majorUpdateFtranFinal: finish iFn uses `multi_choice[CH[iFn]]`.
    /// Finish 2 skips j = 0 in the column buffer and j = 1 in the DSE buffer
    /// (`fabs(pivotX) > kHighsTiny` false), and its j = 1 update reads the column of
    /// finish 1 as updated at iFn = 1.
    const CH: [usize; 3] = [1, 0, 2];
    const ROW_OUT: [usize; 3] = [0, 2, 1];
    const ALPHA: [f64; 3] = [0.5, -0.25, 2.0];

    fn final_cols() -> [RefVec; 3] {
        // indexed by choice
        [
            RefVec::new(&[(0, 0.4), (2, 0.9)]),
            RefVec::new(&[(0, 0.5), (1, 0.3)]),
            RefVec::new(&[(1, 1.25), (2, -0.7)]),
        ]
    }

    fn final_rows() -> [RefVec; 3] {
        [
            RefVec::new(&[(0, -0.6), (1, 0.2)]),
            RefVec::new(&[(0, 0.1), (2, 0.3)]),
            RefVec::new(&[(0, 1.5)]),
        ]
    }

    fn final_setup(work_count: HighsInt) -> HEkkDual {
        let mut dual = HEkkDual::default();
        dual.solver_num_row = 3;
        dual.dual_rhs.work_count = work_count;
        let cols = final_cols();
        let rows = final_rows();
        for ch in 0..3 {
            dual.multi_choice[ch].col_aq = to_hvector(&cols[ch]);
            dual.multi_choice[ch].row_ep = to_hvector(&rows[ch]);
        }
        for i_fn in 0..3 {
            let f = &mut dual.multi_finish[i_fn];
            f.choice = Some(CH[i_fn]);
            f.row_out = ROW_OUT[i_fn] as HighsInt;
            f.alpha_row = ALPHA[i_fn];
        }
        dual.multi_n_finish = 3;
        dual
    }

    /// majorUpdateFtranFinal dense branch (HEkkDualMulti.cpp 696-736).
    #[test]
    fn w5_pami_c_major_update_ftran_final_dense() {
        let mut ekk = HEkk::default();
        let mut dual = final_setup(-1);

        // Expected, serially in C++ order
        let mut cols = final_cols().map(|v| v.a);
        let mut rows = final_rows().map(|v| v.a);
        for i_fn in 0..3 {
            for j_fn in 0..i_fn {
                let pivot_row = ROW_OUT[j_fn];
                let pivot_alpha = ALPHA[j_fn];
                let pivot_array = cols[CH[j_fn]];
                let pivot_x1 = cols[CH[i_fn]][pivot_row];
                let pivot_x2 = rows[CH[i_fn]][pivot_row];
                if pivot_x1.abs() > kHighsTiny {
                    let pivot = pivot_x1 / pivot_alpha;
                    for i in 0..3 {
                        cols[CH[i_fn]][i] -= pivot * pivot_array[i];
                    }
                    cols[CH[i_fn]][pivot_row] = pivot;
                }
                if pivot_x2.abs() > kHighsTiny {
                    let pivot = pivot_x2 / pivot_alpha;
                    for i in 0..3 {
                        rows[CH[i_fn]][i] -= pivot * pivot_array[i];
                    }
                    rows[CH[i_fn]][pivot_row] = pivot;
                }
            }
        }

        dual.major_update_ftran_final(&mut ekk);
        for ch in 0..3 {
            assert_eq!(dual.multi_choice[ch].col_aq.count, -1);
            assert_eq!(dual.multi_choice[ch].row_ep.count, -1);
            for r in 0..3 {
                assert_eq!(
                    dual.multi_choice[ch].col_aq.array[r].to_bits(),
                    cols[ch][r].to_bits()
                );
                assert_eq!(
                    dual.multi_choice[ch].row_ep.array[r].to_bits(),
                    rows[ch][r].to_bits()
                );
            }
        }
        // Finish 2 took the j = 1 column update and the j = 0 DSE update only
        assert_eq!(cols[2][2].to_bits(), (-0.7f64 / -0.25).to_bits());
        assert_eq!(rows[2][0].to_bits(), (1.5f64 / 0.5).to_bits());
    }

    /// majorUpdateFtranFinal sparse branch (HEkkDualMulti.cpp 737-761).
    #[test]
    fn w5_pami_c_major_update_ftran_final_sparse() {
        let mut ekk = HEkk::default();
        let mut dual = final_setup(0);

        // Expected, serially in C++ order
        let mut cols = final_cols();
        let mut rows = final_rows();
        for i_fn in 0..3 {
            for j_fn in 0..i_fn {
                let pivot_row = ROW_OUT[j_fn];
                let j_col = cols[CH[j_fn]].clone();
                let mut pivot_x1 = cols[CH[i_fn]].a[pivot_row];
                if pivot_x1.abs() > kHighsTiny {
                    pivot_x1 /= ALPHA[j_fn];
                    saxpy_ref(&mut cols[CH[i_fn]], -pivot_x1, &j_col);
                    cols[CH[i_fn]].a[pivot_row] = pivot_x1;
                }
                let mut pivot_x2 = rows[CH[i_fn]].a[pivot_row];
                if pivot_x2.abs() > kHighsTiny {
                    pivot_x2 /= ALPHA[j_fn];
                    saxpy_ref(&mut rows[CH[i_fn]], -pivot_x2, &j_col);
                    rows[CH[i_fn]].a[pivot_row] = pivot_x2;
                }
            }
        }

        dual.major_update_ftran_final(&mut ekk);
        for ch in 0..3 {
            assert_matches(&dual.multi_choice[ch].col_aq, &cols[ch]);
            assert_matches(&dual.multi_choice[ch].row_ep, &rows[ch]);
        }
    }
}
