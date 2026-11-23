//! HiGHS source: ipm/hipo/factorhighs/DgemmParallel.h, ipm/hipo/factorhighs/DgemmParallel.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Blocked Schur-complement update `Q -= P * R^T`. Only reached with
//! `parnode == true` (`DenseFactHybrid.cpp` lines 175, 208), which a serial run never
//! sets; ported anyway. The C++ `highs::parallel::for_each(0, row,
//! ..., kBlockGrainSize * jb)` is the real `highs_parallel::for_each`;
//! at one worker it visits the leaf ranges of the HiGHS binary split in
//! increasing `start` order, as the serial loop does, and every entry of `Q` gets the same
//! operations whatever the split. `P_` / `R_` are slices starting at the C
//! pointers; `Q_` is a raw slice pointer (C++ `double* Q_`) so
//! that concurrent `run` calls each reborrow only their own rows of `Q`.

use super::super::Int;
use super::call_and_time_blas::callAndTime_dgemm;
use super::data_collector::DataCollector;
use super::{kBlockGrainSize, kBlockParallelThreshold};
use std::marker::PhantomData;

/// HiGHS: DgemmParallel.h `class dgemmParalleliser`
pub struct dgemmParalleliser<'a> {
    pub P_: &'a [f64],
    pub R_: &'a [f64],
    /// C++ `double* Q_`: the whole `Q` slice passed to `new`, written through by
    /// the `const` `run`.
    pub Q_: *mut [f64],
    pub col_: Int,
    pub jb_: Int,
    pub data_: &'a DataCollector,
    _Q: PhantomData<&'a mut [f64]>,
}

// SAFETY: `run` is the only access to `Q_`, and each call writes only rows
// `[start, end)` of `Q` (`Q_[start * col_ .. end * col_]`); `for_each` hands
// disjoint `[start, end)` ranges to its tasks, so concurrent `run` calls touch
// disjoint parts of `Q` (C++ tasks sharing `pt` and calling the `const` run).
// `P_`, `R_` and `data_` are shared references.
unsafe impl Sync for dgemmParalleliser<'_> {}

impl<'a> dgemmParalleliser<'a> {
    /// HiGHS: DgemmParallel.cpp::dgemmParalleliser::dgemmParalleliser (lines 9-12)
    pub fn new(
        P: &'a [f64],
        R: &'a [f64],
        Q: &'a mut [f64],
        col: Int,
        jb: Int,
        data: &'a DataCollector,
    ) -> Self {
        Self {
            P_: P,
            R_: R,
            Q_: Q as *mut [f64],
            col_: col,
            jb_: jb,
            data_: data,
            _Q: PhantomData,
        }
    }

    /// HiGHS: DgemmParallel.cpp::dgemmParalleliser::run (lines 14-18)
    /// The C++ `const` method writes through `Q_`; the C argument
    /// `&Q_[start * col_]` is the slice of rows `[start, end)` (all that the
    /// `m = col_`, `n = end - start`, `ldc = col_` dgemm touches).
    pub fn run(&self, start: Int, end: Int, beta: f64) {
        let lo = (start * self.col_) as usize;
        let hi = (end * self.col_) as usize;
        assert!(lo <= hi && hi <= self.Q_.len());
        // SAFETY: `Q_` comes from the `&'a mut [f64]` given to `new` (alive and
        // otherwise unused for `'a`), `lo..hi` is in bounds, and concurrent `run`
        // calls get disjoint `[start, end)` from `for_each` (see `impl Sync`).
        let Q = unsafe { std::slice::from_raw_parts_mut((self.Q_ as *mut f64).add(lo), hi - lo) };
        callAndTime_dgemm(
            b'T',
            b'N',
            self.col_,
            end - start,
            self.jb_,
            -1.0,
            self.P_,
            self.jb_,
            &self.R_[(start * self.jb_) as usize..],
            self.jb_,
            beta,
            Q,
            self.col_,
            self.data_,
        );
    }
}

/// HiGHS: DgemmParallel.cpp::dgemmParallel (lines 20-36)
pub fn dgemmParallel(
    P: &[f64],
    R: &[f64],
    Q: &mut [f64],
    col: Int,
    jb: Int,
    row: Int,
    nb: Int,
    beta: f64,
    data: &DataCollector,
) {
    // if there is enough work to be done, parallelise
    if col >= nb / 2 && jb >= nb / 2 && row >= kBlockParallelThreshold * nb {
        let gemmP = dgemmParalleliser::new(P, R, Q, col, jb, data);
        let pt = &gemmP;

        // I need to use an object to call gemm, otherwise the task is too large and
        // static_assert in the parallel deque fails.
        crate::parallel::highs_parallel::for_each(
            0,
            row,
            &|start: Int, end: Int| pt.run(start, end, beta),
            kBlockGrainSize * jb,
        );
    } else {
        callAndTime_dgemm(
            b'T', b'N', col, row, jb, -1.0, P, jb, R, jb, beta, Q, col, data,
        );
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests: the parallel branch matches one serial dgemm.
    use super::*;

    #[test]
    fn leaves_increasing_and_cover_range() {
        let leaves = std::sync::Mutex::new(Vec::<(Int, Int)>::new());
        std::thread::scope(|s| {
            s.spawn(|| {
                crate::parallel::highs_parallel::initialize_scheduler(1);
                crate::parallel::highs_parallel::for_each(
                    0,
                    11,
                    &|s, e| {
                        leaves.lock().unwrap().push((s, e));
                    },
                    2,
                );
                crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            });
        });
        let leaves = leaves.into_inner().unwrap();
        let mut next = 0;
        for &(s, e) in &leaves {
            assert_eq!(s, next);
            assert!(e - s <= 2 && e > s);
            next = e;
        }
        assert_eq!(next, 11);
    }

    #[test]
    fn parallel_branch_equals_serial_gemm() {
        // col = 2, jb = 2, row = 11, nb = 2: 11 >= 5 * 2 takes the split branch.
        let (col, jb, row, nb) = (2, 2, 11, 2);
        let P: Vec<f64> = (0..(jb * col)).map(|k| 0.3 * k as f64 + 0.1).collect();
        let R: Vec<f64> = (0..(jb * row)).map(|k| 1.0 / (k as f64 + 1.0)).collect();
        let Q0: Vec<f64> = (0..(col * row)).map(|k| 0.7 * k as f64).collect();
        let d = DataCollector::new();
        let mut Q1 = Q0.clone();
        dgemmParallel(&P, &R, &mut Q1, col, jb, row, nb, 1.0, &d);
        let mut Q2 = Q0.clone();
        callAndTime_dgemm(
            b'T', b'N', col, row, jb, -1.0, &P, jb, &R, jb, 1.0, &mut Q2, col, &d,
        );
        assert_eq!(Q1, Q2);
    }

    // Unit tests.

    #[test]
    fn w5_hipo_dgemm_for_each_blocks() {
        let (col, jb, row, nb) = (2, 2, 11, 2);
        let P: Vec<f64> = (0..(jb * col)).map(|k| 0.3 * k as f64 + 0.1).collect();
        let R: Vec<f64> = (0..(jb * row)).map(|k| 1.0 / (k as f64 + 1.0)).collect();
        let Q0: Vec<f64> = (0..(col * row)).map(|k| 0.7 * k as f64).collect();
        let d = DataCollector::new();

        let mut Q_serial = Q0.clone();
        callAndTime_dgemm(
            b'T',
            b'N',
            col,
            row,
            jb,
            -1.0,
            &P,
            jb,
            &R,
            jb,
            1.0,
            &mut Q_serial,
            col,
            &d,
        );

        for &n_threads in &[1, 2, 4] {
            let P = P.clone();
            let R = R.clone();
            let Q0 = Q0.clone();
            let d = d.clone();
            let handle = std::thread::spawn(move || {
                crate::parallel::highs_parallel::initialize_scheduler(n_threads);
                let mut Q = Q0.clone();
                dgemmParallel(&P, &R, &mut Q, col, jb, row, nb, 1.0, &d);
                crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
                Q
            });
            let Q = handle.join().unwrap();
            assert_eq!(Q, Q_serial);
        }
    }
}
