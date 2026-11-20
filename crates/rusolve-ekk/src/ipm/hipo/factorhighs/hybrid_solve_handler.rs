//! HiGHS source: ipm/hipo/factorhighs/HybridSolveHandler.h, ipm/hipo/factorhighs/HybridSolveHandler.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Forward, diagonal and backward solves with the hybrid-format factor.

use super::super::{Int, Int64};
use super::call_and_time_blas::{callAndTime_dgemv, callAndTime_dtrsv};
use super::data_collector::DataCollector;
use super::format_handler::extra_space_frontal;
use super::solve_handler::{SolveHandler, SolveHandlerOps};
use super::swaps::permuteWithSwaps;
use super::symbolic::Symbolic;
use crate::util::highs_utils::cpp_min;

/// HiGHS: HybridSolveHandler.h `class HybridSolveHandler`
pub struct HybridSolveHandler<'a> {
    pub base: SolveHandler<'a>,
    pub swaps_: &'a [Vec<Int>],
    pub pivot_2x2_: &'a [Vec<f64>],
}

impl<'a> HybridSolveHandler<'a> {
    /// HiGHS: HybridSolveHandler.cpp::HybridSolveHandler::HybridSolveHandler (lines 12-16)
    pub fn new(
        S: &'a Symbolic,
        sn_columns: &'a [Vec<f64>],
        swaps: &'a [Vec<Int>],
        pivot_2x2: &'a [Vec<f64>],
        data: &'a DataCollector,
    ) -> Self {
        HybridSolveHandler {
            base: SolveHandler::new(S, sn_columns, data),
            swaps_: swaps,
            pivot_2x2_: pivot_2x2,
        }
    }
}

impl<'a> SolveHandlerOps<'a> for HybridSolveHandler<'a> {
    fn base(&self) -> &SolveHandler<'a> {
        &self.base
    }

    /// HiGHS: HybridSolveHandler.cpp::HybridSolveHandler::forwardSolve (lines 18-139)
    fn forwardSolve(&self, x: &mut Vec<f64>) {
        // HIPO_CLOCK_CREATE; (line 24) — level-0 no-op, dropped.

        let nb = self.base.S_.blockSize();

        for sn in 0..self.base.S_.sn() {
            // leading size of supernode
            let ldSn: Int = (self.base.S_.ptr(sn + 1) - self.base.S_.ptr(sn)) as Int;

            // number of columns in the supernode
            let sn_size = self.base.S_.snStart(sn + 1) - self.base.S_.snStart(sn);

            // first columns of the supernode
            let sn_start = self.base.S_.snStart(sn);

            // index to access S->rows for this supernode
            let start_row: Int64 = self.base.S_.ptr(sn);

            // number of blocks of columns
            let n_blocks = (sn_size - 1) / nb + 1;

            // index to access snColumns[sn]
            let mut SnCol_ind: Int64 = 0;

            if sn_size < nb {
                // Fast solve
                // If supernode is small, avoid making BLAS calls
                let jb = sn_size;
                let x_start = sn_start;

                // HIPO_CLOCK_START(2); (line 54) — level-0 no-op, dropped.
                // apply swaps to portion of rhs that is affected
                {
                    let current_swaps = &self.swaps_[sn as usize];
                    permuteWithSwaps(&mut x[x_start as usize..], current_swaps, jb, false);
                }
                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 58) — dropped.

                // HIPO_CLOCK_START(2); (line 61) — level-0 no-op, dropped.
                for row in 0..jb {
                    for col in 0..row {
                        x[x_start as usize + row as usize] -= self.base.sn_columns_[sn as usize]
                            [col as usize + (jb * row) as usize]
                            * x[x_start as usize + col as usize];
                    }
                }

                for row in jb..ldSn {
                    for col in 0..jb {
                        x[self.base.S_.rows(start_row + row as Int64) as usize] -=
                            self.base.sn_columns_[sn as usize][col as usize + (jb * row) as usize]
                                * x[x_start as usize + col as usize];
                    }
                }
                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense); (line 75) — dropped.

                // HIPO_CLOCK_START(2); (line 78) — level-0 no-op, dropped.
                // apply inverse swaps
                {
                    let current_swaps = &self.swaps_[sn as usize];
                    permuteWithSwaps(&mut x[x_start as usize..], current_swaps, jb, true);
                }
                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 81) — dropped.
            } else {
                // go through blocks of columns for this supernode
                for j in 0..n_blocks {
                    // number of columns in the block
                    let jb = cpp_min(nb, sn_size - nb * j);

                    // number of entries in diagonal part
                    let diag_entries = jb * jb;

                    // index to access vector x
                    let x_start = sn_start + nb * j;

                    // HIPO_CLOCK_START(2); (line 97) — level-0 no-op, dropped.
                    // apply swaps to portion of rhs that is affected
                    {
                        let current_swaps = &self.swaps_[sn as usize][(nb * j) as usize..];
                        permuteWithSwaps(&mut x[x_start as usize..], current_swaps, jb, false);
                    }
                    // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 101) — dropped.

                    // HIPO_CLOCK_START(2); (line 104) — level-0 no-op, dropped.
                    callAndTime_dtrsv(
                        b'U',
                        b'T',
                        b'U',
                        jb,
                        &self.base.sn_columns_[sn as usize][SnCol_ind as usize..],
                        jb,
                        &mut x[x_start as usize..],
                        1,
                        self.base.data_,
                    );

                    SnCol_ind += diag_entries as Int64;

                    // temporary space for gemv
                    let gemv_space: Int = ldSn - nb * j - jb;
                    let mut y = vec![0.0f64; gemv_space as usize];
                    if gemv_space > 0 {
                        callAndTime_dgemv(
                            b'T',
                            jb,
                            gemv_space,
                            1.0,
                            &self.base.sn_columns_[sn as usize][SnCol_ind as usize..],
                            jb,
                            &x[x_start as usize..],
                            1,
                            0.0,
                            &mut y,
                            1,
                            self.base.data_,
                        );

                        SnCol_ind += (jb * gemv_space) as Int64;
                        // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense); (line 119) — dropped.

                        // HIPO_CLOCK_START(2); (line 121) — level-0 no-op, dropped.
                        // scatter solution of gemv
                        for i in 0..gemv_space {
                            let row = self
                                .base
                                .S_
                                .rows(start_row + (nb * j) as Int64 + jb as Int64 + i as Int64);
                            x[row as usize] -= y[i as usize];
                        }
                        // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_sparse); (line 127) — dropped.
                    }

                    // HIPO_CLOCK_START(2); (line 131) — level-0 no-op, dropped.
                    // apply inverse swaps
                    {
                        let current_swaps = &self.swaps_[sn as usize][(nb * j) as usize..];
                        permuteWithSwaps(&mut x[x_start as usize..], current_swaps, jb, true);
                    }
                    // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 134) — dropped.
                }
            }
        }
    }

    /// HiGHS: HybridSolveHandler.cpp::HybridSolveHandler::backwardSolve (lines 141-264)
    fn backwardSolve(&self, x: &mut Vec<f64>) {
        // HIPO_CLOCK_CREATE; (line 147) — level-0 no-op, dropped.

        let nb = self.base.S_.blockSize();

        // go through the sn in reverse order
        let mut sn: Int = self.base.S_.sn() - 1;
        while sn >= 0 {
            // leading size of supernode
            let ldSn: Int = (self.base.S_.ptr(sn + 1) - self.base.S_.ptr(sn)) as Int;

            // number of columns in the supernode
            let sn_size = self.base.S_.snStart(sn + 1) - self.base.S_.snStart(sn);

            // first columns of the supernode
            let sn_start = self.base.S_.snStart(sn);

            // index to access S->rows for this supernode
            let start_row: Int64 = self.base.S_.ptr(sn);

            // number of blocks of columns
            let n_blocks = (sn_size - 1) / nb + 1;

            // index to access snColumns[sn]
            // initialised with the total number of entries of snColumns[sn]
            let mut SnCol_ind: Int64 =
                self.base.sn_columns_[sn as usize].len() as Int64 - extra_space_frontal as Int64;

            if sn_size < nb {
                // Fast solve
                // If supernode is small, avoid making BLAS calls
                let jb = sn_size;
                let x_start = sn_start;

                // HIPO_CLOCK_START(2); (line 179) — level-0 no-op, dropped.
                // apply swaps to portion of rhs that is affected
                {
                    let current_swaps = &self.swaps_[sn as usize];
                    permuteWithSwaps(&mut x[x_start as usize..], current_swaps, jb, false);
                }
                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 183) — dropped.

                // HIPO_CLOCK_START(2); (line 186) — level-0 no-op, dropped.
                let mut row: Int = ldSn - 1;
                while row >= jb {
                    let mut col: Int = jb - 1;
                    while col >= 0 {
                        x[x_start as usize + col as usize] -= self.base.sn_columns_[sn as usize]
                            [col as usize + (row * jb) as usize]
                            * x[self.base.S_.rows(start_row + row as Int64) as usize];
                        col -= 1;
                    }
                    row -= 1;
                }

                let mut row: Int = jb - 1;
                while row >= 0 {
                    let mut col: Int = row - 1;
                    while col >= 0 {
                        x[x_start as usize + col as usize] -= self.base.sn_columns_[sn as usize]
                            [col as usize + (row * jb) as usize]
                            * x[x_start as usize + row as usize];
                        col -= 1;
                    }
                    row -= 1;
                }
                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense); (line 200) — dropped.

                // HIPO_CLOCK_START(2); (line 203) — level-0 no-op, dropped.
                // apply inverse swaps
                {
                    let current_swaps = &self.swaps_[sn as usize];
                    permuteWithSwaps(&mut x[x_start as usize..], current_swaps, jb, true);
                }
                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 206) — dropped.
            } else {
                // go through blocks of columns for this supernode in reverse order
                let mut j: Int = n_blocks - 1;
                while j >= 0 {
                    // number of columns in the block
                    let jb = cpp_min(nb, sn_size - nb * j);

                    // number of entries in diagonal part
                    let diag_entries = jb * jb;

                    // index to access vector x
                    let x_start = sn_start + nb * j;

                    // HIPO_CLOCK_START(2); (line 222) — level-0 no-op, dropped.
                    // apply swaps to portion of rhs that is affected
                    {
                        let current_swaps = &self.swaps_[sn as usize][(nb * j) as usize..];
                        permuteWithSwaps(&mut x[x_start as usize..], current_swaps, jb, false);
                    }
                    // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 226) — dropped.

                    // temporary space for gemv
                    let gemv_space: Int = ldSn - nb * j - jb;
                    let mut y = vec![0.0f64; gemv_space as usize];
                    if gemv_space > 0 {
                        // HIPO_CLOCK_START(2); (line 233) — level-0 no-op, dropped.
                        // scatter entries into y
                        for i in 0..gemv_space {
                            let row = self
                                .base
                                .S_
                                .rows(start_row + (nb * j) as Int64 + jb as Int64 + i as Int64);
                            y[i as usize] = x[row as usize];
                        }
                        // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_sparse); (line 239) — dropped.

                        // HIPO_CLOCK_START(2); (line 241) — level-0 no-op, dropped.
                        SnCol_ind -= (jb * gemv_space) as Int64;
                        callAndTime_dgemv(
                            b'N',
                            jb,
                            gemv_space,
                            -1.0,
                            &self.base.sn_columns_[sn as usize][SnCol_ind as usize..],
                            jb,
                            &y,
                            1,
                            1.0,
                            &mut x[x_start as usize..],
                            1,
                            self.base.data_,
                        );
                        // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense); (line 246) — dropped.
                    }

                    // HIPO_CLOCK_START(2); (line 249) — level-0 no-op, dropped.
                    SnCol_ind -= diag_entries as Int64;
                    callAndTime_dtrsv(
                        b'U',
                        b'N',
                        b'U',
                        jb,
                        &self.base.sn_columns_[sn as usize][SnCol_ind as usize..],
                        jb,
                        &mut x[x_start as usize..],
                        1,
                        self.base.data_,
                    );
                    // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense); (line 253) — dropped.

                    // HIPO_CLOCK_START(2); (line 256) — level-0 no-op, dropped.
                    // apply inverse swaps
                    {
                        let current_swaps = &self.swaps_[sn as usize][(nb * j) as usize..];
                        permuteWithSwaps(&mut x[x_start as usize..], current_swaps, jb, true);
                    }
                    // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 259) — dropped.

                    j -= 1;
                }
            }

            sn -= 1;
        }
    }

    /// HiGHS: HybridSolveHandler.cpp::HybridSolveHandler::diagSolve (lines 266-350)
    fn diagSolve(&self, x: &mut Vec<f64>) {
        // HIPO_CLOCK_CREATE; (line 271) — level-0 no-op, dropped.

        let nb = self.base.S_.blockSize();

        for sn in 0..self.base.S_.sn() {
            // leading size of supernode
            let ldSn: Int = (self.base.S_.ptr(sn + 1) - self.base.S_.ptr(sn)) as Int;

            // number of columns in the supernode
            let sn_size = self.base.S_.snStart(sn + 1) - self.base.S_.snStart(sn);

            // first columns of the supernode
            let sn_start = self.base.S_.snStart(sn);

            // number of blocks of columns
            let n_blocks = (sn_size - 1) / nb + 1;

            // index to access diagonal part of block
            let mut diag_start: Int = 0;

            // go through blocks of columns for this supernode
            for j in 0..n_blocks {
                // number of columns in the block
                let jb = cpp_min(nb, sn_size - nb * j);

                // HIPO_CLOCK_START(2); (line 297) — level-0 no-op, dropped.
                // apply swaps to portion of rhs that is affected
                {
                    let current_swaps = &self.swaps_[sn as usize][(nb * j) as usize..];
                    permuteWithSwaps(
                        &mut x[(sn_start + nb * j) as usize..],
                        current_swaps,
                        jb,
                        false,
                    );
                }
                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 301) — dropped.

                // HIPO_CLOCK_START(2); (line 304) — level-0 no-op, dropped.

                let current_2x2 = &self.pivot_2x2_[sn as usize][(nb * j) as usize..];
                // C++ `Int step = 1;` — the initial value is never read (every branch sets it).
                let mut step: Int;

                // go through columns of block
                let mut col: Int = 0;
                while col < jb {
                    if current_2x2[col as usize] == 0.0 {
                        // 1x1 pivots
                        step = 1;
                        let inv_d = self.base.sn_columns_[sn as usize]
                            [diag_start as usize + col as usize + (jb * col) as usize];
                        x[(sn_start + nb * j + col) as usize] *= inv_d;
                    } else {
                        // 2x2 pivots
                        step = 2;

                        // inverse of 2x2 pivot
                        let i_d1 = self.base.sn_columns_[sn as usize]
                            [diag_start as usize + col as usize + (jb * col) as usize];
                        let i_d2 = self.base.sn_columns_[sn as usize]
                            [diag_start as usize + (col + 1) as usize + (jb * (col + 1)) as usize];
                        let i_off = current_2x2[col as usize];

                        let x1 = x[(sn_start + nb * j + col) as usize];
                        let x2 = x[(sn_start + nb * j + col + 1) as usize];

                        x[(sn_start + nb * j + col) as usize] = i_d1 * x1 + i_off * x2;
                        x[(sn_start + nb * j + col + 1) as usize] = i_d2 * x2 + i_off * x1;
                    }
                    col += step;
                }

                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense); (line 334) — dropped.

                // HIPO_CLOCK_START(2); (line 337) — level-0 no-op, dropped.
                // apply inverse swaps
                {
                    let current_swaps = &self.swaps_[sn as usize][(nb * j) as usize..];
                    permuteWithSwaps(
                        &mut x[(sn_start + nb * j) as usize..],
                        current_swaps,
                        jb,
                        true,
                    );
                }
                // HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap); (line 340) — dropped.

                // move diag_start forward by number of diagonal entries in block
                diag_start += jb * jb;

                // move diag_start forward by number of sub-diagonal entries in block
                diag_start += (ldSn - nb * j - jb) * jb;
            }
        }
    }
}
