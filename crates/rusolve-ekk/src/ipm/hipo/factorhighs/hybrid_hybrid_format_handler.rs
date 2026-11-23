//! HiGHS source: ipm/hipo/factorhighs/HybridHybridFormatHandler.h, ipm/hipo/factorhighs/HybridHybridFormatHandler.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The hybrid (blocked lower-packed) frontal and clique format, the only
//! `FormatHandler` HiPO instantiates.

use super::super::{Int, Int64};
use super::Regul;
use super::data_collector::DataCollector;
use super::dense_fact_hybrid::{denseFactFH, denseFactFP2FH};
use super::format_handler::{FormatHandler, FormatHandlerOps, extra_space_frontal};
use super::kMinConsecutiveSums;

/// rusolve: `y[k] += x[k]` over a unit-stride run, what `daxpy(n, 1.0, x, 1, y, 1)`
/// computes (`y + 1.0 * x` is `y + x` exactly) without the strided reference loop, so the compiler
/// vectorises it. Bit-identical to the BLAS call.
#[inline]
fn add_run(y: &mut [f64], x: &[f64], n: usize) {
    for (yv, xv) in y[..n].iter_mut().zip(&x[..n]) {
        *yv += *xv;
    }
}

use super::symbolic::Symbolic;
use crate::ipm::hipo::auxiliary::getDiagStart;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// HiGHS: HybridHybridFormatHandler.h `class HybridHybridFormatHandler`
pub struct HybridHybridFormatHandler<'a> {
    pub base: FormatHandler<'a>,
    pub diag_start_: Vec<Int64>,
    pub data_: &'a DataCollector,
}

impl<'a> HybridHybridFormatHandler<'a> {
    /// HiGHS: HybridHybridFormatHandler.cpp::HybridHybridFormatHandler::HybridHybridFormatHandler (lines 13-23)
    pub fn new(
        S: &'a Symbolic,
        sn: Int,
        regul: &Regul,
        data: &'a DataCollector,
        frontal: &'a mut Vec<f64>,
        clique_ptr: Option<Int64>,
    ) -> Self {
        let base = FormatHandler::new(S, sn, regul, frontal, clique_ptr);
        let mut h = HybridHybridFormatHandler {
            base,
            diag_start_: Vec::new(),
            data_: data,
        };

        // initialise frontal and clique
        h.initFrontal();

        // if CliqueStack is used, clique_ptr already points to a valid region of
        // memory for the clique. Otherwise, allocate it locally.
        if h.base.clique_ptr_.is_none() {
            h.initClique();
        }

        h
    }
}

impl<'a> FormatHandlerOps<'a> for HybridHybridFormatHandler<'a> {
    fn base(&self) -> &FormatHandler<'a> {
        &self.base
    }
    fn base_mut(&mut self) -> &mut FormatHandler<'a> {
        &mut self.base
    }

    /// HiGHS: HybridHybridFormatHandler.cpp::HybridHybridFormatHandler::initFrontal (lines 25-39)
    fn initFrontal(&mut self) {
        let n_blocks = (self.base.sn_size_ - 1) / self.base.nb_ + 1;
        self.diag_start_.resize(n_blocks as usize, 0 as Int64);
        let frontal_size = getDiagStart(
            self.base.ldf_,
            self.base.sn_size_,
            self.base.nb_,
            n_blocks,
            &mut self.diag_start_,
            false,
        ) + extra_space_frontal as Int64;
        self.base.frontal_.resize(frontal_size as usize, 0.0);
        // C++ memset(frontal_.data(), 0, frontal_size * sizeof(double)):
        // explicit zero fill of the whole frontal_.
        self.base.frontal_[..frontal_size as usize].fill(0.0);

        // NB: extra_space_frontal is not strictly needed. However, it removes some
        // weird problem on windows in debug. Who knows what's happening...

        // frontal_ is actually allocated just the first time, then the memory is
        // reused from the previous factorisations and just initialised.
    }

    /// HiGHS: HybridHybridFormatHandler.cpp::HybridHybridFormatHandler::initClique (lines 41-48)
    fn initClique(&mut self) {
        self.base
            .clique_
            .resize(self.base.S_.cliqueSize(self.base.sn_) as usize, 0.0);

        // If the clique size is zero, do not access the underlying pointer. This
        // causes strange issues on windows. It's not a problem if clique_ptr_ remains
        // null, because it will never be used in that case.
        // C++: if (!clique_.empty()) clique_ptr_ = clique_.data();
        // In Rust the pointer is not representable; clique_ptr_ stays None.
    }

    /// HiGHS: HybridHybridFormatHandler.cpp::HybridHybridFormatHandler::assembleFrontal (lines 50-56)
    fn assembleFrontal(&mut self, i: Int, j: Int, val: f64) {
        let block = j / self.base.nb_;
        let ldb = self.base.ldf_ - block * self.base.nb_;
        let ii = i - block * self.base.nb_;
        let jj = j - block * self.base.nb_;
        let idx = (self.diag_start_[block as usize] + ii as Int64 + (ldb * jj) as Int64) as usize;
        self.base.frontal_[idx] = val;
    }

    /// HiGHS: HybridHybridFormatHandler.cpp::HybridHybridFormatHandler::assembleFrontalMultiple (lines 58-82)
    fn assembleFrontalMultiple(
        &mut self,
        num: &mut Int,
        child: &[f64],
        nc: Int,
        child_sn: Int,
        row: Int,
        col: Int,
        i: Int,
        j: Int,
    ) {
        let jblock = col / self.base.nb_;
        let jb = cpp_min(self.base.nb_, nc - self.base.nb_ * jblock);
        let row_ = row - jblock * self.base.nb_;
        let col_ = col - jblock * self.base.nb_;
        let start_block = self.base.S_.cliqueBlockStart(child_sn, jblock);

        let block = j / self.base.nb_;
        let ldb = self.base.ldf_ - block * self.base.nb_;
        let ii = i - block * self.base.nb_;
        let jj = j - block * self.base.nb_;

        if *num > kMinConsecutiveSums {
            let child_idx = (start_block + col_ as Int64 + (jb * row_) as Int64) as usize;
            let frontal_idx =
                (self.diag_start_[block as usize] + ii as Int64 + (ldb * jj) as Int64) as usize;
            // strided child (row-major block of width jb), unit-stride frontal; as daxpy(.., 1.0, ..)
            let y = &mut self.base.frontal_[frontal_idx..frontal_idx + *num as usize];
            for (k, yv) in y.iter_mut().enumerate() {
                *yv += child[child_idx + k * jb as usize];
            }
        } else {
            let frontal_idx =
                (self.diag_start_[block as usize] + ii as Int64 + (ldb * jj) as Int64) as usize;
            let child_idx = (start_block + col_ as Int64 + (jb * row_) as Int64) as usize;
            self.base.frontal_[frontal_idx] += child[child_idx];
            *num = 1;
        }
    }

    /// HiGHS: HybridHybridFormatHandler.cpp::HybridHybridFormatHandler::denseFactorise (lines 84-103)
    fn denseFactorise(&mut self, stack: &mut [f64], reg_thresh: f64) -> Int {
        // either clique is valid, or clique is not needed
        // C++ clique_ptr_ is non-null for a stack clique (Some) or a non-empty
        // local clique_ (initClique sets it to clique_.data()).
        debug_assert!(
            self.base.clique_ptr_.is_some()
                || !self.base.clique_.is_empty()
                || self.base.ldf_ == self.base.sn_size_
        );

        let mut status = denseFactFP2FH(
            self.base.frontal_.as_mut_slice(),
            self.base.ldf_,
            self.base.sn_size_,
            self.base.nb_,
            self.data_,
        );
        if status != 0 {
            return status;
        }

        // find the position within pivot_sign corresponding to this supernode
        let sn_start = self.base.S_.snStart(self.base.sn_);
        let pivot_sign = &self.base.S_.pivotSign()[sn_start as usize..];

        // clique_ptr_ = Some(o): the clique lives in stack[o..]; None (C++ nullptr):
        // the clique is the local clique_.
        let clique: &mut [f64] = match self.base.clique_ptr_ {
            Some(o) => &mut stack[o as usize..],
            None => &mut self.base.clique_[..],
        };

        status = denseFactFH(
            b'H',
            self.base.ldf_,
            self.base.sn_size_,
            self.base.S_.blockSize(),
            self.base.frontal_.as_mut_slice(),
            clique,
            pivot_sign,
            reg_thresh,
            &self.base.regul_,
            self.base.local_reg_.as_mut_slice(),
            self.base.swaps_.as_mut_slice(),
            self.base.pivot_2x2_.as_mut_slice(),
            self.base.S_.parNode(),
            self.data_,
        );
        status
    }

    /// HiGHS: HybridHybridFormatHandler.cpp::HybridHybridFormatHandler::assembleClique (lines 105-186)
    fn assembleClique(&mut self, stack: &mut [f64], child_off: Int64, nc: Int, child_sn: Int) {
        // assemble the child clique into the current clique by blocks of columns.
        // within a block, assemble by rows.

        let nb = self.base.nb_;
        let ldc = self.base.ldc_;
        let sn_size = self.base.sn_size_;
        let sn = self.base.sn_;
        let S = self.base.S_;

        let n_blocks = (nc - 1) / nb + 1;

        let mut row_start: Int = 0;

        // go through the blocks of columns of the child sn
        for b in 0..n_blocks {
            let b_start = S.cliqueBlockStart(child_sn, b);

            let col_start = row_start;
            let col_end = cpp_min(col_start + nb, nc);

            // go through the rows within this block
            for row in row_start..nc {
                let i = S.relindClique(child_sn, row) - sn_size;

                // already assembled into frontal
                if i < 0 {
                    continue;
                }

                // go through the columns of the block
                let mut col = col_start;
                while col < col_end {
                    let mut j = S.relindClique(child_sn, col);
                    if j < sn_size {
                        col += 1;
                        continue;
                    }
                    j -= sn_size;

                    // information and sizes of child sn
                    let jblock_c = b;
                    let jb_c = cpp_min(nb, nc - nb * jblock_c);
                    let row_ = row - jblock_c * nb;
                    let col_ = col - jblock_c * nb;
                    let start_block_c = b_start;

                    // sum consecutive entries in a row.
                    // consecutive need to be reduced, to account for edge of the block
                    let zeros_stored_row = cpp_max(0 as Int, jb_c - (row - row_start) - 1);
                    let mut consecutive = S.consecutiveSums(child_sn, col);
                    let left_in_child = col_end - col - zeros_stored_row;
                    consecutive = cpp_min(consecutive, left_in_child);

                    // consecutive need to account also for edge of block in parent
                    let block_in_parent = j / nb;
                    let col_end_parent = cpp_min((block_in_parent + 1) * nb, ldc);
                    let left_in_parent = col_end_parent - j;
                    consecutive = cpp_min(consecutive, left_in_parent);

                    // needed to deal with zeros stored in upper right part of block
                    if consecutive == 0 {
                        break;
                    }

                    // information and sizes of current sn
                    let jblock = block_in_parent;
                    let jb = cpp_min(nb, ldc - nb * jblock);
                    let i_ = i - jblock * nb;
                    let j_ = j - jblock * nb;
                    let start_block = S.cliqueBlockStart(sn, jblock);

                    if consecutive > kMinConsecutiveSums {
                        // child index in stack: child_off + start_block_c + col_ + jb_c * row_
                        let c_idx =
                            (child_off + start_block_c + col_ as Int64 + (jb_c * row_) as Int64)
                                as usize;

                        match self.base.clique_ptr_ {
                            Some(o) => {
                                let p_idx =
                                    (o + start_block + j_ as Int64 + (jb * i_) as Int64) as usize;
                                // PORT NOTE: child and parent cliques are disjoint regions
                                // of the stack; split it to borrow both.
                                let (child, parent): (&[f64], &mut [f64]) = if c_idx < p_idx {
                                    let (lo, hi) = stack.split_at_mut(p_idx);
                                    (&lo[c_idx..], hi)
                                } else {
                                    let (lo, hi) = stack.split_at_mut(c_idx);
                                    (hi, &mut lo[p_idx..])
                                };
                                add_run(parent, child, consecutive as usize);
                            }
                            None => {
                                // clique_ptr_ is the local clique_.
                                let p_idx =
                                    (start_block + j_ as Int64 + (jb * i_) as Int64) as usize;
                                add_run(
                                    &mut self.base.clique_[p_idx..],
                                    &stack[c_idx..],
                                    consecutive as usize,
                                );
                            }
                        }

                        col += consecutive;
                    } else {
                        // scalar += case
                        let c_idx =
                            (child_off + start_block_c + col_ as Int64 + (jb_c * row_) as Int64)
                                as usize;

                        match self.base.clique_ptr_ {
                            Some(o) => {
                                let p_idx =
                                    (o + start_block + j_ as Int64 + (jb * i_) as Int64) as usize;
                                stack[p_idx] += stack[c_idx];
                            }
                            None => {
                                let p_idx =
                                    (start_block + j_ as Int64 + (jb * i_) as Int64) as usize;
                                self.base.clique_[p_idx] += stack[c_idx];
                            }
                        }

                        col += 1;
                    }
                }
            }

            row_start += nb;
        }
    }

    /// HiGHS: HybridHybridFormatHandler.cpp::HybridHybridFormatHandler::extremeEntries (lines 188-236)
    /// Body is `#ifdef HIPO_COLLECT_EXPENSIVE_DATA`: empty in the HiGHS build.
    fn extremeEntries(&mut self) {
        // HIPO_COLLECT_EXPENSIVE_DATA is not defined in the HiGHS build.
        // The entire body (lines 189-235) is dropped.
    }
}
