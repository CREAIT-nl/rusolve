//! HiGHS source: ipm/ipx/splitted_normal_matrix.h, splitted_normal_matrix.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::basis::{BASIC, BASIC_FREE, Basis, NONBASIC};
use super::linear_operator::LinearOperator;
use super::model::Model;
use super::normal_matrix::{gather_dot, scatter_axpy};
use super::sparse_matrix::{
    SparseMatrix, backward_solve, copy_columns, forward_solve, permute_rows, scale_column,
};
use super::timer::Timer;
use super::utils::{dot, inverse_perm};
use super::*;
use crate::util::highs_cdouble::HighsCDouble;
use std::rc::Rc;

/// HiGHS 1.16 (#3513): ipm/ipx/splitted_normal_matrix.cc::PreciseSolve (lines 13-32)
///
/// Accumulate each triangular solve in double-double precision.
fn precise_solve(
    a: &SparseMatrix,
    x: &mut [HighsCDouble],
    transpose: bool,
    lower: bool,
    unit: bool,
) {
    let forward = transpose != lower;
    for step in 0..a.cols() {
        let j = if forward { step } else { a.cols() - 1 - step };
        let begin = a.begin(j) + (lower && !unit) as Int;
        let end = a.end(j) - (!lower && !unit) as Int;
        let diagonal = if lower { begin - 1 } else { end };
        if transpose {
            let mut sum = HighsCDouble::from(0.0);
            for k in begin..end {
                sum += x[a.index(k) as usize] * a.value(k);
            }
            x[j as usize] -= sum;
        }
        if !unit {
            x[j as usize] /= a.value(diagonal);
        }
        if !transpose {
            for k in begin..end {
                let t = x[j as usize] * a.value(k);
                x[a.index(k) as usize] -= t;
            }
        }
    }
}

/// HiGHS: ipx/splitted_normal_matrix.h::SplittedNormalMatrix
pub struct SplittedNormalMatrix {
    model_: Rc<Model>,
    L_: SparseMatrix,
    U_: SparseMatrix,
    N_: SparseMatrix,
    free_positions_: Vec<Int>,
    colperm_: Vec<Int>,
    rowperm_inv_: Vec<Int>,
    work_: Vector,
    // HiGHS 1.16 #3513 (h16_ipx_scaled): retain precision after roundoff failure; false
    // unless useExtendedPrecision() is called on the flag-on path.
    precise_: bool,
    prepared_: bool,
    /// Port-only: every `N_` row index lies in `0..rows()` (checked in `prepare`).
    N_rowidx_in_range_: bool,
    time_B_: f64,
    time_Bt_: f64,
    time_NNt_: f64,
}

impl SplittedNormalMatrix {
    /// HiGHS: ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::SplittedNormalMatrix (lines 9-14)
    pub fn new(model: Rc<Model>) -> Self {
        let m = model.rows();
        Self {
            model_: model,
            L_: SparseMatrix::new(),
            U_: SparseMatrix::new(),
            N_: SparseMatrix::new(),
            free_positions_: Vec::new(),
            // C++ colperm_.resize(m), rowperm_inv_.resize(m): zeros
            colperm_: vec![0; m as usize],
            rowperm_inv_: vec![0; m as usize],
            work_: Vector::new(m),
            precise_: false,
            prepared_: false,
            N_rowidx_in_range_: false,
            time_B_: 0.0,
            time_Bt_: 0.0,
            time_NNt_: 0.0,
        }
    }
    /// HiGHS: ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::Prepare (lines 16-64)
    pub fn prepare(&mut self, basis: &Basis, colscale: &[f64]) {
        let model = self.model_.clone();
        let m = model.rows();
        let n = model.cols();
        let AI = model.AI();
        // C++ `assert(colscale)` (non-null pointer): a slice is never null.
        self.prepared_ = false;
        self.N_.clear(); // deallocate old memory

        basis.get_lu_factors(
            &mut self.L_,
            &mut self.U_,
            &mut self.rowperm_inv_,
            &mut self.colperm_,
        );
        self.rowperm_inv_ = inverse_perm(&self.rowperm_inv_);

        // Scale columns of U.
        for k in 0..m {
            let p = self.colperm_[k as usize];
            let j = basis.at(p);
            // Nothing to do for BASIC_FREE variables.
            if basis.status_of(j) == BASIC {
                let d = colscale[j as usize];
                debug_assert!(d.is_finite() && d > 0.0);
                scale_column(&mut self.U_, k, d);
            }
        }

        // Build N with permuted row indices.
        let mut nonbasic_vars: Vec<Int> = Vec::new();
        for j in 0..n + m {
            if basis.status_of(j) == NONBASIC {
                nonbasic_vars.push(j);
            }
        }
        self.N_ = copy_columns(AI, &nonbasic_vars);
        permute_rows(&mut self.N_, &self.rowperm_inv_);

        // Scale columns of N.
        for (k, &j) in nonbasic_vars.iter().enumerate() {
            let d = colscale[j as usize];
            debug_assert!(d.is_finite());
            scale_column(&mut self.N_, k as Int, d);
        }

        // Build list of free variables.
        self.free_positions_.clear();
        for k in 0..m {
            let p = self.colperm_[k as usize];
            let j = basis.at(p);
            if basis.status_of(j) == BASIC_FREE {
                self.free_positions_.push(k);
            }
        }
        // PORT NOTE (perf): one range check of `N_`'s row indices per
        // `prepare`, so `_apply` can gather/scatter through them unchecked.
        self.N_rowidx_in_range_ = self.N_.rowidx().iter().all(|&i| 0 <= i && i < m);
        self.prepared_ = true;
    }
    /// HiGHS: ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::colperm (lines 66-68)
    pub fn colperm(&self) -> &[Int] {
        &self.colperm_
    }
    /// HiGHS: ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::time_B (lines 70-72)
    pub fn time_B(&self) -> f64 {
        self.time_B_
    }
    /// HiGHS: ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::time_Bt (lines 74-76)
    pub fn time_Bt(&self) -> f64 {
        self.time_Bt_
    }
    /// HiGHS: ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::time_NNt (lines 78-80)
    pub fn time_NNt(&self) -> f64 {
        self.time_NNt_
    }
    /// HiGHS: ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::reset_time (lines 82-86)
    pub fn reset_time(&mut self) {
        self.time_B_ = 0.0;
        self.time_Bt_ = 0.0;
        self.time_NNt_ = 0.0;
    }
    /// HiGHS 1.16 (#3513): ipm/ipx/splitted_normal_matrix.h::SplittedNormalMatrix::useExtendedPrecision (lines 34-38)
    ///
    /// Returns true only when upgrading, so the caller retries at most once.
    pub fn use_extended_precision(&mut self) -> bool {
        if self.precise_ {
            return false;
        }
        self.precise_ = true;
        true
    }
    /// HiGHS 1.16 (#3513): ipm/ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::ApplyPrecise (lines 144-178)
    fn apply_precise(&mut self, rhs: &Vector, lhs: &mut Vector, rhs_dot_lhs: Option<&mut f64>) {
        debug_assert!(self.prepared_);
        let mut timer = Timer::default();
        let size = rhs.size() as usize;
        let mut w: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); size];
        let mut z: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); size];
        for i in 0..size {
            w[i] = HighsCDouble::from(rhs[i]);
        }
        timer.reset(false);
        precise_solve(&self.U_, &mut w, true, false, false);
        precise_solve(&self.L_, &mut w, true, true, true);
        self.time_Bt_ += timer.elapsed();
        timer.reset(false);
        let mut energy = HighsCDouble::from(0.0);
        for j in 0..self.N_.cols() {
            let ju = j as usize;
            let mut d = HighsCDouble::from(0.0);
            for k in self.N_.colptr()[ju]..self.N_.colptr()[ju + 1] {
                let ku = k as usize;
                d += w[self.N_.rowidx()[ku] as usize] * self.N_.values()[ku];
            }
            // Accumulation can leave cancelling high/low parts; normalize
            // before multiplying two HighsCDouble values.
            d.renormalize();
            energy += d * d;
            for k in self.N_.colptr()[ju]..self.N_.colptr()[ju + 1] {
                let ku = k as usize;
                z[self.N_.rowidx()[ku] as usize] += d * self.N_.values()[ku];
            }
        }
        self.time_NNt_ += timer.elapsed();
        timer.reset(false);
        precise_solve(&self.L_, &mut z, false, true, true);
        precise_solve(&self.U_, &mut z, false, false, false);
        self.time_B_ += timer.elapsed();
        for i in 0..size {
            lhs[i] = f64::from(z[i] + rhs[i]);
            energy += HighsCDouble::from(rhs[i]) * rhs[i];
        }
        for &i in self.free_positions_.iter() {
            lhs[i as usize] = 0.0;
        }
        if let Some(out) = rhs_dot_lhs {
            *out = f64::from(energy);
        }
    }
}

impl LinearOperator for SplittedNormalMatrix {
    /// HiGHS: ipx/splitted_normal_matrix.cc::SplittedNormalMatrix::_Apply (lines 88-115)
    fn _apply(&mut self, rhs: &Vector, lhs: &mut Vector, rhs_dot_lhs: Option<&mut f64>) {
        // HiGHS 1.16 #3513 (h16_ipx_scaled): ipm/ipx/splitted_normal_matrix.cc line 115 (1.16.0)
        if self.precise_ {
            self.apply_precise(rhs, lhs, rhs_dot_lhs);
            return;
        }
        debug_assert!(self.prepared_);
        let mut timer = Timer::default();

        // Compute work = inverse(B') * rhs.
        self.work_.assign(rhs);
        timer.reset(false);
        backward_solve(&self.L_, &self.U_, &mut self.work_);
        self.time_Bt_ += timer.elapsed();

        // Compute lhs = N*N' * work.
        lhs.fill(0.0);
        timer.reset(false);
        // PORT NOTE (perf): `AddNormalProduct(N_, nullptr, work_, lhs)`
        // with the normal_matrix.rs unchecked gather/scatter (index range checked in `prepare`);
        // same column order and accumulation order, so `lhs` is bit-identical.
        {
            assert!(self.N_rowidx_in_range_);
            let m = self.N_.rows() as usize;
            let (Np, Ni, Nx) = (self.N_.colptr(), self.N_.rowidx(), self.N_.values());
            let work = &self.work_[..m];
            let lhs = &mut lhs[..m];
            for j in 0..self.N_.cols() as usize {
                let (begin, end) = (Np[j] as usize, Np[j + 1] as usize);
                let (ai, ax) = (&Ni[begin..end], &Nx[begin..end]);
                let d = gather_dot(work, ai, ax);
                scatter_axpy(lhs, ai, ax, d);
            }
        }
        self.time_NNt_ += timer.elapsed();

        // Compute lhs := inverse(B) * lhs.
        timer.reset(false);
        forward_solve(&self.L_, &self.U_, lhs);
        self.time_B_ += timer.elapsed();

        // lhs += rhs
        for i in 0..lhs.size() as usize {
            lhs[i] = lhs[i] + rhs[i];
        }
        for &i in self.free_positions_.iter() {
            lhs[i as usize] = 0.0;
        }
        if let Some(out) = rhs_dot_lhs {
            *out = dot(rhs, lhs);
        }
    }
}
