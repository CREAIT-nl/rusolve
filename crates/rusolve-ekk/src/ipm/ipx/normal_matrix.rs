//! HiGHS source: ipm/ipx/normal_matrix.h, normal_matrix.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::linear_operator::LinearOperator;
use super::model::Model;
use super::timer::Timer;
use super::utils::dot;
use super::*;
use std::rc::Rc;

// HiGHS normal_matrix.cc: `#define MATVECMETHOD 1`. Only the one-pass variant
// is translated; methods 2 and 3 (which use `work_`) are dead code in HiGHS
// and are not reproduced, so `work_` is never sized.

/// HiGHS: ipx/normal_matrix.h::NormalMatrix
///
/// PORT NOTE: C++ stores the pointer `const double* W_` passed to
/// `Prepare` (it points into `KKTSolverDiag::W_`, which is not changed
/// until the next `Prepare`). Rust copies the `n+m` values into `W_`;
/// `W_is_null_` records `Prepare(nullptr)`.
pub struct NormalMatrix {
    model_: Rc<Model>,
    W_: Vec<f64>,
    W_is_null_: bool,
    prepared_: bool,
    work_: Vector,
    time_: f64,
    /// Every `AI.rowidx()` entry lies in `0..rows()`, checked once in `new`
    /// (the `Rc<Model>` is immutable for this operator's lifetime); it makes
    /// the unchecked gather/scatter in `_apply` sound.
    rowidx_in_range_: bool,
}

impl NormalMatrix {
    /// HiGHS: ipx/normal_matrix.cc::NormalMatrix::NormalMatrix (lines 22-28)
    pub fn new(model: Rc<Model>) -> Self {
        let m = model.rows();
        let rowidx_in_range_ = model.AI().rowidx().iter().all(|&i| 0 <= i && i < m);
        Self {
            model_: model,
            // C++ `const double* W_{nullptr}`
            W_: Vec::new(),
            W_is_null_: true,
            prepared_: false,
            // MATVECMETHOD 1: no workspace
            work_: Vector::default(),
            time_: 0.0,
            rowidx_in_range_,
        }
    }
    /// HiGHS: ipx/normal_matrix.cc::NormalMatrix::Prepare (lines 30-33)
    pub fn prepare(&mut self, W: Option<&[f64]>) {
        match W {
            Some(W) => {
                self.W_ = W.to_vec();
                self.W_is_null_ = false;
            }
            None => {
                self.W_ = Vec::new();
                self.W_is_null_ = true;
            }
        }
        self.prepared_ = true;
    }
    /// HiGHS: ipx/normal_matrix.cc::NormalMatrix::time (lines 35-37)
    pub fn time(&self) -> f64 {
        self.time_
    }
    /// HiGHS: ipx/normal_matrix.cc::NormalMatrix::reset_time (lines 39-41)
    pub fn reset_time(&mut self) {
        self.time_ = 0.0;
    }
}

impl LinearOperator for NormalMatrix {
    /// HiGHS: ipx/normal_matrix.cc::NormalMatrix::_Apply (lines 43-124)
    fn _apply(&mut self, rhs: &Vector, lhs: &mut Vector, rhs_dot_lhs: Option<&mut f64>) {
        let m = self.model_.rows();
        let n = self.model_.cols();
        let AI = self.model_.AI();
        let Ap = AI.colptr();
        let Ai = AI.rowidx();
        let Ax = AI.values();
        let timer = Timer::default();

        debug_assert!(self.prepared_);
        debug_assert_eq!(lhs.size(), m);
        debug_assert_eq!(rhs.size(), m);

        // PORT NOTE (perf): column slices `Ai[begin..end]`/`Ax[begin..end]` are
        // taken once per column instead of indexing `Ai[p]`, `Ax[p]` per nonzero,
        // and the gather/scatter through `Ai` is unchecked (see `gather_dot` /
        // `scatter_axpy`); the accumulation order is that of normal_matrix.cc,
        // so `lhs` is bit-identical.
        assert!(self.rowidx_in_range_);
        let rhs = &rhs[..m as usize];
        let lhs = &mut lhs[..m as usize];
        if !self.W_is_null_ {
            let W = &self.W_;
            for ((l, &r), &w) in lhs
                .iter_mut()
                .zip(rhs)
                .zip(&W[n as usize..(n + m) as usize])
            {
                *l = r * w;
            }
            for j in 0..n {
                let begin = Ap[j as usize] as usize;
                let end = Ap[(j + 1) as usize] as usize;
                let (ai, ax) = (&Ai[begin..end], &Ax[begin..end]);
                let mut d = gather_dot(rhs, ai, ax);
                d *= W[j as usize];
                scatter_axpy(lhs, ai, ax, d);
            }
        } else {
            lhs.fill(0.0);
            for j in 0..n {
                let begin = Ap[j as usize] as usize;
                let end = Ap[(j + 1) as usize] as usize;
                let (ai, ax) = (&Ai[begin..end], &Ax[begin..end]);
                let d = gather_dot(rhs, ai, ax);
                scatter_axpy(lhs, ai, ax, d);
            }
        }
        if let Some(out) = rhs_dot_lhs {
            *out = dot(rhs, lhs);
        }
        self.time_ += timer.elapsed();
    }
}

/// `sum_k x[ai[k]] * ax[k]`, accumulated left to right from `0.0`.
///
/// Callers guarantee every `ai[k]` lies in `0..x.len()`: `ai` is a column
/// slice of the model's `AI.rowidx()`, range-checked once in
/// `NormalMatrix::new` (`rowidx_in_range_`), and `x` has length `rows()`
/// (`SplittedNormalMatrix::_apply`: a column slice of `N_`, range-checked once
/// per `prepare` into `N_rowidx_in_range_`).
#[inline(always)]
pub(super) fn gather_dot(x: &[f64], ai: &[Int], ax: &[f64]) -> f64 {
    let mut d = 0.0;
    for (&i, &a) in ai.iter().zip(ax) {
        debug_assert!((i as usize) < x.len());
        // SAFETY: `i` is in `0..x.len()` (see the function comment).
        d += unsafe { *x.get_unchecked(i as usize) } * a;
    }
    d
}

/// `y[ai[k]] += d * ax[k]` for each k in order; same index contract as
/// [`gather_dot`].
#[inline(always)]
pub(super) fn scatter_axpy(y: &mut [f64], ai: &[Int], ax: &[f64], d: f64) {
    for (&i, &a) in ai.iter().zip(ax) {
        debug_assert!((i as usize) < y.len());
        // SAFETY: `i` is in `0..y.len()` (see `gather_dot`).
        unsafe { *y.get_unchecked_mut(i as usize) += d * a };
    }
}
