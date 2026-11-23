//! HiGHS source: ipm/ipx/diagonal_precond.h, diagonal_precond.cc
//!
//! Copyright (c) 2018-2019 ERGO-Code. See license.txt for license.
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::linear_operator::LinearOperator;
use super::model::Model;
use super::timer::Timer;
use super::*;
use std::rc::Rc;

/// HiGHS: ipx/diagonal_precond.h::DiagonalPrecond
pub struct DiagonalPrecond {
    model_: Rc<Model>,
    factorized_: bool,
    diagonal_: Vector,
    time_: f64,
}

impl DiagonalPrecond {
    /// HiGHS: ipx/diagonal_precond.cc::DiagonalPrecond::DiagonalPrecond (lines 11-14)
    pub fn new(model: Rc<Model>) -> Self {
        // C++: diagonal_.resize(m)
        let m = model.rows();
        Self {
            model_: model,
            factorized_: false,
            diagonal_: Vector::new(m),
            time_: 0.0,
        }
    }
    /// HiGHS: ipx/diagonal_precond.cc::DiagonalPrecond::Factorize (lines 16-41)
    ///
    /// `W` may be `nullptr`.
    pub fn factorize(&mut self, W: Option<&[f64]>, info: &mut Info) {
        let m = self.model_.rows();
        let n = self.model_.cols();
        let AI = self.model_.AI();

        self.factorized_ = false;

        // Build diagonal of normal matrix.
        match W {
            Some(W) => {
                for i in 0..m {
                    self.diagonal_[i as usize] = W[(n + i) as usize];
                }
                for j in 0..n {
                    let w = W[j as usize];
                    for p in AI.begin(j)..AI.end(j) {
                        let i = AI.index(p) as usize;
                        self.diagonal_[i] += AI.value(p) * w * AI.value(p);
                    }
                }
            }
            None => {
                // rightmost m columns have weight zero
                self.diagonal_.fill(0.0);
                for j in 0..n {
                    for p in AI.begin(j)..AI.end(j) {
                        let i = AI.index(p) as usize;
                        self.diagonal_[i] += AI.value(p) * AI.value(p);
                    }
                }
            }
        }

        self.factorized_ = true;
        // C++ does not touch `info` here.
        let _ = info;
    }
    /// HiGHS: ipx/diagonal_precond.cc::DiagonalPrecond::time (lines 43-45)
    pub fn time(&self) -> f64 {
        self.time_
    }
    /// HiGHS: ipx/diagonal_precond.cc::DiagonalPrecond::reset_time (lines 47-49)
    pub fn reset_time(&mut self) {
        self.time_ = 0.0;
    }
}

impl LinearOperator for DiagonalPrecond {
    /// HiGHS: ipx/diagonal_precond.cc::DiagonalPrecond::_Apply (lines 51-68)
    fn _apply(&mut self, rhs: &Vector, lhs: &mut Vector, rhs_dot_lhs: Option<&mut f64>) {
        let m = self.model_.rows();
        let mut rldot = 0.0;
        let timer = Timer::default();

        debug_assert!(self.factorized_);
        debug_assert_eq!(lhs.size(), m);
        debug_assert_eq!(rhs.size(), m);

        // PORT NOTE (perf): zipped slices instead of three bounds-checked
        // indexings per element; same division and accumulation order.
        let mu = m as usize;
        for ((l, &r), &d) in lhs[..mu]
            .iter_mut()
            .zip(&rhs[..mu])
            .zip(&self.diagonal_[..mu])
        {
            *l = r / d;
            rldot += *l * r;
        }
        if let Some(out) = rhs_dot_lhs {
            *out = rldot;
        }
        self.time_ += timer.elapsed();
    }
}
