//! HiGHS source: ipm/ipx/kkt_solver_diag.h, kkt_solver_diag.cc
//!
//! Copyright (c) 2018-2019 ERGO-Code. See license.txt for license.
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::conjugate_residuals::ConjugateResiduals;
use super::control::Control;
use super::diagonal_precond::DiagonalPrecond;
use super::iterate::Iterate;
use super::kkt_solver::KKTSolver;
use super::model::Model;
use super::normal_matrix::NormalMatrix;
use super::sparse_matrix::{dot_column, scatter_column};
use super::*;
use std::rc::Rc;

/// HiGHS: ipx/kkt_solver_diag.h::KKTSolverDiag
pub struct KKTSolverDiag {
    control_: Rc<Control>,
    model_: Rc<Model>,
    normal_matrix_: NormalMatrix,
    precond_: DiagonalPrecond,
    W_: Vector,
    resscale_: Vector,
    factorized_: bool,
    maxiter_: Int,
    iter_: Int,
}

impl KKTSolverDiag {
    /// HiGHS: ipx/kkt_solver_diag.cc::KKTSolverDiag::KKTSolverDiag (lines 10-16)
    pub fn new(control: Rc<Control>, model: Rc<Model>) -> Self {
        let m = model.rows();
        let n = model.cols();
        Self {
            control_: control,
            normal_matrix_: NormalMatrix::new(model.clone()),
            precond_: DiagonalPrecond::new(model.clone()),
            model_: model,
            W_: Vector::new(m + n),
            resscale_: Vector::new(m),
            factorized_: false,
            maxiter_: -1,
            iter_: 0,
        }
    }
    pub fn maxiter(&self) -> Int {
        self.maxiter_
    }
    /// HiGHS: ipx/kkt_solver_diag.h::KKTSolverDiag::maxiter(Int)
    pub fn set_maxiter(&mut self, new_maxiter: Int) {
        self.maxiter_ = new_maxiter;
    }
}

impl KKTSolver for KKTSolverDiag {
    /// HiGHS: ipx/kkt_solver_diag.cc::KKTSolverDiag::_Factorize (lines 18-65)
    fn _factorize(&mut self, iterate: Option<&mut Iterate>, info: &mut Info) {
        let m = self.model_.rows();
        let n = self.model_.cols();
        self.iter_ = 0;
        self.factorized_ = false;

        if let Some(pt) = iterate {
            let xl = pt.xl();
            let xu = pt.xu();
            let zl = pt.zl();
            let zu = pt.zu();

            // Build matrix W for AI*W*AI'. For free variables set W[j] to
            // 1.0/regval, where regval is a regularization value. regval is chosen
            // as the minimum of the complementarity measure mu and the smallest
            // nonzero diagonal entry of the (1,1) block of the KKT matrix.
            let mut regval = pt.mu();
            for j in 0..n + m {
                let ju = j as usize;
                debug_assert!(xl[ju] > 0.0);
                debug_assert!(xu[ju] > 0.0);
                let g = zl[ju] / xl[ju] + zu[ju] / xu[ju];
                debug_assert!(g.is_finite());
                if g != 0.0 && g < regval {
                    regval = g;
                }
                self.W_[ju] = 1.0 / g; // infinity if g is zero
            }
            for j in 0..n + m {
                let ju = j as usize;
                if self.W_[ju].is_infinite() {
                    self.W_[ju] = 1.0 / regval;
                }
                debug_assert!(self.W_[ju].is_finite());
                debug_assert!(self.W_[ju] > 0.0);
            }
        } else {
            self.W_.fill(1.0);
        }

        // Residual scaling factors for termination test of CR method (see below).
        for i in 0..m {
            self.resscale_[i as usize] = 1.0 / self.W_[(n + i) as usize].sqrt();
        }

        // Build normal matrix and preconditioner.
        self.normal_matrix_.prepare(Some(&self.W_[..]));
        self.precond_.factorize(Some(&self.W_[..]), info);
        if info.errflag != 0 {
            return;
        }

        self.factorized_ = true;
    }

    /// HiGHS: ipx/kkt_solver_diag.cc::KKTSolverDiag::_Solve (lines 82-118)
    ///
    /// Reduces the KKT system
    ///
    ///   [ W^{-1}  AI' ] (x) = (a) + (res)
    ///   [ AI       0  ] (y)   (b)   ( 0 )
    ///
    /// to normal equations
    ///
    ///   C * y := (AI*W*AI') * y = -b + AI*W*(a+res)
    ///
    /// and solves by the CR method. The solution to the KKT system is recovered so
    /// that the first n entries of res are zero.
    fn _solve(
        &mut self,
        a: &Vector,
        b: &Vector,
        tol: f64,
        x: &mut Vector,
        y: &mut Vector,
        info: &mut Info,
    ) {
        let m = self.model_.rows();
        let n = self.model_.cols();
        let AI = self.model_.AI();
        debug_assert!(self.factorized_);

        // Compose right-hand side AI*W*a-b.
        let mut rhs = Vector::new(b.size());
        for i in 0..b.size() as usize {
            rhs[i] = -b[i];
        }
        for j in 0..n + m {
            scatter_column(AI, j, self.W_[j as usize] * a[j as usize], &mut rhs);
        }

        // Solve normal equations.
        y.fill(0.0);
        self.normal_matrix_.reset_time();
        self.precond_.reset_time();
        let mut cr = ConjugateResiduals::new(self.control_.clone());
        cr.solve_precond(
            &mut self.normal_matrix_,
            &mut self.precond_,
            &rhs,
            tol,
            Some(&self.resscale_[..]),
            self.maxiter_,
            y,
        );
        info.errflag = cr.errflag();
        info.kktiter1 += cr.iter();
        info.time_cr1 += cr.time();
        info.time_cr1_AAt += self.normal_matrix_.time();
        info.time_cr1_pre += self.precond_.time();
        self.iter_ += cr.iter();

        // Recover solution to KKT system.
        for i in 0..m {
            x[(n + i) as usize] = b[i as usize];
        }
        for j in 0..n {
            let ju = j as usize;
            let aty = dot_column(AI, j, y);
            x[ju] = self.W_[ju] * (a[ju] - aty);
            for p in AI.begin(j)..AI.end(j) {
                let i = (n as usize) + AI.index(p) as usize;
                x[i] -= x[ju] * AI.value(p);
            }
        }
    }

    fn _iter(&self) -> Int {
        self.iter_
    }
}
