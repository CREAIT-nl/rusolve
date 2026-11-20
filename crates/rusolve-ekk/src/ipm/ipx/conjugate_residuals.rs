//! HiGHS source: ipm/ipx/conjugate_residuals.h, conjugate_residuals.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::control::{Control, sci2};
use super::linear_operator::LinearOperator;
use super::timer::Timer;
use super::utils::{dot, infnorm};
use super::*;
use crate::util::highs_utils::cpp_max;
use std::rc::Rc;

/// HiGHS: ipx/conjugate_residuals.h::ConjugateResiduals
pub struct ConjugateResiduals {
    control_: Rc<Control>,
    errflag_: Int,
    iter_: Int,
    time_: f64,
}

impl ConjugateResiduals {
    /// HiGHS: ipx/conjugate_residuals.cc::ConjugateResiduals::ConjugateResiduals (line 9)
    pub fn new(control: Rc<Control>) -> Self {
        Self {
            control_: control,
            errflag_: 0,
            iter_: 0,
            time_: 0.0,
        }
    }
    /// HiGHS: ipx/conjugate_residuals.cc::ConjugateResiduals::Solve(C, rhs, ...) (lines 12-86)
    ///
    /// `resscale` may be `nullptr`.
    pub fn solve(
        &mut self,
        C: &mut dyn LinearOperator,
        rhs: &Vector,
        tol: f64,
        resscale: Option<&[f64]>,
        maxiter: Int,
        lhs: &mut Vector,
    ) {
        let m = rhs.size();
        let mut residual = Vector::new(m); // rhs - C*lhs
        let mut step = Vector::new(m); // update to lhs
        let mut Cresidual = Vector::new(m); // C * residual
        let mut Cstep = Vector::new(m); // C * step
        let mut cdot = 0.0; // dot product from C.Apply
        let timer = Timer::default();

        self.errflag_ = 0;
        self.iter_ = 0;
        self.time_ = 0.0;
        let mut maxiter = maxiter;
        if maxiter < 0 {
            maxiter = m + 100;
        }

        // Initialize residual, step and Cstep.
        if infnorm(lhs) == 0.0 {
            residual.assign(rhs); // saves a matrix-vector op
        } else {
            C.apply(lhs, &mut residual, None);
            for i in 0..m as usize {
                residual[i] = rhs[i] - residual[i];
            }
        }
        C.apply(&residual, &mut Cresidual, Some(&mut cdot));
        step.assign(&residual);
        Cstep.assign(&Cresidual);

        loop {
            // Termination check.
            let mut resnorm = 0.0;
            match resscale {
                Some(resscale) => {
                    for i in 0..m as usize {
                        resnorm = cpp_max(resnorm, (resscale[i] * residual[i]).abs());
                    }
                }
                None => resnorm = infnorm(&residual),
            }
            if resnorm <= tol {
                break;
            }
            if self.iter_ == maxiter {
                self.control_.debug(
                    3,
                    &format!(
                        " CR method not converged in {} iterations. residual = {}, tolerance = {}\n",
                        maxiter,
                        sci2(resnorm),
                        sci2(tol)
                    ),
                );
                self.errflag_ = IPX_ERROR_cr_iter_limit;
                break;
            }
            if cdot <= 0.0 {
                self.errflag_ = IPX_ERROR_cr_matrix_not_posdef;
                break;
            }

            // Update lhs, residual and Cresidual.
            let denom = dot(&Cstep, &Cstep);
            let alpha = cdot / denom;
            if !alpha.is_finite() {
                self.errflag_ = IPX_ERROR_cr_inf_or_nan;
                break;
            }
            for i in 0..m as usize {
                lhs[i] = lhs[i] + alpha * step[i];
            }
            for i in 0..m as usize {
                residual[i] = residual[i] - alpha * Cstep[i];
            }
            let mut cdotnew = 0.0;
            C.apply(&residual, &mut Cresidual, Some(&mut cdotnew));

            // Update step and Cstep.
            let beta = cdotnew / cdot;
            for i in 0..m as usize {
                step[i] = residual[i] + beta * step[i];
            }
            for i in 0..m as usize {
                Cstep[i] = Cresidual[i] + beta * Cstep[i];
            }
            cdot = cdotnew;

            self.iter_ += 1;
            self.errflag_ = self.control_.interrupt_check(-1);
            if self.errflag_ != 0 {
                break;
            }
        }
        self.time_ = timer.elapsed();
    }

    /// HiGHS: ipx/conjugate_residuals.cc::ConjugateResiduals::Solve(C, P, rhs, ...) (lines 88-211)
    pub fn solve_precond(
        &mut self,
        C: &mut dyn LinearOperator,
        P: &mut dyn LinearOperator,
        rhs: &Vector,
        tol: f64,
        resscale: Option<&[f64]>,
        maxiter: Int,
        lhs: &mut Vector,
    ) {
        let m = rhs.size();
        let mut residual = Vector::new(m); // rhs - C*lhs
        let mut sresidual = Vector::new(m); // preconditioned residual
        let mut step = Vector::new(m); // update to lhs
        let mut Csresidual = Vector::new(m); // C * sresidual
        let mut Cstep = Vector::new(m); // C * step
        let mut cdot = 0.0; // dot product from C.Apply
        let timer = Timer::default();

        // resnorm_precond_system = residual' * P * residual.
        // This quantity is minimized by the preconditioned CR method over a Krylov
        // subspace. Hence (in theory) must decrease strictly monotonically. If it
        // does not, then the method stagnated due to round-off errors. This
        // happened in a few cases with augmented diagonal preconditioning (i.e.
        // diagonal preconditioning with dense columns treated as low-rank update)
        // because operations with P were not sufficiently accurate.
        let mut resnorm_precond_system = 0.0;

        self.errflag_ = 0;
        self.iter_ = 0;
        self.time_ = 0.0;
        let mut maxiter = maxiter;
        if maxiter < 0 {
            maxiter = m + 100;
        }

        // Initialize residual, sresidual, step and Cstep.
        if infnorm(lhs) == 0.0 {
            residual.assign(rhs); // saves a matrix-vector op
        } else {
            C.apply(lhs, &mut residual, None);
            for i in 0..m as usize {
                residual[i] = rhs[i] - residual[i];
            }
        }
        P.apply(&residual, &mut sresidual, Some(&mut resnorm_precond_system));
        C.apply(&sresidual, &mut Csresidual, Some(&mut cdot));
        step.assign(&sresidual);
        Cstep.assign(&Csresidual);

        loop {
            // Termination check.
            let mut resnorm = 0.0;
            match resscale {
                Some(resscale) => {
                    for (&rs, &ri) in resscale[..m as usize].iter().zip(residual.iter()) {
                        resnorm = cpp_max(resnorm, (rs * ri).abs());
                    }
                }
                None => resnorm = infnorm(&residual),
            }
            if resnorm <= tol {
                break;
            }
            if self.iter_ == maxiter {
                self.control_.debug(
                    3,
                    &format!(
                        " PCR method not converged in {} iterations. residual = {}, tolerance = {}\n",
                        maxiter,
                        sci2(resnorm),
                        sci2(tol)
                    ),
                );
                self.errflag_ = IPX_ERROR_cr_iter_limit;
                break;
            }
            if cdot <= 0.0 {
                self.control_.debug(
                    3,
                    &format!(
                        " matrix in PCR method not posdef. cdot = {}, infnorm(sresidual) = {}, infnorm(residual) = {}\n",
                        sci2(cdot),
                        sci2(infnorm(&sresidual)),
                        sci2(infnorm(&residual))
                    ),
                );
                self.errflag_ = IPX_ERROR_cr_matrix_not_posdef;
                break;
            }

            // Update lhs, residual, sresidual and Csresidual.
            let mut cdotnew = 0.0;
            {
                // Uses Csresidual as storage for preconditioned Cstep.
                // C++ `Vector& precond_Cstep = Csresidual;`
                let mut pdot = 0.0;
                P.apply(&Cstep, &mut Csresidual, Some(&mut pdot));
                if pdot <= 0.0 {
                    self.errflag_ = IPX_ERROR_cr_precond_not_posdef;
                    break;
                }
                let alpha = cdot / pdot;
                if !alpha.is_finite() {
                    self.errflag_ = IPX_ERROR_cr_inf_or_nan;
                    break;
                }
                axpy_in_place(&mut lhs[..m as usize], alpha, &step);
                axmy_in_place(&mut residual, alpha, &Cstep);
                axmy_in_place(&mut sresidual, alpha, &Csresidual);
                C.apply(&sresidual, &mut Csresidual, Some(&mut cdotnew));
                // Now Csresidual is restored and alias goes out of scope.
            }

            // Update step and Cstep.
            let beta = cdotnew / cdot;
            xpby_in_place(&mut step, &sresidual, beta);
            xpby_in_place(&mut Cstep, &Csresidual, beta);
            cdot = cdotnew;

            self.iter_ += 1;
            if self.iter_ % 5 == 0 {
                // We recompute the preconditioned residual from its definition
                // all 5 iterations. Otherwise it happened that the preconditioned
                // residual approached zero but the true residual stagnated. This
                // was caused by the update to the preconditioned residual being not
                // sufficiently accurate.
                // As a second safeguard, we check that resnorm_precond_system
                // decreased during the last 5 iterations.
                let mut rsdot = 0.0;
                P.apply(&residual, &mut sresidual, Some(&mut rsdot));
                if rsdot >= resnorm_precond_system {
                    self.control_.debug(
                        3,
                        &format!(
                            " resnorm_precond_system old = {}\n resnorm_precond_system new = {}\n",
                            sci2(resnorm_precond_system),
                            sci2(rsdot)
                        ),
                    );
                    self.errflag_ = IPX_ERROR_cr_no_progress;
                    break;
                }
                resnorm_precond_system = rsdot;
            }

            self.errflag_ = self.control_.interrupt_check(-1);
            if self.errflag_ != 0 {
                break;
            }
        }
        self.time_ = timer.elapsed();
    }

    /// HiGHS: ipx/conjugate_residuals.cc::ConjugateResiduals::errflag (line 213)
    pub fn errflag(&self) -> Int {
        self.errflag_
    }
    /// HiGHS: ipx/conjugate_residuals.cc::ConjugateResiduals::iter (line 214)
    pub fn iter(&self) -> Int {
        self.iter_
    }
    /// HiGHS: ipx/conjugate_residuals.cc::ConjugateResiduals::time (line 215)
    pub fn time(&self) -> f64 {
        self.time_
    }
}

// PORT NOTE (perf): the valarray expressions of ConjugateResiduals::Solve
// (`lhs += alpha*step`, `residual -= alpha*Cstep`, `step = sresidual + beta*step`)
// as zipped slice loops. Indexed `for i in 0..m` loops paid a bounds check per
// operand per element and did not vectorise; each element is still computed
// with the same two operations in the same order (no contraction), so the
// results are bit-identical.

/// `x += a*y` (element-wise `x[i] + a*y[i]`).
#[inline]
fn axpy_in_place(x: &mut [f64], a: f64, y: &[f64]) {
    assert_eq!(x.len(), y.len());
    for (xi, &yi) in x.iter_mut().zip(y) {
        *xi = *xi + a * yi;
    }
}

/// `x -= a*y` (element-wise `x[i] - a*y[i]`).
#[inline]
fn axmy_in_place(x: &mut [f64], a: f64, y: &[f64]) {
    assert_eq!(x.len(), y.len());
    for (xi, &yi) in x.iter_mut().zip(y) {
        *xi = *xi - a * yi;
    }
}

/// `y = x + b*y` (element-wise `x[i] + b*y[i]`).
#[inline]
fn xpby_in_place(y: &mut [f64], x: &[f64], b: f64) {
    assert_eq!(x.len(), y.len());
    for (yi, &xi) in y.iter_mut().zip(x) {
        *yi = xi + b * *yi;
    }
}
