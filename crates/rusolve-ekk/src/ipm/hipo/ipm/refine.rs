//! HiGHS source: ipm/hipo/ipm/Refine.cpp, ipm/hipo/ipm/Solver.cpp (lines 322-338)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Iterative refinement of the 6x6 Newton system (`Solver::refine`) and its
//! backward-error estimate `computeOmega` (Arioli-Demmel-Duff
//! componentwise / normwise splitting), together with their only caller
//! `Solver::solveNewtonSystem` (Solver.cpp 322-338). It sits here, not in
//! `solver/newton.rs`, next to `refine` and `solve6x6`, which it calls.

use super::iterate::NewtonDir;
use super::solver::Solver;
use crate::ipm::hipo::auxiliary::Clock;
use crate::ipm::hipo::auxiliary::vector_operations::infNorm;
use crate::ipm::hipo::ipm::parameters::{kMaxIterRefine, kTolRefine};
use crate::ipm::hipo::ipm::status::kStatusError;
use crate::util::highs_utils::cpp_max;

impl Solver {
    /// HiGHS: Solver.cpp::Solver::solveNewtonSystem (lines 322-338)
    pub fn solveNewtonSystem(&mut self, delta: &mut NewtonDir) -> bool {
        // HiGHS: Solver.cpp line 323: solve6x6(delta, it_->res);
        // Rule (i): take res out to pass as rhs, then put it back.
        {
            let res = std::mem::take(&mut self.it_.as_mut().unwrap().res);
            let _g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(25);
            self.solve6x6(delta, &res);
            drop(_g);
            self.it_.as_mut().unwrap().res = res;
        }

        // HiGHS: Solver.cpp line 324: refine(delta);
        self.refine(delta);

        // HiGHS: Solver.cpp lines 326-335: Check for NaN of Inf
        {
            let it = self.it_.as_ref().unwrap();
            if it.isDirNan(delta) {
                // HiGHS: Solver.cpp line 328
                self.logger_.printInfo("Direction is nan\n");
                self.info_.status = kStatusError;
                return true;
            } else if it.isDirInf(delta) {
                // HiGHS: Solver.cpp line 332
                self.logger_.printInfo("Direction is inf\n");
                self.info_.status = kStatusError;
                return true;
            }
        }

        false
    }

    /// HiGHS: Refine.cpp::Solver::refine (lines 7-53)
    pub fn refine(&mut self, delta: &mut NewtonDir) {
        // HiGHS: Refine.cpp line 8: Clock clock;
        let mut clock = Clock::new();

        // HiGHS: Refine.cpp lines 10-11
        let mut correction = NewtonDir::new(self.m_, self.n_);
        let mut temp: NewtonDir;

        // HiGHS: Refine.cpp lines 13-16: compute the residuals of the linear system in it_->ires
        clock.start();
        {
            // HiGHS: Refine.cpp line 15: it_->residuals6x6(delta);
            // model_/regul_ are passed explicitly (Iterate holds no back-references).
            let it = self.it_.as_mut().unwrap();
            let _g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(26);
            it.residuals6x6(&self.model_, &self.regul_, delta);
        }
        // HiGHS: Refine.cpp line 16: info_.residual_time += clock.stop();
        self.info_.residual_time += clock.stop();

        // HiGHS: Refine.cpp lines 18-20
        clock.start();
        let _g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(27);
        let mut omega = self.computeOmega(delta);
        drop(_g);
        self.info_.omega_time += clock.stop();

        // HiGHS: Refine.cpp line 22: double old_omega{};
        let mut old_omega: f64;

        // HiGHS: Refine.cpp line 24: for (Int iter = 0; iter < kMaxIterRefine; ++iter)
        let mut iter: crate::ipm::hipo::Int = 0;
        while iter < kMaxIterRefine {
            // HiGHS: Refine.cpp line 25: if (omega < kTolRefine) break;
            if omega < kTolRefine {
                break;
            }

            // HiGHS: Refine.cpp line 27: correction.clear();
            correction.clear();

            // HiGHS: Refine.cpp line 28: solve6x6(correction, it_->ires);
            // Rule (i): take ires out to pass as rhs, then put it back.
            {
                let ires = std::mem::take(&mut self.it_.as_mut().unwrap().ires);
                let _g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(25);
                self.solve6x6(&mut correction, &ires);
                drop(_g);
                self.it_.as_mut().unwrap().ires = ires;
            }

            // HiGHS: Refine.cpp line 30: temp = delta;
            let _g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(28);
            temp = delta.clone();

            // HiGHS: Refine.cpp line 31: temp.add(correction);
            temp.add(&correction);
            drop(_g);

            // HiGHS: Refine.cpp lines 33-35
            clock.start();
            {
                // HiGHS: Refine.cpp line 34: it_->residuals6x6(temp);
                let it = self.it_.as_mut().unwrap();
                let _g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(26);
                it.residuals6x6(&self.model_, &self.regul_, &temp);
            }
            self.info_.residual_time += clock.stop();

            // HiGHS: Refine.cpp line 37: old_omega = omega;
            old_omega = omega;

            // HiGHS: Refine.cpp lines 39-41
            clock.start();
            let _g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(27);
            omega = self.computeOmega(&temp);
            drop(_g);
            self.info_.omega_time += clock.stop();

            // HiGHS: Refine.cpp lines 43-48
            if omega < old_omega {
                // HiGHS: Refine.cpp line 44: delta = temp;
                *delta = temp.clone();
            } else {
                // HiGHS: Refine.cpp line 46: omega = old_omega;
                omega = old_omega;
                // HiGHS: Refine.cpp line 47: break;
                break;
            }

            // HiGHS: Refine.cpp line 24: ++iter
            iter += 1;
        }

        // HiGHS: Refine.cpp line 52: it_->data.back().omega = std::max(it_->data.back().omega, omega);
        {
            let it = self.it_.as_mut().unwrap();
            let back = it.data.back();
            back.omega = cpp_max(back.omega, omega);
        }
    }

    /// HiGHS: Refine.cpp::Solver::computeOmega (lines 66-212)
    pub fn computeOmega(&self, delta: &NewtonDir) -> f64 {
        // HiGHS: Refine.cpp lines 67-68: comment (dropped)

        // HiGHS: Refine.cpp lines 70-72: Compute |A| * |dx| and |A^T| * |dy|
        let m = self.m_;
        let n = self.n_;
        let mut abs_prod_A: Vec<f64> = vec![0.0; m as usize];
        let mut abs_prod_At: Vec<f64> = vec![0.0; n as usize];

        // HiGHS: Refine.cpp lines 73-81
        let a = self.model_.A();
        for col in 0..n {
            let start_col = a.start_[col as usize] as usize;
            let start_col1 = a.start_[(col + 1) as usize] as usize;
            let mut el = start_col;
            while el < start_col1 {
                let row = a.index_[el] as usize;
                let val = a.value_[el];
                abs_prod_A[row] += val.abs() * delta.x[col as usize].abs();
                abs_prod_At[col as usize] += val.abs() * delta.y[row].abs();
                el += 1;
            }
        }

        // HiGHS: Refine.cpp lines 83-89
        let mut inf_norm_delta: f64 = 0.0;
        inf_norm_delta = cpp_max(inf_norm_delta, infNorm(&delta.x));
        inf_norm_delta = cpp_max(inf_norm_delta, infNorm(&delta.xl));
        inf_norm_delta = cpp_max(inf_norm_delta, infNorm(&delta.xu));
        inf_norm_delta = cpp_max(inf_norm_delta, infNorm(&delta.y));
        inf_norm_delta = cpp_max(inf_norm_delta, infNorm(&delta.zl));
        inf_norm_delta = cpp_max(inf_norm_delta, infNorm(&delta.zu));

        // HiGHS: Refine.cpp lines 91-95: comments (dropped)

        // HiGHS: Refine.cpp line 97: const double tau_const = 1000.0 * (5 * n_ + m_) * 1e-16;
        // C++: 5 * n_ is int arithmetic, then + m_ is int, then promoted to double.
        let tau_const: f64 = 1000.0 * ((5 * self.n_ + self.m_) as f64) * 1e-16;

        // HiGHS: Refine.cpp line 98: double omega_1{}, omega_2{};
        let mut omega_1: f64 = 0.0;
        let mut omega_2: f64 = 0.0;

        // HiGHS: Refine.cpp line 100: assert(it_->Rd || m_ == 0);
        let it = self.it_.as_ref().unwrap();
        debug_assert!(it.Rd.is_some() || m == 0);

        // HiGHS: Refine.cpp lines 102-116: First block
        let rd = it.Rd.unwrap_or(0);
        for i in 0..m {
            let i_us = i as usize;
            // HiGHS: Refine.cpp lines 104-107
            let tau: f64 = tau_const
                * (cpp_max(self.model_.infNormRows(i), it.total_reg[rd + i_us].abs())
                    * inf_norm_delta
                    + it.res.r1[i_us].abs());

            // HiGHS: Refine.cpp lines 108-115
            updateOmega(
                tau,
                &mut omega_1,
                &mut omega_2,
                // numerator: std::abs(it_->ires.r1[i])
                it.ires.r1[i_us].abs(),
                // den1: abs_prod_A[i] + std::abs(it_->Rd[i]) * std::abs(delta.y[i])
                abs_prod_A[i_us] + it.total_reg[rd + i_us].abs() * delta.y[i_us].abs(),
                // den2: std::abs(it_->res.r1[i])
                it.res.r1[i_us].abs(),
                // den3: (model_.oneNormRows(i) + std::abs(it_->Rd[i])) * inf_norm_delta
                (self.model_.oneNormRows(i) + it.total_reg[rd + i_us].abs()) * inf_norm_delta,
            );
        }

        // HiGHS: Refine.cpp lines 118-132: Second block
        for i in 0..n {
            if self.model_.hasLb(i) {
                let i_us = i as usize;
                // HiGHS: Refine.cpp lines 121-122
                let tau: f64 = tau_const * (1.0 * inf_norm_delta + it.res.r2[i_us].abs());

                // HiGHS: Refine.cpp lines 123-130
                updateOmega(
                    tau,
                    &mut omega_1,
                    &mut omega_2,
                    // numerator
                    it.ires.r2[i_us].abs(),
                    // den1: std::abs(delta.x[i]) + std::abs(delta.xl[i])
                    delta.x[i_us].abs() + delta.xl[i_us].abs(),
                    // den2: std::abs(it_->res.r2[i])
                    it.res.r2[i_us].abs(),
                    // den3: 2 * inf_norm_delta
                    2.0 * inf_norm_delta,
                );
            }
        }

        // HiGHS: Refine.cpp lines 134-146: Third block
        for i in 0..n {
            if self.model_.hasUb(i) {
                let i_us = i as usize;
                // HiGHS: Refine.cpp lines 137-138
                let tau: f64 = tau_const * (1.0 * inf_norm_delta + it.res.r3[i_us].abs());

                // HiGHS: Refine.cpp lines 139-144
                updateOmega(
                    tau,
                    &mut omega_1,
                    &mut omega_2,
                    // numerator
                    it.ires.r3[i_us].abs(),
                    // den1: std::abs(delta.x[i]) + std::abs(delta.xu[i])
                    delta.x[i_us].abs() + delta.xu[i_us].abs(),
                    // den2: std::abs(it_->res.r3[i])
                    it.res.r3[i_us].abs(),
                    // den3: 2 * inf_norm_delta
                    2.0 * inf_norm_delta,
                );
            }
        }

        // HiGHS: Refine.cpp lines 148-169: Fourth block
        let rp = it.Rp;
        for i in 0..n {
            let i_us = i as usize;
            // HiGHS: Refine.cpp line 150: double reg_p = it_->Rp ? it_->Rp[i] : regul_.primal;
            let reg_p: f64 = match rp {
                Some(rp_off) => it.total_reg[rp_off + i_us],
                None => self.regul_.primal,
            };

            // HiGHS: Refine.cpp lines 152-156
            let tau: f64 = tau_const
                * (cpp_max(cpp_max(self.model_.infNormCols(i), 1.0), reg_p.abs()) * inf_norm_delta
                    + it.res.r4[i_us].abs());

            // HiGHS: Refine.cpp lines 158-161
            let mut denom1: f64 = abs_prod_At[i_us];
            if self.model_.hasLb(i) {
                denom1 += delta.zl[i_us].abs();
            }
            if self.model_.hasUb(i) {
                denom1 += delta.zu[i_us].abs();
            }
            denom1 += reg_p.abs() * delta.x[i_us].abs();

            // HiGHS: Refine.cpp lines 163-168
            updateOmega(
                tau,
                &mut omega_1,
                &mut omega_2,
                // numerator
                it.ires.r4[i_us].abs(),
                // den1
                denom1,
                // den2: std::abs(it_->res.r4[i])
                it.res.r4[i_us].abs(),
                // den3: (model_.oneNormCols(i) + 2 + std::abs(reg_p)) * inf_norm_delta
                (self.model_.oneNormCols(i) + 2.0 + reg_p.abs()) * inf_norm_delta,
            );
        }

        // HiGHS: Refine.cpp lines 171-189: Fifth block
        for i in 0..n {
            if self.model_.hasLb(i) {
                let i_us = i as usize;
                // HiGHS: Refine.cpp lines 174-177
                let tau: f64 = tau_const
                    * (cpp_max(it.xl[i_us].abs(), it.zl[i_us].abs()) * inf_norm_delta
                        + it.res.r5[i_us].abs());

                // HiGHS: Refine.cpp lines 179-187
                updateOmega(
                    tau,
                    &mut omega_1,
                    &mut omega_2,
                    // numerator
                    it.ires.r5[i_us].abs(),
                    // den1: std::abs(it_->zl[i]) * std::abs(delta.xl[i]) + std::abs(it_->xl[i]) * std::abs(delta.zl[i])
                    it.zl[i_us].abs() * delta.xl[i_us].abs()
                        + it.xl[i_us].abs() * delta.zl[i_us].abs(),
                    // den2: std::abs(it_->res.r5[i])
                    it.res.r5[i_us].abs(),
                    // den3: (std::abs(it_->zl[i]) + std::abs(it_->xl[i])) * inf_norm_delta
                    (it.zl[i_us].abs() + it.xl[i_us].abs()) * inf_norm_delta,
                );
            }
        }

        // HiGHS: Refine.cpp lines 191-209: Sixth block
        for i in 0..n {
            if self.model_.hasUb(i) {
                let i_us = i as usize;
                // HiGHS: Refine.cpp lines 194-197
                let tau: f64 = tau_const
                    * (cpp_max(it.xu[i_us].abs(), it.zu[i_us].abs()) * inf_norm_delta
                        + it.res.r6[i_us].abs());

                // HiGHS: Refine.cpp lines 199-207
                updateOmega(
                    tau,
                    &mut omega_1,
                    &mut omega_2,
                    // numerator
                    it.ires.r6[i_us].abs(),
                    // den1: std::abs(it_->zu[i]) * std::abs(delta.xu[i]) + std::abs(it_->xu[i]) * std::abs(delta.zu[i])
                    it.zu[i_us].abs() * delta.xu[i_us].abs()
                        + it.xu[i_us].abs() * delta.zu[i_us].abs(),
                    // den2: std::abs(it_->res.r6[i])
                    it.res.r6[i_us].abs(),
                    // den3: (std::abs(it_->zu[i]) + std::abs(it_->xu[i])) * inf_norm_delta
                    (it.zu[i_us].abs() + it.xu[i_us].abs()) * inf_norm_delta,
                );
            }
        }

        // HiGHS: Refine.cpp line 211: return omega_1 + omega_2;
        omega_1 + omega_2
    }
}

/// HiGHS: Refine.cpp::updateOmega (lines 55-64)
pub fn updateOmega(
    tau: f64,
    omega1: &mut f64,
    omega2: &mut f64,
    num: f64,
    den1: f64,
    den2: f64,
    den3: f64,
) {
    // HiGHS: Refine.cpp line 57: if (den1 + den2 > tau)
    if den1 + den2 > tau {
        // HiGHS: Refine.cpp line 58: double omega = num / (den1 + den2);
        let omega = num / (den1 + den2);
        // HiGHS: Refine.cpp line 59: omega1 = std::max(omega1, omega);
        *omega1 = cpp_max(*omega1, omega);
    } else {
        // HiGHS: Refine.cpp line 61: double omega = num / (den1 + den3);
        let omega = num / (den1 + den3);
        // HiGHS: Refine.cpp line 62: omega2 = std::max(omega2, omega);
        *omega2 = cpp_max(*omega2, omega);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn update_omega_branches() {
        let (mut o1, mut o2) = (0.0, 0.0);
        updateOmega(1.0, &mut o1, &mut o2, 4.0, 1.0, 1.0, 3.0);
        assert_eq!((o1, o2), (2.0, 0.0));
        updateOmega(10.0, &mut o1, &mut o2, 4.0, 1.0, 1.0, 3.0);
        assert_eq!((o1, o2), (2.0, 1.0));
    }
}
