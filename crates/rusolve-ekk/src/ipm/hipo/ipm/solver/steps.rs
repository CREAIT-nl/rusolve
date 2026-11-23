//! HiGHS source: ipm/hipo/ipm/Solver.cpp (lines 491-808)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Step sizes (Mehrotra), the step itself and the starting point.

use super::Solver;
use crate::ipm::hipo::auxiliary::vector_operations::{infNorm, norm2, vectorAdd, vectorAdd_scalar};
use crate::ipm::hipo::ipm::linear_solver::LsCtx;
use crate::lp_data::highs_options::{kHipoAugmentedString, kHipoNormalEqString};
use crate::util::highs_utils::{cpp_max, cpp_min};

impl Solver {
    /// HiGHS: Solver.cpp::Solver::stepSizes (lines 491-575)
    pub fn stepSizes(&mut self) {
        // Compute primal and dual stepsizes.
        // std::vector<double>& xl = it_->xl;
        // std::vector<double>& xu = it_->xu;
        // std::vector<double>& zl = it_->zl;
        // std::vector<double>& zu = it_->zu;
        // std::vector<double>& dxl = it_->delta.xl;
        // std::vector<double>& dxu = it_->delta.xu;
        // std::vector<double>& dzl = it_->delta.zl;
        // std::vector<double>& dzu = it_->delta.zu;

        // parameters for Mehrotra heuristic
        let gamma_f: f64 = 0.9;
        let gamma_a: f64 = 1.0 / (1.0 - gamma_f);

        // compute stepsizes and blocking components
        let mut block_xl: i32 = 0;
        let mut block_xu: i32 = 0;
        let mut block_zl: i32 = 0;
        let mut block_zu: i32 = 0;

        let mut it = self.it_.take().expect("it_ is created by solve");

        let alpha_xl =
            self.stepToBoundary(&it.xl, &it.delta.xl, None, 0.0, true, Some(&mut block_xl));
        let alpha_xu =
            self.stepToBoundary(&it.xu, &it.delta.xu, None, 0.0, false, Some(&mut block_xu));
        let alpha_zl =
            self.stepToBoundary(&it.zl, &it.delta.zl, None, 0.0, true, Some(&mut block_zl));
        let alpha_zu =
            self.stepToBoundary(&it.zu, &it.delta.zu, None, 0.0, false, Some(&mut block_zu));

        let max_p = cpp_min(alpha_xl, alpha_xu);
        let max_d = cpp_min(alpha_zl, alpha_zu);

        // compute mu with current stepsizes
        let mut mu_full: f64 = 0.0;
        let mut num_finite: i32 = 0;
        for i in 0..self.n_ {
            if self.model_.hasLb(i) {
                mu_full += (it.xl[i as usize] + max_p * it.delta.xl[i as usize])
                    * (it.zl[i as usize] + max_d * it.delta.zl[i as usize]);
                num_finite += 1;
            }
            if self.model_.hasUb(i) {
                mu_full += (it.xu[i as usize] + max_p * it.delta.xu[i as usize])
                    * (it.zu[i as usize] + max_d * it.delta.zu[i as usize]);
                num_finite += 1;
            }
        }
        mu_full /= num_finite as f64;
        mu_full /= gamma_a;

        // compute new stepsizes based on Mehrotra heuristic

        // primal
        let mut alpha_p: f64 = 1.0;
        if max_p < 1.0 {
            let block_p: i32;
            if alpha_xl <= alpha_xu {
                block_p = block_xl;
                let temp =
                    mu_full / (it.zl[block_p as usize] + max_d * it.delta.zl[block_p as usize]);
                alpha_p = (temp - it.xl[block_p as usize]) / it.delta.xl[block_p as usize];
            } else {
                block_p = block_xu;
                let temp =
                    mu_full / (it.zu[block_p as usize] + max_d * it.delta.zu[block_p as usize]);
                alpha_p = (temp - it.xu[block_p as usize]) / it.delta.xu[block_p as usize];
            }
            alpha_p = cpp_max(alpha_p, gamma_f * max_p);
            alpha_p = cpp_min(alpha_p, 1.0);
            debug_assert!(block_p >= 0);
        }

        // dual
        let mut alpha_d: f64 = 1.0;
        if max_d < 1.0 {
            let block_d: i32;
            if alpha_zl <= alpha_zu {
                block_d = block_zl;
                let temp =
                    mu_full / (it.xl[block_d as usize] + max_p * it.delta.xl[block_d as usize]);
                alpha_d = (temp - it.zl[block_d as usize]) / it.delta.zl[block_d as usize];
            } else {
                block_d = block_zu;
                let temp =
                    mu_full / (it.xu[block_d as usize] + max_p * it.delta.xu[block_d as usize]);
                alpha_d = (temp - it.zu[block_d as usize]) / it.delta.zu[block_d as usize];
            }
            alpha_d = cpp_max(alpha_d, gamma_f * max_d);
            alpha_d = cpp_min(alpha_d, 1.0);
            debug_assert!(block_d >= 0);
        }

        self.alpha_primal_ = cpp_min(alpha_p, 1.0 - 1e-4);
        self.alpha_dual_ = cpp_min(alpha_d, 1.0 - 1e-4);

        debug_assert!(
            self.alpha_primal_ > 0.0
                && self.alpha_primal_ < 1.0
                && self.alpha_dual_ > 0.0
                && self.alpha_dual_ < 1.0
        );

        self.it_ = Some(it);
    }

    /// HiGHS: Solver.cpp::Solver::makeStep (lines 577-584)
    pub fn makeStep(&mut self) {
        self.stepSizes();
        let mut it = self.it_.take().expect("it_ is created by solve");
        it.makeStep(self.alpha_primal_, self.alpha_dual_);
        it.residual1234(&self.model_);
        it.computeMu(&self.model_);
        it.indicators(&self.model_);
        self.it_ = Some(it);
        self.printOutput();
    }

    /// HiGHS: Solver.cpp::Solver::startingPoint (lines 586-808)
    pub fn startingPoint(&mut self) -> bool {
        // Return true if an error occurred

        let mut it = self.it_.take().expect("it_ is created by solve");

        // *********************************************************************
        // x starting point
        // *********************************************************************
        // compute feasible x
        for i in 0..self.n_ {
            it.x[i as usize] = 0.0;
            it.x[i as usize] = cpp_max(it.x[i as usize], self.model_.lb(i));
            it.x[i as usize] = cpp_min(it.x[i as usize], self.model_.ub(i));
        }

        let empty_scaling: Vec<f64> = Vec::new();
        let mut temp_m: Vec<f64> = vec![0.0; self.m_ as usize];

        // store b in y
        it.y = self.model_.b().to_vec();
        if norm2(&it.y) < 1e-6 {
            vectorAdd_scalar(&mut it.y, 1e-3);
        }

        if self.options_.nla == kHipoNormalEqString {
            // use y to store b-A*x
            self.model_
                .A()
                .alpha_product_plus_y(-1.0, &it.x, &mut it.y, false);

            // solve A*A^T * dx = b-A*x with factorisation and store the result in
            // temp_m

            // factorise A*A^T
            // PORT NOTE: `split_ls` reads `it_`, which is moved out here, so the
            // context borrows the local iterate's `data` instead.
            let ls = self
                .LS_
                .as_deref_mut()
                .expect("LS_ is created by initialise");
            let status = {
                let mut ctx = LsCtx {
                    kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                    options: &mut self.options_,
                    model: &self.model_,
                    regul: &self.regul_,
                    info: &mut self.info_,
                    data: &mut it.data,
                    logger: &self.logger_,
                };
                ls.factorNE(&mut ctx, &empty_scaling)
            };
            if status != 0 {
                self.logger_
                    .printe("Error while factorising normal equations\n");
                self.info_.status = status;
                self.it_ = Some(it);
                return true;
            }

            let ls = self
                .LS_
                .as_deref_mut()
                .expect("LS_ is created by initialise");
            let status = {
                let mut ctx = LsCtx {
                    kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                    options: &mut self.options_,
                    model: &self.model_,
                    regul: &self.regul_,
                    info: &mut self.info_,
                    data: &mut it.data,
                    logger: &self.logger_,
                };
                ls.solveNE(&mut ctx, &it.y, &mut temp_m)
            };
            if status != 0 {
                self.logger_
                    .printe("Error while solving normal equations\n");
                self.info_.status = status;
                self.it_ = Some(it);
                return true;
            }
        } else if self.options_.nla == kHipoAugmentedString {
            // obtain solution of A*A^T * dx = b-A*x by solving
            // [ -I  A^T] [...] = [ -x]
            // [  A   0 ] [ dx] = [ b ]

            let ls = self
                .LS_
                .as_deref_mut()
                .expect("LS_ is created by initialise");
            let status = {
                let mut ctx = LsCtx {
                    kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                    options: &mut self.options_,
                    model: &self.model_,
                    regul: &self.regul_,
                    info: &mut self.info_,
                    data: &mut it.data,
                    logger: &self.logger_,
                };
                ls.factorAS(&mut ctx, &empty_scaling)
            };
            if status != 0 {
                self.logger_
                    .printe("Error while factorising augmented system\n");
                self.info_.status = status;
                self.it_ = Some(it);
                return true;
            }

            let mut rhs_x: Vec<f64> = vec![0.0; self.n_ as usize];
            for i in 0..self.n_ {
                rhs_x[i as usize] = -it.x[i as usize];
            }
            let mut lhs_x: Vec<f64> = vec![0.0; self.n_ as usize];

            let ls = self
                .LS_
                .as_deref_mut()
                .expect("LS_ is created by initialise");
            let status = {
                let mut ctx = LsCtx {
                    kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                    options: &mut self.options_,
                    model: &self.model_,
                    regul: &self.regul_,
                    info: &mut self.info_,
                    data: &mut it.data,
                    logger: &self.logger_,
                };
                ls.solveAS(&mut ctx, &rhs_x, &it.y, &mut lhs_x, &mut temp_m)
            };
            if status != 0 {
                self.logger_
                    .printe("Error while solving augmented system\n");
                self.info_.status = status;
                self.it_ = Some(it);
                return true;
            }
        }

        // compute dx = A^T * (A*A^T)^{-1} * (b-A*x) and store the result in xl
        it.xl.resize(self.n_ as usize, 0.0);
        for i in 0..self.n_ {
            it.xl[i as usize] = 0.0;
        }
        self.model_
            .A()
            .alpha_product_plus_y(1.0, &temp_m, &mut it.xl, true);

        // x += dx;
        vectorAdd(&mut it.x, &it.xl, 1.0, 1.0);
        // *********************************************************************

        // *********************************************************************
        // xl, xu starting point
        // *********************************************************************
        // compute xl, xu that satisfy linear constraints
        {
            let mut violation: f64 = 0.0;
            for i in 0..self.n_ {
                if self.model_.hasLb(i) {
                    it.xl[i as usize] = it.x[i as usize] - self.model_.lb(i);
                    violation = cpp_min(violation, it.xl[i as usize]);
                } else {
                    it.xl[i as usize] = 0.0;
                }
                if self.model_.hasUb(i) {
                    it.xu[i as usize] = self.model_.ub(i) - it.x[i as usize];
                    violation = cpp_min(violation, it.xu[i as usize]);
                } else {
                    it.xu[i as usize] = 0.0;
                }
            }

            // shift to be positive
            violation = 1.0 + cpp_max(0.0, -1.5 * violation);
            vectorAdd_scalar(&mut it.xl, violation);
            vectorAdd_scalar(&mut it.xu, violation);
        }
        // *********************************************************************

        // *********************************************************************
        // y starting point
        // *********************************************************************

        let mut cPlusQx: Vec<f64> = self.model_.c().to_vec();
        if norm2(&cPlusQx) < 1e-6 {
            vectorAdd_scalar(&mut cPlusQx, 1e-3);
        }
        if self.model_.qp() {
            self.model_
                .Q()
                .alpha_product_plus_y(self.model_.sense(), &it.x, &mut cPlusQx);
        }

        let max_norm_x = cpp_max(infNorm(&it.x), cpp_max(infNorm(&it.xl), infNorm(&it.xu)));
        if infNorm(&cPlusQx) > max_norm_x * 1e3 {
            cPlusQx = self.model_.c().to_vec();
        }

        if self.options_.nla == kHipoNormalEqString {
            // compute A*(c+Qx)
            for i in 0..self.m_ {
                temp_m[i as usize] = 0.0;
            }
            self.model_
                .A()
                .alpha_product_plus_y(1.0, &cPlusQx, &mut temp_m, false);

            let ls = self
                .LS_
                .as_deref_mut()
                .expect("LS_ is created by initialise");
            let status = {
                let mut ctx = LsCtx {
                    kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                    options: &mut self.options_,
                    model: &self.model_,
                    regul: &self.regul_,
                    info: &mut self.info_,
                    data: &mut it.data,
                    logger: &self.logger_,
                };
                ls.solveNE(&mut ctx, &temp_m, &mut it.y)
            };
            if status != 0 {
                self.logger_
                    .printe("Error while solvingF normal equations\n");
                self.info_.status = status;
                self.it_ = Some(it);
                return true;
            }
        } else if self.options_.nla == kHipoAugmentedString {
            // obtain solution of A*A^T * y = A*(c+Qx) by solving
            // [ -I  A^T] [...] = [ c+Qx ]
            // [  A   0 ] [ y ] = [   0  ]

            let rhs_y: Vec<f64> = vec![0.0; self.m_ as usize];
            let mut lhs_x: Vec<f64> = vec![0.0; self.n_ as usize];

            let ls = self
                .LS_
                .as_deref_mut()
                .expect("LS_ is created by initialise");
            let status = {
                let mut ctx = LsCtx {
                    kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                    options: &mut self.options_,
                    model: &self.model_,
                    regul: &self.regul_,
                    info: &mut self.info_,
                    data: &mut it.data,
                    logger: &self.logger_,
                };
                ls.solveAS(&mut ctx, &cPlusQx, &rhs_y, &mut lhs_x, &mut it.y)
            };
            if status != 0 {
                self.logger_
                    .printe("Error while solving augmented system\n");
                self.info_.status = status;
                self.it_ = Some(it);
                return true;
            }
        }
        // *********************************************************************

        // *********************************************************************
        // zl, zu starting point
        // *********************************************************************
        // compute c + Q * x - A^T * y and store in zl
        it.zl = cPlusQx.clone();
        self.model_
            .A()
            .alpha_product_plus_y(-1.0, &it.y, &mut it.zl, true);

        // split result between zl and zu
        {
            let mut violation: f64 = 0.0;
            for i in 0..self.n_ {
                let val = it.zl[i as usize];
                it.zl[i as usize] = 0.0;
                it.zu[i as usize] = 0.0;

                if self.model_.hasLb(i) && self.model_.hasUb(i) {
                    it.zl[i as usize] = 0.5 * val;
                    it.zu[i as usize] = -0.5 * val;
                } else if self.model_.hasLb(i) {
                    it.zl[i as usize] = val;
                } else if self.model_.hasUb(i) {
                    it.zu[i as usize] = -val;
                }

                violation = cpp_min(violation, it.zl[i as usize]);
                violation = cpp_min(violation, it.zu[i as usize]);
            }

            // shift to be positive

            violation = 1.0 + cpp_max(0.0, -1.5 * violation);
            for i in 0..self.n_ {
                if self.model_.hasLb(i) {
                    it.zl[i as usize] += violation;
                }
                if self.model_.hasUb(i) {
                    it.zu[i as usize] += violation;
                }
            }
        }
        // *********************************************************************

        // *********************************************************************
        // improve centrality
        // *********************************************************************
        {
            let mut xsum: f64 = 1.0;
            let mut zsum: f64 = 1.0;
            let mut mu: f64 = 1.0;

            for i in 0..self.n_ {
                if self.model_.hasLb(i) {
                    xsum += it.xl[i as usize];
                    zsum += it.zl[i as usize];
                    mu += it.xl[i as usize] * it.zl[i as usize];
                }
                if self.model_.hasUb(i) {
                    xsum += it.xu[i as usize];
                    zsum += it.zu[i as usize];
                    mu += it.xu[i as usize] * it.zu[i as usize];
                }
            }

            let dx = 0.5 * mu / zsum;
            let dz = 0.5 * mu / xsum;

            vectorAdd_scalar(&mut it.xl, dx);
            vectorAdd_scalar(&mut it.xu, dx);
            for i in 0..self.n_ {
                if self.model_.hasLb(i) {
                    it.zl[i as usize] += dz;
                }
                if self.model_.hasUb(i) {
                    it.zu[i as usize] += dz;
                }
            }
        }
        // *********************************************************************

        self.it_ = Some(it);
        false
    }
}
