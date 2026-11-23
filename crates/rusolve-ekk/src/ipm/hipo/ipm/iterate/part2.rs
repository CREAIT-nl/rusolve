//! HiGHS source: ipm/hipo/ipm/Iterate.cpp (lines 299-634)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complementary point, final residuals, regularisation extraction, the 6x6
//! residual for iterative refinement, the step and the stagnation test.
//! `std::stringstream& log_stream` is `&mut String`.

use super::{Iterate, NewtonDir};
use crate::ipm::hipo::Int;
use crate::ipm::hipo::auxiliary::logger::sci;
use crate::ipm::hipo::auxiliary::vector_operations::{infNorm, vectorAdd, vectorDivide};
use crate::ipm::hipo::ipm::info::Info;
use crate::ipm::hipo::ipm::linear_solver::LinearSolver;
use crate::ipm::hipo::ipm::model::Model;
use crate::ipm::hipo::ipm::parameters::Regularisation;
use crate::ipm::hipo::ipm::status::kStatusOk;
use crate::lp_data::highs_options::kHipoNormalEqString;
use crate::util::highs_utils::{cpp_max, cpp_min};

impl Iterate {
    /// HiGHS: Iterate.cpp::Iterate::dropToComplementarity (lines 299-378)
    pub fn dropToComplementarity(
        &self,
        model: &Model,
        x_cmp: &mut Vec<f64>,
        y_cmp: &mut Vec<f64>,
        z_cmp: &mut Vec<f64>,
    ) {
        let n = model.n();
        *x_cmp = vec![0.0; n as usize];
        *z_cmp = vec![0.0; n as usize];
        *y_cmp = self.y.clone();

        for j in 0..n {
            // value of x_[j] within bounds
            let mut xj = cpp_max(self.x[j as usize], model.lb(j));
            xj = cpp_min(xj, model.ub(j));

            // FIXED VARIABLE
            if model.lb(j) == model.ub(j) {
                x_cmp[j as usize] = model.lb(j);
                z_cmp[j as usize] = self.zl[j as usize] - self.zu[j as usize];
            }
            // BOTH BOUNDS FINITE
            else if model.hasLb(j) && model.hasUb(j) {
                if self.zl[j as usize] * self.xu[j as usize]
                    >= self.zu[j as usize] * self.xl[j as usize]
                {
                    // xlj/zlj <= xuj/zuj
                    // Primal lower is smaller than primal upper, wrt respective duals
                    if self.zl[j as usize] >= self.xl[j as usize] {
                        // drop x to lower bound, set z positive
                        x_cmp[j as usize] = model.lb(j);
                        z_cmp[j as usize] = cpp_max(0.0, self.zl[j as usize] - self.zu[j as usize]);
                    } else {
                        // drop z to zero, set x within bounds
                        x_cmp[j as usize] = xj;
                        z_cmp[j as usize] = 0.0;
                    }
                } else {
                    // xuj/zuj < xlj/zlj
                    // Primal upper is smaller than primal lower, wrt respective duals
                    if self.zu[j as usize] >= self.xu[j as usize] {
                        // drop x to upper bound, set z negative
                        x_cmp[j as usize] = model.ub(j);
                        z_cmp[j as usize] = cpp_min(0.0, self.zl[j as usize] - self.zu[j as usize]);
                    } else {
                        // drop z to zero, set x within bounds
                        x_cmp[j as usize] = xj;
                        z_cmp[j as usize] = 0.0;
                    }
                }
            }
            // LOWER BOUND FINITE
            else if model.hasLb(j) {
                if self.zl[j as usize] >= self.xl[j as usize] {
                    // drop x to lower bound, set z positive
                    x_cmp[j as usize] = model.lb(j);
                    z_cmp[j as usize] = cpp_max(0.0, self.zl[j as usize] - self.zu[j as usize]);
                } else {
                    // drop z to zero, set x within bounds
                    x_cmp[j as usize] = xj;
                    z_cmp[j as usize] = 0.0;
                }
            }
            // UPPER BOUND FINITE
            else if model.hasUb(j) {
                if self.zu[j as usize] >= self.xu[j as usize] {
                    // drop x to upper bound, set z negative
                    x_cmp[j as usize] = model.ub(j);
                    z_cmp[j as usize] = cpp_min(0.0, self.zl[j as usize] - self.zu[j as usize]);
                } else {
                    // drop z to zero, set x within bounds
                    x_cmp[j as usize] = xj;
                    z_cmp[j as usize] = 0.0;
                }
            }
            // NO BOUNDS
            else {
                x_cmp[j as usize] = xj;
                z_cmp[j as usize] = 0.0;
            }
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::infeasAfterDropping (lines 380-439)
    pub fn infeasAfterDropping(&self, model: &Model) -> f64 {
        // Compute estimate of residuals after dropping to complementarity (taken from
        // ipx).

        let mut pinf_max = 0.0;
        let mut dinf_max = 0.0;

        let n = model.n();
        for j in 0..n {
            let mut xdrop = 0.0;
            let mut zdrop = 0.0;

            if model.hasLb(j) && model.hasUb(j) {
                // BOTH BOUNDS FINITE
                if self.zl[j as usize] * self.xu[j as usize]
                    >= self.zu[j as usize] * self.xl[j as usize]
                {
                    if self.zl[j as usize] >= self.xl[j as usize] {
                        xdrop = self.x[j as usize] - model.lb(j);
                    } else {
                        zdrop = self.zl[j as usize] - self.zu[j as usize];
                    }
                } else {
                    if self.zu[j as usize] >= self.xu[j as usize] {
                        xdrop = self.x[j as usize] - model.ub(j);
                    } else {
                        zdrop = self.zl[j as usize] - self.zu[j as usize];
                    }
                }
            }
            // LOWER BOUND FINITE
            else if model.hasLb(j) {
                if self.zl[j as usize] >= self.xl[j as usize] {
                    xdrop = self.x[j as usize] - model.lb(j);
                } else {
                    // zu[j] is zero
                    zdrop = self.zl[j as usize];
                }
            }
            // UPPER BOUND FINITE
            else if model.hasUb(j) {
                if self.zu[j as usize] >= self.xu[j as usize] {
                    xdrop = self.x[j as usize] - model.ub(j);
                } else {
                    // zl[j] is zero
                    zdrop = -self.zu[j as usize];
                }
            }

            // largest entry in column j of A
            let mut amax = 0.0;
            let a = model.A();
            let start_j = a.start_[j as usize] as usize;
            let start_j1 = a.start_[(j + 1) as usize] as usize;
            for el in start_j..start_j1 {
                amax = cpp_max(amax, a.value_[el].abs());
            }

            pinf_max = cpp_max(pinf_max, xdrop.abs() * amax);
            dinf_max = cpp_max(dinf_max, zdrop.abs());
        }

        pinf_max /= 1.0 + model.normScaledRhs();
        dinf_max /= 1.0 + model.normScaledObj();

        cpp_max(pinf_max, dinf_max)
    }

    /// HiGHS: Iterate.cpp::Iterate::finalResiduals (lines 441-486)
    pub fn finalResiduals(&self, model: &Model, info: &mut Info) -> Int {
        // If ipx has been used, the information is already available, otherwise,
        // compute it.

        if !info.ipx_used {
            info.p_obj = self.pobj;
            info.d_obj = self.dobj;
            info.pd_gap = self.pdgap;

            let mut pinf = info.p_res_abs;
            let m = model.m();
            for i in 0..m {
                let mut val = self.res.r1[i as usize].abs();
                if model.scaled() {
                    val /= model.rowScale(i);
                }
                pinf = cpp_max(pinf, val);
            }
            let n = model.n();
            for i in 0..n {
                let mut val = self.res.r2[i as usize].abs();
                if model.scaled() {
                    val *= model.colScale(i);
                }
                pinf = cpp_max(pinf, val);

                val = self.res.r3[i as usize].abs();
                if model.scaled() {
                    val *= model.colScale(i);
                }
                pinf = cpp_max(pinf, val);
            }
            info.p_res_abs = pinf;
            info.p_res_rel = pinf / (1.0 + model.normUnscaledRhs());

            let mut dinf = info.d_res_abs;
            for i in 0..n {
                let mut val = self.res.r4[i as usize].abs();
                if model.scaled() {
                    val /= model.colScale(i);
                }
                dinf = cpp_max(dinf, val);
            }
            info.d_res_abs = dinf;
            info.d_res_rel = dinf / (1.0 + model.normUnscaledObj());
        } else {
            info.p_res_abs = info.ipx_info.abs_presidual;
            info.p_res_rel = info.ipx_info.rel_presidual;
            info.d_res_abs = info.ipx_info.abs_dresidual;
            info.d_res_rel = info.ipx_info.rel_dresidual;
            info.p_obj = info.ipx_info.pobjval;
            info.d_obj = info.ipx_info.dobjval;
            info.pd_gap = info.ipx_info.rel_objgap.abs();
        }

        kStatusOk
    }

    /// HiGHS: Iterate.cpp::Iterate::getReg (lines 488-500)
    /// `Rp = None, Rd = Some(0)` for NE; `Rp = Some(0)`, `Rd = Some(n)` if
    /// `m > 0` (else `None`) for AS.
    pub fn getReg(&mut self, model: &Model, LS: &mut dyn LinearSolver, nla: &str) {
        // extract regularisation
        LS.getReg(&mut self.total_reg);

        // easy access to primal/dual regularisation
        if nla == kHipoNormalEqString {
            self.Rp = None;
            self.Rd = Some(0);
        } else {
            self.Rp = Some(0);
            if model.m() > 0 {
                self.Rd = Some(model.n() as usize);
            } else {
                self.Rd = None;
            }
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::residuals6x6 (lines 502-564)
    pub fn residuals6x6(&mut self, model: &Model, regul: &Regularisation, d: &NewtonDir) {
        let dx = &d.x;
        let dy = &d.y;
        let dxl = &d.xl;
        let dxu = &d.xu;
        let dzl = &d.zl;
        let dzu = &d.zu;
        let m = model.m();
        let n = model.n();
        debug_assert!(self.Rd.is_some() || m == 0);

        // res1,2,3,4,5,6 contain the rhs of the linear system

        // ires1 = res1 - A * dx - Rd * dy
        self.ires.r1 = self.res.r1.clone();
        model
            .A()
            .alpha_product_plus_y(-1.0, dx, &mut self.ires.r1, false);
        let rd = self.Rd.unwrap_or(0);
        for i in 0..m {
            self.ires.r1[i as usize] -= self.total_reg[rd + i as usize] * dy[i as usize];
        }

        // ires2 = res2 - dx + dxl
        for i in 0..n {
            if model.hasLb(i) {
                self.ires.r2[i as usize] =
                    self.res.r2[i as usize] - dx[i as usize] + dxl[i as usize];
            } else {
                self.ires.r2[i as usize] = 0.0;
            }
        }

        // ires3 = res3 - dx - dxu
        for i in 0..n {
            if model.hasUb(i) {
                self.ires.r3[i as usize] =
                    self.res.r3[i as usize] - dx[i as usize] - dxu[i as usize];
            } else {
                self.ires.r3[i as usize] = 0.0;
            }
        }

        // ires4 = res4 - A^T * dy - dzl + dzu + Q * dx + Rp * dx
        self.ires.r4 = self.res.r4.clone();
        for i in 0..n {
            if model.hasLb(i) {
                self.ires.r4[i as usize] -= dzl[i as usize];
            }
            if model.hasUb(i) {
                self.ires.r4[i as usize] += dzu[i as usize];
            }
        }
        model
            .A()
            .alpha_product_plus_y(-1.0, dy, &mut self.ires.r4, true);
        for i in 0..n {
            let reg_p = match self.Rp {
                Some(rp) => self.total_reg[rp + i as usize],
                None => regul.primal,
            };
            self.ires.r4[i as usize] += reg_p * dx[i as usize];
        }
        if model.qp() {
            model
                .Q()
                .alpha_product_plus_y(model.sense(), dx, &mut self.ires.r4);
        }

        // ires5 = res5 - zl * dxl - xl * dzl
        for i in 0..n {
            if model.hasLb(i) {
                self.ires.r5[i as usize] = self.res.r5[i as usize]
                    - self.zl[i as usize] * dxl[i as usize]
                    - self.xl[i as usize] * dzl[i as usize];
            } else {
                self.ires.r5[i as usize] = 0.0;
            }
        }

        // ires6 = res6 - zu * dxu - xu * dzu
        for i in 0..n {
            if model.hasUb(i) {
                self.ires.r6[i as usize] = self.res.r6[i as usize]
                    - self.zu[i as usize] * dxu[i as usize]
                    - self.xu[i as usize] * dzu[i as usize];
            } else {
                self.ires.r6[i as usize] = 0.0;
            }
        }
    }

    /// HiGHS: Iterate.cpp::Iterate::assertConsistency (lines 566-573)
    pub fn assertConsistency(&self, n: Int, m: Int) {
        debug_assert!(self.x.len() as Int == n);
        debug_assert!(self.xl.len() as Int == n);
        debug_assert!(self.xu.len() as Int == n);
        debug_assert!(self.y.len() as Int == m);
        debug_assert!(self.zl.len() as Int == n);
        debug_assert!(self.zu.len() as Int == n);
    }

    /// HiGHS: Iterate.cpp::Iterate::makeStep (lines 575-587)
    pub fn makeStep(&mut self, alpha_primal: f64, alpha_dual: f64) {
        if cpp_min(alpha_primal, alpha_dual) < 0.05 {
            self.bad_iter_ += 1;
        } else {
            self.bad_iter_ = 0;
        }

        vectorAdd(&mut self.x, &self.delta.x, alpha_primal, 1.0);
        vectorAdd(&mut self.xl, &self.delta.xl, alpha_primal, 1.0);
        vectorAdd(&mut self.xu, &self.delta.xu, alpha_primal, 1.0);
        vectorAdd(&mut self.y, &self.delta.y, alpha_dual, 1.0);
        vectorAdd(&mut self.zl, &self.delta.zl, alpha_dual, 1.0);
        vectorAdd(&mut self.zu, &self.delta.zu, alpha_dual, 1.0);
    }

    /// HiGHS: Iterate.cpp::Iterate::stagnation (lines 589-634)
    pub fn stagnation(&mut self, log_stream: &mut String) -> bool {
        // too many iterations in a row with small stepsize
        let mut stagnation = self.bad_iter_ >= crate::ipm::hipo::ipm::parameters::kMaxBadIter;

        // the next tests are aimed at problems where things are going
        // catastrophically bad, which may continue iterating forever otherwise.

        // dx / x or dy / y becoming way too small
        let thresh_dxy_xy: f64 = 1e-30;

        let mut temp = self.delta.x.clone();
        vectorDivide(&mut temp, &self.x);
        let dx_x_max = infNorm(&temp);

        temp = self.delta.y.clone();
        vectorDivide(&mut temp, &self.y);
        let dy_y_max = infNorm(&temp);

        self.largest_dx_x_ = cpp_max(self.largest_dx_x_, dx_x_max);
        self.largest_dy_y_ = cpp_max(self.largest_dy_y_, dy_y_max);

        if dx_x_max < self.largest_dx_x_ * thresh_dxy_xy
            || dy_y_max < self.largest_dy_y_ * thresh_dxy_xy
        {
            stagnation = true;
            log_stream.push_str(&format!(
                "Bad direction ratios, dx_x {} (max {}), dy_y {} (max {})\n",
                sci(dx_x_max, 0, 1),
                sci(self.largest_dx_x_, 0, 1),
                sci(dy_y_max, 0, 1),
                sci(self.largest_dy_y_, 0, 1)
            ));
        }

        // infeasibilities jumping back up
        let thresh_inf_to_best: f64 = 1e12;

        if self.pinf > 0.0 {
            self.best_pinf_ = cpp_min(self.best_pinf_, self.pinf);
        }
        if self.dinf > 0.0 {
            self.best_dinf_ = cpp_min(self.best_dinf_, self.dinf);
        }

        if self.pinf > thresh_inf_to_best * self.best_pinf_
            || self.dinf > thresh_inf_to_best * self.best_dinf_
        {
            stagnation = true;
            log_stream.push_str(&format!(
                "Bad infeasibility, pinf {} (best {}), dinf {} (best {})\n",
                sci(self.pinf, 0, 1),
                sci(self.best_pinf_, 0, 1),
                sci(self.dinf, 0, 1),
                sci(self.best_dinf_, 0, 1)
            ));
        }

        stagnation
    }
}
