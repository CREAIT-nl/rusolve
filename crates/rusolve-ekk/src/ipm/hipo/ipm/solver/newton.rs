//! HiGHS source: ipm/hipo/ipm/Solver.cpp (lines 237-321, 340-489)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! IPX refinement hand-off, the Newton system (2x2 NE/AS reduction and the
//! 6x6 recovery; the driver `solveNewtonSystem` is in `ipm/refine.rs`)
//! and the step-to-boundary ratio tests. `prepareIpx` moves
//! the user callback from `control_` into `ipx_lps_` (C++:
//! `ipx_lps_.SetCallback(control_.callback())`); `refineWithIpx` moves it
//! back once IPX returns.

use super::Solver;
use crate::ipm::hipo::Int;
use crate::ipm::hipo::auxiliary::vector_operations::vectorScale;
use crate::ipm::hipo::ipm::iterate::{NewtonDir, Residuals};
use crate::ipm::hipo::ipm::linear_solver::LsCtx;
use crate::ipm::hipo::ipm::status::{IpxToHipoStatus, kStatusBasic};
use crate::ipm::ipx::info::status_string;
use crate::ipm::ipx::{IPX_STATUS_not_run, IPX_STATUS_optimal, Parameters};
use crate::ipm::ipx_wrapper::vec_data;
use crate::lp_data::hconst::{kHighsOffString, kHighsOnString};
use crate::lp_data::highs_options::kHipoNormalEqString;
use crate::util::highs_utils::cpp_min;

impl Solver {
    /// HiGHS: Solver.cpp::Solver::prepareIpx (lines 237-285)
    pub fn prepareIpx(&mut self, restart: bool) -> bool {
        // Return true if an error occurred;

        debug_assert!(!self.model_.qp());

        let mut ipx_param = Parameters::default();
        ipx_param.display = self.options_.display as Int;
        ipx_param.highs_logging = true;
        ipx_param.log_options = self.options_.log_options.clone();
        ipx_param.dualize = 0;
        // rusolve: a cold restart is the stand-alone IPX run, so
        // IPX picks its own formulation (dualized on neos-5251015, where the primal one is 40x slower)
        let cold = restart && self.Hoptions_.lp_ipm_failure_policy;
        if cold {
            ipx_param.dualize = -1;
        }
        // rusolve (HiGHS 1.16 backport, default off): inert flags from the HiPO options
        ipx_param.h16_ipx_stall = self.options_.h16_ipx_stall;
        ipx_param.h16_ipx_scaled = self.options_.h16_ipx_scaled;
        ipx_param.h16_ipm_gap_offset = self.options_.h16_ipm_gap_offset;
        ipx_param.stall_to_crossover = self.Hoptions_.lp_ipm_failure_policy;
        ipx_param.crossover_ipm_crash = self.Hoptions_.lp_crossover_ipm_crash;
        ipx_param.crossover_batch_push = self.Hoptions_.lp_crossover_batch_push;
        // rusolve lp_ipx_no_drop
        if self.Hoptions_.lp_ipx_no_drop {
            ipx_param.ipm_drop_primal = 0.0;
            ipx_param.ipm_drop_dual = 0.0;
        }
        ipx_param.ipm_drop_dual_soft =
            self.Hoptions_.lp_ipx_soft_drop && !self.Hoptions_.run_centring;

        if self.options_.crossover == kHighsOnString {
            ipx_param.run_crossover = 1;
        } else if self.options_.crossover == kHighsOffString {
            ipx_param.run_crossover = 0;
        } else {
            ipx_param.run_crossover = -1;
        }

        ipx_param.ipm_feasibility_tol = self.options_.feasibility_tol;
        ipx_param.ipm_optimality_tol = self.options_.optimality_tol;
        // rusolve lp_crossover_ipx_refine: IPX iterates on HiPO's
        // point down to a 100x tighter gap before crossover (Andersen and Andersen 2000; Bixby
        // and Saltzman 1994: crossover starts from a well-centred, nearly complementary point)
        if self.Hoptions_.lp_crossover_ipx_refine && !cold {
            ipx_param.ipm_optimality_tol = 1e-10;
        }
        if self.support_reduced_ {
            // the support-reduced point is only dual-feasible/optimal to the perturbation size
            ipx_param.ipm_optimality_tol = 1e-5;
            ipx_param.ipm_feasibility_tol = 1e-5;
        }
        ipx_param.start_crossover_tol = self.options_.crossover_tol;
        ipx_param.time_limit = self.options_.time_limit;
        ipx_param.timeless_log = self.options_.timeless_log;
        ipx_param.ipm_maxiter = self.options_.max_iter - self.iter_;
        self.ipx_lps_.set_parameters(ipx_param);

        let cb = self.control_.take_callback();
        self.ipx_lps_.set_callback(cb);

        let load_status = self.model_.loadIntoIpx(&mut self.ipx_lps_);

        if load_status != 0 {
            self.logger_.printInfo("Error loading model into IPX\n");
            return true;
        }

        // rusolve: the policy restarts IPX from its own
        // starting point, not from the HiPO iterate that just failed
        if cold {
            return false;
        }

        let mut x = Vec::new();
        let mut xl = Vec::new();
        let mut xu = Vec::new();
        let mut slack = Vec::new();
        let mut y = Vec::new();
        let mut zl = Vec::new();
        let mut zu = Vec::new();
        self.getInteriorSolution(
            &mut x, &mut xl, &mut xu, &mut slack, &mut y, &mut zl, &mut zu,
        );

        // libc++ `std::vector::data()` of an empty vector is `nullptr`, which
        // IPX rejects (`IPX_ERROR_argument_null`): pass `vec_data`.
        let start_point_status = self.ipx_lps_.load_ipm_starting_point(
            vec_data(&x),
            vec_data(&xl),
            vec_data(&xu),
            vec_data(&slack),
            vec_data(&y),
            vec_data(&zl),
            vec_data(&zu),
        );

        if start_point_status != 0 {
            self.logger_
                .printInfo("Error loading starting point into IPX\n");
            return true;
        }

        false
    }

    /// HiGHS: Solver.cpp::Solver::refineWithIpx (lines 287-320)
    pub fn refineWithIpx(&mut self) {
        if self.checkInterrupt() {
            return;
        }

        let restart = self.statusNeedsRefinement() && self.refinementIsOn();
        if restart {
            self.logger_.print("\nRestarting with IPX\n");
        } else if self.statusAllowsCrossover() && self.crossoverIsOn() {
            self.logger_.print("\nRunning crossover with IPX\n");
        } else {
            return;
        }

        if self.prepareIpx(restart) {
            // PORT NOTE: C++ `control_` keeps its pointer; move the callback
            // back from `ipx_lps_` so `solveHipo` can return it.
            self.control_.callback_ = self.ipx_lps_.take_callback();
            return;
        }

        self.ipx_lps_.solve();
        // PORT NOTE: move the callback back (see above); IPX is done with it.
        self.control_.callback_ = self.ipx_lps_.take_callback();
        self.info_.ipx_used = true;

        self.info_.ipx_info = self.ipx_lps_.get_info();

        // Convert between ipx and hipo status
        self.info_.status = IpxToHipoStatus(self.info_.ipx_info.status_ipm);

        let mut log_str = format!(
            "IPX reports: ipm {}",
            status_string(self.info_.ipx_info.status_ipm)
        );
        if self.info_.ipx_info.status_crossover != IPX_STATUS_not_run {
            log_str.push_str(&format!(
                ", crossover {}",
                status_string(self.info_.ipx_info.status_crossover)
            ));
        }
        log_str.push('\n');
        self.logger_.print(&log_str);

        if self.info_.ipx_info.status_crossover == IPX_STATUS_optimal {
            self.info_.status = kStatusBasic;
        }
    }

    /// HiGHS: Solver.cpp::Solver::solve2x2 (lines 340-402)
    pub fn solve2x2(&mut self, delta: &mut NewtonDir, rhs: &Residuals) -> bool {
        let theta_inv: Vec<f64> = self.it_.as_ref().unwrap().scaling.clone();

        let res7 = self.it_.as_ref().unwrap().residual7(&self.model_, rhs);

        // NORMAL EQUATIONS
        if self.options_.nla == kHipoNormalEqString {
            let res8 = self
                .it_
                .as_ref()
                .unwrap()
                .residual8(&self.model_, &self.regul_, rhs, &res7);

            // factorise normal equations, if not yet done
            if !self.LS_.as_ref().unwrap().valid() {
                let ls = self.LS_.as_deref_mut().expect("LS_");
                let status = {
                    let mut ctx = LsCtx {
                        kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                        options: &mut self.options_,
                        model: &self.model_,
                        regul: &self.regul_,
                        info: &mut self.info_,
                        data: &mut self.it_.as_mut().expect("it_").data,
                        logger: &self.logger_,
                    };
                    ls.factorNE(&mut ctx, &theta_inv)
                };
                if status != 0 {
                    self.logger_
                        .printe("Error while factorising normal equations\n");
                    self.info_.status = status;
                    return true;
                }
                let ls = self.LS_.as_deref_mut().expect("LS_");
                self.it_
                    .as_mut()
                    .unwrap()
                    .getReg(&self.model_, ls, &self.options_.nla);
            }

            // solve with normal equations
            let ls = self.LS_.as_deref_mut().expect("LS_");
            let status = {
                let mut ctx = LsCtx {
                    kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                    options: &mut self.options_,
                    model: &self.model_,
                    regul: &self.regul_,
                    info: &mut self.info_,
                    data: &mut self.it_.as_mut().expect("it_").data,
                    logger: &self.logger_,
                };
                ls.solveNE(&mut ctx, &res8, &mut delta.y)
            };
            if status != 0 {
                self.logger_
                    .printe("Error while solving normal equations\n");
                self.info_.status = status;
                return true;
            }

            // Compute delta.x
            // Deltax = A^T * Deltay - res7;
            delta.x = res7.clone();
            self.model_
                .A()
                .alpha_product_plus_y(-1.0, &delta.y, &mut delta.x, true);
            vectorScale(&mut delta.x, -1.0);

            // Deltax = (Theta^-1+Rp+Q)^-1 * Deltax
            for i in 0..self.n_ {
                let mut denom = theta_inv[i as usize] + self.regul_.primal;
                if self.model_.qp() {
                    denom += self.model_.sense() * self.model_.Q().diag(i);
                }
                delta.x[i as usize] /= denom;
            }
        }
        // AUGMENTED SYSTEM
        else {
            // factorise augmented system, if not yet done
            if !self.LS_.as_ref().unwrap().valid() {
                let ls = self.LS_.as_deref_mut().expect("LS_");
                let status = {
                    let mut ctx = LsCtx {
                        kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                        options: &mut self.options_,
                        model: &self.model_,
                        regul: &self.regul_,
                        info: &mut self.info_,
                        data: &mut self.it_.as_mut().expect("it_").data,
                        logger: &self.logger_,
                    };
                    ls.factorAS(&mut ctx, &theta_inv)
                };
                if status != 0 {
                    self.logger_
                        .printe("Error while factorising augmented system\n");
                    self.info_.status = status;
                    return true;
                }
                let ls = self.LS_.as_deref_mut().expect("LS_");
                self.it_
                    .as_mut()
                    .unwrap()
                    .getReg(&self.model_, ls, &self.options_.nla);
            }

            // solve with augmented system
            let ls = self.LS_.as_deref_mut().expect("LS_");
            let status = {
                let mut ctx = LsCtx {
                    kkt: self.kkt_.as_deref_mut().expect("kkt_"),
                    options: &mut self.options_,
                    model: &self.model_,
                    regul: &self.regul_,
                    info: &mut self.info_,
                    data: &mut self.it_.as_mut().expect("it_").data,
                    logger: &self.logger_,
                };
                ls.solveAS(&mut ctx, &res7, &rhs.r1, &mut delta.x, &mut delta.y)
            };
            if status != 0 {
                self.logger_
                    .printe("Error while solving augmented system\n");
                self.info_.status = status;
                return true;
            }
        }

        false
    }

    /// HiGHS: Solver.cpp::Solver::solve6x6 (lines 404-408)
    pub fn solve6x6(&mut self, delta: &mut NewtonDir, rhs: &Residuals) -> bool {
        if self.solve2x2(delta, rhs) {
            return true;
        }
        self.recoverDirection(delta, rhs);
        false
    }

    /// HiGHS: Solver.cpp::Solver::recoverDirection (lines 410-440)
    pub fn recoverDirection(&self, delta: &mut NewtonDir, rhs: &Residuals) {
        // Recover components xl, xu, zl, zu of partial direction delta.
        let xl = &self.it_.as_ref().unwrap().xl;
        let xu = &self.it_.as_ref().unwrap().xu;
        let zl = &self.it_.as_ref().unwrap().zl;
        let zu = &self.it_.as_ref().unwrap().zu;
        let res2 = &rhs.r2;
        let res3 = &rhs.r3;
        let res4 = &rhs.r4;
        let res5 = &rhs.r5;
        let res6 = &rhs.r6;

        for i in 0..self.n_ {
            if self.model_.hasLb(i) || self.model_.hasUb(i) {
                delta.xl[i as usize] = delta.x[i as usize] - res2[i as usize];
                delta.zl[i as usize] =
                    (res5[i as usize] - zl[i as usize] * delta.xl[i as usize]) / xl[i as usize];
            } else {
                delta.xl[i as usize] = 0.0;
                delta.zl[i as usize] = 0.0;
            }
        }
        for i in 0..self.n_ {
            if self.model_.hasLb(i) || self.model_.hasUb(i) {
                delta.xu[i as usize] = res3[i as usize] - delta.x[i as usize];
                delta.zu[i as usize] =
                    (res6[i as usize] - zu[i as usize] * delta.xu[i as usize]) / xu[i as usize];
            } else {
                delta.xu[i as usize] = 0.0;
                delta.zu[i as usize] = 0.0;
            }
        }
    }

    /// HiGHS: Solver.cpp::Solver::stepToBoundary (lines 442-467)
    /// C++ default argument `block = nullptr`; `cor` may be nullptr.
    pub fn stepToBoundary(
        &self,
        x: &[f64],
        dx: &[f64],
        cor: Option<&[f64]>,
        weight: f64,
        lo: bool,
        block: Option<&mut Int>,
    ) -> f64 {
        // Compute the largest alpha s.t. x + alpha * dx >= 0.
        // If cor is valid, consider x + alpha * (dx + w * cor) instead.
        // Use lo=1 for xl and zl, lo=0 for xu and zu.
        // Return the blocking index in block.

        let damp = 1.0 - 1e3 * f64::EPSILON;

        let mut alpha = 1.0;
        let mut bl: Int = -1;

        for i in 0..self.n_ {
            if (lo && self.model_.hasLb(i)) || (!lo && self.model_.hasUb(i)) {
                let c = if let Some(c) = cor {
                    c[i as usize] * weight
                } else {
                    0.0
                };
                if x[i as usize] + alpha * (dx[i as usize] + c) < 0.0 {
                    alpha = -(x[i as usize] * damp) / (dx[i as usize] + c);
                    bl = i;
                }
            }
        }
        if let Some(b) = block {
            *b = bl;
        }
        alpha
    }

    /// HiGHS: Solver.cpp::Solver::stepsToBoundary (lines 469-489)
    /// C++ default arguments `cor = nullptr`, `weight = 1.0`.
    pub fn stepsToBoundary(
        &self,
        alpha_primal: &mut f64,
        alpha_dual: &mut f64,
        delta: &NewtonDir,
        cor: Option<&NewtonDir>,
        weight: f64,
    ) {
        // compute primal and dual steps to boundary, given direction, corrector and
        // weight.

        let axl = self.stepToBoundary(
            &self.it_.as_ref().unwrap().xl,
            &delta.xl,
            cor.map(|c| c.xl.as_slice()),
            weight,
            true,
            None,
        );
        let axu = self.stepToBoundary(
            &self.it_.as_ref().unwrap().xu,
            &delta.xu,
            cor.map(|c| c.xu.as_slice()),
            weight,
            false,
            None,
        );
        let azl = self.stepToBoundary(
            &self.it_.as_ref().unwrap().zl,
            &delta.zl,
            cor.map(|c| c.zl.as_slice()),
            weight,
            true,
            None,
        );
        let azu = self.stepToBoundary(
            &self.it_.as_ref().unwrap().zu,
            &delta.zu,
            cor.map(|c| c.zu.as_slice()),
            weight,
            false,
            None,
        );

        *alpha_primal = cpp_min(axl, axu);
        *alpha_primal = cpp_min(*alpha_primal, 1.0);
        *alpha_dual = cpp_min(azl, azu);
        *alpha_dual = cpp_min(*alpha_dual, 1.0);
    }
}
