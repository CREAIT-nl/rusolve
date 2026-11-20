//! HiGHS source: ipm/hipo/ipm/Solver.cpp (lines 810-1146)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Centering parameter, multiple centrality correctors (Gondzio) and the
//! iterate / stagnation / termination checks.

use super::Solver;
use crate::ipm::hipo::Int;
use crate::ipm::hipo::auxiliary::vector_operations::vectorAdd;
use crate::ipm::hipo::ipm::iterate::NewtonDir;
use crate::ipm::hipo::ipm::parameters::{
    kDivergeTol, kGammaCorrector, kMccIncreaseAlpha, kMccIncreaseMin, kSigmaAffine,
};
use crate::ipm::hipo::ipm::status::{
    kStatusDualInfeasible, kStatusError, kStatusNoProgress, kStatusPDFeas, kStatusPrimalInfeasible,
};
use crate::ipm::hipo_wrapper::report::get_hipo_non_vertex_solution;
use crate::lp_data::hconst::HighsModelStatus;
use crate::lp_data::highs_info::HighsInfo;
use crate::lp_data::hstruct::{HighsBasis, HighsSolution};
use crate::lp_data::presolve_support::lp_kkt_check;
use crate::util::highs_utils::cpp_max;

impl Solver {
    /// HiGHS: Solver.cpp::Solver::sigmaAffine (lines 810-814)
    pub fn sigmaAffine(&mut self) {
        self.sigma_ = kSigmaAffine;

        self.it_.as_deref_mut().unwrap().data.back().sigma_aff = self.sigma_;
    }

    /// HiGHS: Solver.cpp::Solver::sigmaCorrectors (lines 816-830)
    pub fn sigmaCorrectors(&mut self) {
        if (self.alpha_primal_ > 0.5 && self.alpha_dual_ > 0.5) || self.iter_ == 1 {
            self.sigma_ = 0.01;
        } else if self.alpha_primal_ > 0.2 && self.alpha_dual_ > 0.2 {
            self.sigma_ = 0.1;
        } else if self.alpha_primal_ > 0.1 && self.alpha_dual_ > 0.1 {
            self.sigma_ = 0.25;
        } else if self.alpha_primal_ > 0.05 && self.alpha_dual_ > 0.05 {
            self.sigma_ = 0.5;
        } else {
            self.sigma_ = 0.9;
        }

        self.it_.as_deref_mut().unwrap().data.back().sigma = self.sigma_;
    }

    /// HiGHS: Solver.cpp::Solver::residualsMcc (lines 832-907)
    pub fn residualsMcc(&mut self) {
        // compute right-hand side for multiple centrality correctors

        // clear existing residuals
        self.it_.as_deref_mut().unwrap().clearRes(&self.model_);

        // stepsizes of current direction
        let mut alpha_p: f64 = 0.0;
        let mut alpha_d: f64 = 0.0;
        {
            let delta = &self.it_.as_deref().unwrap().delta;
            self.stepsToBoundary(&mut alpha_p, &mut alpha_d, delta, None, 1.0);
        }

        // compute increased stepsizes
        alpha_p = cpp_max(1.0, alpha_p + kMccIncreaseAlpha);
        alpha_d = cpp_max(1.0, alpha_d + kMccIncreaseAlpha);

        // compute trial point
        let mut xlt = self.it_.as_deref().unwrap().xl.clone();
        let mut xut = self.it_.as_deref().unwrap().xu.clone();
        let mut zlt = self.it_.as_deref().unwrap().zl.clone();
        let mut zut = self.it_.as_deref().unwrap().zu.clone();

        {
            let it = self.it_.as_deref().unwrap();
            vectorAdd(&mut xlt, &it.delta.xl, alpha_p, 1.0);
            vectorAdd(&mut xut, &it.delta.xu, alpha_p, 1.0);
            vectorAdd(&mut zlt, &it.delta.zl, alpha_d, 1.0);
            vectorAdd(&mut zut, &it.delta.zu, alpha_d, 1.0);
        }

        let sigma = self.sigma_;
        let n = self.n_;
        let it = self.it_.as_deref_mut().unwrap();
        let mu = it.mu;

        // compute right-hand side for mcc
        for i in 0..n {
            // res5
            if self.model_.hasLb(i) {
                let prod = xlt[i as usize] * zlt[i as usize];
                if prod < sigma * mu * kGammaCorrector {
                    // prod is small, we add something positive to res5

                    let temp = sigma * mu * kGammaCorrector - prod;
                    it.res.r5[i as usize] += temp;
                } else if prod > sigma * mu / kGammaCorrector {
                    // prod is large, we may subtract something large from res5.
                    // limit the amount to subtract to -sigma*mu/gamma

                    let mut temp = sigma * mu / kGammaCorrector - prod;
                    temp = cpp_max(temp, -sigma * mu / kGammaCorrector);
                    it.res.r5[i as usize] += temp;
                }
            } else {
                it.res.r5[i as usize] = 0.0;
            }

            // res6
            if self.model_.hasUb(i) {
                let prod = xut[i as usize] * zut[i as usize];
                if prod < sigma * mu * kGammaCorrector {
                    // prod is small, we add something positive to res6

                    let temp = sigma * mu * kGammaCorrector - prod;
                    it.res.r6[i as usize] += temp;
                } else if prod > sigma * mu / kGammaCorrector {
                    // prod is large, we may subtract something large from res6.
                    // limit the amount to subtract to -sigma*mu/gamma

                    let mut temp = sigma * mu / kGammaCorrector - prod;
                    temp = cpp_max(temp, -sigma * mu / kGammaCorrector);
                    it.res.r6[i as usize] += temp;
                }
            } else {
                it.res.r6[i as usize] = 0.0;
            }
        }
    }

    /// HiGHS: Solver.cpp::Solver::centralityCorrectors (lines 909-961)
    pub fn centralityCorrectors(&mut self) -> bool {
        // compute stepsizes of current direction
        let mut alpha_p_old: f64 = 0.0;
        let mut alpha_d_old: f64 = 0.0;
        {
            let delta = &self.it_.as_deref().unwrap().delta;
            self.stepsToBoundary(&mut alpha_p_old, &mut alpha_d_old, delta, None, 1.0);
        }

        let mut cor: Int = 0;
        while cor < self.info_.correctors {
            // compute rhs for corrector
            self.residualsMcc();

            // compute corrector
            let mut corr = NewtonDir::new(self.m_, self.n_);
            if self.solveNewtonSystem(&mut corr) {
                return true;
            }

            let mut alpha_p: f64 = 0.0;
            let mut alpha_d: f64 = 0.0;
            let mut wp = alpha_p_old * alpha_d_old;
            let mut wd = wp;
            {
                let delta = &self.it_.as_deref().unwrap().delta;
                self.bestWeight(delta, &corr, &mut wp, &mut wd, &mut alpha_p, &mut alpha_d);
            }

            if alpha_p < alpha_p_old + kMccIncreaseAlpha * kMccIncreaseMin
                && alpha_d < alpha_d_old + kMccIncreaseAlpha * kMccIncreaseMin
            {
                // reject corrector
                break;
            }

            if alpha_p >= alpha_p_old + kMccIncreaseAlpha * kMccIncreaseMin {
                // accept primal corrector
                {
                    let it = self.it_.as_deref_mut().unwrap();
                    vectorAdd(&mut it.delta.x, &corr.x, wp, 1.0);
                    vectorAdd(&mut it.delta.xl, &corr.xl, wp, 1.0);
                    vectorAdd(&mut it.delta.xu, &corr.xu, wp, 1.0);
                }
                alpha_p_old = alpha_p;
            }
            if alpha_d >= alpha_d_old + kMccIncreaseAlpha * kMccIncreaseMin {
                // accept dual corrector
                {
                    let it = self.it_.as_deref_mut().unwrap();
                    vectorAdd(&mut it.delta.y, &corr.y, wd, 1.0);
                    vectorAdd(&mut it.delta.zl, &corr.zl, wd, 1.0);
                    vectorAdd(&mut it.delta.zu, &corr.zu, wd, 1.0);
                }
                alpha_d_old = alpha_d;
            }

            if alpha_p > 0.95 && alpha_d > 0.95 {
                // stepsizes are large enough, stop
                cor += 1;
                break;
            }

            // else, keep computing correctors
            cor += 1;
        }

        self.it_.as_deref_mut().unwrap().data.back().correctors = cor;

        false
    }

    /// HiGHS: Solver.cpp::Solver::bestWeight (lines 963-996)
    pub fn bestWeight(
        &self,
        delta: &NewtonDir,
        corrector: &NewtonDir,
        wp: &mut f64,
        wd: &mut f64,
        alpha_p: &mut f64,
        alpha_d: &mut f64,
    ) {
        // Find the best primal and dual weights for the corrector in the interval
        // [alpha_p_old * alpha_d_old, 1].
        // Upon return, wp and wd are the optimal weights, alpha_p and alpha_d are the
        // corresponding stepsizes.

        // keep track of best stepsizes
        *alpha_p = 0.0;
        *alpha_d = 0.0;

        // initial weight
        let mut w = *wp;

        // divide interval into 9 points
        let step = (1.0 - w) / 8.0;

        // for each weight, compute stepsizes and save best ones
        while w <= 1.0 {
            let mut ap: f64 = 0.0;
            let mut ad: f64 = 0.0;
            self.stepsToBoundary(&mut ap, &mut ad, delta, Some(corrector), w);
            if ap > *alpha_p {
                *alpha_p = ap;
                *wp = w;
            }
            if ad > *alpha_d {
                *alpha_d = ad;
                *wd = w;
            }

            if step == 0.0 {
                break;
            }

            w += step;
        }
    }

    /// HiGHS: Solver.cpp::Solver::checkIterate (lines 998-1022)
    pub fn checkIterate(&mut self) -> bool {
        // Check that iterate is not NaN or Inf
        if self.it_.as_deref().unwrap().isNan() {
            self.logger_.printInfo("\nIterate is nan\n");
            self.info_.status = kStatusError;
            return true;
        } else if self.it_.as_deref().unwrap().isInf() {
            self.logger_.printInfo("\nIterate is inf\n");
            self.info_.status = kStatusError;
            return true;
        }

        // check that no component is negative
        for i in 0..self.n_ {
            if (self.model_.hasLb(i) && self.it_.as_deref().unwrap().xl[i as usize] < 0.0)
                || (self.model_.hasLb(i) && self.it_.as_deref().unwrap().zl[i as usize] < 0.0)
                || (self.model_.hasUb(i) && self.it_.as_deref().unwrap().xu[i as usize] < 0.0)
                || (self.model_.hasUb(i) && self.it_.as_deref().unwrap().zu[i as usize] < 0.0)
            {
                self.logger_.printInfo("\nIterate has negative component\n");
                return true;
            }
        }

        false
    }

    /// HiGHS: Solver.cpp::Solver::checkStagnation (lines 1024-1029)
    pub fn checkStagnation(&mut self) -> bool {
        let mut log_stream = String::new();
        let stagnation = self.it_.as_deref_mut().unwrap().stagnation(&mut log_stream);
        self.logger_.printInfo(&log_stream);
        stagnation
    }

    /// rusolve `ipm_hipo_stall_detect`: dual infeasibility flat while mu
    /// collapses. A healthy infeasible IPM cuts dinf by (1 - alpha_d) per step, so a 1e4 fall of
    /// mu in K steps with dinf above tolerance and not 10x lower is a stall (ns1688926, bdry2).
    pub fn checkDualStall(&mut self) -> bool {
        const K: usize = 5;
        let (mu, dinf) = {
            let it = self.it_.as_deref().unwrap();
            (it.mu, it.dinf)
        };
        self.stall_hist_.push((mu, dinf));
        let n = self.stall_hist_.len();
        if n <= K {
            return false;
        }
        let (mu0, dinf0) = self.stall_hist_[n - 1 - K];
        let stalled = mu > 0.0
            && mu0 > 1e4 * mu
            && dinf > self.options_.feasibility_tol
            && dinf > 0.1 * dinf0;
        if stalled {
            self.logger_.printInfo(&format!(
                "Dual stall: mu {:.1e} -> {:.1e}, dinf {:.1e} -> {:.1e} over {} iterations\n",
                mu0, mu, dinf0, dinf, K
            ));
        }
        stalled
    }

    /// HiGHS: Solver.cpp::Solver::checkBadIter (lines 1031-1077)
    pub fn checkBadIter(&mut self) -> bool {
        let mut terminate = false;
        let mut stagnation = if self.iter_ > 0 {
            self.checkStagnation()
        } else {
            false
        };
        if self.options_.stall_detect && self.checkDualStall() {
            stagnation = true;
        }

        // check for infeasibility
        let mu_is_large = if self.it_.as_deref().unwrap().best_mu > 0.0 {
            self.it_.as_deref().unwrap().mu > self.it_.as_deref().unwrap().best_mu * kDivergeTol
        } else {
            false
        };
        let pobj_is_very_large = self.it_.as_deref().unwrap().pobj
            < -cpp_max(
                self.it_.as_deref().unwrap().dobj.abs() * kDivergeTol * kDivergeTol,
                1.0,
            );
        let dobj_is_very_large = self.it_.as_deref().unwrap().dobj
            > cpp_max(
                self.it_.as_deref().unwrap().pobj.abs() * kDivergeTol * kDivergeTol,
                1.0,
            );
        let clearly_infeasible = self.iter_ > 5 && (pobj_is_very_large || dobj_is_very_large);

        if stagnation || mu_is_large || clearly_infeasible {
            let pobj_is_larger = self.it_.as_deref().unwrap().pobj
                < -cpp_max(self.it_.as_deref().unwrap().dobj.abs() * kDivergeTol, 1.0);
            let dobj_is_larger = self.it_.as_deref().unwrap().dobj
                > cpp_max(self.it_.as_deref().unwrap().pobj.abs() * kDivergeTol, 1.0);

            if pobj_is_larger {
                // problem is likely to be primal unbounded, i.e. dual infeasible
                self.logger_.print("=== Dual infeasible\n");
                self.info_.status = kStatusDualInfeasible;
                terminate = true;
            } else if dobj_is_larger {
                // problem is likely to be dual unbounded, i.e. primal infeasible
                self.logger_.print("=== Primal infeasible\n");
                self.info_.status = kStatusPrimalInfeasible;
                terminate = true;
            } else if stagnation {
                // stagnation detected, solution may still be good for highs kktCheck
                if self.info_.status != kStatusPDFeas && self.checkTerminationKkt() {
                    self.logger_
                        .printw("HiPO stagnated but HiGHS considers the solution acceptable\n");
                    self.logger_.print("=== Primal-dual feasible point found\n");
                    self.info_.status = kStatusPDFeas;
                } else {
                    self.info_.status = kStatusNoProgress;
                }
                terminate = true;
            }
        }

        terminate
    }

    /// HiGHS: Solver.cpp::Solver::checkTermination (lines 1079-1109)
    pub fn checkTermination(&mut self) -> bool {
        let feasible = self.it_.as_deref().unwrap().pinf < self.options_.feasibility_tol
            && self.it_.as_deref().unwrap().dinf < self.options_.feasibility_tol;
        let optimal = self.it_.as_deref().unwrap().pdgap < self.options_.optimality_tol;

        let mut terminate = false;

        if feasible && optimal {
            if self.crossoverIsOn() {
                if self.info_.status != kStatusPDFeas {
                    self.logger_.print("=== Primal-dual feasible point found\n");
                }
                self.info_.status = kStatusPDFeas;
                let ready_for_crossover = self
                    .it_
                    .as_deref()
                    .unwrap()
                    .infeasAfterDropping(&self.model_)
                    < self.options_.crossover_tol;
                if ready_for_crossover {
                    self.logger_.print("=== Ready for crossover\n");
                    terminate = true;
                }
            } else {
                terminate = if self.model_.qp() {
                    true
                } else {
                    self.checkTerminationKkt()
                };
                if terminate {
                    debug_assert!(self.info_.status != kStatusPDFeas);
                    self.logger_.print("=== Primal-dual feasible point found\n");
                    self.info_.status = kStatusPDFeas;
                }
            }
        }

        terminate
    }

    /// HiGHS: Solver.cpp::Solver::checkTerminationKkt (lines 1111-1146)
    pub fn checkTerminationKkt(&mut self) -> bool {
        if self.model_.qp() {
            // Not yet implemented for QP
            self.logger_.printInfo("kktCheck skipped for QP\n");
            return false;
        } else {
            // LP check
            self.logger_
                .printInfo("Solution may be optimal, perform kktCheck\n");
            let mut model_status = HighsModelStatus::kOptimal;
            let mut highs_info = HighsInfo::default();
            let mut highs_solution = HighsSolution::default();

            // Allow kktCheck to print only in debug mode (this is a copy of highs
            // options, not the original)
            self.Hoptions_.output_flag = self.logger_.debug(2);
            // C++ log_options.output_flag points at output_flag.
            self.Hoptions_.sync_log_options();

            if self.model_.lpOrig().is_none() {
                return false;
            }

            {
                let lp = self.model_.lpOrig().unwrap();
                let rhs = self.model_.rhsOrig();
                let constraints = self.model_.constraintsOrig();
                get_hipo_non_vertex_solution(
                    &self.Hoptions_,
                    lp,
                    self.model_.n_orig(),
                    self.model_.m_orig(),
                    rhs,
                    constraints,
                    self,
                    model_status,
                    &mut highs_solution,
                );
            }

            // HiGHS: HighsSolution.cpp::lpNoBasisKktCheck (lines 1035-1041),
            // inlined: an empty HighsBasis passed to lpKktCheck.
            {
                let basis = HighsBasis::default();
                lp_kkt_check(
                    &mut model_status,
                    &mut highs_info,
                    self.model_.lpOrig().unwrap(),
                    &highs_solution,
                    &basis,
                    &self.Hoptions_,
                    "During HiPO solve",
                );
            }

            if model_status == HighsModelStatus::kOptimal {
                self.logger_.printInfo("Check successfull\n");
                return true;
            } else {
                self.logger_.printInfo("Check failed\n");
            }
        }

        false
    }
}
