//! HiGHS source: ipm/hipo/ipm/Solver.h, ipm/hipo/ipm/Solver.cpp (lines 1-235)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiPO IPM driver. `LS_`, `kkt_`, `it_` are the C++ `unique_ptr`s
//! (`None` before `initialise` / `solve`). Borrowing conventions:
//!
//! * [`Solver::split_ls`] hands out the linear solver and its [`LsCtx`]
//!   (disjoint field borrows, including `it_->data`).
//! * A C++ argument that aliases a member (`solveNewtonSystem(it_->delta)`,
//!   `refine(it_->delta)`) is moved out with `std::mem::take` before the
//!   call and moved back after it; the arithmetic is unchanged.
//! * `start_time_` is uninitialised in C++ and unused in 1.15.1; zero.
//!
//! The out-of-line methods are split over `newton.rs` (237-489), `steps.rs`
//! (491-808), `correctors.rs` (810-1146) and `output.rs` (1148-1342); the
//! iterative refinement of `Refine.cpp` is `super::refine`.

pub mod correctors;
pub mod newton;
pub mod output;
pub mod steps;
pub mod support;

use crate::ipm::hipo::Int;
use crate::ipm::hipo::auxiliary::logger::Logger;
use crate::ipm::hipo::ipm::control::Control;
use crate::ipm::hipo::ipm::factor_highs_solver::FactorHiGHSSolver;
use crate::ipm::hipo::ipm::info::Info;
use crate::ipm::hipo::ipm::iterate::Iterate;
use crate::ipm::hipo::ipm::kkt_matrix::KktMatrix;
use crate::ipm::hipo::ipm::linear_solver::{LinearSolver, LsCtx};
use crate::ipm::hipo::ipm::model::Model;
use crate::ipm::hipo::ipm::options::Options;
use crate::ipm::hipo::ipm::parameters::Regularisation;
use crate::ipm::hipo::ipm::status::*;
use crate::ipm::hipo::oracle_trace;
use crate::ipm::hipo::oracle_trace::IterRecord;
use crate::ipm::ipx::lp_solver::LpSolver;
use crate::lp_data::hconst::kDefaultKktTolerance;
use crate::lp_data::highs_callback::HighsCallback;
use crate::lp_data::highs_lp::HighsLp;
use crate::lp_data::highs_options::HighsOptions;
use crate::model::highs_hessian::HighsHessian;
use crate::util::highs_timer::HighsTimer;
use crate::util::highs_utils::cpp_min;

/// HiGHS: Solver.h `class Solver`
#[derive(Default)]
pub struct Solver {
    pub model_: Model,
    pub LS_: Option<Box<dyn LinearSolver>>,
    pub kkt_: Option<Box<KktMatrix>>,
    pub it_: Option<Box<Iterate>>,
    pub m_: Int,
    pub n_: Int,
    pub iter_: Int,
    pub alpha_primal_: f64,
    pub alpha_dual_: f64,
    pub sigma_: f64,
    pub regul_: Regularisation,
    pub info_: Info,
    pub control_: Control,
    pub options_: Options,
    pub options_orig_: Options,
    pub Hoptions_: HighsOptions,
    pub ipx_lps_: LpSolver,
    pub logger_: Logger,
    pub start_time_: f64,
    /// rusolve: (mu, dinf) per iteration for the stall detector
    pub stall_hist_: Vec<(f64, f64)>,
    /// rusolve (dualize-guard): stop after the symbolic analysis and
    /// keep the predicted factor flops and nnz of the chosen Newton system
    pub analyse_only_: bool,
    pub analysed_flops_: f64,
    pub analysed_nz_: f64,
    /// dualize guard: the dual's factor was predicted costlier than the primal's; solve stopped
    pub guard_aborted_: bool,
    /// rusolve lp_ipm_support_reduction ran and kept its point
    pub support_reduced_: bool,
}

impl Solver {
    /// Rust-only: the linear solver and the C++ references it caches.
    /// Panics if `initialise` has not created `LS_`, `kkt_` and `it_`.
    pub fn split_ls(&mut self) -> (&mut dyn LinearSolver, LsCtx<'_>) {
        let ls = self
            .LS_
            .as_deref_mut()
            .expect("LS_ is created by initialise");
        let ctx = LsCtx {
            kkt: self
                .kkt_
                .as_deref_mut()
                .expect("kkt_ is created by initialise"),
            options: &mut self.options_,
            model: &self.model_,
            regul: &self.regul_,
            info: &mut self.info_,
            data: &mut self
                .it_
                .as_deref_mut()
                .expect("it_ is created by solve")
                .data,
            logger: &self.logger_,
        };
        (ls, ctx)
    }

    /// Rust-only: move the user callback back out of `control_`.
    pub fn take_callback(&mut self) -> Option<HighsCallback> {
        self.control_.take_callback()
    }

    /// HiGHS: Solver.h::Solver::setIpxTimerOffset (inline)
    pub fn setIpxTimerOffset(&mut self, offset: f64) {
        self.ipx_lps_.set_timer_offset(offset);
    }

    /// HiGHS: Solver.cpp::Solver::load (lines 14-20)
    pub fn load(&mut self, lp: &HighsLp, Q: &HighsHessian) -> Int {
        // rusolve (HiGHS 1.16 backport, default off): set before init, which runs preprocessing
        self.model_.h16_ipx_free_rows_ = self.options_.h16_ipx_free_rows;
        self.model_.h16_ipm_gap_offset_ = self.options_.h16_ipm_gap_offset;
        self.model_.h16_hipo_empty_rows_ = self.options_.h16_hipo_empty_rows;
        if self.model_.init(lp, Q) != 0 {
            self.logger_.printInfo("Error with model\n");
            return kStatusBadModel;
        }
        kStatusOk
    }

    /// HiGHS: Solver.cpp::Solver::setOptions (lines 22-59)
    /// HiGHS 1.16 (#3592): ipm/hipo/ipm/Solver.cpp::Solver::setOptions (lines 69-104; changed line 70)
    pub fn setOptions(&mut self, highs_options: &HighsOptions) {
        // HiGHS 1.16 #3592 (h16_hipo_display): ipm/hipo/ipm/Solver.cpp lines 69-70 (1.16.0)
        if highs_options.h16_hipo_display {
            self.options_.display = highs_options.output_flag;
        } else {
            self.options_.display = true;
            if !highs_options.output_flag | !highs_options.log_to_console {
                self.options_.display = false;
            }
        }

        self.options_.log_options = Some(highs_options.log_options.clone());

        // Debug option is already considered through log_options.log_dev_level in
        // hipo::Logger::debug

        self.options_.timeless_log = highs_options.timeless_log;
        self.options_.feasibility_tol = cpp_min(
            highs_options.primal_feasibility_tolerance,
            highs_options.dual_feasibility_tolerance,
        );
        self.options_.optimality_tol = highs_options.ipm_optimality_tolerance;
        self.options_.crossover_tol = highs_options.start_crossover_tolerance;

        if highs_options.kkt_tolerance != kDefaultKktTolerance {
            self.options_.feasibility_tol = highs_options.kkt_tolerance;
            self.options_.optimality_tol = 1e-1 * highs_options.kkt_tolerance;
            self.options_.crossover_tol = 1e-1 * highs_options.kkt_tolerance;
        }

        // hipo uses same timer as highs, so it is fine to pass the same time limit
        self.options_.time_limit = highs_options.time_limit;

        self.options_.max_iter = highs_options.ipm_iteration_limit;
        self.options_.crossover = highs_options.run_crossover.clone();
        self.options_.parallel = highs_options.parallel.clone();
        self.options_.parallel_type = highs_options.hipo_parallel_type.clone();
        self.options_.nla = highs_options.hipo_system.clone();
        self.options_.ordering = highs_options.hipo_ordering.clone();
        self.options_.block_size = highs_options.hipo_block_size;
        // rusolve: `lp_ipm_failure_policy` lets HiPO hand an
        // oversize factorisation to IPX at analyse time
        self.options_.flops_bailout = highs_options.lp_ipm_failure_policy;
        self.options_.stall_detect = highs_options.ipm_hipo_stall_detect;
        self.options_.amd_first = highs_options.ipm_hipo_amd_first;
        // rusolve: packed BLAS kernels, process-wide for the solve
        crate::ipm::hipo::blas::FAST_BLAS.store(
            highs_options.ipm_hipo_fast_blas,
            std::sync::atomic::Ordering::Relaxed,
        );
        crate::ipm::hipo::blas::CACHE_PERM.store(
            highs_options.ipm_hipo_cache_perm,
            std::sync::atomic::Ordering::Relaxed,
        );
        crate::ipm::hipo::blas::DENSE_COLS.store(
            highs_options.ipm_hipo_dense_cols || highs_options.ipm_hipo_dense_schur,
            std::sync::atomic::Ordering::Relaxed,
        );
        crate::ipm::hipo::blas::DENSE_SCHUR.store(
            highs_options.ipm_hipo_dense_schur,
            std::sync::atomic::Ordering::Relaxed,
        );
        crate::ipm::hipo::blas::PAR_NE.store(
            highs_options.ipm_hipo_par_ne,
            std::sync::atomic::Ordering::Relaxed,
        );
        // rusolve (HiGHS 1.16 backport, default off): inert flags, set from HighsOptions
        self.options_.h16_ipx_stall = highs_options.h16_ipx_stall;
        self.options_.h16_ipx_scaled = highs_options.h16_ipx_scaled;
        self.options_.h16_ipm_gap_offset = highs_options.h16_ipm_gap_offset;
        self.options_.h16_ipx_free_rows = highs_options.h16_ipx_free_rows;
        self.options_.h16_hipo_empty_rows = highs_options.h16_hipo_empty_rows;

        self.options_orig_ = self.options_.clone();
        self.Hoptions_ = highs_options.clone();
        self.resetOptions();
    }

    /// HiGHS: Solver.cpp::Solver::resetOptions (lines 61-65)
    pub fn resetOptions(&mut self) {
        self.options_ = self.options_orig_.clone();
        if self.options_.display {
            self.logger_.setOptions(self.options_.log_options.as_ref());
        }
        self.control_.setOptions(&self.options_);
    }

    /// HiGHS: Solver.cpp::Solver::setCallback (lines 66-68)
    /// The callback is moved in; `solveHipo` takes it back with
    /// [`Solver::take_callback`].
    pub fn setCallback(&mut self, callback: HighsCallback) {
        self.control_.setCallback(callback);
    }

    /// HiGHS: Solver.cpp::Solver::setTimer (lines 69-69)
    pub fn setTimer(&mut self, timer: &HighsTimer) {
        self.control_.setTimer(timer);
    }

    /// HiGHS: Solver.cpp::Solver::reset (lines 71-93)
    pub fn reset(&mut self) {
        self.m_ = self.model_.m();
        self.n_ = self.model_.n();

        self.info_ = Info::default();
        self.stall_hist_.clear();

        self.info_.ipx_used = false;
        self.info_.m_solver = self.m_;
        self.info_.n_solver = self.n_;
        self.info_.m_original = self.model_.m_orig();
        self.info_.n_original = self.model_.n_orig();

        self.resetOptions();

        self.iter_ = 0;
        self.alpha_primal_ = 0.0;
        self.alpha_dual_ = 0.0;
        self.sigma_ = 0.0;
        self.regul_ = Regularisation::default();

        self.LS_ = None;
        self.it_ = None;
    }

    /// HiGHS: Solver.cpp::Solver::solve (lines 95-118)
    pub fn solve(&mut self) {
        if !self.model_.ready() {
            self.info_.status = kStatusBadModel;
            return;
        }

        self.reset();

        // iterate object needs to be initialised before potentially interrupting
        self.it_ = Some(Box::new(Iterate::new(&self.model_)));
        if self.it_.is_none() {
            self.info_.status = kStatusError;
            return;
        }

        if self.checkInterrupt() {
            return;
        }

        self.printInfo();

        self.runIpm();
        if self.analyse_only_ || self.guard_aborted_ {
            return;
        }
        self.reduceSupport();
        self.refineWithIpx();
        let it = self.it_.as_ref().unwrap();
        let model = &self.model_;
        let info = &mut self.info_;
        it.finalResiduals(model, info);
        self.printSummary();
    }

    /// HiGHS: Solver.cpp::Solver::runIpm (lines 120-131)
    /// Oracle trace hook: one `oracle_trace::hipo_iter` record after each `makeStep()`
    /// (`tools/highs-oracle/FORMAT.md`, HiPO records).
    pub fn runIpm(&mut self) {
        if self.initialise() {
            return;
        }

        while self.iter_ < self.options_.max_iter {
            let g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(21);
            let stop = self.prepareIter();
            drop(g);
            if stop {
                break;
            }
            let g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(22);
            let stop = self.predictor();
            drop(g);
            if stop {
                break;
            }
            let g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(23);
            let stop = self.correctors();
            drop(g);
            if stop {
                break;
            }
            let g = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(24);
            self.makeStep();
            drop(g);

            // Oracle hook: hipo_iter right after makeStep()
            {
                let it = self.it_.as_ref().unwrap();
                let rec = IterRecord {
                    iter: self.iter_,
                    pobj: it.pobj,
                    dobj: it.dobj,
                    pinf: it.pinf,
                    dinf: it.dinf,
                    pdgap: it.pdgap,
                    mu: it.mu,
                    alpha_p: self.alpha_primal_,
                    alpha_d: self.alpha_dual_,
                    sigma: self.sigma_,
                    correctors: self.info_.correctors,
                };
                oracle_trace::hipo_iter(&rec, &it.x, &it.xl, &it.xu, &it.y, &it.zl, &it.zu);
            }
        }

        self.terminate();
    }

    /// HiGHS: Solver.cpp::Solver::initialise (lines 133-176)
    pub fn initialise(&mut self) -> bool {
        // Prepare ipm for execution.
        // Return true if an error occurred.

        self.start_time_ = self.control_.elapsed();
        // rusolve: the factorisation polls this deadline (it did not poll the interrupt)
        self.options_.deadline_ = if self.start_time_ >= 0.0
            && self.options_.time_limit >= 0.0
            && self.options_.time_limit.is_finite()
        {
            Some(
                std::time::Instant::now()
                    + std::time::Duration::from_secs_f64(
                        (self.options_.time_limit - self.start_time_).max(0.0),
                    ),
            )
        } else {
            None
        };

        self.kkt_ = Some(Box::new(KktMatrix::new()));
        if self.kkt_.is_none() {
            self.info_.status = kStatusError;
            return true;
        }

        // initialise linear solver
        self.LS_ = Some(Box::new(FactorHiGHSSolver::new(
            &self.options_,
            &self.logger_,
        )));
        if self.LS_.is_none() {
            self.info_.status = kStatusError;
            return true;
        }
        let analyse_only = self.analyse_only_;
        let mut analysed = (0.0, 0.0);
        let mut setup_status = 0;
        let guard_flops;
        {
            let (ls, mut ctx) = self.split_ls();
            let status = ls.setup(&mut ctx);
            guard_flops = ls.flops();
            if analyse_only {
                analysed = (ls.flops(), ls.nz());
                setup_status = status;
            } else if status != 0 {
                self.info_.status = status;
                return true;
            }
        }
        let pf = f64::from_bits(
            crate::ipm::hipo::blas::GUARD_PRIMAL_FLOPS.load(std::sync::atomic::Ordering::Relaxed),
        );
        if pf > 0.0 && !analyse_only {
            // dualize guard: below 1e9 flops the comparison is moot;
            // a dual that HiPO itself refuses as too expensive goes to IPX, keep it
            let o = &self.options_;
            let refused = o.flops_bailout
                && o.time_limit > 0.0
                && o.time_limit.is_finite()
                && 20.0 * guard_flops / 5.0e10 > o.time_limit;
            if !refused && guard_flops > 1.0e9 && guard_flops > 1.5 * pf {
                self.guard_aborted_ = true;
                self.info_.status = kStatusNoProgress;
                return true;
            }
        }
        if analyse_only {
            self.analysed_flops_ = analysed.0;
            self.analysed_nz_ = analysed.1;
            self.info_.status = setup_status;
            return true;
        }
        {
            let (ls, _ctx) = self.split_ls();
            ls.clear();
        }

        {
            let it = self.it_.as_mut().unwrap();
            it.data.append();
        }

        if self.checkInterrupt() {
            return true;
        }

        // decide number of correctors to use
        self.maxCorrectors();

        if self.startingPoint() {
            return true;
        }

        {
            let it = self.it_.as_mut().unwrap();
            let model = &self.model_;
            it.residual1234(model);
            it.computeMu(model);
            it.indicators(model);
        }

        self.printOutput();

        if self.checkInterrupt() {
            return true;
        }

        false
    }

    /// HiGHS: Solver.cpp::Solver::terminate (lines 178-181)
    pub fn terminate(&mut self) {
        self.info_.ipm_iter = self.iter_;
        if self.info_.status == kStatusNotRun {
            self.info_.status = kStatusMaxIter;
        }
    }

    /// HiGHS: Solver.cpp::Solver::prepareIter (lines 183-208)
    pub fn prepareIter(&mut self) -> bool {
        // Prepare next iteration.
        // Return true if Ipm main loop should be stopped

        if self.checkIterate() {
            return true;
        }
        if self.checkBadIter() {
            return true;
        }
        if self.checkTermination() {
            return true;
        }
        if self.checkInterrupt() {
            return true;
        }

        self.iter_ += 1;

        {
            let it = self.it_.as_mut().unwrap();
            let model = &mut self.model_;
            let logger = &self.logger_;
            model.adjustFreeVars(&mut it.x, &mut it.xl, &mut it.xu, logger);
        }

        // Clear Newton direction
        {
            let it = self.it_.as_mut().unwrap();
            it.clearDir();
        }

        // Clear any existing data in the linear solver
        {
            let (ls, _ctx) = self.split_ls();
            ls.clear();
        }

        {
            let it = self.it_.as_mut().unwrap();
            it.data.append();
        }

        // compute theta inverse
        {
            let it = self.it_.as_mut().unwrap();
            let model = &self.model_;
            it.computeScaling(model);
        }

        false
    }

    /// HiGHS: Solver.cpp::Solver::predictor (lines 210-223)
    pub fn predictor(&mut self) -> bool {
        // Compute affine scaling direction.
        // Return true if an error occurred.

        if self.checkInterrupt() {
            return true;
        }

        // compute sigma and residuals for affine scaling direction
        self.sigmaAffine();

        {
            let it = self.it_.as_mut().unwrap();
            let model = &self.model_;
            let sigma = self.sigma_;
            it.residual56(model, sigma);
        }

        {
            let mut delta = std::mem::take(&mut self.it_.as_mut().unwrap().delta);
            let result = self.solveNewtonSystem(&mut delta);
            self.it_.as_mut().unwrap().delta = delta;
            if result {
                return true;
            }
        }

        false
    }

    /// HiGHS: Solver.cpp::Solver::correctors (lines 225-235)
    pub fn correctors(&mut self) -> bool {
        // Compute multiple centrality correctors.
        // Return true if an error occurred.

        if self.checkInterrupt() {
            return true;
        }

        self.sigmaCorrectors();
        if self.centralityCorrectors() {
            return true;
        }

        false
    }
}
