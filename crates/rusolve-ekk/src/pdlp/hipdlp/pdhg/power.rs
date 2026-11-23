//! HiGHS source: pdlp/hipdlp/pdhg.cc (lines 1529-1670, 1783-1901)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The power method for `||A||^2`, solver setup from `HighsOptions`, and
//! the scale / unscale / summary entry points. `powerMethodGpu` (1673-1780)
//! is not ported.

use super::super::defs::*;
use super::super::linalg;
use super::super::logger::LogLevel;
use super::PDLPSolver;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g_width;

impl<'a> PDLPSolver<'a> {
    /// HiGHS: pdhg.cc::PDLPSolver::powerMethod (lines 1529-1670)
    ///
    /// `power_method` is the constant `kCuPdlpAATPowerMethod` (2): only that
    /// branch is live. The `std::mt19937(12345)` / `uniform_real_distribution`
    /// draws into `x_vec` and the `normalize` that follows are overwritten by
    /// `x_vec.assign(num_row_, 1)` before any use, so they are not ported
    /// (no observable effect; `std::random_device rd` is unused). Loop of
    /// `max_iter = 20`: `y = A' x` (`a_matrix_.product_transpose`), `z = A
    /// y` (`product`), `z /= sqrt(dot(z, z))`, `y = A' z`, `lambda = dot(y,
    /// y)`, `x = z`. Returns `lambda` (or `1.0` when `num_col_ == 0 ||
    /// num_row_ == 0`). The `log_iters` output is print only.
    pub(crate) fn powerMethod(&self) -> f64 {
        // HiGHS: pdhg.cc::PDLPSolver::powerMethod (lines 1529-1670)
        let lp = &self.lp_;
        let mut op_norm_sq = 1.0;
        if lp.num_col_ == 0 || lp.num_row_ == 0 {
            return op_norm_sq;
        }

        // Parameters for the power method
        let max_iter: HighsInt = 20;
        let tol: f64 = 1e-6;

        // initialize a random vector x
        // PORT NOTE: std::mt19937(12345) draws and linalg::normalize (lines 1540-1546)
        // are overwritten by x_vec.assign(num_row_, 1) (line 1592) before any use; not ported.
        let mut x_vec = vec![0.0; lp.num_col_ as usize];

        let k_yanyu_power_method: HighsInt = 0;
        let k_yanyu_power_method_dev: HighsInt = 1;
        let k_cu_pdlp_aat_power_method: HighsInt = 2;

        let power_method: HighsInt =
            // kYanyuPowerMethod;
            // kYanyuPowerMethodDev;
            k_cu_pdlp_aat_power_method;
        crate::highs_log_dev!(
            &self.params_.log_options_,
            HighsLogType::kInfo,
            "Power method: {}\n",
            if power_method == k_yanyu_power_method {
                "Yanyu"
            } else if power_method == k_yanyu_power_method_dev {
                "Yanyu dev"
            } else {
                "CuPdlp-C"
            }
        );
        // Dev version of Yanyu power method (based on A'A) has
        //
        // * First iterate as vector of ones
        //
        // * No "convergence" check

        // cuPDLP-C version of power method is based on AA'
        //
        // Allocate memory for matrix-vector products
        let mut y_vec: Vec<f64>;
        let mut z_vec: Vec<f64>;
        if power_method < k_cu_pdlp_aat_power_method {
            y_vec = vec![0.0; lp.num_row_ as usize];
            z_vec = vec![0.0; lp.num_col_ as usize];
        } else {
            y_vec = vec![0.0; lp.num_col_ as usize];
            z_vec = vec![0.0; lp.num_row_ as usize];
        }

        let mut op_norm_sq_old: f64 = 0.0;
        let log_level: LogLevel = self.logger_.getLogLevel();
        let log_iters: HighsInt =
            (log_level == LogLevel::kVerbose || log_level == LogLevel::kDebug) as HighsInt;

        if log_iters != 0 {
            crate::highs_log_user!(
                &self.params_.log_options_,
                HighsLogType::kInfo,
                "It       lambda   dl_lambda\n"
            );
        }

        if power_method == k_yanyu_power_method_dev {
            // Start from a vector
            x_vec = vec![1.0; lp.num_col_ as usize];
        } else if power_method == k_cu_pdlp_aat_power_method {
            x_vec = vec![1.0; lp.num_row_ as usize];
        }
        let matrix = &lp.a_matrix_;
        let mut lambda: f64 = 0.0;
        let mut previous_lambda: f64 = lambda;
        let mut iter: HighsInt = 0;
        while iter < max_iter {
            if power_method == k_yanyu_power_method {
                // Original Yanyu power method
                //
                // Compute A'ax = A'(ax)
                linalg::ax(lp, &x_vec, &mut y_vec);
                linalg::aTy(lp, &y_vec, &mut z_vec);

                // Estimate the squared operator norm (largest eigenvalue of A'A)
                op_norm_sq = linalg::dot(&x_vec, &z_vec);

                // check for convergence
                if (op_norm_sq - op_norm_sq_old).abs() < tol * op_norm_sq {
                    return op_norm_sq;
                }
                let dl_op_norm_sq = (op_norm_sq - op_norm_sq_old).abs();
                op_norm_sq_old = op_norm_sq;

                // Prepare for the next iteration
                linalg::normalize(&mut z_vec);
                x_vec = z_vec.clone();
                if log_iters != 0 {
                    crate::highs_log_user!(
                        &self.params_.log_options_,
                        HighsLogType::kInfo,
                        "{:2} {} {}\n",
                        iter,
                        c_fmt_g_width(op_norm_sq, 6, 12),
                        c_fmt_g_width(dl_op_norm_sq, 4, 11)
                    );
                }
            } else {
                if power_method == k_yanyu_power_method_dev {
                    // Yanyu power method without "convergence" check
                    //
                    // Compute z = A'ax = A'(ax)
                    matrix.product(&mut y_vec, &x_vec);
                    matrix.product_transpose(&mut z_vec, &y_vec);
                    // Compute the Rayleigh quotient of x which, since x'x=1, is
                    //
                    // lambda = x'A'ax = x'(A'ax) = x'z
                    lambda = linalg::dot(&x_vec, &z_vec);
                    // x = z / norm(z)
                    let z_norm = linalg::dot(&z_vec, &z_vec).sqrt();
                    let mut i_col: HighsInt = 0;
                    while i_col < lp.num_col_ {
                        x_vec[i_col as usize] = z_vec[i_col as usize] / z_norm;
                        i_col += 1;
                    }
                } else {
                    // cuPDLP-C power method
                    //
                    // Compute z = AA'x = A(A'x)
                    matrix.product_transpose(&mut y_vec, &x_vec);
                    matrix.product(&mut z_vec, &y_vec);
                    // q = z / norm(z)
                    let z_norm = linalg::dot(&z_vec, &z_vec).sqrt();
                    let mut i_row: HighsInt = 0;
                    while i_row < lp.num_row_ {
                        z_vec[i_row as usize] /= z_norm;
                        i_row += 1;
                    }
                    // Compute the Rayleigh quotient of q which, since q'q=1, is
                    //
                    // lambda = q'AA'q = w^Tw = ||w||_2
                    //
                    // where w = A'q
                    //
                    // y_vec is no longer needed, so w is stored in it
                    matrix.product_transpose(&mut y_vec, &z_vec);
                    lambda = linalg::dot(&y_vec, &y_vec);
                    x_vec = z_vec.clone();
                }
                let dl_lambda = (lambda - previous_lambda).abs();
                previous_lambda = lambda;
                if log_iters != 0 {
                    crate::highs_log_user!(
                        &self.params_.log_options_,
                        HighsLogType::kInfo,
                        "{:2} {} {}\n",
                        iter,
                        c_fmt_g_width(lambda, 6, 12),
                        c_fmt_g_width(dl_lambda, 4, 11)
                    );
                }
            }
            iter += 1;
        }
        if power_method != k_yanyu_power_method {
            op_norm_sq = lambda;
        }
        // If the method did not converge within max_iter
        op_norm_sq
    }

    /// HiGHS: pdhg.cc::PDLPSolver::setup (lines 1783-1874)
    ///
    /// `logger_.initialise(options.log_dev_level, &options.log_options)`;
    /// `highs_timer_p_ = Some(timer)`; the `highsLogUser` banner and the
    /// `logger_` messages are print only. `params_.initialise()`, then
    /// `tolerance = pdlp_optimality_tolerance` (`kkt_tolerance` when it is not
    /// `kDefaultKktTolerance`), `max_iterations = pdlp_iteration_limit`,
    /// `device_type = CPU`, `time_limit = options.time_limit`; scaling flags
    /// from `pdlp_scaling_mode` bits unless `pdlp_features_off &
    /// kPdlpScalingOff`; `ruiz_iterations`; `restart_strategy =
    /// ADAPTIVE_RESTART`, `use_halpern_restart = true`;
    /// `step_size_strategy` = `PID` when `pdlp_step_size_strategy` is PID or
    /// anything other than Fixed, else `FIXED`; `params_.log_options_ =
    /// options.log_options`. `passParams` calls are not stored.
    pub fn setup(&mut self, options: &HighsOptions, timer: &'a HighsTimer) {
        // HiGHS: pdhg.cc::PDLPSolver::setup (lines 1783-1874)
        // PORT NOTE: the `&timer` argument of logger_.initialise is dropped
        // (Logger keeps no timer, logger.rs).
        self.logger_
            .initialise(options.log_dev_level, &options.log_options);
        self.logger_.printHeader();
        self.highs_timer_p_ = Some(timer);
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "Using HiPDLP first order PDLP solver on a {}\n",
            "CPU: performance may be disappointing!"
        );
        // #ifdef CUPDLP_GPU block (lines 1795-1810): not ported.
        // #if PDLP_PROFILE block (lines 1812-1815): not ported.
        // #if PDLP_DEBUG_LOG block (lines 1817-1819): not ported.

        self.params_.initialise();
        // params.eta = 0; Not set in parse_options_file
        // params.omega = 0; Not set in parse_options_file
        self.params_.tolerance = options.pdlp_optimality_tolerance;
        if options.kkt_tolerance != kDefaultKktTolerance {
            self.params_.tolerance = options.kkt_tolerance;
        }
        self.params_.max_iterations = options.pdlp_iteration_limit;
        self.params_.device_type = Device::CPU;
        self.params_.time_limit = options.time_limit;
        // rusolve (HiGHS 1.16 backport, default off): read by the PDHG solve
        self.params_.h16_hipdlp_max = options.h16_hipdlp_max;

        self.params_.use_ruiz_scaling = false;
        self.params_.use_pc_scaling = false;
        self.params_.use_l2_scaling = false;
        if (options.pdlp_features_off & kPdlpScalingOff) == 0 {
            // Use scaling: now see which
            self.params_.use_ruiz_scaling = (options.pdlp_scaling_mode & kPdlpScalingRuiz) != 0;
            self.params_.use_pc_scaling = (options.pdlp_scaling_mode & kPdlpScalingPC) != 0;
            self.params_.use_l2_scaling = (options.pdlp_scaling_mode & kPdlpScalingL2) != 0;
        }
        self.params_.ruiz_iterations = options.pdlp_ruiz_iterations;
        // params_.ruiz_norm = INFINITY; Not set in parse_options_file
        // params_.pc_alpha = 1.0; Not set in parse_options_file

        // HiPDLP solve loop is Halpern-only.
        self.params_.restart_strategy = RestartStrategy::ADAPTIVE_RESTART;
        self.params_.use_halpern_restart = true;
        if (options.pdlp_features_off & kPdlpRestartOff) != 0
            || options.pdlp_restart_strategy != kPdlpRestartStrategyHalpern
        {
            self.logger_.info(
                "HiPDLP uses Halpern restart only; ignoring pdlp_restart_strategy \
                 and restart-off feature flag.",
            );
        }
        // params_.fixed_restart_interval = 0; Not set in parse_options_file

        self.params_.step_size_strategy = StepSizeStrategy::FIXED;
        if options.pdlp_step_size_strategy == kPdlpStepSizeStrategyPid {
            self.params_.step_size_strategy = StepSizeStrategy::PID;
        } else if options.pdlp_step_size_strategy != kPdlpStepSizeStrategyFixed {
            self.logger_.info(
                "HiPDLP only supports pdlp_step_size_strategy in {0: fixed, \
                 3: PID}; using PID for unsupported values.",
            );
            self.params_.step_size_strategy = StepSizeStrategy::PID;
        }

        // params_.malitsky_pock_params.initialise(); Not set in parse_options_file
        // params_.adaptive_linesearch_params.initialise(); Not set in
        // parse_options_file
        // scaling_.passParams(&params_); (line 1868): not stored.
        // restart_scheme_.passParams(&params_); (line 1869): not stored.
        // Copy what's needed to use HiGHS logging
        self.params_.log_options_ = options.log_options.clone();
        // log the options
        self.logger_.printParams(&self.params_);
    }

    /// HiGHS: pdhg.cc::PDLPSolver::scaleProblem (lines 1876-1881)
    ///
    /// `scaling_.initialize(&lp_)` then `scaling_.scaleProblem(&mut lp_,
    /// &params_)`.
    pub fn scaleProblem(&mut self) {
        // HiGHS: pdhg.cc::PDLPSolver::scaleProblem (lines 1876-1881)
        // scaling_.passLp(&lp_); (line 1877): not stored.
        // scaling_.passParams(&params_); (line 1878): not stored.
        let lp = &self.lp_;
        self.scaling_.initialize(lp);
        let lp = &mut self.lp_;
        let params = &self.params_;
        self.scaling_.scaleProblem(lp, params);
    }

    /// HiGHS: pdhg.cc::PDLPSolver::unscaleSolution (lines 1883-1897)
    pub fn unscaleSolution(&mut self, x: &mut Vec<f64>, y: &mut Vec<f64>) {
        // HiGHS: pdhg.cc::PDLPSolver::unscaleSolution (lines 1883-1897)
        self.scaling_.unscaleSolution(x, y);

        self.x_current_ = x.clone();
        self.y_current_ = y.clone();

        let col_scale = self.scaling_.getColScaling();
        if !self.dSlackPos_.is_empty() && col_scale.len() == self.dSlackPos_.len() {
            let mut i: usize = 0;
            while i < self.dSlackPos_.len() {
                self.dSlackPos_[i] *= col_scale[i];
                self.dSlackNeg_[i] *= col_scale[i];
                i += 1;
            }
        }
    }

    /// HiGHS: pdhg.cc::PDLPSolver::logSummary (lines 1899-1901). Print only.
    pub fn logSummary(&self) {
        let _ = &self.logger_;
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn power_method_on_diagonal() {
        // A = diag(2, 1): ||A||^2 = 4.
        let mut s = PDLPSolver::default();
        s.lp_.num_col_ = 2;
        s.lp_.num_row_ = 2;
        s.lp_.a_matrix_.num_col_ = 2;
        s.lp_.a_matrix_.num_row_ = 2;
        s.lp_.a_matrix_.start_ = vec![0, 1, 2];
        s.lp_.a_matrix_.index_ = vec![0, 1];
        s.lp_.a_matrix_.value_ = vec![2.0, 1.0];
        let l = s.powerMethod();
        assert!((l - 4.0).abs() < 1e-9);
        assert_eq!(PDLPSolver::default().powerMethod(), 1.0);
    }

    #[test]
    fn setup_from_default_options() {
        let options = HighsOptions::default();
        let timer = HighsTimer::default();
        let mut s = PDLPSolver::default();
        s.setup(&options, &timer);
        assert_eq!(s.params_.tolerance, options.pdlp_optimality_tolerance);
        assert_eq!(s.params_.max_iterations, options.pdlp_iteration_limit);
        assert!(s.params_.use_ruiz_scaling && s.params_.use_pc_scaling);
        assert!(!s.params_.use_l2_scaling);
        assert!(s.params_.use_halpern_restart);
        // Default pdlp_step_size_strategy is Adaptive: not Fixed -> PID.
        assert_eq!(s.params_.step_size_strategy, StepSizeStrategy::PID);
        let _ = linalg::norm2(&[]);
    }

    #[test]
    fn setup_scaling_off_kkt_tolerance_fixed_step() {
        // pdhg.cc 1825-1826, 1834, 1860-1862.
        let mut options = HighsOptions::default();
        options.pdlp_features_off = kPdlpScalingOff;
        options.kkt_tolerance = 1e-5;
        options.pdlp_step_size_strategy = kPdlpStepSizeStrategyFixed;
        let timer = HighsTimer::default();
        let mut s = PDLPSolver::default();
        s.setup(&options, &timer);
        assert_eq!(s.params_.tolerance, 1e-5);
        assert!(!s.params_.use_ruiz_scaling);
        assert!(!s.params_.use_pc_scaling);
        assert!(!s.params_.use_l2_scaling);
        assert_eq!(s.params_.step_size_strategy, StepSizeStrategy::FIXED);
    }

    #[test]
    fn unscale_solution_scales_dual_slacks() {
        // pdhg.cc 1883-1897: unscaled problem leaves x / y, copies them to the
        // current iterates, and multiplies the slacks by col_scale (all 1).
        let mut s = PDLPSolver::default();
        s.lp_.num_col_ = 2;
        s.lp_.num_row_ = 1;
        s.lp_.col_cost_ = vec![1.0, 2.0];
        s.lp_.row_lower_ = vec![1.0];
        s.lp_.row_upper_ = vec![1.0];
        s.scaling_.initialize(&s.lp_);
        s.dSlackPos_ = vec![3.0, 4.0];
        s.dSlackNeg_ = vec![5.0, 6.0];
        let mut x = vec![1.5, 2.5];
        let mut y = vec![0.5];
        s.unscaleSolution(&mut x, &mut y);
        assert_eq!(x, vec![1.5, 2.5]);
        assert_eq!(s.x_current_, x);
        assert_eq!(s.y_current_, y);
        assert_eq!(s.dSlackPos_, vec![3.0, 4.0]);
        assert_eq!(s.dSlackNeg_, vec![5.0, 6.0]);
    }
}
