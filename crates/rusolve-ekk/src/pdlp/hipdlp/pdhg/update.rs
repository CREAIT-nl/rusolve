//! HiGHS source: pdlp/hipdlp/pdhg.cc (lines 1944-2359)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Step sizes, the PID primal-weight update, the plain PDHG updates and the
//! stability helpers. `updateX` / `updateY` / `updateIterates*` /
//! `checkNumericalStability` / `computeMovement` / `computeNonlinearity`
//! are reached only through `performPdhgStep`, which the CPU `solve` never
//! calls; they are ported for completeness. `hipdlpTimer*` calls are
//! dropped.

use super::super::defs::*;
use super::super::kDivergentMovement;
use super::super::linalg;
use super::super::solver_results::SolverResults;
use super::PDLPSolver;
use crate::prelude::*;
use crate::util::c_printf::{c_fmt_f, c_fmt_g};

impl<'a> PDLPSolver<'a> {
    /// HiGHS: pdhg.cc::PDLPSolver::initializeStepSizes (lines 1944-1977)
    ///
    /// `omega = (unscaled_c_norm_ + 1) / (unscaled_rhs_norm_ + 1)`;
    /// strategies other than FIXED / PID fall back to FIXED;
    /// `power_method_lambda = powerMethod()`; `base_step = 0.998 /
    /// sqrt(lambda)`; `eta = base_step`, `primal_step = base_step / omega`,
    /// `dual_step = base_step * omega`.
    pub(crate) fn initializeStepSizes(&mut self) {
        self.params_.omega = (self.unscaled_c_norm_ + 1.0) / (self.unscaled_rhs_norm_ + 1.0);
        self.primal_weight_ = self.params_.omega;
        self.best_primal_weight_ = self.primal_weight_;
        self.stepsize_.beta = self.primal_weight_ * self.primal_weight_;

        if self.params_.step_size_strategy != StepSizeStrategy::FIXED
            && self.params_.step_size_strategy != StepSizeStrategy::PID
        {
            self.logger_
                .info("HiPDLP only supports fixed/PID step size; falling back to fixed.");
            self.params_.step_size_strategy = StepSizeStrategy::FIXED;
        }

        // #ifdef CUPDLP_GPU branch dropped
        let op_norm_sq = self.powerMethod();

        self.stepsize_.power_method_lambda = op_norm_sq;

        let safety_factor: f64 = 0.998;
        let base_step = safety_factor / op_norm_sq.sqrt();

        self.params_.eta = base_step;
        self.stepsize_.primal_step = base_step / self.params_.omega;
        self.stepsize_.dual_step = base_step * self.params_.omega;

        highs_log_dev(
            &self.params_.log_options_,
            HighsLogType::kInfo,
            &format!(
                "Initial step sizes from power method lambda = {}: primal step= {}; dual step = {}, eta = {}, omega = {}\n",
                c_fmt_g(op_norm_sq, 6),
                c_fmt_g(self.stepsize_.primal_step, 6),
                c_fmt_g(self.stepsize_.dual_step, 6),
                c_fmt_g(self.params_.eta, 6),
                c_fmt_g(self.params_.omega, 6)
            ),
        );
    }

    /// HiGHS: pdhg.cc::PDLPSolver::updatePrimalWeightAtRestart (lines 1979-2050)
    ///
    /// Called as `updatePrimalWeightAtRestart(results_)`: the caller passes
    /// a clone of `results_`. `std::log` is `ln`, `std::log10`
    /// is `log10`, `std::exp` is `exp`. Ends with
    /// `restart_scheme_.update_beta(stepsize_.beta, &mut params_)`.
    pub(crate) fn updatePrimalWeightAtRestart(&mut self, results: &SolverResults) {
        // Compute delta = pdhg_iterate - anchor
        let mut primal_dist: f64 = 0.0;
        let mut dual_dist: f64 = 0.0;

        // #ifdef CUPDLP_GPU branch dropped
        // CPU version
        for i in 0..self.lp_.num_col_ {
            let diff = self.x_next_[i as usize] - self.x_anchor_[i as usize];
            primal_dist += diff * diff;
        }
        for j in 0..self.lp_.num_row_ {
            let diff = self.y_next_[j as usize] - self.y_anchor_[j as usize];
            dual_dist += diff * diff;
        }
        primal_dist = primal_dist.sqrt();
        dual_dist = dual_dist.sqrt();

        let rel_primal = results.primal_feasibility / (1.0 + self.unscaled_rhs_norm_);
        let rel_dual = results.dual_feasibility / (1.0 + self.unscaled_c_norm_);
        let ratio_infeas = if rel_primal > 0.0 {
            rel_dual / rel_primal
        } else {
            1e300
        };

        let silent = true;
        // cuPDLPx-style weight update (PID control)
        if primal_dist > 1e-16
            && dual_dist > 1e-16
            && primal_dist < 1e12
            && dual_dist < 1e12
            && ratio_infeas > 1e-8
            && ratio_infeas < 1e8
        {
            let error = dual_dist.ln() - primal_dist.ln() - self.primal_weight_.ln();

            self.primal_weight_error_sum_ =
                self.params_.i_smooth * self.primal_weight_error_sum_ + error;
            let delta_error = error - self.primal_weight_last_error_;

            self.primal_weight_ *= (self.params_.k_p * error
                + self.params_.k_i * self.primal_weight_error_sum_
                + self.params_.k_d * delta_error)
                .exp();

            self.primal_weight_last_error_ = error;
            if !silent {
                self.logger_.info(&format!(
                    "Primal weight updated: {}",
                    c_fmt_f(self.primal_weight_, 6)
                ));
            }
        } else {
            // Revert to best known weight
            self.primal_weight_ = self.best_primal_weight_;
            self.primal_weight_error_sum_ = 0.0;
            self.primal_weight_last_error_ = 0.0;
            if !silent {
                self.logger_.info(&format!(
                    "Weight update failed (bad norms/ratio), reverted to best: {}",
                    c_fmt_f(self.primal_weight_, 6)
                ));
            }
        }

        // Track best weight
        let gap = if rel_primal > 0.0 && rel_dual > 0.0 {
            (rel_dual / rel_primal).log10().abs()
        } else {
            self.best_primal_dual_residual_gap_
        };
        if gap < self.best_primal_dual_residual_gap_ {
            self.best_primal_dual_residual_gap_ = gap;
            self.best_primal_weight_ = self.primal_weight_;
        }

        let eta = (self.stepsize_.primal_step * self.stepsize_.dual_step).sqrt(); // invariant
        self.stepsize_.beta = self.primal_weight_ * self.primal_weight_;
        self.stepsize_.primal_step = eta / self.primal_weight_;
        self.stepsize_.dual_step = eta * self.primal_weight_;
        self.params_.omega = self.primal_weight_;
        self.restart_scheme_
            .update_beta(self.stepsize_.beta, &mut self.params_);
    }

    /// HiGHS: pdhg.cc::PDLPSolver::updateX (lines 2052-2062)
    pub(crate) fn updateX(&self, x: &[f64], aty: &[f64], primal_step: f64) -> Vec<f64> {
        let mut x_new = vec![0.0; self.lp_.num_col_ as usize];
        for i in 0..self.lp_.num_col_ {
            let gradient = self.lp_.col_cost_[i as usize] - aty[i as usize];
            x_new[i as usize] = linalg::projectBox(
                x[i as usize] - primal_step * gradient,
                self.lp_.col_lower_[i as usize],
                self.lp_.col_upper_[i as usize],
            );
        }
        x_new
    }

    /// HiGHS: pdhg.cc::PDLPSolver::updateY (lines 2064-2078)
    pub(crate) fn updateY(
        &self,
        y: &[f64],
        ax: &[f64],
        ax_next: &[f64],
        dual_step: f64,
    ) -> Vec<f64> {
        let mut y_new = vec![0.0; self.lp_.num_row_ as usize];
        for j in 0..self.lp_.num_row_ {
            let extr_ax = 2.0 * ax_next[j as usize] - ax[j as usize];
            let is_equality = self.lp_.row_lower_[j as usize] == self.lp_.row_upper_[j as usize];
            let q = self.lp_.row_lower_[j as usize];
            let dual_update = y[j as usize] + dual_step * (q - extr_ax);
            y_new[j as usize] = if is_equality {
                dual_update
            } else {
                linalg::projectNonNegative(dual_update)
            };
        }
        y_new
    }

    /// HiGHS: pdhg.cc::PDLPSolver::updateIteratesFixed (lines 2080-2134)
    pub(crate) fn updateIteratesFixed(&mut self) {
        // #if PDLP_PROFILE dropped
        let x_next = self.updateX(
            &self.x_current_,
            &self.ATy_cache_,
            self.stepsize_.primal_step,
        );
        self.x_next_ = x_next;

        // #if PDLP_PROFILE dropped
        linalg::ax(&self.lp_, &self.x_next_, &mut self.Ax_next_);

        // #if PDLP_PROFILE dropped
        let y_next = self.updateY(
            &self.y_current_,
            &self.Ax_cache_,
            &self.Ax_next_,
            self.stepsize_.dual_step,
        );
        self.y_next_ = y_next;

        // #if PDLP_PROFILE dropped
        linalg::aTy(&self.lp_, &self.y_next_, &mut self.ATy_next_);

        // #ifdef CUPDLP_GPU dropped
    }

    /// HiGHS: pdhg.cc::PDLPSolver::updateIteratesAdaptive (lines 2136-2305)
    ///
    /// The adaptive line search (`adaptive_linesearch_params`,
    /// `num_rejected_steps_`, `ratio_last_two_step_sizes_`,
    /// `stepsize_.step_size_iter`, `current_eta_`). The stability check and
    /// the iteration cap are commented out upstream and stay out.
    pub(crate) fn updateIteratesAdaptive(&mut self) {
        // #if PDLP_PROFILE dropped

        let _min_eta: f64 = 1e-6;
        let _max_eta: f64 = 1.0;

        let mut accepted_step = false;
        let mut inner_iterations: HighsInt = 0;
        let mut num_rejected_steps: HighsInt = 0;

        let mut d_step_size_update = (self.stepsize_.primal_step * self.stepsize_.dual_step).sqrt();

        // Compute candidate solution
        let x_candidate = self.x_current_.clone(); // Start from current x
        let y_candidate = self.y_current_.clone(); // Start from current y
        let ax_candidate = self.Ax_cache_.clone(); // Start from current ax
        let aty_candidate = self.ATy_cache_.clone(); // Start from current aTy
        let mut xupdate = self.x_next_.clone();
        let mut yupdate = self.y_next_.clone();
        let mut axupdate = self.Ax_next_.clone();
        let mut atyupdate = self.ATy_next_.clone();

        while !accepted_step {
            self.stepsize_.step_size_iter += 1; // nStepSizeIter
            inner_iterations += 1;
            // commented-out max-iteration check dropped

            // Calculate step sizes for this iteration
            let primal_step_update = d_step_size_update / self.stepsize_.beta.sqrt();
            let dual_step_update = d_step_size_update * self.stepsize_.beta.sqrt();

            // #ifdef CUPDLP_GPU dropped
            // Primal update
            // #if PDLP_PROFILE dropped
            xupdate = self.updateX(&x_candidate, &aty_candidate, primal_step_update);
            // #if PDLP_PROFILE dropped

            // #if PDLP_PROFILE dropped
            linalg::ax(&self.lp_, &xupdate, &mut axupdate);
            // #if PDLP_PROFILE dropped

            // Dual update with timing
            // #if PDLP_PROFILE dropped
            yupdate = self.updateY(&y_candidate, &ax_candidate, &axupdate, dual_step_update);
            // #if PDLP_PROFILE dropped

            // #if PDLP_PROFILE dropped
            linalg::aTy(&self.lp_, &yupdate, &mut atyupdate);
            // #if PDLP_PROFILE dropped

            // Compute deltas
            let mut delta_x = vec![0.0; self.lp_.num_col_ as usize];
            let mut delta_y = vec![0.0; self.lp_.num_row_ as usize];
            let mut delta_aty = vec![0.0; self.lp_.num_col_ as usize];

            for i in 0..x_candidate.len() {
                delta_x[i] = xupdate[i] - x_candidate[i];
            }
            for i in 0..y_candidate.len() {
                delta_y[i] = yupdate[i] - y_candidate[i];
            }
            for i in 0..aty_candidate.len() {
                delta_aty[i] = atyupdate[i] - aty_candidate[i];
            }

            // commented-out checkNumericalStability dropped

            // Compute movement and nonlinearity
            let movement = self.computeMovement(&delta_x, &delta_y);
            let nonlinearity = self.computeNonlinearity(&delta_x, &delta_aty);

            // Compute step size limit
            let step_size_limit = if nonlinearity != 0.0 {
                movement / nonlinearity.abs()
            } else {
                f64::INFINITY
            };

            // commented-out highsLogDev dropped

            if d_step_size_update <= step_size_limit {
                accepted_step = true;
            } else {
                num_rejected_steps += 1;
            }

            // Compute new step size
            let first_term = if step_size_limit.is_infinite() {
                step_size_limit
            } else {
                (1.0 - (self.stepsize_.step_size_iter as f64 + 1.0).powf(
                    -self
                        .params_
                        .adaptive_linesearch_params
                        .step_size_reduction_exponent,
                )) * step_size_limit
            };

            let second_term =
                (1.0 + (self.stepsize_.step_size_iter as f64 + 1.0).powf(
                    -self
                        .params_
                        .adaptive_linesearch_params
                        .step_size_growth_exponent,
                )) * d_step_size_update;

            d_step_size_update = cpp_min(first_term, second_term);
            // commented-out current_eta_ clamp dropped
        }

        let _ = (inner_iterations, num_rejected_steps); // C++ counters, read only by commented-out code (pdhg.cc 2163-2170)

        self.x_next_ = xupdate;
        self.y_next_ = yupdate;
        self.Ax_next_ = axupdate;
        self.ATy_next_ = atyupdate;

        self.current_eta_ = d_step_size_update;
        self.stepsize_.primal_step = d_step_size_update / self.stepsize_.beta.sqrt();
        self.stepsize_.dual_step = d_step_size_update * self.stepsize_.beta.sqrt();

        // #if PDLP_PROFILE dropped
    }

    /// HiGHS: pdhg.cc::PDLPSolver::updateIteratesMalitskyPock (lines 2307-2310)
    ///
    /// Upstream body is a stub; translate it as written.
    pub(crate) fn updateIteratesMalitskyPock(&mut self, first_malitsky_iteration: bool) -> bool {
        let _ = first_malitsky_iteration;
        let _x_candidate = vec![0.0; self.lp_.num_col_ as usize];
        true
    }

    /// HiGHS: pdhg.cc::PDLPSolver::checkNumericalStability (lines 2316-2332)
    ///
    /// Compares the movement with [`kDivergentMovement`].
    pub(crate) fn checkNumericalStability(&self, delta_x: &[f64], delta_y: &[f64]) -> bool {
        let movement = self.computeMovement(delta_x, delta_y);

        if movement == 0.0 {
            highs_log_user(
                &self.log_options_,
                HighsLogType::kInfo,
                "Warning: Zero movement detected - numerical termination\n",
            );
            return false;
        }

        if movement > kDivergentMovement {
            highs_log_user(
                &self.log_options_,
                HighsLogType::kInfo,
                &format!(
                    "Warning: Divergent movement detected: {}\n",
                    c_fmt_g(movement, 6)
                ),
            );
            return false;
        }
        true // Placeholder
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeMovement (lines 2334-2349)
    pub(crate) fn computeMovement(&self, delta_primal: &[f64], delta_dual: &[f64]) -> f64 {
        let mut primal_squared_norm: f64 = 0.0;
        for val in delta_primal {
            primal_squared_norm += val * val;
        }

        let mut dual_squared_norm: f64 = 0.0;
        for val in delta_dual {
            dual_squared_norm += val * val;
        }

        let primal_weight = self.stepsize_.beta.sqrt();
        (0.5 * primal_weight * primal_squared_norm) + (0.5 / primal_weight) * dual_squared_norm
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeNonlinearity (lines 2351-2359)
    pub(crate) fn computeNonlinearity(&self, delta_primal: &[f64], delta_aty: &[f64]) -> f64 {
        // Nonlinearity = |Dx' * D(A'y)|
        let mut nonlinearity: f64 = 0.0;
        for i in 0..delta_primal.len() {
            nonlinearity += delta_primal[i] * delta_aty[i];
        }
        nonlinearity // cupdlp does not take absolute value
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn step_sizes_from_norms() {
        let mut s = PDLPSolver::default();
        s.unscaled_c_norm_ = 3.0;
        s.unscaled_rhs_norm_ = 1.0;
        s.params_.step_size_strategy = StepSizeStrategy::ADAPTIVE;
        // Empty LP: powerMethod returns 1.
        s.initializeStepSizes();
        assert_eq!(s.params_.omega, 2.0);
        assert_eq!(s.stepsize_.beta, 4.0);
        assert_eq!(s.params_.step_size_strategy, StepSizeStrategy::FIXED);
        assert_eq!(s.params_.eta, 0.998);
        assert_eq!(s.stepsize_.primal_step, 0.998 / 2.0);
        assert_eq!(s.stepsize_.dual_step, 0.998 * 2.0);
    }

    #[test]
    fn primal_weight_reverts_to_best_on_bad_norms() {
        let mut s = PDLPSolver::default();
        s.primal_weight_ = 5.0;
        s.best_primal_weight_ = 2.0;
        s.stepsize_.primal_step = 1.0;
        s.stepsize_.dual_step = 4.0;
        // Empty LP: both distances are 0 -> revert branch.
        s.updatePrimalWeightAtRestart(&SolverResults::default());
        assert_eq!(s.primal_weight_, 2.0);
        assert_eq!(s.params_.omega, 2.0);
        assert_eq!(s.stepsize_.primal_step, 1.0);
        assert_eq!(s.stepsize_.dual_step, 4.0);
        assert_eq!(s.restart_scheme_.getBeta(), 4.0);
    }

    #[test]
    fn movement_and_nonlinearity() {
        let s = PDLPSolver::default();
        let _ = (
            s.computeMovement(&[], &[]),
            linalg::norm2(&[]),
            kDivergentMovement,
        );
        assert!(s.checkNumericalStability(&[], &[]));
    }

    #[test]
    fn update_x_y_1col_1row() {
        // 1 column, 1 row LP.
        // col_cost = [2.0], col_lower = [-1.0], col_upper = [1.0]
        // row_lower = [0.0], row_upper = [1.0] (inequality, not equality)
        // x = [0.5], aty = [1.0], primal_step = 0.5
        // gradient = 2.0 - 1.0 = 1.0
        // x_new = projectBox(0.5 - 0.5*1.0, -1.0, 1.0) = projectBox(0.0, -1.0, 1.0) = 0.0
        let mut s = PDLPSolver::default();
        s.lp_.num_col_ = 1;
        s.lp_.num_row_ = 1;
        s.lp_.col_cost_ = vec![2.0];
        s.lp_.col_lower_ = vec![-1.0];
        s.lp_.col_upper_ = vec![1.0];
        s.lp_.row_lower_ = vec![0.0];
        s.lp_.row_upper_ = vec![1.0];

        let x = vec![0.5];
        let aty = vec![1.0];
        let x_new = s.updateX(&x, &aty, 0.5);
        assert_eq!(x_new, vec![0.0]);

        // updateY: y = [1.0], ax = [2.0], ax_next = [3.0], dual_step = 0.5
        // extr_ax = 2*3.0 - 2.0 = 4.0
        // is_equality = (0.0 == 1.0) = false
        // q = 0.0
        // dual_update = 1.0 + 0.5*(0.0 - 4.0) = 1.0 - 2.0 = -1.0
        // projectNonNegative(-1.0) = 0.0
        let y = vec![1.0];
        let ax = vec![2.0];
        let ax_next = vec![3.0];
        let y_new = s.updateY(&y, &ax, &ax_next, 0.5);
        assert_eq!(y_new, vec![0.0]);

        // Equality row: row_lower == row_upper
        s.lp_.row_lower_ = vec![1.0];
        s.lp_.row_upper_ = vec![1.0];
        // is_equality = true
        // q = 1.0
        // dual_update = 1.0 + 0.5*(1.0 - 4.0) = 1.0 - 1.5 = -0.5
        // y_new = -0.5 (no projection)
        let y_new_eq = s.updateY(&y, &ax, &ax_next, 0.5);
        assert_eq!(y_new_eq, vec![-0.5]);
    }

    #[test]
    fn compute_movement_and_nonlinearity_values() {
        let mut s = PDLPSolver::default();
        s.stepsize_.beta = 4.0; // sqrt(beta) = 2.0

        // delta_primal = [1.0, 2.0], delta_dual = [3.0, 4.0]
        // primal_squared_norm = 1 + 4 = 5
        // dual_squared_norm = 9 + 16 = 25
        // primal_weight = 2.0
        // movement = 0.5*2.0*5 + 0.5/2.0*25 = 5.0 + 6.25 = 11.25
        let dp = vec![1.0, 2.0];
        let dd = vec![3.0, 4.0];
        let mv = s.computeMovement(&dp, &dd);
        assert_eq!(mv, 11.25);

        // delta_primal = [1.0, 2.0], delta_aty = [3.0, 4.0]
        // nonlinearity = 1*3 + 2*4 = 3 + 8 = 11
        let da = vec![3.0, 4.0];
        let nl = s.computeNonlinearity(&dp, &da);
        assert_eq!(nl, 11.0);
    }

    #[test]
    fn check_numerical_stability_false_cases() {
        let mut s = PDLPSolver::default();
        s.stepsize_.beta = 1.0; // sqrt(beta) = 1.0

        // Zero movement: delta_x = [], delta_y = []
        // movement = 0.0 -> false
        assert!(!s.checkNumericalStability(&[], &[]));

        // Divergent movement: need movement > kDivergentMovement
        // movement = 0.5*1.0*sum(dp^2) + 0.5/1.0*sum(dd^2)
        // Use large values to exceed kDivergentMovement
        let big = vec![1e200; 1];
        assert!(!s.checkNumericalStability(&big, &[]));
    }
}
