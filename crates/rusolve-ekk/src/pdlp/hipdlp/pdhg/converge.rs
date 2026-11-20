//! HiGHS source: pdlp/hipdlp/pdhg.cc (lines 1138-1527)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Vector setup, averaging, and the residual / objective / convergence
//! computations. `hipdlpTimerStart` / `Stop` calls are dropped
//! (`PDLP_PROFILE == 0`).

use super::super::defs::*;
use super::super::linalg;
use super::super::solver_results::*;
use super::PDLPSolver;
use crate::prelude::*;

impl<'a> PDLPSolver<'a> {
    /// HiGHS: pdhg.cc::PDLPSolver::initialize (lines 1138-1177)
    ///
    /// `resize(n, 0.0)` keeps existing entries (as `std::vector::resize`):
    /// use `Vec::resize`, not `vec![0.0; n]`.
    pub(crate) fn initialize(&mut self) {
        // initialize x and y based on the LP problem
        self.x_current_.resize(self.lp_.num_col_ as usize, 0.0);
        self.y_current_.resize(self.lp_.num_row_ as usize, 0.0);
        self.x_next_.resize(self.lp_.num_col_ as usize, 0.0);
        self.y_next_.resize(self.lp_.num_row_ as usize, 0.0);

        self.x_at_last_restart_ = self.x_current_.clone();
        self.y_at_last_restart_ = self.y_current_.clone();
        self.x_sum_.resize(self.lp_.num_col_ as usize, 0.0);
        self.y_sum_.resize(self.lp_.num_row_ as usize, 0.0);

        self.sum_weights_ = 0.0;

        self.x_avg_.resize(self.lp_.num_col_ as usize, 0.0);
        self.y_avg_.resize(self.lp_.num_row_ as usize, 0.0);

        self.reflected_x_.resize(self.lp_.num_col_ as usize, 0.0);
        self.reflected_y_.resize(self.lp_.num_row_ as usize, 0.0);
        self.x_anchor_.resize(self.lp_.num_col_ as usize, 0.0);
        self.y_anchor_.resize(self.lp_.num_row_ as usize, 0.0);
        self.halpern_iteration_ = 0;

        self.Ax_cache_.resize(self.lp_.num_row_ as usize, 0.0);
        self.ATy_cache_.resize(self.lp_.num_col_ as usize, 0.0);
        self.Ax_next_.resize(self.lp_.num_row_ as usize, 0.0);
        self.ATy_next_.resize(self.lp_.num_col_ as usize, 0.0);

        // Member variables added back in HPP:
        self.Ax_avg_.resize(self.lp_.num_row_ as usize, 0.0);
        self.ATy_avg_.resize(self.lp_.num_col_ as usize, 0.0);
        self.K_times_x_diff_.resize(self.lp_.num_row_ as usize, 0.0);

        self.dSlackPos_.resize(self.lp_.num_col_ as usize, 0.0);
        self.dSlackNeg_.resize(self.lp_.num_col_ as usize, 0.0);
        self.dSlackPosAvg_.resize(self.lp_.num_col_ as usize, 0.0);
        self.dSlackNegAvg_.resize(self.lp_.num_col_ as usize, 0.0);
        self.halpern_dual_slack_next_
            .resize(self.lp_.num_col_ as usize, 0.0);
        self.halpern_dual_slack_next_valid_ = false;
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeStepSizeRatio (lines 1180-1203)
    ///
    /// Called with `params_` only (by the GPU path); `working_params` is
    /// therefore `self.params_` and the method takes no parameter. The
    /// final `restart_scheme_.updateBeta` writes `params_.omega` again.
    pub(crate) fn computeStepSizeRatio(&mut self) {
        // 1. Calculate the L2 norm of the difference between current and last-restart
        // iterates.
        let primal_diff_norm = linalg::diffTwoNorm(&self.x_at_last_restart_, &self.x_current_);
        let dual_diff_norm = linalg::diffTwoNorm(&self.y_at_last_restart_, &self.y_current_);

        let dMeanStepSize = (self.stepsize_.primal_step * self.stepsize_.dual_step).sqrt();

        // 2. Update the primal weight (beta = omega^2) if movements are significant.
        if cpp_min(primal_diff_norm, dual_diff_norm) > 1e-10 {
            let beta_update_ratio = dual_diff_norm / primal_diff_norm;
            let old_beta = self.stepsize_.beta;

            let dLogBetaUpdate = 0.5 * beta_update_ratio.ln() + 0.5 * old_beta.sqrt().ln();
            self.stepsize_.beta = (2.0 * dLogBetaUpdate).exp();
        }

        self.stepsize_.primal_step = dMeanStepSize / self.stepsize_.beta.sqrt();
        self.stepsize_.dual_step = self.stepsize_.primal_step * self.stepsize_.beta;
        self.params_.eta = (self.stepsize_.primal_step * self.stepsize_.dual_step).sqrt();
        self.params_.omega = self.stepsize_.beta.sqrt();
        self.restart_scheme_
            .update_beta(self.stepsize_.beta, &mut self.params_);
    }

    /// HiGHS: pdhg.cc::PDLPSolver::updateAverageIterates (lines 1205-1227)
    ///
    /// `x` / `y` alias `x_current_` / `x_next_` at the call site: the caller
    /// takes them.
    pub(crate) fn updateAverageIterates(&mut self, x: &[f64], y: &[f64], inner_iter: HighsInt) {
        let dMeanStepSize = (self.stepsize_.primal_step * self.stepsize_.dual_step).sqrt();

        // hipdlpTimerStart/Stop dropped (PDLP_PROFILE == 0), lines 1210-1216.
        for i in 0..x.len() {
            self.x_sum_[i] += x[i] * dMeanStepSize;
        }

        // hipdlpTimerStart/Stop dropped (PDLP_PROFILE == 0), lines 1218-1224.
        for i in 0..y.len() {
            self.y_sum_[i] += y[i] * dMeanStepSize;
        }

        self.sum_weights_ += dMeanStepSize;
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeAverageIterate (lines 1229-1275)
    ///
    /// Called as `computeAverageIterate(Ax_avg_, ATy_avg_)`: the caller
    /// takes the two members.
    pub(crate) fn computeAverageIterate(&mut self, ax_avg: &mut [f64], aty_avg: &mut [f64]) {
        let dPrimalScale = if self.sum_weights_ > 1e-10 {
            1.0 / self.sum_weights_
        } else {
            1.0
        };
        let dDualScale = if self.sum_weights_ > 1e-10 {
            1.0 / self.sum_weights_
        } else {
            1.0
        };

        // hipdlpTimerStart/Stop dropped (PDLP_PROFILE == 0), lines 1234-1241.
        for i in 0..self.x_avg_.len() {
            self.x_avg_[i] = self.x_sum_[i] * dPrimalScale;
        }

        // hipdlpTimerStart/Stop dropped (PDLP_PROFILE == 0), lines 1243-1249.
        for i in 0..self.y_avg_.len() {
            self.y_avg_[i] = self.y_sum_[i] * dDualScale;
        }

        // #if PDLP_DEBUG_LOG dropped (lines 1251-1253).

        // hipdlpTimerStart/Stop dropped (PDLP_PROFILE == 0), lines 1255-1261.
        linalg::ax(&self.lp_, &self.x_avg_, ax_avg);

        // hipdlpTimerStart/Stop dropped (PDLP_PROFILE == 0), lines 1263-1269.
        linalg::aTy(&self.lp_, &self.y_avg_, aty_avg);

        // #if PDLP_DEBUG_LOG dropped (lines 1271-1274).
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeLambda (lines 1278-1295)
    pub(crate) fn computeLambda(&self, y: &[f64], ATy_vector: &[f64]) -> Vec<f64> {
        let mut lambda = vec![0.0; self.lp_.num_col_ as usize];
        for i in 0..self.lp_.num_col_ {
            let residual = self.lp_.col_cost_[i as usize] - ATy_vector[i as usize];

            if self.lp_.col_lower_[i as usize] <= -kHighsInf
                && self.lp_.col_upper_[i as usize] >= kHighsInf
            {
                lambda[i as usize] = 0.0;
            } else if self.lp_.col_lower_[i as usize] <= -kHighsInf {
                lambda[i as usize] = cpp_min(0.0, residual);
            } else if self.lp_.col_upper_[i as usize] >= kHighsInf {
                lambda[i as usize] = cpp_max(0.0, residual);
            } else {
                lambda[i as usize] = residual;
            }
        }
        lambda
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computePrimalFeasibility (lines 1297-1320)
    pub(crate) fn computePrimalFeasibility(&self, Ax_vector: &[f64]) -> f64 {
        let mut primal_residual = vec![0.0; self.lp_.num_row_ as usize];

        // Compute ax - rhs
        for i in 0..self.lp_.num_row_ {
            primal_residual[i as usize] = Ax_vector[i as usize] - self.lp_.row_lower_[i as usize];

            // For inequality constraints, project to negative part
            if !self.is_equality_row_[i as usize] {
                primal_residual[i as usize] = cpp_min(0.0, primal_residual[i as usize]);
            }
        }

        // Apply scaling if needed
        if self.scaling_.isScaled() {
            let row_scale = self.scaling_.getRowScaling();
            for i in 0..self.lp_.num_row_ {
                primal_residual[i as usize] *= row_scale[i as usize];
            }
        }

        linalg::vectorNorm(&primal_residual, 2.0)
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeDualSlacks (lines 1322-1378)
    ///
    /// `slacks_are_members` is the C++ pointer test `(&dSlackPos ==
    /// &dSlackPos_) && (&dSlackNeg == &dSlackNeg_)` (see [`super`]); the
    /// Halpern cached slack is used only when it is true,
    /// `halpern_dual_slack_next_valid_` and the sizes match.
    pub(crate) fn computeDualSlacks(
        &self,
        dualResidual: &[f64],
        dSlackPos: &mut Vec<f64>,
        dSlackNeg: &mut Vec<f64>,
        slacks_are_members: bool,
    ) {
        // Ensure vectors are correctly sized
        if dSlackPos.len() != self.lp_.num_col_ as usize {
            dSlackPos.resize(self.lp_.num_col_ as usize, 0.0);
        }
        if dSlackNeg.len() != self.lp_.num_col_ as usize {
            dSlackNeg.resize(self.lp_.num_col_ as usize, 0.0);
        }

        for i in 0..self.lp_.num_col_ {
            if self.params_.use_halpern_restart {
                // For Halpern current iterate checks, use major-step dual slack:
                // dual_slack = (proj(temp) - temp) / primal_step.
                // For other cases (e.g. average iterate), fall back to sign/bounds
                // projection.
                let use_cached_halpern_slack = slacks_are_members
                    && self.halpern_dual_slack_next_valid_
                    && (self.halpern_dual_slack_next_.len() == self.lp_.num_col_ as usize);

                let mut dual_slack = 0.0;
                if use_cached_halpern_slack {
                    dual_slack = self.halpern_dual_slack_next_[i as usize];
                } else {
                    let has_lower = self.lp_.col_lower_[i as usize] > -kHighsInf;
                    let has_upper = self.lp_.col_upper_[i as usize] < kHighsInf;
                    if has_lower && has_upper {
                        dual_slack = dualResidual[i as usize];
                    } else if has_lower {
                        dual_slack = cpp_max(0.0, dualResidual[i as usize]);
                    } else if has_upper {
                        dual_slack = cpp_min(0.0, dualResidual[i as usize]);
                    }
                }

                dSlackPos[i as usize] = cpp_max(0.0, dual_slack);
                dSlackNeg[i as usize] = cpp_max(0.0, -dual_slack);
            } else {
                // Compute positive slack (for lower bounds)
                // CUPDLP: max(dual_residual, 0) * hasLower
                if self.lp_.col_lower_[i as usize] > -kHighsInf {
                    dSlackPos[i as usize] = cpp_max(0.0, dualResidual[i as usize]);
                } else {
                    dSlackPos[i as usize] = 0.0;
                }

                // Compute negative slack (for upper bounds)
                // CUPDLP: -min(dual_residual, 0) * hasUpper
                if self.lp_.col_upper_[i as usize] < kHighsInf {
                    dSlackNeg[i as usize] = cpp_max(0.0, -dualResidual[i as usize]);
                } else {
                    dSlackNeg[i as usize] = 0.0;
                }
            }
        }
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeDualFeasibility (lines 1380-1408)
    ///
    /// `slacks_are_members`: see [`Self::computeDualSlacks`].
    pub(crate) fn computeDualFeasibility(
        &self,
        ATy_vector: &[f64],
        dSlackPos: &mut Vec<f64>,
        dSlackNeg: &mut Vec<f64>,
        slacks_are_members: bool,
    ) -> f64 {
        // dualResidual = c-A'y
        let dualResidual = linalg::vectorSubtract(&self.lp_.col_cost_, ATy_vector);

        // Call the refactored function to populate dSlackPos and dSlackNeg
        self.computeDualSlacks(&dualResidual, dSlackPos, dSlackNeg, slacks_are_members);

        let mut dual_residual = vec![0.0; self.lp_.num_col_ as usize];

        for i in 0..self.lp_.num_col_ {
            // Matching CUPDLP: c - A'y - dSlackPos + dSlackNeg
            dual_residual[i as usize] =
                dualResidual[i as usize] - dSlackPos[i as usize] + dSlackNeg[i as usize];
        }

        // Apply scaling if needed
        if self.scaling_.isScaled() {
            let col_scale = self.scaling_.getColScaling();
            for i in 0..self.lp_.num_col_ {
                dual_residual[i as usize] *= col_scale[i as usize];
            }
        }

        linalg::vectorNorm(&dual_residual, 2.0)
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeDualityGap (lines 1410-1445)
    ///
    /// Returns the C++ `std::tuple<double, double, double, double, double>`
    /// in the same order.
    pub(crate) fn computeDualityGap(
        &self,
        x: &[f64],
        y: &[f64],
        lambda: &[f64],
    ) -> (f64, f64, f64, f64, f64) {
        let mut qTy = 0.0;
        for i in 0..self.lp_.num_row_ {
            if self.lp_.row_lower_[i as usize] > -kHighsInf {
                qTy += self.lp_.row_lower_[i as usize] * y[i as usize];
            } else {
                qTy += self.lp_.row_upper_[i as usize] * y[i as usize];
            }
        }

        let mut lT_lambda_plus = 0.0;
        let mut uT_lambda_minus = 0.0;

        for i in 0..self.lp_.num_col_ {
            if self.lp_.col_lower_[i as usize] > -kHighsInf {
                lT_lambda_plus +=
                    self.lp_.col_lower_[i as usize] * cpp_max(0.0, lambda[i as usize]);
            }
            if self.lp_.col_upper_[i as usize] < kHighsInf {
                uT_lambda_minus +=
                    self.lp_.col_upper_[i as usize] * cpp_max(0.0, -lambda[i as usize]);
            }
        }

        let dual_objective = qTy + lT_lambda_plus - uT_lambda_minus;

        let mut cTx = 0.0;
        for i in 0..self.lp_.num_col_ {
            cTx += self.lp_.col_cost_[i as usize] * x[i as usize];
        }

        let duality_gap = (dual_objective - cTx).abs();
        (duality_gap, qTy, lT_lambda_plus, uT_lambda_minus, cTx)
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeDualObjective (lines 1447-1472)
    pub(crate) fn computeDualObjective(
        &self,
        y: &[f64],
        dSlackPos: &[f64],
        dSlackNeg: &[f64],
    ) -> f64 {
        let mut dual_obj = self.lp_.offset_;

        // Compute b'y (or rhs'y in cuPDLP notation)
        for i in 0..self.lp_.num_row_ {
            dual_obj += self.lp_.row_lower_[i as usize] * y[i as usize];
        }

        // Add contribution from lower bounds: l'*slackPos
        for i in 0..self.lp_.num_col_ {
            if self.lp_.col_lower_[i as usize] > -kHighsInf {
                dual_obj += self.lp_.col_lower_[i as usize] * dSlackPos[i as usize];
            }
        }

        // Subtract contribution from upper bounds: u'*slackNeg
        for i in 0..self.lp_.num_col_ {
            if self.lp_.col_upper_[i as usize] < kHighsInf {
                dual_obj -= self.lp_.col_upper_[i as usize] * dSlackNeg[i as usize];
            }
        }

        dual_obj
    }

    /// HiGHS: pdhg.cc::PDLPSolver::checkConvergence (lines 1474-1527)
    ///
    /// `type_` is the C++ `const char* type` (`"[L]"` / `"[A]"`, logging
    /// only). `slacks_are_members`: see [`Self::computeDualSlacks`].
    pub(crate) fn checkConvergence(
        &self,
        iter: HighsInt,
        x: &[f64],
        y: &[f64],
        ax_vector: &[f64],
        aty_vector: &[f64],
        epsilon: f64,
        results: &mut SolverResults,
        type_: &str,
        dSlackPos: &mut Vec<f64>,
        dSlackNeg: &mut Vec<f64>,
        slacks_are_members: bool,
    ) -> bool {
        // computeDualSlacks is now called inside computeDualFeasibility

        // Compute primal feasibility
        let primal_feasibility = self.computePrimalFeasibility(ax_vector);
        results.primal_feasibility = primal_feasibility;

        // Compute dual feasibility
        // This will populate dSlackPos and dSlackNeg
        let dual_feasibility =
            self.computeDualFeasibility(aty_vector, dSlackPos, dSlackNeg, slacks_are_members);
        results.dual_feasibility = dual_feasibility;

        // Compute objectives
        let mut primal_obj = self.lp_.offset_;
        for i in 0..self.lp_.num_col_ {
            primal_obj += self.lp_.col_cost_[i as usize] * x[i as usize];
        }
        results.primal_obj = primal_obj;

        // Pass the now-populated slack vectors to computeDualObjective
        let dual_obj = self.computeDualObjective(y, dSlackPos, dSlackNeg);
        results.dual_obj = dual_obj;

        // Compute duality gap
        let duality_gap = primal_obj - dual_obj;
        results.duality_gap = duality_gap.abs();

        // Compute relative gap (matching cuPDLP formula)
        let relative_obj_gap = duality_gap.abs() / (1.0 + primal_obj.abs() + dual_obj.abs());
        results.relative_obj_gap = relative_obj_gap;

        // #if PDLP_DEBUG_LOG dropped (lines 1513-1518).

        // check convergence criteria (matching cuPDLP)
        let primal_feasible = primal_feasibility < epsilon * (1.0 + self.unscaled_rhs_norm_);
        let dual_feasible = dual_feasibility < epsilon * (1.0 + self.unscaled_c_norm_);
        let gap_small = relative_obj_gap < epsilon;

        primal_feasible && dual_feasible && gap_small
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    fn solver(lp: &HighsLp) -> PDLPSolver<'_> {
        let mut s = PDLPSolver::default();
        s.passLp(lp);
        s.preprocessLp();
        s.initialize();
        s
    }

    #[test]
    fn initialize_sizes_every_vector() {
        let lp = super::super::tests::four_row_lp();
        let s = solver(&lp);
        let (n, m) = (s.lp_.num_col_ as usize, s.lp_.num_row_ as usize);
        assert_eq!((s.x_current_.len(), s.y_current_.len()), (n, m));
        assert_eq!((s.Ax_cache_.len(), s.ATy_cache_.len()), (m, n));
        assert_eq!(s.halpern_dual_slack_next_.len(), n);
        assert!(!s.halpern_dual_slack_next_valid_);
    }

    #[test]
    fn lambda_respects_bound_types() {
        let lp = super::super::tests::four_row_lp();
        let s = solver(&lp);
        let n = s.lp_.num_col_ as usize;
        let lambda = s.computeLambda(&vec![0.0; s.lp_.num_row_ as usize], &vec![0.0; n]);
        // Boxed columns: lambda = c - A'y = c.
        assert_eq!(lambda[0], 1.0);
        assert_eq!(lambda[1], 1.0);
    }

    #[test]
    fn zero_point_is_not_converged() {
        let lp = super::super::tests::four_row_lp();
        let s = solver(&lp);
        let (n, m) = (s.lp_.num_col_ as usize, s.lp_.num_row_ as usize);
        let mut r = SolverResults::default();
        let (mut sp, mut sn) = (vec![0.0; n], vec![0.0; n]);
        let conv = s.checkConvergence(
            0,
            &vec![0.0; n],
            &vec![0.0; m],
            &vec![0.0; m],
            &vec![0.0; n],
            1e-4,
            &mut r,
            "[L]",
            &mut sp,
            &mut sn,
            false,
        );
        assert!(!conv);
        assert!(r.primal_feasibility > 0.0);
        let _ = (linalg::norm2(&sp), StepSizeStrategy::FIXED);
    }

    /// Two columns, two rows, set directly on `lp_` (no preprocessing):
    /// col0 in [0, 10], col1 in (-inf, 5]; r0 >= 1, r1 <= 2.
    /// A = [[1, 1], [1, -1]].
    fn raw_solver() -> PDLPSolver<'static> {
        let mut s = PDLPSolver::default();
        let lp = &mut s.lp_;
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![1.0, 2.0];
        lp.col_lower_ = vec![0.0, -kHighsInf];
        lp.col_upper_ = vec![10.0, 5.0];
        lp.row_lower_ = vec![1.0, -kHighsInf];
        lp.row_upper_ = vec![kHighsInf, 2.0];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, -1.0];
        s.is_equality_row_ = vec![false, false];
        s.initialize();
        s
    }

    #[test]
    fn initialize_keeps_existing_entries() {
        // pdhg.cc:1140 `x_current_.resize(n, 0.0)` keeps entries, and
        // x_at_last_restart_ copies x_current_ (1145).
        let mut s = raw_solver();
        s.x_current_ = vec![3.0];
        s.initialize();
        assert_eq!(s.x_current_, vec![3.0, 0.0]);
        assert_eq!(s.x_at_last_restart_, vec![3.0, 0.0]);
        assert_eq!(s.sum_weights_, 0.0);
    }

    #[test]
    fn lambda_one_sided_and_free_columns() {
        let mut s = raw_solver();
        // c - A'y with A'y = 0 is c = [1, 2]. col0 boxed -> residual;
        // col1 lower -inf, upper 5 -> min(0, 2) = 0.
        let lambda = s.computeLambda(&[0.0, 0.0], &[0.0, 0.0]);
        assert_eq!(lambda, vec![1.0, 0.0]);
        s.lp_.col_upper_ = vec![kHighsInf, 5.0];
        let lambda = s.computeLambda(&[0.0, 0.0], &[-3.0, 4.0]);
        // col0 lower 0, upper inf: max(0, 1 + 3) = 4; col1: min(0, 2 - 4) = -2.
        assert_eq!(lambda, vec![4.0, -2.0]);
        // col0 free -> 0.
        s.lp_.col_lower_ = vec![-kHighsInf, -kHighsInf];
        let lambda = s.computeLambda(&[0.0, 0.0], &[3.0, 1.0]);
        assert_eq!(lambda, vec![0.0, 0.0]);
    }

    #[test]
    fn dual_slacks_non_halpern_bounds() {
        let s = raw_solver();
        let (mut sp, mut sn) = (Vec::new(), Vec::new());
        s.computeDualSlacks(&[-1.5, 2.0], &mut sp, &mut sn, false);
        // col0 has both bounds: pos = max(0, -1.5) = 0, neg = max(0, 1.5).
        // col1 has no lower: pos = 0; upper 5: neg = max(0, -2) = 0.
        assert_eq!(sp, vec![0.0, 0.0]);
        assert_eq!(sn, vec![1.5, 0.0]);
    }

    #[test]
    fn dual_slacks_halpern_cached_only_for_members() {
        let mut s = raw_solver();
        s.params_.use_halpern_restart = true;
        s.halpern_dual_slack_next_ = vec![2.0, -3.0];
        s.halpern_dual_slack_next_valid_ = true;
        let (mut sp, mut sn) = (vec![0.0; 2], vec![0.0; 2]);
        // Members: the cached slack is used.
        s.computeDualSlacks(&[-1.5, 4.0], &mut sp, &mut sn, true);
        assert_eq!(sp, vec![2.0, 0.0]);
        assert_eq!(sn, vec![0.0, 3.0]);
        // Not members: projection. col0 boxed -> -1.5; col1 upper only ->
        // min(0, 4) = 0.
        s.computeDualSlacks(&[-1.5, 4.0], &mut sp, &mut sn, false);
        assert_eq!(sp, vec![0.0, 0.0]);
        assert_eq!(sn, vec![1.5, 0.0]);
        // Members but cache invalid: projection again.
        s.halpern_dual_slack_next_valid_ = false;
        s.computeDualSlacks(&[-1.5, -4.0], &mut sp, &mut sn, true);
        assert_eq!(sp, vec![0.0, 0.0]);
        assert_eq!(sn, vec![1.5, 4.0]);
    }

    #[test]
    fn duality_gap_tuple() {
        let s = raw_solver();
        // qTy = 1*3 + 2*4 = 11 (r1 has no lower: uses row_upper_).
        // lT_lambda_plus = 0*1.5 (col1 has no lower) = 0.
        // uT_lambda_minus = 10*max(0,-1.5) + 5*max(0,2) = 10.
        // dual = 11 + 0 - 10 = 1; cTx = 1 + 2 = 3; gap = |1 - 3| = 2.
        let t = s.computeDualityGap(&[1.0, 1.0], &[3.0, 4.0], &[1.5, -2.0]);
        assert_eq!(t, (2.0, 11.0, 0.0, 10.0, 3.0));
    }

    #[test]
    fn dual_objective_uses_row_lower_and_bounds() {
        let mut s = raw_solver();
        s.lp_.offset_ = 0.5;
        s.lp_.row_lower_ = vec![1.0, -1.0];
        // 0.5 + 1*3 + (-1)*2 = 1.5; + l'pos (col0 only) 0*2 = 1.5;
        // - u'neg 10*1 - 5*2 = -18.5.
        let d = s.computeDualObjective(&[3.0, 2.0], &[2.0, 7.0], &[1.0, 2.0]);
        assert_eq!(d, -18.5);
    }

    #[test]
    fn averages_weight_by_mean_step() {
        let mut s = raw_solver();
        s.stepsize_.primal_step = 4.0;
        s.stepsize_.dual_step = 1.0;
        s.updateAverageIterates(&[1.0, 3.0], &[2.0, -1.0], 0);
        assert_eq!(s.x_sum_, vec![2.0, 6.0]);
        assert_eq!(s.y_sum_, vec![4.0, -2.0]);
        assert_eq!(s.sum_weights_, 2.0);
        let (mut ax, mut aty) = (vec![9.0; 2], vec![9.0; 2]);
        s.computeAverageIterate(&mut ax, &mut aty);
        assert_eq!(s.x_avg_, vec![1.0, 3.0]);
        assert_eq!(s.y_avg_, vec![2.0, -1.0]);
        // A x = [1 + 3, 1 - 3]; A'y = [2 - 1, 2 + 1].
        assert_eq!(ax, vec![4.0, -2.0]);
        assert_eq!(aty, vec![1.0, 3.0]);
    }

    #[test]
    fn step_size_ratio_updates_beta_when_both_moved() {
        let mut s = raw_solver();
        s.stepsize_.primal_step = 1.0;
        s.stepsize_.dual_step = 1.0;
        s.stepsize_.beta = 1.0;
        s.x_current_ = vec![1.0, 0.0];
        s.y_current_ = vec![0.0, 4.0];
        s.computeStepSizeRatio();
        // ratio 4: beta = exp(2 * (0.5 ln 4 + 0.5 ln 1)) = 4.
        let tol = 1e-12;
        assert!((s.stepsize_.beta - 4.0).abs() < tol);
        assert!((s.stepsize_.primal_step - 0.5).abs() < tol);
        assert!((s.stepsize_.dual_step - 2.0).abs() < tol);
        assert!((s.params_.eta - 1.0).abs() < tol);
        assert!((s.params_.omega - 2.0).abs() < tol);
        assert_eq!(s.restart_scheme_.getBeta(), s.stepsize_.beta);
    }

    #[test]
    fn step_size_ratio_keeps_beta_without_movement() {
        let mut s = raw_solver();
        s.stepsize_.primal_step = 1.0;
        s.stepsize_.dual_step = 4.0;
        s.stepsize_.beta = 4.0;
        s.y_current_ = vec![0.0, 4.0]; // primal did not move
        s.computeStepSizeRatio();
        assert_eq!(s.stepsize_.beta, 4.0);
        // mean step sqrt(4) = 2; primal = 2 / sqrt(4).
        assert_eq!(s.stepsize_.primal_step, 1.0);
        assert_eq!(s.stepsize_.dual_step, 4.0);
        assert_eq!(s.params_.omega, 2.0);
    }

    #[test]
    fn convergence_at_optimum() {
        // min x s.t. x = 1, 0 <= x <= 10: x = 1, y = 1 is optimal.
        let mut s = PDLPSolver::default();
        let lp = &mut s.lp_;
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![1.0];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        s.is_equality_row_ = vec![true];
        s.initialize();
        let mut r = SolverResults::default();
        let (mut sp, mut sn) = (Vec::new(), Vec::new());
        let conv = s.checkConvergence(
            0,
            &[1.0],
            &[1.0],
            &[1.0],
            &[1.0],
            1e-4,
            &mut r,
            "[A]",
            &mut sp,
            &mut sn,
            false,
        );
        assert!(conv);
        assert_eq!(r.primal_feasibility, 0.0);
        assert_eq!(r.dual_feasibility, 0.0);
        assert_eq!((r.primal_obj, r.dual_obj), (1.0, 1.0));
        assert_eq!((r.duality_gap, r.relative_obj_gap), (0.0, 0.0));
        assert_eq!((sp, sn), (vec![0.0], vec![0.0]));
    }
}
