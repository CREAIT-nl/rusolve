//! HiGHS source: pdlp/hipdlp/pdhg.cc (lines 784-1136, CPU paths)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Convergence / restart checks and the PDHG step. GPU blocks 742-781 and
//! 1087-1125 are not ported. `performPdhgStep`, `accumulateAverages` and
//! `prepareNextIteration` are compiled but not called by the CPU `solve`
//! (it always runs the Halpern step); they are ported for completeness and
//! covered by unit tests only.

use super::super::PDHG_CHECK_INTERVAL;
use super::super::defs::*;
use super::super::linalg;
use super::super::solver_results::*;
use super::PDLPSolver;
use crate::prelude::*;

impl<'a> PDLPSolver<'a> {
    /// HiGHS: pdhg.cc::PDLPSolver::runConvergenceCheck (lines 784-899)
    ///
    /// `computeAverageIterate(Ax_avg_, ATy_avg_)` (take / put back, rule
    /// (i)); with `use_halpern_restart && iter > 0` the `"[L]"` check runs
    /// on `x_next_` / `y_next_` with fresh local `Ax_pdhg` / `ATy_pdhg`,
    /// otherwise on `x_current_` / `y_current_` with the caches; both pass
    /// `dSlackPos_` / `dSlackNeg_` (`slacks_are_members = true`). The
    /// `"[A]"` check (only without Halpern) uses the averages and
    /// `dSlackPosAvg_` / `dSlackNegAvg_` (`false`). On convergence: copy the
    /// chosen iterate to `output_x` / `output_y` (and the average slacks
    /// into `dSlackPos_` / `dSlackNeg_` when the average won),
    /// `final_iter_count_ = iter`, `results_ =` the chosen results,
    /// `status = OPTIMAL`, return true. Otherwise `results_ =
    /// current_results`, return false. `iter` is `size_t` upstream.
    pub(crate) fn runConvergenceCheck(
        &mut self,
        iter: usize,
        output_x: &mut Vec<f64>,
        output_y: &mut Vec<f64>,
        status: &mut TerminationStatus,
    ) -> bool {
        // 1. Compute Average Iterate (GPU or CPU)
        // GPU branch dropped (lines 788-789).
        {
            let mut ax_avg = std::mem::take(&mut self.Ax_avg_);
            let mut aty_avg = std::mem::take(&mut self.ATy_avg_);
            self.computeAverageIterate(&mut ax_avg, &mut aty_avg);
            self.Ax_avg_ = ax_avg;
            self.ATy_avg_ = aty_avg;
        }

        let mut current_results = SolverResults::default();
        let mut average_results = SolverResults::default();
        let current_converged: bool;
        let mut average_converged = false;

        // 2. Check Convergence (GPU or CPU)
        // GPU branch dropped (lines 800-818).
        // CPU branch (lines 820-838):
        if self.params_.use_halpern_restart && iter > 0 {
            let mut Ax_pdhg = vec![0.0; self.lp_.num_row_ as usize];
            let mut ATy_pdhg = vec![0.0; self.lp_.num_col_ as usize];

            linalg::ax(&self.lp_, &self.x_next_, &mut Ax_pdhg);
            linalg::aTy(&self.lp_, &self.y_next_, &mut ATy_pdhg);

            let mut dsp = std::mem::take(&mut self.dSlackPos_);
            let mut dsn = std::mem::take(&mut self.dSlackNeg_);
            current_converged = self.checkConvergence(
                iter as HighsInt,
                &self.x_next_,
                &self.y_next_,
                &Ax_pdhg,
                &ATy_pdhg,
                self.params_.tolerance,
                &mut current_results,
                "[L]",
                &mut dsp,
                &mut dsn,
                true,
            );
            self.dSlackPos_ = dsp;
            self.dSlackNeg_ = dsn;
        } else {
            let mut dsp = std::mem::take(&mut self.dSlackPos_);
            let mut dsn = std::mem::take(&mut self.dSlackNeg_);
            current_converged = self.checkConvergence(
                iter as HighsInt,
                &self.x_current_,
                &self.y_current_,
                &self.Ax_cache_,
                &self.ATy_cache_,
                self.params_.tolerance,
                &mut current_results,
                "[L]",
                &mut dsp,
                &mut dsn,
                true,
            );
            self.dSlackPos_ = dsp;
            self.dSlackNeg_ = dsn;
        }

        if !self.params_.use_halpern_restart {
            let mut dsp_avg = std::mem::take(&mut self.dSlackPosAvg_);
            let mut dsn_avg = std::mem::take(&mut self.dSlackNegAvg_);
            average_converged = self.checkConvergence(
                iter as HighsInt,
                &self.x_avg_,
                &self.y_avg_,
                &self.Ax_avg_,
                &self.ATy_avg_,
                self.params_.tolerance,
                &mut average_results,
                "[A]",
                &mut dsp_avg,
                &mut dsn_avg,
                false,
            );
            self.dSlackPosAvg_ = dsp_avg;
            self.dSlackNegAvg_ = dsn_avg;
        }

        // 3. Handle Convergence Success
        if current_converged || average_converged {
            let prefer_avg = average_converged;

            // Copy result to output vectors
            // GPU branch dropped (lines 846-865).
            // CPU branch (lines 866-878):
            if prefer_avg {
                *output_x = self.x_avg_.clone();
                *output_y = self.y_avg_.clone();
                self.dSlackPos_ = self.dSlackPosAvg_.clone();
                self.dSlackNeg_ = self.dSlackNegAvg_.clone();
            } else if self.params_.use_halpern_restart && iter > 0 {
                *output_x = self.x_next_.clone();
                *output_y = self.y_next_.clone();
            } else {
                *output_x = self.x_current_.clone();
                *output_y = self.y_current_.clone();
            }

            self.final_iter_count_ = iter as HighsInt;
            self.results_ = if prefer_avg {
                average_results
            } else {
                current_results
            };

            // Final logging on optimality
            let forced = true;
            self.logger_.printIterationStats(
                self.final_iter_count_,
                &self.results_,
                self.primal_weight_,
                forced,
            );
            self.logger_.info(if prefer_avg {
                "Average solution converged"
            } else {
                "Current solution converged"
            });

            *status = TerminationStatus::OPTIMAL;
            return true; // Stop
        } else {
            // Log the iteration if not converged
            self.logger_.printIterationStats(
                iter as HighsInt,
                &current_results,
                self.primal_weight_,
                false,
            );
        }

        self.results_ = current_results;
        false
    }

    /// HiGHS: pdhg.cc::PDLPSolver::checkRestartCriteria (lines 901-927)
    ///
    /// `halpern_iteration` and `final_iter_count` are `int` upstream; the
    /// comparison with `getArtificialRestartThreshold() * final_iter_count`
    /// is in `double`.
    pub(crate) fn checkRestartCriteria(
        &self,
        current_fpe: f64,
        initial_fpe: f64,
        last_trial_fpe: f64,
        halpern_iteration: i32,
        final_iter_count: i32,
    ) -> bool {
        if final_iter_count == PDHG_CHECK_INTERVAL {
            return true;
        }

        if final_iter_count > PDHG_CHECK_INTERVAL {
            if current_fpe <= self.restart_scheme_.getSufficientDecayFactor() * initial_fpe {
                return true;
            }
            if current_fpe <= self.restart_scheme_.getNecessaryDecayFactor() * initial_fpe {
                if current_fpe > last_trial_fpe {
                    return true;
                }
            }
            if (halpern_iteration as f64)
                >= self.restart_scheme_.getArtificialRestartThreshold() * (final_iter_count as f64)
            {
                return true;
            }
        }

        false
    }

    /// HiGHS: pdhg.cc::PDLPSolver::performPdhgStep (lines 932-956)
    ///
    /// `PID` hits `assert(111 == 999)`: `debug_assert!(false)` and nothing.
    pub(crate) fn performPdhgStep(&mut self) {
        // #if PDLP_PROFILE dropped (lines 933-935, 953-955)

        // CPU/GPU agnostic logic for step strategy
        match self.params_.step_size_strategy {
            StepSizeStrategy::FIXED => {
                self.updateIteratesFixed();
            }
            StepSizeStrategy::ADAPTIVE => {
                self.updateIteratesAdaptive();
            }
            StepSizeStrategy::MALITSKY_POCK => {
                // Note: Error handling for MP failure simplified for brevity
                self.updateIteratesMalitskyPock(false);
            }
            StepSizeStrategy::PID => {
                debug_assert!(false);
            }
        }
    }

    /// HiGHS: pdhg.cc::PDLPSolver::performHalpernPdhgStep (lines 961-1018)
    ///
    /// `current_k = halpern_iteration_ + k_offset` (`int`), `w = current_k
    /// / (current_k + 1.0)`; the primal and dual projections, reflections
    /// and the anchored blend exactly as written; `is_major` additionally
    /// stores `x_next_`, `y_next_` and `halpern_dual_slack_next_` (resized
    /// to `num_col_` first if needed) and sets
    /// `halpern_dual_slack_next_valid_`. The dual clamp is `cpp_max(lower,
    /// cpp_min(temp, upper))`.
    pub(crate) fn performHalpernPdhgStep(&mut self, is_major: bool, k_offset: i32) {
        let primal_step = self.stepsize_.primal_step;
        let dual_step = self.stepsize_.dual_step;
        let rho = self.params_.halpern_gamma;

        let current_k = self.halpern_iteration_ + k_offset;
        let w = (current_k as f64) / (current_k as f64 + 1.0);

        // --- STEP 1: Primal projection + reflection ---
        if is_major {
            if self.halpern_dual_slack_next_.len() != self.lp_.num_col_ as usize {
                self.halpern_dual_slack_next_
                    .resize(self.lp_.num_col_ as usize, 0.0);
            }
            self.halpern_dual_slack_next_valid_ = true;
        }
        for i in 0..self.lp_.num_col_ {
            let idx = i as usize;
            let temp = self.x_current_[idx]
                - primal_step * (self.lp_.col_cost_[idx] - self.ATy_cache_[idx]);
            let temp_proj =
                linalg::projectBox(temp, self.lp_.col_lower_[idx], self.lp_.col_upper_[idx]);
            if is_major {
                self.x_next_[idx] = temp_proj; // pdhg_primal (for convergence check)
                self.halpern_dual_slack_next_[idx] = (temp_proj - temp) / primal_step;
            }
            self.reflected_x_[idx] = 2.0 * temp_proj - self.x_current_[idx];
        }

        // --- STEP 2: SpMV on reflected primal ---
        linalg::ax(&self.lp_, &self.reflected_x_, &mut self.Ax_next_);

        // --- STEP 3: Dual projection + reflection ---
        for j in 0..self.lp_.num_row_ {
            let idx = j as usize;
            let temp = self.y_current_[idx] / dual_step - self.Ax_next_[idx];
            let lower_proj = -self.lp_.row_upper_[idx];
            let upper_proj = -self.lp_.row_lower_[idx];
            let temp_proj = cpp_max(lower_proj, cpp_min(temp, upper_proj));

            let pdhg_dual_j = (temp - temp_proj) * dual_step;

            if is_major {
                self.y_next_[idx] = pdhg_dual_j;
            }
            self.reflected_y_[idx] = 2.0 * pdhg_dual_j - self.y_current_[idx];
        }

        // --- STEP 4: Halpern blend ---
        for i in 0..self.lp_.num_col_ {
            let idx = i as usize;
            let blended = rho * self.reflected_x_[idx] + (1.0 - rho) * self.x_current_[idx];
            self.x_current_[idx] = w * blended + (1.0 - w) * self.x_anchor_[idx];
        }
        for j in 0..self.lp_.num_row_ {
            let idx = j as usize;
            let blended = rho * self.reflected_y_[idx] + (1.0 - rho) * self.y_current_[idx];
            self.y_current_[idx] = w * blended + (1.0 - w) * self.y_anchor_[idx];
        }

        // --- STEP 5: Recompute aTy ---
        linalg::aTy(&self.lp_, &self.y_current_, &mut self.ATy_cache_);
    }

    /// HiGHS: pdhg.cc::PDLPSolver::accumulateAverages (lines 1023-1053)
    ///
    /// `inner_iter = iter - restart_scheme_.getLastRestartIter()` (size_t
    /// minus HighsInt, converted to HighsInt). `dMeanStepSize` is computed
    /// and unused upstream (kept).
    pub(crate) fn accumulateAverages(&mut self, iter: usize) {
        // #if PDLP_PROFILE dropped (lines 1024-1026, 1050-1052)

        // size_t minus HighsInt (converted to size_t, wrapping), then narrowed.
        let inner_iter =
            iter.wrapping_sub(self.restart_scheme_.getLastRestartIter() as usize) as HighsInt;
        let _d_mean_step_size = (self.stepsize_.primal_step * self.stepsize_.dual_step).sqrt();

        // GPU branch dropped (lines 1031-1041).
        // CPU branch (lines 1042-1048):
        if self.params_.use_halpern_restart {
            let x_cur = std::mem::take(&mut self.x_current_);
            let y_cur = std::mem::take(&mut self.y_current_);
            self.updateAverageIterates(&x_cur, &y_cur, inner_iter);
            self.x_current_ = x_cur;
            self.y_current_ = y_cur;
        } else {
            let x_next = std::mem::take(&mut self.x_next_);
            let y_next = std::mem::take(&mut self.y_next_);
            self.updateAverageIterates(&x_next, &y_next, inner_iter);
            self.x_next_ = x_next;
            self.y_next_ = y_next;
        }
    }

    /// HiGHS: pdhg.cc::PDLPSolver::prepareNextIteration (lines 1058-1085)
    #[allow(clippy::needless_return)] // C++ early return kept (line 1062)
    pub(crate) fn prepareNextIteration(&mut self) {
        if self.params_.use_halpern_restart {
            // Halpern update already wrote into x_current_, no swap needed
            // Just ensure caches are correct if not done in PerformHalpernStep
            return;
        } else {
            // Standard PDHG: x_current becomes x_next
            // GPU branch dropped (lines 1065-1077).
            // CPU branch (lines 1078-1083):
            self.x_current_.clone_from(&self.x_next_);
            self.y_current_.clone_from(&self.y_next_);
            self.Ax_cache_.clone_from(&self.Ax_next_);
            self.ATy_cache_.clone_from(&self.ATy_next_);
        }
    }

    /// HiGHS: pdhg.cc::PDLPSolver::solveReturn (lines 1128-1136)
    ///
    /// `results_.term_code = term_code`, then the oracle hook
    /// `hipdlp_end(self)` ([`crate::pdlp::oracle_trace`]).
    pub(crate) fn solveReturn(&mut self, term_code: TerminationStatus) {
        // GPU cleanup dropped (lines 1129-1131).
        self.results_.term_code = term_code;
        // #if PDLP_PROFILE dropped (lines 1133-1135)
        crate::pdlp::oracle_trace::hipdlp_end(self);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn restart_criteria() {
        let s = PDLPSolver::default();
        // First check interval always restarts.
        assert!(s.checkRestartCriteria(1.0, 1.0, 1.0, 40, PDHG_CHECK_INTERVAL));
        // Sufficient decay 0.2.
        assert!(s.checkRestartCriteria(0.1, 1.0, f64::INFINITY, 1, 80));
        // Necessary decay 0.8 and no progress since last trial.
        assert!(s.checkRestartCriteria(0.5, 1.0, 0.4, 1, 80));
        assert!(!s.checkRestartCriteria(0.5, 1.0, 0.6, 1, 80));
        // Artificial: 0.36 * 100 = 36 <= 40.
        assert!(s.checkRestartCriteria(0.9, 1.0, 0.1, 40, 100));
        assert!(!s.checkRestartCriteria(0.9, 1.0, 0.1, 40, 200));
    }

    #[test]
    fn solve_return_sets_term_code() {
        let mut s = PDLPSolver::default();
        s.solveReturn(TerminationStatus::MAXITER);
        assert_eq!(s.getTerminationCode(), TerminationStatus::MAXITER);
    }

    #[test]
    fn halpern_major_step_fills_next_and_slack() {
        let lp = super::super::tests::four_row_lp();
        let mut s = PDLPSolver::default();
        s.passLp(&lp);
        s.preprocessLp();
        s.initialize();
        s.stepsize_.primal_step = 0.5;
        s.stepsize_.dual_step = 0.5;
        s.params_.halpern_gamma = 1.0;
        s.performHalpernPdhgStep(true, 1);
        assert!(s.halpern_dual_slack_next_valid_);
        assert_eq!(s.halpern_dual_slack_next_.len(), s.lp_.num_col_ as usize);
        // x = 0 - 0.5 * cost = -0.5 projected to the lower bound 0.
        assert_eq!(s.x_next_[0], 0.0);
        let _ = linalg::norm2(&s.x_current_);
        let _ = StepSizeStrategy::FIXED;
    }

    fn solver(lp: &HighsLp) -> PDLPSolver<'_> {
        let mut s = PDLPSolver::default();
        s.passLp(lp);
        s.preprocessLp();
        s.initialize();
        s
    }

    #[test]
    fn halpern_weight_and_minor_step() {
        // w = k / (k + 1.0) with k = halpern_iteration_ + k_offset; a minor
        // step leaves x_next_ / y_next_ alone.
        let lp = super::super::tests::four_row_lp();
        let mut s = solver(&lp);
        s.stepsize_.primal_step = 0.5;
        s.stepsize_.dual_step = 0.5;
        s.params_.halpern_gamma = 1.0;
        s.halpern_iteration_ = 0;
        s.x_next_.iter_mut().for_each(|v| *v = 7.0);
        let anchor = s.x_anchor_.clone();
        // k = 0: w = 0, the blend returns the anchor.
        s.performHalpernPdhgStep(false, 0);
        assert!(!s.halpern_dual_slack_next_valid_);
        assert!(s.x_next_.iter().all(|&v| v == 7.0));
        assert_eq!(s.x_current_, anchor);
    }

    #[test]
    fn prepare_next_iteration_copies_only_without_halpern() {
        let lp = super::super::tests::four_row_lp();
        let mut s = solver(&lp);
        let n = s.lp_.num_col_ as usize;
        let m = s.lp_.num_row_ as usize;
        s.x_next_ = vec![1.0; n];
        s.y_next_ = vec![2.0; m];
        s.Ax_next_ = vec![3.0; m];
        s.ATy_next_ = vec![4.0; n];
        s.params_.use_halpern_restart = true;
        s.prepareNextIteration();
        assert!(s.x_current_.iter().all(|&v| v != 1.0));
        s.params_.use_halpern_restart = false;
        s.prepareNextIteration();
        assert_eq!(s.x_current_, vec![1.0; n]);
        assert_eq!(s.y_current_, vec![2.0; m]);
        assert_eq!(s.Ax_cache_, vec![3.0; m]);
        assert_eq!(s.ATy_cache_, vec![4.0; n]);
    }

    #[test]
    fn accumulate_averages_picks_iterate_by_halpern_flag() {
        let lp = super::super::tests::four_row_lp();
        let mut s = solver(&lp);
        let n = s.lp_.num_col_ as usize;
        s.stepsize_.primal_step = 4.0;
        s.stepsize_.dual_step = 1.0;
        s.x_current_ = vec![1.0; n];
        s.x_next_ = vec![10.0; n];
        s.x_sum_ = vec![0.0; n];
        s.sum_weights_ = 0.0;
        s.params_.use_halpern_restart = true;
        s.accumulateAverages(5);
        assert_eq!(s.x_sum_, vec![2.0; n]);
        assert_eq!(s.sum_weights_, 2.0);
        s.params_.use_halpern_restart = false;
        s.accumulateAverages(6);
        assert_eq!(s.x_sum_, vec![22.0; n]);
        assert_eq!(s.sum_weights_, 4.0);
    }

    #[test]
    fn pdhg_step_fixed_fills_next() {
        let lp = super::super::tests::four_row_lp();
        let mut s = solver(&lp);
        s.params_.step_size_strategy = StepSizeStrategy::FIXED;
        s.stepsize_.primal_step = 0.5;
        s.stepsize_.dual_step = 0.5;
        s.performPdhgStep();
        assert_eq!(s.x_next_.len(), s.lp_.num_col_ as usize);
        assert_eq!(s.y_next_.len(), s.lp_.num_row_ as usize);
        let mut ax = vec![0.0; s.lp_.num_row_ as usize];
        linalg::ax(&s.lp_, &s.x_next_, &mut ax);
        assert_eq!(ax, s.Ax_next_);
    }

    #[test]
    fn convergence_check_not_converged_keeps_current_results() {
        // Zero point of four_row_lp violates row 0 (>= 1): no convergence on
        // either path; results_ is the "[L]" result and status is untouched.
        for halpern in [false, true] {
            let lp = super::super::tests::four_row_lp();
            let mut s = solver(&lp);
            s.params_.use_halpern_restart = halpern;
            s.final_iter_count_ = -1;
            let mut out_x = vec![9.0];
            let mut out_y = vec![9.0];
            let mut status = TerminationStatus::NOTSET;
            let stop = s.runConvergenceCheck(3, &mut out_x, &mut out_y, &mut status);
            assert!(!stop);
            assert_eq!(status, TerminationStatus::NOTSET);
            assert_eq!(s.final_iter_count_, -1);
            assert_eq!((out_x, out_y), (vec![9.0], vec![9.0]));
            assert!(s.results_.primal_feasibility > 0.0);
        }
    }
}
