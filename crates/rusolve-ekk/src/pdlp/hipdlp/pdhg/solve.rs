//! HiGHS source: pdlp/hipdlp/pdhg.cc (lines 494-739)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiPDLP main loop (CPU, Halpern restarts) and the fixed-point error.

use super::super::PDHG_CHECK_INTERVAL;
use super::super::defs::*;
use super::super::linalg;
use super::super::solver_results::*;
use super::PDLPSolver;
use crate::prelude::*;

impl<'a> PDLPSolver<'a> {
    /// HiGHS: pdhg.cc::PDLPSolver::solve (lines 494-707)
    ///
    /// CPU path (the `#ifdef CUPDLP_GPU` blocks 523-548, 584-595, 610-627,
    /// 633-637, 668-680, 693-699 are dropped):
    /// `initializeStepSizes`; `initialize`; `best_primal_weight_ =
    /// primal_weight_`; `best_primal_dual_residual_gap_ = inf`;
    /// `restart_scheme_.initialize(results_, &params_)`; `x_current_ = x`,
    /// `y_current_ = y`, `projectBounds`; Halpern anchors and
    /// `halpern_iteration_ = 0`; `Ax_cache_` / `ATy_cache_`;
    /// `final_iter_count_ = 0`; `runConvergenceCheck(0, x, y, status)`
    /// (return through `solveReturn` on convergence). Then while
    /// `final_iter_count_ < params_.max_iterations`:
    /// time check `highs_timer_p_.read() > params_.time_limit` -> TIMEOUT;
    /// `performHalpernPdhgStep(true, 1)`; if `do_restart` then `fpe_ =
    /// computeFixedPointError()`, `initial_fpe_ = fpe_`, `do_restart =
    /// false`; `performHalpernPdhgStep(false, i)` for `i` in
    /// `2..=PDHG_CHECK_INTERVAL - 1`; `performHalpernPdhgStep(true,
    /// PDHG_CHECK_INTERVAL)`; `fpe_ = computeFixedPointError()`;
    /// `halpern_iteration_ += 40`, `final_iter_count_ += 40`;
    /// `runConvergenceCheck` (break on convergence); `do_restart =
    /// checkRestartCriteria(fpe_, initial_fpe_, last_trial_fpe,
    /// halpern_iteration_, final_iter_count_)`; `last_trial_fpe = fpe_`; on
    /// restart: PID -> `updatePrimalWeightAtRestart(&results_)` (take a
    /// copy), anchors and current := next, recompute the caches,
    /// `halpern_iteration_ = 0`, `last_trial_fpe = inf`. After the loop
    /// NOTSET -> MAXITER; `solveReturn`.
    ///
    /// Note: `x` / `y` are written only by a converging
    /// `runConvergenceCheck`; on MAXITER / TIMEOUT they keep their input.
    ///
    /// Oracle hooks ([`crate::pdlp::oracle_trace`]): `hipdlp_check(self,
    /// stop)` right after each of the two `runConvergenceCheck` calls
    /// (`stop` = its return value); `hipdlp_restart(self)` at the end of the
    /// restart block (after `last_trial_fpe = inf`); `hipdlp_end` is in
    /// `solveReturn`.
    pub fn solve(&mut self, x: &mut Vec<f64>, y: &mut Vec<f64>) {
        // 1. Initialization and setup
        // #ifdef CUPDLP_GPU setupGpu() dropped (line 501)
        self.initializeStepSizes();
        self.initialize(); // Resets vectors, caches, and sets initial x_current, y_current
        self.best_primal_weight_ = self.primal_weight_;
        self.best_primal_dual_residual_gap_ = f64::INFINITY;

        // initialize internal solver state
        // C++: restart_scheme_.passParams(&params_); restart_scheme_.initialize(results_);
        // Rust: passParams is a no-op; initialize takes params as a parameter.
        self.restart_scheme_
            .initialize(&self.results_, &self.params_);

        // Copy initial input x,y to internal state
        self.x_current_ = x.clone();
        self.y_current_ = y.clone();
        linalg::projectBounds(&self.lp_, &mut self.x_current_);

        if self.params_.use_halpern_restart {
            self.x_anchor_ = self.x_current_.clone();
            self.y_anchor_ = self.y_current_.clone();
            self.halpern_iteration_ = 0;
        }

        // #ifdef CUPDLP_GPU blocks 524-548 dropped
        // #else (line 549-552):
        linalg::ax(&self.lp_, &self.x_current_, &mut self.Ax_cache_);
        linalg::aTy(&self.lp_, &self.y_current_, &mut self.ATy_cache_);

        let mut termination_status = TerminationStatus::NOTSET;
        let mut do_restart = false;
        let mut last_trial_fpe = f64::INFINITY;
        self.final_iter_count_ = 0;

        // =========================================================================
        // INITIAL CONVERGENCE CHECK (ITERATION 0)
        // =========================================================================
        let should_terminate = self.runConvergenceCheck(
            self.final_iter_count_ as usize,
            x,
            y,
            &mut termination_status,
        );
        crate::pdlp::oracle_trace::hipdlp_check(self, should_terminate);

        if should_terminate {
            self.solveReturn(termination_status);
            return;
        }

        // #ifdef CUPDLP_GPU graph_created / graphExec declarations dropped (lines 570-573)

        // 2. Main cuPDLPx-style Loop
        while self.final_iter_count_ < self.params_.max_iterations {
            // check global time limit
            if self.highs_timer_p_.unwrap().read() > self.params_.time_limit {
                self.logger_.info("Time limit reached.");
                termination_status = TerminationStatus::TIMEOUT;
                break;
            }

            // -- Step 1 (Major, isolated for restart-FPE check) --
            // #ifdef CUPDLP_GPU blocks 585-595 dropped
            // #else (line 597):
            self.performHalpernPdhgStep(true, 1);

            if do_restart {
                // #ifdef CUPDLP_GPU computeFixedPointErrorGpu dropped (line 602)
                // #else (line 604):
                self.fpe_ = self.computeFixedPointError();
                self.initial_fpe_ = self.fpe_;
                #[allow(unused_assignments)]
                {
                    do_restart = false;
                }
            }
            // -- Steps 2 to PDHG_CHECK_INTERVAL - 1 (Minor) --
            // #ifdef CUPDLP_GPU graph capture blocks 610-627 dropped
            // #else (lines 629-632):
            for i in 2..=(PDHG_CHECK_INTERVAL - 1) {
                self.performHalpernPdhgStep(false, i);
            }
            self.performHalpernPdhgStep(true, PDHG_CHECK_INTERVAL);

            // Compute Error for Restart Check
            // #ifdef CUPDLP_GPU computeFixedPointErrorGpu dropped (line 637)
            // #else (line 639):
            self.fpe_ = self.computeFixedPointError();
            self.halpern_iteration_ += PDHG_CHECK_INTERVAL;
            self.final_iter_count_ += PDHG_CHECK_INTERVAL;

            // =======================================================================
            // CONVERGENCE CHECK
            // =======================================================================
            let mut check_status = TerminationStatus::NOTSET;
            let should_terminate =
                self.runConvergenceCheck(self.final_iter_count_ as usize, x, y, &mut check_status);
            crate::pdlp::oracle_trace::hipdlp_check(self, should_terminate);
            if should_terminate {
                termination_status = check_status;
                break;
            }

            // =======================================================================
            // RESTART CHECK
            // =======================================================================
            do_restart = self.checkRestartCriteria(
                self.fpe_,
                self.initial_fpe_,
                last_trial_fpe,
                self.halpern_iteration_,
                self.final_iter_count_,
            );
            last_trial_fpe = self.fpe_;

            // -- Perform Restart --
            if do_restart {
                if self.params_.step_size_strategy == StepSizeStrategy::PID {
                    let results_copy = self.results_.clone();
                    self.updatePrimalWeightAtRestart(&results_copy);
                }
                // #ifdef CUPDLP_GPU blocks 667-679 dropped
                // #else (lines 681-686):
                self.x_anchor_ = self.x_next_.clone();
                self.y_anchor_ = self.y_next_.clone();
                self.x_current_ = self.x_next_.clone();
                self.y_current_ = self.y_next_.clone();
                linalg::ax(&self.lp_, &self.x_current_, &mut self.Ax_cache_);
                linalg::aTy(&self.lp_, &self.y_current_, &mut self.ATy_cache_);

                self.halpern_iteration_ = 0;
                last_trial_fpe = f64::INFINITY;
                crate::pdlp::oracle_trace::hipdlp_restart(self);
            }
        }

        // #ifdef CUPDLP_GPU graphExec destroy dropped (lines 694-698)

        // 3. Loop Finished
        if termination_status == TerminationStatus::NOTSET {
            self.logger_
                .info("Iteration limit reached without convergence");
            termination_status = TerminationStatus::MAXITER;
        }

        self.solveReturn(termination_status);
    }

    /// HiGHS: pdhg.cc::PDLPSolver::computeFixedPointError (lines 709-739)
    ///
    /// `delta = next - reflected`; `sqrt(max(0, |dx|^2 * omega + |dy|^2 /
    /// omega + 2 * eta * dx'A'dy))` with `params_.omega` / `params_.eta`.
    pub(crate) fn computeFixedPointError(&self) -> f64 {
        let mut primal_norm_sq: f64 = 0.0;
        let mut dual_norm_sq: f64 = 0.0;
        let mut cross_term: f64 = 0.0;

        let num_col = self.lp_.num_col_ as usize;
        let num_row = self.lp_.num_row_ as usize;

        let mut delta_x = vec![0.0f64; num_col];
        let mut delta_y = vec![0.0f64; num_row];

        for i in 0..num_col {
            delta_x[i] = self.x_next_[i] - self.reflected_x_[i];
            primal_norm_sq += delta_x[i] * delta_x[i];
        }

        for i in 0..num_row {
            delta_y[i] = self.y_next_[i] - self.reflected_y_[i];
            dual_norm_sq += delta_y[i] * delta_y[i];
        }

        let mut AT_delta_y = vec![0.0f64; num_col];
        linalg::aTy(&self.lp_, &delta_y, &mut AT_delta_y);

        for i in 0..num_col {
            cross_term += delta_x[i] * AT_delta_y[i];
        }

        let movement = primal_norm_sq * self.params_.omega + dual_norm_sq / self.params_.omega;
        let interaction = 2.0 * self.params_.eta * cross_term;

        // std::max(0.0, ...) -> cpp_max (first argument wins on ties and NaN)
        cpp_max(0.0, movement + interaction).sqrt()
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests. The loop itself is covered by
    //! the oracle test `tests/t_w3_hipdlp.rs`.
    use super::*;

    #[test]
    fn fixed_point_error_of_identical_points_is_zero() {
        let lp = super::super::tests::four_row_lp();
        let mut s = PDLPSolver::default();
        s.passLp(&lp);
        s.preprocessLp();
        s.initialize();
        s.params_.omega = 1.0;
        s.params_.eta = 1.0;
        assert_eq!(s.computeFixedPointError(), 0.0);
        s.x_next_[0] = 3.0;
        s.reflected_x_[0] = 1.0;
        // |dx|^2 * omega = 4, no dual movement.
        assert_eq!(s.computeFixedPointError(), 2.0);
    }

    #[test]
    fn zero_iteration_limit_is_maxiter() {
        let lp = super::super::tests::four_row_lp();
        let timer = HighsTimer::default();
        let options = HighsOptions::default();
        let mut s = PDLPSolver::default();
        s.setup(&options, &timer);
        s.params_.max_iterations = 0;
        s.passLp(&lp);
        s.preprocessLp();
        s.scaleProblem();
        let mut x = vec![0.0; s.getNumCol() as usize];
        let mut y = vec![0.0; s.getNumRow() as usize];
        s.solve(&mut x, &mut y);
        assert_eq!(s.getTerminationCode(), TerminationStatus::MAXITER);
        assert_eq!(s.getIterationCount(), 0);
        let _ = (linalg::norm2(&x), PDHG_CHECK_INTERVAL);
    }

    /// pdhg.cc 578-581: the time check runs before any PDHG step, so a
    /// negative time limit stops at iteration 0 with TIMEOUT and leaves
    /// `x` / `y` at their input.
    #[test]
    fn negative_time_limit_is_timeout_before_any_step() {
        let lp = super::super::tests::four_row_lp();
        let timer = HighsTimer::default();
        let options = HighsOptions::default();
        let mut s = PDLPSolver::default();
        s.setup(&options, &timer);
        s.params_.max_iterations = 1000;
        s.params_.time_limit = -1.0;
        s.passLp(&lp);
        s.preprocessLp();
        s.scaleProblem();
        let mut x = vec![0.0; s.getNumCol() as usize];
        let mut y = vec![0.0; s.getNumRow() as usize];
        s.solve(&mut x, &mut y);
        assert_eq!(s.getTerminationCode(), TerminationStatus::TIMEOUT);
        assert_eq!(s.getIterationCount(), 0);
        assert!(x.iter().all(|&v| v == 0.0) && y.iter().all(|&v| v == 0.0));
    }

    /// pdhg.cc 576-692: one pass of the loop advances `final_iter_count_`
    /// by PDHG_CHECK_INTERVAL; with `max_iterations = PDHG_CHECK_INTERVAL`
    /// the loop runs exactly once, and the first in-loop restart check
    /// (`final_iter_count == PDHG_CHECK_INTERVAL`) always restarts, which
    /// resets `halpern_iteration_` to 0 and copies x_next_ into the anchor.
    #[test]
    fn one_loop_pass_counts_check_interval_and_restarts() {
        let lp = super::super::tests::four_row_lp();
        let timer = HighsTimer::default();
        let options = HighsOptions::default();
        let mut s = PDLPSolver::default();
        s.setup(&options, &timer);
        s.params_.max_iterations = PDHG_CHECK_INTERVAL;
        s.passLp(&lp);
        s.preprocessLp();
        s.scaleProblem();
        let mut x = vec![0.0; s.getNumCol() as usize];
        let mut y = vec![0.0; s.getNumRow() as usize];
        s.solve(&mut x, &mut y);
        assert_eq!(s.getIterationCount(), PDHG_CHECK_INTERVAL);
        let status = s.getTerminationCode();
        if status == TerminationStatus::MAXITER {
            if s.params_.use_halpern_restart {
                assert_eq!(s.halpern_iteration_, 0);
                assert_eq!(s.x_anchor_, s.x_next_);
                assert_eq!(s.y_anchor_, s.y_next_);
            }
            assert_eq!(s.x_current_, s.x_next_);
            assert_eq!(s.y_current_, s.y_next_);
        } else {
            assert_eq!(status, TerminationStatus::OPTIMAL);
        }
        assert!(s.fpe_.is_finite() && s.fpe_ >= 0.0);
    }
}
