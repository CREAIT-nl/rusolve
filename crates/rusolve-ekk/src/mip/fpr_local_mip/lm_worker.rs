//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! LocalMipWorker (local_mip_worker.h/cpp), perturb_solution.

use super::context::LeastViolatedPoint;
use super::lm_search::{infeasible_step, select_lift_move};
use super::rng::{
    S5Rng, is_integer, s5_bernoulli, s5_uniform_int_i64, s5_uniform_int_u64, s5_uniform_real,
};
use super::sink::S5Sink;
use super::*;
use crate::prelude::*;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// mip-heuristics (36f7d9d): local_mip_worker.h:18 (kPerturbBinaryFraction)
pub(crate) const kPerturbBinaryFraction: f64 = 0.2;

/// mip-heuristics (36f7d9d): local_mip_worker.cpp::perturb_solution (lines 20-80)
pub(crate) fn perturb_solution(
    solution: &mut [f64],
    binary: &[u8],
    integrality: &[HighsVarType],
    col_lb: &[f64],
    col_ub: &[f64],
    ncol: HighsInt,
    rng: &mut S5Rng,
) {
    for j in 0..ncol {
        if !is_integer(integrality, j) {
            continue;
        }
        if s5_uniform_real(rng, 0.0, 1.0) > kPerturbBinaryFraction {
            continue;
        }
        let j = j as usize;
        if binary[j] != 0u8 {
            solution[j] = if solution[j] < 0.5 { 1.0 } else { 0.0 };
        } else {
            if !solution[j].is_finite() {
                continue;
            }
            let mut lo = col_lb[j].ceil();
            let mut hi = col_ub[j].floor();
            let current = solution[j].round();
            if !lo.is_finite() || !hi.is_finite() || hi - lo > K_SAFE_INT64_DOUBLE_RANGE {
                lo = current - K_INF_BOUND_SHIFT_WINDOW;
                hi = current + K_INF_BOUND_SHIFT_WINDOW;
            }
            if hi <= lo {
                continue;
            }
            let irange = (hi - lo) as i64;
            if irange < 1 {
                continue;
            }
            let shift = s5_uniform_int_i64(rng, 1, irange);
            solution[j] = lo + (current - lo + (shift as f64)) % ((irange as f64) + 1.0);
            solution[j] = cpp_max(col_lb[j], cpp_min(col_ub[j], solution[j]));
        }
    }
}

/// mip-heuristics (36f7d9d): local_mip_worker.h::LocalMipWorker (lines 34-116)
///
/// `SolutionSink& sink_` (local_mip_worker.h:75) is not a field: the slot
/// that owns this worker across attempts (local_mip.cpp LmState) also hands
/// the sink to the runner between attempts, so the sink is passed to
/// `run_attempt` per call instead (lib.rs ownership rule). Same sink, same
/// calls, same order.
pub(crate) struct LocalMipWorker<'a, 'l> {
    exec_: ExecutionContext,
    least_violated: Option<&'l mut LeastViolatedPoint>,
    rng_: S5Rng,
    base_: WorkerBudgetState,
    trace_: WorkerTrace,
    ctx_: WorkerCtx<'a>,
    costed_vars_: Vec<HighsInt>,
    binary_vars_: Vec<HighsInt>,
    best_feasible_: bool,
    best_objective_: f64,
    best_solution_: Vec<f64>,
    steps_since_improvement_: HighsInt,
    restart_count_: HighsInt,
    step_: HighsInt,
    feasible_random_walks_done_: HighsInt,
    retained_violation_: f64,
    retained_: Vec<f64>,
    retained_unpublished_: bool,
}

impl<'a, 'l> LocalMipWorker<'a, 'l> {
    /// mip-heuristics (36f7d9d): local_mip_worker.cpp::LocalMipWorker::LocalMipWorker (lines 82-128)
    pub(crate) fn new(
        problem: &'a S5Problem,
        exec: ExecutionContext,
        total_budget: usize,
        stale_budget: usize,
        seed: u32,
        initial_solution: Option<&[f64]>,
        trace: WorkerTrace,
        least_violated: Option<&'l mut LeastViolatedPoint>,
    ) -> Self {
        let track_violation = least_violated.is_some();
        let mut ctx_ = WorkerCtx::new(problem);
        let mut base_ = WorkerBudgetState::default();
        base_.total_budget = total_budget;
        base_.stale_budget = stale_budget;
        // heuristic_context.h:228 (problem.degenerate())
        base_.finished = problem.ncol == 0 || problem.nrow == 0;
        let ncol = ctx_.ncol;

        let mut costed_vars_ = Vec::new();
        let mut binary_vars_ = Vec::new();
        for j in 0..ncol {
            if ctx_.col_cost[j as usize].abs() >= kEpsZero {
                costed_vars_.push(j);
            }
            if ctx_.is_binary(j) {
                binary_vars_.push(j);
            }
        }
        ctx_.lift.costed_vars = Some(Rc::new(costed_vars_.clone()));

        if let Some(src) = initial_solution {
            for j in 0..ncol {
                let j = j as usize;
                let mut v = src[j];
                if ctx_.is_int(j as HighsInt) {
                    v = v.round();
                }
                ctx_.solution[j] = cpp_max(ctx_.col_lb[j], cpp_min(ctx_.col_ub[j], v));
            }
        } else {
            for j in 0..ncol {
                let j = j as usize;
                let lb = ctx_.col_lb[j];
                let ub = ctx_.col_ub[j];
                let mut v = if 0.0 < lb {
                    lb
                } else if ub < 0.0 {
                    ub
                } else {
                    0.0
                };
                if ctx_.is_int(j as HighsInt) {
                    v = v.round();
                }
                ctx_.solution[j] = v;
            }
        }

        ctx_.track_violation = track_violation;
        // rusolve GUB swap (`gub_swap.rs`): off unless mip_ls_gub_swap.
        if exec.gub_swap {
            ctx_.gub = super::gub_swap::GubSwap::build(&ctx_);
        }
        ctx_.rebuild_state();

        LocalMipWorker {
            exec_: exec,
            least_violated,
            rng_: S5Rng::new(seed as u64),
            base_,
            trace_: trace,
            ctx_,
            costed_vars_,
            binary_vars_,
            best_feasible_: false,
            best_objective_: f64::INFINITY,
            best_solution_: vec![0.0; ncol as usize],
            steps_since_improvement_: 0,
            restart_count_: 0,
            step_: 0,
            feasible_random_walks_done_: 0,
            retained_violation_: f64::INFINITY,
            retained_: Vec::new(),
            retained_unpublished_: false,
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_worker.h::LocalMipWorker::finished (line 49)
    pub(crate) fn finished(&self) -> bool {
        self.base_.finished
    }

    /// mip-heuristics (36f7d9d): local_mip_worker.h::LocalMipWorker::traced_effort (line 56)
    pub(crate) fn traced_effort(&self) -> usize {
        self.trace_.at(self.ctx_.effort)
    }

    /// mip-heuristics (36f7d9d): local_mip_worker.h::LocalMipWorker::randomize_weights (lines 65, 130-139)
    pub(crate) fn randomize_weights(&mut self, rng: &mut S5Rng) {
        if !s5_bernoulli(rng, 0.5) {
            return;
        }
        for w in self.ctx_.weight.iter_mut() {
            *w = s5_uniform_int_u64(rng, 1, 2);
        }
        self.ctx_.obj_weight = s5_uniform_int_u64(rng, 1, 2);
    }

    /// mip-heuristics (36f7d9d): local_mip_worker.h::LocalMipWorker::weights_for_test (line 67)
    pub(crate) fn weights_for_test(&self) -> &[u64] {
        &self.ctx_.weight
    }

    /// mip-heuristics (36f7d9d): local_mip_worker.h::LocalMipWorker::obj_weight_for_test (line 68)
    pub(crate) fn obj_weight_for_test(&self) -> u64 {
        self.ctx_.obj_weight
    }

    /// mip-heuristics (36f7d9d): local_mip_worker.cpp::LocalMipWorker::run_attempt (lines 146-502)
    ///
    /// kInstrumented blocks dropped (local_mip_worker.cpp lines 191-196, 202-205,
    /// 212-215, 481-487 at 36f7d9d; production build has kInstrumented == false,
    /// local_mip.h:45).
    pub(crate) fn run_attempt(
        &mut self,
        attempt_budget: usize,
        sink: &mut S5Sink<'_>,
    ) -> AttemptResult {
        if self.base_.finished {
            return AttemptResult::default();
        }

        let ncol = self.ctx_.ncol;

        let mut attempt = AttemptResult::default();
        let effort_start = self.ctx_.effort;
        let mut effort_at_last_improvement = effort_start;
        let (evals_start, commits_start) = (self.ctx_.ls_evals, self.ctx_.ls_commits);

        let mut next_deadline_poll: usize = 0;
        let mut steps_since_poll: usize = 0;

        loop {
            let spent = self.ctx_.effort - effort_start;
            if spent >= attempt_budget || self.base_.exhausted_with(spent) {
                break;
            }

            self.retain_if_least_violated();

            let spent = self.ctx_.effort - effort_start;
            if self.base_.stale_with(spent) {
                self.base_.finished = true;
                break;
            }

            let spent = self.ctx_.effort - effort_start;
            let poll_now = if spent >= next_deadline_poll {
                true
            } else {
                steps_since_poll += 1;
                steps_since_poll > kTermCheckWork
            };
            if poll_now {
                self.publish_retained();
                if self.exec_.past_deadline() {
                    break;
                }
                next_deadline_poll = (self.ctx_.effort - effort_start) + kTermCheckWork;
                steps_since_poll = 0;
            }

            let feasible_mode = self.ctx_.violated.empty();

            if feasible_mode {
                let need_full_recheck = self.ctx_.was_infeasible
                    || (self.ctx_.feasible_recheck_counter % kFeasibleRecheckPeriod == 0);
                self.ctx_.was_infeasible = false;
                self.ctx_.feasible_recheck_counter += 1;

                let mut truly_feasible = true;
                if need_full_recheck {
                    truly_feasible = self.ctx_.full_recheck(true, false);
                }
                if !truly_feasible {
                    self.step_ += 1;
                    continue;
                }

                let obj = self.ctx_.current_obj;
                let improved = if !self.best_feasible_ {
                    true
                } else {
                    obj < self.best_objective_ - self.ctx_.epsilon
                };

                if improved {
                    if !need_full_recheck {
                        if !self.ctx_.full_recheck(false, true) {
                            self.ctx_.rebuild_state();
                            self.step_ += 1;
                            continue;
                        }
                    }
                    self.best_feasible_ = true;
                    self.best_objective_ = obj;
                    self.best_solution_ = self.ctx_.solution.clone();
                    self.steps_since_improvement_ = 0;
                    self.feasible_random_walks_done_ = 0;

                    let trace = self.trace_;
                    let effort_at = trace.at(self.ctx_.effort);
                    let result = sink.offer(obj, &self.ctx_.solution, &trace, effort_at);
                    if result.improved_incumbent {
                        attempt.found_improvement = true;
                        self.base_.reset_staleness();
                        effort_at_last_improvement = self.ctx_.effort;
                    }
                }

                let mut lift = std::mem::take(&mut self.ctx_.lift);
                lift.recompute_all(&mut self.ctx_);
                self.ctx_.lift = lift;
                let lift_best = select_lift_move(&mut self.ctx_);

                if lift_best.var_idx != -1 {
                    self.ctx_.apply_move_with_tabu(
                        lift_best.var_idx,
                        lift_best.new_val,
                        self.step_,
                        &mut self.rng_,
                    );
                    self.steps_since_improvement_ += 1;
                } else {
                    let cand = infeasible_step(
                        &mut self.ctx_,
                        &mut self.rng_,
                        self.step_,
                        self.best_feasible_,
                        self.best_objective_,
                        &self.costed_vars_,
                        &self.binary_vars_,
                    );
                    if cand.var_idx != -1 {
                        self.ctx_.apply_move_with_tabu(
                            cand.var_idx,
                            cand.new_val,
                            self.step_,
                            &mut self.rng_,
                        );
                        self.steps_since_improvement_ += 1;
                        if !self.ctx_.violated.empty() {
                            self.ctx_.lift.mark_all_dirty();
                        }
                    }
                }

                if self.steps_since_improvement_ >= kFeasiblePlateau {
                    if self.feasible_random_walks_done_ < kFeasibleMaxRandomWalks {
                        perturb_solution(
                            &mut self.ctx_.solution,
                            self.ctx_.binary,
                            self.ctx_.integrality,
                            self.ctx_.col_lb,
                            self.ctx_.col_ub,
                            self.ctx_.ncol,
                            &mut self.rng_,
                        );
                        self.ctx_.reset_weights();
                        self.ctx_.rebuild_state();
                        self.feasible_random_walks_done_ += 1;
                        self.steps_since_improvement_ = 0;
                    } else {
                        self.base_.finished = true;
                        break;
                    }
                }
            } else {
                let cand = infeasible_step(
                    &mut self.ctx_,
                    &mut self.rng_,
                    self.step_,
                    self.best_feasible_,
                    self.best_objective_,
                    &self.costed_vars_,
                    &self.binary_vars_,
                );

                if cand.var_idx != -1 {
                    self.ctx_.apply_move_with_tabu(
                        cand.var_idx,
                        cand.new_val,
                        self.step_,
                        &mut self.rng_,
                    );
                    self.steps_since_improvement_ += 1;
                    if self.ctx_.violated.empty() {
                        self.steps_since_improvement_ = 0;
                    }
                }
            }

            if self.step_ % kActivityPeriod == 0 && self.step_ > 0 {
                self.ctx_.rebuild_state();
            }

            if self.steps_since_improvement_ >= kRestartInterval {
                self.steps_since_improvement_ = 0;
                self.restart_count_ += 1;

                if self.best_feasible_ && (self.restart_count_ % 2 == 1) {
                    self.ctx_.solution = self.best_solution_.clone();
                } else {
                    for j in 0..ncol {
                        let j = j as usize;
                        if self.ctx_.is_binary(j as HighsInt) {
                            self.ctx_.solution[j] = if self.rng_.next_u64().is_multiple_of(2) {
                                0.0
                            } else {
                                1.0
                            };
                        } else if self.ctx_.is_int(j as HighsInt) {
                            let lo = cpp_max(self.ctx_.col_lb[j], -1e8);
                            let hi = cpp_min(self.ctx_.col_ub[j], lo + 100.0);
                            self.ctx_.solution[j] = cpp_max(
                                self.ctx_.col_lb[j],
                                cpp_min(
                                    self.ctx_.col_ub[j],
                                    s5_uniform_real(&mut self.rng_, lo, hi).round(),
                                ),
                            );
                        } else {
                            let lo = if self.ctx_.col_lb[j] > -kHighsInf {
                                self.ctx_.col_lb[j]
                            } else {
                                -1e6
                            };
                            let hi = if self.ctx_.col_ub[j] < kHighsInf {
                                self.ctx_.col_ub[j]
                            } else {
                                lo + 1e6
                            };
                            if hi > lo {
                                self.ctx_.solution[j] = s5_uniform_real(&mut self.rng_, lo, hi);
                            } else {
                                self.ctx_.solution[j] = lo;
                            }
                        }
                    }
                }

                self.ctx_.rebuild_state();
                self.ctx_.tabu_inc_until.fill(0);
                self.ctx_.tabu_dec_until.fill(0);
            }

            self.step_ += 1;
        }

        self.retain_if_least_violated();
        self.publish_retained();

        let attempt_effort = self.ctx_.effort - effort_start;
        self.base_.total_effort += attempt_effort;
        self.base_.effort_since_improvement += self.ctx_.effort - effort_at_last_improvement;
        if self.base_.exhausted() || self.base_.stale() {
            self.base_.finished = true;
        }
        attempt.effort = attempt_effort;
        super::stats::LOCAL_MIP.add_moves(
            self.ctx_.ls_evals - evals_start,
            self.ctx_.ls_commits - commits_start,
        );

        attempt
    }

    /// mip-heuristics (36f7d9d): local_mip_worker.cpp::LocalMipWorker::retain_if_least_violated (lines 509-520)
    fn retain_if_least_violated(&mut self) {
        let Some(holder) = self.least_violated.as_deref() else {
            return;
        };
        if self.ctx_.violated.elements.is_empty() {
            return;
        }
        let violation = self.ctx_.total_violation;
        if violation < self.retained_violation_ && holder.would_take(violation) {
            self.retained_violation_ = violation;
            self.retained_ = self.ctx_.solution.clone();
            self.retained_unpublished_ = true;
        }
    }

    /// mip-heuristics (36f7d9d): local_mip_worker.cpp::LocalMipWorker::publish_retained (lines 522-527)
    fn publish_retained(&mut self) {
        if self.retained_unpublished_ {
            if let Some(holder) = self.least_violated.as_deref_mut() {
                let _ = holder.offer(self.retained_violation_, &self.retained_);
            }
            self.retained_unpublished_ = false;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_problem_2x2_binary() -> S5Problem {
        S5Problem {
            col_lower: vec![0.0, 0.0],
            col_upper: vec![1.0, 1.0],
            row_lower: vec![-kHighsInf, 1.0],
            row_upper: vec![2.0, kHighsInf],
            col_cost: vec![1.0, 1.0],
            offset: 0.0,
            integrality: vec![HighsVarType::kInteger, HighsVarType::kInteger],
            ar_start: vec![0, 2, 4],
            ar_index: vec![0, 1, 0, 1],
            ar_value: vec![1.0, 1.0, 1.0, 1.0],
            csc: CscMatrix {
                col_start: vec![0, 2, 4],
                col_row: vec![0, 1, 0, 1],
                col_val: vec![1.0, 1.0, 1.0, 1.0],
            },
            uplocks: vec![0, 0],
            downlocks: vec![0, 0],
            feastol: 1e-6,
            epsilon: 1e-6,
            ncol: 2,
            nrow: 2,
            nnz: 4,
            incumbent: vec![0.0, 0.0],
            binary: vec![1, 1],
        }
    }

    #[test]
    fn s5_lmworker_run_attempt_reaches_feasible_and_offers() {
        let problem = make_problem_2x2_binary();
        let mut sink = S5Sink::new(&problem.integrality, 0);
        let exec = ExecutionContext::default();
        let trace = WorkerTrace::default();
        {
            let mut worker =
                LocalMipWorker::new(&problem, exec, 100_000, 100_000, 42, None, trace, None);
            worker.run_attempt(100_000, &mut sink);
            assert!(worker.best_feasible_);
        }
        assert!(sink.accepted() >= 1);
    }

    #[test]
    fn s5_lmworker_perturb_flips_about_20_percent_of_binaries() {
        let n = 1000;
        let mut solution = vec![0.0; n];
        let binary = vec![1u8; n];
        let integrality = vec![HighsVarType::kInteger; n];
        let col_lb = vec![0.0; n];
        let col_ub = vec![1.0; n];
        let mut rng = S5Rng::new(42);
        perturb_solution(
            &mut solution,
            &binary,
            &integrality,
            &col_lb,
            &col_ub,
            n as HighsInt,
            &mut rng,
        );
        let count = solution.iter().filter(|&&v| v == 1.0).count();
        assert!(count >= 150 && count <= 250, "count = {count}");
    }

    #[test]
    fn s5_lmworker_randomize_weights_draws_one_or_two() {
        let problem = make_problem_2x2_binary();
        let exec = ExecutionContext::default();
        let trace = WorkerTrace::default();

        let mut worker = LocalMipWorker::new(&problem, exec, 1000, 1000, 42, None, trace, None);
        let mut rng = S5Rng::new(7);
        worker.randomize_weights(&mut rng);
        for &w in worker.weights_for_test() {
            assert!(w == 1 || w == 2);
        }
        assert!(worker.obj_weight_for_test() == 1 || worker.obj_weight_for_test() == 2);

        let mut worker2 = LocalMipWorker::new(&problem, exec, 1000, 1000, 42, None, trace, None);
        let mut rng2 = S5Rng::new(7);
        worker2.randomize_weights(&mut rng2);
        assert_eq!(worker.weights_for_test(), worker2.weights_for_test());
        assert_eq!(worker.obj_weight_for_test(), worker2.obj_weight_for_test());
    }

    #[test]
    fn s5_lmworker_retain_keeps_least_violated_point() {
        let problem = make_problem_2x2_binary();
        let mut sink = S5Sink::new(&problem.integrality, 0);
        let exec = ExecutionContext::default();
        let trace = WorkerTrace::default();
        let mut holder = LeastViolatedPoint::default();
        {
            let mut worker = LocalMipWorker::new(
                &problem,
                exec,
                100000,
                100000,
                42,
                None,
                trace,
                Some(&mut holder),
            );
            worker.run_attempt(1000, &mut sink);
        }
        assert!(holder.best().is_some());
    }
}
