//! HiGHS source: pdlp/hipdlp/restart.hpp, pdlp/hipdlp/restart.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The adaptive restart scheme of HiPDLP. Back-pointers are dropped:
//! `initialize` and `update_beta` take the
//! `PrimalDualParams` they read or write; `passParams`, `passLogOptions`,
//! `passDebugPdlpLogFile`, `passDebugPdlpData` are not ported, and
//! `debugPdlpRestartLog` (line 71) is a no-op under `PDLP_DEBUG_LOG == 0`.

use super::defs::*;
use super::solver_results::SolverResults;
use crate::prelude::*;

/// HiGHS: restart.hpp `struct RestartInfo`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RestartInfo {
    pub should_restart: bool,
    /// If true, restart to average; otherwise, to current
    pub restart_to_average: bool,
}

impl Default for RestartInfo {
    fn default() -> Self {
        Self {
            should_restart: false,
            restart_to_average: true,
        }
    }
}

impl RestartInfo {
    /// HiGHS: restart.hpp `RestartInfo(bool, bool)`
    pub fn new(should_restart_: bool, restart_to_average_: bool) -> Self {
        Self {
            should_restart: should_restart_,
            restart_to_average: restart_to_average_,
        }
    }
}

/// HiGHS: restart.cc::compute_score (lines 21-25, file static)
pub fn compute_score(beta: f64, primal_feas: f64, dual_feas: f64, duality_gap: f64) -> f64 {
    (beta * primal_feas * primal_feas + dual_feas * dual_feas / beta + duality_gap * duality_gap)
        .sqrt()
}

/// HiGHS: restart.hpp `class RestartScheme`
#[derive(Clone, Debug)]
pub struct RestartScheme {
    pub last_restart_iter_: HighsInt,
    // State for adaptive restart
    pub primal_feas_last_restart_: f64,
    pub dual_feas_last_restart_: f64,
    pub duality_gap_last_restart_: f64,

    pub(crate) strategy_: RestartStrategy,
    pub(crate) fixed_restart_interval_: HighsInt,
    pub(crate) beta_: f64,
    pub(crate) primal_feas_last_candidate_: f64,
    pub(crate) dual_feas_last_candidate_: f64,
    pub(crate) duality_gap_last_candidate_: f64,
    pub(crate) sufficient_decay_factor_: f64,
    pub(crate) necessary_decay_factor_: f64,
    pub(crate) artificial_restart_threshold_: f64,
}

impl Default for RestartScheme {
    /// `beta_` has no initialiser in C++; `initialize` sets it.
    fn default() -> Self {
        Self {
            last_restart_iter_: 0,
            primal_feas_last_restart_: 0.0,
            dual_feas_last_restart_: 0.0,
            duality_gap_last_restart_: 0.0,
            strategy_: RestartStrategy::NO_RESTART,
            fixed_restart_interval_: 100,
            beta_: 0.0,
            primal_feas_last_candidate_: 0.0,
            dual_feas_last_candidate_: 0.0,
            duality_gap_last_candidate_: 0.0,
            sufficient_decay_factor_: 0.2,
            necessary_decay_factor_: 0.8,
            artificial_restart_threshold_: 0.36,
        }
    }
}

impl RestartScheme {
    /// HiGHS: restart.cc::RestartScheme::initialize (lines 28-38)
    ///
    /// `params` is `*params_`. `results` is unused upstream (kept).
    pub fn initialize(&mut self, results: &SolverResults, params: &PrimalDualParams) {
        let _ = results;
        self.strategy_ = params.restart_strategy;
        self.fixed_restart_interval_ = params.fixed_restart_interval;
        self.last_restart_iter_ = 0;

        // Stepsize ratio,
        // β = dBeta = dDualStep / dPrimalStep,
        // in the paper, primal weight is the ω:
        // ω = √β
        self.beta_ = params.omega * params.omega;
    }

    /// HiGHS: restart.cc::RestartScheme::check (lines 40-132)
    pub fn check(
        &mut self,
        current_iter: HighsInt,
        current_results: &SolverResults,
        average_results: &SolverResults,
    ) -> RestartInfo {
        if self.strategy_ != RestartStrategy::ADAPTIVE_RESTART {
            if self.strategy_ == RestartStrategy::FIXED_RESTART
                && current_iter - self.last_restart_iter_ >= self.fixed_restart_interval_
            {
                return RestartInfo::new(true, true);
            }
            return RestartInfo::new(false, false);
        }

        // Base Case: First check after a restart. Initialize state and exit.
        if current_iter == self.last_restart_iter_ {
            self.primal_feas_last_restart_ = current_results.primal_feasibility;
            self.dual_feas_last_restart_ = current_results.dual_feasibility;
            self.duality_gap_last_restart_ = current_results.duality_gap;

            self.primal_feas_last_candidate_ = current_results.primal_feasibility;
            self.dual_feas_last_candidate_ = current_results.dual_feasibility;
            self.duality_gap_last_candidate_ = current_results.duality_gap;

            return RestartInfo::new(false, false);
        }

        // 1. Calculate scores and determine the best candidate
        let mu_current = compute_score(
            self.beta_,
            current_results.primal_feasibility,
            current_results.dual_feasibility,
            current_results.duality_gap,
        );
        let mu_average = compute_score(
            self.beta_,
            average_results.primal_feasibility,
            average_results.dual_feasibility,
            average_results.duality_gap,
        );
        // debugPdlpRestartLog (line 71) is a no-op under PDLP_DEBUG_LOG == 0.
        // C++ `bool restart_to_average = false;` is assigned on both branches below.
        let (candidate_score, restart_to_average) = if mu_current < mu_average {
            (mu_current, false)
        } else {
            (mu_average, true)
        };

        // 2. Check the three restart conditions in order
        let mut should_restart = false;
        // PDHG_restart_choice restart_choice_for_logic = PDHG_NO_RESTART; // unused

        if ((current_iter - self.last_restart_iter_) as f64)
            >= self.artificial_restart_threshold_ * (current_iter as f64)
        {
            // Condition 1: Artificial Restart
            should_restart = true;
        } else {
            let mu_last_restart = compute_score(
                self.beta_,
                self.primal_feas_last_restart_,
                self.dual_feas_last_restart_,
                self.duality_gap_last_restart_,
            );

            if candidate_score < self.sufficient_decay_factor_ * mu_last_restart {
                // Condition 2: Sufficient Decay
                should_restart = true;
            } else {
                let mu_last_candidate = compute_score(
                    self.beta_,
                    self.primal_feas_last_candidate_,
                    self.dual_feas_last_candidate_,
                    self.duality_gap_last_candidate_,
                );

                if candidate_score < self.necessary_decay_factor_ * mu_last_restart
                    && candidate_score > mu_last_candidate
                {
                    // Condition 3: Necessary Decay
                    should_restart = true;
                }
            }
        }

        // 3. ALWAYS update the "last candidate" metrics for the next check
        if restart_to_average {
            self.primal_feas_last_candidate_ = average_results.primal_feasibility;
            self.dual_feas_last_candidate_ = average_results.dual_feasibility;
            self.duality_gap_last_candidate_ = average_results.duality_gap;
        } else {
            self.primal_feas_last_candidate_ = current_results.primal_feasibility;
            self.dual_feas_last_candidate_ = current_results.dual_feasibility;
            self.duality_gap_last_candidate_ = current_results.duality_gap;
        }
        // Lines 122-130: commented-out debug print block, dropped.
        RestartInfo::new(should_restart, restart_to_average)
    }

    /// HiGHS: restart.hpp::RestartScheme::getLastRestartIter
    pub fn getLastRestartIter(&self) -> HighsInt {
        self.last_restart_iter_
    }

    /// HiGHS: restart.hpp::RestartScheme::getBeta
    pub fn getBeta(&self) -> f64 {
        self.beta_
    }

    /// HiGHS: restart.hpp::RestartScheme::updateBeta
    ///
    /// `params` is `*params_`: `params_->omega = std::sqrt(beta)`.
    pub fn update_beta(&mut self, beta: f64, params: &mut PrimalDualParams) {
        self.beta_ = beta;
        params.omega = beta.sqrt();
    }

    /// HiGHS: restart.hpp::RestartScheme::setLastRestartIter (takes `int`)
    pub fn setLastRestartIter(&mut self, iter: i32) {
        self.last_restart_iter_ = iter;
    }

    /// HiGHS: restart.hpp::RestartScheme::getSufficientDecayFactor
    pub fn getSufficientDecayFactor(&self) -> f64 {
        self.sufficient_decay_factor_
    }

    /// HiGHS: restart.hpp::RestartScheme::getNecessaryDecayFactor
    pub fn getNecessaryDecayFactor(&self) -> f64 {
        self.necessary_decay_factor_
    }

    /// HiGHS: restart.hpp::RestartScheme::getArtificialRestartThreshold
    pub fn getArtificialRestartThreshold(&self) -> f64 {
        self.artificial_restart_threshold_
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn score_and_beta() {
        // sqrt(4*1 + 4/4 + 0) = sqrt(5)
        assert_eq!(compute_score(4.0, 1.0, 2.0, 0.0), 5.0f64.sqrt());
        let mut r = RestartScheme::default();
        let mut p = PrimalDualParams {
            omega: 3.0,
            restart_strategy: RestartStrategy::ADAPTIVE_RESTART,
            ..Default::default()
        };
        r.initialize(&SolverResults::default(), &p);
        assert_eq!(r.getBeta(), 9.0);
        r.update_beta(4.0, &mut p);
        assert_eq!(p.omega, 2.0);
    }

    #[test]
    fn first_check_records_state_without_restart() {
        let mut r = RestartScheme::default();
        let p = PrimalDualParams {
            omega: 1.0,
            restart_strategy: RestartStrategy::ADAPTIVE_RESTART,
            ..Default::default()
        };
        r.initialize(&SolverResults::default(), &p);
        let cur = SolverResults {
            primal_feasibility: 1.5,
            duality_gap: 2.5,
            ..Default::default()
        };
        let info = r.check(0, &cur, &SolverResults::default());
        assert_eq!(info, RestartInfo::new(false, false));
        assert_eq!(r.primal_feas_last_restart_, 1.5);
        assert_eq!(r.duality_gap_last_candidate_, 2.5);
    }
}
