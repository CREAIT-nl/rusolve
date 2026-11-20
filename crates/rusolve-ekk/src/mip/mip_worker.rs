//! HiGHS source: mip/HighsMipWorker.h, mip/HighsMipWorker.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsMipWorker.h::HighsMipWorker (lines 23-189)
///
/// Removed C++ members:
/// * C++: `const HighsMipSolver & mipsolver_`: removed; methods that use it take `ms: &mut HighsMipSolver`.
/// * C++: `const HighsMipSolverData & mipdata_`: removed; methods that use it take `ms: &mut HighsMipSolver`.
#[derive(Clone, Default, Debug)]
pub struct HighsMipWorker {
    /// C++ `HighsLpRelaxation* lp_`: index into `ms.mipdata_.lps`.
    pub lp_: usize,
    /// C++ `HighsDomain* globaldom_`: index into `ms.mipdata_.domains`.
    pub globaldom_: usize,
    /// C++ `HighsCutPool* cutpool_`: index into `ms.mipdata_.cutpools`.
    pub cutpool_: usize,
    /// C++ `HighsConflictPool* conflictpool_`: index into `ms.mipdata_.conflictpools`.
    pub conflictpool_: usize,
    /// C++ `HighsPseudocost* pseudocost_`: index into `ms.mipdata_.pseudocosts`.
    pub pseudocost_: usize,
    /// C++: `bool heuristics_allowed`
    pub heuristics_allowed: bool,
    /// C++: `HeurStatistics heur_stats`
    pub heur_stats: HighsMipWorker_HeurStatistics,
    /// C++: `SepaStatistics sepa_stats`
    pub sepa_stats: HighsMipWorker_SepaStatistics,
    /// C++: `std::unique_ptr<HighsSearch> search_ptr_`
    pub search_ptr_: Option<Box<HighsSearch>>,
    /// C++: `std::unique_ptr<HighsSeparation> sepa_ptr_`
    pub sepa_ptr_: Option<Box<HighsSeparation>>,
    /// C++: `HighsNodeQueue nodequeue`
    pub nodequeue: HighsNodeQueue,
    /// C++: `double upper_bound`
    pub upper_bound: f64,
    /// C++: `double upper_limit`
    pub upper_limit: f64,
    /// C++: `double optimality_limit`
    pub optimality_limit: f64,
    /// C++: `std::vector<std::tuple<std::vector<double>, double, int>> solutions_`
    pub solutions_: Vec<(Vec<f64>, f64, i32)>,
    /// C++: `HighsRandom randgen`
    pub randgen: HighsRandom,
    /// HiGHS 1.16 (#3192): HighsMipWorker.h::processedNodes (line 73): the nodes the
    /// worker's search stashed in this batch, with the search's `countTreeWeight` at the
    /// stash. Used only under `h16_mip_workers`; empty otherwise.
    pub processed_nodes: Vec<(HighsNodeQueue_OpenNode, bool)>,
    /// HiGHS 1.16 (#3192): HighsMipWorker.h::preparedNodes (line 74): the nodes assigned
    /// to the worker for this batch. Used only under `h16_mip_workers`; empty otherwise.
    pub prepared_nodes: Vec<HighsNodeQueue_OpenNode>,
    /// HiGHS 1.16 (#3192): HighsMipWorker.h::prepNodeIdx (line 75). C++ leaves it
    /// uninitialised until `processWorkerNodes` sets it to 0; here it starts at 0.
    pub prep_node_idx: usize,
    /// HiGHS 1.16 (#3124): HighsMipWorker.h::early_termination (line 80). Only written
    /// under `h16_mip_workers`.
    pub early_termination: bool,
}

/// HiGHS: HighsMipWorker.h::HighsMipWorker::SepaStatistics (lines 25-30)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsMipWorker_SepaStatistics {
    /// C++: `int64_t numNeighbourhoodQueries`
    pub num_neighbourhood_queries: i64,
    /// C++: `int64_t sepa_lp_iterations`
    pub sepa_lp_iterations: i64,
}

/// HiGHS: HighsMipWorker.h::HighsMipWorker::HeurStatistics (lines 31-56)
#[derive(Clone, Default, Debug)]
pub struct HighsMipWorker_HeurStatistics {
    /// C++: `int64_t total_repair_lp`
    pub total_repair_lp: i64,
    /// C++: `int64_t total_repair_lp_feasible`
    pub total_repair_lp_feasible: i64,
    /// C++: `int64_t total_repair_lp_iterations`
    pub total_repair_lp_iterations: i64,
    /// C++: `int64_t lp_iterations`
    pub lp_iterations: i64,
    /// C++: `double successObservations`
    pub success_observations: f64,
    /// C++: `HighsInt numSuccessObservations`
    pub num_success_observations: HighsInt,
    /// C++: `double infeasObservations`
    pub infeas_observations: f64,
    /// C++: `HighsInt numInfeasObservations`
    pub num_infeas_observations: HighsInt,
    /// C++: `HighsInt max_submip_level`
    pub max_submip_level: HighsInt,
    /// C++: `HighsModelStatus termination_status_`
    pub termination_status_: HighsModelStatus,
}

impl HighsMipWorker {
    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::HighsMipWorker (lines 13-35)
    /// C++: `HighsMipWorker(const HighsMipSolver & mipsolver, HighsLpRelaxation * lp, HighsDomain * domain, HighsCutPool * cutpool, HighsConflictPool * conflictpool, HighsPseudocost * pseudocost)`
    // PORT NOTE: the pointers are indices into `ms.mipdata_.{lps, domains, cutpools,
    // conflictpools, pseudocosts}`. `ms` is `&mut` because `HighsSearch::new` copies the
    // global domain (registers with the pools). The search's `setLpRelaxation(lp_)` /
    // the separation's `setLpRelaxation(lp_)` map to the worker's `lp_` (search: no-op,
    // the LP is `SearchCtx::lp`; separation: `Some(lp_)`).
    pub fn new(
        mipsolver: &mut HighsMipSolver,
        lp: usize,
        domain: usize,
        cutpool: usize,
        conflictpool: usize,
        pseudocost: usize,
    ) -> Self {
        let mut this = HighsMipWorker {
            lp_: lp,
            globaldom_: domain,
            cutpool_: cutpool,
            conflictpool_: conflictpool,
            pseudocost_: pseudocost,
            heuristics_allowed: true,
            heur_stats: HighsMipWorker_HeurStatistics::new(),
            sepa_stats: HighsMipWorker_SepaStatistics::new(),
            search_ptr_: None,
            sepa_ptr_: None,
            nodequeue: HighsNodeQueue::default(),
            upper_bound: mipsolver.mipdata_.upper_bound,
            upper_limit: mipsolver.mipdata_.upper_limit,
            optimality_limit: mipsolver.mipdata_.optimality_limit,
            solutions_: Vec::new(),
            randgen: HighsRandom::new(mipsolver.options_mip_.random_seed as HighsUInt),
            processed_nodes: Vec::new(),
            prepared_nodes: Vec::new(),
            prep_node_idx: 0,
            // HiGHS 1.16 (#3124): HighsMipWorker.cpp line 30
            early_termination: false,
        };
        this.reset_search(mipsolver);
        this.reset_sepa(mipsolver);
        this
    }

    /// HiGHS 1.16 (#3192): HighsMipWorker.cpp::HighsMipWorker::getOptimalityLimit (lines 152-158)
    /// C++: `double getOptimalityLimit() const`
    ///
    /// Moved here from `HighsSearch::getOptimalityLimit` (same body); called only on the
    /// `h16_mip_workers` path.
    pub fn get_optimality_limit(&self, ms: &HighsMipSolver) -> f64 {
        if !ms.parallel_lock_active() {
            ms.mipdata_.optimality_limit
        } else {
            self.optimality_limit
        }
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::getMipSolver (lines 37-39)
    /// C++: `const HighsMipSolver & getMipSolver() const`
    // PORT NOTE: the solver is the context `ms` (M1); identity function kept for 1:1 call sites.
    pub fn get_mip_solver(ms: &HighsMipSolver) -> &HighsMipSolver {
        ms
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::resetSearch (lines 41-46)
    /// C++: `void resetSearch()`
    // PORT NOTE: `self` may be detached, or not yet in `ms.mipdata_.workers` (constructor):
    // `HighsSearch::new` reads only the worker's `globaldom_`. `search_ptr_.reset()`
    // destroys the old search, so its local domain is released first.
    pub fn reset_search(&mut self, ms: &mut HighsMipSolver) {
        if let Some(mut old) = self.search_ptr_.take() {
            old.release(ms);
        }
        let search = HighsSearch::new(ms, self.globaldom_);
        self.search_ptr_ = Some(Box::new(search));
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::resetSepa (lines 48-52)
    /// C++: `void resetSepa()`
    pub fn reset_sepa(&mut self, ms: &mut HighsMipSolver) {
        let mut sepa = HighsSeparation::new(ms);
        sepa.set_lp_relaxation(Some(self.lp_));
        self.sepa_ptr_ = Some(Box::new(sepa));
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::addIncumbent (lines 54-75)
    /// C++: `bool addIncumbent(const std::vector<double> & sol, double solobj, int solution_source)`
    pub fn add_incumbent(
        &mut self,
        ms: &mut HighsMipSolver,
        sol: &[f64],
        solobj: f64,
        solution_source: i32,
    ) -> bool {
        if solobj < self.upper_bound {
            // Get the transformed objective and solution if required
            let transformed_solobj = self.transform_new_integer_feasible_solution(ms, sol);
            if transformed_solobj.0 && transformed_solobj.1 < self.upper_bound {
                self.upper_bound = transformed_solobj.1;
                let new_upper_limit = ms.compute_new_upper_limit(self.upper_bound, 0.0, 0.0);
                if new_upper_limit < self.upper_limit {
                    self.upper_limit = new_upper_limit;
                    self.optimality_limit = ms.compute_new_upper_limit(
                        self.upper_bound,
                        ms.options_mip_.mip_abs_gap,
                        ms.options_mip_.mip_rel_gap,
                    );
                }
            }
            // Can't repair solutions locally, so also buffer infeasible ones
            self.solutions_
                .push((sol.to_vec(), solobj, solution_source));
        }
        true
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::transformNewIntegerFeasibleSolution (lines 77-111)
    /// C++: `std::pair<bool, double> transformNewIntegerFeasibleSolution(const std::vector<double> & sol) const`
    pub fn transform_new_integer_feasible_solution(
        &self,
        ms: &mut HighsMipSolver,
        sol: &[f64],
    ) -> (bool, f64) {
        let mut solution = HighsSolution {
            col_value: sol.to_vec(),
            value_valid: true,
            ..Default::default()
        };

        // Perform primal postsolve to get the original column values
        ms.mipdata_
            .post_solve_stack
            .undo_primal(&ms.options_mip_, &mut solution, -1, true);

        // Determine the row values, as they aren't computed in primal postsolve
        let return_status = crate::lp_data::presolve_support::calculate_row_values_quad(
            &ms.orig_model_,
            &mut solution,
            -1,
        );
        if kAllowDeveloperAssert {
            debug_assert!(return_status == HighsStatus::kOk);
        }

        // compute the objective value in the original space
        let mut bound_violation: f64 = 0.0;
        let mut row_violation: f64 = 0.0;
        let mut integrality_violation: f64 = 0.0;

        let mut mipsolver_quad_objective_value = HighsCDouble::from(0.0);

        let feasible = ms.solution_feasible(
            Some(&ms.orig_model_),
            &solution.col_value,
            Some(&solution.row_value),
            &mut bound_violation,
            &mut row_violation,
            &mut integrality_violation,
            &mut mipsolver_quad_objective_value,
        );

        let transformed_solobj = f64::from(
            ((ms.orig_model_.sense_ as HighsInt) as f64) * mipsolver_quad_objective_value
                - mip_model!(ms).offset_,
        );

        (feasible, transformed_solobj)
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::trySolution (lines 113-146)
    /// C++: `bool trySolution(const std::vector<double> & solution, int solution_source)`
    pub fn try_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        solution: &[f64],
        solution_source: i32,
    ) -> bool {
        if solution.len() as HighsInt != ms.num_col() {
            return false;
        }

        let mut obj = HighsCDouble::from(0.0);
        // HiGHS 1.16 #3149 (h16_cdouble): MipWorker.cpp::trySolution lines 133 (1.16.0)
        let h16_cdouble_obj = ms.options_mip_.h16_cdouble;

        for i in 0..ms.num_col() {
            let iu = i as usize;
            if solution[iu] < mip_model!(ms).col_lower_[iu] - ms.mipdata_.feastol {
                return false;
            }
            if solution[iu] > mip_model!(ms).col_upper_[iu] + ms.mipdata_.feastol {
                return false;
            }
            if ms.variable_type_2(i) == HighsVarType::kInteger
                && fractionality(solution[iu], None) > ms.mipdata_.feastol
            {
                return false;
            }

            if h16_cdouble_obj {
                obj += HighsCDouble::from(ms.col_cost_2(i)) * solution[iu];
            } else {
                obj += ms.col_cost_2(i) * solution[iu];
            }
        }

        for i in 0..ms.num_row() {
            let iu = i as usize;
            let mut rowactivity: f64 = 0.0;

            let start = ms.mipdata_.a_rstart_[iu];
            let end = ms.mipdata_.a_rstart_[iu + 1];

            for j in start..end {
                let ju = j as usize;
                rowactivity +=
                    solution[ms.mipdata_.a_rindex_[ju] as usize] * ms.mipdata_.a_rvalue_[ju];
            }

            if rowactivity > ms.row_upper_2(i) + ms.mipdata_.feastol {
                return false;
            }
            if rowactivity < ms.row_lower_2(i) - ms.mipdata_.feastol {
                return false;
            }
        }

        self.add_incumbent(ms, solution, f64::from(obj), solution_source)
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::resetSepaStats (lines 148-151)
    /// C++: `void resetSepaStats()`
    pub fn reset_sepa_stats(&mut self) {
        self.sepa_stats.num_neighbourhood_queries = 0;
        self.sepa_stats.sepa_lp_iterations = 0;
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::updateHeurStatsLpIters (lines 153-161)
    /// C++: `void updateHeurStatsLpIters(int64_t lp_iters, int64_t total_repair_lp, int64_t total_repair_lp_feasible, int64_t total_repair_lp_iters)`
    pub fn update_heur_stats_lp_iters(
        &mut self,
        lp_iters: i64,
        total_repair_lp: i64,
        total_repair_lp_feasible: i64,
        total_repair_lp_iters: i64,
    ) {
        self.heur_stats.lp_iterations += lp_iters;
        self.heur_stats.total_repair_lp += total_repair_lp;
        self.heur_stats.total_repair_lp_feasible += total_repair_lp_feasible;
        self.heur_stats.total_repair_lp_iterations += total_repair_lp_iters;
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::updateHeurStatsInfeasObservations (lines 163-166)
    /// C++: `void updateHeurStatsInfeasObservations(double fixingRate)`
    pub fn update_heur_stats_infeas_observations(&mut self, fixing_rate: f64) {
        self.heur_stats.infeas_observations += fixing_rate;
        self.heur_stats.num_infeas_observations += 1;
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::updateHeurStatsSuccessObservations (lines 168-171)
    /// C++: `void updateHeurStatsSuccessObservations(double fixingRate)`
    pub fn update_heur_stats_success_observations(&mut self, fixing_rate: f64) {
        self.heur_stats.success_observations += fixing_rate;
        self.heur_stats.num_success_observations += 1;
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::getHeurStatsValues (lines 173-189)
    /// C++: `void getHeurStatsValues(int64_t & total_repair_lp, int64_t & total_repair_lp_feasible, int64_t & total_repair_lp_iterations, int64_t & lp_iterations, double & successObservations, HighsInt & numSuccessObservations, double & infeasObservations, HighsInt & numInfeasObservations, HighsInt & max_submip_level, HighsModelStatus & termination_status) const`
    pub fn get_heur_stats_values(
        &self,
        total_repair_lp: &mut i64,
        total_repair_lp_feasible: &mut i64,
        total_repair_lp_iterations: &mut i64,
        lp_iterations: &mut i64,
        success_observations: &mut f64,
        num_success_observations: &mut HighsInt,
        infeas_observations: &mut f64,
        num_infeas_observations: &mut HighsInt,
        max_submip_level: &mut HighsInt,
        termination_status: &mut HighsModelStatus,
    ) {
        *total_repair_lp = self.heur_stats.total_repair_lp;
        *total_repair_lp_feasible = self.heur_stats.total_repair_lp_feasible;
        *total_repair_lp_iterations = self.heur_stats.total_repair_lp_iterations;
        *lp_iterations = self.heur_stats.lp_iterations;
        *success_observations = self.heur_stats.success_observations;
        *num_success_observations = self.heur_stats.num_success_observations;
        *infeas_observations = self.heur_stats.infeas_observations;
        *num_infeas_observations = self.heur_stats.num_infeas_observations;
        *max_submip_level = self.heur_stats.max_submip_level;
        *termination_status = self.heur_stats.termination_status_;
    }

    /// HiGHS: HighsMipWorker.cpp::HighsMipWorker::resetHeurStats (lines 191-202)
    /// C++: `void resetHeurStats()`
    pub fn reset_heur_stats(&mut self) {
        self.heur_stats.total_repair_lp = 0;
        self.heur_stats.total_repair_lp_feasible = 0;
        self.heur_stats.total_repair_lp_iterations = 0;
        self.heur_stats.lp_iterations = 0;
        self.heur_stats.max_submip_level = 0;
        self.heur_stats.termination_status_ = HighsModelStatus::kNotset;
        self.heur_stats.success_observations = 0.0;
        self.heur_stats.num_success_observations = 0;
        self.heur_stats.infeas_observations = 0.0;
        self.heur_stats.num_infeas_observations = 0;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::~HighsMipWorker (lines 89-92)
    /// C++: `~HighsMipWorker()`
    // PORT NOTE: takes `ms` to release the search's local domain (`search_ptr_.reset()`).
    pub fn release(&mut self, ms: &mut HighsMipSolver) {
        if let Some(mut search) = self.search_ptr_.take() {
            search.release(ms);
        }
        self.sepa_ptr_ = None;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getGlobalDomain (lines 98-98)
    /// C++: `HighsDomain & getGlobalDomain() const`
    // PORT NOTE: index into `ms.mipdata_` (the C++ returns `*globaldom_`).
    pub fn get_global_domain(&self) -> usize {
        self.globaldom_
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::setGlobalDomain (lines 100-100)
    /// C++: `void setGlobalDomain(HighsDomain * globaldom)`
    // PORT NOTE: index into `ms.mipdata_` (the C++ pointer).
    pub fn set_global_domain(&mut self, globaldom: usize) {
        self.globaldom_ = globaldom;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getPseudocost (lines 102-102)
    /// C++: `HighsPseudocost & getPseudocost() const`
    // PORT NOTE: index into `ms.mipdata_` (the C++ returns `*pseudocost_`).
    pub fn get_pseudocost(&self) -> usize {
        self.pseudocost_
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::setPseudocost (lines 104-104)
    /// C++: `void setPseudocost(HighsPseudocost * pseudocost)`
    // PORT NOTE: index into `ms.mipdata_` (the C++ pointer).
    pub fn set_pseudocost(&mut self, pseudocost: usize) {
        self.pseudocost_ = pseudocost;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getConflictPool (lines 106-106)
    /// C++: `HighsConflictPool & getConflictPool() const`
    // PORT NOTE: index into `ms.mipdata_` (the C++ returns `*conflictpool_`).
    pub fn get_conflict_pool(&self) -> usize {
        self.conflictpool_
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::setConflictPool (lines 108-110)
    /// C++: `void setConflictPool(HighsConflictPool * conflictpool)`
    // PORT NOTE: index into `ms.mipdata_` (the C++ pointer).
    pub fn set_conflict_pool(&mut self, conflictpool: usize) {
        self.conflictpool_ = conflictpool;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getCutPool (lines 112-112)
    /// C++: `HighsCutPool & getCutPool() const`
    // PORT NOTE: index into `ms.mipdata_` (the C++ returns `*cutpool_`).
    pub fn get_cut_pool(&self) -> usize {
        self.cutpool_
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::setCutPool (lines 114-114)
    /// C++: `void setCutPool(HighsCutPool * cutpool)`
    // PORT NOTE: index into `ms.mipdata_` (the C++ pointer).
    pub fn set_cut_pool(&mut self, cutpool: usize) {
        self.cutpool_ = cutpool;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getLpRelaxation (lines 116-116)
    /// C++: `HighsLpRelaxation & getLpRelaxation() const`
    // PORT NOTE: index into `ms.mipdata_` (the C++ returns `*lp_`).
    pub fn get_lp_relaxation(&self) -> usize {
        self.lp_
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::setLpRelaxation (lines 118-118)
    /// C++: `void setLpRelaxation(HighsLpRelaxation * lp)`
    // PORT NOTE: index into `ms.mipdata_` (the C++ pointer).
    pub fn set_lp_relaxation(&mut self, lp: usize) {
        self.lp_ = lp;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::setAllowHeuristics (lines 128-128)
    /// C++: `void setAllowHeuristics(const bool allowed)`
    pub fn set_allow_heuristics(&mut self, allowed: bool) {
        self.heuristics_allowed = allowed;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getAllowHeuristics (lines 130-130)
    /// C++: `bool getAllowHeuristics() const`
    pub fn get_allow_heuristics(&self) -> bool {
        self.heuristics_allowed
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getNumNeighbourhoodQueries (lines 132-134)
    /// C++: `int64_t & getNumNeighbourhoodQueries()`
    pub fn get_num_neighbourhood_queries(&mut self) -> &mut i64 {
        &mut self.sepa_stats.num_neighbourhood_queries
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getSepaLpIterations (lines 136-136)
    /// C++: `int64_t & getSepaLpIterations()`
    pub fn get_sepa_lp_iterations(&mut self) -> &mut i64 {
        &mut self.sepa_stats.sepa_lp_iterations
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getHeurLpIterations (lines 148-148)
    /// C++: `int64_t & getHeurLpIterations()`
    pub fn get_heur_lp_iterations(&mut self) -> &mut i64 {
        &mut self.heur_stats.lp_iterations
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getHeurNumSuccessObservations (lines 150-152)
    /// C++: `HighsInt & getHeurNumSuccessObservations()`
    pub fn get_heur_num_success_observations(&mut self) -> &mut HighsInt {
        &mut self.heur_stats.num_success_observations
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getHeurNumInfeasObservations (lines 154-156)
    /// C++: `HighsInt & getHeurNumInfeasObservations()`
    pub fn get_heur_num_infeas_observations(&mut self) -> &mut HighsInt {
        &mut self.heur_stats.num_infeas_observations
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getHeurSuccessObservations (lines 158-160)
    /// C++: `double & getHeurSuccessObservations()`
    pub fn get_heur_success_observations(&mut self) -> &mut f64 {
        &mut self.heur_stats.success_observations
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getHeurInfeasObservations (lines 162-162)
    /// C++: `double & getHeurInfeasObservations()`
    pub fn get_heur_infeas_observations(&mut self) -> &mut f64 {
        &mut self.heur_stats.infeas_observations
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::setHeurTerminationStatus (lines 164-166)
    /// C++: `void setHeurTerminationStatus(HighsModelStatus status)`
    pub fn set_heur_termination_status(&mut self, status: HighsModelStatus) {
        self.heur_stats.termination_status_ = status;
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::getHeurTerminationStatus (lines 168-170)
    /// C++: `HighsModelStatus getHeurTerminationStatus() const`
    pub fn get_heur_termination_status(&self) -> HighsModelStatus {
        self.heur_stats.termination_status_
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::terminatorTerminated (lines 172-174)
    /// C++: `bool terminatorTerminated() const`
    pub fn terminator_terminated(&self) -> bool {
        self.heur_stats.termination_status_ != HighsModelStatus::kNotset
    }

    /// HiGHS: HighsMipWorker.h::HighsMipWorker::updateHeurStatsMaxSubMipLevel (lines 176-179)
    /// C++: `void updateHeurStatsMaxSubMipLevel(HighsInt new_level)`
    pub fn update_heur_stats_max_sub_mip_level(&mut self, new_level: HighsInt) {
        self.heur_stats.max_submip_level = cpp_max(new_level, self.heur_stats.max_submip_level);
    }
}

impl HighsMipWorker_SepaStatistics {
    /// HiGHS: HighsMipWorker.h::HighsMipWorker::SepaStatistics::SepaStatistics (lines 26-26)
    /// C++: `SepaStatistics()`
    pub fn new() -> Self {
        HighsMipWorker_SepaStatistics {
            num_neighbourhood_queries: 0,
            sepa_lp_iterations: 0,
        }
    }
}

impl HighsMipWorker_HeurStatistics {
    /// HiGHS: HighsMipWorker.h::HighsMipWorker::HeurStatistics::HeurStatistics (lines 32-43)
    /// C++: `HeurStatistics()`
    pub fn new() -> Self {
        HighsMipWorker_HeurStatistics {
            total_repair_lp: 0,
            total_repair_lp_feasible: 0,
            total_repair_lp_iterations: 0,
            lp_iterations: 0,
            success_observations: 0.0,
            num_success_observations: 0,
            infeas_observations: 0.0,
            num_infeas_observations: 0,
            max_submip_level: 0,
            termination_status_: HighsModelStatus::kNotset,
        }
    }
}
