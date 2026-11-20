//! HiGHS source: mip/HighsMipSolverData.h (class definitions and header-inline methods)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;
use std::sync::Arc;

mod part1;
pub use part1::*;
mod part2;
pub use part2::*;
mod part3;
pub use part3::*;

/// HiGHS: HighsMipSolverData.h::HighsPrimaDualIntegral (lines 34-41)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsPrimaDualIntegral {
    /// C++: `double value`
    pub value: f64,
    /// C++: `double prev_lb`
    pub prev_lb: f64,
    /// C++: `double prev_ub`
    pub prev_ub: f64,
    /// C++: `double prev_gap`
    pub prev_gap: f64,
    /// C++: `double prev_time`
    pub prev_time: f64,
}

/// HiGHS 1.16 (#3124): HighsMipSolverData.h::worker_lp_iterations_stop (line 81)
/// C++: `std::atomic<int64_t> worker_lp_iterations_stop`
///
/// The smallest local LP iteration count at which a worker asked for early termination
/// in the current batch (`i64::MAX`: none). Read and written only under
/// `h16_mip_workers`. Workers of one batch update it concurrently, so it is an atomic, as
/// in C++; `Clone` copies the current value and `Default` is `i64::MAX` (the constructor
/// value, HighsMipSolverData.cpp line 86).
#[derive(Debug)]
pub struct H16WorkerLpIterationsStop(pub std::sync::atomic::AtomicI64);

impl Default for H16WorkerLpIterationsStop {
    fn default() -> Self {
        Self(std::sync::atomic::AtomicI64::new(i64::MAX))
    }
}

impl Clone for H16WorkerLpIterationsStop {
    fn clone(&self) -> Self {
        Self(std::sync::atomic::AtomicI64::new(
            self.0.load(std::sync::atomic::Ordering::Relaxed),
        ))
    }
}

/// HiGHS: HighsMipSolverData.h::HighsMipSolverData (lines 70-277)
///
/// Removed C++ members:
/// * C++ `HighsMipSolver& mipsolver`: removed (M1: HighsMipSolverData methods are `impl HighsMipSolver`).
#[derive(Clone, Default, Debug)]
pub struct HighsMipSolverData {
    /// C++ `std::deque<HighsLpRelaxation> lps` (M2: each element is a detachable `Slot`).
    pub lps: Vec<Slot<HighsLpRelaxation>>,
    /// C++ `std::deque<HighsCutPool> cutpools` (M2: each element is a detachable `Slot`).
    pub cutpools: Vec<Slot<HighsCutPool>>,
    /// C++ `std::deque<HighsConflictPool> conflictpools` (M2: each element is a detachable `Slot`).
    pub conflictpools: Vec<Slot<HighsConflictPool>>,
    /// C++ `std::deque<HighsDomain> domains` (M2: each element is a detachable `Slot`).
    pub domains: Vec<Slot<HighsDomain>>,
    /// C++ `std::deque<HighsPseudocost> pseudocosts` (M2: each element is a detachable `Slot`).
    pub pseudocosts: Vec<Slot<HighsPseudocost>>,
    /// C++ `std::deque<HighsMipWorker> workers` (M2: each element is a detachable `Slot`).
    pub workers: Vec<Slot<HighsMipWorker>>,
    /// C++: `bool parallel_lock`
    pub parallel_lock: bool,
    /// C++ `HighsPrimalHeuristics heuristics` (M2: detachable `Slot`).
    pub heuristics: Slot<HighsPrimalHeuristics>,
    /// C++ `HighsCliqueTable cliquetable` (M2: detachable `Slot`).
    pub cliquetable: Slot<HighsCliqueTable>,
    /// C++ `HighsImplications implications` (M2: detachable `Slot`).
    pub implications: Slot<HighsImplications>,
    /// C++ `HighsRedcostFixing redcostfixing` (M2: detachable `Slot`).
    pub redcostfixing: Slot<HighsRedcostFixing>,
    /// C++ `HighsObjectiveFunction objectiveFunction` (M2: detachable `Slot`).
    pub objective_function: Slot<HighsObjectiveFunction>,
    /// C++: `presolve::HighsPostsolveStack postSolveStack`
    pub post_solve_stack: HighsPostsolveStack,
    /// C++: `HighsPresolveStatus presolve_status`
    pub presolve_status: HighsPresolveStatus,
    /// C++: `HighsLp presolvedModel`
    pub presolved_model: HighsLp,
    /// C++: `bool cliquesExtracted`
    pub cliques_extracted: bool,
    /// C++: `bool rowMatrixSet`
    pub row_matrix_set: bool,
    /// C++: `bool analyticCenterComputed`
    pub analytic_center_computed: bool,
    /// C++: `HighsModelStatus analyticCenterStatus`
    pub analytic_center_status: HighsModelStatus,
    /// C++: `bool detectSymmetries`
    pub detect_symmetries: bool,
    /// C++: `HighsInt numRestarts`
    pub num_restarts: HighsInt,
    /// C++: `HighsInt numRestartsRoot`
    pub num_restarts_root: HighsInt,
    /// rusolve `mip_cutoff_bigm` (Savelsbergh 1994 coefficient tightening): the cutoff at the last
    /// cutoff-driven restart and how many such restarts ran
    pub cutoff_bigm_last: Option<f64>,
    pub cutoff_bigm_restarts: HighsInt,
    /// C++: `HighsInt numCliqueEntriesAfterPresolve`
    pub num_clique_entries_after_presolve: HighsInt,
    /// C++: `HighsInt numCliqueEntriesAfterFirstPresolve`
    pub num_clique_entries_after_first_presolve: HighsInt,
    /// C++: `std::vector<HighsInt> ARstart_`
    pub a_rstart_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ARindex_`
    pub a_rindex_: Vec<HighsInt>,
    /// C++: `std::vector<double> ARvalue_`
    pub a_rvalue_: Vec<f64>,
    /// C++: `std::vector<double> maxAbsRowCoef`
    pub max_abs_row_coef: Vec<f64>,
    /// C++: `std::vector<uint8_t> rowintegral`
    pub rowintegral: Vec<u8>,
    /// C++: `std::vector<HighsInt> uplocks`
    pub uplocks: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> downlocks`
    pub downlocks: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> integer_cols`
    pub integer_cols: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> implint_cols`
    pub implint_cols: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> integral_cols`
    pub integral_cols: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> continuous_cols`
    pub continuous_cols: Vec<HighsInt>,
    /// C++: `HighsSymmetries symmetries`
    ///
    /// PORT NOTE (architecture): not a duplicate of
    /// [`HighsMipSolverData_SymmetryDetectionData::symmetries`]; the C++ has both
    /// (HighsMipSolverData.h:111 and :175). `finishSymmetryDetection` moves the
    /// detection result into this one.
    pub symmetries: HighsSymmetries,
    /// C++: `std::shared_ptr<const StabilizerOrbits> globalOrbits`
    pub global_orbits: Option<Arc<StabilizerOrbits>>,
    /// C++: `double feastol`
    pub feastol: f64,
    /// C++: `double epsilon`
    pub epsilon: f64,
    /// C++: `double heuristic_effort`
    pub heuristic_effort: f64,
    /// C++: `int64_t dispfreq`
    pub dispfreq: i64,
    /// C++: `std::vector<double> analyticCenter`
    pub analytic_center: Vec<f64>,
    /// C++: `std::vector<double> firstlpsol`
    pub firstlpsol: Vec<f64>,
    /// C++: `std::vector<double> rootlpsol`
    pub rootlpsol: Vec<f64>,
    /// C++: `double firstlpsolobj`
    pub firstlpsolobj: f64,
    /// C++: `HighsBasis firstrootbasis`
    pub firstrootbasis: HighsBasis,
    /// C++: `double rootlpsolobj`
    pub rootlpsolobj: f64,
    /// C++: `HighsInt numintegercols`
    pub numintegercols: HighsInt,
    /// C++: `HighsInt maxTreeSizeLog2`
    pub max_tree_size_log2: HighsInt,
    /// C++: `HighsCDouble pruned_treeweight`
    pub pruned_treeweight: HighsCDouble,
    /// C++: `double avgrootlpiters`
    pub avgrootlpiters: f64,
    /// C++: `double disptime`
    pub disptime: f64,
    /// C++: `double last_disptime`
    pub last_disptime: f64,
    /// C++: `int64_t firstrootlpiters`
    pub firstrootlpiters: i64,
    /// C++: `int64_t num_nodes`
    pub num_nodes: i64,
    /// C++: `int64_t num_leaves`
    pub num_leaves: i64,
    /// C++: `int64_t num_leaves_before_run`
    pub num_leaves_before_run: i64,
    /// C++: `int64_t num_nodes_before_run`
    pub num_nodes_before_run: i64,
    /// C++: `int64_t total_repair_lp`
    pub total_repair_lp: i64,
    /// C++: `int64_t total_repair_lp_feasible`
    pub total_repair_lp_feasible: i64,
    /// C++: `int64_t total_repair_lp_iterations`
    pub total_repair_lp_iterations: i64,
    /// C++: `int64_t total_lp_iterations`
    pub total_lp_iterations: i64,
    /// C++: `int64_t heuristic_lp_iterations`
    pub heuristic_lp_iterations: i64,
    /// C++: `int64_t sepa_lp_iterations`
    pub sepa_lp_iterations: i64,
    /// C++: `int64_t sb_lp_iterations`
    pub sb_lp_iterations: i64,
    /// C++: `int64_t total_lp_iterations_before_run`
    pub total_lp_iterations_before_run: i64,
    /// C++: `int64_t heuristic_lp_iterations_before_run`
    pub heuristic_lp_iterations_before_run: i64,
    /// C++: `int64_t sepa_lp_iterations_before_run`
    pub sepa_lp_iterations_before_run: i64,
    /// C++: `int64_t sb_lp_iterations_before_run`
    pub sb_lp_iterations_before_run: i64,
    /// C++: `int64_t num_disp_lines`
    pub num_disp_lines: i64,
    /// C++: `HighsInt numImprovingSols`
    pub num_improving_sols: HighsInt,
    /// C++: `double lower_bound`
    pub lower_bound: f64,
    /// C++: `double upper_bound`
    pub upper_bound: f64,
    /// C++: `double upper_limit`
    pub upper_limit: f64,
    /// C++: `double optimality_limit`
    pub optimality_limit: f64,
    /// C++: `std::vector<double> incumbent`
    pub incumbent: Vec<f64>,
    /// C++ `HighsNodeQueue nodequeue` (M2: detachable `Slot`).
    pub nodequeue: Slot<HighsNodeQueue>,
    /// C++: `HighsPrimaDualIntegral primal_dual_integral`
    pub primal_dual_integral: HighsPrimaDualIntegral,
    /// C++: `HighsDebugSol debugSolution`
    pub debug_solution: HighsDebugSol,
    /// HiGHS 1.16 (#3124): `std::atomic<int64_t> worker_lp_iterations_stop` (h16_mip_workers only).
    pub worker_lp_iterations_stop: H16WorkerLpIterationsStop,
    /// mip_conflict_resolution: SCIP conflict_resolution.c workspaces, one per conflict pool; inert on the off path.
    pub conflict_resolution: crate::mip::conflict_resolution::S2Workspaces,
    /// mip_symmetry_dynamic: dynamic symmetry data (SCIP symmetry_orbitopal.c / symmetry_lexred.c); None when the option is off.
    pub(crate) symmetry_dynamic: Option<Box<crate::mip::symmetry_dynamic::S3SymData>>,
    /// mip_heuristic_run_oneopt / mip_heuristic_run_crossover: SCIP heur_oneopt.c / heur_crossover.c state and primal.c solution pool; None when both options are off.
    pub(crate) s6: Option<Box<crate::mip::scip_heuristics::S6State>>,
    /// rusolve (mip_alns): bandit LNS state (Hendel 2022); None while the option is off.
    pub(crate) alns: Option<Box<crate::mip::scip_heuristics::alns::AlnsState>>,
    /// rusolve (mip_proximity): proximity search state (Fischetti, Monaci 2014); None while the option is off.
    pub(crate) proximity: Option<Box<crate::mip::scip_heuristics::proximity::ProximityState>>,
    /// rusolve (mip_local_branching): local branching state (Fischetti, Lodi 2003); None while the option is off.
    pub(crate) local_branching:
        Option<Box<crate::mip::scip_heuristics::local_branching::LocalBranchingState>>,
    /// mip_heuristic_run_fpr / mip_heuristic_run_local_mip: mip-heuristics patch field (apply_patch.cmake:654-664); written only on the path, never read by HiGHS code.
    pub heuristic_effort_used: usize,
    /// FPR / Local-MIP: accepted offers and per-heuristic charges of the presolve chain (test-observable); None when both options are off.
    pub(crate) s5_report: Option<Box<crate::mip::fpr_local_mip::S5Report>>,
    /// rusolve (mip_ls_int_lp_repair): LS end-point LP repair counters; None while the option
    /// is off.
    pub(crate) ls_int_lp_repair: Option<Box<crate::mip::ls_int_lp_repair::LsIntLpRepair>>,
}

/// HiGHS: HighsMipSolverData.h::HighsMipSolverData::SymmetryDetectionData (lines 173-177)
#[derive(Clone, Debug)]
pub struct HighsMipSolverData_SymmetryDetectionData {
    /// C++: `HighsSymmetryDetection symDetection`
    pub sym_detection: HighsSymmetryDetection,
    /// C++: `HighsSymmetries symmetries`
    pub symmetries: HighsSymmetries,
    /// C++: `double detectionTime = 0.0`
    pub detection_time: f64,
}

impl Default for HighsMipSolverData_SymmetryDetectionData {
    fn default() -> Self {
        Self {
            sym_detection: Default::default(),
            symmetries: Default::default(),
            detection_time: 0.0,
        }
    }
}

/// C++ `highs::parallel::TaskGroup` as `evaluateRootNode` declares it (HighsMipSolverData.cpp
/// 1995) and passes it to `start/finishAnalyticCenterComputation` and
/// `start/finishSymmetryDetection` (`const highs::parallel::TaskGroup&`).
///
/// The real scoped group ([`crate::parallel::highs_parallel::TaskGroup`], rule (g)):
/// `evaluateRootNode` runs its whole body inside `highs_parallel::task_group(|tg| ..)`, so
/// the group lives from line 1995 to the end of the function across every `goto restart`,
/// and every `return` runs `~TaskGroup` (`cancel`, `taskWait`). At one worker nothing
/// steals: `spawn` only pushes, `sync` pops the top (LIFO) and runs it, a cancelled task
/// never runs. At more workers a thief may run the analytic-centre or symmetry task before
/// the owner reaches `sync` / `cancel`.
pub type HighsMipTaskGroup<'scope, 'env> = crate::parallel::highs_parallel::TaskGroup<'scope, 'env>;

/// HiGHS: HighsMipSolverData.h::MipSolutionSource (plain C enum: module consts).
pub type MipSolutionSource = i32;
pub const kSolutionSourceNone: i32 = -1;
pub const kSolutionSourceMin: i32 = -1;
pub const kSolutionSourceBranching: i32 = 0;
pub const kSolutionSourceCentralRounding: i32 = 1;
pub const kSolutionSourceFeasibilityPump: i32 = 2;
pub const kSolutionSourceHeuristic: i32 = 3;
pub const kSolutionSourceShifting: i32 = 4;
pub const kSolutionSourceFeasibilityJump: i32 = 5;
pub const kSolutionSourceSubMip: i32 = 6;
pub const kSolutionSourceEmptyMip: i32 = 7;
pub const kSolutionSourceRandomizedRounding: i32 = 8;
pub const kSolutionSourceSolveLp: i32 = 9;
pub const kSolutionSourceEvaluateNode: i32 = 10;
pub const kSolutionSourceUnbounded: i32 = 11;
pub const kSolutionSourceUserSolution: i32 = 12;
pub const kSolutionSourceHighsSolution: i32 = 13;
pub const kSolutionSourceZiRound: i32 = 14;
pub const kSolutionSourceTrivialL: i32 = 15;
pub const kSolutionSourceTrivialP: i32 = 16;
pub const kSolutionSourceTrivialU: i32 = 17;
pub const kSolutionSourceTrivialZ: i32 = 18;
pub const kSolutionSourceCleanup: i32 = 19;
pub const kSolutionSourceCount: i32 = 20;
/// mip_heuristic_run_fpr / mip_heuristic_run_local_mip: mip-heuristics apply_patch.cmake lines 654-703 at 36f7d9d
pub const kSolutionSourceFPR: i32 = 20;
/// mip_heuristic_run_fpr / mip_heuristic_run_local_mip: mip-heuristics apply_patch.cmake lines 654-703 at 36f7d9d
pub const kSolutionSourceLocalMIP: i32 = 21;
/// rusolve mip_alns (Hendel 2022): sub-MIP solutions of the ALNS bandit, log letter `N`.
pub const kSolutionSourceAlns: i32 = 22;
/// rusolve mip_shift_and_propagate (Berthold, Hendel 2015): log letter `Q` (SCIP's
/// `P` is taken by Empty MIP).
pub const kSolutionSourceShiftProp: i32 = 23;
/// rusolve mip_proximity (Fischetti, Monaci 2014): proximity-search sub-MIP solutions, log letter `O`.
pub const kSolutionSourceProximity: i32 = 24;
/// rusolve mip_local_branching (Fischetti, Lodi 2003): local-branching sub-MIP solutions, log letter `K`.
pub const kSolutionSourceLocalBranching: i32 = 25;

impl HighsMipSolver {
    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getRow (lines 232-238)
    /// C++: `void getRow(HighsInt row, HighsInt & rowlen, const HighsInt *& rowinds, const double *& rowvals) const`
    pub fn get_row<'s>(
        &'s self,
        row: HighsInt,
        rowlen: &mut HighsInt,
        rowinds: &mut &'s [HighsInt],
        rowvals: &mut &'s [f64],
    ) {
        let start = self.mipdata_.a_rstart_[row as usize];
        *rowlen = self.mipdata_.a_rstart_[(row + 1) as usize] - start;
        let (s, e) = (start as usize, (start + *rowlen) as usize);
        *rowinds = &self.mipdata_.a_rindex_[s..e];
        *rowvals = &self.mipdata_.a_rvalue_[s..e];
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::terminatorActive (lines 255-255)
    /// C++: `bool terminatorActive() const`
    pub fn terminator_active(&self) -> bool {
        self.terminator_concurrency() > 0
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::parallelLockActive (lines 261-263)
    /// C++: `bool parallelLockActive() const`
    pub fn parallel_lock_active(&self) -> bool {
        self.mipdata_.parallel_lock && self.has_multiple_workers()
    }

    /// HiGHS 1.16 (#3124): HighsMipSolverData.h::updateWorkerEarlyTermination (lines 268-276)
    /// C++: `void updateWorkerEarlyTermination(HighsMipWorker & worker)`
    ///
    /// Lowers `worker_lp_iterations_stop` to worker `w`'s local LP iteration count (a
    /// relaxed CAS loop, as in C++: concurrent workers may race here) and marks the worker
    /// early-terminated. Called only under `h16_mip_workers`, from the task of worker `w`.
    pub fn update_worker_early_termination(&mut self, w: usize) {
        use std::sync::atomic::Ordering;
        let candidate = self.mipdata_.workers[w]
            .search_ptr_
            .as_deref()
            .expect("updateWorkerEarlyTermination: worker search is running (re-entrant call)")
            .lpiterations;
        let stop = &self.mipdata_.worker_lp_iterations_stop.0;
        let mut current = stop.load(Ordering::Relaxed);
        while candidate < current {
            match stop.compare_exchange_weak(
                current,
                candidate,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => current = actual,
            }
        }
        self.mipdata_.workers[w].early_termination = true;
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::hasMultipleWorkers (lines 265-265)
    /// C++: `bool hasMultipleWorkers() const`
    pub fn has_multiple_workers(&self) -> bool {
        self.mipdata_.workers.len() > 1
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getDomain (lines 267-267)
    /// C++: `HighsDomain & getDomain()`
    pub fn get_domain(&mut self) -> &mut HighsDomain {
        &mut self.mipdata_.domains[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getConflictPool (lines 268-268)
    /// C++: `HighsConflictPool & getConflictPool()`
    pub fn get_conflict_pool(&mut self) -> &mut HighsConflictPool {
        &mut self.mipdata_.conflictpools[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getCutPool (lines 269-269)
    /// C++: `HighsCutPool & getCutPool()`
    pub fn get_cut_pool(&mut self) -> &mut HighsCutPool {
        &mut self.mipdata_.cutpools[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getLp (lines 270-270)
    /// C++: `HighsLpRelaxation & getLp()`
    pub fn get_lp(&mut self) -> &mut HighsLpRelaxation {
        &mut self.mipdata_.lps[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getPseudoCost (lines 271-271)
    /// C++: `HighsPseudocost & getPseudoCost()`
    pub fn get_pseudo_cost(&mut self) -> &mut HighsPseudocost {
        &mut self.mipdata_.pseudocosts[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getDomain (lines 272-272)
    /// C++: `const HighsDomain & getDomain() const`
    pub fn get_domain_2(&self) -> &HighsDomain {
        &self.mipdata_.domains[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getConflictPool (lines 273-273)
    /// C++: `const HighsConflictPool & getConflictPool() const`
    pub fn get_conflict_pool_2(&self) -> &HighsConflictPool {
        &self.mipdata_.conflictpools[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getCutPool (lines 274-274)
    /// C++: `const HighsCutPool & getCutPool() const`
    pub fn get_cut_pool_2(&self) -> &HighsCutPool {
        &self.mipdata_.cutpools[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getLp (lines 275-275)
    /// C++: `const HighsLpRelaxation & getLp() const`
    pub fn get_lp_2(&self) -> &HighsLpRelaxation {
        &self.mipdata_.lps[0]
    }

    /// HiGHS: HighsMipSolverData.h::HighsMipSolverData::getPseudoCost (lines 276-276)
    /// C++: `const HighsPseudocost & getPseudoCost() const`
    pub fn get_pseudo_cost_2(&self) -> &HighsPseudocost {
        &self.mipdata_.pseudocosts[0]
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use super::*;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    #[test]
    fn s5_hook_source_constants() {
        assert_eq!(kSolutionSourceCount, 20);
        assert_eq!(kSolutionSourceCleanup, 19);
        assert_eq!(kSolutionSourceFPR, 20);
        assert_eq!(kSolutionSourceLocalMIP, 21);
        let ms = HighsMipSolver::default();
        assert_eq!(ms.solution_source_to_string(kSolutionSourceFPR, true), "A");
        assert_eq!(
            ms.solution_source_to_string(kSolutionSourceLocalMIP, false),
            "Local MIP"
        );
    }

    #[test]
    fn s5_hook_default_fields_off() {
        let data = HighsMipSolverData::default();
        assert_eq!(data.heuristic_effort_used, 0);
        assert!(data.s5_report.is_none());
    }

    #[test]
    fn s3_datafield_default_none() {
        let data = HighsMipSolverData::default();
        assert!(data.symmetry_dynamic.is_none());
    }

    #[test]
    fn s6_datafield_default_none() {
        let data = HighsMipSolverData::default();
        assert!(data.s6.is_none());
    }

    /// Items (weight row 1, weight row 2, value) of a two-row binary knapsack.
    type Items = &'static [(f64, f64, f64)];

    const ITEMS: Items = &[
        (23.0, 12.0, 92.0),
        (31.0, 45.0, 57.0),
        (29.0, 33.0, 49.0),
        (44.0, 21.0, 68.0),
        (53.0, 67.0, 60.0),
        (38.0, 14.0, 43.0),
        (63.0, 28.0, 67.0),
        (85.0, 35.0, 84.0),
        (89.0, 16.0, 87.0),
        (82.0, 40.0, 72.0),
        (17.0, 52.0, 31.0),
        (41.0, 27.0, 55.0),
    ];

    /// Six item types, each twice: the column pairs are symmetric, so symmetry detection
    /// spawns its task (HighsMipSolverData.cpp:553).
    const PAIRED_ITEMS: Items = &[
        (23.0, 12.0, 92.0),
        (23.0, 12.0, 92.0),
        (31.0, 45.0, 57.0),
        (31.0, 45.0, 57.0),
        (44.0, 21.0, 68.0),
        (44.0, 21.0, 68.0),
        (38.0, 14.0, 43.0),
        (38.0, 14.0, 43.0),
        (63.0, 28.0, 67.0),
        (63.0, 28.0, 67.0),
        (89.0, 16.0, 87.0),
        (89.0, 16.0, 87.0),
    ];

    const CAP: [f64; 2] = [165.0, 150.0];

    /// min -value.x, row1.x <= 165, row2.x <= 150, x binary.
    fn knapsack_lp(items: Items) -> HighsLp {
        let n = items.len();
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = 2;
        lp.col_cost_ = items.iter().map(|it| -it.2).collect();
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf; 2];
        lp.row_upper_ = CAP.to_vec();
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = (0..=n).map(|j| (2 * j) as HighsInt).collect();
        lp.a_matrix_.index_ = (0..n).flat_map(|_| [0, 1]).collect();
        lp.a_matrix_.value_ = items.iter().flat_map(|it| [it.0, it.1]).collect();
        lp
    }

    fn brute_force_optimum(items: Items) -> f64 {
        let mut best = 0.0f64;
        for mask in 0u32..(1 << items.len()) {
            let (mut a, mut b, mut v) = (0.0, 0.0, 0.0);
            for (j, it) in items.iter().enumerate() {
                if mask & (1 << j) != 0 {
                    a += it.0;
                    b += it.1;
                    v += it.2;
                }
            }
            if a <= CAP[0] && b <= CAP[1] && v > best {
                best = v;
            }
        }
        -best
    }

    /// What the root `TaskGroup` can change, bit for bit: status, objective, dual bound,
    /// nodes, the incumbent, `analyticCenterComputed`, the analytic centre and the number of
    /// symmetry generators.
    type Outcome = (
        HighsModelStatus,
        u64,
        u64,
        i64,
        Vec<u64>,
        bool,
        Vec<u64>,
        HighsInt,
    );

    /// Solve on a fresh thread with `threads` scheduler workers (the executor is per thread).
    fn solve_on_fresh_thread(threads: i32, items: Items, presolve: &'static str) -> Outcome {
        std::thread::spawn(move || {
            initialize_scheduler(threads);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.presolve = presolve.to_string();
            let lp = knapsack_lp(items);
            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let out = (
                ms.modelstatus_,
                ms.solution_objective_.to_bits(),
                ms.dual_bound_.to_bits(),
                ms.node_count_,
                ms.solution_.iter().map(|x| x.to_bits()).collect(),
                ms.mipdata_.analytic_center_computed,
                // `analyticCenter` is read only once `analyticCenterComputed` is set
                // (centralRounding after finishAnalyticCenterComputation). Before that a
                // thief may have written it (a stolen task that an early return's
                // ~TaskGroup waited for), as in C++: that value is unobservable.
                if ms.mipdata_.analytic_center_computed {
                    ms.mipdata_
                        .analytic_center
                        .iter()
                        .map(|x| x.to_bits())
                        .collect()
                } else {
                    Vec::new()
                },
                ms.mipdata_.symmetries.num_generators,
            );
            drop(ms);
            HighsTaskExecutor::shutdown(true);
            out
        })
        .join()
        .unwrap()
    }

    #[test]
    fn w5_mip_root_root_task_group() {
        // `HighsMipTaskGroup` is the real `highs::parallel::TaskGroup` (HighsMipSolverData.cpp:1995):
        // the AC (411) and symmetry (553) tasks are spawned on it, `sync` (497, 565) and the
        // restart `cancel` + `taskWait` (2131-2135, 2497-2502) act on it. At threads = 1 nothing
        // steals (serial LIFO order); at threads > 1 a thief may run the AC / symmetry task
        // early, and the result must not change (verdict T): each count runs twice.
        for items in [ITEMS, PAIRED_ITEMS] {
            let optimum = brute_force_optimum(items);
            for presolve in [kHighsOffString, kHighsOnString] {
                let case = format!("{} items, presolve {presolve}", items.len());
                let reference = solve_on_fresh_thread(1, items, presolve);
                assert_eq!(reference.0, HighsModelStatus::kOptimal, "{case}");
                assert_eq!(f64::from_bits(reference.1), optimum, "{case}");
                if presolve == kHighsOffString && items == ITEMS {
                    // Without presolve the distinct-item root reaches
                    // finishAnalyticCenterComputation (the paired set closes its gap and
                    // returns earlier, so its AC task is only cancelled and awaited).
                    assert!(reference.5, "{case}: analytic centre not computed");
                    assert_eq!(reference.6.len(), items.len(), "{case}");
                }
                for threads in [1, 2, 4, 8] {
                    for run in 0..2 {
                        let out = solve_on_fresh_thread(threads, items, presolve);
                        assert_eq!(out, reference, "{case}, threads {threads}, run {run}");
                    }
                }
            }
        }
    }
}
