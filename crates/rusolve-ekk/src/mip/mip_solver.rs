//! HiGHS source: mip/HighsMipSolver.h, mip/HighsMipSolver.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::lp_data::presolve_support::is_basis_consistent;
use crate::mip::search::plunge_levers_active;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g_width;
use std::sync::Arc;

/// HiGHS: HighsMipSolver.h::HighsTerminator (lines 22-34)
#[derive(Clone, Default, Debug)]
pub struct HighsTerminator {
    /// C++: `HighsInt num_instance`
    pub num_instance: HighsInt,
    /// C++: `HighsInt my_instance`
    pub my_instance: HighsInt,
    /// C++ `HighsModelStatus* record`: concurrency record (serial: one entry).
    pub record: Vec<HighsModelStatus>,
}

/// HiGHS: HighsMipSolver.h::HighsMipSolver (lines 36-169)
#[derive(Clone, Default, Debug)]
pub struct HighsMipSolver {
    /// C++ `HighsCallback* callback_`: owned; `Highs::callSolveMip` moves `Highs::callback_` in and back out (rule f).
    pub callback_: HighsCallback,
    /// C++ `const HighsOptions* options_mip_`: owned clone (never written through the pointer).
    pub options_mip_: HighsOptions,
    /// C++ `const HighsLp* model_`: selector, read through [`mip_model!`] (`&orig_model_` or `&mipdata_.presolvedModel`).
    pub model_: HighsMipModel,
    /// C++ `const HighsLp* orig_model_`: owned clone of the `lp` passed to [`HighsMipSolver::new`] (never written; `model_` selects it or `mipdata_.presolved_model`).
    pub orig_model_: HighsLp,
    /// C++: `HighsModelStatus modelstatus_`
    pub modelstatus_: HighsModelStatus,
    /// C++: `std::vector<double> solution_`
    pub solution_: Vec<f64>,
    /// C++: `double solution_objective_`
    pub solution_objective_: f64,
    /// C++: `double bound_violation_`
    pub bound_violation_: f64,
    /// C++: `double integrality_violation_`
    pub integrality_violation_: f64,
    /// C++: `double row_violation_`
    pub row_violation_: f64,
    /// C++: `double dual_bound_`
    pub dual_bound_: f64,
    /// C++: `double primal_bound_`
    pub primal_bound_: f64,
    /// C++: `double gap_`
    pub gap_: f64,
    /// C++: `int64_t node_count_`
    pub node_count_: i64,
    /// C++: `int64_t total_lp_iterations_`
    pub total_lp_iterations_: i64,
    /// C++: `double primal_dual_integral_`
    pub primal_dual_integral_: f64,
    /// C++ `FILE* improving_solution_file_` (`None` = nullptr). An `Arc` so that
    /// `Clone` copies the handle as C++ copies the pointer; `&File` is `Write`.
    pub improving_solution_file_: Option<std::sync::Arc<std::fs::File>>,
    /// C++: `std::vector<HighsObjectiveSolution> saved_objective_and_solution_`
    pub saved_objective_and_solution_: Vec<HighsObjectiveSolution>,
    /// C++: `bool submip`
    pub submip: bool,
    /// C++: `HighsInt submip_level`
    pub submip_level: HighsInt,
    /// C++: `HighsInt max_submip_level`
    pub max_submip_level: HighsInt,
    /// C++ `const HighsBasis* rootbasis` (sub-MIPs): owned; `solveSubMip` builds the
    /// basis as a local and moves it in (`None` = nullptr).
    pub rootbasis: Option<Box<HighsBasis>>,
    /// C++ `const HighsPseudocostInitialization* pscostinit` (sub-MIPs): owned; built by
    /// `solveSubMip` from the parent's pseudocost and moved in.
    pub pscostinit: Option<Box<HighsPseudocostInitialization>>,
    /// rusolve: stop a sub-MIP once this much work (`src/work.rs`) has been spent since it
    /// started; the LP-iteration limit of SCIP's proximity sub-MIP, in work units.
    /// `None` outside the proximity sub-MIP.
    pub work_limit: Option<u64>,
    /// mip_proximity: the work ledger and stall clock that survive MIP restarts
    pub(crate) prox: crate::mip::scip_heuristics::proximity::ProxLedger,
    /// C++ `const HighsCliqueTable* clqtableinit` (sub-MIPs): points at the parent's
    /// `mipdata_.cliquetable`, which stays in its slot (under parallel workers it is
    /// pinned and several sub-MIPs may read it at once); set by `solveSubMip` for the
    /// sub-MIP's `run()` and cleared right after it. The sub-MIP only reads it.
    pub clqtableinit: Option<ParentPtr<HighsCliqueTable>>,
    /// C++ `const HighsImplications* implicinit` (sub-MIPs): like `clqtableinit`.
    pub implicinit: Option<ParentPtr<HighsImplications>>,
    /// C++ `std::unique_ptr<HighsMipSolverData> mipdata_` (always set while solving).
    pub mipdata_: Box<HighsMipSolverData>,
    /// C++: `HighsModelStatus termination_status_`
    pub termination_status_: HighsModelStatus,
    /// C++: `HighsTerminator terminator_`
    pub terminator_: HighsTerminator,
    /// rusolve racing (Shinano et al. 2018): the racing store, `None` outside a race and in every sub-MIP.
    pub race: Option<Arc<crate::mip::race::RaceShared>>,
    /// rusolve: the racing store version this solver last imported from.
    pub race_seen: u64,
    /// C++: `HighsTimer timer_`
    pub timer_: HighsTimer,
    /// C++ `HighsProfiling* profiling_ = nullptr`: owned no-op stub (timers are stubbed).
    pub profiling_: HighsProfiling,
    /// reporting-only per-phase statistics (`mip_stats.rs`); the search never reads it.
    pub stats_: crate::mip::mip_stats::MipRunStats,
}

/// Conflict budget of one `mip_cdcl_dive` run.
const CDCL_DIVE_CONFLICTS: u64 = 50_000;

impl HighsMipSolver {
    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::HighsMipSolver (lines 31-70)
    /// C++: `HighsMipSolver(HighsCallback & callback, const HighsOptions & options, const HighsLp & lp, const HighsSolution & solution, bool submip = false, HighsInt submip_level = 0)`
    /// C++ default arguments: `submip = false`, `submip_level = 0`
    ///
    /// PORT NOTE (architecture): C++ stores `&callback`; here the callback is *moved*
    /// in (`mem::take(callback)`, crate rule f) and every caller moves it back after
    /// `run()`: `*callback = mem::take(&mut solver.callback_)` (`Highs::callSolveMip`,
    /// `solveSubMip`). `options` and `lp` are cloned (the C++ only reads them).
    pub fn new(
        callback: &mut HighsCallback,
        options: &HighsOptions,
        lp: &HighsLp,
        solution: &HighsSolution,
        submip: bool,
        submip_level: HighsInt,
    ) -> Self {
        let mut this = HighsMipSolver {
            callback_: std::mem::take(callback),
            options_mip_: options.clone(),
            model_: HighsMipModel::Orig,
            orig_model_: lp.clone(),
            solution_objective_: kHighsInf,
            submip,
            submip_level,
            rootbasis: None,
            pscostinit: None,
            work_limit: None,
            prox: Default::default(),
            clqtableinit: None,
            implicinit: None,
            ..Default::default()
        };
        // timer_.setPrintfFlag(output_flag, log_to_console): the stub timer prints nothing
        debug_assert!(!submip || submip_level > 0);
        this.max_submip_level = 0;
        // Initialise empty terminator
        this.initialise_terminator(0, kNoThreadInstance, None);
        debug_assert!(this.termination_status_ == HighsModelStatus::kNotset);
        if solution.value_valid {
            // (the NDEBUG-only assessLpPrimalSolution assert is not ported)
            // Initial solution can be infeasible, but need to set values for violation
            // and objective
            let mut quad_solution_objective_ = HighsCDouble::default();
            let (mut bv, mut rv, mut iv) = (0.0, 0.0, 0.0);
            this.solution_feasible(
                None,
                &solution.col_value,
                Some(&solution.row_value),
                &mut bv,
                &mut rv,
                &mut iv,
                &mut quad_solution_objective_,
            );
            this.bound_violation_ = bv;
            this.row_violation_ = rv;
            this.integrality_violation_ = iv;
            this.solution_objective_ = f64::from(quad_solution_objective_);
            this.solution_ = solution.col_value.clone();
        }
        this
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::~HighsMipSolver (lines 72-72)
    /// C++: `~HighsMipSolver()`
    ///
    /// `= default`: nothing to undo. The solver's pools, domains and LPs register only
    /// with each other, inside `mipdata_`, and go with it. The caller must have moved
    /// the callback back first (see [`HighsMipSolver::new`]).
    pub fn release(&mut self) {}

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::runTask (lines 74-93)
    /// C++: `template <class F> void runTask(F&& f, TaskGroup& tg, bool parallel_lock, bool force_serial, const std::vector<HighsInt>& indices = {0})`
    ///
    /// PORT NOTE: `tg` is the `TaskGroup` that `run()` declares at C++ 248 and uses only
    /// here; every `runTask` ends with `tg.taskWait()`, so its deque head is the head at
    /// entry to each call. The port opens a scoped group (rule (g)) on this thread's deque
    /// for the spawn branch instead: the same head, the same spawns and the same
    /// `taskWait`. Its destructor's `cancel` finds no task left.
    ///
    /// `f` is `Fn + Sync`: C++ calls the one lambda `f` from every spawned task at once
    /// (`[&f, i] { f(i); }`), so the callers write per-index state through `SharedMut`
    /// slots (rule (s)) rather than through a `FnMut` capture.
    pub fn run_task(
        &mut self,
        f: &(dyn Fn(&mut HighsMipSolver, HighsInt) + Sync),
        parallel_lock: bool,
        force_serial: bool,
        indices: &[HighsInt],
    ) {
        if indices.is_empty() {
            return;
        }
        self.set_parallel_lock(parallel_lock);
        let spawn_tasks = !force_serial
            && indices.len() > 1
            && !self.options_mip_.mip_search_simulate_concurrency;
        // C++ 81-91: one loop with the `spawn_tasks` test inside, then `taskWait` when
        // spawning. `spawn_tasks` is loop-invariant, so the loop is split by branch; the
        // calls and spawns are in index order either way.
        if spawn_tasks {
            self.pin_shared_slots(true);
            // C++ captures `this` and `f` by reference in each spawned lambda (cpp 84).
            let ms_shared = crate::parallel::shared_mut::SharedMut::new(&mut *self);
            crate::parallel::highs_parallel::task_group(|tg| {
                for &i in indices {
                    let ms_ref = &ms_shared;
                    tg.spawn(move || {
                        // SAFETY: MIP worker protocol under `setParallelLock` (C++ 80):
                        // task `i` touches worker `i`'s own LP, search, domain and pools,
                        // and shared solver state only as the C++ does under the parallel
                        // lock. The scoped group (rule (g)) joins every task before
                        // `ms_shared` is dropped.
                        let ms: &mut HighsMipSolver = unsafe { ms_ref.get_mut_unchecked() };
                        f(ms, i);
                    });
                }
                tg.task_wait();
            });
            self.pin_shared_slots(false);
        } else {
            for &i in indices {
                f(self, i);
            }
        }
        self.set_parallel_lock(false);
    }

    /// Parallel glue (no C++ counterpart): pin, for the concurrent tasks of [`Self::run_task`],
    /// the slots of the objects all C++ workers reach through `mipdata_` (the singletons
    /// and the global `getDomain()` / `getCutPool()` / `getConflictPool()` /
    /// `getPseudoCost()`), so that every task uses the one object in place instead of
    /// moving it out of its slot ([`Slot::set_pinned`] has the C++ protocol). The other
    /// deque entries belong to one worker each (`createNewWorkers` /
    /// `constructAdditionalWorkerData`) and stay move-detached by their worker.
    ///
    /// Before pinning, the global domain replays its pending pool events (M3) serially:
    /// replay writes the domain, and the tasks may all enter `getDomain()` methods. The
    /// replay is exact at any point between two domain methods, and no task appends to
    /// the global pools' logs (workers add cuts and conflicts to their own pools), so
    /// the pinned domain has nothing left to replay while the tasks run.
    fn pin_shared_slots(&mut self, pinned: bool) {
        if pinned
            && !self.mipdata_.domains.is_empty()
            && self.mipdata_.domains[0].has_pending_pool_events(self)
        {
            detach!(self.mipdata_.domains[0], |d| d.replay_pool_events(self));
        }
        let data = &mut *self.mipdata_;
        data.heuristics.set_pinned(pinned);
        data.cliquetable.set_pinned(pinned);
        data.implications.set_pinned(pinned);
        data.redcostfixing.set_pinned(pinned);
        data.objective_function.set_pinned(pinned);
        data.nodequeue.set_pinned(pinned);
        if let Some(slot) = data.domains.first_mut() {
            slot.set_pinned(pinned);
        }
        if let Some(slot) = data.cutpools.first_mut() {
            slot.set_pinned(pinned);
        }
        if let Some(slot) = data.conflictpools.first_mut() {
            slot.set_pinned(pinned);
        }
        if let Some(slot) = data.pseudocosts.first_mut() {
            slot.set_pinned(pinned);
        }
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::run (lines 95-1023)
    /// C++: `void run()`
    ///
    /// PORT NOTE: the C++ lambdas are private methods (`run_*`) taking the worker index
    /// (`usize`) and the [`RunState`] of the restart iteration. `goto restart` is
    /// `continue 'restart`: every per-restart local is declared inside the loop body.
    /// Only the presolve/solve/postsolve clocks are live (`profiling_.mip_()` is false, so the
    /// timing logs stay silent), and the C++ `assert` calls are omitted (the oracle is NDEBUG).
    pub fn run(&mut self) {
        self.modelstatus_ = HighsModelStatus::kNotset;
        self.max_submip_level = std::cmp::max(self.submip_level, self.max_submip_level);
        // Start the timer local to HighsMipSolver
        self.timer_.start();
        // reporting-only.
        // Per-solve counters: this solve's scope (or the enclosing solve's, for a sub-MIP).
        let work_scope = crate::work::SolveScope::enter();
        self.stats_.work_counters = Some(work_scope.counters().clone());
        self.stats_.work_start = work_scope.counters().snapshot();
        // `fopen(.., "w")`: create or truncate; a failed open is nullptr (no file output)
        self.improving_solution_file_ =
            if !self.submip && !self.options_mip_.mip_improving_solution_file.is_empty() {
                std::fs::File::create(&self.options_mip_.mip_improving_solution_file)
                    .ok()
                    .map(std::sync::Arc::new)
            } else {
                None
            };

        // mipdata_ = new HighsMipSolverData(*this); init(); runMipPresolve(...) (inline in C++)
        self.new_data();
        self.profiling_.start(kPresolveTime, false);
        self.init();
        // HIGHS_DEBUGSOL (mip_debug_solution_file): HighsMipSolver.cpp lines 117-121
        HighsDebugSol::activate(self);
        self.mipdata_.debug_solution.suspend();
        self.run_mip_presolve_data(self.options_mip_.presolve_reduction_limit);
        self.profiling_.stop(kPresolveTime);
        self.stats_.presolve_end_s = self.timer_.read(); // C5
        if self.profiling_.mip_() && !self.submip {
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "MIP-Timing: {} - completed presolve\n",
                c_fmt_g_width(self.timer_.read(), 2, 11)
            );
        }
        // HiGHS 1.16 #3578 (h16_mip_workers): HighsMipSolver.cpp run lines 136-156 (1.16.0)
        if self.options_mip_.h16_mip_workers {
            if self.modelstatus_ == HighsModelStatus::kNotset {
                let presolve_status = self.mipdata_.presolve_status;
                // Time limit has been reached (in presolve)
                if presolve_status == HighsPresolveStatus::kTimeout {
                    self.modelstatus_ = HighsModelStatus::kTimeLimit;
                } else if presolve_status == HighsPresolveStatus::kOutOfMemory {
                    // Handle failures in presolve
                    // PORT NOTE: 1.16 also maps `kException` and `kMatrixValueError` to
                    // `kPresolveError` here; neither presolve status is in the port's
                    // enum (the port's presolve never raises them), so only the
                    // out-of-memory branch remains.
                    self.modelstatus_ = HighsModelStatus::kMemoryLimit;
                }
            }
        } else if self.modelstatus_ == HighsModelStatus::kNotset
            && self.timer_.read() >= self.options_mip_.time_limit
        {
            // Identify whether time limit has been reached (in presolve)
            self.modelstatus_ = HighsModelStatus::kTimeLimit;
        }

        if self.modelstatus_ != HighsModelStatus::kNotset {
            // HiGHS 1.16 #3135 (h16_mip_start): read the option once; off = the 1.15.1 statements
            let h16_mip_start = self.options_mip_.h16_mip_start;
            if h16_mip_start {
                // HiGHS 1.16 #3135: HighsMipSolver.cpp run lines 158-168 (1.16.0)
                let mut presolve_status = self.modelstatus_;
                if presolve_status == HighsModelStatus::kInfeasible
                    && self.solution_objective_ != kHighsInf
                    && self.bound_violation_ <= self.options_mip_.mip_feasibility_tolerance
                    && self.integrality_violation_ <= self.options_mip_.mip_feasibility_tolerance
                    && self.row_violation_ <= self.options_mip_.mip_feasibility_tolerance
                {
                    presolve_status = HighsModelStatus::kOptimal;
                }
                crate::highs_log_user!(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "Presolve: {}\n",
                    util_model_status_to_string(presolve_status)
                );
            } else {
                crate::highs_log_user!(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "Presolve: {}\n",
                    util_model_status_to_string(self.modelstatus_)
                );
            }
            if self.modelstatus_ == HighsModelStatus::kOptimal {
                self.mipdata_.lower_bound = 0.0;
                self.mipdata_.upper_bound = 0.0;
                if h16_mip_start {
                    // HiGHS 1.16 #3272 (h16_mip_start): HighsMipSolver.cpp run lines 148-190 (1.16.0)
                    // Model has been reduced to empty in presolve, so optimality
                    // has been deduced
                    let mut solution_objective =
                        self.transform_new_integer_feasible_solution(&[], true);
                    if solution_objective == kHighsInf {
                        // Solution is not feasible (due to an error in presolve) so
                        // don't return a model status of HighsModelStatus::kOptimal!
                        self.modelstatus_ = HighsModelStatus::kSolveError;
                    } else {
                        // Set solution_objective = -kHighsInf, which ensures that the
                        // incumbent is viewed as an improving solution in
                        // saveReportMipSolution
                        solution_objective = -kHighsInf;
                    }
                    // Pass solution_objective as new_upper_limit
                    self.save_report_mip_solution(solution_objective);
                } else {
                    self.transform_new_integer_feasible_solution(&[], true);
                    self.save_report_mip_solution(-kHighsInf);
                }
            }
            self.cleanup_solve();
            return;
        }

        if self.profiling_.mip_() && !self.submip {
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "MIP-Timing: {} - starting  setup\n",
                c_fmt_g_width(self.timer_.read(), 2, 11)
            );
        }
        // HIGHS_DEBUGSOL (mip_debug_solution_file): HighsMipSolver.cpp lines 153-155
        self.mipdata_.debug_solution.resume();
        self.profiling_.start(kSolveTime, false);
        self.run_setup();
        self.stats_.setup_end_s = self.timer_.read(); // C5
        if self.profiling_.mip_() && !self.submip {
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "MIP-Timing: {} - completed setup\n",
                c_fmt_g_width(self.timer_.read(), 2, 11)
            );
        }

        if self.get_domain_2().infeasible() {
            self.cleanup_solve();
            return;
        }
        // Initialise master worker.
        let master = HighsMipWorker::new(self, 0, 0, 0, 0, 0);
        self.mipdata_.workers.push(Slot::new(master));

        'restart: loop {
            if self.modelstatus_ == HighsModelStatus::kNotset {
                // Check limits have not been reached before evaluating root node
                if self.check_limits(0) {
                    self.cleanup_solve();
                    return;
                }
                // Possibly query existence of an external solution
                if !self.submip {
                    self.query_external_solution(
                        self.solution_objective_,
                        kExternalMipSolutionQueryOriginAfterSetup,
                    );
                }

                // Apply the trivial heuristics
                let returned_model_status = self.trivial_heuristics();
                if self.modelstatus_ == HighsModelStatus::kNotset
                    && returned_model_status == HighsModelStatus::kInfeasible
                {
                    // trivialHeuristics can spot trivial infeasibility, so act on it
                    self.modelstatus_ = returned_model_status;
                    self.cleanup_solve();
                    return;
                }
                // Apply the feasibility jump heuristic (if enabled)
                if self.options_mip_.mip_heuristic_run_feasibility_jump {
                    let returned_model_status = self.feasibility_jump();
                    if self.modelstatus_ == HighsModelStatus::kNotset
                        && returned_model_status == HighsModelStatus::kInfeasible
                    {
                        // feasibilityJump can spot trivial infeasibility, so act on it
                        self.modelstatus_ = returned_model_status;
                        self.cleanup_solve();
                        return;
                    }
                }
                // rusolve mip_first_inc: no incumbent after the FJ
                if self.options_mip_.mip_first_inc
                    && !self.submip
                    && self.mipdata_.upper_limit == kHighsInf
                {
                    self.first_incumbent_fallback();
                }
                // rusolve mip_shift_and_propagate (Berthold, Hendel 2015): SCIP runs it
                // before the root node (HEUR_TIMING BEFORENODE, freq 0), only while no solution is known
                if self.options_mip_.mip_shift_and_propagate && !self.submip {
                    crate::mip::scip_heuristics::shiftprop::shift_and_propagate(self);
                }
                // rusolve mip_zero_obj_cdcl
                self.cdcl_zero_obj();
                // rusolve mip_cdcl_dive: LP-free conflict-driven DFS
                if (self.options_mip_.mip_cdcl_dive || self.options_mip_.mip_cdcl_uip)
                    && self.mipdata_.upper_limit == kHighsInf
                {
                    self.cdcl_dive(CDCL_DIVE_CONFLICTS);
                }
                // mip_heuristic_run_fpr / mip_heuristic_run_local_mip: mip-heuristics apply_patch.cmake lines 1380-1382 at 36f7d9d
                let s5_on = self.options_mip_.mip_heuristic_run_fpr
                    || self.options_mip_.mip_heuristic_run_local_mip;
                if s5_on && crate::mip::fpr_local_mip::driver::s5_run_presolve(self) {
                    self.modelstatus_ = HighsModelStatus::kInfeasible;
                    self.cleanup_solve();
                    return;
                }
                // rusolve mip_improve_early: the improvement chain on the
                // feasibility-jump incumbent, before the root LP
                if self.options_mip_.mip_improve_early && !self.submip {
                    crate::mip::scip_heuristics::proximity::improve_early(self);
                }
                // End of pre-root-node heuristics
                if self.profiling_.mip_() && !self.submip {
                    crate::highs_log_user!(
                        &self.options_mip_.log_options,
                        HighsLogType::kInfo,
                        "MIP-Timing: {} - starting evaluate root node\n",
                        c_fmt_g_width(self.timer_.read(), 2, 11)
                    );
                }
                let c5_root_start = self.timer_.read(); // C5
                self.evaluate_root_node(0);
                self.stats_.root_end_s = self.timer_.read();
                self.stats_.root_node_s += self.stats_.root_end_s - c5_root_start;
                if self.terminate() {
                    self.modelstatus_ = self.termination_status();
                    self.cleanup_solve();
                    return;
                }
                if (self.options_mip_.mip_cdcl_dive || self.options_mip_.mip_cdcl_uip)
                    && self.mipdata_.upper_limit == kHighsInf
                {
                    self.cdcl_dive(CDCL_DIVE_CONFLICTS);
                }
                // mip_heuristic_run_oneopt: SCIP heur_oneopt.c lines 422-437 at 32dd425
                if self.options_mip_.mip_heuristic_run_oneopt {
                    crate::mip::scip_heuristics::driver::s6_run_root_oneopt(self);
                }
                // rusolve mip_alns: the one call right after the root
                if self.options_mip_.mip_alns {
                    crate::mip::scip_heuristics::alns::alns_step(self, true);
                }
                crate::mip::scip_heuristics::proximity::fischetti_step(self);
                if self.profiling_.mip_() && !self.submip {
                    crate::highs_log_user!(
                        &self.options_mip_.log_options,
                        HighsLogType::kInfo,
                        "MIP-Timing: {} - completed evaluate root node\n",
                        c_fmt_g_width(self.timer_.read(), 2, 11)
                    );
                }
                // age 5 times to remove stored but never violated cuts after root separation
                for _ in 0..5 {
                    self.get_cut_pool().perform_aging();
                }
            }
            if self.mipdata_.nodequeue.empty() || self.check_limits(0) {
                self.cleanup_solve();
                return;
            }

            let best_lower_bound = self.mipdata_.nodequeue.get_best_lower_bound();
            self.update_lower_bound(best_lower_bound, true, true);
            self.print_display_line(kSolutionSourceNone);

            let max_num_workers = self.get_max_num_workers();
            let mut num_workers: HighsInt = 1;

            // Ensure the master worker's search is set up from the current global state
            {
                let ub = self.mipdata_.upper_bound;
                let ul = self.mipdata_.upper_limit;
                let ol = self.mipdata_.optimality_limit;
                let nc = self.num_col();
                detach!(self.mipdata_.workers[0], |w| {
                    w.reset_search(self);
                    w.reset_sepa(self);
                });
                self.mipdata_.workers[0].nodequeue.clear();
                self.mipdata_.workers[0].nodequeue.set_num_col(nc);
                self.mipdata_.workers[0].upper_bound = ub;
                self.mipdata_.workers[0].upper_limit = ul;
                self.mipdata_.workers[0].optimality_limit = ol;
            }
            self.run_destroy_old_workers();
            {
                let local = self.run_worker_search(0).get_local_domain_2();
                self.mipdata_.debug_solution.register_domain(local);
            }

            // HiGHS 1.16 #3124 #3192 #3336 #3406 (h16_mip_workers): HighsMipSolver.cpp run lines
            // 497-1201 (1.16.0): the 1.16 search loop. Off = the 1.15.1 loop below.
            if self.options_mip_.h16_mip_workers {
                if self.run_search_h16(max_num_workers) {
                    continue 'restart;
                }
                self.run_sync_solutions();
                self.cleanup_solve();
                return;
            }

            let mut st = RunState::new(self);
            let mut search_indices: Vec<usize> = vec![0];
            let mut worker_restart_votes: Vec<RestartVote> =
                vec![RestartVote::NoCheck; max_num_workers as usize];
            let mut root_node = true; // Don't separate the root node again
            // for ramp up
            let mut node_lim: HighsInt = if max_num_workers > 1 { 1 } else { kHighsIInf };

            while !self.mipdata_.nodequeue.empty() {
                // Possibly query existence of an external solution
                if !self.submip {
                    self.query_external_solution(
                        self.solution_objective_,
                        kExternalMipSolutionQueryOriginBeforeDive,
                    );
                }
                self.race_import();

                // S6 (mip_heuristic_run_oneopt|crossover): SCIP heur_oneopt.c lines 472-889, heur_crossover.c lines 903-1044 at 32dd425
                let s6_on = self.options_mip_.mip_heuristic_run_oneopt
                    || self.options_mip_.mip_heuristic_run_crossover;
                if s6_on {
                    crate::mip::scip_heuristics::driver::s6_run_tree_heuristics(self);
                }
                // rusolve mip_alns (Hendel 2022): sibling of the SCIP tree-heuristics hook
                if self.options_mip_.mip_alns {
                    crate::mip::scip_heuristics::alns::alns_step(self, false);
                }
                crate::mip::scip_heuristics::proximity::fischetti_step(self);

                // Update global pseudo-cost with worker information
                self.run_sync_global_pseudo_cost();

                // Get new candidate worker search indices
                self.run_get_search_indices_with_no_nodes(
                    &mut search_indices,
                    &mut worker_restart_votes,
                    num_workers,
                );

                // Only update worker's pseudo-costs that have been assigned a node
                self.run_reset_worker_pseudo_costs(&search_indices);

                // Assign nodes to workers
                let mut limit_reached = false;
                if root_node {
                    let node = self.mipdata_.nodequeue.pop_best_bound_node();
                    HighsSearch::with_worker_search(self, 0, |s, ctx| s.install_node(ctx, node));
                } else {
                    self.run_install_nodes(&search_indices, &mut limit_reached, &mut st);
                }
                if limit_reached {
                    break;
                }

                if node_lim != kHighsIInf {
                    if num_workers >= max_num_workers {
                        node_lim = std::cmp::max(20, 2 * node_lim);
                    }
                    if node_lim > 100 {
                        node_lim = kHighsIInf;
                    }
                }

                // Process nodes (separation / heuristics / dives)
                let avgiter = self.get_lp().get_avg_solve_iters();
                self.run_process_nodes(
                    &search_indices,
                    &mut worker_restart_votes,
                    root_node,
                    node_lim,
                    100,
                    avgiter,
                    &st,
                );

                root_node = false;

                // Sync statistics, check infeasibility, and flush nodes from worker queues
                let mut infeasible = false;
                for &w in &search_indices {
                    let gd = self.mipdata_.workers[w].get_global_domain();
                    if self.mipdata_.domains[gd].infeasible() {
                        infeasible = true;
                    }
                    HighsSearch::with_worker_search(self, w, |s, ctx| {
                        s.open_nodes_to_queue(ctx, HighsSearch_NodeQueueRef::Global)
                    });
                    // rusolve mip_symmetry_dynamic: carry the branching path on the queue transfer
                    while self.mipdata_.workers[w].nodequeue.num_nodes() > 0 {
                        let node = self.mipdata_.workers[w].nodequeue.pop_best_node();
                        // rusolve mip_node_install_basis: and the parent basis
                        if self.options_mip_.mip_symmetry_dynamic || node.basis.is_some() {
                            let _ = self.mipdata_.nodequeue.emplace_node_ext(
                                node.domchgstack,
                                node.branchings,
                                node.lower_bound,
                                node.estimate,
                                node.depth,
                                node.s3_path,
                                node.basis,
                            );
                        } else {
                            let _ = self.mipdata_.nodequeue.emplace_node(
                                node.domchgstack,
                                node.branchings,
                                node.lower_bound,
                                node.estimate,
                                node.depth,
                            );
                        }
                    }
                    self.run_with_search_only(w, |s, ms| s.flush_statistics(ms));
                    self.run_sync_sepa_stats(w);
                    detach!(self.mipdata_.heuristics, |h| h.flush_statistics(self, w));
                }

                if infeasible {
                    self.mipdata_.nodequeue.clear();
                    self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                    let ub = self.mipdata_.upper_bound;
                    self.update_lower_bound(cpp_min(kHighsInf, ub), true, true);
                    break;
                }

                let best = self.mipdata_.nodequeue.get_best_lower_bound();
                let ub = self.mipdata_.upper_bound;
                self.update_lower_bound(cpp_min(ub, best), true, true);

                self.run_sync_solutions();

                limit_reached = self.check_limits(0);
                if limit_reached {
                    self.print_display_line(kSolutionSourceNone);
                    break;
                }

                // Sync global information
                self.run_sync_pools(&search_indices);
                self.run_sync_global_domain(&search_indices);
                detach!(self.mipdata_.domains[0], |d| {
                    d.propagate(self);
                });

                let feastol = self.mipdata_.feastol;
                let pruned = detach!(self.mipdata_.nodequeue, |nq| {
                    detach!(self.mipdata_.domains[0], |d| nq
                        .prune_infeasible_nodes(self, d, feastol))
                });
                self.mipdata_.pruned_treeweight += pruned;

                if self.get_domain_2().infeasible() {
                    self.mipdata_.nodequeue.clear();
                    self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                    let ub = self.mipdata_.upper_bound;
                    self.update_lower_bound(cpp_min(kHighsInf, ub), true, true);
                    self.print_display_line(kSolutionSourceNone);
                    break;
                }

                let best = self.mipdata_.nodequeue.get_best_lower_bound();
                let ub = self.mipdata_.upper_bound;
                self.update_lower_bound(cpp_min(ub, best), true, true);
                self.print_display_line(kSolutionSourceNone);

                if self.mipdata_.nodequeue.empty() {
                    break;
                }

                // Reset global domain and sync worker's global domains.
                let spawn_more_workers = num_workers < max_num_workers
                    && self.mipdata_.nodequeue.num_nodes() > num_workers as i64;
                let reset_workers = self.has_multiple_workers();
                self.run_reset_global_domain(spawn_more_workers, reset_workers, num_workers);

                if node_lim == kHighsIInf
                    && self.run_check_worker_restart_votes(&worker_restart_votes, &mut st)
                {
                    continue 'restart;
                }

                if spawn_more_workers {
                    let new_max_num_workers = std::cmp::min(
                        self.mipdata_.nodequeue.num_nodes() as HighsInt,
                        max_num_workers,
                    );
                    self.get_pseudo_cost().remove_changed();
                    if num_workers == 1 {
                        self.run_construct_additional_worker_data();
                    }
                    self.run_create_new_workers(new_max_num_workers - num_workers);
                    num_workers = new_max_num_workers;
                }
            }
            self.run_sync_solutions();

            self.cleanup_solve();
            return;
        }
    }

    /// Worker `w`'s search, for non-context calls (C++ `*workers[w].search_ptr_`).
    fn run_worker_search(&self, w: usize) -> &HighsSearch {
        self.mipdata_.workers[w]
            .search_ptr_
            .as_deref()
            .expect("HighsMipSolver::run: worker search is running (re-entrant call)")
    }

    /// Mutable form of [`HighsMipSolver::run_worker_search`].
    fn run_worker_search_mut(&mut self, w: usize) -> &mut HighsSearch {
        self.mipdata_.workers[w]
            .search_ptr_
            .as_deref_mut()
            .expect("HighsMipSolver::run: worker search is running (re-entrant call)")
    }

    /// Run `f` on worker `w`'s search, which is taken out of the worker for the call and
    /// put back afterwards (for the `flushStatistics(*this)` style calls that take `ms`).
    fn run_with_search_only<R>(
        &mut self,
        w: usize,
        f: impl FnOnce(&mut HighsSearch, &mut HighsMipSolver) -> R,
    ) -> R {
        let mut search = self.mipdata_.workers[w]
            .search_ptr_
            .take()
            .expect("HighsMipSolver::run: worker search is running (re-entrant call)");
        let result = f(&mut search, self);
        self.mipdata_.workers[w].search_ptr_ = Some(search);
        result
    }

    /// C++ `parallelLockActive() ? workers[i].nodequeue : mipdata_->nodequeue` (the
    /// `globalqueue` choice in `evaluateNode`, `separateAndStoreBasis`, `backtrackPlunge`).
    fn run_queue_ref(&self) -> HighsSearch_NodeQueueRef {
        if self.parallel_lock_active() {
            HighsSearch_NodeQueueRef::Worker
        } else {
            HighsSearch_NodeQueueRef::Global
        }
    }

    /// C++ `destroyOldWorkers` (lambda in run()). Popped domains and workers are released
    /// first (their destructors unregister from the pools, which `Slot` drop does not do).
    fn run_destroy_old_workers(&mut self) {
        if self.mipdata_.workers.len() <= 1 {
            return;
        }
        while self.mipdata_.domains.len() > 1 {
            let mut slot = self.mipdata_.domains.pop().expect("domains non-empty");
            let mut domain = slot.take();
            domain.release(self);
        }
        while self.mipdata_.lps.len() > 1 {
            self.mipdata_.lps.pop();
        }
        while self.mipdata_.pseudocosts.len() > 1 {
            self.mipdata_.pseudocosts.pop();
        }
        while self.mipdata_.workers.len() > 1 {
            let mut slot = self.mipdata_.workers.pop().expect("workers non-empty");
            let mut worker = slot.take();
            worker.release(self);
        }
        while self.mipdata_.cutpools.len() > 1 {
            self.mipdata_.cutpools.pop();
        }
        while self.mipdata_.conflictpools.len() > 1 {
            self.mipdata_.conflictpools.pop();
        }
    }

    /// C++ `createNewWorkers` (lambda in run()).
    fn run_create_new_workers(&mut self, num_new_workers: HighsInt) {
        if num_new_workers <= 0 {
            return;
        }
        let nc = self.num_col();
        let age = self.options_mip_.mip_pool_age_limit;
        let soft = self.options_mip_.mip_pool_soft_limit;
        // HiGHS 1.16 #3574 (h16_cutpool_age): HighsMipSolver.cpp run_create_new_workers (1.16.0)
        let h16_cutpool_age = self.options_mip_.h16_cutpool_age;
        // Remove all cuts from non-global pool for copied LP
        let lp_copy = HighsLpRelaxation::copy_from(self, &self.mipdata_.lps[0]);
        self.mipdata_.lps.push(Slot::new(lp_copy));
        let last_lp = self.mipdata_.lps.len() - 1;
        detach!(self.mipdata_.lps[last_lp], |lp| {
            lp.remove_worker_specific_rows(self);
        });
        for i in 0..num_new_workers {
            if i != 0 {
                let last = self.mipdata_.lps.len() - 1;
                let copy = HighsLpRelaxation::copy_from(self, &self.mipdata_.lps[last]);
                self.mipdata_.lps.push(Slot::new(copy));
            }
            let domain_copy = detach!(self.mipdata_.domains[0], |g| HighsDomain::copy_from(
                self, g
            ));
            self.mipdata_.domains.push(Slot::new(domain_copy));
            let cp = self.mipdata_.cutpools.len();
            self.mipdata_
                .cutpools
                .push(Slot::new(HighsCutPool::new_with_age_flag(
                    nc,
                    age,
                    soft,
                    cp as HighsInt,
                    h16_cutpool_age,
                )));
            let cf = self.mipdata_.conflictpools.len();
            self.mipdata_
                .conflictpools
                .push(Slot::new(HighsConflictPool::new(5 * age, soft)));
            let dl = self.mipdata_.domains.len() - 1;
            detach!(self.mipdata_.domains[dl], |d| {
                d.add_cutpool(self, cp);
                d.add_conflict_pool(self, cf);
            });
            let pc = self.mipdata_.pseudocosts.len();
            let pseudocost = HighsPseudocost::new(self);
            self.mipdata_.pseudocosts.push(Slot::new(pseudocost));
            let lp = self.mipdata_.lps.len() - 1;
            let worker = HighsMipWorker::new(self, lp, dl, cp, cf, pc);
            self.mipdata_.workers.push(Slot::new(worker));
            let wi = self.mipdata_.workers.len() - 1;
            self.mipdata_.lps[lp].set_mip_worker(wi);
            detach!(self.mipdata_.lps[lp], |l| {
                l.notify_cut_pools_lp_copied(self, 1);
            });
            // C++ `random_seed + workers.size() - 1`
            let seed = (self.options_mip_.random_seed as HighsUInt).wrapping_add(wi as HighsUInt);
            self.mipdata_.workers[wi].randgen.initialise(seed);
            self.mipdata_.workers[wi].nodequeue.set_num_col(nc);
            let local = self.run_worker_search(wi).get_local_domain_2();
            self.mipdata_.debug_solution.register_domain(local);
        }
    }

    /// C++ `constructAdditionalWorkerData` (lambda in run()): gives master worker 0
    /// its own cut pool, conflict pool, global domain and pseudocost.
    fn run_construct_additional_worker_data(&mut self) {
        let nc = self.num_col();
        let age = self.options_mip_.mip_pool_age_limit;
        let soft = self.options_mip_.mip_pool_soft_limit;
        // HiGHS 1.16 #3574 (h16_cutpool_age): HighsMipSolver.cpp run_construct_additional_worker_data (1.16.0)
        let h16_cutpool_age = self.options_mip_.h16_cutpool_age;
        let cp = self.mipdata_.cutpools.len();
        self.mipdata_
            .cutpools
            .push(Slot::new(HighsCutPool::new_with_age_flag(
                nc,
                age,
                soft,
                1,
                h16_cutpool_age,
            )));
        self.mipdata_.workers[0].set_cut_pool(cp);
        let cf = self.mipdata_.conflictpools.len();
        self.mipdata_
            .conflictpools
            .push(Slot::new(HighsConflictPool::new(5 * age, soft)));
        self.mipdata_.workers[0].set_conflict_pool(cf);
        let domain_copy = detach!(self.mipdata_.domains[0], |g| HighsDomain::copy_from(
            self, g
        ));
        self.mipdata_.domains.push(Slot::new(domain_copy));
        let dl = self.mipdata_.domains.len() - 1;
        self.mipdata_.workers[0].set_global_domain(dl);
        detach!(self.mipdata_.domains[dl], |d| {
            d.add_cutpool(self, cp);
            d.add_conflict_pool(self, cf);
        });
        let pc = self.mipdata_.pseudocosts.len();
        let pseudocost = HighsPseudocost::new(self);
        self.mipdata_.pseudocosts.push(Slot::new(pseudocost));
        self.mipdata_.workers[0].set_pseudocost(pc);
        let lp = self.mipdata_.workers[0].get_lp_relaxation();
        self.mipdata_.lps[lp].set_mip_worker(0);
        detach!(self.mipdata_.workers[0], |w| {
            w.reset_search(self);
            w.reset_sepa(self);
        });
        self.mipdata_.workers[0].nodequeue.clear();
        self.mipdata_.workers[0].nodequeue.set_num_col(nc);
    }

    /// C++ `syncSolutions` (lambda in run()).
    fn run_sync_solutions(&mut self) {
        for w in 0..self.mipdata_.workers.len() {
            // Note: Upper bound / limit of workers updated via addIncumbent
            let solutions = std::mem::take(&mut self.mipdata_.workers[w].solutions_);
            for (sol, solobj, solution_source) in solutions {
                self.add_incumbent(&sol, solobj, solution_source, true, false);
            }
        }
    }

    /// C++ `syncPools` (lambda in run()).
    fn run_sync_pools(&mut self, indices: &[usize]) {
        if !self.has_multiple_workers() || self.parallel_lock_active() {
            return;
        }
        for &w in indices {
            let cf = self.mipdata_.workers[w].get_conflict_pool();
            detach!(self.mipdata_.conflictpools[cf], |pool| {
                detach!(self.mipdata_.conflictpools[0], |global| {
                    pool.sync_conflict_pool(global);
                });
            });
            let cp = self.mipdata_.workers[w].get_cut_pool();
            detach!(self.mipdata_.cutpools[cp], |pool| {
                detach!(self.mipdata_.cutpools[0], |global| {
                    pool.sync_cut_pool(self, global);
                });
            });
        }
        self.mipdata_.cutpools[0].perform_aging();
        self.mipdata_.conflictpools[0].perform_aging(false);
    }

    /// C++ `syncGlobalDomain` (lambda in run()).
    fn run_sync_global_domain(&mut self, indices: &[usize]) {
        if !self.has_multiple_workers() {
            return;
        }
        for &w in indices {
            let gd = self.mipdata_.workers[w].get_global_domain();
            // PORT NOTE: the worker's stack is cloned; changeBound below mutates the global domain.
            let stack: Vec<HighsDomainChange> =
                self.mipdata_.domains[gd].get_domain_change_stack().to_vec();
            for domchg in stack {
                let col = domchg.column as usize;
                let global = self.get_domain_2();
                if (domchg.boundtype == HighsBoundType::kLower
                    && domchg.boundval > global.col_lower_[col])
                    || (domchg.boundtype == HighsBoundType::kUpper
                        && domchg.boundval < global.col_upper_[col])
                {
                    detach!(self.mipdata_.domains[0], |d| {
                        d.change_bound(self, domchg, HighsDomain_Reason::unspecified());
                    });
                }
            }
        }
    }

    /// C++ `doResetWorkerDomain` (lambda in run()), run as a task on worker `i`.
    fn run_do_reset_worker_domain(&mut self, i: usize) {
        let gd = self.mipdata_.workers[i].get_global_domain();
        let stack: Vec<HighsDomainChange> = self.get_domain_2().get_domain_change_stack().to_vec();
        // HiGHS 1.16 #3406 (h16_mip_workers): HighsMipSolver.cpp doResetWorkerDomain lines 421-427 (1.16.0)
        let h16_mip_workers = self.options_mip_.h16_mip_workers;
        for domchg in stack {
            detach!(self.mipdata_.domains[gd], |d| {
                d.change_bound(self, domchg, HighsDomain_Reason::unspecified());
            });
            if h16_mip_workers && self.mipdata_.domains[gd].infeasible() {
                return;
            }
        }
        detach!(self.mipdata_.domains[gd], |d| {
            d.set_domain_change_stack(self, &[]);
        });
        // Warning: Resetting local domain cannot be done in parallel (changes
        // propagationDomains of main pool)
        HighsSearch::with_worker_search(self, i, |s, ctx| s.reset_local_domain(ctx));
        self.mipdata_.domains[gd].clear_changed_cols();
    }

    /// C++ `resetGlobalDomain` (lambda in run()).
    fn run_reset_global_domain(&mut self, force: bool, reset_workers: bool, num_workers: HighsInt) {
        // if global propagation found bound changes, we update the domain
        let changed_nonempty = !self.get_domain_2().get_changed_cols().is_empty();
        if changed_nonempty || force {
            let num_changed = self.get_domain_2().get_changed_cols().len();
            crate::highs_log_dev!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "added {} global bound changes\n",
                num_changed as HighsInt
            );
            detach!(self.mipdata_.cliquetable, |ct| {
                detach!(self.mipdata_.domains[0], |global| {
                    ct.cleanup_fixed(self, global);
                });
            });
            if self.has_multiple_workers() && reset_workers {
                // Sync worker domains here. cleanupFixed might have found extra changes
                let indices: Vec<HighsInt> = (0..num_workers).collect();
                self.run_task(
                    &|ms: &mut HighsMipSolver, i: HighsInt| {
                        ms.run_do_reset_worker_domain(i as usize)
                    },
                    false,
                    true,
                    &indices,
                );
            }
            let changed: Vec<HighsInt> = self.get_domain_2().get_changed_cols().to_vec();
            for col in changed {
                detach!(self.mipdata_.implications, |imp| {
                    imp.cleanup_varbounds(self, col);
                });
            }

            detach!(self.mipdata_.domains[0], |d| {
                d.set_domain_change_stack(self, &[]);
            });
            if !self.has_multiple_workers() {
                HighsSearch::with_worker_search(self, 0, |s, ctx| s.reset_local_domain(ctx));
            }
            self.get_domain().clear_changed_cols();
            self.remove_fixed_indices();
        }
    }

    /// C++ `syncGlobalPseudoCost` (lambda in run()).
    fn run_sync_global_pseudo_cost(&mut self) {
        if !self.has_multiple_workers() {
            return;
        }
        for w in 0..self.mipdata_.workers.len() {
            let pc = self.mipdata_.workers[w].get_pseudocost();
            detach!(self.mipdata_.pseudocosts[0], |global| {
                detach!(self.mipdata_.pseudocosts[pc], |worker_pc| {
                    global.flush_pseudo_cost(worker_pc);
                });
            });
        }
    }

    /// C++ `resetWorkerPseudoCosts` (lambda in run()).
    fn run_reset_worker_pseudo_costs(&mut self, indices: &[usize]) {
        if !self.has_multiple_workers() {
            return;
        }
        let task_indices: Vec<HighsInt> = indices.iter().map(|&w| w as HighsInt).collect();
        self.run_task(
            &|ms: &mut HighsMipSolver, i: HighsInt| {
                let pc = ms.mipdata_.workers[i as usize].get_pseudocost();
                detach!(ms.mipdata_.pseudocosts[0], |global| {
                    detach!(ms.mipdata_.pseudocosts[pc], |worker_pc| {
                        global.sync_pseudo_cost(worker_pc);
                    });
                });
            },
            false,
            false,
            &task_indices,
        );
    }

    /// C++ `getSearchIndicesWithNoNodes` (lambda in run()).
    fn run_get_search_indices_with_no_nodes(
        &mut self,
        indices: &mut Vec<usize>,
        restarts: &mut [RestartVote],
        num_workers: HighsInt,
    ) {
        indices.clear();
        for r in restarts.iter_mut() {
            *r = RestartVote::NoCheck;
        }
        let num_search_workers = std::cmp::min(
            num_workers,
            self.mipdata_.nodequeue.num_active_nodes() as HighsInt,
        );
        // HiGHS 1.16 #3192 (h16_mip_workers): HighsMipSolver.cpp getSearchIndicesWithNoNodes
        // lines 643-656 (1.16.0): every search worker may run heuristics.
        if self.options_mip_.h16_mip_workers {
            for i in 0..num_search_workers {
                indices.push(i as usize);
                self.mipdata_.workers[i as usize].set_allow_heuristics(true);
            }
            return;
        }
        let num_heuristic_workers = if self.mipdata_.upper_bound < kHighsInf {
            std::cmp::max(1, (num_search_workers + 3) / 4)
        } else {
            std::cmp::max(1, (num_search_workers + 1) / 2)
        };

        for i in 0..num_search_workers {
            indices.push(i as usize);
            self.mipdata_.workers[i as usize].set_allow_heuristics(i < num_heuristic_workers);
        }
    }

    /// C++ `installNodes` (lambda in run()).
    fn run_install_nodes(
        &mut self,
        indices: &[usize],
        limit_reached: &mut bool,
        st: &mut RunState,
    ) {
        // rusolve mip_plunge_bestbound_freq: the best-bound pick period (HiGHS: every 10th leaf)
        let best_bound_freq = self.run_best_bound_freq();
        for &w in indices {
            if (indices.len() == 1 && st.num_queue_leaves - st.last_lb_leave >= best_bound_freq)
                || (indices.len() > 1 && w == 0)
            {
                let node = self.mipdata_.nodequeue.pop_best_bound_node();
                if crate::mip::search::tree_trace() {
                    eprintln!(
                        "TQ bb {} {} {}",
                        node.lower_bound, node.estimate, node.depth
                    );
                }
                HighsSearch::with_worker_search(self, w, |s, ctx| s.install_node(ctx, node));
                st.last_lb_leave = st.num_queue_leaves;
            } else {
                let best_bound_node_stack_size =
                    self.mipdata_.nodequeue.get_best_bound_domchg_stack_size();
                let best_bound_node_lb = self.mipdata_.nodequeue.get_best_lower_bound();
                let next_node = self.pop_best_node_hooked();
                if next_node.lower_bound == best_bound_node_lb
                    && next_node.domchgstack.len() as HighsInt == best_bound_node_stack_size
                {
                    st.last_lb_leave = st.num_queue_leaves;
                }
                if crate::mip::search::tree_trace() {
                    eprintln!(
                        "TQ est {} {} {} {}",
                        next_node.lower_bound,
                        next_node.estimate,
                        next_node.depth,
                        best_bound_node_lb
                    );
                }
                HighsSearch::with_worker_search(self, w, |s, ctx| s.install_node(ctx, next_node));
            }

            st.num_queue_leaves += 1;

            if self.run_worker_search(w).get_current_estimate() >= self.mipdata_.upper_limit {
                st.num_stall_nodes += 1;
                if self.options_mip_.mip_max_stall_nodes != kHighsIInf
                    && st.num_stall_nodes >= self.options_mip_.mip_max_stall_nodes as i64
                {
                    *limit_reached = true;
                    self.modelstatus_ = HighsModelStatus::kSolutionLimit;
                    break;
                }
            } else {
                st.num_stall_nodes = 0;
            }
        }
    }

    /// C++ `evaluateNode` (lambda in run()).
    fn run_evaluate_node(&mut self, w: usize) -> bool {
        let queue = self.run_queue_ref();
        let result = HighsSearch::with_worker_search(self, w, |s, ctx| s.evaluate_node(ctx));
        if result == HighsSearch_NodeResult::kSubOptimal {
            if self.options_mip_.h16_mip_workers {
                // HiGHS 1.16 #3192 (h16_mip_workers): HighsMipSolver.cpp evaluateNode lines 690-704 (1.16.0)
                HighsSearch::with_worker_search(self, w, |s, ctx| s.stash_current_node(ctx));
            } else {
                HighsSearch::with_worker_search(self, w, |s, ctx| {
                    s.current_node_to_queue(ctx, queue)
                });
            }
            return true;
        }
        false
    }

    /// C++ `pruneNode` (lambda in run()).
    fn run_prune_node(&mut self, w: usize) -> bool {
        let mut pruned = false;
        if self.run_worker_search(w).current_node_pruned() {
            HighsSearch::with_worker_search(self, w, |s, ctx| s.backtrack(ctx, true));
            let gd = self.mipdata_.workers[w].get_global_domain();
            detach!(self.mipdata_.domains[gd], |d| {
                d.propagate(self);
            });
            pruned = true;
            let search = self.run_worker_search_mut(w);
            *search.get_local_nodes() += 1;
            *search.get_local_leaves() += 1;
        }
        let gd = self.mipdata_.workers[w].get_global_domain();
        self.mipdata_.domains[gd].infeasible() || pruned
    }

    /// C++ `separateAndStoreBasis` (lambda in run()). Returns true when the node is done.
    fn run_separate_and_store_basis(&mut self, w: usize) -> bool {
        // rusolve mip_node_sepa_maxbounddist / mip_node_sepa_freq: `NodeSepa::Full` unless set
        let gate = if self.options_mip_.mip_allow_cut_separation_at_nodes {
            self.node_sepa_gate(w)
        } else {
            NodeSepa::Skip
        };
        if gate == NodeSepa::Full {
            self.separate_at_node(w);
        } else if gate == NodeSepa::PoolOnly {
            self.mipdata_.workers[w]
                .sepa_ptr_
                .as_mut()
                .unwrap()
                .pool_only = true;
            self.separate_at_node(w);
            self.mipdata_.workers[w]
                .sepa_ptr_
                .as_mut()
                .unwrap()
                .pool_only = false;
        } else {
            let cp = self.mipdata_.workers[w].get_cut_pool();
            self.mipdata_.cutpools[cp].perform_aging();
        }

        let gd = self.mipdata_.workers[w].get_global_domain();
        if self.mipdata_.domains[gd].infeasible() {
            self.run_worker_search_mut(w).cutoff_node();
            if self.options_mip_.h16_mip_workers {
                // HiGHS 1.16 #3192 #3124 (h16_mip_workers): HighsMipSolver.cpp separateAndStoreBasis
                // lines 746-752 (1.16.0); true here is the caller's `infeasible = true`.
                HighsSearch::with_worker_search(self, w, |s, ctx| s.stash_open_nodes(ctx));
                self.run_assign_early_termination(w);
            } else {
                let queue = self.run_queue_ref();
                HighsSearch::with_worker_search(self, w, |s, ctx| {
                    s.open_nodes_to_queue(ctx, queue)
                });
            }
            return true;
        }

        let lp = self.mipdata_.workers[w].get_lp_relaxation();
        let status = self.mipdata_.lps[lp].get_status();
        if status != HighsLpRelaxation_Status::kError && status != HighsLpRelaxation_Status::kNotSet
        {
            self.mipdata_.lps[lp].store_basis();
        }

        let basis = self.mipdata_.lps[lp].get_stored_basis();
        let consistent = match &basis {
            Some(b) => is_basis_consistent(self.mipdata_.lps[lp].get_lp(), b),
            None => false,
        };
        if !consistent {
            let mut b = self.mipdata_.firstrootbasis.clone();
            let num_rows = self.mipdata_.lps[lp].num_rows() as usize;
            b.row_status.resize(num_rows, HighsBasisStatus::kBasic);
            self.mipdata_.lps[lp].set_stored_basis(Some(Arc::new(b)));
        }

        false
    }

    /// rusolve: whether the tree node of worker `w` separates cuts (mip_node_sepa_maxbounddist,
    /// mip_node_sepa_freq, mip_node_sepa_pool_only). With both gates at their defaults (1.0 and
    /// 1) every node separates, as in HiGHS, and nothing below is read.
    ///
    /// The two rules are SCIP's per-separator gates in `SCIPsepaExecLP` (src/scip/sepa.c,
    /// Apache-2.0): a separator runs at depth `d` only if `freq > 0 && d % freq == 0`, and only if
    /// `bounddist <= maxbounddist`, where `priceAndCutLoop` (src/scip/solve.c) computes
    /// `bounddist = (node lower bound - global lower bound) / (cutoff bound - global lower
    /// bound)`. See also Achterberg, "Constraint Integer Programming" (PhD thesis, TU Berlin,
    /// 2007), on separation in the tree. Here the gate applies to the whole node round, the
    /// node bound is the node's LP bound after `evaluateNode`, the global bound is
    /// `mipdata_.lower_bound` and the cutoff is `mipdata_.upper_limit`. Without an incumbent the
    /// distance is 0, so every node passes the bound rule until a solution is known.
    fn node_sepa_gate(&self, w: usize) -> NodeSepa {
        let o = &self.options_mip_;
        if o.mip_node_sepa_maxbounddist >= 1.0 && o.mip_node_sepa_freq <= 1 {
            return NodeSepa::Full;
        }
        let search = self.run_worker_search(w);
        let mut pass = true;
        if o.mip_node_sepa_freq > 1 && search.get_current_depth() % o.mip_node_sepa_freq != 0 {
            pass = false;
        }
        if pass && o.mip_node_sepa_maxbounddist < 1.0 {
            let glb = self.mipdata_.lower_bound;
            let cutoff = self.mipdata_.upper_limit;
            let denom = cutoff - glb;
            if cutoff < kHighsInf && denom > self.mipdata_.feastol {
                let bounddist = (search.get_current_lower_bound() - glb) / denom;
                pass = bounddist <= o.mip_node_sepa_maxbounddist;
            }
        }
        if pass {
            NodeSepa::Full
        } else if o.mip_node_sepa_pool_only {
            NodeSepa::PoolOnly
        } else {
            NodeSepa::Skip
        }
    }

    /// rusolve `mip_plunge_max_nodes`: the node limit of one plunge (HiGHS: `plungeLimit`, 100),
    /// lowered to the option value while `plunge_levers_active` holds. 0 = HiGHS.
    fn run_plunge_limit(&self, plunge_limit: HighsInt) -> HighsInt {
        match self.options_mip_.mip_plunge_max_nodes {
            m if m > 0 && plunge_levers_active(self, false) => m,
            _ => plunge_limit,
        }
    }

    /// rusolve `mip_plunge_bestbound_freq`: a single worker takes the best-bound queue node once
    /// this many queue picks have passed without one (HiGHS: 10), lowered to the option value
    /// while `plunge_levers_active` holds. 0 = HiGHS.
    fn run_best_bound_freq(&self) -> i64 {
        match self.options_mip_.mip_plunge_bestbound_freq {
            f if f > 0 && plunge_levers_active(self, false) => f as i64,
            _ => 10,
        }
    }

    /// C++ `backtrackPlunge` (lambda in run()). Returns true when the plunge must stop.
    fn run_backtrack_plunge(&mut self, w: usize) -> bool {
        let queue = self.run_queue_ref();
        let backtrack_plunge =
            HighsSearch::with_worker_search(self, w, |s, ctx| s.backtrack_plunge(ctx, queue));
        if !backtrack_plunge {
            return true;
        }

        let cf = self.mipdata_.workers[w].get_conflict_pool();
        if self.mipdata_.conflictpools[cf].get_num_conflicts()
            > self.options_mip_.mip_pool_soft_limit
        {
            self.mipdata_.conflictpools[cf].perform_aging(false);
        }
        false
    }

    /// C++ `runHeuristics` (lambda in run()). Returns true when the node is done.
    fn run_heuristics(&mut self, w: usize) -> bool {
        match self.run_heuristics_body(w) {
            RunHeuristicsOutcome::SubOptimal => true,
            RunHeuristicsOutcome::Pruned => false,
            RunHeuristicsOutcome::Ran => {
                let gd = self.mipdata_.workers[w].get_global_domain();
                self.mipdata_.domains[gd].infeasible()
            }
        }
    }

    /// HiGHS 1.16 (#3192): HighsMipSolver.cpp::run lambda `runHeuristics` (lines 791-847)
    /// C++: `auto runHeuristics = [&](const HighsInt i, bool& infeasible) -> bool`
    ///
    /// Used only under `h16_mip_workers`: a global-domain infeasibility sets `infeasible`
    /// and returns false, where the 1.15.1 lambda returns true.
    fn run_heuristics_h16(&mut self, w: usize, infeasible: &mut bool) -> bool {
        match self.run_heuristics_body(w) {
            RunHeuristicsOutcome::SubOptimal => true,
            RunHeuristicsOutcome::Pruned => false,
            RunHeuristicsOutcome::Ran => {
                let gd = self.mipdata_.workers[w].get_global_domain();
                if self.mipdata_.domains[gd].infeasible() {
                    *infeasible = true;
                }
                false
            }
        }
    }

    /// The `runHeuristics` body shared by the 1.15.1 and 1.16 lambdas, up to their
    /// different final global-domain test.
    fn run_heuristics_body(&mut self, w: usize) -> RunHeuristicsOutcome {
        let evaluate_node_result =
            HighsSearch::with_worker_search(self, w, |s, ctx| s.evaluate_node(ctx));

        if evaluate_node_result == HighsSearch_NodeResult::kSubOptimal {
            return RunHeuristicsOutcome::SubOptimal;
        }

        if self.run_worker_search(w).current_node_pruned() {
            *self.run_worker_search_mut(w).get_local_leaves() += 1;
            return RunHeuristicsOutcome::Pruned;
        }

        let lp = self.mipdata_.workers[w].get_lp_relaxation();
        if self.mipdata_.incumbent.is_empty() {
            // PORT NOTE: the LP solution is cloned (the heuristics may re-solve the LP).
            let lpsol = self.mipdata_.lps[lp]
                .get_lp_solver()
                .get_solution()
                .col_value
                .clone();
            detach!(self.mipdata_.heuristics, |h| {
                h.randomized_rounding(self, w, &lpsol);
            });
        }
        if self.mipdata_.incumbent.is_empty() {
            if self.options_mip_.mip_heuristic_run_rens {
                let lpsol = self.mipdata_.lps[lp]
                    .get_lp_solver()
                    .get_solution()
                    .col_value
                    .clone();
                detach!(self.mipdata_.heuristics, |h| {
                    h.rens(self, w, &lpsol);
                });
            }
        } else if self.options_mip_.mip_heuristic_run_rins && !self.options_mip_.mip_alns {
            let lpsol = self.mipdata_.lps[lp]
                .get_lp_solver()
                .get_solution()
                .col_value
                .clone();
            detach!(self.mipdata_.heuristics, |h| {
                h.rins(self, w, &lpsol);
            });
        }

        // rusolve mip_guided_dive (SCIP heur_guideddiving.c, freq 10 offset 7): dive toward the incumbent
        if self.options_mip_.mip_guided_dive
            && !self.mipdata_.incumbent.is_empty()
            && self.mipdata_.num_nodes % 10 == 7
        {
            detach!(self.mipdata_.heuristics, |h| {
                h.guided_dive(self, w);
            });
        }

        RunHeuristicsOutcome::Ran
    }

    /// C++ `dive` (lambda in run()). Returns true when the plunge must stop.
    fn run_dive(&mut self, w: usize, node_lim: HighsInt) -> bool {
        let dive_node_lim: i64 = if node_lim == kHighsIInf {
            i64::MAX
        } else {
            self.run_worker_search(w).nnodes + node_lim as i64
        };
        if !self.run_worker_search(w).current_node_pruned() {
            let search_dive_result =
                HighsSearch::with_worker_search(self, w, |s, ctx| s.dive(ctx, dive_node_lim));
            if search_dive_result == HighsSearch_NodeResult::kSubOptimal {
                return true;
            }
            // rusolve mip_plunge_quot: the dive stopped at a child above the plunge bound, not at
            // a leaf; end the plunge so the open nodes go to the queue.
            if std::mem::take(&mut self.run_worker_search_mut(w).plunge_bound_stop) {
                return true;
            }
            *self.run_worker_search_mut(w).get_local_leaves() += 1;
        }
        node_lim != kHighsIInf && self.run_worker_search(w).nnodes >= dive_node_lim
    }

    /// C++ `checkRestart` (lambda in run()). Reads the restart state only.
    fn run_check_restart(
        &self,
        w: usize,
        num_worker_votes: HighsInt,
        st: &RunState,
    ) -> RestartVote {
        let treeweight = self.run_worker_search(w).treeweight;
        self.run_check_restart_tw(w, num_worker_votes, treeweight, st)
    }

    /// HiGHS 1.16 (#3192): HighsMipSolver.cpp::run lambda `checkRestart` (lines 521-586)
    /// C++: `auto checkRestart = [&](const HighsMipWorker& worker, HighsInt numWorkerVotes, HighsCDouble currTreeWeight) -> RestartVote`
    ///
    /// The 1.15.1 lambda is this body with `currTreeWeight = worker.search_ptr_->treeweight`
    /// ([`Self::run_check_restart`]).
    fn run_check_restart_tw(
        &self,
        w: usize,
        mut num_worker_votes: HighsInt,
        curr_tree_weight: HighsCDouble,
        st: &RunState,
    ) -> RestartVote {
        let worker = &self.mipdata_.workers[w];
        let search = worker
            .search_ptr_
            .as_deref()
            .expect("HighsMipSolver::run: worker search is running (re-entrant call)");
        let n_nodes: i64 = search.nnodes + self.mipdata_.num_nodes;
        if !self.submip && n_nodes >= st.next_check && self.options_mip_.mip_allow_restart {
            let n_tree_restarts: HighsInt =
                self.mipdata_.num_restarts - self.mipdata_.num_restarts_root;
            let tree_weight: HighsCDouble = curr_tree_weight + self.mipdata_.pruned_treeweight;
            let one_minus_tree_weight: f64 = f64::from(1.0 - tree_weight);
            let tree_weight_delta: f64 = f64::from(tree_weight - st.treeweight_last_check);
            let curr_node_estim: f64 = (st.num_nodes_last_check
                - self.mipdata_.num_nodes_before_run) as f64
                + ((n_nodes - st.num_nodes_last_check) as f64) * one_minus_tree_weight
                    / cpp_max(tree_weight_delta, self.mipdata_.epsilon);

            let mut do_restart = false;

            let mut active_integer_ratio = 1.0 - self.percentage_inactive_integers() / 100.0;
            active_integer_ratio *= active_integer_ratio;
            let mut huge_tree = false;

            if !do_restart {
                let mut gap_reduction = 1.0;
                if worker.upper_limit != kHighsInf {
                    let old_gap = st.upper_lim_last_check - st.lower_bound_last_check;
                    let new_gap = worker.upper_limit - self.mipdata_.lower_bound;
                    gap_reduction = old_gap / new_gap;
                }

                if gap_reduction < 1.0 + (0.05 / active_integer_ratio)
                    && curr_node_estim
                        >= active_integer_ratio
                            * 20.0
                            * ((n_nodes - self.mipdata_.num_nodes_before_run) as f64)
                {
                    huge_tree = true;
                } else {
                    num_worker_votes = -kHighsIInf;
                }

                let min_huge_tree_offset: i64 = (self.mipdata_.num_leaves + search.nleaves
                    - self.mipdata_.num_leaves_before_run)
                    / 1000;
                // PORT NOTE: std::pow(1.5, int) is mapped to powf (see the doubts list).
                let min_huge_tree_estim: i64 = HighsIntegers::nearest_integer(
                    active_integer_ratio
                        * ((10 + min_huge_tree_offset) as f64)
                        * 1.5f64.powf(n_tree_restarts as f64),
                );

                do_restart = (st.num_huge_tree_estim as i64 + num_worker_votes as i64)
                    >= min_huge_tree_estim;
            }

            if do_restart {
                return RestartVote::WouldRestart;
            }
            if huge_tree {
                return RestartVote::HugeTree;
            }
            return RestartVote::NoHugeTree;
        }
        RestartVote::NoCheck
    }

    /// C++ `performRestart` (lambda in run()).
    fn run_perform_restart(&mut self) {
        crate::highs_log_user!(
            &self.options_mip_.log_options,
            HighsLogType::kInfo,
            "\nRestarting search from the root node\n"
        );
        if crate::mip::search::tree_trace() {
            eprintln!("TR {:.3} {}", self.timer_.read(), self.mipdata_.num_nodes);
        }
        self.perform_restart();
    }

    /// C++ `checkWorkerRestartVotes` (lambda in run()). Returns true when the search restarts.
    fn run_check_worker_restart_votes(&mut self, votes: &[RestartVote], st: &mut RunState) -> bool {
        let mut num_voters: HighsInt = 0;
        let mut num_huge_tree_votes: HighsInt = 0;
        let mut num_restart_votes: HighsInt = 0;
        for &vote in votes {
            if vote != RestartVote::NoCheck {
                num_voters += 1;
                if vote == RestartVote::WouldRestart {
                    num_restart_votes += 1;
                } else if vote == RestartVote::HugeTree {
                    num_huge_tree_votes += 1;
                }
            }
        }
        // rusolve mip_cutoff_bigm: a cutoff that fell enough to tighten big-M rows
        if num_voters >= 1 && self.cutoff_bigm_restart_due() {
            self.run_perform_restart();
            return true;
        }
        let forced_restart = self.options_mip_.mip_allow_restart
            && num_voters >= 2
            && (num_restart_votes as f64 / num_voters as f64) >= 0.25;
        // Force a restart if enough individual workers vote for it
        if forced_restart {
            self.run_perform_restart();
            return true;
        }
        // Using joint information after workers are synced, query a restart
        let joint_votes = (num_restart_votes + num_huge_tree_votes + 2) / 2;
        let vote = self.run_check_restart(0, joint_votes, st);
        if vote == RestartVote::WouldRestart {
            self.run_perform_restart();
            return true;
        }

        // Update the appropriate values in case of no restart
        if vote == RestartVote::HugeTree {
            st.next_check = self.mipdata_.num_nodes + 100;
            st.num_huge_tree_estim += joint_votes;
        } else if vote == RestartVote::NoHugeTree {
            st.num_huge_tree_estim = 0;
            st.treeweight_last_check = f64::from(self.mipdata_.pruned_treeweight);
            st.num_nodes_last_check = self.mipdata_.num_nodes;
            st.upper_lim_last_check = self.mipdata_.upper_limit;
            st.lower_bound_last_check = self.mipdata_.lower_bound;
        }
        false
    }

    /// C++ `evaluateNode`/`pruneNode`/... chain of `processNodes`'s `processNode` lambda:
    /// returns `Some(vote)` when `restarts[i] = checkRestart(worker, 1)` applies.
    fn run_process_node(
        &mut self,
        w: usize,
        skip_separation: bool,
        node_lim: HighsInt,
        plunge_limit: HighsInt,
        avgiter: f64,
        st: &RunState,
    ) -> Option<RestartVote> {
        if !self.run_process_node_search(w, skip_separation, node_lim, plunge_limit, avgiter) {
            return None;
        }
        if node_lim == kHighsIInf {
            return Some(self.run_check_restart(w, 1, st));
        }
        None
    }

    /// The body of C++ `processNode` before its final restart check. Returns false on each
    /// early `return` of the lambda.
    fn run_process_node_search(
        &mut self,
        w: usize,
        skip_separation: bool,
        node_lim: HighsInt,
        plunge_limit: HighsInt,
        avgiter: f64,
    ) -> bool {
        let mut nodes_explored: i64 = 0;
        if !skip_separation {
            self.run_evaluate_node(w);
            if self.run_prune_node(w) {
                return false;
            }
            if !self.parallel_lock_active() {
                if self.check_limits(0) {
                    return false;
                }
                self.print_display_line(kSolutionSourceNone);
            } else if HighsSearch::with_worker_search(self, w, |s, ctx| s.check_local_limits(ctx)) {
                return false;
            }
            if self.run_separate_and_store_basis(w) {
                return false;
            }
        }
        let cf = self.mipdata_.workers[w].get_conflict_pool();
        self.mipdata_.conflictpools[cf].perform_aging(false);

        let mut iterlimit: HighsInt =
            (10.0 * cpp_max(avgiter, self.mipdata_.avgrootlpiters)) as HighsInt;
        iterlimit = std::cmp::max(
            std::cmp::max(10000, iterlimit),
            ((3 * self.mipdata_.firstrootlpiters) / 2) as HighsInt,
        );
        let lp = self.mipdata_.workers[w].get_lp_relaxation();
        self.mipdata_.lps[lp].set_iteration_limit(iterlimit);

        let mut consider_heuristics = true;
        loop {
            if consider_heuristics
                && (skip_separation || node_lim == kHighsIInf)
                && self.mipdata_.workers[w].get_allow_heuristics()
                && self.more_heuristics_allowed()
                && self.run_heuristics(w)
            {
                break;
            }
            consider_heuristics = false;
            let gd = self.mipdata_.workers[w].get_global_domain();
            if self.mipdata_.domains[gd].infeasible() {
                break;
            }
            if self.run_dive(w, node_lim) {
                break;
            }
            let local_nodes = *self.run_worker_search_mut(w).get_local_nodes();
            if HighsSearch::with_worker_search(self, w, |s, ctx| s.check_limits(ctx, local_nodes)) {
                break;
            }

            if *self.run_worker_search_mut(w).get_local_nodes() + nodes_explored
                >= self.run_plunge_limit(plunge_limit) as i64
            {
                break;
            }
            if !self.parallel_lock_active() {
                nodes_explored += *self.run_worker_search_mut(w).get_local_nodes();
            }
            if self.run_backtrack_plunge(w) {
                break;
            }
            if !self.parallel_lock_active() {
                self.run_with_search_only(w, |s, ms| s.flush_statistics(ms));
                self.print_display_line(kSolutionSourceNone);
            }
        }
        true
    }

    /// C++ `processNodes` (lambda in run()): runs `processNode` on each index as a task.
    fn run_process_nodes(
        &mut self,
        indices: &[usize],
        restarts: &mut [RestartVote],
        skip_separation: bool,
        node_lim: HighsInt,
        plunge_limit: HighsInt,
        avgiter: f64,
        st: &RunState,
    ) {
        let task_indices: Vec<HighsInt> = indices.iter().map(|&w| w as HighsInt).collect();
        let restarts = crate::parallel::shared_mut::SharedMut::new(restarts);
        let process_node = |ms: &mut HighsMipSolver, i: HighsInt| {
            let w = i as usize;
            if let Some(vote) =
                ms.run_process_node(w, skip_separation, node_lim, plunge_limit, avgiter, st)
            {
                // SAFETY: disjoint index: task `i` writes only `restarts[i]` (C++ 870),
                // and the slots are read after `taskWait` (C++ 873, rule (o)).
                unsafe { restarts.get_mut_unchecked()[w] = vote };
            }
        };
        self.run_task(&process_node, true, false, &task_indices);
    }

    /// C++ `syncSepaStats` (lambda in run()).
    fn run_sync_sepa_stats(&mut self, w: usize) {
        let neighbourhood_queries = *self.mipdata_.workers[w].get_num_neighbourhood_queries();
        *self.mipdata_.cliquetable.get_num_neighbourhood_queries() += neighbourhood_queries;
        let sepa_lp_iterations = *self.mipdata_.workers[w].get_sepa_lp_iterations();
        self.mipdata_.sepa_lp_iterations += sepa_lp_iterations;
        self.mipdata_.total_lp_iterations += sepa_lp_iterations;
        self.mipdata_.workers[w].reset_sepa_stats();
    }

    /// HiGHS 1.16 (#3124 #3192 #3336 #3406): HighsMipSolver.cpp::run main solve loop (lines 976-1197)
    /// C++: `while (!mipdata_->nodequeue.empty()) { ... }` after `restart:` set-up
    ///
    /// The 1.16 search loop, run under `h16_mip_workers`. Workers no longer own node
    /// queues: each gets a batch of prepared nodes, stashes the open nodes it produces in
    /// `processed_nodes`, and the main thread moves them to the global queue after the
    /// batch. Returns true when the search restarts (`goto restart`).
    fn run_search_h16(&mut self, max_num_workers: HighsInt) -> bool {
        // Maximum node processing budget for each worker in a batch during parallel
        const K_MAX_NODES_PER_WORKER: HighsInt = 100;
        let max_nodes_per_worker_lim: HighsInt = if max_num_workers > 1 {
            K_MAX_NODES_PER_WORKER
        } else {
            10000
        };
        let mut st = RunState::new(self);
        let mut num_workers: HighsInt = 1;
        let mut search_indices: Vec<usize> = vec![0];
        let mut worker_restart_votes: Vec<RestartVote> =
            vec![RestartVote::NoCheck; max_num_workers as usize];
        let mut stall_nodes_per_worker: Vec<HighsInt> = vec![0; max_num_workers as usize];
        let mut root_node = true; // Don't separate the root node again
        // ramp-up
        let mut node_lim: HighsInt = if max_num_workers > 1 {
            1
        } else {
            max_nodes_per_worker_lim
        };
        while !self.mipdata_.nodequeue.empty() {
            // Possibly query existence of an external solution
            if !self.submip {
                self.query_external_solution(
                    self.solution_objective_,
                    kExternalMipSolutionQueryOriginBeforeDive,
                );
            }
            self.race_import();

            // S6 (mip_heuristic_run_oneopt|crossover): SCIP heur_oneopt.c lines 472-889, heur_crossover.c lines 903-1044 at 32dd425
            let s6_on = self.options_mip_.mip_heuristic_run_oneopt
                || self.options_mip_.mip_heuristic_run_crossover;
            if s6_on {
                crate::mip::scip_heuristics::driver::s6_run_tree_heuristics(self);
            }
            // rusolve mip_alns (Hendel 2022): sibling of the SCIP tree-heuristics hook
            if self.options_mip_.mip_alns {
                crate::mip::scip_heuristics::alns::alns_step(self, false);
            }
            crate::mip::scip_heuristics::proximity::fischetti_step(self);

            // Update global pseudo-cost with worker information
            self.run_sync_global_pseudo_cost();

            // Get new candidate worker search indices
            self.run_get_search_indices_with_no_nodes(
                &mut search_indices,
                &mut worker_restart_votes,
                num_workers,
            );

            // Only update worker's pseudo-costs that have been assigned a node
            self.run_reset_worker_pseudo_costs(&search_indices);

            // Assign nodes to workers
            if root_node {
                let node = self.mipdata_.nodequeue.pop_best_bound_node();
                self.mipdata_.workers[0].prepared_nodes.push(node);
            } else {
                self.run_prepare_nodes(&search_indices, &mut st);
            }

            if node_lim != max_nodes_per_worker_lim {
                if num_workers >= max_num_workers {
                    node_lim = std::cmp::max(20, 2 * node_lim);
                }
                if node_lim >= max_nodes_per_worker_lim {
                    node_lim = max_nodes_per_worker_lim;
                }
            }

            // Process nodes (separation / heuristics / dives)
            let avgiter = self.get_lp().get_avg_solve_iters();
            self.run_process_nodes_h16(
                &search_indices,
                &mut worker_restart_votes,
                &mut stall_nodes_per_worker,
                root_node,
                node_lim,
                100,
                avgiter,
                max_nodes_per_worker_lim,
                &st,
            );

            root_node = false;

            // Check if one worker sent an early global termination signal
            let early_termination_lp_iterations = self
                .mipdata_
                .worker_lp_iterations_stop
                .0
                .load(std::sync::atomic::Ordering::Relaxed);
            let mut early_terminated_worker: Option<usize> = None;
            if early_termination_lp_iterations != i64::MAX {
                for &i in &search_indices {
                    if self.mipdata_.workers[i].early_termination
                        && self.run_worker_search(i).lpiterations == early_termination_lp_iterations
                    {
                        early_terminated_worker = Some(i);
                        break;
                    }
                }
                debug_assert!(early_terminated_worker.is_some());
                match early_terminated_worker {
                    None => {
                        self.mipdata_
                            .worker_lp_iterations_stop
                            .0
                            .store(i64::MAX, std::sync::atomic::Ordering::Relaxed);
                        for &i in &search_indices {
                            self.mipdata_.workers[i].early_termination = false;
                        }
                    }
                    Some(e) => {
                        // Clear all buffered solve information from workers that are not
                        // the earliest to terminate
                        for &i in &search_indices {
                            if i == e {
                                continue;
                            }
                            self.mipdata_.workers[i].processed_nodes.clear();
                            self.run_worker_search_mut(i).reset_statistics();
                            self.mipdata_.workers[i].reset_heur_stats();
                            self.mipdata_.workers[i].reset_sepa_stats();
                            self.mipdata_.workers[i].solutions_.clear();
                        }
                        search_indices.clear();
                        search_indices.push(e);
                    }
                }
            }

            // Sync statistics, check infeasibility, and flush nodes from worker queues
            let mut infeasible = false;
            for &w in &search_indices {
                let gd = self.mipdata_.workers[w].get_global_domain();
                if self.mipdata_.domains[gd].infeasible() {
                    infeasible = true;
                }
                debug_assert!(!self.run_worker_search(w).has_node());
                let processed = std::mem::take(&mut self.mipdata_.workers[w].processed_nodes);
                for (node, count_tree_weight) in processed {
                    // rusolve mip_symmetry_dynamic: carry the branching path on the queue transfer
                    // rusolve mip_node_install_basis: and the parent basis
                    let tmp_tree_weight =
                        if self.options_mip_.mip_symmetry_dynamic || node.basis.is_some() {
                            self.mipdata_.nodequeue.emplace_node_ext(
                                node.domchgstack,
                                node.branchings,
                                node.lower_bound,
                                node.estimate,
                                node.depth,
                                node.s3_path,
                                node.basis,
                            )
                        } else {
                            self.mipdata_.nodequeue.emplace_node(
                                node.domchgstack,
                                node.branchings,
                                node.lower_bound,
                                node.estimate,
                                node.depth,
                            )
                        };
                    if count_tree_weight {
                        self.run_worker_search_mut(w).treeweight += tmp_tree_weight;
                    }
                }
                self.run_with_search_only(w, |s, ms| s.flush_statistics(ms));
                self.run_sync_sepa_stats(w);
                detach!(self.mipdata_.heuristics, |h| h.flush_statistics(self, w));
            }

            if infeasible {
                self.mipdata_.nodequeue.clear();
                self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                let ub = self.mipdata_.upper_bound;
                self.update_lower_bound(cpp_min(kHighsInf, ub), true, true);
                break;
            }

            self.run_sync_solutions();

            if early_terminated_worker.is_none() {
                let best = self.mipdata_.nodequeue.get_best_lower_bound();
                let ub = self.mipdata_.upper_bound;
                self.update_lower_bound(cpp_min(ub, best), true, true);
            }

            let mut limit_reached = self.check_limits(0);
            if limit_reached || early_terminated_worker.is_some() {
                self.print_display_line(kSolutionSourceNone);
                break;
            }

            for &i in &search_indices {
                if stall_nodes_per_worker[i] == 0 {
                    st.num_stall_nodes = 0;
                } else {
                    st.num_stall_nodes += stall_nodes_per_worker[i] as i64;
                }
                stall_nodes_per_worker[i] = 0;
                if self.options_mip_.mip_max_stall_nodes != kHighsIInf
                    && st.num_stall_nodes
                        >= std::cmp::max(1, self.options_mip_.mip_max_stall_nodes) as i64
                {
                    limit_reached = true;
                    self.modelstatus_ = HighsModelStatus::kSolutionLimit;
                    break;
                }
            }
            if limit_reached {
                break;
            }

            // Sync global information
            self.run_sync_pools(&search_indices);
            self.run_sync_global_domain(&search_indices);
            detach!(self.mipdata_.domains[0], |d| {
                d.propagate(self);
            });

            let feastol = self.mipdata_.feastol;
            let pruned = detach!(self.mipdata_.nodequeue, |nq| {
                detach!(self.mipdata_.domains[0], |d| nq
                    .prune_infeasible_nodes(self, d, feastol))
            });
            self.mipdata_.pruned_treeweight += pruned;

            if self.get_domain_2().infeasible() {
                self.mipdata_.nodequeue.clear();
                self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                let ub = self.mipdata_.upper_bound;
                self.update_lower_bound(cpp_min(kHighsInf, ub), true, true);
                self.print_display_line(kSolutionSourceNone);
                break;
            }

            let best = self.mipdata_.nodequeue.get_best_lower_bound();
            let ub = self.mipdata_.upper_bound;
            self.update_lower_bound(cpp_min(ub, best), true, true);
            self.print_display_line(kSolutionSourceNone);

            if self.mipdata_.nodequeue.empty() {
                break;
            }

            // Reset global domain and sync worker's global domains.
            let spawn_more_workers = num_workers < max_num_workers
                && self.mipdata_.nodequeue.num_nodes() > num_workers as i64;
            let reset_workers = self.has_multiple_workers();
            self.run_reset_global_domain(spawn_more_workers, reset_workers, num_workers);

            for w in 0..self.mipdata_.workers.len() {
                let gd = self.mipdata_.workers[w].get_global_domain();
                if self.mipdata_.domains[gd].infeasible() {
                    infeasible = true;
                    break;
                }
            }
            if infeasible {
                self.mipdata_.nodequeue.clear();
                self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                let ub = self.mipdata_.upper_bound;
                self.update_lower_bound(cpp_min(kHighsInf, ub), true, true);
                self.print_display_line(kSolutionSourceNone);
                break;
            }

            if node_lim == max_nodes_per_worker_lim
                && self.run_check_worker_restart_votes(&worker_restart_votes, &mut st)
            {
                return true;
            }

            if spawn_more_workers {
                let new_max_num_workers = std::cmp::min(
                    self.mipdata_.nodequeue.num_nodes() as HighsInt,
                    max_num_workers,
                );
                self.get_pseudo_cost().remove_changed();
                if num_workers == 1 {
                    self.run_construct_additional_worker_data();
                }
                self.run_create_new_workers(new_max_num_workers - num_workers);
                num_workers = new_max_num_workers;
            }
        }
        false
    }

    /// HiGHS 1.16 (#3192): HighsMipSolver.cpp::run lambda `prepareNodes` (lines 658-688)
    /// C++: `auto prepareNodes = [&](const std::vector<HighsInt>& indices)`
    fn run_prepare_nodes(&mut self, indices: &[usize], st: &mut RunState) {
        let num_nodes_per_worker: HighsInt = if indices.len() > 1
            && self.mipdata_.nodequeue.num_active_nodes() as usize >= 2 * indices.len()
        {
            2
        } else {
            1
        };
        // rusolve mip_plunge_bestbound_freq: the best-bound pick period (HiGHS: every 10th leaf)
        let best_bound_freq = self.run_best_bound_freq();
        for j in 0..num_nodes_per_worker {
            for &i in indices {
                if (indices.len() == 1 && st.num_queue_leaves - st.last_lb_leave >= best_bound_freq)
                    || (indices.len() > 1 && i == 0 && j == 0)
                {
                    let node = self.mipdata_.nodequeue.pop_best_bound_node();
                    if crate::mip::search::tree_trace() {
                        eprintln!(
                            "TQ bb {} {} {}",
                            node.lower_bound, node.estimate, node.depth
                        );
                    }
                    self.mipdata_.workers[i].prepared_nodes.push(node);
                    st.last_lb_leave = st.num_queue_leaves;
                } else {
                    let best_bound_node_stack_size =
                        self.mipdata_.nodequeue.get_best_bound_domchg_stack_size();
                    let best_bound_node_lb = self.mipdata_.nodequeue.get_best_lower_bound();
                    let next_node = self.pop_best_node_hooked();
                    if next_node.lower_bound == best_bound_node_lb
                        && next_node.domchgstack.len() as HighsInt == best_bound_node_stack_size
                    {
                        st.last_lb_leave = st.num_queue_leaves;
                    }
                    if crate::mip::search::tree_trace() {
                        eprintln!(
                            "TQ est {} {} {} {}",
                            next_node.lower_bound,
                            next_node.estimate,
                            next_node.depth,
                            best_bound_node_lb
                        );
                    }
                    self.mipdata_.workers[i].prepared_nodes.push(next_node);
                }
                st.num_queue_leaves += 1;
            }
        }
    }

    /// HiGHS 1.16 (#3124): HighsMipSolver.cpp::run lambda `assignEarlyTermination` (lines 706-711)
    /// C++: `auto assignEarlyTermination = [&](HighsMipWorker& worker)`
    fn run_assign_early_termination(&mut self, w: usize) {
        let gd = self.mipdata_.workers[w].get_global_domain();
        if self.mipdata_.domains[gd].infeasible()
            || self.mipdata_.lower_bound > self.mipdata_.workers[w].get_optimality_limit(self)
        {
            self.update_worker_early_termination(w);
        }
    }

    /// HiGHS 1.16 (#3124 #3192): HighsMipSolver.cpp::run lambda `pruneNode` (lines 713-732)
    /// C++: `auto pruneNode = [&](const HighsInt i, bool& infeasible) -> bool`
    fn run_prune_node_h16(&mut self, w: usize, infeasible: &mut bool) -> bool {
        debug_assert!(!*infeasible);
        let mut pruned = false;
        if self.run_worker_search(w).current_node_pruned() {
            HighsSearch::with_worker_search(self, w, |s, ctx| s.backtrack(ctx, true));
            let gd = self.mipdata_.workers[w].get_global_domain();
            detach!(self.mipdata_.domains[gd], |d| {
                d.propagate(self);
            });
            pruned = true;
            let search = self.run_worker_search_mut(w);
            *search.get_local_nodes() += 1;
            *search.get_local_leaves() += 1;
            if self.mipdata_.domains[gd].infeasible() {
                self.run_assign_early_termination(w);
                *infeasible = true;
            }
        }
        pruned
    }

    /// HiGHS 1.16 (#3192 #3336): HighsMipSolver.cpp::run lambda `processNodes` (lines 865-978)
    /// C++: `auto processNodes = [&](const std::vector<HighsInt>& indices, std::vector<RestartVote>& restarts, std::vector<HighsInt>& stallNodes, const bool skip_separation, const HighsInt nodeLim, const HighsInt plungeLimit, double avgiter)`
    #[allow(clippy::too_many_arguments)]
    fn run_process_nodes_h16(
        &mut self,
        indices: &[usize],
        restarts: &mut [RestartVote],
        stall_nodes: &mut [HighsInt],
        skip_separation: bool,
        node_lim: HighsInt,
        plunge_limit: HighsInt,
        avgiter: f64,
        max_nodes_per_worker_lim: HighsInt,
        st: &RunState,
    ) {
        let task_indices: Vec<HighsInt> = indices.iter().map(|&w| w as HighsInt).collect();
        let restarts = crate::parallel::shared_mut::SharedMut::new(restarts);
        let stall_nodes = crate::parallel::shared_mut::SharedMut::new(stall_nodes);
        let process_worker_nodes = |ms: &mut HighsMipSolver, i: HighsInt| {
            let w = i as usize;
            let (stall, vote) = ms.run_process_worker_nodes_h16(
                w,
                skip_separation,
                node_lim,
                plunge_limit,
                avgiter,
                max_nodes_per_worker_lim,
                st,
            );
            // SAFETY: disjoint index: task `i` writes only `stallNodes[i]` and
            // `restarts[i]` (C++ 881, 973), read after `taskWait` (rule (o)).
            unsafe {
                stall_nodes.get_mut_unchecked()[w] += stall;
                if let Some(vote) = vote {
                    restarts.get_mut_unchecked()[w] = vote;
                }
            }
        };
        self.run_task(&process_worker_nodes, true, false, &task_indices);
    }

    /// HiGHS 1.16 (#3124 #3192 #3336): HighsMipSolver.cpp::run lambda `processWorkerNodes` (lines 872-976)
    /// C++: `auto processWorkerNodes = [&](const HighsInt i)`
    ///
    /// Returns the worker's stall-node increment (C++ `stallNodes[i]++`) and its restart
    /// vote when `nodeLim == maxNodesPerWorkerLim` (C++ `restarts[i] = ...`).
    #[allow(clippy::too_many_arguments)]
    fn run_process_worker_nodes_h16(
        &mut self,
        w: usize,
        skip_separation: bool,
        node_lim: HighsInt,
        plunge_limit: HighsInt,
        avgiter: f64,
        max_nodes_per_worker_lim: HighsInt,
        st: &RunState,
    ) -> (HighsInt, Option<RestartVote>) {
        let mut stall: HighsInt = 0;
        self.mipdata_.workers[w].prep_node_idx = 0;
        let mut total_nodes_explored: i64 = 0;
        let mut consider_heuristics = true;
        let mut stop = false;
        while self.mipdata_.workers[w].prep_node_idx
            != self.mipdata_.workers[w].prepared_nodes.len()
        {
            let mut nodes_explored: i64 = 0;
            let mut last_local_node_count: i64 = *self.run_worker_search_mut(w).get_local_nodes();
            let idx = self.mipdata_.workers[w].prep_node_idx;
            let node = std::mem::take(&mut self.mipdata_.workers[w].prepared_nodes[idx]);
            HighsSearch::with_worker_search(self, w, |s, ctx| s.install_node(ctx, node));
            self.mipdata_.workers[w].prep_node_idx += 1;
            if self.run_worker_search(w).get_current_estimate()
                >= self.mipdata_.workers[w].upper_limit
            {
                stall += 1;
            }
            if !skip_separation {
                self.run_evaluate_node(w);
                self.run_assign_early_termination(w);
                if self.run_prune_node_h16(w, &mut stop) {
                    total_nodes_explored += 1;
                    if self.mipdata_.workers[w].early_termination {
                        stop = true;
                    }
                    if stop || total_nodes_explored >= node_lim as i64 {
                        break;
                    }
                    continue;
                }
                if stop || HighsSearch::with_worker_search(self, w, |s, ctx| s.check_limits(ctx, 0))
                {
                    // C++ sets `stop = true` here; it is not read after the break.
                    break;
                }
                if !self.parallel_lock_active() {
                    self.print_display_line(kSolutionSourceNone);
                }
                if self.run_separate_and_store_basis(w) {
                    stop = true;
                }
                if stop {
                    break;
                }
            }
            let cf = self.mipdata_.workers[w].get_conflict_pool();
            self.mipdata_.conflictpools[cf].perform_aging(false);
            let mut iterlimit: HighsInt =
                (10.0 * cpp_max(avgiter, self.mipdata_.avgrootlpiters)) as HighsInt;
            iterlimit = std::cmp::max(
                std::cmp::max(10000, iterlimit),
                ((3 * self.mipdata_.firstrootlpiters) / 2) as HighsInt,
            );
            let lp = self.mipdata_.workers[w].get_lp_relaxation();
            self.mipdata_.lps[lp].set_iteration_limit(iterlimit);
            loop {
                if consider_heuristics
                    && (skip_separation || node_lim == max_nodes_per_worker_lim)
                    && self.mipdata_.workers[w].get_allow_heuristics()
                    && self.more_heuristics_allowed()
                {
                    if self.run_heuristics_h16(w, &mut stop) {
                        break;
                    }
                    if stop {
                        break;
                    }
                }
                consider_heuristics = false;
                let gd = self.mipdata_.workers[w].get_global_domain();
                if self.mipdata_.domains[gd].infeasible() {
                    stop = true;
                    break;
                }
                if self.run_dive(w, node_lim) {
                    break;
                }
                let local_nodes = *self.run_worker_search_mut(w).get_local_nodes();
                if HighsSearch::with_worker_search(self, w, |s, ctx| {
                    s.check_limits(ctx, local_nodes)
                }) {
                    stop = true;
                    break;
                }

                let local_nodes = *self.run_worker_search_mut(w).get_local_nodes();
                nodes_explored += local_nodes - last_local_node_count;
                last_local_node_count = local_nodes;
                if nodes_explored >= self.run_plunge_limit(plunge_limit) as i64 {
                    break;
                }
                if self.run_backtrack_plunge(w) {
                    break;
                }
                if !self.parallel_lock_active() {
                    self.run_with_search_only(w, |s, ms| s.flush_statistics(ms));
                    self.print_display_line(kSolutionSourceNone);
                    last_local_node_count = 0;
                }
            }
            let local_nodes = *self.run_worker_search_mut(w).get_local_nodes();
            nodes_explored += local_nodes - last_local_node_count;
            total_nodes_explored += nodes_explored;
            self.run_assign_early_termination(w);
            if self.mipdata_.workers[w].early_termination {
                stop = true;
            }
            HighsSearch::with_worker_search(self, w, |s, ctx| s.stash_open_nodes(ctx));
            if total_nodes_explored >= node_lim as i64 || stop {
                break;
            }
        }
        if self.run_worker_search(w).has_node() {
            HighsSearch::with_worker_search(self, w, |s, ctx| s.stash_open_nodes(ctx));
        }
        let count_tree_weight = self.run_worker_search(w).count_tree_weight;
        let idx = self.mipdata_.workers[w].prep_node_idx;
        let worker = &mut self.mipdata_.workers[w];
        let rest: Vec<HighsNodeQueue_OpenNode> = worker.prepared_nodes.drain(idx..).collect();
        for node in rest {
            worker.processed_nodes.push((node, count_tree_weight));
        }
        worker.prepared_nodes.clear();
        worker.prep_node_idx = 0;
        let mut vote = None;
        if node_lim == max_nodes_per_worker_lim {
            // Need to spoof tree weight
            let worker = &self.mipdata_.workers[w];
            let mut curr_tree_weight: HighsCDouble = self.run_worker_search(w).treeweight;
            for (node, counted) in &worker.processed_nodes {
                if *counted && node.lower_bound > worker.optimality_limit {
                    curr_tree_weight += crate::util::highs_cdouble::ldexp_f64(1.0, 1 - node.depth);
                }
            }
            vote = Some(self.run_check_restart_tw(w, 1, curr_tree_weight, st));
        }
        (stall, vote)
    }

    /// Glue (ours): C++ `worker.sepa_ptr_->separate(worker.search_ptr_->getLocalDomain())`
    /// in the `separateAndStoreBasis` lambda of `run()` (cpp 689). Detaches
    /// `workers[worker]`, takes its `sepa_ptr_` and `search_ptr_` out, detaches the
    /// worker's LP `lps[lp_]` (the one `sepa_ptr_` was built on), runs
    /// `separate` with the search's local domain as `propdomain`, and puts everything
    /// back.
    pub fn separate_at_node(&mut self, worker: usize) {
        let pending = detach!(self.mipdata_.workers[worker], |wk| {
            let mut sepa = wk.sepa_ptr_.take().expect("worker without sepa_ptr_");
            let mut search = wk.search_ptr_.take().expect("worker without search_ptr_");
            let lp_index = wk.lp_;
            detach!(self.mipdata_.lps[lp_index], |lp| {
                let mut ctx = SepaCtx {
                    ms: &mut *self,
                    worker: &mut *wk,
                    lp: &mut *lp,
                };
                sepa.separate(&mut ctx, search.get_local_domain());
            });
            // rusolve mip_ti_list_sched / mip_netdesign_heur: schedules built inside a separator.
            // Under the parallel lock they are buffered in the worker and merged at the sync
            // point (the same path as every other worker-side heuristic), never added to the
            // shared incumbent from a task.
            let mut pending = sepa.take_pending_solutions();
            wk.search_ptr_ = Some(search);
            wk.sepa_ptr_ = Some(sepa);
            if self.parallel_lock_active() {
                for x in pending.drain(..) {
                    wk.try_solution(self, &x, kSolutionSourceHeuristic);
                }
            }
            pending
        });
        // Serial path: `add_incumbent` visits every worker, so it has to run after this
        // worker's slot is back (inside the detach it would find the slot empty).
        for x in pending {
            self.try_solution(&x, kSolutionSourceHeuristic);
        }
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::runMipPresolve (lines 1213-1217)
    /// C++: `void runMipPresolve(const HighsInt presolve_reduction_limit)`
    pub fn run_mip_presolve(&mut self, presolve_reduction_limit: HighsInt) {
        self.new_data();
        self.init();
        self.run_mip_presolve_data(presolve_reduction_limit);
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::getPresolvedModel (lines 1219-1221)
    /// C++: `const HighsLp & getPresolvedModel() const`
    pub fn get_presolved_model(&self) -> &HighsLp {
        &self.mipdata_.presolved_model
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::getPresolveStatus (lines 1223-1225)
    /// C++: `HighsPresolveStatus getPresolveStatus() const`
    pub fn get_presolve_status(&self) -> HighsPresolveStatus {
        self.mipdata_.presolve_status
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::getPostsolveStack (lines 1227-1229)
    /// C++: `presolve::HighsPostsolveStack getPostsolveStack() const`
    pub fn get_postsolve_stack(&self) -> HighsPostsolveStack {
        self.mipdata_.post_solve_stack.clone()
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::initialiseTerminatorRecord (lines 1365-1369)
    /// C++: `std::vector<HighsModelStatus> initialiseTerminatorRecord(HighsInt num_instance) const`
    pub fn initialise_terminator_record(&self, num_instance: HighsInt) -> Vec<HighsModelStatus> {
        vec![HighsModelStatus::kNotset; num_instance as usize]
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::initialiseTerminator (lines 1371-1377)
    /// C++: `void initialiseTerminator(HighsInt num_instance_ = 0, HighsInt my_instance_ = kNoThreadInstance, HighsModelStatus * record_ = nullptr)`
    /// C++ default arguments: `num_instance_ = 0`, `my_instance_ = kNoThreadInstance`, `record_ = nullptr`
    pub fn initialise_terminator(
        &mut self,
        num_instance_: HighsInt,
        my_instance_: HighsInt,
        record_: Option<Vec<HighsModelStatus>>,
    ) {
        self.termination_status_ = HighsModelStatus::kNotset;
        self.terminator_.clear();
        self.terminator_
            .initialise(num_instance_, my_instance_, record_);
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::initialiseTerminator (lines 1379-1386)
    /// C++: `void initialiseTerminator(const HighsMipSolver & mip_solver)`
    ///
    /// PORT NOTE: the record is copied, not shared (see `HighsTerminator::initialise`).
    pub fn initialise_terminator_2(&mut self, mip_solver: &HighsMipSolver) {
        self.terminator_.clear();
        if !mip_solver.terminator_active() {
            return;
        }
        debug_assert!(mip_solver.terminator_concurrency() > 0);
        self.initialise_terminator(
            mip_solver.terminator_concurrency(),
            mip_solver.terminator_my_instance(),
            Some(mip_solver.terminator_.record.clone()),
        );
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::setParallelLock (lines 1388-1395)
    /// C++: `void setParallelLock(bool lock) const`
    ///
    /// `&mut self`: C++ `const`, writes `mipdata_` through the pointer.
    pub fn set_parallel_lock(&mut self, lock: bool) {
        if !self.has_multiple_workers() {
            return;
        }
        self.mipdata_.parallel_lock = lock;
        for conflictpool in self.mipdata_.conflictpools.iter_mut() {
            conflictpool.set_age_lock(lock);
        }
        let allow = !lock && !self.submip;
        self.mipdata_.cliquetable.set_allow_parallel(allow);
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::setProfiling (lines 1397-1400)
    /// C++: `void setProfiling(HighsProfiling * profiling)`
    ///
    /// Profiling is a stub (`profiling_` is an owned no-op `HighsProfiling`), so this
    /// takes no argument, like `HighsLpRelaxation::set_profiling`.
    pub fn set_profiling(&mut self) {}

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::numCol (lines 75-75)
    /// C++: `HighsInt numCol() const`
    pub fn num_col(&self) -> HighsInt {
        mip_model!(self).num_col_
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::numRow (lines 77-77)
    /// C++: `HighsInt numRow() const`
    pub fn num_row(&self) -> HighsInt {
        mip_model!(self).num_row_
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::numNonzero (lines 79-79)
    /// C++: `HighsInt numNonzero() const`
    pub fn num_nonzero(&self) -> HighsInt {
        mip_model!(self).a_matrix_.num_nz()
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::colCost (lines 81-81)
    /// C++: `const double * colCost() const`
    pub fn col_cost(&self) -> &[f64] {
        &mip_model!(self).col_cost_
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::colCost (lines 83-83)
    /// C++: `double colCost(HighsInt col) const`
    pub fn col_cost_2(&self, col: HighsInt) -> f64 {
        mip_model!(self).col_cost_[col as usize]
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::rowLower (lines 85-85)
    /// C++: `const double * rowLower() const`
    pub fn row_lower(&self) -> &[f64] {
        &mip_model!(self).row_lower_
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::rowLower (lines 87-87)
    /// C++: `double rowLower(HighsInt col) const`
    pub fn row_lower_2(&self, col: HighsInt) -> f64 {
        mip_model!(self).row_lower_[col as usize]
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::rowUpper (lines 89-89)
    /// C++: `const double * rowUpper() const`
    pub fn row_upper(&self) -> &[f64] {
        &mip_model!(self).row_upper_
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::rowUpper (lines 91-91)
    /// C++: `double rowUpper(HighsInt col) const`
    pub fn row_upper_2(&self, col: HighsInt) -> f64 {
        mip_model!(self).row_upper_[col as usize]
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::variableType (lines 93-95)
    /// C++: `const HighsVarType * variableType() const`
    pub fn variable_type(&self) -> &[HighsVarType] {
        &mip_model!(self).integrality_
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::variableType (lines 97-99)
    /// C++: `HighsVarType variableType(HighsInt col) const`
    // PORT NOTE (perf): `#[inline]`, as the C++ `variableType` header
    // inline it translates.
    #[inline]
    pub fn variable_type_2(&self, col: HighsInt) -> HighsVarType {
        mip_model!(self).integrality_[col as usize]
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::isColIntegral (lines 101-103)
    /// C++: `bool isColIntegral(const HighsInt col) const`
    // PORT NOTE (perf): `#[inline]`, as the C++ header inline.
    #[inline]
    pub fn is_col_integral(&self, col: HighsInt) -> bool {
        self.variable_type_2(col) != HighsVarType::kContinuous
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::isColInteger (lines 105-107)
    /// C++: `bool isColInteger(const HighsInt col) const`
    pub fn is_col_integer(&self, col: HighsInt) -> bool {
        self.variable_type_2(col) == HighsVarType::kInteger
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::isColContinuous (lines 109-111)
    /// C++: `bool isColContinuous(const HighsInt col) const`
    pub fn is_col_continuous(&self, col: HighsInt) -> bool {
        self.variable_type_2(col) == HighsVarType::kContinuous
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::setModel (lines 125-128)
    /// C++: `void setModel(const HighsLp & model)`
    ///
    /// `model_ = &model` -> the selector: `HighsMipModel::Orig` for `*orig_model_`,
    /// `HighsMipModel::Presolved` for `mipdata_->presolvedModel` (the only two models the
    /// C++ ever points it at).
    pub fn set_model(&mut self, model: HighsMipModel) {
        self.model_ = model;
        self.solution_objective_ = kHighsInf;
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::terminate (lines 153-155)
    /// C++: `bool terminate() const`
    pub fn terminate(&self) -> bool {
        self.termination_status_ != HighsModelStatus::kNotset
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::terminationStatus (lines 156-158)
    /// C++: `HighsModelStatus terminationStatus() const`
    pub fn termination_status(&self) -> HighsModelStatus {
        self.termination_status_
    }

    /// HiGHS: HighsMipSolver.h::HighsMipSolver::getMaxNumWorkers (lines 161-168)
    /// C++: `HighsInt getMaxNumWorkers() const`
    ///
    pub fn get_max_num_workers(&self) -> HighsInt {
        // HiGHS 1.16 #3192 #3406 (h16_mip_workers): HighsMipSolver.h lines 177-182 (1.16.0)
        if self.options_mip_.h16_mip_workers {
            if self.options_mip_.parallel == kHighsOffString || self.submip {
                return 1;
            }
            return crate::parallel::highs_parallel::num_threads();
        }
        if crate::parallel::highs_parallel::num_threads() == 1
            || self.options_mip_.parallel != kHighsOnString
            || self.submip
        {
            return 1;
        }
        (1.7 * crate::parallel::highs_parallel::num_threads() as f64).ceil() as HighsInt
    }
}

/// HiGHS: HStruct.h::HighsObjectiveSolution (lines 37-41)
#[derive(Clone, Debug, Default)]
pub struct HighsObjectiveSolution {
    /// C++: `double objective`
    pub objective: f64,
    /// C++: `std::vector<double> col_value`
    pub col_value: Vec<f64>,
}

impl HighsObjectiveSolution {
    /// HiGHS: HighsSolution.cpp::HighsObjectiveSolution::clear (line 2192)
    pub fn clear(&mut self) {
        self.col_value.clear();
    }
}

/// rusolve: what the tree-node separation gate (`HighsMipSolver::node_sepa_gate`) decides.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum NodeSepa {
    /// The full HiGHS separation round.
    Full,
    /// mip_node_sepa_pool_only: propagation and violated cut-pool cuts only.
    PoolOnly,
    /// No separation; the cut pool ages, as with mip_allow_cut_separation_at_nodes=false.
    Skip,
}

/// Outcome of the `runHeuristics` body before its final global-domain test.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RunHeuristicsOutcome {
    /// `evaluateNode` returned `kSubOptimal` (the lambda returns true).
    SubOptimal,
    /// The node was pruned (the lambda returns false).
    Pruned,
    /// The heuristics ran; the lambda's result depends on the global domain.
    Ran,
}

/// HiGHS: HighsMipSolver.cpp::HighsMipSolver::run, local `enum class RestartVote`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RestartVote {
    NoCheck,
    NoHugeTree,
    HugeTree,
    WouldRestart,
}

/// HiGHS: HighsMipSolver.cpp::HighsMipSolver::run, the local counters declared after
/// `restart:` (reset on every restart).
struct RunState {
    num_stall_nodes: i64,
    last_lb_leave: i64,
    num_queue_leaves: i64,
    num_huge_tree_estim: HighsInt,
    num_nodes_last_check: i64,
    next_check: i64,
    treeweight_last_check: f64,
    upper_lim_last_check: f64,
    lower_bound_last_check: f64,
}

impl RunState {
    /// C++ initialisation of the `run()` locals, after `destroyOldWorkers()`.
    fn new(ms: &HighsMipSolver) -> Self {
        RunState {
            num_stall_nodes: 0,
            last_lb_leave: 0,
            num_queue_leaves: 0,
            num_huge_tree_estim: 0,
            num_nodes_last_check: ms.mipdata_.num_nodes,
            next_check: ms.mipdata_.num_nodes,
            treeweight_last_check: 0.0,
            upper_lim_last_check: ms.mipdata_.upper_limit,
            lower_bound_last_check: ms.mipdata_.lower_bound,
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use crate::prelude::*;

    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;
    use std::sync::Mutex;
    const S6_GATE_BASE_OBJ: f64 = -6.0;
    const S6_GATE_BASE_NODES: i64 = 1;
    const S6_GATE_BASE_LP_ITERS: i64 = 1;

    fn s6_gate_knapsack_lp() -> HighsLp {
        use crate::lp_data::hconst::HighsVarType;
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 2;
        lp.sense_ = ObjSense::kMaximize;
        lp.offset_ = 0.0;
        lp.col_cost_ = vec![3.0, 2.0, 4.0, 1.0];
        lp.col_lower_ = vec![0.0; 4];
        lp.col_upper_ = vec![1.0; 4];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![4.0, 5.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 4, 6, 8];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![2.0, 1.0, 1.0, 2.0, 3.0, 1.0, 1.0, 2.0];
        lp.integrality_ = vec![HighsVarType::kInteger; 4];
        lp
    }

    fn s6_gate_solve(oneopt: bool, crossover: bool, submip: bool) -> (f64, i64, i64, i64) {
        let handle = std::thread::spawn(move || {
            crate::parallel::highs_parallel::initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.mip_heuristic_run_oneopt = oneopt;
            options.mip_heuristic_run_crossover = crossover;
            let lp = s6_gate_knapsack_lp();
            let mut callback = HighsCallback::default();
            let solution = HighsSolution::default();
            // C++ asserts `!submip || submip_level > 0` (HighsMipSolver.cpp:47).
            let submip_level = HighsInt::from(submip);
            let mut ms = HighsMipSolver::new(
                &mut callback,
                &options,
                &lp,
                &solution,
                submip,
                submip_level,
            );
            ms.run();
            let xover_ncalls = ms.mipdata_.s6.as_ref().map_or(0, |s| s.xover_ncalls);
            let result = (
                ms.mipdata_.upper_bound,
                ms.mipdata_.num_nodes,
                ms.mipdata_.total_lp_iterations,
                xover_ncalls,
            );
            drop(ms);
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            result
        });
        handle.join().unwrap()
    }

    #[test]
    fn s6_gate_off_matches_on_optimum() {
        let (off, _, _, _) = s6_gate_solve(false, false, false);
        let (o1, _, _, _) = s6_gate_solve(true, false, false);
        let (c1, _, _, _) = s6_gate_solve(false, true, false);
        let (both, _, _, _) = s6_gate_solve(true, true, false);
        assert!(off.is_finite(), "off optimum is not finite");
        assert!(o1.is_finite(), "oneopt optimum is not finite");
        assert!(c1.is_finite(), "crossover optimum is not finite");
        assert!(both.is_finite(), "both optimum is not finite");
        assert!((off - o1).abs() < 1e-9, "off {} != oneopt {}", off, o1);
        assert!((off - c1).abs() < 1e-9, "off {} != crossover {}", off, c1);
        assert!((off - both).abs() < 1e-9, "off {} != both {}", off, both);
    }

    /// A solution a separator built at a node (mip_ti_list_sched / mip_netdesign_heur) that
    /// improves the incumbent is added serially at threads = 1. It used to be added while the
    /// worker's slot was detached, which panicked ("Slot: value is detached"); the planted
    /// solution is the optimum of a finished solve whose bound is reset.
    #[test]
    fn serial_node_separator_solution_improves_incumbent() {
        std::thread::spawn(|| {
            initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.threads = 1;
            let lp = s6_gate_knapsack_lp();
            let mut callback = HighsCallback::default();
            let solution = HighsSolution::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let best = ms.mipdata_.upper_bound;
            let x = ms.mipdata_.incumbent.clone();
            assert!(best.is_finite() && !x.is_empty());
            ms.mipdata_.upper_bound = kHighsInf;
            // keep the primal-dual integral's bookkeeping consistent with the reset bound
            ms.mipdata_.primal_dual_integral.prev_ub = kHighsInf;
            ms.mipdata_.primal_dual_integral.prev_gap = kHighsInf;
            let planted = {
                let wk = &mut ms.mipdata_.workers[0];
                let sepa = wk.sepa_ptr_.as_mut().expect("worker sepa");
                let mut done = false;
                for b in &mut sepa.separators {
                    if let crate::mip::HighsSeparatorBox::Netdesign(n) = b {
                        n.pending.push((0, x.clone()));
                        done = true;
                    }
                }
                done
            };
            assert!(planted, "no netdesign separator in the worker");
            ms.separate_at_node(0);
            assert!(
                (ms.mipdata_.upper_bound - best).abs() < 1e-9,
                "incumbent {} != {}",
                ms.mipdata_.upper_bound,
                best
            );
            drop(ms);
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }

    /// A 3-row multi-dimensional knapsack with a real tree (the 4-column model solves at the
    /// root and never creates a second worker).
    fn tree_knapsack_lp(n: usize) -> HighsLp {
        use crate::lp_data::hconst::HighsVarType;
        let m = 3;
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as f64
        };
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = m as HighsInt;
        let mut weights = vec![vec![0.0; n]; m];
        for j in 0..n {
            let mut sum = 0.0;
            for row in weights.iter_mut() {
                row[j] = 10.0 + (next() % 40.0).floor();
                sum += row[j];
            }
            lp.col_cost_
                .push(-(sum / m as f64 + (next() % 10.0).floor()).floor());
        }
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf; m];
        lp.row_upper_ = weights
            .iter()
            .map(|r| (r.iter().sum::<f64>() / 3.0).floor())
            .collect();
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = m as HighsInt;
        lp.a_matrix_.start_ = (0..=n).map(|j| (j * m) as HighsInt).collect();
        for j in 0..n {
            for (i, row) in weights.iter().enumerate() {
                lp.a_matrix_.index_.push(i as HighsInt);
                lp.a_matrix_.value_.push(row[j]);
            }
        }
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp
    }

    /// The same planted separator solution under the parallel lock (2 workers): it is buffered
    /// in the worker, the shared incumbent is untouched, and `run_sync_solutions` adds it
    /// serially. Debug-friendly stand-in for the beasleyC3 solve at 8 threads.
    #[test]
    fn locked_node_separator_solution_is_buffered_then_synced() {
        std::thread::spawn(|| {
            initialize_scheduler(2);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.threads = 2;
            options.parallel = "on".to_string();
            options.h16_mip_workers = true;
            let lp = tree_knapsack_lp(48);
            let mut callback = HighsCallback::default();
            let solution = HighsSolution::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let best = ms.mipdata_.upper_bound;
            let x = ms.mipdata_.incumbent.clone();
            assert!(best.is_finite() && !x.is_empty());
            assert!(ms.has_multiple_workers(), "needs two workers");
            ms.mipdata_.upper_bound = kHighsInf;
            ms.mipdata_.primal_dual_integral.prev_ub = kHighsInf;
            ms.mipdata_.primal_dual_integral.prev_gap = kHighsInf;
            for w in 0..2 {
                ms.mipdata_.workers[w].upper_bound = kHighsInf;
            }
            {
                let sepa = ms.mipdata_.workers[0]
                    .sepa_ptr_
                    .as_mut()
                    .expect("worker sepa");
                let n = sepa
                    .separators
                    .iter_mut()
                    .find_map(|b| match b {
                        crate::mip::HighsSeparatorBox::Netdesign(n) => Some(n),
                        _ => None,
                    })
                    .expect("no netdesign separator in the worker");
                n.pending.push((0, x.clone()));
            }
            ms.set_parallel_lock(true);
            assert!(ms.parallel_lock_active());
            ms.separate_at_node(0);
            assert_eq!(ms.mipdata_.upper_bound, kHighsInf, "added under the lock");
            assert_eq!(ms.mipdata_.workers[0].solutions_.len(), 1, "not buffered");
            ms.set_parallel_lock(false);
            ms.run_sync_solutions();
            assert!((ms.mipdata_.upper_bound - best).abs() < 1e-9);
            drop(ms);
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }

    /// Solutions from the two separator kinds come back in production order, as the single
    /// `ti_pending_solutions` queue they replace gave them.
    #[test]
    fn pending_solutions_keep_production_order_across_separators() {
        std::thread::spawn(|| {
            initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            let lp = s6_gate_knapsack_lp();
            let mut callback = HighsCallback::default();
            let solution = HighsSolution::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let sepa = ms.mipdata_.workers[0]
                .sepa_ptr_
                .as_mut()
                .expect("worker sepa");
            for b in &mut sepa.separators {
                match b {
                    crate::mip::HighsSeparatorBox::TiPrecedence(t) => {
                        t.pending.push((20, vec![2.0]));
                        t.pending.push((30, vec![3.0]));
                    }
                    crate::mip::HighsSeparatorBox::Netdesign(n) => {
                        n.pending.push((10, vec![1.0]));
                        n.pending.push((40, vec![4.0]));
                    }
                    _ => {}
                }
            }
            let got: Vec<f64> = sepa.take_pending_solutions().iter().map(|v| v[0]).collect();
            assert_eq!(got, vec![1.0, 2.0, 3.0, 4.0]);
            drop(ms);
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }

    #[test]
    fn s6_gate_off_is_deterministic() {
        let (_, n1, lp1, _) = s6_gate_solve(false, false, false);
        let (_, n2, lp2, _) = s6_gate_solve(false, false, false);
        assert_eq!(n1, n2, "node count differs between off runs");
        assert_eq!(lp1, lp2, "lp iterations differ between off runs");
    }

    #[test]
    fn s6_gate_off_matches_baseline() {
        // Baseline hard-coded from an off run (both flags false); upper_bound is the minimised obj.
        let (obj, nodes, lp_iters, _) = s6_gate_solve(false, false, false);
        assert_eq!(
            (obj, nodes, lp_iters),
            (S6_GATE_BASE_OBJ, S6_GATE_BASE_NODES, S6_GATE_BASE_LP_ITERS),
            "off run differs from the hard-coded baseline"
        );
    }

    #[test]
    fn s6_gate_crossover_skipped_in_submip() {
        let (_, _, _, xover) = s6_gate_solve(false, true, true);
        assert_eq!(xover, 0, "crossover ran in submip");
    }

    /// Drives `s6_run_tree_heuristics` directly with a pool of 3 sols, so crossover passes its
    /// first throttle (it sets `xover_nextnodenumber` when the 3rd sol is new) and then stops at
    /// nwaitingnodes. In a sub-MIP the driver must not enter crossover at all.
    #[test]
    fn s6_gate_driver_skips_crossover_in_submip() {
        use crate::mip::scip_heuristics::S6Origin;
        fn probe(submip: bool) -> (i64, S6Origin) {
            let handle = std::thread::spawn(move || {
                crate::parallel::highs_parallel::initialize_scheduler(1);
                let mut options = HighsOptions::default();
                options.output_flag = false;
                options.mip_heuristic_run_crossover = true;
                let lp = s6_gate_knapsack_lp();
                let mut callback = HighsCallback::default();
                let solution = HighsSolution::default();
                let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
                ms.run();
                let nc = ms.num_col() as usize;
                {
                    let s6 = ms
                        .mipdata_
                        .s6
                        .as_mut()
                        .expect("s6 state exists when crossover is on");
                    s6.pool.sols.clear();
                    let mut next = s6.next_sol_index;
                    for k in 0..3 {
                        let mut x = vec![0.0; nc];
                        if k > 0 {
                            x[k - 1] = 1.0;
                        }
                        s6.pool
                            .insert(-(k as f64), x, S6Origin::Oneopt, 0, &mut next)
                            .expect("pool insert");
                    }
                    s6.next_sol_index = next;
                    s6.xover_prevlastsol = None;
                    s6.xover_nextnodenumber = i64::MIN;
                }
                ms.submip = submip;
                crate::mip::scip_heuristics::driver::s6_run_tree_heuristics(&mut ms);
                let s6 = ms.mipdata_.s6.as_ref().unwrap();
                let result = (s6.xover_nextnodenumber, s6.cur_origin);
                drop(ms);
                crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
                result
            });
            handle.join().unwrap()
        }
        let none = S6Origin::Highs(crate::mip::mip_solver_data::kSolutionSourceNone);
        let (sub_next, sub_origin) = probe(true);
        assert_eq!(sub_next, i64::MIN, "crossover was entered in a sub-MIP");
        assert_eq!(sub_origin, none);
        let (top_next, top_origin) = probe(false);
        assert_ne!(
            top_next,
            i64::MIN,
            "crossover was not entered outside a sub-MIP"
        );
        assert_eq!(top_origin, none, "cur_origin not reset after crossover");
    }

    #[test]
    fn s3_gate3_queue_transfer_keeps_path() {
        // the worker-to-main transfer and the root emplace carry the path only through emplace_node_s3
        use crate::mip::node_queue::HighsNodeQueue;
        use crate::mip::symmetry_dynamic::S3PathNode;
        use std::sync::Arc;
        let mut ms = HighsMipSolver::default();
        assert!(!ms.options_mip_.mip_symmetry_dynamic);
        let path = Arc::new(S3PathNode {
            parent: None,
            branch_col: 2,
            depth: 1,
            orbi: Vec::new(),
        });
        let mut worker = HighsNodeQueue::default();
        worker.emplace_node_s3(Vec::new(), vec![2], 1.0, 2.0, 1, Some(path.clone()));
        let node = worker.pop_best_node();
        let mut main = HighsNodeQueue::default();
        main.emplace_node_s3(
            node.domchgstack.clone(),
            node.branchings.clone(),
            node.lower_bound,
            node.estimate,
            node.depth,
            node.s3_path_handle(),
        );
        let moved = main.pop_best_node();
        assert!(Arc::ptr_eq(&moved.s3_path_handle().unwrap(), &path));
        main.emplace_node_s3(Vec::new(), Vec::new(), 1.0, 2.0, 1, None);
        assert!(main.pop_best_node().s3_path_handle().is_none());
    }

    /// Runs `f` on a fresh thread that owns a `threads`-worker executor (rule "Executor per
    /// thread").
    fn on_fresh_thread(threads: i32, f: impl FnOnce() + Send + 'static) {
        std::thread::spawn(move || {
            initialize_scheduler(threads);
            f();
            HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }

    /// Calls `run_task` and returns the indices `f` saw, in call order.
    fn run_and_log(
        ms: &mut HighsMipSolver,
        force_serial: bool,
        indices: &[HighsInt],
    ) -> Vec<HighsInt> {
        let log: Mutex<Vec<HighsInt>> = Mutex::new(Vec::new());
        ms.run_task(
            &|_, i| log.lock().unwrap().push(i),
            false,
            force_serial,
            indices,
        );
        log.into_inner().unwrap()
    }

    #[test]
    fn w5_mip_runtask_run_task_and_worker_count() {
        // getMaxNumWorkers (HighsMipSolver.h 161-168): ceil(1.7 * threads) only with
        // threads > 1, parallel = "on" and not a sub-MIP; otherwise 1.
        for (threads, expect_on) in [(1, 1), (2, 4), (4, 7), (8, 14)] {
            on_fresh_thread(threads, move || {
                let mut ms = HighsMipSolver::default();
                ms.options_mip_.parallel = kHighsOnString.to_string();
                ms.submip = false;
                assert_eq!(ms.get_max_num_workers(), expect_on);
                ms.submip = true;
                assert_eq!(ms.get_max_num_workers(), 1);
                ms.submip = false;
                ms.options_mip_.parallel = kHighsChooseString.to_string();
                assert_eq!(ms.get_max_num_workers(), 1);
                ms.options_mip_.parallel = kHighsOffString.to_string();
                assert_eq!(ms.get_max_num_workers(), 1);
            });
        }

        // runTask (HighsMipSolver.cpp 74-93)
        for threads in [1, 2, 4, 8] {
            for _repeat in 0..2 {
                on_fresh_thread(threads, move || {
                    let mut ms = HighsMipSolver::default();
                    // empty indices: early return, `f` never called
                    assert!(run_and_log(&mut ms, false, &[]).is_empty());
                    // force_serial: `f(i)` in index order on this thread
                    assert_eq!(run_and_log(&mut ms, true, &[0, 1, 2]), vec![0, 1, 2]);
                    // one index: never spawned
                    assert_eq!(run_and_log(&mut ms, false, &[5]), vec![5]);
                    // simulate concurrency: serial in index order
                    ms.options_mip_.mip_search_simulate_concurrency = true;
                    assert_eq!(run_and_log(&mut ms, false, &[0, 1, 2]), vec![0, 1, 2]);
                    ms.options_mip_.mip_search_simulate_concurrency = false;
                    // spawn path: every index once
                    let mut seen = run_and_log(&mut ms, false, &[0, 1, 2, 3]);
                    if threads == 1 {
                        // one worker: `taskWait` syncs the deque LIFO
                        assert_eq!(seen, vec![3, 2, 1, 0]);
                    }
                    seen.sort();
                    assert_eq!(seen, vec![0, 1, 2, 3]);
                });
            }
        }
    }

    /// a 12-item binary knapsack with a side row, min -value.x.
    fn s5_gate_lp() -> crate::lp_data::highs_lp::HighsLp {
        use crate::lp_data::highs_lp::HighsLp;
        let weights: [f64; 12] = [
            23.0, 31.0, 29.0, 44.0, 53.0, 38.0, 63.0, 85.0, 89.0, 82.0, 41.0, 57.0,
        ];
        let values: [f64; 12] = [
            92.0, 57.0, 49.0, 68.0, 60.0, 43.0, 67.0, 84.0, 87.0, 72.0, 45.0, 61.0,
        ];
        let n = weights.len();
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = 2;
        lp.col_cost_ = values.iter().map(|v| -v).collect();
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![265.0, 6.0];
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = 2;
        let mut start = Vec::new();
        let mut index = Vec::new();
        let mut value = Vec::new();
        for &w in weights.iter() {
            start.push(index.len() as HighsInt);
            index.push(0);
            value.push(w);
            index.push(1);
            value.push(1.0);
        }
        start.push(index.len() as HighsInt);
        lp.a_matrix_.start_ = start;
        lp.a_matrix_.index_ = index;
        lp.a_matrix_.value_ = value;
        lp
    }

    /// Brute-force optimum of [`s5_gate_lp`].
    fn s5_gate_brute_force() -> f64 {
        let lp = s5_gate_lp();
        let n = lp.num_col_ as usize;
        let mut best = 0.0f64;
        for mask in 0u32..(1 << n) {
            let (mut w, mut k, mut v) = (0.0f64, 0.0f64, 0.0f64);
            for j in 0..n {
                if mask & (1 << j) != 0 {
                    w += lp.a_matrix_.value_[2 * j];
                    k += 1.0;
                    v -= lp.col_cost_[j];
                }
            }
            if w <= 265.0 && k <= 6.0 && v > best {
                best = v;
            }
        }
        -best
    }

    /// One solve of [`s5_gate_lp`] with the given flags (None: default options).
    /// Returns (objective, num_nodes, total_lp_iterations, report heur
    /// (name, effort) list or None, heuristic_effort_used).
    #[allow(clippy::type_complexity)]
    fn s5_gate_run(
        flags: Option<(bool, bool)>,
    ) -> (f64, i64, i64, Option<Vec<(&'static str, usize)>>, usize) {
        use crate::lp_data::highs_callback::HighsCallback;
        use crate::lp_data::highs_options::HighsOptions;
        use crate::lp_data::hstruct::HighsSolution;
        on_fresh_thread_ret(move || {
            let lp = s5_gate_lp();
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.log_options.output_flag = false;
            options.presolve = kHighsOffString.to_string();
            options.mip_heuristic_run_feasibility_jump = false;
            if let Some((fpr, lm)) = flags {
                options.mip_heuristic_run_fpr = fpr;
                options.mip_heuristic_run_local_mip = lm;
            }
            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let report = ms
                .mipdata_
                .s5_report
                .as_ref()
                .map(|r| r.heur.iter().map(|e| (e.name, e.effort)).collect());
            (
                ms.solution_objective_,
                ms.mipdata_.num_nodes,
                ms.mipdata_.total_lp_iterations,
                report,
                ms.mipdata_.heuristic_effort_used,
            )
        })
    }

    /// Runs `f` on a fresh single-worker executor thread and returns its value.
    fn on_fresh_thread_ret<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
        std::thread::spawn(move || {
            initialize_scheduler(1);
            let out = f();
            HighsTaskExecutor::shutdown(true);
            out
        })
        .join()
        .unwrap()
    }

    #[test]
    fn s5_gate_each_flag_same_optimum() {
        let expected = s5_gate_brute_force();
        for flags in [
            None,
            Some((false, false)),
            Some((true, false)),
            Some((false, true)),
            Some((true, true)),
        ] {
            let (obj, _, _, _, _) = s5_gate_run(flags);
            assert!(
                (obj - expected).abs() < 1e-6,
                "flags {:?}: objective {} != brute force {}",
                flags,
                obj,
                expected
            );
        }
    }

    #[test]
    fn s5_gate_off_matches_baseline() {
        let (obj_d, nodes_d, iters_d, report_d, used_d) = s5_gate_run(None);
        let (obj_o, nodes_o, iters_o, report_o, used_o) = s5_gate_run(Some((false, false)));
        assert_eq!(obj_d.to_bits(), obj_o.to_bits());
        assert_eq!((nodes_d, iters_d), (nodes_o, iters_o));
        assert_eq!(
            (nodes_o, iters_o),
            (S5_GATE_BASELINE_NODES, S5_GATE_BASELINE_ITERS)
        );
        assert!(report_d.is_none() && report_o.is_none());
        assert_eq!((used_d, used_o), (0, 0));
    }

    /// Node count and LP iterations of an off run (s5-int at 4e929c0).
    const S5_GATE_BASELINE_NODES: i64 = 1;
    const S5_GATE_BASELINE_ITERS: i64 = 3;

    #[test]
    fn s5_gate_fpr_report() {
        let (_, _, _, report, used) = s5_gate_run(Some((true, false)));
        let report = report.expect("s5_report is created by the driver");
        assert!(!report.is_empty());
        assert_eq!(report[0].0, "fpr");
        assert!(report.iter().all(|e| e.0 == "fpr"));
        assert_eq!(used, report.iter().map(|e| e.1).sum::<usize>());
    }

    #[test]
    fn s5_gate_lm_report() {
        let (_, _, _, report, used) = s5_gate_run(Some((false, true)));
        let report = report.expect("s5_report is created by the driver");
        assert_eq!(report[0].0, "local_mip");
        assert_eq!(used, report.iter().map(|e| e.1).sum::<usize>());
        let (_, _, _, report, used) = s5_gate_run(Some((true, true)));
        let report = report.expect("s5_report is created by the driver");
        assert_eq!((report[0].0, report[1].0), ("fpr", "local_mip"));
        assert_eq!(used, report.iter().map(|e| e.1).sum::<usize>());
    }
}

/// HiGHS 1.16 #3135 / #3272 / #3523 (h16_mip_start) tests. The solve tests check the optimum of each
/// hand-built model, which both paths must reach; the `_log` tests check the user log text of
/// 1.16.0 checkAddSolution / solvingReport (on) against 1.15.1 runSetup / solvingReport (off).
#[cfg(test)]
mod h16_mip_start_tests {
    use crate::lp_data::highs_options::HighsOptions;
    use crate::lp_data::highs_run::Highs;
    use crate::prelude::*;

    /// min x + 2y, x integer in [0,3], y in [0,5], 2 <= x + y <= 4. Optimum 2 at (2, 0).
    fn two_col_model() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 2.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![3.0, 5.0];
        lp.row_lower_ = vec![2.0];
        lp.row_upper_ = vec![4.0];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kContinuous];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        lp
    }

    /// min x, x integer in [0,3], x >= 1: presolve reduces the model to empty.
    fn singleton_model() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![3.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![f64::INFINITY];
        lp.integrality_ = vec![HighsVarType::kInteger];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        lp
    }

    /// max 4a + 5b + 6c + 7d + 8e (as min of the negation), binaries a..e,
    /// 3a + 4b + 5c + 6d + 7e <= 12. Optimum -15 at (1, 1, 1, 0, 0); presolve does not reduce
    /// it, so the MIP-start code runs in `run_setup`.
    fn knapsack_model() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 5;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![-4.0, -5.0, -6.0, -7.0, -8.0];
        lp.col_lower_ = vec![0.0; 5];
        lp.col_upper_ = vec![1.0; 5];
        lp.row_lower_ = vec![f64::NEG_INFINITY];
        lp.row_upper_ = vec![12.0];
        lp.integrality_ = vec![HighsVarType::kInteger; 5];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 5;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2, 3, 4, 5];
        lp.a_matrix_.index_ = vec![0, 0, 0, 0, 0];
        lp.a_matrix_.value_ = vec![3.0, 4.0, 5.0, 6.0, 7.0];
        lp
    }

    fn start(col_value: Vec<f64>, row_value: Vec<f64>) -> HighsSolution {
        let mut s = HighsSolution::default();
        s.col_value = col_value;
        s.row_value = row_value;
        s.value_valid = true;
        s
    }

    /// Solve `lp` with `h16_mip_start` set to `h16`, optionally with a MIP start.
    fn solve(
        lp: HighsLp,
        h16: bool,
        mip_start: Option<HighsSolution>,
        output: bool,
    ) -> (HighsModelStatus, f64, Vec<f64>) {
        let mut options = HighsOptions::default();
        options.output_flag = output;
        options.h16_mip_start = h16;
        let mut highs = Highs::with_options(options);
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        if let Some(s) = mip_start {
            assert_eq!(highs.set_solution(&s), HighsStatus::kOk);
        }
        highs.run();
        (
            highs.get_model_status(),
            highs.get_info().objective_function_value,
            highs.get_solution().col_value.clone(),
        )
    }

    /// Solve `lp` with `h16_mip_start` set to `h16` and a MIP start, capturing the user log
    /// through the logging callback (console output off).
    fn solve_logged(
        lp: HighsLp,
        h16: bool,
        mip_start: HighsSolution,
    ) -> (HighsModelStatus, String) {
        use crate::lp_data::highs_callback::{
            HighsCallbackFunctionType, HighsCallbackInput, HighsCallbackOutput,
        };
        use std::sync::{Arc, Mutex};
        let mut options = HighsOptions::default();
        options.output_flag = true;
        options.log_to_console = false;
        options.h16_mip_start = h16;
        let mut highs = Highs::with_options(options);
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        assert_eq!(highs.set_solution(&mip_start), HighsStatus::kOk);
        let log = Arc::new(Mutex::new(String::new()));
        let sink = Arc::clone(&log);
        highs.set_callback(HighsCallbackFunctionType::new(
            move |kind: HighsInt,
                  message: &str,
                  _out: &HighsCallbackOutput,
                  _cb_in: &mut HighsCallbackInput| {
                if kind == kCallbackLogging {
                    sink.lock().unwrap().push_str(message);
                }
            },
        ));
        highs.start_callback(kCallbackLogging);
        highs.run();
        let text = log.lock().unwrap().clone();
        (highs.get_model_status(), text)
    }

    #[test]
    fn h16_mip_start_on_knapsack_feasible_start_log() {
        // 1.16 checkAddSolution: "MIP start solution is feasible: objective value is %.12g";
        // the incumbent objective equals the supplied one, so no "reduced space" line.
        // 1.16 solvingReport: zero violations are not reported.
        let (status, log) = solve_logged(knapsack_model(), true, start(vec![0.0; 5], vec![0.0]));
        assert_eq!(status, HighsModelStatus::kOptimal);
        assert!(
            log.contains("\nMIP start solution is feasible: objective value is 0\n"),
            "{log}"
        );
        assert!(
            !log.contains("MIP start solution is feasible, objective"),
            "{log}"
        );
        assert!(!log.contains("in the reduced space"), "{log}");
        assert!(
            log.contains("                    -15 (objective)\n"),
            "{log}"
        );
        assert!(!log.contains("(bound viol.)"), "{log}");
        assert!(!log.contains("violation)"), "{log}");
    }

    #[test]
    fn h16_mip_start_off_knapsack_feasible_start_log() {
        // 1.15.1 runSetup / solvingReport text
        let (status, log) = solve_logged(knapsack_model(), false, start(vec![0.0; 5], vec![0.0]));
        assert_eq!(status, HighsModelStatus::kOptimal);
        assert!(
            log.contains("\nMIP start solution is feasible, objective value is 0\n"),
            "{log}"
        );
        assert!(
            log.contains("                    0 (bound viol.)\n"),
            "{log}"
        );
        assert!(!log.contains("(variable bound violation)"), "{log}");
    }

    fn assert_optimum(result: (HighsModelStatus, f64, Vec<f64>), obj: f64) {
        assert_eq!(result.0, HighsModelStatus::kOptimal);
        assert!((result.1 - obj).abs() < 1e-9, "objective {}", result.1);
    }

    #[test]
    fn h16_mip_start_on_feasible_start_reaches_optimum() {
        let s = start(vec![2.0, 0.0], vec![2.0]);
        let r = solve(two_col_model(), true, Some(s), false);
        assert_optimum(r.clone(), 2.0);
        assert_eq!(r.2[0], 2.0);
    }

    #[test]
    fn h16_mip_start_on_infeasible_start_reaches_optimum() {
        // (0, 0) violates the row lower bound 2: the start is logged as infeasible
        let s = start(vec![0.0, 0.0], vec![0.0]);
        assert_optimum(solve(two_col_model(), true, Some(s), true), 2.0);
    }

    #[test]
    fn h16_mip_start_on_presolve_empty_model() {
        let s = start(vec![1.0], vec![1.0]);
        assert_optimum(solve(singleton_model(), true, Some(s), false), 1.0);
    }

    #[test]
    fn h16_mip_start_on_knapsack_infeasible_start_reaches_optimum() {
        // (2, 2) violates 2x + 2y <= 3: the start is logged as infeasible (on path)
        let s = start(vec![1.0; 5], vec![25.0]);
        assert_optimum(solve(knapsack_model(), true, Some(s), true), -15.0);
    }

    #[test]
    fn h16_mip_start_on_knapsack_feasible_start_reaches_optimum() {
        let s = start(vec![0.0; 5], vec![0.0]);
        assert_optimum(solve(knapsack_model(), true, Some(s), true), -15.0);
    }

    #[test]
    fn h16_mip_start_off_knapsack_infeasible_start_reaches_optimum() {
        let s = start(vec![1.0; 5], vec![25.0]);
        assert_optimum(solve(knapsack_model(), false, Some(s), true), -15.0);
    }

    #[test]
    fn h16_mip_start_off_feasible_start_reaches_optimum() {
        let s = start(vec![2.0, 0.0], vec![2.0]);
        let r = solve(two_col_model(), false, Some(s), false);
        assert_optimum(r.clone(), 2.0);
        assert_eq!(r.2[0], 2.0);
    }

    #[test]
    fn h16_mip_start_off_infeasible_start_reaches_optimum() {
        let s = start(vec![0.0, 0.0], vec![0.0]);
        assert_optimum(solve(two_col_model(), false, Some(s), true), 2.0);
    }

    #[test]
    fn h16_mip_start_off_presolve_empty_model() {
        let s = start(vec![1.0], vec![1.0]);
        assert_optimum(solve(singleton_model(), false, Some(s), false), 1.0);
    }
}

/// HiGHS 1.16-MIP-WORKERS: the 1.16.0 worker search loop (#3192 #3406 #3124) behind
/// `h16_mip_workers`. A 2-row binary knapsack with presolve off branches, so the loop in
/// `run_search_h16` runs; serial and 4-thread solves must reach the brute-force optimum,
/// and the off path must match the default options bit for bit.
#[cfg(test)]
mod h16_mip_workers_tests {
    use crate::lp_data::highs_callback::HighsCallback;
    use crate::lp_data::highs_lp::HighsLp;
    use crate::lp_data::highs_options::HighsOptions;
    use crate::lp_data::hstruct::HighsSolution;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;
    use crate::prelude::*;

    const N: usize = 40;
    const M: usize = 2;

    /// Strongly correlated 2-row knapsack from a fixed LCG: weights in [10, 59], value of
    /// item j = w0[j] + w1[j] + a term in [5, 14], each capacity 40 % of its row sum
    /// (floored). The root LP plus cuts does not close it, so the solve branches.
    #[allow(clippy::type_complexity)]
    fn data() -> ([[f64; N]; M], [f64; N], [f64; M]) {
        let mut s: u64 = 12345;
        let mut next = |lo: u64, span: u64| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (lo + (s >> 33) % span) as f64
        };
        let mut w = [[0.0; N]; M];
        let mut v = [0.0; N];
        for j in 0..N {
            for row in w.iter_mut() {
                row[j] = next(10, 50);
            }
            v[j] = w[0][j] + w[1][j] + next(5, 10);
        }
        let mut cap = [0.0; M];
        for (cap_i, row) in cap.iter_mut().zip(w.iter()) {
            *cap_i = (0.4 * row.iter().sum::<f64>()).floor();
        }
        (w, v, cap)
    }

    fn model() -> HighsLp {
        let (w, v, cap) = data();
        let mut lp = HighsLp::default();
        lp.num_col_ = N as HighsInt;
        lp.num_row_ = M as HighsInt;
        lp.col_cost_ = v.iter().map(|x| -x).collect();
        lp.col_lower_ = vec![0.0; N];
        lp.col_upper_ = vec![1.0; N];
        lp.row_lower_ = vec![-kHighsInf; M];
        lp.row_upper_ = cap.to_vec();
        lp.integrality_ = vec![HighsVarType::kInteger; N];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = N as HighsInt;
        lp.a_matrix_.num_row_ = M as HighsInt;
        let (mut start, mut index, mut value) = (Vec::new(), Vec::new(), Vec::new());
        for j in 0..N {
            start.push(index.len() as HighsInt);
            for (i, row) in w.iter().enumerate() {
                index.push(i as HighsInt);
                value.push(row[j]);
            }
        }
        start.push(index.len() as HighsInt);
        lp.a_matrix_.start_ = start;
        lp.a_matrix_.index_ = index;
        lp.a_matrix_.value_ = value;
        lp
    }

    /// Exact optimum by 0/1 dynamic programming over both (integer) capacities.
    fn brute_force() -> f64 {
        let (w, v, cap) = data();
        let (c0, c1) = (cap[0] as usize, cap[1] as usize);
        let mut best = vec![0.0f64; (c0 + 1) * (c1 + 1)];
        for j in 0..N {
            let (w0, w1) = (w[0][j] as usize, w[1][j] as usize);
            for a in (w0..=c0).rev() {
                for b in (w1..=c1).rev() {
                    let cand = best[(a - w0) * (c1 + 1) + (b - w1)] + v[j];
                    let cell = &mut best[a * (c1 + 1) + b];
                    if cand > *cell {
                        *cell = cand;
                    }
                }
            }
        }
        -best[c0 * (c1 + 1) + c1]
    }

    /// One solve with `h16` (None: default options) on a fresh executor of `threads`
    /// threads; `parallel=on` when `threads > 1`. Returns (objective, nodes, lp iterations,
    /// max_num_workers).
    fn run(h16: Option<bool>, threads: i32) -> (f64, i64, i64, HighsInt) {
        std::thread::spawn(move || {
            initialize_scheduler(threads);
            let lp = model();
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.log_options.output_flag = false;
            options.presolve = kHighsOffString.to_string();
            options.threads = threads as HighsInt;
            if threads > 1 {
                options.parallel = kHighsOnString.to_string();
            }
            if let Some(on) = h16 {
                options.h16_mip_workers = on;
            }
            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            let workers = ms.get_max_num_workers();
            ms.run();
            let out = (
                ms.solution_objective_,
                ms.mipdata_.num_nodes,
                ms.mipdata_.total_lp_iterations,
                workers,
            );
            HighsTaskExecutor::shutdown(true);
            out
        })
        .join()
        .unwrap()
    }

    #[test]
    fn h16_mip_workers_on_serial_reaches_optimum() {
        let (obj, nodes, _, workers) = run(Some(true), 1);
        assert_eq!(workers, 1);
        assert!(nodes > 1, "model must branch to reach the search loop");
        assert!(
            (obj - brute_force()).abs() < 1e-6,
            "{obj} != {}",
            brute_force()
        );
    }

    #[test]
    fn h16_mip_workers_on_parallel_reaches_optimum() {
        let expected = brute_force();
        for _ in 0..3 {
            let (obj, _, _, workers) = run(Some(true), 4);
            assert_eq!(workers, 4, "1.16: one worker per thread");
            assert!((obj - expected).abs() < 1e-6, "{obj} != {expected}");
        }
    }

    #[test]
    fn h16_mip_workers_off_matches_default() {
        let (obj_d, nodes_d, iters_d, workers_d) = run(None, 1);
        let (obj_o, nodes_o, iters_o, workers_o) = run(Some(false), 1);
        assert_eq!(obj_d.to_bits(), obj_o.to_bits());
        assert_eq!((nodes_d, iters_d, workers_d), (nodes_o, iters_o, workers_o));
        assert!((obj_o - brute_force()).abs() < 1e-6);
    }

    #[test]
    fn h16_mip_workers_off_parallel_worker_count() {
        let (obj, _, _, workers) = run(Some(false), 4);
        assert_eq!(workers, 7, "1.15.1: ceil(1.7 * threads) workers");
        assert!((obj - brute_force()).abs() < 1e-6);
    }

    /// rusolve racing: four racers with different seeds and settings share incumbents and
    /// return the brute-force optimum, and `mip_racing_threads = 0` is the plain solve.
    #[test]
    fn racing_reaches_optimum() {
        let expected = brute_force();
        let out = std::thread::spawn(move || {
            let lp = model();
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.log_options.output_flag = false;
            options.presolve = kHighsOffString.to_string();
            let mut callback = HighsCallback::default();
            let solver = crate::lp_data::highs_run::race_solve_mip(
                &mut callback,
                &options,
                &lp,
                &HighsSolution::default(),
                4,
            );
            (solver.solution_objective_, solver.modelstatus_)
        })
        .join()
        .unwrap();
        assert_eq!(out.1, HighsModelStatus::kOptimal);
        assert!((out.0 - expected).abs() < 1e-6, "{} != {expected}", out.0);
    }

    /// rusolve mip_root_lp_concurrent: the root LP race hands an optimal basis to the node
    /// loop and the solve reaches the brute-force optimum.
    #[test]
    fn root_lp_concurrent_reaches_optimum() {
        let expected = brute_force();
        let obj = std::thread::spawn(move || {
            let lp = model();
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.log_options.output_flag = false;
            options.presolve = kHighsOffString.to_string();
            options.threads = 2;
            options.mip_root_lp_concurrent = true;
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(
                &mut callback,
                &options,
                &lp,
                &HighsSolution::default(),
                false,
                0,
            );
            ms.run();
            ms.solution_objective_
        })
        .join()
        .unwrap();
        assert!((obj - expected).abs() < 1e-6, "{obj} != {expected}");
    }
}
