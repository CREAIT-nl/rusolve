//! HiGHS source: mip/HighsMipSolverData.cpp (lines 1842-3017)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;
use crate::util::c_printf::c_fmt_g;

impl HighsMipSolver {
    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::rootSeparationRound (lines 1842-1868)
    /// C++: `bool rootSeparationRound(HighsMipWorker & worker, HighsSeparation & sepa, HighsInt & ncuts, HighsLpRelaxation::Status & status)`
    ///
    /// `worker` is the index into `mipdata_.workers` (attached). The C++
    /// `sepa.separationRound(getDomain(), status)` is
    /// [`HighsMipSolver::separation_round_global`].
    pub fn root_separation_round(
        &mut self,
        worker: usize,
        sepa: &mut HighsSeparation,
        ncuts: &mut HighsInt,
        status: &mut HighsLpRelaxation_Status,
    ) -> bool {
        let mut tmp_lp_iters: i64 = -self.get_lp_2().get_num_lp_iterations();
        let lp_src_guard = crate::work::lp_kind(crate::work::LP_ROOT_SEPA);
        *ncuts = self.separation_round_global(worker, sepa, status);
        drop(lp_src_guard);
        tmp_lp_iters += self.get_lp_2().get_num_lp_iterations();
        self.mipdata_.avgrootlpiters = self.get_lp().get_avg_solve_iters();
        self.mipdata_.total_lp_iterations += tmp_lp_iters;
        self.mipdata_.sepa_lp_iterations += tmp_lp_iters;

        *status = self.evaluate_root_lp(worker);
        if *status == HighsLpRelaxation_Status::kInfeasible {
            return true;
        }

        // C++ holds `solvals` by reference while the heuristics run on `this`; a copy is
        // read-only equivalent.
        let solvals = self
            .get_lp_2()
            .get_lp_solver_2()
            .get_solution()
            .col_value
            .clone();

        if self.submip || self.mipdata_.incumbent.is_empty() {
            detach!(self.mipdata_.heuristics, |h| h
                .randomized_rounding(self, worker, &solvals));
            if self.options_mip_.mip_heuristic_run_shifting {
                detach!(self.mipdata_.heuristics, |h| h
                    .shifting(self, worker, &solvals));
            }
            detach!(self.mipdata_.heuristics, |h| h
                .flush_statistics(self, worker));
            *status = self.evaluate_root_lp(worker);
            if *status == HighsLpRelaxation_Status::kInfeasible {
                return true;
            }
        }

        false
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::evaluateRootLp (lines 1870-1965)
    /// C++: `HighsLpRelaxation::Status evaluateRootLp(HighsMipWorker & worker)`
    ///
    /// `worker`: index into `mipdata_.workers` (attached; heuristics take it as is).
    pub fn evaluate_root_lp(&mut self, worker: usize) -> HighsLpRelaxation_Status {
        let _lp_src_guard = crate::work::lp_kind(crate::work::LP_ROOT);
        loop {
            detach!(self.mipdata_.domains[0], |d| d.propagate(self));

            // rusolve mip_symmetry_dynamic: skip root global orbital fixing under the flag
            if !self.options_mip_.mip_symmetry_dynamic {
                if let Some(orbits) = self.mipdata_.global_orbits.clone() {
                    if !self.get_domain_2().infeasible() {
                        detach!(self.mipdata_.domains[0], |d| {
                            self.with_symmetries_ref(|ms, sym| orbits.orbital_fixing(ms, sym, d))
                        });
                    }
                }
            }

            if self.get_domain_2().infeasible() {
                let ub = cpp_min(kHighsInf, self.mipdata_.upper_bound);
                self.update_lower_bound(ub, true, true);
                self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                self.mipdata_.num_nodes += 1;
                self.mipdata_.num_leaves += 1;
                return HighsLpRelaxation_Status::kInfeasible;
            }

            let mut lp_bounds_changed = false;
            if !self.get_domain_2().get_changed_cols().is_empty() {
                lp_bounds_changed = true;
                self.remove_fixed_indices();
                detach!(self.mipdata_.lps[0], |lp| detach!(
                    self.mipdata_.domains[0],
                    |d| lp.flush_domain(self, d, false)
                ));
            }

            let mut lp_was_solved = false;
            let status: HighsLpRelaxation_Status;
            if lp_bounds_changed
                || self.get_lp_2().get_lp_solver_2().get_model_status() == HighsModelStatus::kNotset
            {
                let mut lp_iters: i64 = -self.get_lp_2().get_num_lp_iterations();
                let lp_src_guard = crate::work::lp_kind(crate::work::LP_ROOT);
                status = detach!(self.mipdata_.lps[0], |lp| detach!(
                    self.mipdata_.domains[0],
                    |d| lp.resolve_lp(self, Some(d))
                ));
                drop(lp_src_guard);
                lp_iters += self.get_lp_2().get_num_lp_iterations();
                self.mipdata_.total_lp_iterations += lp_iters;
                self.mipdata_.avgrootlpiters = self.get_lp().get_avg_solve_iters();
                lp_was_solved = true;

                if status == HighsLpRelaxation_Status::kUnbounded {
                    if self.solution_.is_empty() {
                        self.modelstatus_ = HighsModelStatus::kUnboundedOrInfeasible;
                    } else {
                        self.modelstatus_ = HighsModelStatus::kUnbounded;
                    }
                    self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                    self.mipdata_.num_nodes += 1;
                    self.mipdata_.num_leaves += 1;
                    return status;
                }

                if status == HighsLpRelaxation_Status::kOptimal
                    && self.get_lp_2().get_fractional_integers().is_empty()
                {
                    let sol = self
                        .get_lp_2()
                        .get_lp_solver_2()
                        .get_solution()
                        .col_value
                        .clone();
                    let obj = self.get_lp_2().get_objective();
                    if self.add_incumbent(&sol, obj, kSolutionSourceEvaluateNode, true, false) {
                        self.modelstatus_ = HighsModelStatus::kOptimal;
                        let ub = self.mipdata_.upper_bound;
                        self.update_lower_bound(ub, true, true);
                        self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                        self.mipdata_.num_nodes += 1;
                        self.mipdata_.num_leaves += 1;
                        return HighsLpRelaxation_Status::kInfeasible;
                    }
                }

                if status == HighsLpRelaxation_Status::kOptimal
                    && self.options_mip_.mip_heuristic_run_zi_round
                {
                    let sol = self
                        .get_lp_2()
                        .get_lp_solver_2()
                        .get_solution()
                        .col_value
                        .clone();
                    detach!(self.mipdata_.heuristics, |h| h.zi_round(self, worker, &sol));
                }
            } else {
                status = self.get_lp_2().get_status();
            }

            if status == HighsLpRelaxation_Status::kInfeasible {
                let ub = cpp_min(kHighsInf, self.mipdata_.upper_bound);
                self.update_lower_bound(ub, true, true);
                self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                self.mipdata_.num_nodes += 1;
                self.mipdata_.num_leaves += 1;
                return status;
            }

            if HighsLpRelaxation::unscaled_dual_feasible(self.get_lp_2().get_status()) {
                let obj = self.get_lp_2().get_objective();
                let lb = self.mipdata_.lower_bound;
                self.update_lower_bound(cpp_max(obj, lb), true, true);

                if lp_was_solved {
                    let col_dual = self
                        .get_lp_2()
                        .get_lp_solver_2()
                        .get_solution()
                        .col_dual
                        .clone();
                    let obj = self.get_lp_2().get_objective();
                    detach!(self.mipdata_.redcostfixing, |rc| rc
                        .add_root_redcost(self, &col_dual, obj));
                    if self.mipdata_.upper_limit != kHighsInf {
                        detach!(self.mipdata_.redcostfixing, |rc| rc
                            .propagate_root_redcost(self));
                    }
                }
            }

            if self.mipdata_.lower_bound > self.mipdata_.optimality_limit {
                self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                self.mipdata_.num_nodes += 1;
                self.mipdata_.num_leaves += 1;
                return HighsLpRelaxation_Status::kInfeasible;
            }

            if self.get_domain_2().get_changed_cols().is_empty() {
                return status;
            }
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::evaluateRootNode (lines 1987-2532)
    /// C++: `void evaluateRootNode(HighsMipWorker & worker)`
    ///
    /// Locals: `std::unique_ptr<SymmetryDetectionData> symData` ->
    /// `let mut sym_data: Option<Box<HighsMipSolverData_SymmetryDetectionData>> = None;`,
    /// `highs::parallel::TaskGroup tg` -> the body from line 1996 on runs inside
    /// `highs_parallel::task_group(|tg| ..)` (rule (g); declared before `restart:`, so one
    /// group across restarts, as in C++). `symData` is declared before `tg`, so it outlives
    /// the group, as in C++. Every `return` leaves the closure, and `task_group` then runs
    /// `~TaskGroup` (`cancel`, `taskWait`: a pending task never runs, a stolen one is
    /// waited for).
    /// `HighsSeparation sepa(worker); sepa.setLpRelaxation(&getLp());` ->
    /// `let mut sepa = HighsSeparation::new(self); sepa.set_lp_relaxation(Some(0));`.
    /// `goto restart` -> `loop { .. continue .. }` around the body from `restart:`.
    ///
    /// Profiling is a stub (mip_timer.rs): the `profiling->start/stop` calls, the
    /// `if (profiling->mip_)` timing messages and `clockOff` (a no-op) are omitted.
    /// `HIGHS_DEBUGSOL` is not defined in the reference build, so its block is omitted.
    /// `compute_analytic_centre` is `true` unless `mip_skip_analytic_centre` (rusolve, off), so its `printf` branch is dead.
    #[allow(clippy::needless_range_loop)]
    pub fn evaluate_root_node(&mut self, worker: usize) {
        let compute_analytic_centre = !self.options_mip_.mip_skip_analytic_centre;
        let mut max_sepa_rounds: HighsInt = if self.submip { 5 } else { kHighsIInf };
        if self.mipdata_.num_restarts == 0 {
            max_sepa_rounds = cpp_min(
                (2.0 * f64::from(self.mipdata_.max_tree_size_log2).sqrt()) as HighsInt,
                max_sepa_rounds,
            );
        }
        let mut sym_data: Option<Box<HighsMipSolverData_SymmetryDetectionData>> = None;
        crate::parallel::highs_parallel::task_group(|tg| {
            'restart: loop {
                let mut sym_early = false;
                let mut sym_stash: Option<SymmetryStash> = None;
                let sym_lazy = (self.options_mip_.mip_orbitope_root_fix
                    || self.options_mip_.mip_orbitope_full_fix)
                    && !self.submip;
                if self.mipdata_.detect_symmetries && !sym_lazy {
                    self.start_symmetry_detection(tg, &mut sym_data);
                }
                if compute_analytic_centre && !self.mipdata_.analytic_center_computed {
                    self.start_analytic_center_computation(tg);
                }

                // lp.getLpSolver().setOptionValue("dual_simplex_cost_perturbation_multiplier", 10.0);
                self.get_lp().set_iteration_limit(kHighsIInf);
                detach!(self.mipdata_.lps[0], |lp| lp.load_model(self));
                self.get_domain().clear_changed_cols();
                let upper_limit = self.mipdata_.upper_limit;
                detach!(self.mipdata_.lps[0], |lp| lp
                    .set_objective_limit(self, upper_limit));

                let lb = cpp_max(
                    self.mipdata_.lower_bound,
                    self.get_domain_2().get_objective_lower_bound(),
                );
                self.update_lower_bound(lb, true, true);

                self.print_display_line(kSolutionSourceNone);

                // Possibly query existence of an external solution
                if !self.submip {
                    let obj = self.solution_objective_;
                    self.query_external_solution(
                        obj,
                        kExternalMipSolutionQueryOriginEvaluateRootNode0,
                    );
                }

                // check if only root presolve is allowed
                if self.mipdata_.firstrootbasis.valid {
                    let basis = self.mipdata_.firstrootbasis.clone();
                    self.get_lp()
                        .get_lp_solver()
                        .set_basis(&basis, "HighsMipSolverData::evaluateRootNode");
                } else if self.options_mip_.mip_root_presolve_only {
                    self.get_lp()
                        .get_lp_solver()
                        .set_option_value_str("presolve", kHighsOffString);
                } else {
                    self.get_lp()
                        .get_lp_solver()
                        .set_option_value_str("presolve", kHighsOnString);
                }
                if self.options_mip_.highs_debug_level != 0 {
                    let output_flag = self.options_mip_.output_flag;
                    self.get_lp()
                        .get_lp_solver()
                        .set_option_value_bool("output_flag", output_flag);
                }

                let mut status = self.evaluate_root_lp(worker);
                if self.mipdata_.num_restarts == 0 {
                    self.mipdata_.firstrootlpiters = self.mipdata_.total_lp_iterations;
                }

                self.get_lp()
                    .get_lp_solver()
                    .set_option_value_bool("output_flag", false);
                self.get_lp()
                    .get_lp_solver()
                    .set_option_value_str("presolve", kHighsOffString);
                self.get_lp()
                    .get_lp_solver()
                    .set_option_value_str("parallel", kHighsOffString);

                if status == HighsLpRelaxation_Status::kInfeasible
                    || status == HighsLpRelaxation_Status::kUnbounded
                {
                    return;
                }

                self.mipdata_.firstlpsol = self.get_lp_2().get_solution().col_value.clone();
                self.mipdata_.firstlpsolobj = self.get_lp_2().get_objective();
                self.mipdata_.rootlpsolobj = self.mipdata_.firstlpsolobj;

                let lp_basis_valid = self.get_lp_2().get_lp_solver_2().get_basis().valid;
                if lp_basis_valid && self.get_lp_2().num_rows() == self.num_row() {
                    self.mipdata_.firstrootbasis =
                        self.get_lp_2().get_lp_solver_2().get_basis().clone();
                } else {
                    // the root basis is later expected to be consistent for the model without
                    // cuts so set it to the slack basis if the current basis already includes
                    // cuts, e.g. due to a restart
                    let ncol = self.num_col() as usize;
                    let nrow = self.num_row() as usize;
                    let fb = &mut self.mipdata_.firstrootbasis;
                    fb.col_status = vec![HighsBasisStatus::kNonbasic; ncol];
                    fb.row_status = vec![HighsBasisStatus::kBasic; nrow];
                    fb.valid = true;
                    fb.useful = true;
                }

                // rusolve mip_orbitope_root_fix: detection runs
                // inline and lazily, so it is finished here, after the first root LP and only
                // when that LP is fractional (models solved at the LP never pay for it).
                if sym_lazy
                    && self.mipdata_.detect_symmetries
                    && !self.get_lp_2().get_fractional_integers().is_empty()
                {
                    self.start_symmetry_detection(tg, &mut sym_data);
                }
                if sym_lazy && self.mipdata_.detect_symmetries && sym_data.is_some() {
                    self.finish_symmetry_detection(tg, &mut sym_data);
                    sym_early = true;
                    let mut nfixed = 0;
                    if self.options_mip_.mip_orbitope_full_fix {
                        nfixed += self.orbitope_full_fix();
                    }
                    if self.options_mip_.mip_orbitope_root_fix
                        && !self.mipdata_.domains[0].infeasible()
                    {
                        nfixed += self.orbitope_root_fix();
                    }
                    if nfixed == 0 && self.mipdata_.detect_symmetries {
                        // nothing fixed: put the symmetry data back until the point where HiGHS
                        // finishes detection, so models without a fixable orbitope keep the
                        // HiGHS search bit for bit (identity.sh)
                        sym_early = false;
                        sym_stash = Some(self.stash_symmetries());
                    }
                    if self.mipdata_.domains[0].infeasible() {
                        self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                        self.modelstatus_ = HighsModelStatus::kInfeasible;
                        return;
                    }
                    if nfixed > 0 {
                        status = self.evaluate_root_lp(worker);
                        if status == HighsLpRelaxation_Status::kInfeasible {
                            return;
                        }
                        self.mipdata_.firstlpsol = self.get_lp_2().get_solution().col_value.clone();
                        self.mipdata_.firstlpsolobj = self.get_lp_2().get_objective();
                        self.mipdata_.rootlpsolobj = self.mipdata_.firstlpsolobj;
                    }
                }

                if self.get_cut_pool_2().get_num_cuts() != 0 {
                    debug_assert!(self.mipdata_.num_restarts != 0);
                    let mut cutset = HighsCutSet::default();
                    self.mipdata_.cutpools[0].separate_lp_cuts_after_restart(&mut cutset);
                    // HIGHS_DEBUGSOL block (debugSolution.checkCut per cut) omitted: not defined.
                    self.mipdata_.lps[0].add_cuts(&mut cutset);
                    status = self.evaluate_root_lp(worker);
                    detach!(self.mipdata_.lps[0], |lp| lp
                        .remove_obsolete_rows(self, true));
                    if status == HighsLpRelaxation_Status::kInfeasible {
                        return;
                    }
                }

                let iter_limit = cpp_max(10000, (10.0 * self.mipdata_.avgrootlpiters) as HighsInt);
                self.get_lp().set_iteration_limit(iter_limit);

                // make sure first line after solving root LP is printed
                self.mipdata_.last_disptime = -kHighsInf;
                self.mipdata_.disptime = 0.0;

                let firstlpsol = self.mipdata_.firstlpsol.clone();
                if self.options_mip_.mip_heuristic_run_zi_round {
                    detach!(self.mipdata_.heuristics, |h| h.zi_round(
                        self,
                        worker,
                        &firstlpsol
                    ));
                }
                detach!(self.mipdata_.heuristics, |h| h.randomized_rounding(
                    self,
                    worker,
                    &firstlpsol
                ));
                if self.options_mip_.mip_heuristic_run_shifting {
                    detach!(self.mipdata_.heuristics, |h| h.shifting(
                        self,
                        worker,
                        &firstlpsol
                    ));
                }
                detach!(self.mipdata_.heuristics, |h| h
                    .flush_statistics(self, worker));

                status = self.evaluate_root_lp(worker);
                if status == HighsLpRelaxation_Status::kInfeasible {
                    return;
                }

                // rusolve mip_first_inc_root: still no incumbent after the
                // first root LP and its rounding heuristics, before the cut loop (which can take
                // the whole time limit on instances without an incumbent)
                if self.options_mip_.mip_first_inc_root
                    && !self.submip
                    && self.mipdata_.upper_limit == kHighsInf
                    && !self.check_limits(0)
                {
                    self.first_incumbent_fallback();
                }

                self.mipdata_.rootlpsolobj = self.mipdata_.firstlpsolobj;
                self.remove_fixed_indices();
                if self.options_mip_.mip_allow_restart
                    && self.options_mip_.presolve != kHighsOffString
                {
                    let fixing_rate = self.percentage_inactive_integers();
                    if fixing_rate >= 10.0 || self.cutoff_bigm_restart_due() {
                        tg.cancel();
                        highs_log_user(
                            &self.options_mip_.log_options,
                            HighsLogType::kInfo,
                            &format!(
                                "\n{:.1}% inactive integer columns, restarting\n",
                                fixing_rate
                            ),
                        );
                        tg.task_wait();
                        self.perform_restart();
                        self.mipdata_.num_restarts_root += 1;
                        if self.modelstatus_ == HighsModelStatus::kNotset {
                            continue 'restart;
                        }
                        return;
                    }
                }

                // begin separation (the MIP-Timing messages are profiling-only and omitted)
                let ncol = self.num_col() as usize;
                let mut avgdirection: Vec<f64> = vec![0.0; ncol];
                let mut curdirection: Vec<f64> = vec![0.0; ncol];

                let mut stall: HighsInt = 0;
                let mut smoothprogress: f64 = 0.0;
                let mut nseparounds: HighsInt = 0;
                // rusolve: flat root bound with an incumbent
                let mut bound_flat: HighsInt = 0;
                let mut bound_prev: f64 = self.get_lp_2().get_objective();
                let mut sepa = HighsSeparation::new(self);
                sepa.set_lp_relaxation(Some(0));

                while HighsLpRelaxation::scaled_optimal(status)
                    && !self.get_lp_2().get_fractional_integers().is_empty()
                    && stall
                        < (if self.options_mip_.mip_root_stall_scip {
                            10
                        } else {
                            3
                        })
                {
                    self.print_display_line(kSolutionSourceNone);

                    if self.check_limits(0) {
                        return;
                    }

                    if nseparounds == max_sepa_rounds {
                        break;
                    }

                    self.remove_fixed_indices();

                    if !self.submip && self.options_mip_.presolve != kHighsOffString {
                        let fixing_rate = self.percentage_inactive_integers();
                        if fixing_rate >= 10.0 {
                            stall = -1;
                            break;
                        }
                    }

                    nseparounds += 1;

                    let mut ncuts: HighsInt = 0;
                    let root_separation_round_result =
                        self.root_separation_round(worker, &mut sepa, &mut ncuts, &mut status);
                    if root_separation_round_result {
                        return;
                    }

                    if nseparounds >= 5
                        && !self.submip
                        && !self.mipdata_.analytic_center_computed
                        && compute_analytic_centre
                    {
                        if self.check_limits(0) {
                            return;
                        }
                        self.finish_analytic_center_computation(tg);

                        detach!(self.mipdata_.heuristics, |h| h
                            .central_rounding(self, worker));
                        detach!(self.mipdata_.heuristics, |h| h
                            .flush_statistics(self, worker));

                        if self.check_limits(0) {
                            return;
                        }
                        status = self.evaluate_root_lp(worker);
                        if status == HighsLpRelaxation_Status::kInfeasible {
                            return;
                        }
                    }

                    let mut sqrnorm = HighsCDouble::from(0.0);
                    let solvals = self.get_lp_2().get_solution().col_value.clone();

                    for i in 0..ncol {
                        curdirection[i] = self.mipdata_.firstlpsol[i] - solvals[i];

                        // if (mip.integrality_[i] == 2 && lp.getObjective() > firstobj &&
                        //    std::abs(curdirection[i]) > 1e-6)
                        //  pseudocost.addObservation(i, -curdirection[i],
                        //                            lp.getObjective() - firstobj);

                        sqrnorm += curdirection[i] * curdirection[i];
                    }

                    let scale = f64::from(1.0 / sqrnorm.sqrt());
                    sqrnorm = HighsCDouble::from(0.0);
                    let mut dotproduct = HighsCDouble::from(0.0);
                    for i in 0..ncol {
                        avgdirection[i] =
                            (scale * curdirection[i] - avgdirection[i]) / f64::from(nseparounds);
                        sqrnorm += avgdirection[i] * avgdirection[i];
                        dotproduct += avgdirection[i] * curdirection[i];
                    }

                    // HiGHS quirk: Apple clang (arm64, -ffp-contract=on) may fuse the
                    // a*b+c below into an FMA; this port does not (no mul_add).
                    let progress = f64::from(dotproduct / sqrnorm.sqrt());

                    if nseparounds == 1 {
                        smoothprogress = progress;
                    } else {
                        let alpha = 1.0 / 3.0;
                        let nextprogress = (1.0 - alpha) * smoothprogress + alpha * progress;

                        if nextprogress < smoothprogress * 1.01
                            && (self.get_lp_2().get_objective() - self.mipdata_.firstlpsolobj)
                                <= (self.mipdata_.rootlpsolobj - self.mipdata_.firstlpsolobj)
                                    * 1.001
                        {
                            stall += 1;
                        } else {
                            stall = 0;
                        }
                        smoothprogress = nextprogress;
                    }

                    self.mipdata_.rootlpsolobj = self.get_lp_2().get_objective();
                    let iter_limit =
                        cpp_max(10000, (10.0 * self.mipdata_.avgrootlpiters) as HighsInt);
                    self.get_lp().set_iteration_limit(iter_limit);
                    if ncuts == 0 {
                        break;
                    }
                    // rusolve `mip_root_bound_stall`: with an
                    // incumbent, cuts that leave the root bound
                    // unchanged for 5 rounds cannot close the gap; hand over to the root
                    // heuristics (RENS, sub-MIP) instead of separating on
                    if self.options_mip_.mip_root_bound_stall
                        && !self.submip
                        && self.mipdata_.upper_limit != kHighsInf
                    {
                        let cur = self.mipdata_.rootlpsolobj;
                        // gap to the incumbent must be small (decomp2: 18 % gap, the late
                        // rounds led to the optimum; neos-1171448: 0.3 %)
                        let small_gap = self.mipdata_.upper_bound - cur
                            <= 0.02 * cur.abs().max(self.mipdata_.upper_bound.abs()).max(1e-9);
                        if (cur - bound_prev).abs() <= 1e-9 * (1.0 + cur.abs()) {
                            bound_flat += 1;
                        } else {
                            bound_flat = 0;
                        }
                        bound_prev = cur;
                        if bound_flat >= 5 && small_gap {
                            // Log line only (no state): lets a run say whether the lever fired.
                            highs_log_user(
                                &self.options_mip_.log_options,
                                HighsLogType::kInfo,
                                "rusolve: root bound stall, leaving the cut loop\n",
                            );
                            break;
                        }
                    }

                    // Possibly query existence of an external solution
                    if !self.submip {
                        let obj = self.solution_objective_;
                        self.query_external_solution(
                            obj,
                            kExternalMipSolutionQueryOriginEvaluateRootNode1,
                        );
                    }
                }

                self.get_lp().set_iteration_limit(kHighsIInf);
                status = self.evaluate_root_lp(worker);
                if status == HighsLpRelaxation_Status::kInfeasible {
                    return;
                }

                self.mipdata_.rootlpsol = self
                    .get_lp_2()
                    .get_lp_solver_2()
                    .get_solution()
                    .col_value
                    .clone();
                self.mipdata_.rootlpsolobj = self.get_lp_2().get_objective();
                let iter_limit = cpp_max(10000, (10.0 * self.mipdata_.avgrootlpiters) as HighsInt);
                self.get_lp().set_iteration_limit(iter_limit);

                if self.options_mip_.mip_heuristic_run_zi_round {
                    detach!(self.mipdata_.heuristics, |h| h.zi_round(
                        self,
                        worker,
                        &firstlpsol
                    ));
                    detach!(self.mipdata_.heuristics, |h| h
                        .flush_statistics(self, worker));
                }
                if self.options_mip_.mip_heuristic_run_shifting {
                    let rootlpsol = self.mipdata_.rootlpsol.clone();
                    detach!(self.mipdata_.heuristics, |h| h
                        .shifting(self, worker, &rootlpsol));
                    detach!(self.mipdata_.heuristics, |h| h
                        .flush_statistics(self, worker));
                }

                if !self.mipdata_.analytic_center_computed && compute_analytic_centre {
                    if self.check_limits(0) {
                        return;
                    }

                    self.finish_analytic_center_computation(tg);

                    detach!(self.mipdata_.heuristics, |h| h
                        .central_rounding(self, worker));
                    detach!(self.mipdata_.heuristics, |h| h
                        .flush_statistics(self, worker));

                    // if there are new global bound changes we re-evaluate the LP and do one
                    // more separation round
                    if self.check_limits(0) {
                        return;
                    }
                    let separate = !self.get_domain_2().get_changed_cols().is_empty();
                    status = self.evaluate_root_lp(worker);
                    if status == HighsLpRelaxation_Status::kInfeasible {
                        return;
                    }
                    if separate && HighsLpRelaxation::scaled_optimal(status) {
                        let mut ncuts: HighsInt = 0;
                        let root_separation_round_result =
                            self.root_separation_round(worker, &mut sepa, &mut ncuts, &mut status);
                        if root_separation_round_result {
                            return;
                        }
                        nseparounds += 1;
                        self.print_display_line(kSolutionSourceNone);
                    }
                }

                self.print_display_line(kSolutionSourceNone);
                // Possibly query existence of an external solution
                if !self.submip {
                    let obj = self.solution_objective_;
                    self.query_external_solution(
                        obj,
                        kExternalMipSolutionQueryOriginEvaluateRootNode2,
                    );
                }

                // Possible cut extraction callback
                if !self.submip
                    && self.callback_.user_callback.is_some()
                    && self.callback_.callback_active(kCallbackMipGetCutPool)
                {
                    self.callback_get_cut_pool();
                }
                if self.check_limits(0) {
                    return;
                }

                // C++ `do { ... } while (false);` with `break`: a labelled block.
                'do_while: {
                    if self.mipdata_.rootlpsol.is_empty() {
                        break 'do_while;
                    }
                    if (self.mipdata_.upper_limit != kHighsInf
                        || self.options_mip_.mip_heur_lp_quota)
                        && !self.more_heuristics_allowed()
                    {
                        break 'do_while;
                    }

                    if self.options_mip_.mip_heuristic_run_root_reduced_cost {
                        detach!(self.mipdata_.heuristics, |h| h
                            .root_reduced_cost(self, worker));
                        detach!(self.mipdata_.heuristics, |h| h
                            .flush_statistics(self, worker));
                    }

                    if self.check_limits(0) {
                        return;
                    }

                    // if there are new global bound changes we re-evaluate the LP and do one
                    // more separation round
                    let mut separate = !self.get_domain_2().get_changed_cols().is_empty();
                    status = self.evaluate_root_lp(worker);
                    if status == HighsLpRelaxation_Status::kInfeasible {
                        return;
                    }
                    if separate && HighsLpRelaxation::scaled_optimal(status) {
                        let mut ncuts: HighsInt = 0;
                        let root_separation_round_result =
                            self.root_separation_round(worker, &mut sepa, &mut ncuts, &mut status);
                        if root_separation_round_result {
                            return;
                        }
                        nseparounds += 1;
                        self.print_display_line(kSolutionSourceNone);
                    }

                    if (self.mipdata_.upper_limit != kHighsInf
                        || self.options_mip_.mip_heur_lp_quota)
                        && !self.more_heuristics_allowed()
                    {
                        break 'do_while;
                    }

                    if self.check_limits(0) {
                        return;
                    }
                    if self.options_mip_.mip_heuristic_run_rens {
                        let rootlpsol = self.mipdata_.rootlpsol.clone();
                        detach!(self.mipdata_.heuristics, |h| h
                            .rens(self, worker, &rootlpsol));
                        detach!(self.mipdata_.heuristics, |h| h
                            .flush_statistics(self, worker));
                    }

                    if self.check_limits(0) {
                        return;
                    }
                    // if there are new global bound changes we re-evaluate the LP and do one
                    // more separation round
                    separate = !self.get_domain_2().get_changed_cols().is_empty();
                    status = self.evaluate_root_lp(worker);
                    if status == HighsLpRelaxation_Status::kInfeasible {
                        return;
                    }
                    if separate && HighsLpRelaxation::scaled_optimal(status) {
                        let mut ncuts: HighsInt = 0;
                        let root_separation_round_result =
                            self.root_separation_round(worker, &mut sepa, &mut ncuts, &mut status);
                        if root_separation_round_result {
                            return;
                        }
                        nseparounds += 1;

                        self.print_display_line(kSolutionSourceNone);
                        // Possibly query existence of an external solution
                        if !self.submip {
                            let obj = self.solution_objective_;
                            self.query_external_solution(
                                obj,
                                kExternalMipSolutionQueryOriginEvaluateRootNode3,
                            );
                        }
                    }

                    if self.mipdata_.upper_limit != kHighsInf || self.submip {
                        break 'do_while;
                    }

                    if self.check_limits(0) {
                        return;
                    }
                    detach!(self.mipdata_.heuristics, |h| h
                        .feasibility_pump(self, worker));
                    detach!(self.mipdata_.heuristics, |h| h
                        .flush_statistics(self, worker));

                    if self.check_limits(0) {
                        return;
                    }
                    status = self.evaluate_root_lp(worker);
                    if status == HighsLpRelaxation_Status::kInfeasible {
                        return;
                    }
                }

                if self.mipdata_.lower_bound > self.mipdata_.upper_limit {
                    self.modelstatus_ = HighsModelStatus::kOptimal;
                    self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                    self.mipdata_.num_nodes += 1;
                    self.mipdata_.num_leaves += 1;
                    return;
                }

                // if there are new global bound changes we re-evaluate the LP and do one
                // more separation round
                let separate = !self.get_domain_2().get_changed_cols().is_empty();
                status = self.evaluate_root_lp(worker);
                if status == HighsLpRelaxation_Status::kInfeasible {
                    return;
                }
                if separate && HighsLpRelaxation::scaled_optimal(status) {
                    let mut ncuts: HighsInt = 0;
                    let root_separation_round_result =
                        self.root_separation_round(worker, &mut sepa, &mut ncuts, &mut status);
                    if root_separation_round_result {
                        return;
                    }
                    nseparounds += 1;
                    self.print_display_line(kSolutionSourceNone);
                }

                // Possibly query existence of an external solution
                if !self.submip {
                    let obj = self.solution_objective_;
                    self.query_external_solution(
                        obj,
                        kExternalMipSolutionQueryOriginEvaluateRootNode4,
                    );
                }

                self.remove_fixed_indices();
                if self.get_lp_2().get_lp_solver_2().get_basis().valid {
                    detach!(self.mipdata_.lps[0], |lp| lp
                        .remove_obsolete_rows(self, true));
                }
                self.mipdata_.rootlpsolobj = self.get_lp_2().get_objective();

                self.print_display_line(kSolutionSourceNone);

                if self.mipdata_.lower_bound <= self.mipdata_.upper_limit {
                    if !self.submip
                        && self.options_mip_.mip_allow_restart
                        && self.options_mip_.presolve != kHighsOffString
                    {
                        if !self.mipdata_.analytic_center_computed && compute_analytic_centre {
                            self.finish_analytic_center_computation(tg);
                        }
                        let fixing_rate = self.percentage_inactive_integers();
                        let submip_weight = if self.submip { 1.0 } else { 0.0 };
                        if fixing_rate >= 2.5 + 7.5 * submip_weight
                            || (!self.submip
                                && fixing_rate > 0.0
                                && self.mipdata_.num_restarts == 0)
                            || self.cutoff_bigm_restart_due()
                        {
                            tg.cancel();
                            highs_log_user(
                                &self.options_mip_.log_options,
                                HighsLogType::kInfo,
                                &format!(
                                    "\n{:.1}% inactive integer columns, restarting\n",
                                    fixing_rate
                                ),
                            );
                            if stall != -1 {
                                max_sepa_rounds = cpp_min(max_sepa_rounds, nseparounds);
                            }
                            tg.task_wait();
                            self.perform_restart();
                            if self.terminate() {
                                return;
                            }
                            self.mipdata_.num_restarts_root += 1;
                            if self.modelstatus_ == HighsModelStatus::kNotset {
                                continue 'restart;
                            }
                            return;
                        }
                    }

                    if self.mipdata_.detect_symmetries && !sym_early {
                        if let Some(stash) = sym_stash.take() {
                            self.restore_symmetries(stash);
                        } else {
                            if sym_data.is_none() {
                                self.start_symmetry_detection(tg, &mut sym_data);
                            }
                            self.finish_symmetry_detection(tg, &mut sym_data);
                        }
                        status = self.evaluate_root_lp(worker);
                        if status == HighsLpRelaxation_Status::kInfeasible {
                            return;
                        }
                    }

                    // add the root node to the nodequeue to initialize the search
                    let ps_index = self.mipdata_.workers[worker].get_pseudocost();
                    let estimate = self
                        .get_lp_2()
                        .compute_best_estimate(self, &self.mipdata_.pseudocosts[ps_index]);
                    let lower_bound = self.mipdata_.lower_bound;
                    // rusolve mip_symmetry_dynamic: root node carries no branching path
                    if self.options_mip_.mip_symmetry_dynamic {
                        self.mipdata_.nodequeue.emplace_node_s3(
                            Vec::new(),
                            Vec::new(),
                            lower_bound,
                            estimate,
                            1,
                            None,
                        );
                    } else {
                        self.mipdata_.nodequeue.emplace_node(
                            Vec::new(),
                            Vec::new(),
                            lower_bound,
                            estimate,
                            1,
                        );
                    }
                }
                // End of HighsMipSolverData::evaluateRootNode()
                break 'restart;
            }
        })
    }

    /// Glue (ours): C++ `sepa.separationRound(getDomain(), status)` from the root
    /// (cpp 1846). Detaches `workers[worker]`, the separation LP `lps[0]` (`getLp()`, the
    /// one `sepa.setLpRelaxation(&getLp())` installed) and the propagation domain
    /// `domains[0]` (`getDomain()`), runs `separation_round` with
    /// `SepaCtx { ms, worker, lp }`, and puts all three back. Its body lends them back
    /// with `lend_in!` where it re-enters through `ms`.
    pub fn separation_round_global(
        &mut self,
        worker: usize,
        sepa: &mut HighsSeparation,
        status: &mut HighsLpRelaxation_Status,
    ) -> HighsInt {
        let result = detach!(self.mipdata_.workers[worker], |wk| {
            detach!(self.mipdata_.lps[0], |lp| {
                detach!(self.mipdata_.domains[0], |dom| {
                    let mut ctx = SepaCtx {
                        ms: &mut *self,
                        worker: &mut *wk,
                        lp: &mut *lp,
                    };
                    sepa.separation_round(&mut ctx, dom, status)
                })
            })
        });
        // rusolve mip_ti_list_sched: schedules built inside a separator
        for x in sepa.take_pending_solutions() {
            self.try_solution(&x, kSolutionSourceHeuristic);
        }
        result
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::checkLimits (lines 2534-2630)
    /// C++: `bool checkLimits(int64_t nodeOffset = 0) const`
    /// C++ default arguments: `node_offset = 0`
    ///
    /// PORT NOTE (architecture): the C++ method is `const` but writes
    /// `mipsolver.modelstatus_` and the callback (`clearHighsCallbackOutput`,
    /// `interruptFromCallbackWithData`) through the `HighsMipSolver&` member, so the final
    /// receiver is `&mut self`. HighsSearch::check_limits (search/part2.rs:172) takes
    /// `ctx: &mut SearchCtx` to reach it.
    pub fn check_limits(&mut self, node_offset: i64) -> bool {
        // rusolve racing: another racer finished
        if let Some(race) = &self.race {
            if race.stop.load(std::sync::atomic::Ordering::Acquire) {
                if self.modelstatus_ == HighsModelStatus::kNotset {
                    self.modelstatus_ = HighsModelStatus::kInterrupt;
                }
                return true;
            }
        }
        // This MIP instance may have been terminated
        if self.terminator_active() && self.terminator_terminated() {
            return true;
        }

        // Possible user interrupt
        if !self.submip && !self.parallel_lock_active() && self.callback_.user_callback.is_some() {
            self.callback_.clear_highs_callback_output();
            let obj = self.solution_objective_;
            if self.interrupt_from_callback_with_data(
                kCallbackMipInterrupt,
                obj,
                "MIP check limits",
            ) {
                if self.modelstatus_ == HighsModelStatus::kNotset {
                    highs_log_dev(
                        &self.options_mip_.log_options,
                        HighsLogType::kInfo,
                        "User interrupt\n",
                    );
                    self.modelstatus_ = HighsModelStatus::kInterrupt;
                }
                return true;
            }
        }

        // Possible termination due to objective being at least as good as
        // the target value
        if !self.submip
            && self.solution_objective_ < kHighsInf
            && self.options_mip_.objective_target > -kHighsInf
        {
            // Note:
            //
            // Whether the sense is ObjSense::kMinimize or
            // ObjSense::kMaximize, the undefined value of
            // mipsolver.solution_objective_ is kHighsInf, and the default
            // target value is -kHighsInf, so had to rule out these cases in
            // the conditional statement above.
            //
            // mipsolver.solution_objective_ is the actual objective of the
            // MIP - including the offset, and independent of objective sense
            //
            // The target is reached if the objective is below (above) the
            // target value when minimizing (maximizing).
            let int_sense = self.orig_model_.sense_ as HighsInt;
            let reached_objective_target = (int_sense as f64) * self.solution_objective_
                < (int_sense as f64) * self.options_mip_.objective_target;
            if reached_objective_target {
                if self.modelstatus_ == HighsModelStatus::kNotset {
                    highs_log_dev(
                        &self.options_mip_.log_options,
                        HighsLogType::kInfo,
                        "Reached objective target\n",
                    );
                    self.modelstatus_ = HighsModelStatus::kObjectiveTarget;
                }
                return true;
            }
        }

        if self.options_mip_.mip_max_nodes != kHighsIInf
            && self.mipdata_.num_nodes + node_offset >= self.options_mip_.mip_max_nodes as i64
        {
            if self.modelstatus_ == HighsModelStatus::kNotset {
                highs_log_dev(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "Reached node limit\n",
                );
                self.modelstatus_ = HighsModelStatus::kSolutionLimit;
            }
            return true;
        }

        if self.options_mip_.mip_max_leaves != kHighsIInf
            && self.mipdata_.num_leaves >= self.options_mip_.mip_max_leaves as i64
        {
            if self.modelstatus_ == HighsModelStatus::kNotset {
                highs_log_dev(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "Reached leaf node limit\n",
                );
                self.modelstatus_ = HighsModelStatus::kSolutionLimit;
            }
            return true;
        }

        if self.options_mip_.mip_max_improving_sols != kHighsIInf
            && self.mipdata_.num_improving_sols >= self.options_mip_.mip_max_improving_sols
        {
            if self.modelstatus_ == HighsModelStatus::kNotset {
                highs_log_dev(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "Reached improving solution limit\n",
                );
                self.modelstatus_ = HighsModelStatus::kSolutionLimit;
            }
            return true;
        }

        if let Some(limit) = self.work_limit
            && crate::work::WorkSnapshot::now()
                .since(&self.stats_.work_start)
                .work()
                >= limit
        {
            if self.modelstatus_ == HighsModelStatus::kNotset {
                highs_log_dev(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "Reached work limit\n",
                );
                self.modelstatus_ = HighsModelStatus::kSolutionLimit;
            }
            return true;
        }

        //  const double time = mipsolver.timer_.read();
        //  printf("checkLimits: time = %g\n", time);
        if self.options_mip_.time_limit < kHighsInf
            && self.timer_.read() >= self.options_mip_.time_limit
        {
            if self.modelstatus_ == HighsModelStatus::kNotset {
                highs_log_dev(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "Reached time limit\n",
                );
                self.modelstatus_ = HighsModelStatus::kTimeLimit;
            }
            return true;
        }

        false
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::checkObjIntegrality (lines 2632-2639)
    /// C++: `void checkObjIntegrality()`
    pub fn check_obj_integrality(&mut self) {
        let epsilon = self.mipdata_.epsilon;
        self.mipdata_.objective_function.check_integrality(epsilon);
        if self.mipdata_.objective_function.is_integral() && self.mipdata_.num_restarts == 0 {
            let scale = self.mipdata_.objective_function.integral_scale();
            highs_log_user(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                &format!(
                    "Objective function is integral with scale {}\n",
                    c_fmt_g(scale, 6)
                ),
            );
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::setupDomainPropagation (lines 2641-2664)
    /// C++: `void setupDomainPropagation()`
    #[allow(clippy::needless_range_loop)]
    pub fn setup_domain_propagation(&mut self) {
        let (ar_start, ar_index, ar_value) = {
            let model = mip_model!(self);
            highs_sparse_transpose(
                model.num_row_,
                model.num_col_,
                &model.a_matrix_.start_,
                &model.a_matrix_.index_,
                &model.a_matrix_.value_,
            )
        };
        self.mipdata_.a_rstart_ = ar_start;
        self.mipdata_.a_rindex_ = ar_index;
        self.mipdata_.a_rvalue_ = ar_value;

        let pseudocost = HighsPseudocost::new(self);
        *self.get_pseudo_cost() = pseudocost;

        // compute the maximal absolute coefficients to filter propagation
        let num_row = self.num_row() as usize;
        self.mipdata_.max_abs_row_coef.resize(num_row, 0.0);
        for i in 0..num_row {
            let mut maxabsval = 0.0f64;

            let start = self.mipdata_.a_rstart_[i] as usize;
            let end = self.mipdata_.a_rstart_[i + 1] as usize;
            for j in start..end {
                maxabsval = cpp_max(maxabsval, self.mipdata_.a_rvalue_[j].abs());
            }

            self.mipdata_.max_abs_row_coef[i] = maxabsval;
        }

        // C++ `getDomain() = HighsDomain(mipsolver)`: the copy assignment drops the
        // registrations of the old domain with its pools (`release`), then takes the
        // fresh domain (the temporary registers with no pool).
        let fresh = HighsDomain::new(self);
        detach!(self.mipdata_.domains[0], |d| d.release(self));
        *self.mipdata_.domains[0] = fresh;
        detach!(self.mipdata_.domains[0], |d| d.compute_row_activities(self));
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::saveReportMipSolution (lines 2666-2697)
    /// C++: `void saveReportMipSolution(const double new_upper_limit = - kHighsInf)`
    /// C++ default arguments: `new_upper_limit = - kHighsInf`
    pub fn save_report_mip_solution(&mut self, new_upper_limit: f64) {
        let non_improving = new_upper_limit >= self.mipdata_.upper_limit;
        if self.submip {
            return;
        }
        if non_improving {
            return;
        }

        if self.callback_.user_callback.is_some()
            && self.callback_.active[kCallbackMipImprovingSolution as usize]
        {
            self.callback_.clear_highs_callback_output();
            self.callback_.data_out.mip_solution = self.solution_.clone();
            let obj = self.solution_objective_;
            let interrupt = self.interrupt_from_callback_with_data(
                kCallbackMipImprovingSolution,
                obj,
                "Improving solution",
            );
            debug_assert!(!interrupt);
            // the source belongs to this one incumbent
            self.callback_.last_solution_source = -1;
        }

        if self.options_mip_.mip_improving_solution_save {
            let record = HighsObjectiveSolution {
                objective: self.solution_objective_,
                col_value: self.solution_.clone(),
            };
            self.saved_objective_and_solution_.push(record);
        }
        if let Some(file) = self.improving_solution_file_.clone() {
            // `FILE* file = mipsolver.improving_solution_file_`
            let mut file: &std::fs::File = &file;
            write_lp_objective(
                &mut file,
                &self.options_mip_.log_options,
                &self.orig_model_,
                &self.solution_,
            );
            write_primal_solution(
                &mut file,
                &self.options_mip_.log_options,
                &self.orig_model_,
                &self.solution_,
                self.options_mip_.mip_improving_solution_report_sparse,
            );
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::limitsToBounds (lines 2699-2710)
    /// C++: `void limitsToBounds(double & dual_bound, double & primal_bound, double & mip_rel_gap) const`
    pub fn limits_to_bounds(
        &self,
        dual_bound: &mut f64,
        primal_bound: &mut f64,
        mip_rel_gap: &mut f64,
    ) {
        *mip_rel_gap = self.limits_to_gap(
            self.mipdata_.lower_bound,
            self.mipdata_.upper_bound,
            dual_bound,
            primal_bound,
        );
        *primal_bound = cpp_min(self.options_mip_.objective_bound, *primal_bound);
        // Adjust objective sense in case of maximization problem
        if matches!(self.orig_model_.sense_, ObjSense::kMaximize) {
            *dual_bound = -*dual_bound;
            *primal_bound = -*primal_bound;
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::updateLowerBound (lines 2712-2721)
    /// C++: `void updateLowerBound(double new_lower_bound, const bool check_bound_change = true, const bool check_prev_data = true)`
    /// C++ default arguments: `check_bound_change = true`, `check_prev_data = true`
    pub fn update_lower_bound(
        &mut self,
        new_lower_bound: f64,
        check_bound_change: bool,
        check_prev_data: bool,
    ) {
        // Update lower bound
        let prev_lower_bound = self.mipdata_.lower_bound;
        self.mipdata_.lower_bound = new_lower_bound;
        if !self.submip && self.mipdata_.lower_bound != prev_lower_bound {
            let lower_bound = self.mipdata_.lower_bound;
            let upper_bound = self.mipdata_.upper_bound;
            self.update_primal_dual_integral(
                prev_lower_bound,
                lower_bound,
                upper_bound,
                upper_bound,
                check_bound_change,
                check_prev_data,
            );
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::setCallbackDataOut (lines 2727-2742)
    /// C++: `void setCallbackDataOut(const double mipsolver_objective_value) const`
    ///
    /// `&mut self`: writes `callback_->data_out` (C++ `const`, through the pointer).
    pub fn set_callback_data_out(&mut self, mipsolver_objective_value: f64) {
        let mut dual_bound = 0.0;
        let mut primal_bound = 0.0;
        let mut mip_rel_gap = 0.0;
        self.limits_to_bounds(&mut dual_bound, &mut primal_bound, &mut mip_rel_gap);
        self.callback_.data_out.running_time = self.timer_.read();
        self.callback_.data_out.objective_function_value = mipsolver_objective_value;
        self.callback_.data_out.mip_node_count = self.mipdata_.num_nodes;
        self.callback_.data_out.mip_total_lp_iterations = self.mipdata_.total_lp_iterations;
        self.callback_.data_out.mip_primal_bound = primal_bound;
        self.callback_.data_out.mip_dual_bound = dual_bound;
        self.callback_.data_out.mip_gap = mip_rel_gap;
        // information for the callback only (see `HighsCallbackOutput`)
        self.callback_.data_out.mip_solution_source = self.callback_.last_solution_source;
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::interruptFromCallbackWithData (lines 2744-2751)
    /// C++: `bool interruptFromCallbackWithData(const int callback_type, const double mipsolver_objective_value, const std::string message = "") const`
    /// C++ default arguments: `message = ""`
    ///
    /// `&mut self`: runs `callback_->callbackAction` (C++ `const`, through the pointer).
    pub fn interrupt_from_callback_with_data(
        &mut self,
        callback_type: i32,
        mipsolver_objective_value: f64,
        message: &str,
    ) -> bool {
        if !self.callback_.callback_active(callback_type) {
            return false;
        }
        debug_assert!(!self.submip);
        self.set_callback_data_out(mipsolver_objective_value);
        self.callback_.callback_action(callback_type, message)
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::queryExternalSolution (lines 2753-2817)
    /// HiGHS 1.16 (#3492): HighsMipSolverData.cpp::queryExternalSolution (lines 2813-2877)
    /// C++: `void queryExternalSolution(const double mipsolver_objective_value, const ExternalMipSolutionQueryOrigin external_solution_query_origin)`
    pub fn query_external_solution(
        &mut self,
        mipsolver_objective_value: f64,
        external_solution_query_origin: ExternalMipSolutionQueryOrigin,
    ) {
        debug_assert!(!self.submip);
        let use_callback = self.callback_.user_callback.is_some()
            && self.callback_.active[kCallbackMipUserSolution as usize];
        if use_callback {
            self.set_callback_data_out(mipsolver_objective_value);
            self.callback_.data_out.external_solution_query_origin = external_solution_query_origin;
            // C++ `clearHighsCallbackInput()` uses the owning Highs' numCol, i.e. the
            // original model's column count, not the presolved model's.
            let num_col = self.orig_model_.num_col_ as usize;
            self.callback_.clear_highs_callback_input(num_col);

            let interrupt = self
                .callback_
                .callback_action(kCallbackMipUserSolution, "MIP User solution");
            debug_assert!(!interrupt);
            if self.callback_.data_in.user_has_solution {
                // Objective is assumed to be original_offset +
                // (original_c)^T(original_x), but MIP solver bounds are based on the
                // reduced objective (reduced_c)^T(reduced_x)
                //
                // Now, original_sense*[reduced_offset + (reduced_c)^T(reduced_x)] is an
                // objective in the original space, so
                //
                // f0 + c0^Tx0 = s*(f1 + c1^Tx1)
                //
                // where 0 => original; 1 => reduced
                //
                // This allows the reduced objective value to be deduced as
                //
                // c1^Tx1 = s*(f0 + c0^Tx0) - f1
                //
                // (reduced_c)^T(reduced_x) = original_sense*[original_offset +
                // (original_c)^T(original_x) - reduced_offset]
                let mut bound_violation = 0.0;
                let mut row_violation = 0.0;
                let mut integrality_violation = 0.0;
                let mut user_solution_quad_objective_value = HighsCDouble::from(0.0);
                let feasible = self.solution_feasible(
                    Some(&self.orig_model_),
                    &self.callback_.data_in.user_solution,
                    None,
                    &mut bound_violation,
                    &mut row_violation,
                    &mut integrality_violation,
                    &mut user_solution_quad_objective_value,
                );
                let user_solution_objective_value = f64::from(user_solution_quad_objective_value);
                if !feasible {
                    // count the rejected user solution (mip/hooks.rs)
                    self.callback_
                        .hook_stats
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .user_solutions_rejected += 1;
                    highs_log_user(
                        &self.options_mip_.log_options,
                        HighsLogType::kWarning,
                        &format!(
                            "User-supplied solution has with objective {} has violations: bound = {}; integrality = {}; row = {}\n",
                            c_fmt_g(user_solution_objective_value, 6),
                            c_fmt_g(bound_violation, 4),
                            c_fmt_g(integrality_violation, 4),
                            c_fmt_g(row_violation, 4)
                        ),
                    );
                    return;
                }
                self.callback_
                    .hook_stats
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .user_solutions_checked_feasible += 1;
                let reduced_user_solution = self
                    .mipdata_
                    .post_solve_stack
                    .get_reduced_primal_solution(&self.callback_.data_in.user_solution);
                // HiGHS 1.16 #3492 (h16_user_sol_obj): HighsMipSolverData.cpp lines 2866-2869 (1.16.0)
                let h16_user_sol_obj = self.options_mip_.h16_user_sol_obj;
                let reduced_user_solution_objective_value = if h16_user_sol_obj {
                    // C++ `double(HighsInt(orig_model_->sense_) * quad - model_->offset_)`:
                    // the int-to-double factor is the `double` overload of the friend
                    // `operator*`, which is `quad * sense` (HighsCDouble.h lines 185-187).
                    f64::from(
                        user_solution_quad_objective_value
                            * ((self.orig_model_.sense_ as i32) as f64)
                            - mip_model!(self).offset_,
                    )
                } else {
                    user_solution_objective_value
                };
                let print_display_line = true;
                let is_user_solution = true;
                self.add_incumbent(
                    &reduced_user_solution,
                    reduced_user_solution_objective_value,
                    kSolutionSourceUserSolution,
                    print_display_line,
                    is_user_solution,
                );
            }
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::terminatorConcurrency (lines 2819-2821)
    /// C++: `HighsInt terminatorConcurrency() const`
    pub fn terminator_concurrency(&self) -> HighsInt {
        self.terminator_.num_instance
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::terminatorMyInstance (lines 2823-2825)
    /// C++: `HighsInt terminatorMyInstance() const`
    pub fn terminator_my_instance(&self) -> HighsInt {
        self.terminator_.my_instance
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::terminatorTerminate (lines 2827-2830)
    /// C++: `void terminatorTerminate()`
    pub fn terminator_terminate(&mut self) {
        debug_assert!(self.terminator_active());
        self.terminator_.terminate();
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::terminatorTerminated (lines 2832-2836)
    /// C++: `bool terminatorTerminated() const`
    ///
    /// PORT NOTE (architecture): C++ `const`, but when the terminator is active it assigns
    /// `mipsolver.termination_status_`. In HiGHS 1.15.1 the terminator is never active
    /// (`initialiseTerminator` is only called with the defaults, or copied from a parent
    /// that has the defaults: HighsMipSolver.cpp:50, HighsPrimalHeuristics.cpp:145), so
    /// the write is unreachable and the receiver stays `&self` (parallel callers hold `&ms`).
    pub fn terminator_terminated(&self) -> bool {
        assert!(
            !self.terminator_active(),
            "terminator active: unreachable in HiGHS 1.15.1 (see PORT NOTE)"
        );
        self.termination_status_ != HighsModelStatus::kNotset
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::terminatorReport (lines 2838-2841)
    /// C++: `void terminatorReport() const`
    pub fn terminator_report(&self) {
        if self.terminator_active() {
            self.terminator_
                .report(self.options_mip_.log_options.clone());
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::updatePrimalDualIntegral (lines 2862-2975)
    /// C++: `void updatePrimalDualIntegral(const double from_lower_bound, const double to_lower_bound, const double from_upper_bound, const double to_upper_bound, const bool check_bound_change = true, const bool check_prev_data = true)`
    /// C++ default arguments: `check_bound_change = true`, `check_prev_data = true`
    ///
    /// The C++ `assert`s are `debug_assert!` (active in debug builds, as in C++ without
    /// NDEBUG). The C++ difference variables that only fed the asserts (`lower_bound_difference`,
    /// `upper_bound_difference`) are pure computations and are not kept.
    pub fn update_primal_dual_integral(
        &mut self,
        from_lower_bound: f64,
        to_lower_bound: f64,
        from_upper_bound: f64,
        to_upper_bound: f64,
        check_bound_change: bool,
        check_prev_data: bool,
    ) {
        // Parameters to updatePrimalDualIntegral are lower and upper bounds
        // before/after a change
        //
        // updatePrimalDualIntegral should only be called when there is a
        // change in one of the bounds, except when the final update is
        // made, in which case the bounds must NOT have changed. By default,
        // a check for some bound change is made, unless check_bound_change
        // is false, in which case there is a check for unchanged bounds.
        //
        // HighsPrimaDualIntegral struct contains the following data
        //
        // * value: Current value of the P-D integral
        //
        // * prev_lb: Value of lb that was computed from to_lower_bound in
        //   the previous call. Used as a check that the value of lb
        //   computed from from_lower_bound in this call is equal - to
        //   within bound_change_tolerance. If not true, then a change in lb
        //   has been missed. Only for checking/debugging
        //
        // * prev_ub: Ditto for upper_bound. Only for checking/debugging
        //
        // * prev_gap: Ditto for gap. Only for checking/debugging
        //
        // * prev_time: Used to determine the time spent at the previous gap

        let mut from_lb = 0.0;
        let mut from_ub = 0.0;
        let from_gap = self.limits_to_gap(
            from_lower_bound,
            from_upper_bound,
            &mut from_lb,
            &mut from_ub,
        );
        let mut to_lb = 0.0;
        let mut to_ub = 0.0;
        let to_gap = self.limits_to_gap(to_lower_bound, to_upper_bound, &mut to_lb, &mut to_ub);

        let lb_difference = poss_inf_rel_diff(from_lb, to_lb, to_lb);
        let ub_difference = poss_inf_rel_diff(from_ub, to_ub, to_ub);
        let bound_change_tolerance = 0.0;
        let bound_change =
            lb_difference > bound_change_tolerance || ub_difference > bound_change_tolerance;

        if check_bound_change {
            if !bound_change
                && from_lower_bound == to_lower_bound
                && from_upper_bound == to_upper_bound
            {
                debug_assert!(bound_change);
            }
        } else if bound_change
            && (from_lower_bound != to_lower_bound || from_upper_bound != to_upper_bound)
        {
            debug_assert!(!bound_change);
        }

        // HiGHS 1.16 #3521 (h16_mip_workers): HighsMipSolverData.cpp lines 3024-3062 (1.16.0)
        if self.options_mip_.h16_mip_workers {
            if self.mipdata_.primal_dual_integral.value <= -kHighsInf {
                // First call, so zero the P-D integral and no housekeeping to
                // test
                self.mipdata_.primal_dual_integral.value = 0.0;
                self.mipdata_.primal_dual_integral.prev_gap = 1.0;
            } else if check_prev_data {
                // updatePrimalDualIntegral has been called previously, so can
                // test housekeeping, except after restart
                let pdi = &self.mipdata_.primal_dual_integral;
                let lb_inconsistency = poss_inf_rel_diff(from_lb, pdi.prev_lb, pdi.prev_lb);
                let lb_consistent = lb_inconsistency < 1e-12;
                let ub_inconsistency = poss_inf_rel_diff(from_ub, pdi.prev_ub, pdi.prev_ub);
                let ub_consistent = ub_inconsistency < 1e-12;
                let gap_inconsistency = poss_inf_rel_diff(from_gap, pdi.prev_gap, 1.0);
                let gap_consistent = gap_inconsistency < 1e-12;
                debug_assert!(lb_consistent);
                debug_assert!(ub_consistent);
                debug_assert!(gap_consistent);
            }
            // Update the P-D integral using the gap capped at 1, so an infinite
            // gap (no incumbent, no lower bound, or gap bigger than 1)
            // contributes 1
            let time = self.timer_.read();
            let time_diff = time - self.mipdata_.primal_dual_integral.prev_time;
            debug_assert!(time_diff >= 0.0);
            let prev_gap = self.mipdata_.primal_dual_integral.prev_gap;
            self.mipdata_.primal_dual_integral.value += time_diff * cpp_min(prev_gap, 1.0);
            self.mipdata_.primal_dual_integral.prev_time = time;
            self.mipdata_.primal_dual_integral.prev_lb = to_lb;
            self.mipdata_.primal_dual_integral.prev_ub = to_ub;
            self.mipdata_.primal_dual_integral.prev_gap = to_gap;
            return;
        }
        if self.mipdata_.primal_dual_integral.value > -kHighsInf {
            // updatePrimalDualIntegral has been called previously, so can
            // usually test housekeeping, even if gap is still inf
            //
            // The one case where the checking can't be done comes after restart, where
            // the
            //
            if check_prev_data {
                // These housekeeping tests check that the previous saved
                // lower/upper bounds and gap are very close to the "from"
                // lower/upper bounds and corresponding gap. They are usually
                // identical, but rounding error can occur when passing through
                // reset, when the old/new offsets are added/subtracted from the
                // bounds due to changes in offset during presolve.
                let pdi = &self.mipdata_.primal_dual_integral;
                let lb_inconsistency = poss_inf_rel_diff(from_lb, pdi.prev_lb, pdi.prev_lb);
                let lb_consistent = lb_inconsistency < 1e-12;
                let ub_inconsistency = poss_inf_rel_diff(from_ub, pdi.prev_ub, pdi.prev_ub);
                let ub_consistent = ub_inconsistency < 1e-12;
                let gap_inconsistency = poss_inf_rel_diff(from_gap, pdi.prev_gap, 1.0);
                let gap_consistent = gap_inconsistency < 1e-12;
                debug_assert!(lb_consistent);
                debug_assert!(ub_consistent);
                debug_assert!(gap_consistent);
            }
            if to_gap < kHighsInf {
                let time = self.timer_.read();
                if from_gap < kHighsInf {
                    // Need to update the P-D integral
                    let time_diff = time - self.mipdata_.primal_dual_integral.prev_time;
                    debug_assert!(time_diff >= 0.0);
                    let prev_gap = self.mipdata_.primal_dual_integral.prev_gap;
                    self.mipdata_.primal_dual_integral.value += time_diff * prev_gap;
                }
                self.mipdata_.primal_dual_integral.prev_time = time;
            }
        } else {
            self.mipdata_.primal_dual_integral.value = 0.0;
        }
        self.mipdata_.primal_dual_integral.prev_lb = to_lb;
        self.mipdata_.primal_dual_integral.prev_ub = to_ub;
        self.mipdata_.primal_dual_integral.prev_gap = to_gap;
    }
}

/// HiGHS: HighsMipSolverData.cpp::clockOff (lines 1967-1985)
/// C++: `static void clockOff(HighsProfiling * profiling)`
///
/// The body starts `if (!profiling->mip_) return;`: profiling is a stub (mip_timer.rs), so
/// `mip_` is false and this is a no-op.
pub fn clock_off(_profiling: &mut HighsProfiling) {}

/// HiGHS: HighsMipSolverData.cpp::possInfRelDiff (lines 2843-2860)
/// C++: `static double possInfRelDiff(const double v0, const double v1, const double den)`
pub fn poss_inf_rel_diff(v0: f64, v1: f64, den: f64) -> f64 {
    let rel_diff: f64;
    if v0.abs() == kHighsInf {
        if v1.abs() == kHighsInf {
            rel_diff = 0.0;
        } else {
            rel_diff = kHighsInf;
        }
    } else if v1.abs() == kHighsInf {
        rel_diff = kHighsInf;
    } else {
        rel_diff = (v1 - v0).abs() / 1.0f64.max(den.abs());
    }
    rel_diff
}

impl HighsPrimaDualIntegral {
    /// HiGHS: HighsMipSolverData.cpp::HighsPrimaDualIntegral::initialise (lines 2977-2977)
    /// C++: `void initialise()`
    pub fn initialise(&mut self) {
        self.value = -kHighsInf;
    }
}

impl HighsTerminator {
    /// HiGHS: HighsMipSolverData.cpp::HighsTerminator::clear (lines 2979-2983)
    /// C++: `void clear()`
    pub fn clear(&mut self) {
        self.num_instance = 0;
        self.my_instance = kNoThreadInstance;
        self.record = Vec::new();
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsTerminator::initialise (lines 2985-2991)
    /// C++: `void initialise(HighsInt num_instance_, HighsInt my_instance_, HighsModelStatus * record_)`
    ///
    /// PORT NOTE (architecture): `record_` (a pointer into an array shared by racing
    /// instances) is an owned copy; `None` = nullptr. The terminator is never active in
    /// 1.15.1 (see `terminator_terminated`), so nothing writes a shared record.
    pub fn initialise(
        &mut self,
        num_instance_: HighsInt,
        my_instance_: HighsInt,
        record_: Option<Vec<HighsModelStatus>>,
    ) {
        self.clear();
        self.num_instance = num_instance_;
        self.my_instance = my_instance_;
        self.record = record_.unwrap_or_default();
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsTerminator::concurrency (lines 2993-2993)
    /// C++: `HighsInt concurrency() const`
    pub fn concurrency(&self) -> HighsInt {
        self.num_instance
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsTerminator::terminate (lines 2995-2999)
    /// C++: `void terminate()`
    pub fn terminate(&mut self) {
        debug_assert!(!self.record.is_empty());
        debug_assert!(self.my_instance < self.num_instance);
        self.record[self.my_instance as usize] = HighsModelStatus::kHighsInterrupt;
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsTerminator::terminationStatus (lines 3001-3008)
    /// C++: `HighsModelStatus terminationStatus() const`
    pub fn termination_status(&self) -> HighsModelStatus {
        debug_assert!(!self.record.is_empty());
        for instance in 0..self.num_instance {
            if self.record[instance as usize] != HighsModelStatus::kNotset {
                return self.record[instance as usize];
            }
        }
        HighsModelStatus::kNotset
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsTerminator::report (lines 3010-3016)
    /// C++: `void report(const HighsLogOptions log_options) const`
    pub fn report(&self, log_options: HighsLogOptions) {
        highs_log_user(&log_options, HighsLogType::kInfo, "\nTerminator:        ");
        for instance in 0..self.num_instance {
            highs_log_user(
                &log_options,
                HighsLogType::kInfo,
                &format!(" {:20}", self.record[instance as usize] as i32),
            );
        }
        highs_log_user(&log_options, HighsLogType::kInfo, "\n");
    }
}

/// C++ `highsSparseTranspose` (util/HighsUtils.cpp:60-86): the row-wise copy of a
/// column-wise matrix, returned as `(ARstart, ARindex, ARvalue)`.
#[allow(clippy::needless_range_loop)]
fn highs_sparse_transpose(
    num_row: HighsInt,
    num_col: HighsInt,
    astart: &[HighsInt],
    aindex: &[HighsInt],
    avalue: &[f64],
) -> (Vec<HighsInt>, Vec<HighsInt>, Vec<f64>) {
    // Make a AR copy
    let nrow = num_row as usize;
    let mut iwork: Vec<HighsInt> = vec![0; nrow];
    let mut arstart: Vec<HighsInt> = vec![0; nrow + 1];
    let acount_x = aindex.len();
    let mut arindex: Vec<HighsInt> = vec![0; acount_x];
    let mut arvalue: Vec<f64> = vec![0.0; acount_x];
    for k in 0..acount_x {
        debug_assert!(aindex[k] < num_row);
        iwork[aindex[k] as usize] += 1;
    }
    for i in 1..=nrow {
        arstart[i] = arstart[i - 1] + iwork[i - 1];
    }
    iwork[..nrow].copy_from_slice(&arstart[..nrow]);
    for icol in 0..num_col as usize {
        for k in astart[icol] as usize..astart[icol + 1] as usize {
            let irow = aindex[k] as usize;
            let iput = iwork[irow] as usize;
            iwork[irow] += 1;
            arindex[iput] = icol as HighsInt;
            arvalue[iput] = avalue[k];
        }
    }
    (arstart, arindex, arvalue)
}

impl HighsMipSolver {
    /// rusolve racing (Shinano et al. 2018): offer the incumbent just stored to the other racers.
    pub fn race_publish(&self, solobj: f64) {
        if self.submip {
            return;
        }
        if let Some(race) = &self.race {
            let fp = crate::mip::race::model_fingerprint(&self.mipdata_.presolved_model);
            race.publish(fp, solobj, &self.mipdata_.incumbent);
        }
    }

    /// rusolve: install a better incumbent another racer found, between two nodes.
    pub fn race_import(&mut self) {
        if self.submip {
            return;
        }
        let Some(race) = self.race.clone() else {
            return;
        };
        let version = race.version.load(std::sync::atomic::Ordering::Acquire);
        if version == self.race_seen {
            return;
        }
        self.race_seen = version;
        let fp = crate::mip::race::model_fingerprint(&self.mipdata_.presolved_model);
        if let Some((obj, sol)) = race.fetch(fp, self.mipdata_.upper_bound) {
            self.add_incumbent(&sol, obj, kSolutionSourceUserSolution, false, false);
        }
    }
}

#[cfg(test)]
mod h16_user_sol_obj_tests {
    //! HiGHS 1.16 #3492 (h16_user_sol_obj) tests. A maximisation knapsack with a nonzero
    //! model offset; a user callback offers the optimal solution at the first query
    //! and records the primal bound seen at every query.
    use super::*;
    use crate::lp_data::highs_run::Highs;
    use std::sync::{Arc, Mutex};

    const N: usize = 12;
    const OFFSET: f64 = 10.0;

    /// Values and weights of the knapsack items.
    fn items() -> (Vec<f64>, Vec<f64>, f64) {
        let values: Vec<f64> = (0..N).map(|i| (10 + (i * 37) % 23) as f64).collect();
        let weights: Vec<f64> = (0..N).map(|i| (5 + (i * 53) % 19) as f64).collect();
        let capacity = 60.0;
        (values, weights, capacity)
    }

    /// Brute-force optimum of the knapsack (without the offset) and its 0/1 solution.
    fn brute_force() -> (f64, Vec<f64>) {
        let (values, weights, capacity) = items();
        let mut best = (-1.0, vec![0.0; N]);
        for mask in 0u32..(1 << N) {
            let (mut value, mut weight) = (0.0, 0.0);
            let mut x = vec![0.0; N];
            for i in 0..N {
                if mask & (1 << i) != 0 {
                    value += values[i];
                    weight += weights[i];
                    x[i] = 1.0;
                }
            }
            if weight <= capacity && value > best.0 {
                best = (value, x);
            }
        }
        best
    }

    fn knapsack_model() -> HighsLp {
        let (values, weights, capacity) = items();
        let mut lp = HighsLp::default();
        lp.num_col_ = N as HighsInt;
        lp.num_row_ = 1;
        lp.sense_ = ObjSense::kMaximize;
        lp.offset_ = OFFSET;
        lp.col_cost_ = values;
        lp.col_lower_ = vec![0.0; N];
        lp.col_upper_ = vec![1.0; N];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![capacity];
        lp.integrality_ = vec![HighsVarType::kInteger; N];
        lp.a_matrix_.num_col_ = N as HighsInt;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = (0..=N as HighsInt).collect();
        lp.a_matrix_.index_ = vec![0; N];
        lp.a_matrix_.value_ = weights;
        lp
    }

    /// Runs the knapsack with a user solution offered at the first query. Returns the
    /// solver objective and the primal bound recorded at each query.
    fn run_with_offer(user_sol_obj: bool) -> (f64, Vec<f64>) {
        let (_, offer) = brute_force();
        let mut highs = Highs::default();
        highs.set_option_value_bool("output_flag", false);
        highs.set_option_value_double("mip_heuristic_effort", 0.0);
        highs.set_option_value_bool("mip_heuristic_run_feasibility_jump", false);
        highs.set_option_value_bool("h16_user_sol_obj", user_sol_obj);
        highs.pass_model(knapsack_model());
        let mut start = HighsSolution::default();
        start.col_value = vec![0.0; N];
        start.col_value[0] = 1.0;
        highs.set_solution(&start);

        let bounds = Arc::new(Mutex::new(Vec::new()));
        let recorder = Arc::clone(&bounds);
        let offer_for_callback = offer.clone();
        highs.set_callback(HighsCallbackFunctionType::new(
            move |_, _, data_out, data_in| {
                let mut recorded = recorder.lock().unwrap();
                if recorded.is_empty() {
                    data_in.user_has_solution = true;
                    data_in.user_solution = offer_for_callback.clone();
                }
                recorded.push(data_out.mip_primal_bound);
            },
        ));
        highs.start_callback(kCallbackMipUserSolution);
        assert_eq!(highs.run(), HighsStatus::kOk);
        let objective = highs.get_info().objective_function_value;
        let recorded = bounds.lock().unwrap().clone();
        (objective, recorded)
    }

    #[test]
    fn h16_user_sol_obj_on_offer_reaches_incumbent() {
        let (best, _) = brute_force();
        let optimum = OFFSET + best;
        let (objective, bounds) = run_with_offer(true);
        assert_eq!(objective, optimum);
        assert!(bounds.len() >= 2, "callback queried {} times", bounds.len());
        for &bound in &bounds[1..] {
            assert_eq!(bound, optimum);
        }
    }

    #[test]
    fn h16_user_sol_obj_off_keeps_1_15_1_objective() {
        let (best, _) = brute_force();
        let (objective, bounds) = run_with_offer(false);
        // The final objective is still the optimum (branch and bound finds it).
        assert_eq!(objective, OFFSET + best);
        // The start solution (item 0 only) gives the first primal bound.
        let (values, _, _) = items();
        assert_eq!(bounds[0], OFFSET + values[0]);
        // 1.15.1 passes the offered solution's plain objective to addIncumbent, and
        // the offered optimum is not accepted: every later query still sees the start
        // bound.
        assert!(bounds.len() >= 2, "callback queried {} times", bounds.len());
        for &bound in &bounds[1..] {
            assert_eq!(bound, bounds[0]);
        }
    }
}
