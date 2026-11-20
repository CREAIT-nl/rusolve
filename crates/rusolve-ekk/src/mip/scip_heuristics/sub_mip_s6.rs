//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! sub-MIP copy: solve_sub_mip with the SCIP crossover cutoff (heur_crossover.c:700-720).

use crate::prelude::*;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// Extras of `s6_solve_sub_mip_ext`; the default is the plain sub-MIP.
#[derive(Clone, Copy, Default)]
pub(crate) struct S6SubMipExt {
    /// the sub-LP carries its own objective: `cutoff` is its objective bound as given, not
    /// capped by the parent's incumbent limit
    pub own_objective: bool,
    /// stop the sub-MIP after this many improving solutions
    pub max_improving_sols: Option<HighsInt>,
    /// stop the sub-MIP after this much work (the LP-iteration limit of SCIP's proximity
    /// sub-MIP, in work units)
    pub work_limit: Option<u64>,
    /// SCIP's sub-SCIP settings for proximity (`heur_proximity.c` `SCIPsetCommonSubscipParams`
    /// block): no cut separation, no restarts, no recursive heuristics, light heuristic effort
    pub cheap: bool,
}

impl HighsPrimalHeuristics {
    /// HiGHS: HighsPrimalHeuristics.cpp::solveSubMip (lines 77-247)
    /// SCIP 11.0.0 (32dd425): heur_crossover.c (lines 700-720)
    // Copy of primal_heuristics/part1.rs solve_sub_mip; the only
    // semantic change is objective_bound = min(upper_limit, cutoff), plus the returned
    // sub-MIP node count.
    // PORT NOTE: the profiling calls (`profiling_->start/stop/solveCall/setSubMip`, the
    // `setProfiling` hand-over, and the `sub_solver_` printf) are dropped: `HighsProfiling` is a
    // no-op stub and `HighsMipSolver::set_profiling` is not yet written.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn s6_solve_sub_mip(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        lp: HighsLp,
        basis: HighsBasis,
        fixing_rate: f64,
        col_lower: Vec<f64>,
        col_upper: Vec<f64>,
        maxleaves: HighsInt,
        maxnodes: HighsInt,
        stallnodes: HighsInt,
        cutoff: f64,
    ) -> (bool, i64) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_SCIP_SUBMIP);
        let r = self.s6_solve_sub_mip_src(
            ms,
            worker,
            lp,
            basis,
            fixing_rate,
            col_lower,
            col_upper,
            maxleaves,
            maxnodes,
            stallnodes,
            cutoff,
            kSolutionSourceSubMip,
        );
        (r.0, r.1)
    }

    /// `s6_solve_sub_mip` with the solution source of the point handed to `try_solution` and the
    /// sub-MIP's final model status (rusolve mip_alns: the bandit needs infeasible / optimal /
    /// limit to adapt the fixing rate, and its own source letter in the log).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn s6_solve_sub_mip_src(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        lp: HighsLp,
        basis: HighsBasis,
        fixing_rate: f64,
        col_lower: Vec<f64>,
        col_upper: Vec<f64>,
        maxleaves: HighsInt,
        maxnodes: HighsInt,
        stallnodes: HighsInt,
        cutoff: f64,
        source: i32,
    ) -> (bool, i64, HighsModelStatus) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_SCIP_SUBMIP);
        self.s6_solve_sub_mip_ext(
            ms,
            worker,
            lp,
            basis,
            fixing_rate,
            col_lower,
            col_upper,
            maxleaves,
            maxnodes,
            stallnodes,
            cutoff,
            source,
            S6SubMipExt::default(),
        )
    }

    /// `s6_solve_sub_mip_src` with the two extras of the Fischetti heuristics (proximity search,
    /// local branching): a sub-MIP whose objective is not the
    /// parent's (no `objective_bound` from the parent's incumbent) and an improving-solution limit.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn s6_solve_sub_mip_ext(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        lp: HighsLp,
        basis: HighsBasis,
        fixing_rate: f64,
        col_lower: Vec<f64>,
        col_upper: Vec<f64>,
        maxleaves: HighsInt,
        maxnodes: HighsInt,
        stallnodes: HighsInt,
        cutoff: f64,
        source: i32,
        ext: S6SubMipExt,
    ) -> (bool, i64, HighsModelStatus) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_SCIP_SUBMIP);
        let mut submipoptions: HighsOptions = ms.options_mip_.clone();
        let mut submip: HighsLp = lp;

        // set bounds and restore integrality of the lp relaxation copy
        submip.col_lower_ = col_lower;
        submip.col_upper_ = col_upper;
        submip.integrality_ = mip_model!(ms).integrality_.clone();
        submip.offset_ = 0.0;

        // set limits
        submipoptions.mip_max_leaves = maxleaves;
        submipoptions.output_flag = false;

        // HiGHS quirk: the sub-MIP log is switched on only for one size pair that does not occur
        // for real models. The printf of the profiling branch is dropped (profiling is a stub).
        let allow_submip_log = true;
        if allow_submip_log && submip.num_col_ == -54 && submip.num_row_ == -172 {
            submipoptions.output_flag = true;
        }
        // PORT NOTE: C++ `log_options.output_flag` points at `submipoptions.output_flag`
        // (the HighsOptions copy constructor re-targets it), so the sub-MIP log is silent;
        // the port's copies must be re-synced.
        submipoptions.sync_log_options();

        submipoptions.mip_max_nodes = maxnodes;
        submipoptions.mip_max_stall_nodes = stallnodes;
        submipoptions.mip_pscost_minreliable = 0;
        submipoptions.time_limit -= ms.timer_.read();
        // objective_bound = min(upper_limit, cutoff) (heur_crossover.c:700-720 cutoff use)
        submipoptions.objective_bound = if ext.own_objective {
            cutoff
        } else {
            cpp_min(ms.mipdata_.workers[worker].upper_limit, cutoff)
        };
        if let Some(n) = ext.max_improving_sols {
            submipoptions.mip_max_improving_sols = n;
        }

        if !ms.submip {
            let mut curr_abs_gap =
                ms.mipdata_.workers[worker].upper_limit - ms.mipdata_.lower_bound;

            if curr_abs_gap == kHighsInf {
                curr_abs_gap = ms.mipdata_.lower_bound.abs();
                if curr_abs_gap == kHighsInf {
                    curr_abs_gap = 0.0;
                }
            }

            submipoptions.mip_rel_gap = 0.0;
            submipoptions.mip_abs_gap = ms.mipdata_.feastol * cpp_max(curr_abs_gap, 1000.0);
        }

        // check if only root presolve is allowed
        if submipoptions.mip_root_presolve_only {
            submipoptions.presolve = kHighsOffString.to_string();
        } else {
            submipoptions.presolve = kHighsOnString.to_string();
        }
        submipoptions.mip_detect_symmetry = false;
        // rusolve mip_alns: the lever lives in the parent only
        submipoptions.mip_alns = false;
        submipoptions.mip_proximity = false;
        submipoptions.mip_local_branching = false;
        submipoptions.mip_heuristic_effort = 0.8;
        if ext.cheap {
            // heur_proximity.c: separation off, no recursive LNS heuristics (the sub-SCIP's
            // `SCIPsetSeparating(OFF)` and `SCIPsetSubscipsOff`); no restarts inside the node-limited
            // sub-MIP
            submipoptions.mip_allow_cut_separation_at_nodes = false;
            submipoptions.mip_allow_restart = false;
            submipoptions.mip_heuristic_run_rins = false;
            submipoptions.mip_heuristic_run_rens = false;
        }
        // setup solver and run it

        let solution = HighsSolution {
            value_valid: false,
            dual_valid: false,
            ..HighsSolution::default()
        };

        // Create HighsMipSolver instance for sub-MIP
        //
        // C++ hands the sub-MIP a pointer to the parent's callback (`*mipsolver.callback_`).
        // The port moves the callback in and back out. Under parallel workers several
        // sub-MIPs run at once, so moving the parent's field would race: each worker gets
        // its own clone instead. The sub-MIP never invokes the user callback (every call
        // site is gated on `!submip`), and the clone shares the user closure (`Arc`), so
        // only the sub-MIP's scratch `data_out` writes stay local to the clone.
        let parallel_submip = ms.parallel_lock_active();
        let mut parallel_callback = parallel_submip.then(|| ms.callback_.clone());
        let mut submipsolver = HighsMipSolver::new(
            match parallel_callback.as_mut() {
                Some(cb) => cb,
                None => &mut ms.callback_,
            },
            &submipoptions,
            &submip,
            &solution,
            true,
            ms.submip_level + 1,
        );
        // Initialise termination_status_ and propagate any terminator to the sub-MIP
        submipsolver.initialise_terminator_2(ms);
        // before the root LP (mip_improve_early) the parent has no basis to hand over
        submipsolver.rootbasis = basis.valid.then(|| Box::new(basis));
        submipsolver.work_limit = ext.work_limit;
        let pscost_k = ms.mipdata_.workers[worker].get_pseudocost();
        submipsolver.pscostinit = Some(Box::new(detach!(ms.mipdata_.pseudocosts[pscost_k], |p| {
            HighsPseudocostInitialization::new(p, 1)
        })));
        // C++ `const` pointers to the parent's clique table and implications: the sub-MIP reads
        // them in place (they stay in the parent's slots, pinned under parallel workers).
        submipsolver.clqtableinit = Some(ParentPtr::new(&*ms.mipdata_.cliquetable));
        submipsolver.implicinit = Some(ParentPtr::new(&*ms.mipdata_.implications));
        // Solve the sub-MIP
        submipsolver.run();

        // The pointers must not outlive this call (C++ `submipsolver` is a local).
        submipsolver.clqtableinit = None;
        submipsolver.implicinit = None;
        if !parallel_submip {
            ms.callback_ = std::mem::take(&mut submipsolver.callback_);
        }

        let sub_node_count: i64 = submipsolver.node_count_;
        let sub_status: HighsModelStatus = submipsolver.modelstatus_;

        ms.mipdata_.workers[worker]
            .update_heur_stats_max_sub_mip_level(submipsolver.max_submip_level + 1);

        // C++ `submipsolver` is destroyed at the end of the function, so it is released once
        // the result handling below has run, on every path.
        let result = 'sub: {
            // 22/07/25: Seems impossible for submipsolver.mipdata_ to be a null pointer after
            // calling HighsMipSolver::run(). The C++ null check and assert are not carried over:
            // `mipdata_` is a Box in Rust.
            if submipsolver.termination_status_ != HighsModelStatus::kNotset {
                ms.mipdata_.workers[worker]
                    .set_heur_termination_status(submipsolver.termination_status_);
                break 'sub false;
            }

            let num_unfixed =
                (ms.mipdata_.integral_cols.len() + ms.mipdata_.continuous_cols.len()) as f64;
            let adjustmentfactor = ((1.0 - fixing_rate) * ms.mipdata_.integral_cols.len() as f64
                + ms.mipdata_.continuous_cols.len() as f64)
                / cpp_max(1.0, num_unfixed);
            // (double)mipsolver.orig_model_->a_matrix_.value_.size();
            let adjusted_lp_iterations: i64 = (adjustmentfactor
                * submipsolver.mipdata_.total_lp_iterations as f64)
                as usize as i64;
            ms.mipdata_.workers[worker].update_heur_stats_lp_iters(
                adjusted_lp_iterations,
                submipsolver.mipdata_.total_repair_lp,
                submipsolver.mipdata_.total_repair_lp_feasible,
                submipsolver.mipdata_.total_repair_lp_iterations,
            );
            // Warning: This will not be deterministic if sub-mips are run in parallel
            if ms.submip {
                let sub_nodes = (adjustmentfactor * submipsolver.node_count_ as f64) as i64;
                ms.mipdata_.num_nodes += cpp_max(1, sub_nodes);
            }

            if submipsolver.modelstatus_ == HighsModelStatus::kInfeasible {
                ms.mipdata_.workers[worker].update_heur_stats_infeas_observations(fixing_rate);
            }
            if submipsolver.node_count_ <= 1
                && submipsolver.modelstatus_ == HighsModelStatus::kInfeasible
            {
                break 'sub false;
            }
            let old_upper_limit = ms.mipdata_.workers[worker].upper_limit;
            if submipsolver.modelstatus_ != HighsModelStatus::kInfeasible
                && !submipsolver.solution_.is_empty()
            {
                self.try_solution(ms, &submipsolver.solution_, source, worker);
            }

            if ms.mipdata_.workers[worker].upper_limit < old_upper_limit {
                // remember fixing rate as good
                ms.mipdata_.workers[worker].update_heur_stats_success_observations(fixing_rate);
            }

            true
        };
        submipsolver.release();
        (result, sub_node_count, sub_status)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

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

    /// Solve on a fresh thread with one scheduler worker.
    fn solve_on_fresh_thread(items: Items) -> HighsMipSolver {
        std::thread::spawn(move || {
            initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            let lp = knapsack_lp(items);
            let mut ms = HighsMipSolver::new(
                &mut HighsCallback::default(),
                &options,
                &lp,
                &HighsSolution::default(),
                false,
                0,
            );
            ms.run();
            HighsTaskExecutor::shutdown(true);
            ms
        })
        .join()
        .unwrap()
    }

    #[test]
    fn s6_submip_cutoff_inf_matches_solve_sub_mip() {
        let mut ms1 = solve_on_fresh_thread(ITEMS);
        let mut ms2 = solve_on_fresh_thread(ITEMS);

        let lp1 = mip_model!(ms1).clone();
        let basis1 = ms1.mipdata_.firstrootbasis.clone();
        let cl1 = lp1.col_lower_.clone();
        let cu1 = lp1.col_upper_.clone();

        let lp2 = mip_model!(ms2).clone();
        let basis2 = ms2.mipdata_.firstrootbasis.clone();
        let cl2 = lp2.col_lower_.clone();
        let cu2 = lp2.col_upper_.clone();

        let r1 = detach!(ms1.mipdata_.heuristics, |h| h.solve_sub_mip(
            &mut ms1, 0, lp1, basis1, 0.0, cl1, cu1, kHighsIInf, kHighsIInf, kHighsIInf,
        ));

        let (r2, nodes2) = detach!(ms2.mipdata_.heuristics, |h| h.s6_solve_sub_mip(
            &mut ms2,
            0,
            lp2,
            basis2,
            0.0,
            cl2,
            cu2,
            kHighsIInf,
            kHighsIInf,
            kHighsIInf,
            f64::INFINITY,
        ));

        assert_eq!(r1, r2);
        assert!(nodes2 >= 0);
        assert_eq!(ms1.mipdata_.num_nodes, ms2.mipdata_.num_nodes);
        assert_eq!(ms1.mipdata_.incumbent.len(), ms2.mipdata_.incumbent.len());
        for (a, b) in ms1.mipdata_.incumbent.iter().zip(&ms2.mipdata_.incumbent) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(
            ms1.mipdata_.upper_bound.to_bits(),
            ms2.mipdata_.upper_bound.to_bits()
        );
        assert_eq!(
            ms1.solution_objective_.to_bits(),
            ms2.solution_objective_.to_bits()
        );
        assert_eq!(
            ms1.mipdata_.total_lp_iterations,
            ms2.mipdata_.total_lp_iterations
        );
        assert_eq!(
            ms1.mipdata_.workers[0].upper_limit.to_bits(),
            ms2.mipdata_.workers[0].upper_limit.to_bits()
        );
    }

    #[test]
    fn s6_submip_cutoff_below_optimum_adds_no_incumbent() {
        let mut ms = solve_on_fresh_thread(ITEMS);

        let optimum = ms.solution_objective_;
        let upper_before = ms.mipdata_.workers[0].upper_limit;
        let obj_before = ms.solution_objective_;
        let ub_before = ms.mipdata_.upper_bound;
        let inc_before = ms.mipdata_.incumbent.clone();

        let lp = mip_model!(ms).clone();
        let basis = ms.mipdata_.firstrootbasis.clone();
        let cl = lp.col_lower_.clone();
        let cu = lp.col_upper_.clone();

        let cutoff = optimum - 1.0;

        let (_result, _nodes) = detach!(ms.mipdata_.heuristics, |h| h.s6_solve_sub_mip(
            &mut ms, 0, lp, basis, 0.0, cl, cu, kHighsIInf, kHighsIInf, kHighsIInf, cutoff,
        ));

        assert_eq!(
            ms.mipdata_.workers[0].upper_limit.to_bits(),
            upper_before.to_bits()
        );
        // The result is HiGHS's: the sub-MIP's cleanupSolve reports kOptimal whenever it holds
        // a feasible solution (HighsMipSolver.cpp:1089-1095), which objective_bound does not
        // filter, so solveSubMip returns true here. The parent's trySolution rejects the
        // non-improving solution, so no incumbent appears.
        assert_eq!(ms.solution_objective_.to_bits(), obj_before.to_bits());
        assert_eq!(ms.mipdata_.upper_bound.to_bits(), ub_before.to_bits());
        assert_eq!(ms.mipdata_.incumbent, inc_before);
    }
}
