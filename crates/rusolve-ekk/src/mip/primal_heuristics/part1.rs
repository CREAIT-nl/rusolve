//! HiGHS source: mip/HighsPrimalHeuristics.cpp (lines 1-970)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS latest (#3445, b1282872a): HighsPrimalHeuristics.cpp (line 419, b1282872a)
const K_DEPTH_FIRST_NODE_LIMIT: i64 = 1000;

/// C++ `std::pair<double, uint64_t>::operator<` (libc++): `x.first < y.first ||
/// (!(y.first < x.first) && x.second < y.second)`.
fn pair_lt(x: (f64, u64), y: (f64, u64)) -> bool {
    x.0 < y.0 || (!(y.0 < x.0) && x.1 < y.1)
}

/// C++ `std::tuple<double, uint64_t, HighsInt>::operator<`, element by element, so that a NaN
/// key falls through to the next element as it does in C++.
fn tuple_lt(x: (f64, u64, HighsInt), y: (f64, u64, HighsInt)) -> bool {
    if x.0 < y.0 {
        return true;
    }
    if y.0 < x.0 {
        return false;
    }
    if x.1 < y.1 {
        return true;
    }
    if y.1 < x.1 {
        return false;
    }
    x.2 < y.2
}

/// C++ `HighsHashHelpers::hash((uint64_t(col) << 32) + n)`: the RENS / RINS sort key tie-break,
/// where `n` is `heurlp.getFractionalIntegers().size()`.
fn frac_key_hash(col: HighsInt, n: usize) -> u64 {
    let x = ((col as u64) << 32).wrapping_add(n as u64);
    hash_u64_bytes(x as u32, (x >> 32) as u32)
}

/// C++ `getFixVal` lambda in RENS (HighsPrimalHeuristics.cpp 428-439).
pub(crate) fn rens_get_fix_val(
    ms: &HighsMipSolver,
    localdom: &HighsDomain,
    col: HighsInt,
    fracval: f64,
) -> f64 {
    let ci = col as usize;
    // reinforce direction of this solution away from root
    // solution if the change is at least 0.4
    // otherwise take the direction where the objective gets worse
    // if objective is zero round to nearest integer
    let rootchange = if ms.mipdata_.rootlpsol.is_empty() {
        0.0
    } else {
        fracval - ms.mipdata_.rootlpsol[ci]
    };
    let mut fixval = calc_fix_val(rootchange, fracval, ms.col_cost_2(col));
    // make sure we do not set an infeasible domain
    fixval = cpp_min(localdom.col_upper_[ci], fixval);
    fixval = cpp_max(localdom.col_lower_[ci], fixval);
    fixval
}

/// C++ `getFixVal` lambda in RINS (HighsPrimalHeuristics.cpp 723-737). The lambda captures
/// `fixtolpsol` by reference, so the current value is passed at each call.
fn rins_get_fix_val(
    ms: &HighsMipSolver,
    localdom: &HighsDomain,
    relaxationsol: &[f64],
    fixtolpsol: bool,
    col: HighsInt,
    fracval: f64,
) -> f64 {
    let ci = col as usize;
    let mut fixval = if fixtolpsol {
        // RINS neighbourhood (with extension)
        (relaxationsol[ci] + 0.5).floor()
    } else {
        // reinforce direction of this solution away from root
        // solution if the change is at least 0.4
        // otherwise take the direction where the objective gets worse
        // if objective is zero round to nearest integer
        // HiGHS quirk: no `rootlpsol.empty()` guard here, unlike RENS. C++ indexes an empty
        // vector (undefined behaviour); Rust panics on the index.
        calc_fix_val(
            fracval - ms.mipdata_.rootlpsol[ci],
            fracval,
            ms.col_cost_2(col),
        )
    };
    // make sure we do not set an infeasible domain
    fixval = cpp_min(localdom.col_upper_[ci], fixval);
    fixval = cpp_max(localdom.col_lower_[ci], fixval);
    fixval
}

/// HiGHS 1.16 (#3237): HighsPrimalHeuristics.cpp::subMipShouldRetry (lines 43-66), its limit
/// test `worker.search_ptr_->checkLimits()` (`HighsSearch::checkLimits`, nodeOffset 0): the
/// worker's local limits under the parallel lock, else `mipdata_->checkLimits()`.
fn h16_sub_mip_check_limits(ms: &mut HighsMipSolver, worker: usize) -> bool {
    if ms.parallel_lock_active() {
        let ms_ref: &HighsMipSolver = ms;
        let search = ms_ref.mipdata_.workers[worker]
            .search_ptr_
            .as_deref()
            .expect("subMipShouldRetry: worker search is running (re-entrant call)");
        return search.check_local_limits_w(ms_ref, worker);
    }
    ms.check_limits(0)
}

/// C++ `localdom.conflictAnalysis(worker.getConflictPool(), worker.getGlobalDomain(),
/// worker.getPseudocost())` for a domain that is not in `ms` (the worker's
/// pseudocost is detached for the call).
fn conflict_analysis_worker(ms: &mut HighsMipSolver, localdom: &mut HighsDomain, worker: usize) {
    let _src = crate::work::conflict_src(crate::work::ConflictSrc::Heur);
    let cp = ms.mipdata_.workers[worker].get_conflict_pool();
    let gd = ms.mipdata_.workers[worker].get_global_domain();
    let k = ms.mipdata_.workers[worker].get_pseudocost();
    detach!(ms.mipdata_.pseudocosts[k], |pc| localdom
        .conflict_analysis(ms, cp, gd, pc));
}

impl HighsPrimalHeuristics {
    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::HighsPrimalHeuristics (lines 35-41)
    /// C++: `HighsPrimalHeuristics(HighsMipSolver & mipsolver)`
    pub fn new(mipsolver: &HighsMipSolver) -> Self {
        HighsPrimalHeuristics {
            intcols: Vec::new(),
            success_observations: 0.0,
            num_success_observations: 0,
            infeas_observations: 0.0,
            num_infeas_observations: 0,
            randgen: HighsRandom::new(mipsolver.options_mip_.random_seed as HighsUInt),
            guided_lp_iterations: 0,
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::setupIntCols (lines 43-75)
    /// C++: `void setupIntCols()`
    pub fn setup_int_cols(&mut self, ms: &mut HighsMipSolver) {
        self.intcols = ms.mipdata_.integer_cols.clone();

        pdqsort_by(&mut self.intcols[..], |&c1, &c2| {
            // HiGHS: FP_32BIT_VOLATILE is a no-op on the oracle platform (arm64).
            let lock_score1 = (ms.mipdata_.feastol + ms.mipdata_.uplocks[c1 as usize] as f64)
                * (ms.mipdata_.feastol + ms.mipdata_.downlocks[c1 as usize] as f64);

            let lock_score2 = (ms.mipdata_.feastol + ms.mipdata_.uplocks[c2 as usize] as f64)
                * (ms.mipdata_.feastol + ms.mipdata_.downlocks[c2 as usize] as f64);

            if lock_score1 > lock_score2 {
                return true;
            }
            if lock_score2 > lock_score1 {
                return false;
            }

            let clique_score1 = (ms.mipdata_.feastol
                + ms.mipdata_.cliquetable.get_num_implications_2(c1, true) as f64)
                * (ms.mipdata_.feastol
                    + ms.mipdata_.cliquetable.get_num_implications_2(c1, false) as f64);

            let clique_score2 = (ms.mipdata_.feastol
                + ms.mipdata_.cliquetable.get_num_implications_2(c2, true) as f64)
                * (ms.mipdata_.feastol
                    + ms.mipdata_.cliquetable.get_num_implications_2(c2, false) as f64);

            // C++ `std::make_tuple(..) > std::make_tuple(..)` is `rhs < lhs`.
            let key1 = (
                clique_score1,
                hash_u64_bytes(c1 as u64 as u32, ((c1 as u64) >> 32) as u32),
                c1,
            );
            let key2 = (
                clique_score2,
                hash_u64_bytes(c2 as u64 as u32, ((c2 as u64) >> 32) as u32),
                c2,
            );
            tuple_lt(key2, key1)
        });
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::solveSubMip (lines 77-247)
    /// C++: `bool solveSubMip(HighsMipWorker & worker, const HighsLp & lp, const HighsBasis & basis, double fixingRate, std::vector<double> colLower, std::vector<double> colUpper, HighsInt maxleaves, HighsInt maxnodes, HighsInt stallnodes)`
    // PORT NOTE: `lp` / `basis` are taken by value: call sites pass members of `ms`
    // (`*mipsolver.model_`, `mipdata_->firstrootbasis`) or of the local `heurlp`, and the C++
    // copies both anyway (`HighsLp submip = lp`; the sub-MIP's `rootbasis` is an owned clone).
    // PORT NOTE: the profiling calls (`profiling_->start/stop/solveCall/setSubMip`, the
    // `setProfiling` hand-over, and the `sub_solver_` printf) are dropped: `HighsProfiling` is a
    // no-op stub and `HighsMipSolver::set_profiling` is not yet written.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_sub_mip(
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
    ) -> bool {
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
        submipoptions.objective_bound = ms.mipdata_.workers[worker].upper_limit;

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
        submipsolver.rootbasis = Some(Box::new(basis));
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
                self.try_solution(ms, &submipsolver.solution_, kSolutionSourceSubMip, worker);
            }

            if ms.mipdata_.workers[worker].upper_limit < old_upper_limit {
                // remember fixing rate as good
                ms.mipdata_.workers[worker].update_heur_stats_success_observations(fixing_rate);
            }

            true
        };
        submipsolver.release();
        result
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::determineTargetFixingRate (lines 249-274)
    /// C++: `double determineTargetFixingRate(HighsMipWorker & worker)`
    pub fn determine_target_fixing_rate(&mut self, ms: &mut HighsMipSolver, worker: usize) -> f64 {
        let mut low_fixing_rate = 0.6;
        let mut high_fixing_rate = 0.6;

        {
            let wk: &HighsMipWorker = &ms.mipdata_.workers[worker];
            if self.get_num_infeas_observations(wk) != 0 {
                let infeas_rate =
                    self.get_infeas_observations(wk) / self.get_num_infeas_observations(wk) as f64;
                high_fixing_rate = 0.9 * infeas_rate;
                low_fixing_rate = cpp_min(low_fixing_rate, high_fixing_rate);
            }

            if self.get_num_success_observations(wk) != 0 {
                let success_fixing_rate = self.get_success_observations(wk)
                    / self.get_num_success_observations(wk) as f64;
                low_fixing_rate = cpp_min(low_fixing_rate, 0.9 * success_fixing_rate);
                high_fixing_rate = cpp_max(success_fixing_rate * 1.1, high_fixing_rate);
            }
        }

        // HighsRandom& randgen = parallelLockActive() ? worker.randgen : this->randgen;
        let fixing_rate = if ms.parallel_lock_active() {
            let mut randgen = std::mem::take(&mut ms.mipdata_.workers[worker].randgen);
            let value = randgen.real(low_fixing_rate, high_fixing_rate);
            ms.mipdata_.workers[worker].randgen = randgen;
            value
        } else {
            self.randgen.real(low_fixing_rate, high_fixing_rate)
        };
        // if (!mipsolver.submip) printf("fixing rate: %.2f\n", 100.0 * fixingRate);
        fixing_rate
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::rootReducedCost (lines 314-375)
    /// C++: `void rootReducedCost(HighsMipWorker & worker)`
    pub fn root_reduced_cost(&mut self, ms: &mut HighsMipSolver, worker: usize) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_ROOT_REDCOST);
        let gd = ms.mipdata_.workers[worker].get_global_domain();
        let mut lurking_bounds = ms
            .mipdata_
            .redcostfixing
            .get_lurking_bounds(ms, &ms.mipdata_.domains[gd]);
        if 10 * lurking_bounds.len() < ms.mipdata_.integral_cols.len() {
            return;
        }
        pdqsort_by(&mut lurking_bounds[..], |a, b| a.0 > b.0);

        let mut localdom: HighsDomain =
            detach!(ms.mipdata_.domains[gd], |g| HighsDomain::copy_from(ms, g));

        'body: {
            let mut neighbourhood = HeuristicNeighbourhood::new(ms, &localdom);

            // HiGHS quirk: C++ initialises currCutoff to kHighsInf, but every read follows an
            // assignment inside the loop, so the initial value is never observed.
            let mut curr_cutoff: f64;
            let lower_bound = ms.mipdata_.lower_bound + ms.mipdata_.feastol;

            for domchg in lurking_bounds.iter() {
                curr_cutoff = domchg.0;

                if curr_cutoff <= lower_bound {
                    break;
                }

                if localdom.is_active(&domchg.1) {
                    continue;
                }
                // C++ `changeBound(boundchg)` takes the default reason `Reason::branching()`.
                localdom.change_bound(ms, domchg.1.clone(), HighsDomain_Reason::branching());

                loop {
                    localdom.propagate(ms);
                    if localdom.infeasible() {
                        conflict_analysis_worker(ms, &mut localdom, worker);

                        let new_lower_bound = cpp_max(ms.mipdata_.lower_bound, curr_cutoff);
                        ms.update_lower_bound(new_lower_bound, true, true);

                        localdom.backtrack(ms);
                        if localdom.get_branch_depth() == 0 {
                            break;
                        }
                        neighbourhood.backtracked();
                        continue;
                    }
                    break;
                }
                let fixing_rate = neighbourhood.get_fixing_rate(ms, &localdom);
                if fixing_rate >= 0.5 {
                    break;
                }
                // double gap = (currCutoff - mipsolver.mipdata_->lower_bound) /
                //             std::max(std::abs(mipsolver.mipdata_->lower_bound), 1.0);
                // if (gap < 0.001) break;
            }

            let fixing_rate = neighbourhood.get_fixing_rate(ms, &localdom);
            if fixing_rate < 0.3 {
                break 'body;
            }

            let lp: HighsLp = mip_model!(ms).clone();
            let basis: HighsBasis = ms.mipdata_.firstrootbasis.clone();
            let col_lower = localdom.col_lower_.clone();
            let col_upper = localdom.col_upper_.clone();
            let maxnodes: HighsInt = 200 + (ms.mipdata_.num_nodes / 20) as HighsInt;
            // std::max(50, int(0.05 * (mipsolver.mipdata_->num_leaves))),
            self.solve_sub_mip(
                ms,
                worker,
                lp,
                basis,
                fixing_rate,
                col_lower,
                col_upper,
                500,
                maxnodes,
                12,
            );
        }
        localdom.release(ms);
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::RENS (lines 394-659)
    /// C++: `void RENS(HighsMipWorker & worker, const std::vector<double> & tmp)`
    // HiGHS quirk: C++ takes `tmp` and never reads it.
    pub fn rens(&mut self, ms: &mut HighsMipSolver, worker: usize, tmp: &[f64]) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_RENS);
        let _ = tmp;
        let gd = ms.mipdata_.workers[worker].get_global_domain();
        // return if domain is infeasible
        if ms.mipdata_.domains[gd].infeasible() {
            return;
        }

        let pscost_k = ms.mipdata_.workers[worker].get_pseudocost();
        let mut pscost: HighsPseudocost = detach!(ms.mipdata_.pseudocosts[pscost_k], |p| p.clone());
        let mut heur = HighsSearch::new(ms, gd);
        heur.set_heuristic(true);

        // HiGHS: under the parallel lock `intcols` is a copy; otherwise the member is filtered.
        let par = ms.parallel_lock_active();
        let mut intcols_par: Vec<HighsInt> = if par {
            self.intcols.clone()
        } else {
            Vec::new()
        };
        {
            let intcols_mut: &mut Vec<HighsInt> = if par {
                &mut intcols_par
            } else {
                &mut self.intcols
            };
            intcols_mut.retain(|&i| !ms.mipdata_.domains[gd].is_fixed(i));
        }

        // LP relaxation instantiation
        let lp_k = ms.mipdata_.workers[worker].get_lp_relaxation();
        let mut heurlp = HighsLpRelaxation::copy_from(ms, &ms.mipdata_.lps[lp_k]);
        heurlp.set_mip_worker(worker);
        heurlp.set_profiling();
        // only use the global upper limit as LP limit so that dual proofs are valid
        let worker_upper_limit = ms.mipdata_.workers[worker].upper_limit;
        heurlp.set_objective_limit(ms, worker_upper_limit);
        heurlp.set_adjust_symmetric_branching_col(false);
        // C++ `heur.setLpRelaxation(&heurlp)` is a no-op: the LP is the SearchCtx's.

        let num_col_minus_one = ms.num_col() - 1;
        let _ = heurlp.get_lp_solver().change_cols_bounds_interval(
            0,
            num_col_minus_one,
            &heur.localdom.col_lower_,
            &heur.localdom.col_upper_,
        );
        heur.localdom.clear_changed_cols();
        heur.create_new_node();

        // determine the initial number of unfixed variables fixing rate to decide if the problem
        // is restricted enough to be considered for solving a submip
        let mut max_fixing_rate = self.determine_target_fixing_rate(ms, worker);
        // HiGHS quirk: C++ declares fixingrate = 0.0 and every read follows an assignment.
        let mut fixing_rate: f64;
        // HiGHS quirk: C++ `stop` is never set to true in this function.
        let stop = false;
        let mut target_depth: HighsInt = 1;
        let mut nbacktracks: HighsInt = -1;
        let mut neighbourhood = HeuristicNeighbourhood::new(ms, &heur.localdom);

        'body: {
            // C++ `goto retry` is `continue 'retry`.
            'retry: loop {
                nbacktracks += 1;
                neighbourhood.backtracked();
                if heur.get_current_depth() > target_depth {
                    let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                    if !heur.backtrack_until_depth(&mut ctx, target_depth) {
                        *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                            heur.get_local_lp_iterations();
                        break 'body;
                    }
                }

                debug_assert!(heur.has_node());
                loop {
                    {
                        let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                        heur.evaluate_node(&mut ctx);
                    }
                    if heur.current_node_pruned() {
                        nbacktracks += 1;
                        if ms.mipdata_.domains[gd].infeasible() {
                            *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                                heur.get_local_lp_iterations();
                            break 'body;
                        }

                        let backtracked = {
                            let mut ctx =
                                SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                            heur.backtrack(&mut ctx, true)
                        };
                        if !backtracked {
                            break;
                        }
                        neighbourhood.backtracked();
                        continue;
                    }

                    fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                    if fixing_rate >= max_fixing_rate {
                        break;
                    }
                    if stop {
                        break;
                    }
                    if nbacktracks >= 10 {
                        break;
                    }

                    let mut num_branched: HighsInt = 0;
                    let stop_fixing_rate = cpp_min(
                        1.0 - (1.0 - neighbourhood.get_fixing_rate(ms, &heur.localdom)) * 0.9,
                        max_fixing_rate,
                    );
                    let relaxationsol: Vec<f64> = heurlp.get_solution().col_value.clone();
                    let intcols: &Vec<HighsInt> = if par { &intcols_par } else { &self.intcols };
                    for &i in intcols.iter() {
                        let ci = i as usize;
                        if heur.localdom.col_lower_[ci] == heur.localdom.col_upper_[ci] {
                            continue;
                        }

                        let mut downval = (relaxationsol[ci] + ms.mipdata_.feastol).floor();
                        let mut upval = (relaxationsol[ci] - ms.mipdata_.feastol).ceil();

                        downval = cpp_min(downval, heur.localdom.col_upper_[ci]);
                        upval = cpp_max(upval, heur.localdom.col_lower_[ci]);
                        if heur.localdom.col_lower_[ci] < downval {
                            num_branched += 1;
                            {
                                let mut ctx =
                                    SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                                heur.branch_upwards(&mut ctx, i, downval, downval - 0.5);
                            }
                            heur.localdom.propagate(ms);
                            if heur.localdom.infeasible() {
                                conflict_analysis_worker(ms, &mut heur.localdom, worker);
                                break;
                            }
                        }
                        if heur.localdom.col_upper_[ci] > upval {
                            num_branched += 1;
                            {
                                let mut ctx =
                                    SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                                heur.branch_downwards(&mut ctx, i, upval, upval + 0.5);
                            }
                            heur.localdom.propagate(ms);
                            if heur.localdom.infeasible() {
                                conflict_analysis_worker(ms, &mut heur.localdom, worker);
                                break;
                            }
                        }

                        if neighbourhood.get_fixing_rate(ms, &heur.localdom) >= stop_fixing_rate {
                            break;
                        }
                    }

                    if num_branched == 0 {
                        // sort the variables by their distance towards the value they will be fixed to
                        let frac_len = heurlp.get_fractional_integers().len();
                        let mut fracs = std::mem::take(heurlp.get_fractional_integers_2());
                        {
                            let msr: &HighsMipSolver = &*ms;
                            let localdom_ref: &HighsDomain = &heur.localdom;
                            pdqsort_by(&mut fracs[..], |a, b| {
                                pair_lt(
                                    (
                                        (rens_get_fix_val(msr, localdom_ref, a.0, a.1) - a.1).abs(),
                                        frac_key_hash(a.0, frac_len),
                                    ),
                                    (
                                        (rens_get_fix_val(msr, localdom_ref, b.0, b.1) - b.1).abs(),
                                        frac_key_hash(b.0, frac_len),
                                    ),
                                )
                            });
                        }
                        *heurlp.get_fractional_integers_2() = fracs;

                        let mut change = 0.0;
                        // select a set of fractional variables to fix
                        for k in 0..heurlp.get_fractional_integers().len() {
                            let fracint = heurlp.get_fractional_integers()[k];
                            let fixval = rens_get_fix_val(ms, &heur.localdom, fracint.0, fracint.1);

                            if heur.localdom.col_lower_[fracint.0 as usize] < fixval {
                                num_branched += 1;
                                {
                                    let mut ctx =
                                        SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                                    heur.branch_upwards(&mut ctx, fracint.0, fixval, fracint.1);
                                }
                                heur.localdom.propagate(ms);
                                if heur.localdom.infeasible() {
                                    conflict_analysis_worker(ms, &mut heur.localdom, worker);
                                    break;
                                }

                                fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                            }

                            if heur.localdom.col_upper_[fracint.0 as usize] > fixval {
                                num_branched += 1;
                                {
                                    let mut ctx =
                                        SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                                    heur.branch_downwards(&mut ctx, fracint.0, fixval, fracint.1);
                                }
                                heur.localdom.propagate(ms);
                                if heur.localdom.infeasible() {
                                    conflict_analysis_worker(ms, &mut heur.localdom, worker);
                                    break;
                                }

                                fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                            }

                            if fixing_rate >= max_fixing_rate {
                                break;
                            }

                            change += (fixval - fracint.1).abs();
                            if change >= 0.5 {
                                break;
                            }
                        }
                    }

                    if num_branched == 0 {
                        break;
                    }
                    heurlp.flush_domain(ms, &mut heur.localdom, false);
                }

                // if there is no node left it means we backtracked to the global domain and
                // the subproblem was solved with the dive
                if !heur.has_node() {
                    *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                        heur.get_local_lp_iterations();
                    break 'body;
                }
                // determine the fixing rate to decide if the problem is restricted enough to be
                // considered for solving a submip
                fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                if fixing_rate < 0.1 || (ms.submip && ms.mipdata_.num_improving_sols != 0) {
                    // heur.childselrule = ChildSelectionRule::kBestCost;
                    {
                        // HiGHS #3445 (mip_dfs_node_limit): HighsPrimalHeuristics.cpp lines 639/932 at b1282872a
                        let dfs_node_limit = if ms.options_mip_.mip_dfs_node_limit {
                            K_DEPTH_FIRST_NODE_LIMIT
                        } else {
                            i64::MAX
                        };
                        let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                        heur.set_min_reliable(&mut ctx, 0);
                        heur.solve_depth_first(&mut ctx, 10, dfs_node_limit);
                    }
                    *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                        heur.get_local_lp_iterations();
                    if ms.submip {
                        let local_nodes = *heur.get_local_nodes();
                        ms.mipdata_.num_nodes += local_nodes;
                    }
                    // lpiterations += heur.lpiterations;
                    // pseudocost = heur.pseudocost;
                    break 'body;
                }

                heurlp.remove_obsolete_rows(ms, false);
                let node_reduction_factor: HighsInt = if ms.parallel_lock_active() {
                    cpp_max(1, ms.mipdata_.workers.len() as HighsInt / 4)
                } else {
                    1
                };
                let sub_lp = heurlp.get_lp().clone();
                let sub_basis = heurlp.get_lp_solver().get_basis().clone();
                let sub_col_lower = heur.localdom.col_lower_.clone();
                let sub_col_upper = heur.localdom.col_upper_.clone();
                let maxnodes: HighsInt =
                    (200 + ms.mipdata_.num_nodes / (node_reduction_factor as i64 * 20)) as HighsInt;
                // std::max(50, int(0.05 * (mipsolver.mipdata_->num_leaves))),
                // rusolve `mip_rens_min_fixing_rate`: a RENS sub-MIP with
                // most variables still free is the original problem; SCIP RENS gives up below
                // minfixingrate (Berthold, Math. Prog. Comp. 2014)
                if ms.options_mip_.mip_rens_min_fixing_rate > 0.0
                    && fixing_rate < ms.options_mip_.mip_rens_min_fixing_rate
                {
                    *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                        heur.get_local_lp_iterations();
                    break 'body;
                }
                let solve_sub_mip_return = self.solve_sub_mip(
                    ms,
                    worker,
                    sub_lp,
                    sub_basis,
                    fixing_rate,
                    sub_col_lower,
                    sub_col_upper,
                    500,
                    maxnodes,
                    12,
                );
                if ms.mipdata_.workers[worker].terminator_terminated() {
                    break 'body;
                }
                if !solve_sub_mip_return {
                    let new_lp_iterations = *ms.mipdata_.workers[worker].get_heur_lp_iterations()
                        + heur.get_local_lp_iterations();
                    if new_lp_iterations + ms.mipdata_.heuristic_lp_iterations
                        > 100000
                            + ((ms.mipdata_.total_lp_iterations
                                - ms.mipdata_.heuristic_lp_iterations
                                - ms.mipdata_.sb_lp_iterations)
                                >> 1)
                    {
                        *ms.mipdata_.workers[worker].get_heur_lp_iterations() = new_lp_iterations;
                        break 'body;
                    }

                    target_depth = heur.get_current_depth() / 2;
                    // C++ `||` short-circuits: `checkLimits` (which may run the interrupt
                    // callback) is only evaluated when `targetdepth > 1`.
                    // HiGHS 1.16 #3237 (h16_mip_workers): HighsPrimalHeuristics.cpp subMipShouldRetry
                    // lines 59-63 (1.16.0)
                    if target_depth <= 1
                        || (if ms.options_mip_.h16_mip_workers {
                            h16_sub_mip_check_limits(ms, worker)
                        } else {
                            !ms.parallel_lock_active() && ms.check_limits(0)
                        })
                    {
                        *ms.mipdata_.workers[worker].get_heur_lp_iterations() = new_lp_iterations;
                        break 'body;
                    }
                    max_fixing_rate = fixing_rate * 0.5;
                    // printf("infeasible in root node, trying with lower fixing rate %g\n",
                    //        maxfixingrate);
                    continue 'retry;
                }

                *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                    heur.get_local_lp_iterations();
                break 'body;
            }
        }
        heur.release(ms);
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::RINS (lines 661-969)
    /// C++: `void RINS(HighsMipWorker & worker, const std::vector<double> & relaxationsol)`
    pub fn rins(&mut self, ms: &mut HighsMipSolver, worker: usize, relaxationsol: &[f64]) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_RINS);
        let gd = ms.mipdata_.workers[worker].get_global_domain();
        // return if domain is infeasible
        if ms.mipdata_.domains[gd].infeasible() {
            return;
        }

        if relaxationsol.len() != ms.num_col() as usize {
            return;
        }

        // HiGHS: under the parallel lock `intcols` is a copy; otherwise the member is filtered.
        let par = ms.parallel_lock_active();
        let mut intcols_par: Vec<HighsInt> = if par {
            self.intcols.clone()
        } else {
            Vec::new()
        };
        {
            let intcols_mut: &mut Vec<HighsInt> = if par {
                &mut intcols_par
            } else {
                &mut self.intcols
            };
            intcols_mut.retain(|&i| !ms.mipdata_.domains[gd].is_fixed(i));
        }

        let pscost_k = ms.mipdata_.workers[worker].get_pseudocost();
        let mut pscost: HighsPseudocost = detach!(ms.mipdata_.pseudocosts[pscost_k], |p| p.clone());
        let mut heur = HighsSearch::new(ms, gd);
        heur.set_heuristic(true);

        // LP relaxation instantiation
        let lp_k = ms.mipdata_.workers[worker].get_lp_relaxation();
        let mut heurlp = HighsLpRelaxation::copy_from(ms, &ms.mipdata_.lps[lp_k]);
        heurlp.set_mip_worker(worker);
        heurlp.set_profiling();
        // only use the global upper limit as LP limit so that dual proofs are valid
        let worker_upper_limit = ms.mipdata_.workers[worker].upper_limit;
        heurlp.set_objective_limit(ms, worker_upper_limit);
        heurlp.set_adjust_symmetric_branching_col(false);
        // C++ `heur.setLpRelaxation(&heurlp)` is a no-op: the LP is the SearchCtx's.

        let num_col_minus_one = ms.num_col() - 1;
        let _ = heurlp.get_lp_solver().change_cols_bounds_interval(
            0,
            num_col_minus_one,
            &heur.localdom.col_lower_,
            &heur.localdom.col_upper_,
        );
        heur.localdom.clear_changed_cols();
        heur.create_new_node();

        // determine the initial number of unfixed variables fixing rate to decide if the problem
        // is restricted enough to be considered for solving a submip
        let mut max_fixing_rate = self.determine_target_fixing_rate(ms, worker);
        let min_fixing_rate = 0.25;
        // HiGHS quirk: C++ declares fixingrate = 0.0 and every read follows an assignment.
        let mut fixing_rate: f64;
        // HiGHS quirk: C++ `stop` is never set to true in this function.
        let stop = false;
        let mut nbacktracks: HighsInt = -1;
        let mut target_depth: HighsInt = 1;
        let mut neighbourhood = HeuristicNeighbourhood::new(ms, &heur.localdom);

        'body: {
            // C++ `goto retry` is `continue 'retry`.
            'retry: loop {
                nbacktracks += 1;
                neighbourhood.backtracked();
                // printf("current depth : %" HIGHSINT_FORMAT "   target depth : %" HIGHSINT_FORMAT "\n", heur.getCurrentDepth(), targetdepth);
                if heur.get_current_depth() > target_depth {
                    let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                    if !heur.backtrack_until_depth(&mut ctx, target_depth) {
                        *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                            heur.get_local_lp_iterations();
                        break 'body;
                    }
                }

                debug_assert!(heur.has_node());

                loop {
                    {
                        let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                        heur.evaluate_node(&mut ctx);
                    }
                    if heur.current_node_pruned() {
                        nbacktracks += 1;
                        // printf("backtrack1\n");
                        if ms.mipdata_.domains[gd].infeasible() {
                            *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                                heur.get_local_lp_iterations();
                            break 'body;
                        }

                        let backtracked = {
                            let mut ctx =
                                SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                            heur.backtrack(&mut ctx, true)
                        };
                        if !backtracked {
                            break;
                        }
                        neighbourhood.backtracked();
                        continue;
                    }

                    fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);

                    if stop {
                        break;
                    }
                    if fixing_rate >= max_fixing_rate {
                        break;
                    }
                    if nbacktracks >= 10 {
                        break;
                    }

                    // partition the fractional variables to consider which ones should we fix in
                    // this dive first if there is an incumbent, we dive towards the RINS
                    // neighbourhood
                    let mut fixcand_end: usize = {
                        let incumbent: &[f64] = &ms.mipdata_.incumbent;
                        let feastol = ms.mipdata_.feastol;
                        cpp_partition(
                            heurlp.get_fractional_integers_2().as_mut_slice(),
                            |fracvar: &(HighsInt, f64)| {
                                (relaxationsol[fracvar.0 as usize] - incumbent[fracvar.0 as usize])
                                    .abs()
                                    <= feastol
                            },
                        )
                    };

                    let mut fixtolpsol = true;

                    // no candidates left to fix for getting to the neighbourhood, therefore we
                    // switch to a different diving strategy until the minimal fixing rate is
                    // reached
                    let mut num_branched: HighsInt = 0;
                    if fixcand_end == 0 {
                        fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                        let stop_fixing_rate =
                            cpp_min(max_fixing_rate, 1.0 - (1.0 - fixing_rate) * 0.9);
                        let currlpsol: Vec<f64> = heurlp.get_solution().col_value.clone();
                        let intcols: &Vec<HighsInt> =
                            if par { &intcols_par } else { &self.intcols };
                        for &i in intcols.iter() {
                            let ci = i as usize;
                            if heur.localdom.col_lower_[ci] == heur.localdom.col_upper_[ci] {
                                continue;
                            }

                            if (currlpsol[ci] - ms.mipdata_.incumbent[ci]).abs()
                                <= ms.mipdata_.feastol
                            {
                                let fixval = HighsIntegers::nearest_integer(currlpsol[ci]) as f64;
                                if heur.localdom.col_lower_[ci] < fixval {
                                    num_branched += 1;
                                    {
                                        let mut ctx = SearchCtx::new(
                                            &mut *ms,
                                            worker,
                                            &mut heurlp,
                                            &mut pscost,
                                        );
                                        heur.branch_upwards(&mut ctx, i, fixval, fixval - 0.5);
                                    }
                                    heur.localdom.propagate(ms);
                                    if heur.localdom.infeasible() {
                                        conflict_analysis_worker(ms, &mut heur.localdom, worker);
                                        break;
                                    }

                                    fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                                }
                                if heur.localdom.col_upper_[ci] > fixval {
                                    num_branched += 1;
                                    {
                                        let mut ctx = SearchCtx::new(
                                            &mut *ms,
                                            worker,
                                            &mut heurlp,
                                            &mut pscost,
                                        );
                                        heur.branch_downwards(&mut ctx, i, fixval, fixval + 0.5);
                                    }
                                    heur.localdom.propagate(ms);
                                    if heur.localdom.infeasible() {
                                        conflict_analysis_worker(ms, &mut heur.localdom, worker);
                                        break;
                                    }

                                    fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                                }

                                if fixing_rate >= stop_fixing_rate {
                                    break;
                                }
                            }
                        }

                        if num_branched != 0 {
                            // printf("fixed %d additional cols, old fixing rate: %.2f%%, new fixing rate: %.2f%%\n", ...);
                            heurlp.flush_domain(ms, &mut heur.localdom, false);
                            continue;
                        }

                        // if the RINS neighbourhood achieved a high enough fixing rate by itself we
                        // stop here
                        if fixing_rate >= min_fixing_rate {
                            break;
                        }
                        fixcand_end = heurlp.get_fractional_integers().len();
                        // now sort the variables by their distance towards the value they will be
                        // fixed to
                        fixtolpsol = false;
                    }

                    // now sort the variables by their distance towards the value they will be
                    // fixed to
                    let frac_len = heurlp.get_fractional_integers().len();
                    let fix_to_lp = fixtolpsol;
                    let mut fracs = std::mem::take(heurlp.get_fractional_integers_2());
                    {
                        let msr: &HighsMipSolver = &*ms;
                        let localdom_ref: &HighsDomain = &heur.localdom;
                        pdqsort_by(&mut fracs[..fixcand_end], |a, b| {
                            pair_lt(
                                (
                                    (rins_get_fix_val(
                                        msr,
                                        localdom_ref,
                                        relaxationsol,
                                        fix_to_lp,
                                        a.0,
                                        a.1,
                                    ) - a.1)
                                        .abs(),
                                    frac_key_hash(a.0, frac_len),
                                ),
                                (
                                    (rins_get_fix_val(
                                        msr,
                                        localdom_ref,
                                        relaxationsol,
                                        fix_to_lp,
                                        b.0,
                                        b.1,
                                    ) - b.1)
                                        .abs(),
                                    frac_key_hash(b.0, frac_len),
                                ),
                            )
                        });
                    }
                    *heurlp.get_fractional_integers_2() = fracs;

                    let mut change = 0.0;
                    // select a set of fractional variables to fix
                    for k in 0..fixcand_end {
                        let fracint = heurlp.get_fractional_integers()[k];
                        let fixval = rins_get_fix_val(
                            ms,
                            &heur.localdom,
                            relaxationsol,
                            fixtolpsol,
                            fracint.0,
                            fracint.1,
                        );

                        if heur.localdom.col_lower_[fracint.0 as usize] < fixval {
                            num_branched += 1;
                            {
                                let mut ctx =
                                    SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                                heur.branch_upwards(&mut ctx, fracint.0, fixval, fracint.1);
                            }
                            // HiGHS quirk: this dive does not call `localdom.propagate()` before
                            // the infeasibility test (the RINS neighbourhood loop above does).
                            if heur.localdom.infeasible() {
                                conflict_analysis_worker(ms, &mut heur.localdom, worker);
                                break;
                            }

                            fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                        }

                        if heur.localdom.col_upper_[fracint.0 as usize] > fixval {
                            num_branched += 1;
                            {
                                let mut ctx =
                                    SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                                heur.branch_downwards(&mut ctx, fracint.0, fixval, fracint.1);
                            }
                            // HiGHS quirk: no `localdom.propagate()` here either.
                            if heur.localdom.infeasible() {
                                conflict_analysis_worker(ms, &mut heur.localdom, worker);
                                break;
                            }

                            fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                        }

                        if fixing_rate >= max_fixing_rate {
                            break;
                        }

                        change += (fixval - fracint.1).abs();
                        if change >= 0.5 {
                            break;
                        }
                    }

                    if num_branched == 0 {
                        break;
                    }

                    heurlp.flush_domain(ms, &mut heur.localdom, false);

                    // printf("%d/%d fixed, fixingrate is %g\n", nfixed, ntotal, fixingrate);
                }

                // if there is no node left it means we backtracked to the global domain and the
                // subproblem was solved with the dive
                if !heur.has_node() {
                    *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                        heur.get_local_lp_iterations();
                    break 'body;
                }
                // determine the fixing rate to decide if the problem is restricted enough to be
                // considered for solving a submip

                // printf("fixing rate is %g\n", fixingrate);
                fixing_rate = neighbourhood.get_fixing_rate(ms, &heur.localdom);
                if fixing_rate < 0.1 || (ms.submip && ms.mipdata_.num_improving_sols != 0) {
                    // heur.childselrule = ChildSelectionRule::kBestCost;
                    {
                        // HiGHS #3445 (mip_dfs_node_limit): HighsPrimalHeuristics.cpp lines 639/932 at b1282872a
                        let dfs_node_limit = if ms.options_mip_.mip_dfs_node_limit {
                            K_DEPTH_FIRST_NODE_LIMIT
                        } else {
                            i64::MAX
                        };
                        let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                        heur.set_min_reliable(&mut ctx, 0);
                        heur.solve_depth_first(&mut ctx, 10, dfs_node_limit);
                    }
                    *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                        heur.get_local_lp_iterations();
                    if ms.submip {
                        let local_nodes = *heur.get_local_nodes();
                        ms.mipdata_.num_nodes += local_nodes;
                    }
                    // lpiterations += heur.lpiterations;
                    // pseudocost = heur.pseudocost;
                    break 'body;
                }

                heurlp.remove_obsolete_rows(ms, false);
                let node_reduction_factor: HighsInt = if ms.parallel_lock_active() {
                    cpp_max(1, ms.mipdata_.workers.len() as HighsInt / 4)
                } else {
                    1
                };
                let sub_lp = heurlp.get_lp().clone();
                let sub_basis = heurlp.get_lp_solver().get_basis().clone();
                let sub_col_lower = heur.localdom.col_lower_.clone();
                let sub_col_upper = heur.localdom.col_upper_.clone();
                let maxnodes: HighsInt =
                    (200 + ms.mipdata_.num_nodes / (node_reduction_factor as i64 * 20)) as HighsInt;
                // std::max(50, int(0.05 * (mipsolver.mipdata_->num_leaves))),
                let solve_sub_mip_return = self.solve_sub_mip(
                    ms,
                    worker,
                    sub_lp,
                    sub_basis,
                    fixing_rate,
                    sub_col_lower,
                    sub_col_upper,
                    500,
                    maxnodes,
                    12,
                );
                if ms.mipdata_.workers[worker].terminator_terminated() {
                    break 'body;
                }
                if !solve_sub_mip_return {
                    let new_lp_iterations = *ms.mipdata_.workers[worker].get_heur_lp_iterations()
                        + heur.get_local_lp_iterations();
                    if new_lp_iterations + ms.mipdata_.heuristic_lp_iterations
                        > 100000
                            + ((ms.mipdata_.total_lp_iterations
                                - ms.mipdata_.heuristic_lp_iterations
                                - ms.mipdata_.sb_lp_iterations)
                                >> 1)
                    {
                        *ms.mipdata_.workers[worker].get_heur_lp_iterations() = new_lp_iterations;
                        break 'body;
                    }

                    target_depth = heur.get_current_depth() / 2;
                    // C++ `||` short-circuits: `checkLimits` (which may run the interrupt
                    // callback) is only evaluated when `targetdepth > 1`.
                    // HiGHS 1.16 #3237 (h16_mip_workers): HighsPrimalHeuristics.cpp subMipShouldRetry
                    // lines 59-63 (1.16.0)
                    if target_depth <= 1
                        || (if ms.options_mip_.h16_mip_workers {
                            h16_sub_mip_check_limits(ms, worker)
                        } else {
                            !ms.parallel_lock_active() && ms.check_limits(0)
                        })
                    {
                        *ms.mipdata_.workers[worker].get_heur_lp_iterations() = new_lp_iterations;
                        break 'body;
                    }
                    // printf("infeasible in root node, trying with lower fixing rate\n");
                    max_fixing_rate = fixing_rate * 0.5;
                    continue 'retry;
                }

                *ms.mipdata_.workers[worker].get_heur_lp_iterations() +=
                    heur.get_local_lp_iterations();
                break 'body;
            }
        }
        heur.release(ms);
    }
}

impl HeuristicNeighbourhood {
    /// HiGHS: HighsPrimalHeuristics.cpp::HeuristicNeighbourhood::HeuristicNeighbourhood (lines 285-294)
    /// C++: `HeuristicNeighbourhood(const HighsMipSolver & mipsolver, HighsDomain & localdom)`
    // PORT NOTE: `localdom` is not stored; `get_fixing_rate` takes it (it is the heuristic
    // search's `heur.localdom` or a local domain copy).
    pub fn new(mipsolver: &HighsMipSolver, localdom: &HighsDomain) -> Self {
        let start_checked_changes = localdom.get_domain_change_stack().len();
        let mut num_fixed: HighsInt = 0;
        for &i in mipsolver.mipdata_.integral_cols.iter() {
            if localdom.col_lower_[i as usize] == localdom.col_upper_[i as usize] {
                num_fixed += 1;
            }
        }

        HeuristicNeighbourhood {
            num_fixed,
            fixed_cols: HighsHashTable::new(),
            start_checked_changes,
            n_checked_changes: start_checked_changes,
            num_total: mipsolver.mipdata_.integral_cols.len() as HighsInt - num_fixed,
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HeuristicNeighbourhood::getFixingRate (lines 296-306)
    /// C++: `double getFixingRate()`
    // PORT NOTE: takes `ms` and `localdom`: C++ reads the stored `localdom` reference, and
    // `localdom.variableType(col)` (HighsDomain.h:643) forwards to `mipsolver->variableType`.
    pub fn get_fixing_rate(&mut self, ms: &HighsMipSolver, localdom: &HighsDomain) -> f64 {
        while self.n_checked_changes < localdom.get_domain_change_stack().len() {
            let col = localdom.get_domain_change_stack()[self.n_checked_changes].column;
            self.n_checked_changes += 1;
            if localdom.variable_type(ms, col) == HighsVarType::kContinuous {
                continue;
            }
            if localdom.is_fixed(col) {
                self.fixed_cols.insert(HighsHashTableEntry::new(col, ()));
            }
        }

        if self.num_total != 0 {
            self.fixed_cols.size() as f64 / self.num_total as f64
        } else {
            0.0
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HeuristicNeighbourhood::backtracked (lines 308-311)
    /// C++: `void backtracked()`
    pub fn backtracked(&mut self) {
        self.n_checked_changes = self.start_checked_changes;
        if self.fixed_cols.size() != 0 {
            self.fixed_cols.clear();
        }
    }
}

/// HiGHS: HighsPrimalHeuristics.cpp::calcFixVal (lines 377-392)
/// C++: `static double calcFixVal(double rootchange, double fracval, double cost)`
pub fn calc_fix_val(rootchange: f64, fracval: f64, cost: f64) -> f64 {
    // reinforce direction of this solution away from root
    // solution if the change is at least 0.4
    // otherwise take the direction where the objective gets worse
    // if objective is zero round to nearest integer
    if rootchange >= 0.4 {
        fracval.ceil()
    } else if rootchange <= -0.4 {
        fracval.floor()
    } else if cost > 0.0 {
        fracval.ceil()
    } else if cost < 0.0 {
        fracval.floor()
    } else {
        (fracval + 0.5).floor()
    }
}

#[cfg(test)]
mod tests {
    use super::{calc_fix_val, pair_lt, tuple_lt};

    /// libc++ `std::pair<double, uint64_t>::operator<`: a NaN first element compares neither
    /// less nor greater, so the hash decides (a `partial_cmp`-based comparator would not).
    #[test]
    fn pair_lt_nan_falls_through_to_hash() {
        assert!(pair_lt((f64::NAN, 1), (0.5, 2)));
        assert!(!pair_lt((f64::NAN, 2), (0.5, 1)));
        assert!(pair_lt((0.25, 9), (0.5, 1)));
        assert!(!pair_lt((0.5, 1), (0.5, 1)));
    }

    #[test]
    fn tuple_lt_is_lexicographic_with_nan_fallthrough() {
        assert!(tuple_lt((1.0, 5, 3), (1.0, 5, 4)));
        assert!(!tuple_lt((1.0, 5, 4), (1.0, 5, 4)));
        assert!(tuple_lt((f64::NAN, 4, 9), (2.0, 5, 0)));
        assert!(!tuple_lt((2.0, 4, 9), (1.0, 5, 0)));
    }

    /// RENS / RINS `calcFixVal`: the root-change test wins over the cost sign, and a zero cost
    /// rounds to nearest (`floor(x + 0.5)`).
    #[test]
    fn calc_fix_val_branches() {
        assert_eq!(calc_fix_val(0.4, 2.3, -1.0), 3.0);
        assert_eq!(calc_fix_val(-0.4, 2.7, 1.0), 2.0);
        assert_eq!(calc_fix_val(0.0, 2.3, 1.0), 3.0);
        assert_eq!(calc_fix_val(0.0, 2.7, -1.0), 2.0);
        assert_eq!(calc_fix_val(0.0, 2.5, 0.0), 3.0);
        assert_eq!(calc_fix_val(0.0, 2.4, 0.0), 2.0);
    }
}
