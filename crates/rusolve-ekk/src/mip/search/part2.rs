//! HiGHS source: mip/HighsSearch.cpp (lines 1099-1998)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::mip::symmetry_dynamic::S3PathNode;
use crate::mip::symmetry_dynamic::hook::{s3_child_path, s3_propagate};
use crate::prelude::*;
use std::sync::Arc;

/// C++ `currnode.branchingdecision.boundtype` / `boundval` and `currnode.other_child_lb`
/// assignments of `HighsSearch::branch` (the same three writes in each child-selection arm).
fn set_branch_decision(
    node: &mut HighsSearch_NodeData,
    boundtype: HighsBoundType,
    boundval: f64,
    other_child_lb: f64,
) {
    node.branchingdecision.boundtype = boundtype;
    node.branchingdecision.boundval = boundval;
    node.other_child_lb = other_child_lb;
}

/// C++ `assert`s after `localdom.backtrack()` in `backtrack`, `backtrackPlunge` and
/// `backtrackUntilDepth`: the undone change matches the branching decision of `node`.
fn debug_assert_undone_branching(node: &HighsSearch_NodeData, branchchg: &HighsDomainChange) {
    debug_assert!(
        (branchchg.boundtype == HighsBoundType::kLower
            && branchchg.boundval >= node.branchingdecision.boundval)
            || (branchchg.boundtype == HighsBoundType::kUpper
                && branchchg.boundval <= node.branchingdecision.boundval)
    );
    debug_assert!(branchchg.boundtype == node.branchingdecision.boundtype);
    debug_assert!(branchchg.column == node.branchingdecision.column);
}

/// C++ flip of the exhausted node's decision (`backtrack`, `backtrackPlunge`,
/// `backtrackUntilDepth`): lower becomes upper at `floor(boundval - 0.5)`, upper becomes
/// lower at `ceil(boundval + 0.5)`; a fallback branch (decision bound equal to the branching
/// point) moves the branching point along.
fn flip_branching_decision(node: &mut HighsSearch_NodeData, fallbackbranch: bool) {
    if node.branchingdecision.boundtype == HighsBoundType::kLower {
        node.branchingdecision.boundtype = HighsBoundType::kUpper;
        node.branchingdecision.boundval = (node.branchingdecision.boundval - 0.5).floor();
    } else {
        node.branchingdecision.boundtype = HighsBoundType::kLower;
        node.branchingdecision.boundval = (node.branchingdecision.boundval + 0.5).ceil();
    }

    if fallbackbranch {
        node.branching_point = node.branchingdecision.boundval;
    }
}

impl HighsSearch {
    /// C++ `if (nodestack.back().stabilizerOrbits) orbitalFixing(localdom); else
    /// getSymmetries().propagateOrbitopes(localdom);` (the re-propagation of `backtrack` and
    /// `backtrackPlunge`).
    ///
    /// returns true when `s3_propagate` detected infeasibility (always false
    /// with mip_symmetry_dynamic off; the 1.15.1 calls report through `localdom` only).
    fn propagate_node_symmetries(&mut self, ms: &mut HighsMipSolver) -> bool {
        if ms.options_mip_.mip_symmetry_dynamic {
            // mip_symmetry_dynamic: SCIP sym_lexorbred.c lines 1091-1148 at 32dd425
            // (G14: the node on top of the stack, SCIP depth = current depth - 1)
            let path = self.nodestack.last().unwrap().s3_path.clone();
            let depth = self.get_current_depth() - 1;
            return self.s3_propagate_localdom(ms, path, depth);
        }
        let orbits = self.nodestack.last().unwrap().stabilizer_orbits.clone();
        match orbits {
            Some(orbits) => {
                self.orbital_fixing(ms, &orbits);
            }
            None => {
                self.propagate_orbitopes(ms);
            }
        }
        false
    }

    /// SCIP 11.0.0 (32dd425): sym_lexorbred.c::symhdlrPropLexOrbRed (lines 1091-1148)
    ///
    /// rusolve mip_symmetry_dynamic: `s3_propagate` on `self.localdom` for the node with
    /// branching path `path` at SCIP depth `depth`. Not run in a heuristic search.
    /// Returns true when infeasible. Called only under mip_symmetry_dynamic.
    fn s3_propagate_localdom(
        &mut self,
        ms: &mut HighsMipSolver,
        path: Option<Arc<S3PathNode>>,
        depth: HighsInt,
    ) -> bool {
        if self.inheuristic {
            return false;
        }
        s3_propagate(
            ms,
            &mut self.localdom,
            path.as_ref(),
            depth,
            &mut self.s3_stats,
        )
    }

    /// C++ `localdom` conflict analysis of a pruned branch (`conflictAnalysis(getConflictPool(),
    /// mipworker.getGlobalDomain(), pseudocost)`).
    fn analyze_local_conflict(&mut self, ctx: &mut SearchCtx<'_>) {
        let conflict_pool = self.get_conflict_pool(ctx);
        let globaldom = ctx.worker().get_global_domain();
        self.localdom
            .conflict_analysis(ctx.ms, conflict_pool, globaldom, ctx.pseudocost);
    }

    /// C++ inner `while` body of `backtrack` and `backtrackPlunge` for a node with
    /// `opensubtrees == 0` when `nodestack.size() > 1`: pop it, undo its branching change and
    /// re-propagate the parent node.
    fn pop_exhausted_node(&mut self, ctx: &mut SearchCtx<'_>) {
        self.nodestack.pop();
        let branchchg = self.localdom.backtrack(ctx.ms);

        if self.nodestack.last().unwrap().opensubtrees != 0 {
            self.count_tree_weight = self.nodestack.last().unwrap().skip_depth_count == 0;
            // repropagate the node, as it may have become infeasible due to conflicts
            let old_num_domchgs = self.localdom.get_num_domain_changes();
            let old_num_changed_cols = self.localdom.get_changed_cols().len();
            self.localdom.propagate(ctx.ms);
            let mut s3_infeasible = false;
            if !self.localdom.infeasible()
                && old_num_domchgs != self.localdom.get_num_domain_changes()
            {
                s3_infeasible = self.propagate_node_symmetries(ctx.ms);
            }
            if self.localdom.infeasible() || s3_infeasible {
                self.localdom.clear_changed_cols_2(old_num_changed_cols);
                if self.count_tree_weight {
                    let depth = self.get_current_depth();
                    self.treeweight += ldexp_f64(1.0, -depth);
                }
                self.nodestack.last_mut().unwrap().opensubtrees = 0;
            }
        }

        debug_assert_undone_branching(self.nodestack.last().unwrap(), &branchchg);
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::branch (lines 1099-1501)
    /// C++: `NodeResult branch()`
    pub fn branch(&mut self, ctx: &mut SearchCtx<'_>) -> HighsSearch_NodeResult {
        debug_assert!(self.localdom.get_changed_cols().is_empty());

        debug_assert!(self.nodestack.last().unwrap().opensubtrees == 2);
        self.nodestack.last_mut().unwrap().branchingdecision.column = -1;
        self.inbranching = true;

        let minrel: HighsInt = ctx.pseudocost.get_min_reliable();
        let mut child_lb: f64 = self.get_current_lower_bound();
        let mut result = HighsSearch_NodeResult::kOpen;
        while self.nodestack.last().unwrap().opensubtrees == 2
            && HighsLpRelaxation::scaled_optimal(ctx.lp.get_status())
            && !ctx.lp.get_fractional_integers().is_empty()
        {
            let mut sbmaxiters: i64 = 0;
            if minrel > 0 {
                let sbiters: i64 = self.get_strong_branching_lp_iterations(ctx.ms);
                sbmaxiters = 100000
                    + ((self.get_total_lp_iterations(ctx.ms)
                        - self.get_heuristic_lp_iterations(ctx.ms)
                        - self.get_strong_branching_lp_iterations(ctx.ms))
                        >> 1);
                // HiGHS #3213 (mip_sb_budget_pre_restart): HighsSearch.cpp lines 1118 at 4a638fa93
                if ctx.ms.options_mip_.mip_sb_budget_pre_restart
                    && ctx.ms.mipdata_.num_restarts <= 2
                {
                    sbmaxiters >>= 2;
                }
                if sbiters > sbmaxiters {
                    ctx.pseudocost.set_min_reliable(0);
                } else if sbiters > (sbmaxiters >> 1) {
                    let reductionratio: f64 = ((sbiters - (sbmaxiters >> 1)) as f64)
                        / ((sbmaxiters - (sbmaxiters >> 1)) as f64);
                    let minrelreduced: HighsInt =
                        (minrel as f64 - reductionratio * ((minrel - 1) as f64)) as HighsInt;
                    ctx.pseudocost
                        .set_min_reliable(cpp_min(minrel, minrelreduced));
                }
            }

            let degeneracy_fac: f64 = ctx.lp.compute_lp_degneracy(&self.localdom);
            ctx.pseudocost.set_degeneracy_factor(degeneracy_fac);
            if degeneracy_fac >= 10.0 {
                ctx.pseudocost.set_min_reliable(0);
            }
            let mut down_node_lb: f64 = self.get_current_lower_bound();
            let mut up_node_lb: f64 = self.get_current_lower_bound();
            let branchcand = self.select_branching_candidate(
                ctx,
                sbmaxiters,
                &mut down_node_lb,
                &mut up_node_lb,
            );
            child_lb = self.nodestack.last().unwrap().lower_bound;
            if branchcand != -1 {
                let branching = ctx.lp.get_fractional_integers()[branchcand as usize];
                {
                    let currnode = self.nodestack.last_mut().unwrap();
                    currnode.branchingdecision.column = branching.0;
                    currnode.branching_point = branching.1;
                }

                let col = branching.0;
                let bp: f64 = self.nodestack.last().unwrap().branching_point;
                let feastol: f64 = self.get_feas_tol(ctx);
                let epsilon: f64 = self.get_epsilon(ctx);

                match self.childselrule {
                    HighsSearch_ChildSelectionRule::kUp => {
                        set_branch_decision(
                            self.nodestack.last_mut().unwrap(),
                            HighsBoundType::kLower,
                            bp.ceil(),
                            down_node_lb,
                        );
                        child_lb = up_node_lb;
                    }
                    HighsSearch_ChildSelectionRule::kDown => {
                        set_branch_decision(
                            self.nodestack.last_mut().unwrap(),
                            HighsBoundType::kUpper,
                            bp.floor(),
                            up_node_lb,
                        );
                        child_lb = down_node_lb;
                    }
                    HighsSearch_ChildSelectionRule::kRootSol => {
                        let mut downprio = ctx.pseudocost.get_avg_inferences_down(col) + epsilon;
                        let mut upprio = ctx.pseudocost.get_avg_inferences_up(col) + epsilon;
                        let downval = bp.floor();
                        let upval = bp.ceil();
                        if !self.subrootsol.is_empty() {
                            let mut rootsol = self.subrootsol[col as usize];
                            if rootsol < downval {
                                rootsol = downval;
                            } else if rootsol > upval {
                                rootsol = upval;
                            }

                            upprio *= 1.0 + (bp - rootsol);
                            downprio *= 1.0 + (rootsol - bp);
                        } else {
                            if self.nodestack.last().unwrap().lp_objective != -kHighsInf {
                                self.subrootsol = ctx.lp.get_solution().col_value.clone();
                            }
                            if !self.get_root_lp_sol(ctx).is_empty() {
                                let mut rootsol = self.get_root_lp_sol(ctx)[col as usize];
                                if rootsol < downval {
                                    rootsol = downval;
                                } else if rootsol > upval {
                                    rootsol = upval;
                                }

                                upprio *= 1.0 + (bp - rootsol);
                                downprio *= 1.0 + (rootsol - bp);
                            }
                        }
                        if upprio + epsilon >= downprio {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kLower,
                                upval,
                                down_node_lb,
                            );
                            child_lb = up_node_lb;
                        } else {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kUpper,
                                downval,
                                up_node_lb,
                            );
                            child_lb = down_node_lb;
                        }
                    }
                    HighsSearch_ChildSelectionRule::kObj => {
                        if ctx.ms.col_cost_2(col) >= 0.0 {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kLower,
                                bp.ceil(),
                                down_node_lb,
                            );
                            child_lb = up_node_lb;
                        } else {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kUpper,
                                bp.floor(),
                                up_node_lb,
                            );
                            child_lb = down_node_lb;
                        }
                    }
                    HighsSearch_ChildSelectionRule::kRandom => {
                        if self.random.bit() {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kLower,
                                bp.ceil(),
                                down_node_lb,
                            );
                            child_lb = up_node_lb;
                        } else {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kUpper,
                                bp.floor(),
                                up_node_lb,
                            );
                            child_lb = down_node_lb;
                        }
                    }
                    HighsSearch_ChildSelectionRule::kBestCost => {
                        if ctx.pseudocost.get_pseudocost_up(col, bp, feastol)
                            > ctx.pseudocost.get_pseudocost_down(col, bp, feastol)
                        {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kUpper,
                                bp.floor(),
                                up_node_lb,
                            );
                            child_lb = down_node_lb;
                        } else {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kLower,
                                bp.ceil(),
                                down_node_lb,
                            );
                            child_lb = up_node_lb;
                        }
                    }
                    HighsSearch_ChildSelectionRule::kWorstCost => {
                        // HiGHS quirk: the worst-cost rule compares the 2-argument
                        // pseudocosts (default offset), not the feasibility-tolerance ones.
                        if ctx.pseudocost.get_pseudocost_up_2(col, bp)
                            >= ctx.pseudocost.get_pseudocost_down_2(col, bp)
                        {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kLower,
                                bp.ceil(),
                                down_node_lb,
                            );
                            child_lb = up_node_lb;
                        } else {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kUpper,
                                bp.floor(),
                                up_node_lb,
                            );
                            child_lb = down_node_lb;
                        }
                    }
                    HighsSearch_ChildSelectionRule::kDisjunction => {
                        let numnodesup: i64 = self.get_node_queue(ctx).num_nodes_up(col);
                        let numnodesdown: i64 = self.get_node_queue(ctx).num_nodes_down(col);
                        if numnodesup > numnodesdown {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kLower,
                                bp.ceil(),
                                down_node_lb,
                            );
                            child_lb = up_node_lb;
                        } else if numnodesdown > numnodesup {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kUpper,
                                bp.floor(),
                                up_node_lb,
                            );
                            child_lb = down_node_lb;
                        } else if ctx.ms.col_cost_2(col) >= 0.0 {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kLower,
                                bp.ceil(),
                                down_node_lb,
                            );
                            child_lb = up_node_lb;
                        } else {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kUpper,
                                bp.floor(),
                                up_node_lb,
                            );
                            child_lb = down_node_lb;
                        }
                    }
                    HighsSearch_ChildSelectionRule::kHybridInferenceCost => {
                        let upval = bp.ceil();
                        let downval = bp.floor();
                        let up_score = (1.0 + ctx.pseudocost.get_avg_inferences_up(col))
                            / ctx.pseudocost.get_pseudocost_up(col, bp, feastol);
                        let down_score = (1.0 + ctx.pseudocost.get_avg_inferences_down(col))
                            / ctx.pseudocost.get_pseudocost_down(col, bp, feastol);

                        if up_score >= down_score {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kLower,
                                upval,
                                down_node_lb,
                            );
                            child_lb = up_node_lb;
                        } else {
                            set_branch_decision(
                                self.nodestack.last_mut().unwrap(),
                                HighsBoundType::kUpper,
                                downval,
                                up_node_lb,
                            );
                            child_lb = down_node_lb;
                        }
                    }
                }
                result = HighsSearch_NodeResult::kBranched;
                break;
            }

            debug_assert!(!self.localdom.get_changed_cols().is_empty());
            result = self.evaluate_node(ctx);
            if result == HighsSearch_NodeResult::kSubOptimal {
                break;
            }
        }
        self.inbranching = false;
        ctx.pseudocost.set_min_reliable(minrel);
        ctx.pseudocost.set_degeneracy_factor(1.0);

        debug_assert!(
            self.nodestack.last().unwrap().opensubtrees == 2
                || self.nodestack.last().unwrap().opensubtrees == 0
        );

        if self.nodestack.last().unwrap().opensubtrees != 2
            || result == HighsSearch_NodeResult::kSubOptimal
        {
            return result;
        }

        if self.nodestack.last().unwrap().branchingdecision.column == -1 {
            let mut bestscore: f64 = -1.0;
            // solution branching failed, so choose any integer variable to branch
            // on in case we have a different solution status could happen due to a
            // fail in the LP solution process
            ctx.pseudocost.set_degeneracy_factor(1e6);

            let integral_cols: Vec<HighsInt> = ctx.ms.mipdata_.integral_cols.clone();
            for &i in integral_cols.iter() {
                let iu = i as usize;
                if self.localdom.col_upper_[iu] - self.localdom.col_lower_[iu] < 0.5 {
                    continue;
                }

                let lo: f64 = self.localdom.col_lower_[iu];
                let up: f64 = self.localdom.col_upper_[iu];
                let fracval: f64 = if lo != -kHighsInf && up != kHighsInf {
                    (0.5 * (lo + up + 0.5)).floor() + 0.5
                } else if lo != -kHighsInf {
                    lo + 0.5
                } else if up != kHighsInf {
                    up - 0.5
                } else {
                    0.5
                };

                let score: f64 = ctx.pseudocost.get_score_2(i, fracval);
                debug_assert!(score >= 0.0);

                if score > bestscore {
                    bestscore = score;
                    let cost: f64 =
                        if HighsLpRelaxation::unscaled_dual_feasible(ctx.lp.get_status()) {
                            ctx.lp.get_solution().col_dual[iu]
                        } else {
                            ctx.ms.col_cost_2(i)
                        };
                    let feastol: f64 = self.get_feas_tol(ctx);
                    let branch_upwards: bool =
                        if cost.abs() > feastol && self.get_cutoff_bound(ctx) < kHighsInf {
                            // branch in direction of worsening cost first in case the column has
                            // cost and we do have an upper bound
                            cost > 0.0
                        } else if ctx.pseudocost.get_avg_inferences_up(i)
                            > ctx.pseudocost.get_avg_inferences_down(i) + feastol
                        {
                            // column does not have (reduced) cost above tolerance so branch in
                            // direction of more inferences
                            true
                        } else if ctx.pseudocost.get_avg_inferences_up(i)
                            < ctx.pseudocost.get_avg_inferences_down(i) - feastol
                        {
                            false
                        } else {
                            // number of inferences give a tie, so we branch in the direction that
                            // does have a less recent domain change to avoid branching the same
                            // integer column into the same direction over and over
                            let mut col_lower_pos: HighsInt = 0;
                            let mut col_upper_pos: HighsInt = 0;
                            let num_domchgs = self.localdom.get_num_domain_changes();
                            self.localdom
                                .get_col_lower_pos(i, num_domchgs, &mut col_lower_pos);
                            self.localdom
                                .get_col_upper_pos(i, num_domchgs, &mut col_upper_pos);
                            col_lower_pos <= col_upper_pos
                        };

                    let currnode = self.nodestack.last_mut().unwrap();
                    if branch_upwards {
                        let upval = fracval.ceil();
                        currnode.branching_point = upval;
                        currnode.branchingdecision.boundtype = HighsBoundType::kLower;
                        currnode.branchingdecision.column = i;
                        currnode.branchingdecision.boundval = upval;
                    } else {
                        let downval = fracval.floor();
                        currnode.branching_point = downval;
                        currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
                        currnode.branchingdecision.column = i;
                        currnode.branchingdecision.boundval = downval;
                    }
                }
            }

            ctx.pseudocost.set_degeneracy_factor(1.0);
        }

        if self.nodestack.last().unwrap().branchingdecision.column == -1 {
            if ctx.lp.get_status() == HighsLpRelaxation_Status::kOptimal {
                // if the LP was solved to optimality and all columns are fixed, then this
                // particular assignment is not feasible or has a worse objective in the
                // original space, otherwise the node would not be open. Hence we prune
                // this particular assignment
                self.nodestack.last_mut().unwrap().opensubtrees = 0;
                result = HighsSearch_NodeResult::kLpInfeasible;
                return result;
            }
            ctx.lp.set_iteration_limit(kHighsIInf);

            // create a fresh LP only with model rows since all integer columns are
            // fixed, the cutting planes are not required and the LP could not be solved
            // so we want to make it as easy as possible
            //
            // LP relaxation instantiation
            let mut lp_copy = HighsLpRelaxation::new(&*ctx.ms);
            lp_copy.set_profiling();
            lp_copy.load_model(&*ctx.ms);
            let num_col: HighsInt = ctx.ms.num_col();
            lp_copy.get_lp_solver().change_cols_bounds_interval(
                0,
                num_col - 1,
                &self.localdom.col_lower_,
                &self.localdom.col_upper_,
            );

            // C++ `std::swap(tmpLp, lp)`: the search runs on `lpCopy` for the re-evaluation.
            // The outer `lp` (detached from `lps[lp_]` by the caller) is lent back to its slot
            // for the whole swap: `sub.lend_back` lends only `lpCopy`, which no slot holds, and
            // `addIncumbent` reached from these evaluations reads every `lps[..]`.
            let worker = ctx.worker;
            let ms = &mut *ctx.ms;
            let outer_lp = &mut *ctx.lp;
            let pseudocost = &mut *ctx.pseudocost;
            return lend_in!(ms.mipdata_.lps, outer_lp, {
                let mut sub = SearchCtx::new(&mut *ms, worker, &mut lp_copy, &mut *pseudocost);

                // reevaluate the node with LP presolve enabled
                sub.lp
                    .get_lp_solver()
                    .set_option_value_str("presolve", kHighsOnString);
                result = self.evaluate_node(&mut sub);

                if result == HighsSearch_NodeResult::kOpen {
                    // LP still not solved, reevaluate with primal simplex
                    sub.lp.get_lp_solver().clear_solver();
                    sub.lp
                        .get_lp_solver()
                        .set_option_value_int("simplex_strategy", kSimplexStrategyPrimal);
                    result = self.evaluate_node(&mut sub);
                    sub.lp
                        .get_lp_solver()
                        .set_option_value_int("simplex_strategy", kSimplexStrategyDual);
                    if result == HighsSearch_NodeResult::kOpen {
                        // LP still not solved, reevaluate with IPM instead of simplex
                        sub.lp.get_lp_solver().clear_solver();
                        sub.lp.get_lp_solver().set_option_value_str("solver", "ipm");
                        result = self.evaluate_node(&mut sub);

                        if result == HighsSearch_NodeResult::kOpen {
                            crate::highs_log_user!(
                                &sub.ms.options_mip_.log_options,
                                HighsLogType::kWarning,
                                "Failed to solve node with all integer columns fixed. Declaring node infeasible.\n"
                            );
                            // LP still not solved, give up and declare as infeasible
                            self.nodestack.last_mut().unwrap().opensubtrees = 0;
                            result = HighsSearch_NodeResult::kLpInfeasible;
                        }
                    }
                }

                // restore old lp relaxation (the swap ends with `sub`)
                result
            });
        }

        // finally open a new node with the branching decision added
        // and remember that we have one open subtree left
        let domchg_pos: HighsInt = self.localdom.get_domain_change_stack().len() as HighsInt;

        let branching_decision = self.nodestack.last().unwrap().branchingdecision.clone();
        if ctx.ms.options_mip_.mip_symmetry_dynamic {
            // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at 32dd425
            // (G15: child path built before the branching change_bound, D3; no stabilizer
            // orbits under the flag, Q2/G9)
            let child_path = s3_child_path(
                &*ctx.ms,
                &self.localdom,
                self.nodestack.last().unwrap().s3_path.as_ref(),
                branching_decision.column,
                self.get_current_depth(),
            );
            self.nodestack.last_mut().unwrap().s3_child_path = child_path.clone();
            self.localdom
                .change_bound(ctx.ms, branching_decision, HighsDomain_Reason::branching());
            let currnode = self.nodestack.last_mut().unwrap();
            currnode.opensubtrees = 1;
            let parent_lb = cpp_max(child_lb, currnode.lower_bound);
            let estimate = currnode.estimate;
            let parent_basis = currnode.node_basis.clone();
            self.nodestack.push(HighsSearch_NodeData::new(
                parent_lb,
                estimate,
                parent_basis,
                None,
            ));
            self.nodestack.last_mut().unwrap().s3_path = child_path;
        } else {
            let pass_stabilizer_to_child_node =
                self.orbits_valid_in_child_node(ctx.ms, &branching_decision);
            self.localdom
                .change_bound(ctx.ms, branching_decision, HighsDomain_Reason::branching());
            let currnode = self.nodestack.last_mut().unwrap();
            currnode.opensubtrees = 1;
            let parent_lb = cpp_max(child_lb, currnode.lower_bound);
            let estimate = currnode.estimate;
            let parent_basis = currnode.node_basis.clone();
            let parent_orbits = if pass_stabilizer_to_child_node {
                currnode.stabilizer_orbits.clone()
            } else {
                None
            };
            self.nodestack.push(HighsSearch_NodeData::new(
                parent_lb,
                estimate,
                parent_basis,
                parent_orbits,
            ));
        }
        self.nodestack.last_mut().unwrap().domgchg_stack_pos = domchg_pos;

        HighsSearch_NodeResult::kBranched
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::backtrack (lines 1503-1621)
    /// C++: `bool backtrack(bool recoverBasis = true)`
    /// C++ default arguments: `recover_basis = true`
    pub fn backtrack(&mut self, ctx: &mut SearchCtx<'_>, recover_basis: bool) -> bool {
        if self.nodestack.is_empty() {
            return false;
        }
        debug_assert!(!self.nodestack.is_empty());
        debug_assert!(self.nodestack.last().unwrap().opensubtrees == 0);
        loop {
            while self.nodestack.last().unwrap().opensubtrees == 0 {
                self.count_tree_weight = true;
                self.depthoffset += self.nodestack.last().unwrap().skip_depth_count as HighsInt;
                if self.nodestack.len() == 1 {
                    if recover_basis {
                        if let Some(basis) = self.nodestack.last().unwrap().node_basis.clone() {
                            ctx.lp.set_stored_basis(Some(basis));
                        }
                    }
                    self.nodestack.pop();
                    self.localdom.backtrack_to_global(ctx.ms);
                    ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
                    if recover_basis {
                        ctx.lp.recover_basis();
                    }
                    return false;
                }

                self.pop_exhausted_node(ctx);
            }

            let currnode = self.nodestack.last_mut().unwrap();
            debug_assert!(currnode.opensubtrees == 1);
            currnode.opensubtrees = 0;
            let fallbackbranch = currnode.branchingdecision.boundval == currnode.branching_point;
            let domchg_pos = self.localdom.get_domain_change_stack().len() as HighsInt;
            flip_branching_decision(currnode, fallbackbranch);

            let num_changed_cols = self.localdom.get_changed_cols().len();
            let decision = self.nodestack.last().unwrap().branchingdecision.clone();
            let s3 = ctx.ms.options_mip_.mip_symmetry_dynamic;
            let pass_stabilizer_to_child_node = if s3 {
                // mip_symmetry_dynamic: SCIP sym_lexorbred.c lines 1091-1148 at 32dd425
                // (no stabilizer orbits under the flag, Q2/G9)
                false
            } else {
                self.orbits_valid_in_child_node(&*ctx.ms, &decision)
            };
            // the sibling reuses the parent's child path
            let s3_sibling_path = if s3 {
                self.nodestack.last().unwrap().s3_child_path.clone()
            } else {
                None
            };
            self.localdom
                .change_bound(ctx.ms, decision, HighsDomain_Reason::branching());
            let currnode = self.nodestack.last().unwrap();
            let nodelb = cpp_max(currnode.lower_bound, currnode.other_child_lb);
            let estimate = currnode.estimate;
            let parent_basis = currnode.node_basis.clone();
            let parent_orbits = if pass_stabilizer_to_child_node {
                currnode.stabilizer_orbits.clone()
            } else {
                None
            };

            let mut prune = nodelb > self.get_cutoff_bound(ctx) || self.localdom.infeasible();
            if !prune {
                self.localdom.propagate(ctx.ms);
                prune = self.localdom.infeasible();
                if prune {
                    self.analyze_local_conflict(ctx);
                }
            }
            if !prune {
                if s3 {
                    // mip_symmetry_dynamic: SCIP sym_lexorbred.c lines 1091-1148 at
                    // 32dd425 (sibling child: SCIP depth = current depth)
                    let depth = self.get_current_depth();
                    prune = self.s3_propagate_localdom(ctx.ms, s3_sibling_path.clone(), depth)
                        || self.localdom.infeasible();
                } else {
                    self.propagate_orbitopes(ctx.ms);
                    prune = self.localdom.infeasible();
                }
            }
            if !prune && pass_stabilizer_to_child_node {
                if let Some(orbits) = self.nodestack.last().unwrap().stabilizer_orbits.clone() {
                    self.orbital_fixing(ctx.ms, &orbits);
                    prune = self.localdom.infeasible();
                }
            }
            if prune {
                self.localdom.backtrack(ctx.ms);
                self.localdom.clear_changed_cols_2(num_changed_cols);
                if self.count_tree_weight {
                    let depth = self.get_current_depth();
                    self.treeweight += ldexp_f64(1.0, -depth);
                }
                continue;
            }
            self.nodestack.push(HighsSearch_NodeData::new(
                nodelb,
                estimate,
                parent_basis,
                parent_orbits,
            ));
            if s3 {
                // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at 32dd425
                self.nodestack.last_mut().unwrap().s3_path = s3_sibling_path;
            }

            ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
            self.nodestack.last_mut().unwrap().domgchg_stack_pos = domchg_pos;
            self.eval_origin = crate::work::LP_ORIGIN_SIBLING;
            break;
        }

        if recover_basis {
            if let Some(basis) = self.nodestack.last().unwrap().node_basis.clone() {
                ctx.lp.set_stored_basis(Some(basis));
                ctx.lp.recover_basis();
            }
        }

        true
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::backtrackPlunge (lines 1623-1803)
    /// HiGHS 1.16 (#3192): HighsSearch.cpp::backtrackPlunge (lines 1619-1741), under `h16_mip_workers`
    /// C++: `bool backtrackPlunge(HighsNodeQueue & nodequeue)`
    // PORT NOTE: `nodequeue` is a selector, as in `current_node_to_queue`.
    pub fn backtrack_plunge(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        nodequeue: HighsSearch_NodeQueueRef,
    ) -> bool {
        if self.nodestack.is_empty() {
            return false;
        }
        debug_assert!(!self.nodestack.is_empty());
        debug_assert!(self.nodestack.last().unwrap().opensubtrees == 0);
        // HiGHS 1.16 #3192 (h16_mip_workers): HighsSearch.cpp backtrackPlunge lines 1619-1741 (1.16.0)
        let h16_mip_workers = ctx.ms.options_mip_.h16_mip_workers;

        loop {
            while self.nodestack.last().unwrap().opensubtrees == 0 {
                self.count_tree_weight = true;
                self.depthoffset += self.nodestack.last().unwrap().skip_depth_count as HighsInt;

                if self.nodestack.len() == 1 {
                    if let Some(basis) = self.nodestack.last().unwrap().node_basis.clone() {
                        ctx.lp.set_stored_basis(Some(basis));
                    }
                    self.nodestack.pop();
                    self.localdom.backtrack_to_global(ctx.ms);
                    ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
                    ctx.lp.recover_basis();
                    return false;
                }

                self.pop_exhausted_node(ctx);
            }

            let currnode = self.nodestack.last_mut().unwrap();
            debug_assert!(currnode.opensubtrees == 1);
            currnode.opensubtrees = 0;
            let fallbackbranch = currnode.branchingdecision.boundval == currnode.branching_point;
            let column = currnode.branchingdecision.column;
            let score_branchpoint = if fallbackbranch {
                0.5
            } else {
                currnode.branching_point
            };
            // HiGHS 1.16 #3192 (h16_mip_workers): 1.16 drops the node score (the ancestor check
            // below is gone); the pseudocost score getters are pure, so 0.0 is inert.
            let node_score: f64 = if h16_mip_workers {
                0.0
            } else if currnode.branchingdecision.boundtype == HighsBoundType::kLower {
                ctx.pseudocost.get_score_down(column, score_branchpoint)
            } else {
                ctx.pseudocost.get_score_up(column, score_branchpoint)
            };
            flip_branching_decision(currnode, fallbackbranch);

            let domchg_pos = self.localdom.get_domain_change_stack().len() as HighsInt;
            let num_changed_cols = self.localdom.get_changed_cols().len();
            let decision = self.nodestack.last().unwrap().branchingdecision.clone();
            let s3 = ctx.ms.options_mip_.mip_symmetry_dynamic;
            let pass_stabilizer_to_child_node = if s3 {
                // mip_symmetry_dynamic: SCIP sym_lexorbred.c lines 1091-1148 at 32dd425
                // (no stabilizer orbits under the flag, Q2/G9)
                false
            } else {
                self.orbits_valid_in_child_node(&*ctx.ms, &decision)
            };
            // the sibling reuses the parent's child path
            let s3_sibling_path = if s3 {
                self.nodestack.last().unwrap().s3_child_path.clone()
            } else {
                None
            };
            self.localdom
                .change_bound(ctx.ms, decision, HighsDomain_Reason::branching());
            let currnode = self.nodestack.last().unwrap();
            let mut nodelb = cpp_max(currnode.lower_bound, currnode.other_child_lb);
            let estimate = currnode.estimate;
            let parent_basis = currnode.node_basis.clone();
            let parent_orbits = if pass_stabilizer_to_child_node {
                currnode.stabilizer_orbits.clone()
            } else {
                None
            };

            let mut prune = nodelb > self.get_cutoff_bound(ctx) || self.localdom.infeasible();
            if !prune {
                self.localdom.propagate(ctx.ms);
                prune = self.localdom.infeasible();
                if prune {
                    self.analyze_local_conflict(ctx);
                }
            }
            if !prune {
                if s3 {
                    // mip_symmetry_dynamic: SCIP sym_lexorbred.c lines 1091-1148 at
                    // 32dd425 (sibling child: SCIP depth = current depth)
                    let depth = self.get_current_depth();
                    prune = self.s3_propagate_localdom(ctx.ms, s3_sibling_path.clone(), depth)
                        || self.localdom.infeasible();
                } else {
                    self.propagate_orbitopes(ctx.ms);
                    prune = self.localdom.infeasible();
                }
            }
            if !prune && pass_stabilizer_to_child_node {
                if let Some(orbits) = self.nodestack.last().unwrap().stabilizer_orbits.clone() {
                    self.orbital_fixing(ctx.ms, &orbits);
                    prune = self.localdom.infeasible();
                }
            }
            if prune {
                self.localdom.backtrack(ctx.ms);
                self.localdom.clear_changed_cols_2(num_changed_cols);
                if self.count_tree_weight {
                    let depth = self.get_current_depth();
                    self.treeweight += ldexp_f64(1.0, -depth);
                }
                continue;
            }

            nodelb = cpp_max(nodelb, self.localdom.get_objective_lower_bound());
            if h16_mip_workers {
                // HiGHS 1.16 #3192 (h16_mip_workers): HighsSearch.cpp lines 1715-1725 (1.16.0)
                let node_to_queue = nodelb > ctx.worker().get_optimality_limit(&*ctx.ms);
                if node_to_queue {
                    let queue_estimate = self.nodestack.last().unwrap().estimate;
                    let depth = self.get_current_depth() + 1;
                    // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at
                    // 32dd425 (the queued sibling carries its path, D3)
                    self.stash_node_to_processed(
                        ctx,
                        nodelb,
                        queue_estimate,
                        depth,
                        s3_sibling_path.clone(),
                    );
                    self.localdom.backtrack(ctx.ms);
                    self.localdom.clear_changed_cols_2(num_changed_cols);
                    continue;
                }
            }
            let mut node_to_queue = !h16_mip_workers && nodelb > self.get_optimality_limit(ctx);
            // we check if switching to the other branch of an ancestor yields a higher
            // additive branch score than staying in this node and if so we postpone the
            // node and put it to the queue to backtrack further.
            if !h16_mip_workers && !node_to_queue {
                let mut i: HighsInt = self.nodestack.len() as HighsInt - 2;
                while i >= 0 {
                    let iu = i as usize;
                    if self.nodestack[iu].opensubtrees == 0 {
                        i -= 1;
                        continue;
                    }

                    let ancestor = &self.nodestack[iu];
                    let fallbackbranch =
                        ancestor.branchingdecision.boundval == ancestor.branching_point;
                    let branchpoint = if fallbackbranch {
                        0.5
                    } else {
                        ancestor.branching_point
                    };
                    let ancestor_col = ancestor.branchingdecision.column;
                    let ancestor_score_active: f64;
                    let ancestor_score_inactive: f64;
                    if ancestor.branchingdecision.boundtype == HighsBoundType::kLower {
                        ancestor_score_inactive =
                            ctx.pseudocost.get_score_down(ancestor_col, branchpoint);
                        ancestor_score_active =
                            ctx.pseudocost.get_score_up(ancestor_col, branchpoint);
                    } else {
                        ancestor_score_active =
                            ctx.pseudocost.get_score_down(ancestor_col, branchpoint);
                        ancestor_score_inactive =
                            ctx.pseudocost.get_score_up(ancestor_col, branchpoint);
                    }

                    node_to_queue = ancestor_score_inactive - ancestor_score_active
                        > node_score + self.get_feas_tol(ctx);
                    break;
                }
            }

            if node_to_queue {
                let mut branch_positions: Vec<HighsInt> = Vec::new();
                let domchg_stack = self
                    .localdom
                    .get_reduced_domain_change_stack(&mut branch_positions);
                let queue_estimate = self.nodestack.last().unwrap().estimate;
                let depth = self.get_current_depth() + 1;
                let tmp_tree_weight = if s3 {
                    // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at
                    // 32dd425 (the queued sibling carries its path, D3)
                    ctx.node_queue(nodequeue).emplace_node_s3(
                        domchg_stack,
                        branch_positions,
                        nodelb,
                        queue_estimate,
                        depth,
                        s3_sibling_path.clone(),
                    )
                } else {
                    ctx.node_queue(nodequeue).emplace_node(
                        domchg_stack,
                        branch_positions,
                        nodelb,
                        queue_estimate,
                        depth,
                    )
                };
                if self.count_tree_weight {
                    self.treeweight += tmp_tree_weight;
                }
                self.localdom.backtrack(ctx.ms);
                self.localdom.clear_changed_cols_2(num_changed_cols);
                continue;
            }
            self.nodestack.push(HighsSearch_NodeData::new(
                nodelb,
                estimate,
                parent_basis,
                parent_orbits,
            ));
            if s3 {
                // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at 32dd425
                self.nodestack.last_mut().unwrap().s3_path = s3_sibling_path;
            }

            ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
            self.nodestack.last_mut().unwrap().domgchg_stack_pos = domchg_pos;
            self.eval_origin = crate::work::LP_ORIGIN_SIBLING;
            break;
        }

        if let Some(basis) = self.nodestack.last().unwrap().node_basis.clone() {
            ctx.lp.set_stored_basis(Some(basis));
            ctx.lp.recover_basis();
        }

        true
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::backtrackUntilDepth (lines 1805-1868)
    /// C++: `bool backtrackUntilDepth(HighsInt targetDepth)`
    pub fn backtrack_until_depth(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        target_depth: HighsInt,
    ) -> bool {
        if self.nodestack.is_empty() {
            return false;
        }
        debug_assert!(!self.nodestack.is_empty());
        if self.get_current_depth() >= target_depth {
            self.nodestack.last_mut().unwrap().opensubtrees = 0;
        }

        while self.nodestack.last().unwrap().opensubtrees == 0 {
            self.depthoffset += self.nodestack.last().unwrap().skip_depth_count as HighsInt;
            self.nodestack.pop();

            let branchchg = self.localdom.backtrack(ctx.ms);
            if self.nodestack.is_empty() {
                ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
                return false;
            }
            debug_assert_undone_branching(self.nodestack.last().unwrap(), &branchchg);

            if self.get_current_depth() >= target_depth {
                self.nodestack.last_mut().unwrap().opensubtrees = 0;
            }
        }

        let currnode = self.nodestack.last_mut().unwrap();
        debug_assert!(currnode.opensubtrees == 1);
        currnode.opensubtrees = 0;
        let fallbackbranch = currnode.branchingdecision.boundval == currnode.branching_point;
        flip_branching_decision(currnode, fallbackbranch);

        let domchg_pos = self.localdom.get_domain_change_stack().len() as HighsInt;
        let decision = self.nodestack.last().unwrap().branchingdecision.clone();
        let s3 = ctx.ms.options_mip_.mip_symmetry_dynamic;
        let pass_stabilizer_to_child_node = if s3 {
            // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at 32dd425
            // (no stabilizer orbits under the flag, Q2/G9)
            false
        } else {
            self.orbits_valid_in_child_node(&*ctx.ms, &decision)
        };
        // the sibling reuses the parent's child path
        let s3_sibling_path = if s3 {
            self.nodestack.last().unwrap().s3_child_path.clone()
        } else {
            None
        };
        self.localdom
            .change_bound(ctx.ms, decision, HighsDomain_Reason::branching());
        let currnode = self.nodestack.last().unwrap();
        let lower_bound = currnode.lower_bound;
        let estimate = currnode.estimate;
        let parent_basis = currnode.node_basis.clone();
        let parent_orbits = if pass_stabilizer_to_child_node {
            currnode.stabilizer_orbits.clone()
        } else {
            None
        };
        self.nodestack.push(HighsSearch_NodeData::new(
            lower_bound,
            estimate,
            parent_basis,
            parent_orbits,
        ));
        if s3 {
            // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at 32dd425
            self.nodestack.last_mut().unwrap().s3_path = s3_sibling_path;
        }

        ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
        self.nodestack.last_mut().unwrap().domgchg_stack_pos = domchg_pos;
        if let Some(basis) = self.nodestack.last().unwrap().node_basis.clone() {
            if basis.row_status.len() as HighsInt == ctx.lp.get_lp().num_row_ {
                ctx.lp.set_stored_basis(Some(basis));
            }
        }
        ctx.lp.recover_basis();

        true
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::dive (lines 1870-1885)
    /// C++: `NodeResult dive(int64_t nodeLim = std :: numeric_limits < int64_t > :: max ( ))`
    /// C++ default arguments: `node_lim = std :: numeric_limits < int64_t > :: max ( )`
    pub fn dive(&mut self, ctx: &mut SearchCtx<'_>, node_lim: i64) -> HighsSearch_NodeResult {
        self.reliableatnode.clear();

        loop {
            self.nnodes += 1;
            let result = self.evaluate_node(ctx);

            if self.check_limits(ctx, self.nnodes) {
                return result;
            }

            if result != HighsSearch_NodeResult::kOpen {
                return result;
            }

            let result = self.branch(ctx);
            if result != HighsSearch_NodeResult::kBranched {
                return result;
            }
            if self.nnodes >= node_lim {
                return result;
            }
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::solveDepthFirst (lines 1887-1898)
    /// HiGHS latest (#3445, b1282872a): HighsSearch.cpp::HighsSearch::solveDepthFirst (lines 1809-1822 at b1282872a)
    /// C++: `void solveDepthFirst(int64_t maxbacktracks = 1, int64_t maxnodes = std::numeric_limits<int64_t>::max())`
    /// C++ default arguments: `maxbacktracks = 1`, `maxnodes = i64::MAX`
    pub fn solve_depth_first(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        maxbacktracks: i64,
        maxnodes: i64,
    ) {
        let mut maxbacktracks = maxbacktracks;
        loop {
            if maxbacktracks == 0 {
                break;
            }

            let result = self.dive(ctx, maxnodes);
            // if a limit was reached the result might be open
            if result == HighsSearch_NodeResult::kOpen {
                break;
            }
            if self.nnodes >= maxnodes {
                break;
            }

            maxbacktracks -= 1;

            if !self.backtrack(ctx, true) {
                break;
            }
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getFeasTol (lines 1900-1900)
    /// C++: `double getFeasTol() const`
    pub fn get_feas_tol(&self, ctx: &SearchCtx<'_>) -> f64 {
        ctx.ms.mipdata_.feastol
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getUpperLimit (lines 1902-1908)
    /// C++: `double getUpperLimit() const`
    pub fn get_upper_limit(&self, ctx: &SearchCtx<'_>) -> f64 {
        if !ctx.ms.parallel_lock_active() {
            ctx.ms.mipdata_.upper_limit
        } else {
            ctx.worker().upper_limit
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getEpsilon (lines 1910-1910)
    /// C++: `double getEpsilon() const`
    pub fn get_epsilon(&self, ctx: &SearchCtx<'_>) -> f64 {
        ctx.ms.mipdata_.epsilon
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getOptimalityLimit (lines 1912-1918)
    /// C++: `double getOptimalityLimit() const`
    pub fn get_optimality_limit(&self, ctx: &SearchCtx<'_>) -> f64 {
        if !ctx.ms.parallel_lock_active() {
            ctx.ms.mipdata_.optimality_limit
        } else {
            ctx.worker().optimality_limit
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getRootLpSol (lines 1920-1922)
    /// C++: `const std::vector<double> & getRootLpSol() const`
    pub fn get_root_lp_sol<'c>(&self, ctx: &'c SearchCtx<'_>) -> &'c [f64] {
        &ctx.ms.mipdata_.rootlpsol
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getIntegralCols (lines 1924-1926)
    /// C++: `const std::vector<HighsInt> & getIntegralCols() const`
    pub fn get_integral_cols<'c>(&self, ctx: &'c SearchCtx<'_>) -> &'c [HighsInt] {
        &ctx.ms.mipdata_.integral_cols
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getDomain (lines 1928-1930)
    /// C++: `HighsDomain & getDomain() const`
    // PORT NOTE: returns the worker's global-domain index (`ms.mipdata_.domains[i]`); read it
    // attached, mutate it with `detach!`.
    pub fn get_domain(&self, ctx: &SearchCtx<'_>) -> usize {
        ctx.worker().globaldom_
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getConflictPool (lines 1932-1934)
    /// C++: `HighsConflictPool & getConflictPool() const`
    // PORT NOTE: index into `ms.mipdata_.conflictpools` (the `conflict_pool` parameter of
    // `HighsDomain::conflict_analysis*`).
    pub fn get_conflict_pool(&self, ctx: &SearchCtx<'_>) -> usize {
        ctx.worker().conflictpool_
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getCutPool (lines 1936-1936)
    /// C++: `HighsCutPool & getCutPool() const`
    // PORT NOTE: index into `ms.mipdata_.cutpools` (`HighsCutGeneration::new(ms, lp, cutpool)`).
    pub fn get_cut_pool(&self, ctx: &SearchCtx<'_>) -> usize {
        ctx.worker().cutpool_
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getNodeQueue (lines 1938-1940)
    /// C++: `const HighsNodeQueue & getNodeQueue() const`
    pub fn get_node_queue<'c>(&self, ctx: &'c SearchCtx<'_>) -> &'c HighsNodeQueue {
        &ctx.ms.mipdata_.nodequeue
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::checkLimits (lines 1942-1947)
    /// C++: `bool checkLimits(int64_t nodeOffset = 0) const`
    /// C++ default arguments: `node_offset = 0`
    // PORT NOTE: `ctx` is `&mut`: `HighsMipSolverData::checkLimits` writes `modelstatus_` and
    // runs the interrupt callback.
    pub fn check_limits(&self, ctx: &mut SearchCtx<'_>, node_offset: i64) -> bool {
        if ctx.ms.parallel_lock_active() {
            return self.check_local_limits(ctx);
        };
        ctx.ms.check_limits(node_offset)
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::checkLocalLimits (lines 1949-1982)
    /// HiGHS 1.16 (#3124): HighsSearch.cpp::checkLocalLimits (lines 1867-1906), under `h16_mip_workers`
    /// C++: `bool checkLocalLimits() const`
    pub fn check_local_limits(&self, ctx: &SearchCtx<'_>) -> bool {
        self.check_local_limits_w(&*ctx.ms, ctx.worker)
    }

    /// [`Self::check_local_limits`] for the search of worker `w`, callable without a
    /// `SearchCtx` (the 1.16 `subMipShouldRetry` calls `worker.search_ptr_->checkLimits()`
    /// from inside RENS/RINS, #3237).
    pub fn check_local_limits_w(&self, ms: &HighsMipSolver, w: usize) -> bool {
        if ms.terminator_active() && ms.terminator_terminated() {
            return true;
        }

        // HiGHS 1.16 #3124 (h16_mip_workers): HighsSearch.cpp checkLocalLimits lines 1871-1875 (1.16.0)
        if ms.options_mip_.h16_mip_workers {
            let stop = ms
                .mipdata_
                .worker_lp_iterations_stop
                .0
                .load(std::sync::atomic::Ordering::Relaxed);
            if stop <= self.lpiterations {
                return true;
            }
        }

        let upper_bound = ms.mipdata_.workers[w].upper_bound;
        if !ms.submip && upper_bound < kHighsInf && ms.options_mip_.objective_target > -kHighsInf {
            let internal_target = (ms.orig_model_.sense_ as HighsInt) as f64
                * ms.options_mip_.objective_target
                - mip_model!(ms).offset_;
            if upper_bound < internal_target {
                return true;
            }
        }

        if ms.options_mip_.mip_max_nodes != kHighsIInf
            && ms.mipdata_.num_nodes + self.nnodes >= ms.options_mip_.mip_max_nodes as i64
        {
            return true;
        }

        if ms.options_mip_.mip_max_leaves != kHighsIInf
            && ms.mipdata_.num_leaves + self.nleaves >= ms.options_mip_.mip_max_leaves as i64
        {
            return true;
        }

        if ms.options_mip_.time_limit < kHighsInf && ms.timer_.read() >= ms.options_mip_.time_limit
        {
            return true;
        }

        false
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getSymmetries (lines 1984-1986)
    /// C++: `HighsSymmetries & getSymmetries() const`
    // PORT NOTE: read-only access (`columnPosition`, `numPerms`). Calls that also take `ms`
    // (`propagateOrbitopes`, `orbitalFixing`, `computeStabilizerOrbits`) go through
    // `ms.with_symmetries_ref(..)` / `HighsSearch::propagate_orbitopes` / `orbital_fixing`.
    pub fn get_symmetries<'c>(&self, ctx: &'c SearchCtx<'_>) -> &'c HighsSymmetries {
        &ctx.ms.mipdata_.symmetries
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::addIncumbent (lines 1988-1997)
    /// C++: `bool addIncumbent(const std::vector<double> & sol, double solobj, const int solution_source, const bool print_display_line = true)`
    /// C++ default arguments: `print_display_line = true`
    // PORT NOTE: `sol` cannot borrow from `ctx.lp` (pass a clone of the LP solution). The
    // serial path lends `ctx.lp` / `ctx.pseudocost` back (`HighsMipSolverData::addIncumbent`
    // reads every LP); the worker stays attached, so the C++ loop over all workers sees it.
    pub fn add_incumbent(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        sol: &[f64],
        solobj: f64,
        solution_source: i32,
        print_display_line: bool,
    ) -> bool {
        if ctx.ms.parallel_lock_active() {
            let w = ctx.worker;
            let ms = &mut *ctx.ms;
            detach!(ms.mipdata_.workers[w], |wk| wk.add_incumbent(
                ms,
                sol,
                solobj,
                solution_source
            ))
        } else {
            ctx.lend_back(|ms| {
                ms.add_incumbent(sol, solobj, solution_source, print_display_line, false)
            })
        }
    }
}

#[cfg(test)]
mod tests {
    //! HiGHS backport-SEARCH tests (#3213, #3445).

    use super::*;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    const N_ROWS: usize = 3;

    /// A strongly correlated multi-row binary knapsack (value = mean weight + 10), so the
    /// root LP bound does not close and the tree branches: min -value.x,
    /// sum_j w_ij x_j <= half the row weight, x binary. Weights come from a fixed LCG.
    fn knapsack_lp(n: usize) -> HighsLp {
        let mut state: u64 = 12345;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (10 + (state >> 33) % 91) as f64
        };
        let weights: Vec<[f64; N_ROWS]> = (0..n).map(|_| [next(), next(), next()]).collect();
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = N_ROWS as HighsInt;
        lp.col_cost_ = weights
            .iter()
            .map(|w| -((w[0] + w[1] + w[2]) / 3.0).floor() - 10.0)
            .collect();
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf; N_ROWS];
        lp.row_upper_ = (0..N_ROWS)
            .map(|i| (weights.iter().map(|w| w[i]).sum::<f64>() / 2.0).floor())
            .collect();
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = N_ROWS as HighsInt;
        lp.a_matrix_.start_ = (0..=n).map(|j| (N_ROWS * j) as HighsInt).collect();
        lp.a_matrix_.index_ = (0..n).flat_map(|_| [0, 1, 2]).collect();
        lp.a_matrix_.value_ = weights.iter().flat_map(|w| *w).collect();
        lp
    }

    /// The knapsack optimum by enumeration (for the small instance).
    fn brute_force_optimum(lp: &HighsLp) -> f64 {
        let n = lp.num_col_ as usize;
        let mut best = 0.0f64;
        for mask in 0u32..(1 << n) {
            let mut act = [0.0f64; N_ROWS];
            let mut v = 0.0;
            for j in (0..n).filter(|j| mask & (1 << j) != 0) {
                for (i, a) in act.iter_mut().enumerate() {
                    *a += lp.a_matrix_.value_[N_ROWS * j + i];
                }
                v -= lp.col_cost_[j];
            }
            if act.iter().zip(&lp.row_upper_).all(|(a, u)| a <= u) && v > best {
                best = v;
            }
        }
        -best
    }

    fn base_options() -> HighsOptions {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.presolve = kHighsOffString.to_string();
        // No primal heuristics, so the root does not close the gap and the tree branches.
        options.mip_heuristic_effort = 0.0;
        options.mip_heuristic_run_feasibility_jump = false;
        options.mip_heuristic_run_rens = false;
        options.mip_heuristic_run_rins = false;
        options.mip_heuristic_run_root_reduced_cost = false;
        options.mip_heuristic_run_zi_round = false;
        options.mip_heuristic_run_shifting = false;
        options
    }

    /// #3213 (HighsSearch.cpp:1118 at 4a638fa93): with the option on, `branch` divides the
    /// strong-branching budget by four while `numRestarts <= 2`. On this model the budget
    /// (at least 25000 iterations) never binds, so the solve must reach `branch` and give
    /// the same optimum and tree as with the option off.
    #[test]
    fn x1_search_sb_budget_shift() {
        let solve = |on: bool| {
            std::thread::spawn(move || {
                initialize_scheduler(1);
                let mut options = base_options();
                options.mip_sb_budget_pre_restart = on;
                options.mip_allow_restart = false;
                let lp = knapsack_lp(20);
                let solution = HighsSolution::default();
                let mut callback = HighsCallback::default();
                let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
                ms.run();
                let out = (
                    ms.modelstatus_,
                    ms.solution_objective_.to_bits(),
                    ms.node_count_,
                    ms.mipdata_.num_restarts,
                );
                drop(ms);
                HighsTaskExecutor::shutdown(true);
                out
            })
            .join()
            .unwrap()
        };
        let off = solve(false);
        let on = solve(true);
        assert_eq!(on.0, HighsModelStatus::kOptimal);
        assert_eq!(f64::from_bits(on.1), brute_force_optimum(&knapsack_lp(20)));
        // `branch` (and so the gate) is reached only when the tree has more than the root.
        assert!(on.2 > 1, "node count {}", on.2);
        assert!(on.3 <= 2);
        assert_eq!(on, off);
    }

    /// #3445 (HighsSearch.cpp:1809-1820 at b1282872a): `solveDepthFirst(maxbacktracks,
    /// maxnodes)` passes `maxnodes` to `dive` and stops once `nnodes >= maxnodes`. The
    /// search is set up as `HighsPrimalHeuristics::RENS` sets up `heur` on the
    /// root of the knapsack after a one-node solve.
    #[test]
    fn x1_search_dfs_limit_stops() {
        let run = |maxnodes: i64| {
            std::thread::spawn(move || {
                initialize_scheduler(1);
                let mut options = base_options();
                options.mip_max_nodes = 1;
                let lp = knapsack_lp(40);
                let solution = HighsSolution::default();
                let mut callback = HighsCallback::default();
                let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
                ms.run();
                ms.options_mip_.mip_max_nodes = kHighsIInf;
                // Search the whole root again: no incumbent prunes the depth-first search.
                ms.mipdata_.upper_limit = kHighsInf;
                ms.mipdata_.workers[0].upper_limit = kHighsInf;

                let worker = 0usize;
                let gd = ms.mipdata_.workers[worker].get_global_domain();
                let pscost_k = ms.mipdata_.workers[worker].get_pseudocost();
                let mut pscost: HighsPseudocost =
                    detach!(ms.mipdata_.pseudocosts[pscost_k], |p| p.clone());
                let mut heur = HighsSearch::new(&mut ms, gd);
                heur.set_heuristic(true);
                let lp_k = ms.mipdata_.workers[worker].get_lp_relaxation();
                let mut heurlp = HighsLpRelaxation::copy_from(&ms, &ms.mipdata_.lps[lp_k]);
                heurlp.set_mip_worker(worker);
                heurlp.set_objective_limit(&ms, kHighsInf);
                heurlp.set_adjust_symmetric_branching_col(false);
                let num_col_minus_one = ms.num_col() - 1;
                let _ = heurlp.get_lp_solver().change_cols_bounds_interval(
                    0,
                    num_col_minus_one,
                    &heur.localdom.col_lower_,
                    &heur.localdom.col_upper_,
                );
                heur.localdom.clear_changed_cols();
                heur.create_new_node();
                {
                    let mut ctx = SearchCtx::new(&mut ms, worker, &mut heurlp, &mut pscost);
                    heur.set_min_reliable(&mut ctx, 0);
                    heur.solve_depth_first(&mut ctx, 10, maxnodes);
                }
                let out = (heur.nnodes, heur.has_node());
                heur.release(&mut ms);
                drop(ms);
                HighsTaskExecutor::shutdown(true);
                out
            })
            .join()
            .unwrap()
        };
        let (unlimited, _) = run(i64::MAX);
        assert!(
            unlimited > 3,
            "unlimited depth-first search took {unlimited} nodes"
        );
        for maxnodes in [1i64, 2, 3] {
            let (nnodes, open) = run(maxnodes);
            assert_eq!(nnodes, maxnodes);
            // The search stopped on a limit, not because the tree was exhausted.
            assert!(open);
        }
    }

    /// Knapsack whose columns come in identical pairs (j and j + n), so symmetry
    /// detection finds the column swaps and, with the flag on, data is built.
    fn paired_knapsack_lp(n: usize) -> HighsLp {
        let base = knapsack_lp(n);
        let mut lp = base.clone();
        let m = 2 * n;
        lp.num_col_ = m as HighsInt;
        lp.col_cost_ = [base.col_cost_.clone(), base.col_cost_.clone()].concat();
        lp.col_lower_ = vec![0.0; m];
        lp.col_upper_ = vec![1.0; m];
        lp.integrality_ = vec![HighsVarType::kInteger; m];
        lp.row_upper_ = base.row_upper_.iter().map(|u| u * 2.0).collect();
        lp.a_matrix_.num_col_ = m as HighsInt;
        lp.a_matrix_.start_ = (0..=m).map(|j| (N_ROWS * j) as HighsInt).collect();
        lp.a_matrix_.index_ = (0..m).flat_map(|_| [0, 1, 2]).collect();
        lp.a_matrix_.value_ =
            [base.a_matrix_.value_.clone(), base.a_matrix_.value_.clone()].concat();
        lp
    }

    fn s3_gate2_solve(lp: HighsLp, on: bool) -> (HighsModelStatus, u64, i64, bool) {
        std::thread::spawn(move || {
            initialize_scheduler(1);
            let mut options = base_options();
            options.mip_symmetry_dynamic = on;
            options.mip_allow_restart = false;
            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let out = (
                ms.modelstatus_,
                ms.solution_objective_.to_bits(),
                ms.node_count_,
                ms.mipdata_.symmetry_dynamic.is_some(),
            );
            drop(ms);
            HighsTaskExecutor::shutdown(true);
            out
        })
        .join()
        .unwrap()
    }

    /// without symmetry no data exists, so the gated branch(), backtrack and
    /// backtrack_plunge paths with the flag on give the same optimum and tree as flag off.
    #[test]
    fn s3_gate2_no_symmetry_matches_off() {
        let off = s3_gate2_solve(knapsack_lp(20), false);
        let on = s3_gate2_solve(knapsack_lp(20), true);
        assert_eq!(off.0, HighsModelStatus::kOptimal);
        assert_eq!(f64::from_bits(off.1), brute_force_optimum(&knapsack_lp(20)));
        assert!(off.2 > 1, "node count {}", off.2);
        assert!(!on.3);
        assert_eq!(on, off);
    }

    /// on a model with symmetric column pairs the flag-on search keeps the
    /// brute-force optimum (flag off is the 1.15.1 reference).
    #[test]
    fn s3_gate2_symmetric_keeps_optimum() {
        let lp = paired_knapsack_lp(8);
        let opt = brute_force_optimum(&lp);
        let off = s3_gate2_solve(lp.clone(), false);
        let on = s3_gate2_solve(lp, true);
        assert_eq!(off.0, HighsModelStatus::kOptimal);
        assert_eq!(on.0, HighsModelStatus::kOptimal);
        assert!(
            (f64::from_bits(off.1) - opt).abs() <= 1e-6,
            "off {}",
            f64::from_bits(off.1)
        );
        assert!(
            (f64::from_bits(on.1) - opt).abs() <= 1e-6,
            "on {}",
            f64::from_bits(on.1)
        );
        assert!(!off.3);
        eprintln!(
            "s3_gate2 symmetric: off nodes {} on nodes {} s3 data {}",
            off.2, on.2, on.3
        );
    }
}
