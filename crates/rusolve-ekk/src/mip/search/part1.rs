//! HiGHS source: mip/HighsSearch.cpp (lines 1-1098)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::mip::symmetry_dynamic::hook::{s3_child_path, s3_propagate};
use crate::prelude::*;
use std::sync::Arc;

/// rusolve `mip_child_selection`: the option value (0..8) as HiGHS's child selection rule
/// (`HighsSearch::ChildSelectionRule`, HiGHS 1.16 `mip/HighsSearch.h`). HiGHS itself uses
/// kRootSol in the main search and kHybridInferenceCost in sub-MIPs; this lets the main
/// search use any of the rules HiGHS already implements. -1 (default) keeps HiGHS.
fn child_selection_rule(v: HighsInt) -> HighsSearch_ChildSelectionRule {
    use HighsSearch_ChildSelectionRule::*;
    match v {
        0 => kUp,
        1 => kDown,
        2 => kRootSol,
        3 => kObj,
        4 => kRandom,
        5 => kBestCost,
        6 => kWorstCost,
        7 => kDisjunction,
        _ => kHybridInferenceCost,
    }
}

impl HighsSearch {
    /// rusolve `mip_child_selection_switch_gap`: the child rule for this branching. With
    /// `mip_child_selection` set (main search only), the chosen rule (for example the hybrid
    /// inference/cost rule, which finds incumbents sooner) is used while the relative primal-dual
    /// gap of `limitsToGap` is above the threshold; at or below it the search goes back to the
    /// HiGHS root-solution rule, which closes trees faster once a good incumbent is known. This
    /// split is our own; the evidence is in the pull request. 0 = never switch (the rule stays).
    pub(crate) fn effective_child_rule(
        &self,
        ms: &HighsMipSolver,
    ) -> HighsSearch_ChildSelectionRule {
        let sw = ms.options_mip_.mip_child_selection_switch_gap;
        if sw > 0.0 && !ms.submip && !self.inheuristic && ms.options_mip_.mip_child_selection >= 0 {
            let (mut lb, mut ub) = (0.0, 0.0);
            let gap = ms.limits_to_gap(
                ms.mipdata_.lower_bound,
                ms.mipdata_.upper_bound,
                &mut lb,
                &mut ub,
            );
            if gap <= sw {
                return HighsSearch_ChildSelectionRule::kRootSol;
            }
        }
        self.childselrule
    }

    /// rusolve `mip_plunge_quot`: the largest node lower bound the search still plunges into, or
    /// +inf when the rule is off. Plunging (diving into a child, or backtracking to a sibling on
    /// the node stack) continues only while the node's lower bound is at most
    /// `lb + quot * (cutoff - lb)`, with `lb` the global lower bound and `cutoff` the upper limit;
    /// a node above it goes to the queue, where best-bound and best-estimate picks take over.
    /// This is the `maxplungequot` rule of SCIP's best-first and hybrid best-estimate node
    /// selectors (SCIP, Apache-2.0: `src/scip/nodesel_bfs.c`, `src/scip/nodesel_hybridestim.c`;
    /// Achterberg, "Constraint Integer Programming", PhD thesis, TU Berlin 2007, ch. 6), applied
    /// only when [`plunge_levers_active`] says so.
    pub(crate) fn plunge_max_bound(&self, ms: &HighsMipSolver) -> f64 {
        let quot = ms.options_mip_.mip_plunge_quot;
        if quot <= 0.0 || !plunge_levers_active(ms, self.inheuristic) {
            return kHighsInf;
        }
        let lb = ms.mipdata_.lower_bound;
        lb + quot * (ms.mipdata_.upper_limit - lb)
    }
}

/// rusolve `mip_plunge_gap`: whether the plunge levers (`mip_plunge_quot`,
/// `mip_plunge_max_nodes`, `mip_plunge_bestbound_freq`) act now. They act in the main search
/// (not in sub-MIPs or heuristic dives) once the upper limit is finite and the relative gap between the global lower bound and the upper
/// limit (the incumbent cutoff, or `objective_bound`) is at or below `mip_plunge_gap`. Before
/// that, deep plunging finds incumbents; after it, picking low-bound nodes more often raises the
/// reported dual bound. 0 = off: the HiGHS plunge.
pub(crate) fn plunge_levers_active(ms: &HighsMipSolver, inheuristic: bool) -> bool {
    let max_gap = ms.options_mip_.mip_plunge_gap;
    if max_gap <= 0.0 || ms.submip || inheuristic {
        return false;
    }
    let ul = ms.mipdata_.upper_limit;
    if ul == kHighsInf {
        return false;
    }
    let (mut lb, mut ub) = (0.0, 0.0);
    ms.limits_to_gap(ms.mipdata_.lower_bound, ul, &mut lb, &mut ub) <= max_gap
}

/// C++ `std::pair<A, B>` `operator<` (libc++): `p.first < q.first || (!(q.first < p.first) && p.second < q.second)`.
/// `>` and `>=` are written as `std_pair_lt` with the operands swapped or negated.
fn std_pair_lt<A: PartialOrd, B: PartialOrd>(p: (A, B), q: (A, B)) -> bool {
    p.0 < q.0 || (!(q.0 < p.0) && p.1 < q.1)
}

/// The C++ locals of `selectBranchingCandidate` that its lambdas share.
struct SbState {
    upscore: Vec<f64>,
    downscore: Vec<f64>,
    upscorereliable: Vec<bool>,
    downscorereliable: Vec<bool>,
    upbound: Vec<f64>,
    downbound: Vec<f64>,
    evalqueue: Vec<HighsInt>,
    min_score: f64,
}

/// C++ lambda `selectBestScore(bool finalSelection)` of `selectBranchingCandidate`.
#[allow(clippy::needless_late_init)]
fn sb_select_best_score(
    st: &mut SbState,
    fracints: &[(HighsInt, f64)],
    nodequeue: &HighsNodeQueue,
    pseudocost: &HighsPseudocost,
    feastol: f64,
    final_selection: bool,
) -> HighsInt {
    let mut best: HighsInt = -1;
    let mut bestscore: f64 = -1.0;
    let mut bestnodes: f64 = -1.0;
    let mut bestnumnodes: i64 = 0;

    let oldminscore = st.min_score;
    let evalqueue = std::mem::take(&mut st.evalqueue);
    for &k in evalqueue.iter() {
        let k = k as usize;
        let col = fracints[k].0;
        let fracval = fracints[k].1;

        if st.upscore[k] <= oldminscore {
            st.upscorereliable[k] = true;
        }
        if st.downscore[k] <= oldminscore {
            st.downscorereliable[k] = true;
        }

        let s = 1e-3
            * cpp_min(
                if st.upscorereliable[k] {
                    st.upscore[k]
                } else {
                    0.0
                },
                if st.downscorereliable[k] {
                    st.downscore[k]
                } else {
                    0.0
                },
            );
        st.min_score = cpp_max(s, st.min_score);

        let score: f64;
        if st.upscore[k] <= oldminscore || st.downscore[k] <= oldminscore {
            score = pseudocost.get_score(
                col,
                cpp_min(st.upscore[k], oldminscore),
                cpp_min(st.downscore[k], oldminscore),
            );
        } else if st.upscore[k] == kHighsInf || st.downscore[k] == kHighsInf {
            score = if final_selection {
                pseudocost.get_score_2(col, fracval)
            } else {
                kHighsInf
            };
        } else {
            score = pseudocost.get_score(col, st.upscore[k], st.downscore[k]);
        }

        debug_assert!(score >= 0.0);
        let upnodes = nodequeue.num_nodes_up(col);
        let downnodes = nodequeue.num_nodes_down(col);
        let mut nodes: f64 = 0.0;
        let numnodes: i64 = upnodes + downnodes;
        if upnodes != 0 || downnodes != 0 {
            nodes = (downnodes as f64 / numnodes as f64) * (upnodes as f64 / numnodes as f64);
        }
        if score > bestscore
            || (score > bestscore - feastol
                && std_pair_lt((bestnodes, bestnumnodes), (nodes, numnodes)))
        {
            bestscore = score;
            best = k as HighsInt;
            bestnodes = nodes;
            bestnumnodes = numnodes;
        }
    }
    st.evalqueue = evalqueue;

    best
}

impl HighsSearch {
    /// C++ lambda `analyzeSolution(objdelta, sol)` of `selectBranchingCandidate`.
    #[allow(clippy::needless_range_loop)]
    fn sb_analyze_solution(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        st: &mut SbState,
        fracints: &[(HighsInt, f64)],
        col: HighsInt,
        objdelta: f64,
        sol: &[f64],
    ) {
        let feastol = self.get_feas_tol(ctx);
        let pool = self.get_conflict_pool(ctx);
        let globaldom = ctx.worker().get_global_domain();
        let numchanged = self.localdom.get_changed_cols().len();
        let domchg_stack_size = self.localdom.get_domain_change_stack().len() as HighsInt;
        let numfrac = fracints.len();

        for k in 0..numfrac {
            if fracints[k].0 == col {
                continue;
            }
            let fi = fracints[k].0;
            let fi_u = fi as usize;
            let otherfracval = fracints[k].1;
            let otherdownval = fracints[k].1.floor();
            let otherupval = fracints[k].1.ceil();
            if sol[fi_u] <= otherdownval + feastol {
                if self.localdom.col_upper_[fi_u] > otherdownval + feastol {
                    self.localdom.change_bound_2(
                        ctx.ms,
                        HighsBoundType::kUpper,
                        fi,
                        otherdownval,
                        HighsDomain_Reason::branching(),
                    );
                    if self.localdom.infeasible() {
                        self.localdom
                            .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
                        self.localdom.backtrack(ctx.ms);
                        self.localdom.clear_changed_cols_2(numchanged);
                        continue;
                    }
                    self.localdom.propagate(ctx.ms);
                    if self.localdom.infeasible() {
                        self.localdom
                            .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
                        self.localdom.backtrack(ctx.ms);
                        self.localdom.clear_changed_cols_2(numchanged);
                        continue;
                    }

                    let new_stack_size = self.localdom.get_domain_change_stack().len() as HighsInt;
                    let mut solution_valid = true;
                    for j in (domchg_stack_size + 1)..new_stack_size {
                        let chg = &self.localdom.get_domain_change_stack()[j as usize];
                        if chg.boundtype == HighsBoundType::kLower {
                            if chg.boundval > sol[chg.column as usize] + feastol {
                                solution_valid = false;
                                break;
                            }
                        } else if chg.boundval < sol[chg.column as usize] - feastol {
                            solution_valid = false;
                            break;
                        }
                    }

                    self.localdom.backtrack(ctx.ms);
                    self.localdom.clear_changed_cols_2(numchanged);
                    if !solution_valid {
                        continue;
                    }
                }

                if objdelta <= feastol {
                    ctx.pseudocost
                        .add_observation(fi, otherdownval - otherfracval, objdelta);
                    self.mark_branching_var_down_reliable_at_node(fi);
                }

                st.downscore[k] = cpp_min(st.downscore[k], objdelta);
            } else if sol[fi_u] >= otherupval - feastol {
                if self.localdom.col_lower_[fi_u] < otherupval - feastol {
                    self.localdom.change_bound_2(
                        ctx.ms,
                        HighsBoundType::kLower,
                        fi,
                        otherupval,
                        HighsDomain_Reason::branching(),
                    );

                    if self.localdom.infeasible() {
                        self.localdom
                            .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
                        self.localdom.backtrack(ctx.ms);
                        self.localdom.clear_changed_cols_2(numchanged);
                        continue;
                    }
                    self.localdom.propagate(ctx.ms);
                    if self.localdom.infeasible() {
                        self.localdom
                            .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
                        self.localdom.backtrack(ctx.ms);
                        self.localdom.clear_changed_cols_2(numchanged);
                        continue;
                    }

                    let new_stack_size = self.localdom.get_domain_change_stack().len() as HighsInt;
                    let mut solution_valid = true;
                    for j in (domchg_stack_size + 1)..new_stack_size {
                        let chg = &self.localdom.get_domain_change_stack()[j as usize];
                        if chg.boundtype == HighsBoundType::kLower {
                            if chg.boundval > sol[chg.column as usize] + feastol {
                                solution_valid = false;
                                break;
                            }
                        } else if chg.boundval < sol[chg.column as usize] - feastol {
                            solution_valid = false;
                            break;
                        }
                    }

                    self.localdom.backtrack(ctx.ms);
                    self.localdom.clear_changed_cols_2(numchanged);

                    if !solution_valid {
                        continue;
                    }
                }

                if objdelta <= feastol {
                    ctx.pseudocost
                        .add_observation(fi, otherupval - otherfracval, objdelta);
                    self.mark_branching_var_up_reliable_at_node(fi);
                }

                st.upscore[k] = cpp_min(st.upscore[k], objdelta);
            }
        }
    }

    /// C++ lambda `strongBranch(bool upbranch)` of `selectBranchingCandidate`. `true` is the
    /// C++ `return true` (the node was branched on); `false` is the C++ `return false`.
    fn sb_strong_branch(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        playground: &mut HighsLpRelaxation_Playground,
        st: &mut SbState,
        fracints: &[(HighsInt, f64)],
        candidate: usize,
        upbranch: bool,
    ) -> bool {
        let ci = candidate;
        let col = fracints[ci].0;
        let fracval = fracints[ci].1;
        let upval = fracval.ceil();
        let downval = fracval.floor();

        let mut inferences: HighsInt =
            -(self.localdom.get_domain_change_stack().len() as HighsInt) - 1;
        let boundtype = if upbranch {
            HighsBoundType::kLower
        } else {
            HighsBoundType::kUpper
        };
        let boundval = if upbranch { upval } else { downval };
        let domchg = HighsDomainChange {
            boundval,
            column: col,
            boundtype,
        };

        // mip_symmetry_dynamic: Q6, Q2 — no stabilizer or orbitope calls in strong branching
        let flag = ctx.ms.options_mip_.mip_symmetry_dynamic;
        let orbital_fixing = if flag {
            false
        } else {
            self.nodestack.last().unwrap().stabilizer_orbits.is_some()
                && self.orbits_valid_in_child_node(&*ctx.ms, &domchg)
        };
        self.localdom
            .change_bound(ctx.ms, domchg, HighsDomain_Reason::branching());
        self.localdom.propagate(ctx.ms);

        if !self.localdom.infeasible() {
            if flag {
                // no stabilizer or orbitope calls in strong branching
            } else if orbital_fixing {
                let orbits = self
                    .nodestack
                    .last()
                    .unwrap()
                    .stabilizer_orbits
                    .clone()
                    .unwrap();
                self.orbital_fixing(ctx.ms, &orbits);
            } else {
                self.propagate_orbitopes(ctx.ms);
            }
        }

        inferences += self.localdom.get_domain_change_stack().len() as HighsInt;
        if self.localdom.infeasible() {
            let pool = self.get_conflict_pool(ctx);
            let globaldom = ctx.worker().get_global_domain();
            self.localdom
                .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
            ctx.pseudocost.add_cutoff_observation(col, upbranch);
            self.localdom.backtrack(ctx.ms);
            self.localdom.clear_changed_cols();

            if upbranch {
                self.branch_downwards(ctx, col, downval, fracval);
            } else {
                self.branch_upwards(ctx, col, upval, fracval);
            }
            let len = self.nodestack.len();
            self.nodestack[len - 2].opensubtrees = 0;
            self.nodestack[len - 2].skip_depth_count = 1;
            self.depthoffset -= 1;

            return true;
        }

        ctx.pseudocost
            .add_inference_observation(col, inferences, upbranch);

        let mut numiters = ctx.lp.get_num_lp_iterations();
        let lp_src_guard = crate::work::lp_kind(crate::work::LP_STRONG_BRANCH);
        let status = playground.solve_lp(ctx.ms, ctx.lp, &mut self.localdom);
        drop(lp_src_guard);
        numiters = ctx.lp.get_num_lp_iterations() - numiters;
        self.lpiterations += numiters;
        self.sblpiterations += numiters;

        if HighsLpRelaxation::scaled_optimal(status) {
            ctx.lp.perform_aging(ctx.ms, false);

            let delta = if upbranch {
                upval - fracval
            } else {
                downval - fracval
            };
            // HiGHS: `sol` is a `const&` to the LP solution; cloned here so the LP can be lent
            // to `addIncumbent` below. The LP solution is not changed by the calls in between.
            let sol: Vec<f64> = ctx.lp.get_solution().col_value.clone();
            let mut integerfeasible = false;
            let solobj = self.check_sol(ctx, &sol, &mut integerfeasible);

            let mut objdelta = cpp_max(solobj - ctx.lp.get_objective(), 0.0);
            if objdelta <= self.get_epsilon(ctx) {
                objdelta = 0.0;
            }

            if upbranch {
                st.upscore[ci] = objdelta;
                st.upscorereliable[ci] = true;
                self.mark_branching_var_up_reliable_at_node(col);
            } else {
                st.downscore[ci] = objdelta;
                st.downscorereliable[ci] = true;
                self.mark_branching_var_down_reliable_at_node(col);
            }
            ctx.pseudocost.add_observation(col, delta, objdelta);
            self.sb_analyze_solution(ctx, st, fracints, col, objdelta, &sol);

            if HighsLpRelaxation::unscaled_primal_feasible(status) && integerfeasible {
                let cutoffbnd = self.get_cutoff_bound(ctx);
                let source = if self.inheuristic {
                    kSolutionSourceHeuristic
                } else {
                    kSolutionSourceBranching
                };
                self.add_incumbent(ctx, &sol, solobj, source, true);

                if self.get_upper_limit(ctx) < cutoffbnd {
                    let ul = self.get_upper_limit(ctx);
                    ctx.lp.set_objective_limit(&*ctx.ms, ul);
                }
            }

            if HighsLpRelaxation::unscaled_dual_feasible(status) {
                if upbranch {
                    st.upbound[ci] = solobj;
                } else {
                    st.downbound[ci] = solobj;
                }
                if solobj > self.get_optimality_limit(ctx) {
                    self.add_bound_exceeding_conflict(ctx);

                    let pruned = solobj > self.get_cutoff_bound(ctx);
                    if pruned {
                        ctx.ms.mipdata_.debug_solution.node_pruned(&self.localdom);
                    }

                    self.localdom.backtrack(ctx.ms);
                    ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);

                    if upbranch {
                        self.branch_downwards(ctx, col, downval, fracval);
                    } else {
                        self.branch_upwards(ctx, col, upval, fracval);
                    }
                    let len = self.nodestack.len();
                    self.nodestack[len - 2].opensubtrees = if pruned { 0 } else { 1 };
                    self.nodestack[len - 2].other_child_lb = solobj;
                    self.nodestack[len - 2].skip_depth_count = 1;
                    self.depthoffset -= 1;

                    return true;
                }
            } else if solobj > self.get_cutoff_bound(ctx) {
                self.add_bound_exceeding_conflict(ctx);
                self.localdom.propagate(ctx.ms);
                let infeas = self.localdom.infeasible();
                if infeas {
                    self.localdom.backtrack(ctx.ms);
                    ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);

                    if upbranch {
                        self.branch_downwards(ctx, col, downval, fracval);
                    } else {
                        self.branch_upwards(ctx, col, upval, fracval);
                    }
                    let len = self.nodestack.len();
                    self.nodestack[len - 2].opensubtrees = 0;
                    self.nodestack[len - 2].skip_depth_count = 1;
                    self.depthoffset -= 1;

                    return true;
                }
            }
        } else if status == HighsLpRelaxation_Status::kInfeasible {
            ctx.ms.mipdata_.debug_solution.node_pruned(&self.localdom);
            self.add_infeasible_conflict(ctx);
            ctx.pseudocost.add_cutoff_observation(col, upbranch);
            self.localdom.backtrack(ctx.ms);
            ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);

            if upbranch {
                self.branch_downwards(ctx, col, downval, fracval);
            } else {
                self.branch_upwards(ctx, col, upval, fracval);
            }
            let len = self.nodestack.len();
            self.nodestack[len - 2].opensubtrees = 0;
            self.nodestack[len - 2].skip_depth_count = 1;
            self.depthoffset -= 1;

            return true;
        } else {
            // in case of an LP error we set the score of this variable to zero to
            // avoid choosing it as branching candidate if possible
            st.downscore[ci] = 0.0;
            st.upscore[ci] = 0.0;
            st.downscorereliable[ci] = true;
            st.upscorereliable[ci] = true;
            self.mark_branching_var_up_reliable_at_node(col);
            self.mark_branching_var_down_reliable_at_node(col);
        }

        self.localdom.backtrack(ctx.ms);
        ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
        false
    }

    /// C++ `while (true)` loop of `selectBranchingCandidate`.
    #[allow(clippy::too_many_arguments)]
    fn sb_loop(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        playground: &mut HighsLpRelaxation_Playground,
        st: &mut SbState,
        fracints: &[(HighsInt, f64)],
        max_sb_iters: i64,
        down_node_lb: &mut f64,
        up_node_lb: &mut f64,
    ) -> HighsInt {
        loop {
            let must_stop = self.get_strong_branching_lp_iterations(&*ctx.ms) >= max_sb_iters
                || self.check_limits(ctx, 0);

            let feastol = self.get_feas_tol(ctx);
            let candidate = sb_select_best_score(
                st,
                fracints,
                &ctx.ms.mipdata_.nodequeue,
                &*ctx.pseudocost,
                feastol,
                must_stop,
            );
            let ci = candidate as usize;

            if (st.upscorereliable[ci] && st.downscorereliable[ci]) || must_stop {
                *down_node_lb = st.downbound[ci];
                *up_node_lb = st.upbound[ci];
                return candidate;
            }

            let upper = self.get_upper_limit(ctx);
            ctx.lp.set_objective_limit(&*ctx.ms, upper);

            let col = fracints[ci].0;
            let down_first = !st.downscorereliable[ci]
                && (st.upscorereliable[ci]
                    || !std_pair_lt(
                        (
                            st.downscore[ci],
                            ctx.pseudocost.get_avg_inferences_down(col),
                        ),
                        (st.upscore[ci], ctx.pseudocost.get_avg_inferences_up(col)),
                    ));

            if down_first {
                // evaluate down branch
                if self.sb_strong_branch(ctx, playground, st, fracints, ci, false) {
                    return -1;
                }
            } else if self.sb_strong_branch(ctx, playground, st, fracints, ci, true) {
                // evaluate up branch
                return -1;
            }
        }
    }
}

impl HighsSearch {
    /// HiGHS: HighsSearch.cpp::HighsSearch::HighsSearch (lines 17-40)
    /// C++: `HighsSearch(HighsMipWorker & mipworker, HighsPseudocost & pseudocost)`
    // PORT NOTE: the C++ stores `mipworker` and `pseudocost` (they are `SearchCtx::worker` /
    // `SearchCtx::pseudocost` at every call) and reads only `mipworker.getGlobalDomain()`,
    // passed as its index `globaldom` (`ms.mipdata_.workers[w].globaldom_`; the worker may
    // not be in `workers` yet when its constructor runs this). `localdom(..)` is
    // `detach!(ms.mipdata_.domains[globaldom], |g| HighsDomain::copy_from(ms, g))`.
    // `HighsRandom random` is default-constructed in C++: `HighsRandom::new(0)`, not
    // `Default` (state 0 without `initialise`).
    pub fn new(ms: &mut HighsMipSolver, globaldom: usize) -> Self {
        let localdom = detach!(ms.mipdata_.domains[globaldom], |g| HighsDomain::copy_from(
            ms, g
        ));
        let mut this = HighsSearch {
            localdom,
            random: HighsRandom::new(0),
            nnodes: 0,
            nleaves: 0,
            lpiterations: 0,
            heurlpiterations: 0,
            sblpiterations: 0,
            eval_origin: crate::work::LP_ORIGIN_CHILD,
            upper_limit: kHighsInf,
            treeweight: HighsCDouble::from(0.0),
            inds: Vec::new(),
            vals: Vec::new(),
            depthoffset: 0,
            inbranching: false,
            inheuristic: false,
            count_tree_weight: true,
            childselrule: if ms.submip {
                HighsSearch_ChildSelectionRule::kHybridInferenceCost
            } else if ms.options_mip_.mip_child_selection >= 0 {
                // rusolve mip_child_selection: one of HiGHS's own rules for the main search
                child_selection_rule(ms.options_mip_.mip_child_selection)
            } else {
                HighsSearch_ChildSelectionRule::kRootSol
            },
            subrootsol: Vec::new(),
            nodestack: Vec::new(),
            stabilizer_orbit_workspace: StabilizerOrbitWorkspace::default(),
            reliableatnode: HighsHashTable::default(),
            s3_stats: crate::mip::symmetry_dynamic::S3Stats::default(),
            plunge_bound_stop: false,
        };
        // the infeasibility flag is overwritten and lost when setDomainChangeStack is
        // called. therefore, assert that localdom is not infeasible here.
        debug_assert!(!this.localdom.infeasible());
        this.localdom.set_domain_change_stack(ms, &[]);
        this
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::checkSol (lines 42-58)
    /// C++: `double checkSol(const std::vector<double> & sol, bool & integerfeasible) const`
    pub fn check_sol(&self, ctx: &SearchCtx<'_>, sol: &[f64], integerfeasible: &mut bool) -> f64 {
        let mipsolver = &*ctx.ms;
        let mut objval = HighsCDouble::from(0.0);
        // HiGHS 1.16 #3149 (h16_cdouble): HighsSearch.cpp::checkSol line 48 (1.16.0)
        let h16_cdouble_obj = mipsolver.options_mip_.h16_cdouble;
        *integerfeasible = true;
        for i in 0..mipsolver.num_col() {
            if h16_cdouble_obj {
                objval += HighsCDouble::from(sol[i as usize]) * mipsolver.col_cost_2(i);
            } else {
                objval += sol[i as usize] * mipsolver.col_cost_2(i);
            }
            debug_assert!(sol[i as usize].is_finite());

            if !*integerfeasible || !mipsolver.is_col_integer(i) {
                continue;
            }

            if fractionality(sol[i as usize], None) > self.get_feas_tol(ctx) {
                *integerfeasible = false;
            }
        }

        f64::from(objval)
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::orbitsValidInChildNode (lines 60-77)
    /// C++: `bool orbitsValidInChildNode(const HighsDomainChange & branchChg) const`
    pub fn orbits_valid_in_child_node(
        &self,
        ms: &HighsMipSolver,
        branch_chg: &HighsDomainChange,
    ) -> bool {
        let branch_col = branch_chg.column;
        // if the variable is integral or we are in an up branch the stabilizer only
        // stays valid if the column has been stabilized
        let curr_node = self.nodestack.last().unwrap();
        match &curr_node.stabilizer_orbits {
            None => return true,
            Some(orbits) => {
                if orbits.orbit_cols.is_empty()
                    || orbits.is_stabilized(&ms.mipdata_.symmetries, branch_col)
                {
                    return true;
                }
            }
        }

        // a down branch stays valid if the variable is binary
        if branch_chg.boundtype == HighsBoundType::kUpper
            && self.localdom.is_global_binary(ms, branch_chg.column)
        {
            return true;
        }

        false
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getCutoffBound (lines 79-81)
    /// C++: `double getCutoffBound() const`
    pub fn get_cutoff_bound(&self, ctx: &SearchCtx<'_>) -> f64 {
        cpp_min(self.get_upper_limit(ctx), self.upper_limit)
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::setRINSNeighbourhood (lines 83-101)
    /// C++: `void setRINSNeighbourhood(const std::vector<double> & basesol, const std::vector<double> & relaxsol)`
    pub fn set_rins_neighbourhood(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        basesol: &[f64],
        relaxsol: &[f64],
    ) {
        let numcol = ctx.ms.num_col();
        let feastol = self.get_feas_tol(ctx);
        for i in 0..numcol {
            let iu = i as usize;
            if !ctx.ms.is_col_integer(i) {
                continue;
            }
            if self.localdom.col_lower_[iu] == self.localdom.col_upper_[iu] {
                continue;
            }

            let intval = (basesol[iu] + 0.5).floor();
            if (relaxsol[iu] - intval).abs() < feastol {
                if self.localdom.col_lower_[iu] < intval {
                    let val = cpp_min(intval, self.localdom.col_upper_[iu]);
                    self.localdom.change_bound_2(
                        ctx.ms,
                        HighsBoundType::kLower,
                        i,
                        val,
                        HighsDomain_Reason::unspecified(),
                    );
                }
                if self.localdom.col_upper_[iu] > intval {
                    let val = cpp_max(intval, self.localdom.col_lower_[iu]);
                    self.localdom.change_bound_2(
                        ctx.ms,
                        HighsBoundType::kUpper,
                        i,
                        val,
                        HighsDomain_Reason::unspecified(),
                    );
                }
            }
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::setRENSNeighbourhood (lines 103-124)
    /// C++: `void setRENSNeighbourhood(const std::vector<double> & lpsol)`
    pub fn set_rens_neighbourhood(&mut self, ctx: &mut SearchCtx<'_>, lpsol: &[f64]) {
        let numcol = ctx.ms.num_col();
        let feastol = self.get_feas_tol(ctx);
        for i in 0..numcol {
            let iu = i as usize;
            if !ctx.ms.is_col_integer(i) {
                continue;
            }
            if self.localdom.col_lower_[iu] == self.localdom.col_upper_[iu] {
                continue;
            }

            let downval = (lpsol[iu] + feastol).floor();
            let upval = (lpsol[iu] - feastol).ceil();

            if self.localdom.col_lower_[iu] < downval {
                let val = cpp_min(downval, self.localdom.col_upper_[iu]);
                self.localdom.change_bound_2(
                    ctx.ms,
                    HighsBoundType::kLower,
                    i,
                    val,
                    HighsDomain_Reason::unspecified(),
                );
                if self.localdom.infeasible() {
                    return;
                }
            }
            if self.localdom.col_upper_[iu] > upval {
                let val = cpp_max(upval, self.localdom.col_lower_[iu]);
                self.localdom.change_bound_2(
                    ctx.ms,
                    HighsBoundType::kUpper,
                    i,
                    val,
                    HighsDomain_Reason::unspecified(),
                );
                if self.localdom.infeasible() {
                    return;
                }
            }
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::createNewNode (lines 126-129)
    /// C++: `void createNewNode()`
    pub fn create_new_node(&mut self) {
        self.nodestack.push(HighsSearch_NodeData::new(
            -kHighsInf, -kHighsInf, None, None,
        ));
        self.nodestack.last_mut().unwrap().domgchg_stack_pos =
            self.localdom.get_domain_change_stack().len() as HighsInt;
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::cutoffNode (lines 131-131)
    /// C++: `void cutoffNode()`
    pub fn cutoff_node(&mut self) {
        self.nodestack.last_mut().unwrap().opensubtrees = 0;
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::setMinReliable (lines 133-135)
    /// C++: `void setMinReliable(HighsInt minreliable)`
    pub fn set_min_reliable(&mut self, ctx: &mut SearchCtx<'_>, minreliable: HighsInt) {
        ctx.pseudocost.set_min_reliable(minreliable);
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::branchDownwards (lines 137-158)
    /// C++: `void branchDownwards(HighsInt col, double newub, double branchpoint)`
    pub fn branch_downwards(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        col: HighsInt,
        newub: f64,
        branchpoint: f64,
    ) {
        let ms = &mut *ctx.ms;
        let flag = ms.options_mip_.mip_symmetry_dynamic;
        let currnode = self.nodestack.last_mut().unwrap();

        debug_assert!(currnode.opensubtrees == 2);
        debug_assert!(ms.is_col_integral(col));

        currnode.opensubtrees = 1;
        currnode.branching_point = branchpoint;
        currnode.branchingdecision.column = col;
        currnode.branchingdecision.boundval = newub;
        currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
        let branchingdecision = currnode.branchingdecision.clone();

        // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at 32dd425
        let s3_child = if flag {
            let parent_path = self.nodestack.last().unwrap().s3_path.clone();
            let depth = self.get_current_depth();
            s3_child_path(&*ms, &self.localdom, parent_path.as_ref(), col, depth)
        } else {
            None
        };

        let domchg_pos = self.localdom.get_domain_change_stack().len() as HighsInt;
        let pass_stabilizer_to_child_node = self.orbits_valid_in_child_node(ms, &branchingdecision);
        self.localdom
            .change_bound(ms, branchingdecision, HighsDomain_Reason::branching());
        let currnode = self.nodestack.last().unwrap();
        let child = HighsSearch_NodeData::new(
            currnode.lower_bound,
            currnode.estimate,
            currnode.node_basis.clone(),
            if pass_stabilizer_to_child_node {
                currnode.stabilizer_orbits.clone()
            } else {
                None
            },
        );
        self.nodestack.push(child);
        self.nodestack.last_mut().unwrap().domgchg_stack_pos = domchg_pos;
        if flag {
            let len = self.nodestack.len();
            self.nodestack[len - 2].s3_child_path = s3_child.clone();
            self.nodestack[len - 1].s3_path = s3_child;
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::branchUpwards (lines 160-181)
    /// C++: `void branchUpwards(HighsInt col, double newlb, double branchpoint)`
    pub fn branch_upwards(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        col: HighsInt,
        newlb: f64,
        branchpoint: f64,
    ) {
        let ms = &mut *ctx.ms;
        let flag = ms.options_mip_.mip_symmetry_dynamic;
        let currnode = self.nodestack.last_mut().unwrap();

        debug_assert!(currnode.opensubtrees == 2);
        debug_assert!(ms.is_col_integral(col));

        currnode.opensubtrees = 1;
        currnode.branching_point = branchpoint;
        currnode.branchingdecision.column = col;
        currnode.branchingdecision.boundval = newlb;
        currnode.branchingdecision.boundtype = HighsBoundType::kLower;
        let branchingdecision = currnode.branchingdecision.clone();

        // mip_symmetry_dynamic: SCIP symmetry_orbitopal.c lines 720-962 at 32dd425
        let s3_child = if flag {
            let parent_path = self.nodestack.last().unwrap().s3_path.clone();
            let depth = self.get_current_depth();
            s3_child_path(&*ms, &self.localdom, parent_path.as_ref(), col, depth)
        } else {
            None
        };

        let domchg_pos = self.localdom.get_domain_change_stack().len() as HighsInt;
        let pass_stabilizer_to_child_node = self.orbits_valid_in_child_node(ms, &branchingdecision);
        self.localdom
            .change_bound(ms, branchingdecision, HighsDomain_Reason::branching());
        let currnode = self.nodestack.last().unwrap();
        let child = HighsSearch_NodeData::new(
            currnode.lower_bound,
            currnode.estimate,
            currnode.node_basis.clone(),
            if pass_stabilizer_to_child_node {
                currnode.stabilizer_orbits.clone()
            } else {
                None
            },
        );
        self.nodestack.push(child);
        self.nodestack.last_mut().unwrap().domgchg_stack_pos = domchg_pos;
        if flag {
            let len = self.nodestack.len();
            self.nodestack[len - 2].s3_child_path = s3_child.clone();
            self.nodestack[len - 1].s3_path = s3_child;
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::addBoundExceedingConflict (lines 183-199)
    /// C++: `void addBoundExceedingConflict()`
    pub fn add_bound_exceeding_conflict(&mut self, ctx: &mut SearchCtx<'_>) {
        let _src = crate::work::conflict_src(crate::work::ConflictSrc::BoundLp);
        let upper = self.get_upper_limit(ctx);
        if upper != kHighsInf {
            let mut rhs: f64 = 0.0;
            let globaldom = self.get_domain(ctx);
            if ctx.lp.compute_dual_proof(
                ctx.ms,
                globaldom,
                upper,
                &mut self.inds,
                &mut self.vals,
                &mut rhs,
                true,
            ) {
                if ctx.ms.mipdata_.domains[globaldom].infeasible() {
                    return;
                }
                let pool = self.get_conflict_pool(ctx);
                let len = self.inds.len() as HighsInt;
                self.localdom.conflict_analysis_2(
                    ctx.ms,
                    &self.inds,
                    &self.vals,
                    len,
                    rhs,
                    pool,
                    globaldom,
                    ctx.pseudocost,
                );

                let cutpool = self.get_cut_pool(ctx);
                let mut cut_gen = HighsCutGeneration::new(&*ctx.ms, &*ctx.lp, cutpool);
                ctx.ms
                    .mipdata_
                    .debug_solution
                    .check_cut(&self.inds, &self.vals, len, rhs);
                cut_gen.generate_conflict(
                    ctx.ms,
                    &*ctx.lp,
                    &self.localdom,
                    globaldom,
                    &mut self.inds,
                    &mut self.vals,
                    &mut rhs,
                );
            }
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::addInfeasibleConflict (lines 201-245)
    /// C++: `void addInfeasibleConflict()`
    pub fn add_infeasible_conflict(&mut self, ctx: &mut SearchCtx<'_>) {
        let _src = crate::work::conflict_src(crate::work::ConflictSrc::InfLp);
        let mut rhs: f64 = 0.0;
        if ctx.lp.get_lp_solver_2().get_model_status() == HighsModelStatus::kObjectiveBound {
            ctx.lp.perform_aging(ctx.ms, false);
        }

        let globaldom = self.get_domain(ctx);
        if ctx.lp.compute_dual_inf_proof(
            &ctx.ms.mipdata_.domains[globaldom],
            &mut self.inds,
            &mut self.vals,
            &mut rhs,
        ) {
            if ctx.ms.mipdata_.domains[globaldom].infeasible() {
                return;
            }
            let pool = self.get_conflict_pool(ctx);
            let len = self.inds.len() as HighsInt;
            self.localdom.conflict_analysis_2(
                ctx.ms,
                &self.inds,
                &self.vals,
                len,
                rhs,
                pool,
                globaldom,
                ctx.pseudocost,
            );

            let cutpool = self.get_cut_pool(ctx);
            let mut cut_gen = HighsCutGeneration::new(&*ctx.ms, &*ctx.lp, cutpool);
            ctx.ms
                .mipdata_
                .debug_solution
                .check_cut(&self.inds, &self.vals, len, rhs);
            cut_gen.generate_conflict(
                ctx.ms,
                &*ctx.lp,
                &self.localdom,
                globaldom,
                &mut self.inds,
                &mut self.vals,
                &mut rhs,
            );
        }
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::selectBranchingCandidate (lines 247-711)
    /// C++: `HighsInt selectBranchingCandidate(int64_t maxSbIters, double & downNodeLb, double & upNodeLb)`
    pub fn select_branching_candidate(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        max_sb_iters: i64,
        down_node_lb: &mut f64,
        up_node_lb: &mut f64,
    ) -> HighsInt {
        debug_assert!(!ctx.lp.get_fractional_integers().is_empty());
        // HiGHS: `fracints` is a `const&` to the LP list. `solveLp`, `run` and the Playground
        // never change that list (only `resolveLp` does), so one copy is exact.
        let fracints: Vec<(HighsInt, f64)> = ctx.lp.get_fractional_integers().to_vec();
        let numfrac = fracints.len();
        let current_lb = self.get_current_lower_bound();
        let feastol = self.get_feas_tol(ctx);

        // rusolve branching-score hook (mip/hooks.rs); one `is_some` test when no hooks
        if let Some(h) = ctx.ms.user_hooks() {
            let cands: Vec<crate::mip::hooks::BranchingCandidate> = fracints
                .iter()
                .map(|&(col, value)| crate::mip::hooks::BranchingCandidate { col, value })
                .collect();
            let info = crate::mip::hooks::BranchingInfo {
                depth: self.get_current_depth(),
                node_lower_bound: current_lb,
                upper_limit: ctx.ms.mipdata_.upper_limit,
                candidates: &cands,
                solution: &ctx.lp.get_solution().col_value,
                col_lower: &self.localdom.col_lower_,
                col_upper: &self.localdom.col_upper_,
            };
            if let Some(k) = h.branching_choice(&info) {
                return k as HighsInt;
            }
        }

        let mut st = SbState {
            upscore: vec![kHighsInf; numfrac],
            downscore: vec![kHighsInf; numfrac],
            upscorereliable: vec![false; numfrac],
            downscorereliable: vec![false; numfrac],
            upbound: vec![current_lb; numfrac],
            downbound: vec![current_lb; numfrac],
            evalqueue: (0..numfrac as HighsInt).collect(),
            min_score: feastol,
        };

        // initialize up and down scores of variables that have a
        // reliable pseudocost so that they do not get evaluated
        #[allow(clippy::needless_range_loop)]
        for k in 0..numfrac {
            let col = fracints[k].0;
            let fracval = fracints[k].1;
            let colu = col as usize;

            let lower_residual = (fracval - self.localdom.col_lower_[colu]) - feastol;
            let lower_ok = lower_residual > 0.0;
            if !lower_ok {
                crate::highs_log_user!(
                    &ctx.ms.options_mip_.log_options,
                    HighsLogType::kError,
                    "HighsSearch::selectBranchingCandidate Error fracval = {} <= {} = {} + {} = localdom.col_lower_[col] + getFeasTol(): Residual {}\n",
                    fracval,
                    self.localdom.col_lower_[colu] + feastol,
                    self.localdom.col_lower_[colu],
                    feastol,
                    lower_residual
                );
            }

            let upper_residual = (self.localdom.col_upper_[colu] - fracval) - feastol;
            let upper_ok = upper_residual > 0.0;
            if !upper_ok {
                crate::highs_log_user!(
                    &ctx.ms.options_mip_.log_options,
                    HighsLogType::kError,
                    "HighsSearch::selectBranchingCandidate Error fracval = {} >= {} = {} - {} = localdom.col_upper_[col] - getFeasTol(): Residual {}\n",
                    fracval,
                    self.localdom.col_upper_[colu] - feastol,
                    self.localdom.col_upper_[colu],
                    feastol,
                    upper_residual
                );
            }

            debug_assert!(lower_residual > -1e-12 && upper_residual > -1e-12);

            if ctx.pseudocost.is_reliable(col) {
                st.upscore[k] = ctx.pseudocost.get_pseudocost_up_2(col, fracval);
                st.downscore[k] = ctx.pseudocost.get_pseudocost_down_2(col, fracval);
                st.upscorereliable[k] = true;
                st.downscorereliable[k] = true;
            } else {
                let flags = self.branching_var_reliable_at_node_flags(col);
                if flags & HighsSearch::kUpReliable != 0 {
                    st.upscore[k] = ctx.pseudocost.get_pseudocost_up_2(col, fracval);
                    st.upscorereliable[k] = true;
                }

                if flags & HighsSearch::kDownReliable != 0 {
                    st.downscore[k] = ctx.pseudocost.get_pseudocost_down_2(col, fracval);
                    st.downscorereliable[k] = true;
                }
            }
        }

        let mut playground = ctx.lp.playground();
        let result = self.sb_loop(
            ctx,
            &mut playground,
            &mut st,
            &fracints,
            max_sb_iters,
            down_node_lb,
            up_node_lb,
        );
        // C++ `~Playground()` runs when `selectBranchingCandidate` returns.
        playground.release(ctx.ms, ctx.lp);
        result
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getParentNodeData (lines 713-717)
    /// C++: `const NodeData * getParentNodeData() const`
    pub fn get_parent_node_data(&self) -> Option<&HighsSearch_NodeData> {
        if self.nodestack.len() <= 1 {
            return None;
        }

        Some(&self.nodestack[self.nodestack.len() - 2])
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::currentNodeToQueue (lines 719-744)
    /// C++: `void currentNodeToQueue(HighsNodeQueue & nodequeue)`
    // PORT NOTE: `nodequeue` is a selector ([`HighsSearch_NodeQueueRef`]); reach the queue
    // with `ctx.node_queue(nodequeue).emplace_node(..)` only at the emplace.
    pub fn current_node_to_queue(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        nodequeue: HighsSearch_NodeQueueRef,
    ) {
        let flag = ctx.ms.options_mip_.mip_symmetry_dynamic;
        let oldchangedcols = self.localdom.get_changed_cols().len();
        let mut prune = self.get_current_lower_bound() > self.get_cutoff_bound(ctx);
        if !prune {
            self.localdom.propagate(ctx.ms);
            self.localdom.clear_changed_cols_2(oldchangedcols);
            prune = self.localdom.infeasible();
            if prune {
                let pool = self.get_conflict_pool(ctx);
                let globaldom = ctx.worker().get_global_domain();
                self.localdom
                    .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
            }
        }
        if !prune {
            let mut branch_positions: Vec<HighsInt> = Vec::new();
            let domchg_stack = self
                .localdom
                .get_reduced_domain_change_stack(&mut branch_positions);
            let lower_bound = cpp_max(
                self.get_current_lower_bound(),
                self.localdom.get_objective_lower_bound(),
            );
            let estimate = self.get_current_estimate();
            let depth = self.get_current_depth();
            // mip_symmetry_dynamic: carry s3_path through the node queue
            // mip_node_install_basis: and the parent basis
            let basis = self.queued_node_basis(ctx);
            let tmp_treeweight = if flag || basis.is_some() {
                let s3_path = if flag {
                    self.nodestack.last().unwrap().s3_path.clone()
                } else {
                    None
                };
                ctx.node_queue(nodequeue).emplace_node_ext(
                    domchg_stack,
                    branch_positions,
                    lower_bound,
                    estimate,
                    depth,
                    s3_path,
                    basis,
                )
            } else {
                ctx.node_queue(nodequeue).emplace_node(
                    domchg_stack,
                    branch_positions,
                    lower_bound,
                    estimate,
                    depth,
                )
            };
            if self.count_tree_weight {
                self.treeweight += tmp_treeweight;
            }
        } else {
            ctx.ms.mipdata_.debug_solution.node_pruned(&self.localdom);
            if self.count_tree_weight {
                self.treeweight += ldexp_f64(1.0, 1 - self.get_current_depth());
            }
        }
        self.nodestack.last_mut().unwrap().opensubtrees = 0;
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::openNodesToQueue (lines 746-794)
    /// C++: `void openNodesToQueue(HighsNodeQueue & nodequeue)`
    // PORT NOTE: `nodequeue` is a selector, as in `current_node_to_queue`.
    pub fn open_nodes_to_queue(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        nodequeue: HighsSearch_NodeQueueRef,
    ) {
        if self.nodestack.is_empty() {
            return;
        }

        // get the basis of the node highest up in the tree
        // (mip_node_install_basis: the node keeps it, so that its queued children carry it)
        let keep_basis = ctx.ms.options_mip_.mip_node_install_basis;
        let mut basis: Option<Arc<HighsBasis>> = None;
        for node_data in self.nodestack.iter_mut() {
            let b = if keep_basis {
                node_data.node_basis.clone()
            } else {
                node_data.node_basis.take()
            };
            if let Some(b) = b {
                basis = Some(b);
                break;
            }
        }

        if self.nodestack.last().unwrap().opensubtrees == 0 {
            self.backtrack(ctx, false);
        }

        // HiGHS: the C++ loop body is `currentNodeToQueue(nodequeue)` followed by `backtrack(false)`.
        while !self.nodestack.is_empty() {
            self.current_node_to_queue(ctx, nodequeue);
            self.backtrack(ctx, false);
        }

        ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
        if let Some(basis) = basis {
            if basis.row_status.len() as HighsInt == ctx.lp.num_rows() {
                ctx.lp.set_stored_basis(Some(basis));
            }
            ctx.lp.recover_basis();
        }
    }

    /// HiGHS 1.16 (#3192): HighsSearch.cpp::stashNodeToProcessed (lines 727-736)
    /// C++: `void stashNodeToProcessed(double lb, double estimate, HighsInt depth)`
    ///
    /// Used only under `h16_mip_workers`. `s3_path` (extra argument) is the branching path
    /// the node carries when `mip_symmetry_dynamic` is on, None
    /// otherwise; the 1.15.1 queue emplaces pass it to `emplace_node_s3` the same way.
    pub(crate) fn stash_node_to_processed(
        &mut self,
        ctx: &mut SearchCtx<'_>,
        lb: f64,
        estimate: f64,
        depth: HighsInt,
        s3_path: Option<Arc<crate::mip::symmetry_dynamic::S3PathNode>>,
    ) {
        let mut branch_positions: Vec<HighsInt> = Vec::new();
        let domchg_stack = self
            .localdom
            .get_reduced_domain_change_stack(&mut branch_positions);
        let mut node =
            HighsNodeQueue_OpenNode::new_2(domchg_stack, branch_positions, lb, estimate, depth);
        node.s3_path = s3_path;
        node.basis = self.queued_node_basis(ctx);
        let count_tree_weight = self.count_tree_weight;
        ctx.ms.mipdata_.workers[ctx.worker]
            .processed_nodes
            .push((node, count_tree_weight));
    }

    /// HiGHS 1.16 (#3192): HighsSearch.cpp::stashCurrentNode (lines 738-758)
    /// C++: `void stashCurrentNode()`
    ///
    /// The 1.16 replacement of `currentNodeToQueue`: the node goes to the worker's
    /// `processedNodes` instead of a node queue, and its tree weight is counted when the
    /// main loop moves it to the global queue. Used only under `h16_mip_workers`.
    pub fn stash_current_node(&mut self, ctx: &mut SearchCtx<'_>) {
        let flag = ctx.ms.options_mip_.mip_symmetry_dynamic;
        let oldchangedcols = self.localdom.get_changed_cols().len();
        let mut prune = self.nodestack.last().unwrap().lower_bound > self.get_cutoff_bound(ctx);
        if !prune {
            self.localdom.propagate(ctx.ms);
            self.localdom.clear_changed_cols_2(oldchangedcols);
            prune = self.localdom.infeasible();
            if prune {
                let pool = self.get_conflict_pool(ctx);
                let globaldom = ctx.worker().get_global_domain();
                self.localdom
                    .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
            }
        }
        if !prune {
            let lower_bound = cpp_max(
                self.nodestack.last().unwrap().lower_bound,
                self.localdom.get_objective_lower_bound(),
            );
            let estimate = self.nodestack.last().unwrap().estimate;
            let depth = self.get_current_depth();
            // mip_symmetry_dynamic: carry s3_path through the node queue
            let s3_path = if flag {
                self.nodestack.last().unwrap().s3_path.clone()
            } else {
                None
            };
            self.stash_node_to_processed(ctx, lower_bound, estimate, depth, s3_path);
        } else {
            ctx.ms.mipdata_.debug_solution.node_pruned(&self.localdom);
            if self.count_tree_weight {
                self.treeweight += ldexp_f64(1.0, 1 - self.get_current_depth());
            }
        }
        self.nodestack.last_mut().unwrap().opensubtrees = 0;
    }

    /// HiGHS 1.16 (#3192): HighsSearch.cpp::stashOpenNodes (lines 760-803)
    /// C++: `void stashOpenNodes()`
    ///
    /// The 1.16 replacement of `openNodesToQueue`; the loop body is the body of
    /// `stashCurrentNode` (1.16 writes it inline). Used only under `h16_mip_workers`.
    pub fn stash_open_nodes(&mut self, ctx: &mut SearchCtx<'_>) {
        if self.nodestack.is_empty() {
            return;
        }

        // get the basis of the node highest up in the tree
        // (mip_node_install_basis: the node keeps it, so that its queued children carry it)
        let keep_basis = ctx.ms.options_mip_.mip_node_install_basis;
        let mut basis: Option<Arc<HighsBasis>> = None;
        for node_data in self.nodestack.iter_mut() {
            let b = if keep_basis {
                node_data.node_basis.clone()
            } else {
                node_data.node_basis.take()
            };
            if let Some(b) = b {
                basis = Some(b);
                break;
            }
        }

        if self.nodestack.last().unwrap().opensubtrees == 0 {
            self.backtrack(ctx, false);
        }

        while !self.nodestack.is_empty() {
            self.stash_current_node(ctx);
            self.backtrack(ctx, false);
        }

        ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
        if let Some(basis) = basis {
            if basis.row_status.len() as HighsInt == ctx.lp.num_rows() {
                ctx.lp.set_stored_basis(Some(basis));
            }
            ctx.lp.recover_basis();
        }
    }

    /// HiGHS 1.16 (#3124): HighsSearch.cpp::resetStatistics (lines 815-822)
    /// C++: `void resetStatistics()`
    ///
    /// Used only under `h16_mip_workers` (the 1.16 `flushStatistics` ends with this call;
    /// the port's `flush_statistics` zeroes the same counters inline).
    pub fn reset_statistics(&mut self) {
        self.nnodes = 0;
        self.nleaves = 0;
        self.treeweight = HighsCDouble::from(0.0);
        self.lpiterations = 0;
        self.heurlpiterations = 0;
        self.sblpiterations = 0;
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::flushStatistics (lines 796-814)
    /// C++: `void flushStatistics(HighsMipSolver & mipsolver)`
    // PORT NOTE: takes only the solver (every C++ call passes the search's own `mipsolver`;
    // callable while the search is outside a `SearchCtx`).
    pub fn flush_statistics(&mut self, mipsolver: &mut HighsMipSolver) {
        mipsolver.mipdata_.num_nodes += self.nnodes;
        self.nnodes = 0;

        mipsolver.mipdata_.num_leaves += self.nleaves;
        self.nleaves = 0;

        mipsolver.mipdata_.pruned_treeweight += self.treeweight;
        self.treeweight = HighsCDouble::from(0.0);

        mipsolver.mipdata_.total_lp_iterations += self.lpiterations;
        self.lpiterations = 0;

        mipsolver.mipdata_.heuristic_lp_iterations += self.heurlpiterations;
        self.heurlpiterations = 0;

        mipsolver.mipdata_.sb_lp_iterations += self.sblpiterations;
        self.sblpiterations = 0;
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getHeuristicLpIterations (lines 816-818)
    /// C++: `int64_t getHeuristicLpIterations() const`
    pub fn get_heuristic_lp_iterations(&self, ms: &HighsMipSolver) -> i64 {
        self.heurlpiterations + ms.mipdata_.heuristic_lp_iterations
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getTotalLpIterations (lines 820-822)
    /// C++: `int64_t getTotalLpIterations() const`
    pub fn get_total_lp_iterations(&self, ms: &HighsMipSolver) -> i64 {
        self.lpiterations + ms.mipdata_.total_lp_iterations
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getLocalLpIterations (lines 824-824)
    /// C++: `int64_t getLocalLpIterations() const`
    pub fn get_local_lp_iterations(&self) -> i64 {
        self.lpiterations
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getLocalNodes (lines 826-826)
    /// C++: `int64_t & getLocalNodes()`
    pub fn get_local_nodes(&mut self) -> &mut i64 {
        &mut self.nnodes
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getLocalLeaves (lines 828-828)
    /// C++: `int64_t & getLocalLeaves()`
    pub fn get_local_leaves(&mut self) -> &mut i64 {
        &mut self.nleaves
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::getStrongBranchingLpIterations (lines 830-832)
    /// C++: `int64_t getStrongBranchingLpIterations() const`
    pub fn get_strong_branching_lp_iterations(&self, ms: &HighsMipSolver) -> i64 {
        self.sblpiterations + ms.mipdata_.sb_lp_iterations
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::resetLocalDomain (lines 834-846)
    /// C++: `void resetLocalDomain()`
    pub fn reset_local_domain(&mut self, ctx: &mut SearchCtx<'_>) {
        let globaldom = self.get_domain(ctx);
        ctx.lp
            .reset_to_global_domain(&*ctx.ms, &ctx.ms.mipdata_.domains[globaldom]);
        detach!(ctx.ms.mipdata_.domains[globaldom], |g| self
            .localdom
            .op_assign(ctx.ms, g));

        for i in 0..ctx.ms.num_col() {
            let iu = i as usize;
            debug_assert!(
                ctx.lp.get_lp_solver_2().get_lp().col_lower_[iu] == self.localdom.col_lower_[iu]
                    || ctx.ms.is_col_continuous(i)
            );
            debug_assert!(
                ctx.lp.get_lp_solver_2().get_lp().col_upper_[iu] == self.localdom.col_upper_[iu]
                    || ctx.ms.is_col_continuous(i)
            );
        }
    }

    /// rusolve mip_node_install_basis: the basis a node sent to the queue now carries, the one
    /// stored with the top of the node stack (the node's own basis once its LP is solved, its
    /// parent's before). None when the flag is off, or when the open nodes would hold more
    /// than `NODE_BASIS_QUEUE_BUDGET` basis statuses if each had its own basis (siblings share
    /// one, so the real use is lower).
    pub(crate) fn queued_node_basis(&self, ctx: &SearchCtx<'_>) -> Option<Arc<HighsBasis>> {
        const NODE_BASIS_QUEUE_BUDGET: usize = 1 << 27;
        if !ctx.ms.options_mip_.mip_node_install_basis {
            return None;
        }
        let basis = self.nodestack.last()?.node_basis.as_ref()?;
        let size = basis.col_status.len() + basis.row_status.len();
        let queued = ctx.ms.mipdata_.nodequeue.num_nodes() as usize + 1;
        if queued.saturating_mul(size) > NODE_BASIS_QUEUE_BUDGET {
            return None;
        }
        Some(basis.clone())
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::installNode (lines 848-872)
    /// C++: `void installNode(HighsNodeQueue::OpenNode && node)`
    pub fn install_node(&mut self, ctx: &mut SearchCtx<'_>, node: HighsNodeQueue_OpenNode) {
        let flag = ctx.ms.options_mip_.mip_symmetry_dynamic;
        self.localdom
            .set_domain_change_stack_2(ctx.ms, &node.domchgstack, &node.branchings);
        let mut global_symmetries_valid = true;
        if ctx.ms.mipdata_.global_orbits.is_some() {
            // if global orbits have been computed we check whether they are still valid
            // in this node
            let globaldom = self.get_domain(ctx);
            for &i in self.localdom.get_branching_positions() {
                let iu = i as usize;
                let col = self.localdom.get_domain_change_stack()[iu].column;
                if self.get_symmetries(ctx).column_position[col as usize] == -1 {
                    continue;
                }

                let chg = &self.localdom.get_domain_change_stack()[iu];
                if !ctx.ms.mipdata_.domains[globaldom].is_binary(&*ctx.ms, col)
                    || (chg.boundtype == HighsBoundType::kLower && chg.boundval == 1.0)
                {
                    global_symmetries_valid = false;
                    break;
                }
            }
        }
        let orbits = if global_symmetries_valid {
            ctx.ms.mipdata_.global_orbits.clone()
        } else {
            None
        };
        self.nodestack.push(HighsSearch_NodeData::new(
            node.lower_bound,
            node.estimate,
            None,
            orbits,
        ));
        // mip_symmetry_dynamic: carry s3_path from the queued node
        if flag {
            self.nodestack.last_mut().unwrap().s3_path = node.s3_path_handle();
        }
        // rusolve mip_node_install_basis: start the node's LP from its parent's basis
        // (lp_relaxation/node_basis.rs); HiGHS keeps the basis the LP holds.
        if let Some(basis) = node.basis.as_deref() {
            ctx.lp.restore_node_basis(basis);
        }
        self.subrootsol.clear();
        self.depthoffset = node.depth - 1;
        self.eval_origin = crate::work::LP_ORIGIN_JUMP;
    }

    /// HiGHS: HighsSearch.cpp::HighsSearch::evaluateNode (lines 874-1097)
    /// C++: `NodeResult evaluateNode()`
    pub fn evaluate_node(&mut self, ctx: &mut SearchCtx<'_>) -> HighsSearch_NodeResult {
        debug_assert!(!self.nodestack.is_empty());
        let flag = ctx.ms.options_mip_.mip_symmetry_dynamic;
        let n = self.nodestack.len() - 1;
        // C++ `getParentNodeData()`: nullptr when only the current node is on the stack.
        let parent: Option<usize> = n.checked_sub(1);
        // C++ reads `parent->lp_objective`, `parent->branching_point` and the parent's branching
        // decision at each cutoff site. None of them is changed by this function, so the
        // value the sites share is read once here.
        let parent_cutoff: Option<(HighsInt, bool)> = match parent {
            Some(p) => {
                let pn = &self.nodestack[p];
                if pn.lp_objective != -kHighsInf
                    && pn.branching_point != pn.branchingdecision.boundval
                {
                    Some((
                        pn.branchingdecision.column,
                        pn.branchingdecision.boundtype == HighsBoundType::kLower,
                    ))
                } else {
                    None
                }
            }
            None => None,
        };

        if !self.inheuristic && self.nodestack[n].lower_bound > self.get_optimality_limit(ctx) {
            return HighsSearch_NodeResult::kSubOptimal;
        }

        self.localdom.propagate(ctx.ms);

        if !self.inheuristic && !self.localdom.infeasible() {
            if flag {
                // mip_symmetry_dynamic: SCIP sym_lexorbred.c lines 1091-1148,
                // sym_rowcol.c lines 2530-2570, symmetry_orbitopal.c lines 2139-2206 at 32dd425
                // SCIPgetDepth of the focus node: the HiGHS root is depth 1, SCIP's is 0 (D7:
                // root branching = 1), so SCIP's "do not run in the root" check sees depth - 1.
                let path = self.nodestack[n].s3_path.clone();
                let depth = self.get_current_depth() - 1;
                let infeasible = s3_propagate(
                    ctx.ms,
                    &mut self.localdom,
                    path.as_ref(),
                    depth,
                    &mut self.s3_stats,
                );
                // SCIP_CUTOFF: an infeasibility found without a bound change (D2, peeks) is
                // recorded on the domain, as HighsOrbitopeMatrix::orbital_fixing does.
                if infeasible && !self.localdom.infeasible() {
                    self.localdom
                        .mark_infeasible(HighsDomain_Reason::unspecified());
                }
            } else {
                let parent_allows_orbits = match parent {
                    None => true,
                    Some(p) => match &self.nodestack[p].stabilizer_orbits {
                        None => true,
                        Some(o) => !o.orbit_cols.is_empty(),
                    },
                };
                if ctx.ms.mipdata_.symmetries.num_perms > 0
                    && self.nodestack[n].stabilizer_orbits.is_none()
                    && parent_allows_orbits
                {
                    let orbits = ctx.ms.with_symmetries_ref(|ms, sym| {
                        sym.compute_stabilizer_orbits(
                            ms,
                            &self.localdom,
                            &mut self.stabilizer_orbit_workspace,
                        )
                    });
                    self.nodestack[n].stabilizer_orbits = orbits;
                }

                match self.nodestack[n].stabilizer_orbits.clone() {
                    Some(orbits) => {
                        self.orbital_fixing(ctx.ms, &orbits);
                    }
                    None => {
                        self.propagate_orbitopes(ctx.ms);
                    }
                }
            }
        }

        if let Some(p) = parent {
            // C++ `size_t` arithmetic: `domchgstack.size() - (domgchgStackPos + 1)`.
            let inferences: HighsInt = self.localdom.get_domain_change_stack().len() as HighsInt
                - (self.nodestack[n].domgchg_stack_pos as HighsInt + 1);
            let col = self.nodestack[p].branchingdecision.column;
            let upbranch = self.nodestack[p].branchingdecision.boundtype == HighsBoundType::kLower;
            ctx.pseudocost
                .add_inference_observation(col, inferences, upbranch);
        }

        let mut result = HighsSearch_NodeResult::kOpen;

        if self.localdom.infeasible() {
            result = HighsSearch_NodeResult::kDomainInfeasible;
            self.localdom.clear_changed_cols();
            if let Some((col, upbranch)) = parent_cutoff {
                ctx.pseudocost.add_cutoff_observation(col, upbranch);
            }
            let pool = self.get_conflict_pool(ctx);
            let globaldom = ctx.worker().get_global_domain();
            self.localdom
                .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
        } else {
            ctx.lp.flush_domain(&*ctx.ms, &mut self.localdom, false);
            if ctx.ms.options_mip_.mip_node_bigm {
                // rusolve mip_node_bigm: local big-M tightening
                if self.inheuristic {
                    ctx.lp.restore_local_mods();
                } else {
                    let _ = ctx.lp.apply_local_bigm(&*ctx.ms, &self.localdom);
                }
            }
            let upper = self.get_upper_limit(ctx);
            ctx.lp.set_objective_limit(&*ctx.ms, upper);

            for i in 0..ctx.ms.num_col() {
                let iu = i as usize;
                debug_assert!(
                    ctx.lp.get_lp_solver_2().get_lp().col_lower_[iu]
                        == self.localdom.col_lower_[iu]
                        || ctx.ms.is_col_continuous(i)
                );
                debug_assert!(
                    ctx.lp.get_lp_solver_2().get_lp().col_upper_[iu]
                        == self.localdom.col_upper_[iu]
                        || ctx.ms.is_col_continuous(i)
                );
            }
            let oldnumiters = ctx.lp.get_num_lp_iterations();
            let lp_src_guard = crate::work::lp_kind(crate::work::lp_node_src(
                std::mem::replace(&mut self.eval_origin, crate::work::LP_ORIGIN_CHILD),
                (self.depthoffset.max(0) as usize) + self.nodestack.len(),
            ));
            let status = ctx.lp.resolve_lp(ctx.ms, Some(&mut self.localdom));
            drop(lp_src_guard);
            self.lpiterations += ctx.lp.get_num_lp_iterations() - oldnumiters;

            self.nodestack[n].lower_bound = cpp_max(
                self.localdom.get_objective_lower_bound(),
                self.nodestack[n].lower_bound,
            );

            if self.localdom.infeasible() {
                result = HighsSearch_NodeResult::kDomainInfeasible;
                self.localdom.clear_changed_cols();
                if let Some((col, upbranch)) = parent_cutoff {
                    ctx.pseudocost.add_cutoff_observation(col, upbranch);
                }
                let pool = self.get_conflict_pool(ctx);
                let globaldom = ctx.worker().get_global_domain();
                self.localdom
                    .conflict_analysis(ctx.ms, pool, globaldom, ctx.pseudocost);
            } else if HighsLpRelaxation::scaled_optimal(status) {
                ctx.lp.store_basis();
                ctx.lp.perform_aging(ctx.ms, false);

                self.nodestack[n].node_basis = ctx.lp.get_stored_basis();
                let estimate = ctx.lp.compute_best_estimate(&*ctx.ms, &*ctx.pseudocost);
                self.nodestack[n].estimate = estimate;
                self.nodestack[n].lp_objective = ctx.lp.get_objective();

                if let Some(p) = parent {
                    if self.nodestack[p].lp_objective != -kHighsInf
                        && self.nodestack[p].branching_point
                            != self.nodestack[p].branchingdecision.boundval
                    {
                        let delta = self.nodestack[p].branchingdecision.boundval
                            - self.nodestack[p].branching_point;
                        let objdelta = cpp_max(
                            0.0,
                            self.nodestack[n].lp_objective - self.nodestack[p].lp_objective,
                        );
                        let col = self.nodestack[p].branchingdecision.column;
                        ctx.pseudocost.add_observation(col, delta, objdelta);
                    }
                }

                if HighsLpRelaxation::unscaled_primal_feasible(status)
                    && ctx.lp.get_fractional_integers().is_empty()
                {
                    let cutoffbnd = self.get_cutoff_bound(ctx);
                    // HiGHS: clone, since `addIncumbent` lends the LP back to `ms`.
                    let sol: Vec<f64> = ctx.lp.get_solution().col_value.clone();
                    let objective = ctx.lp.get_objective();
                    let source = if self.inheuristic {
                        kSolutionSourceHeuristic
                    } else {
                        kSolutionSourceEvaluateNode
                    };
                    self.add_incumbent(ctx, &sol, objective, source, true);
                    if self.get_upper_limit(ctx) < cutoffbnd {
                        let ul = self.get_upper_limit(ctx);
                        ctx.lp.set_objective_limit(&*ctx.ms, ul);
                    }

                    if HighsLpRelaxation::unscaled_dual_feasible(status) {
                        self.add_bound_exceeding_conflict(ctx);
                        result = HighsSearch_NodeResult::kBoundExceeding;
                    }
                }

                if result == HighsSearch_NodeResult::kOpen {
                    if HighsLpRelaxation::unscaled_dual_feasible(status) {
                        self.nodestack[n].lower_bound = cpp_max(
                            self.nodestack[n].lp_objective,
                            self.nodestack[n].lower_bound,
                        );

                        if self.nodestack[n].lower_bound > self.get_cutoff_bound(ctx) {
                            result = HighsSearch_NodeResult::kBoundExceeding;
                            self.add_bound_exceeding_conflict(ctx);
                        } else if self.get_upper_limit(ctx) != kHighsInf {
                            if !self.inheuristic {
                                let gap = self.get_upper_limit(ctx) - ctx.lp.get_objective();
                                let threshold = gap
                                    + cpp_max(
                                        10.0 * self.get_feas_tol(ctx),
                                        self.get_epsilon(ctx) * gap,
                                    );
                                let globaldom = self.get_domain(ctx);
                                let conflictpool = self.get_conflict_pool(ctx);
                                ctx.with_worker_pseudocost(|ms, lp, pc| {
                                    lp.compute_basic_degenerate_duals(
                                        ms,
                                        threshold,
                                        &mut self.localdom,
                                        globaldom,
                                        conflictpool,
                                        pc,
                                        true,
                                    );
                                });
                            }
                            let upper = self.get_upper_limit(ctx);
                            let globaldom = self.get_domain(ctx);
                            let conflictpool = self.get_conflict_pool(ctx);
                            ctx.with_worker_pseudocost(|ms, lp, pc| {
                                HighsRedcostFixing::propagate_red_cost(
                                    ms,
                                    &mut self.localdom,
                                    globaldom,
                                    lp,
                                    conflictpool,
                                    pc,
                                    upper,
                                );
                            });
                            self.localdom.propagate(ctx.ms);
                            if self.localdom.infeasible() {
                                result = HighsSearch_NodeResult::kDomainInfeasible;
                                self.localdom.clear_changed_cols();
                                if let Some((col, upbranch)) = parent_cutoff {
                                    ctx.pseudocost.add_cutoff_observation(col, upbranch);
                                }
                                let pool = self.get_conflict_pool(ctx);
                                let globaldom = ctx.worker().get_global_domain();
                                self.localdom.conflict_analysis(
                                    ctx.ms,
                                    pool,
                                    globaldom,
                                    ctx.pseudocost,
                                );
                            } else if !self.localdom.get_changed_cols().is_empty() {
                                return self.evaluate_node(ctx);
                            }
                        } else if !self.inheuristic {
                            let globaldom = self.get_domain(ctx);
                            let conflictpool = self.get_conflict_pool(ctx);
                            ctx.with_worker_pseudocost(|ms, lp, pc| {
                                lp.compute_basic_degenerate_duals(
                                    ms,
                                    kHighsInf,
                                    &mut self.localdom,
                                    globaldom,
                                    conflictpool,
                                    pc,
                                    true,
                                );
                            });
                            self.localdom.propagate(ctx.ms);
                            if self.localdom.infeasible() {
                                result = HighsSearch_NodeResult::kDomainInfeasible;
                                self.localdom.clear_changed_cols();
                                if let Some((col, upbranch)) = parent_cutoff {
                                    ctx.pseudocost.add_cutoff_observation(col, upbranch);
                                }
                                let pool = self.get_conflict_pool(ctx);
                                let globaldom = ctx.worker().get_global_domain();
                                self.localdom.conflict_analysis(
                                    ctx.ms,
                                    pool,
                                    globaldom,
                                    ctx.pseudocost,
                                );
                            } else if !self.localdom.get_changed_cols().is_empty() {
                                return self.evaluate_node(ctx);
                            }
                        }
                    } else if ctx.lp.get_objective() > self.get_cutoff_bound(ctx) {
                        // the LP is not solved to dual feasibility due to scaling/numerics
                        // therefore we compute a conflict constraint as if the LP was bound
                        // exceeding and propagate the local domain again. The lp relaxation
                        // class will take care to consider the dual multipliers with an
                        // increased zero tolerance due to the dual infeasibility when
                        // computing the proof conBoundExceedingstraint.
                        self.add_bound_exceeding_conflict(ctx);
                        self.localdom.propagate(ctx.ms);
                        if self.localdom.infeasible() {
                            result = HighsSearch_NodeResult::kBoundExceeding;
                        }
                    }
                }
            } else if status == HighsLpRelaxation_Status::kInfeasible {
                if ctx.lp.get_lp_solver_2().get_model_status() == HighsModelStatus::kObjectiveBound
                {
                    result = HighsSearch_NodeResult::kBoundExceeding;
                } else {
                    result = HighsSearch_NodeResult::kLpInfeasible;
                }
                self.add_infeasible_conflict(ctx);
                if let Some((col, upbranch)) = parent_cutoff {
                    ctx.pseudocost.add_cutoff_observation(col, upbranch);
                }
            }
        }

        if result != HighsSearch_NodeResult::kOpen {
            ctx.ms.mipdata_.debug_solution.node_pruned(&self.localdom);
            self.treeweight += ldexp_f64(1.0, 1 - self.get_current_depth());
            self.nodestack[n].opensubtrees = 0;
        } else if !self.inheuristic
            && self.nodestack[n].lower_bound > self.get_optimality_limit(ctx)
        {
            result = HighsSearch_NodeResult::kSubOptimal;
            self.add_bound_exceeding_conflict(ctx);
        }

        result
    }
}

#[cfg(test)]
mod s3_gate1_tests {
    //! Tests: the flag gates in part1.rs run the knapsack solve with the option off
    //! and on. The helpers are copies of the search/part2.rs test helpers (private there).

    use super::*;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    const N_ROWS: usize = 3;

    /// Strongly correlated multi-row binary knapsack, as in search/part2.rs tests.
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
        // The HiGHS child rule: the pinned node counts below are the HiGHS path.
        options.mip_child_selection = -1;
        options
    }

    /// knapsack_lp(5) with every item repeated `copies` times: the copies of an item are
    /// interchangeable, so symmetry detection finds permutations.
    fn symmetric_knapsack_lp(copies: usize) -> HighsLp {
        let base = knapsack_lp(5);
        let n = 5 * copies;
        let mut lp = base.clone();
        lp.num_col_ = n as HighsInt;
        lp.col_cost_ = (0..n).map(|j| base.col_cost_[j / copies]).collect();
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_upper_ = base.row_upper_.iter().map(|u| u * copies as f64).collect();
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.start_ = (0..=n).map(|j| (N_ROWS * j) as HighsInt).collect();
        lp.a_matrix_.index_ = (0..n).flat_map(|_| [0, 1, 2]).collect();
        lp.a_matrix_.value_ = (0..n)
            .flat_map(|j| {
                let k = N_ROWS * (j / copies);
                [
                    base.a_matrix_.value_[k],
                    base.a_matrix_.value_[k + 1],
                    base.a_matrix_.value_[k + 2],
                ]
            })
            .collect();
        lp
    }

    /// Solves `lp` with the given flag value. Returns status, objective, node count and
    /// whether the data is absent.
    fn solve(lp: HighsLp, on: bool) -> (HighsModelStatus, f64, i64, bool) {
        std::thread::spawn(move || {
            initialize_scheduler(1);
            let mut options = base_options();
            options.mip_symmetry_dynamic = on;
            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let out = (
                ms.modelstatus_,
                ms.solution_objective_,
                ms.node_count_,
                ms.mipdata_.symmetry_dynamic.is_none(),
            );
            drop(ms);
            HighsTaskExecutor::shutdown(true);
            out
        })
        .join()
        .unwrap()
    }

    #[test]
    fn s3_gate1_flag_off_knapsack() {
        let (status, obj, nodes, s3_absent) = solve(knapsack_lp(20), false);
        assert_eq!(status, HighsModelStatus::kOptimal);
        assert!((obj - brute_force_optimum(&knapsack_lp(20))).abs() < 1e-6);
        // Node count and objective of s3-int (before this job) on the same run: the off
        // path is unchanged (rule F).
        assert_eq!(nodes, 1016, "node count {nodes}");
        assert_eq!(obj, -691.0);
        // Off: no data is built (rule F).
        assert!(s3_absent);
    }

    #[test]
    fn s3_gate1_flag_on_knapsack() {
        let (status, obj, _nodes, _) = solve(knapsack_lp(20), true);
        assert_eq!(status, HighsModelStatus::kOptimal);
        assert!((obj - brute_force_optimum(&knapsack_lp(20))).abs() < 1e-6);
    }

    #[test]
    fn s3_gate1_flag_on_symmetric_knapsack() {
        let lp = symmetric_knapsack_lp(3);
        let opt = brute_force_optimum(&lp);
        let (status_off, obj_off, nodes_off, _) = solve(lp.clone(), false);
        let (status_on, obj_on, nodes_on, s3_absent_on) = solve(lp, true);
        eprintln!(
            "s3_gate1 symmetric: off nodes {nodes_off}, on nodes {nodes_on}, s3 absent {s3_absent_on}"
        );
        assert_eq!(status_off, HighsModelStatus::kOptimal);
        assert_eq!(status_on, HighsModelStatus::kOptimal);
        assert!((obj_off - opt).abs() < 1e-6);
        assert!((obj_on - opt).abs() < 1e-6);
    }
}
