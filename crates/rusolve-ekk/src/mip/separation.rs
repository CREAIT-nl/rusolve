//! HiGHS source: mip/HighsSeparation.h, mip/HighsSeparation.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsSeparation.h::HighsSeparation (lines 23-41)
///
/// Removed C++ members:
/// * C++: `HighsMipWorker & mipworker_`: removed; the
///   methods take [`SepaCtx`], whose `worker` the caller detached from `ms.mipdata_.workers`.
/// * C++ `HighsLpRelaxation* lp`: kept as the index [`HighsSeparation::lp`] into
///   `ms.mipdata_.lps` (every C++ `setLpRelaxation` passes a worker's `lp_` or `&getLp()`,
///   both elements of `mipdata_->lps`). The caller detaches `ms.mipdata_.lps[lp]` into
///   `SepaCtx::lp` for `separate` / `separationRound`.
///
/// Ownership in `separationRound`:
/// * `propdomain` is detached by the caller (the node domain, or `domains[0]` at the root).
///   The C++ `&propdomain == &mipdata.getDomain()` / `== &mipworker_.getGlobalDomain()` tests
///   are `ms.mipdata_.domains[i].holds(propdomain)`.
/// * Around `lp->resolveLp(&propdomain)`, `addRootRedcost` / `propagateRootRedcost` (they read
///   `getLp()` and `getDomain()`) and the clique/implied-bound separators, `ctx.lp`,
///   `ctx.worker` and `propdomain` are lent back with
///   `lend_in!(ms.mipdata_.lps, ctx.lp, ..)` / `lend_in!(ms.mipdata_.workers, ctx.worker, ..)`
///   / `lend_in!(ms.mipdata_.domains, propdomain, ..)` whenever the callee reaches them
///   through `ms`. `resolve_lp` itself takes the LP as `self`, so for that call only the
///   worker and the domain are lent.
/// * `mipworker_.getCutPool()` / `getConflictPool()` / `getPseudocost()` /
///   `getGlobalDomain()` are the worker's indices `cutpool_` / `conflictpool_` /
///   `pseudocost_` / `globaldom_`.
#[derive(Clone, Default, Debug)]
pub struct HighsSeparation {
    /// C++: `HighsInt implBoundClock`
    pub impl_bound_clock: HighsInt,
    /// C++: `HighsInt cliqueClock`
    pub clique_clock: HighsInt,
    /// C++ `std::vector<std::unique_ptr<HighsSeparator>> separators` (virtual dispatch -> enum).
    pub separators: Vec<HighsSeparatorBox>,
    /// C++: `HighsCutSet cutset`
    pub cutset: HighsCutSet,
    /// C++: `HighsLpRelaxation* lp` (index into `ms.mipdata_.lps`; see the struct doc).
    pub lp: Option<usize>,
    /// SCIP cutsel_hybrid.c randnumgen (cutseldata :62-74); None unless mip_cut_selection_hybrid
    pub(crate) s4_randnumgen: Option<crate::mip::cut_selection::S4RandNumGen>,
    /// rusolve mip_node_sepa_pool_only: set by the tree-node gate (`HighsMipSolver::
    /// run_separate_and_store_basis`) for one `separate` call. The round then only propagates
    /// and adds violated cut-pool cuts: implied bounds, cliques, degenerate-dual conflicts, the
    /// separators and the user hook are skipped. Always false otherwise.
    pub(crate) pool_only: bool,
}

/// rusolve: a process-wide increasing number for separator-built solutions; only the relative
/// order of the numbers taken by one worker's separators matters.
pub fn next_pending_seq() -> u64 {
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
}

impl HighsSeparation {
    /// rusolve: solutions the TI / network-design separators built during the last round. They
    /// live in the separator (owned by one worker), never in shared solver state.
    ///
    /// Returned in the order the separators produced them (the single `ti_pending_solutions`
    /// queue this replaces interleaved the two kinds by production order): each entry carries
    /// a sequence number from `next_pending_seq`, which only orders entries within one worker.
    pub fn take_pending_solutions(&mut self) -> Vec<Vec<f64>> {
        let mut out = Vec::new();
        for s in &mut self.separators {
            match s {
                HighsSeparatorBox::TiPrecedence(s) => out.append(&mut s.pending),
                HighsSeparatorBox::Netdesign(s) => out.append(&mut s.pending),
                _ => {}
            }
        }
        out.sort_by_key(|e| e.0);
        out.into_iter().map(|e| e.1).collect()
    }

    /// HiGHS: HighsSeparation.cpp::HighsSeparation::HighsSeparation (lines 25-41)
    /// C++: `HighsSeparation(HighsMipWorker & mipworker)`
    // PORT NOTE: the worker is not stored (it is `SepaCtx::worker` at every call); the body
    // reads only `mipworker.getMipSolver()`, which is `ms`.
    pub fn new(ms: &HighsMipSolver) -> Self {
        // HiGHS #2990 (mip_machine_sched_sepa): HighsSeparation.cpp lines 41-42 at dc83470a9
        let mut separators = vec![
            HighsSeparatorBox::Tableau(HighsTableauSeparator::new(ms)),
            HighsSeparatorBox::Path(HighsPathSeparator::new(ms)),
            HighsSeparatorBox::Modk(HighsModkSeparator::new(ms)),
        ];
        if ms.options_mip_.mip_machine_sched_sepa {
            separators.push(HighsSeparatorBox::MachineSched(
                crate::mip::machine_sched_separator::HighsMachineSchedSeparator::new(ms),
            ));
        }
        if ms.options_mip_.mip_knapsack_cover_sepa {
            separators.push(HighsSeparatorBox::KnapsackCover(
                crate::mip::knapsack_cover_separator::HighsKnapsackCoverSeparator::new(ms),
            ));
        }
        if ms.options_mip_.mip_dff_separator {
            separators.push(HighsSeparatorBox::Dff(
                crate::mip::dff_separator::HighsDffSeparator::new(ms),
            ));
        }
        if ms.options_mip_.mip_ti_precedence_sepa || ms.options_mip_.mip_ti_list_sched {
            separators.push(HighsSeparatorBox::TiPrecedence(
                crate::mip::ti_precedence_separator::HighsTiPrecedenceSeparator::new(ms),
            ));
        }
        if ms.options_mip_.mip_netdesign_sepa
            || ms.options_mip_.mip_netdesign_heur
            || ms.options_mip_.mip_netdesign_cutset
        {
            separators.push(HighsSeparatorBox::Netdesign(
                crate::mip::netdesign_separator::HighsNetdesignSeparator::new(ms),
            ));
        }
        // mip_cut_selection_hybrid: SCIP cutsel_hybrid.c lines 306-317 at 32dd425
        let s4_randnumgen = if ms.options_mip_.mip_cut_selection_hybrid {
            Some(crate::mip::cut_selection::S4RandNumGen::new(
                crate::mip::cut_selection::s4_initialize_random_seed(ms.options_mip_.random_seed),
            ))
        } else {
            None
        };
        HighsSeparation {
            impl_bound_clock: 990,
            clique_clock: 991,
            separators,
            cutset: HighsCutSet::default(),
            lp: None,
            s4_randnumgen,
            pool_only: false,
        }
    }

    /// HiGHS: HighsSeparation.cpp::HighsSeparation::separationRound (lines 43-179)
    /// C++: `HighsInt separationRound(HighsDomain & propdomain, HighsLpRelaxation::Status & status)`
    pub fn separation_round(
        &mut self,
        ctx: &mut SepaCtx<'_>,
        propdomain: &mut HighsDomain,
        status: &mut HighsLpRelaxation_Status,
    ) -> HighsInt {
        // C++ `const HighsSolution& sol = lp->getLpSolver().getSolution()` is a reference into the
        // LP; it is re-read at each use (`ctx.lp.get_solution()`) because the LP changes in between.
        let globaldom = ctx.worker.globaldom_;
        let cutpool = ctx.worker.cutpool_;
        ctx.ms.stats_.sepa_rounds.add(1); // reporting-only

        // rusolve mip_node_sepa_pool_only: propagation and the cut pool only (see `pool_only`)
        if self.pool_only {
            return self.pool_only_round(ctx, propdomain, status);
        }

        if !ctx.ms.parallel_lock_active() {
            ctx.ms.profiling_.start(self.impl_bound_clock, false);
        }
        let feastol = ctx.ms.mipdata_.feastol;
        let thread_safe = ctx.ms.parallel_lock_active();
        // C++ `mipdata.implications.separateImpliedBounds(*lp, sol.col_value, getCutPool(), feastol,
        // getGlobalDomain(), parallelLockActive())`. The global domain may be `propdomain` itself
        // (root), so it is lent back for the call.
        lend_in!(
            ctx.ms.mipdata_.domains,
            propdomain,
            detach!(ctx.ms.mipdata_.implications, |imp| detach!(
                ctx.ms.mipdata_.cutpools[cutpool],
                |cp| imp.separate_implied_bounds(
                    ctx.ms,
                    &*ctx.lp,
                    &ctx.lp.get_solution().col_value,
                    cp,
                    feastol,
                    globaldom,
                    thread_safe
                )
            ))
        );
        if !ctx.ms.parallel_lock_active() {
            ctx.ms.profiling_.stop(self.impl_bound_clock);
        }

        let mut ncuts: HighsInt = 0;
        let mut numboundchgs = propagate_and_resolve(ctx, propdomain, status);
        if numboundchgs == -1 {
            return 0;
        } else {
            ncuts += numboundchgs;
        }

        if !ctx.ms.parallel_lock_active() {
            ctx.ms.profiling_.start(self.clique_clock, false);
        }
        {
            // serial path: the clique table's own randgen / counter (`None`, HighsSeparation.cpp:110-116);
            // parallel path: the worker's.
            let (randgen, local_nq): (Option<&mut HighsRandom>, Option<&mut i64>) =
                if ctx.ms.parallel_lock_active() {
                    (
                        Some(&mut ctx.worker.randgen),
                        Some(&mut ctx.worker.sepa_stats.num_neighbourhood_queries),
                    )
                } else {
                    (None, None)
                };
            let sol = &ctx.lp.get_solution().col_value;
            // `separateCliques` reads `mipdata_->getDomain()` (the clique subsumption), which is
            // `propdomain` at the root: lend it back for the call.
            lend_in!(
                ctx.ms.mipdata_.domains,
                propdomain,
                detach!(ctx.ms.mipdata_.cliquetable, |ct| detach!(
                    ctx.ms.mipdata_.cutpools[cutpool],
                    |cp| ct.separate_cliques(ctx.ms, sol, cp, feastol, randgen, local_nq)
                ))
            );
        }
        if !ctx.ms.parallel_lock_active() {
            ctx.ms.profiling_.stop(self.clique_clock);
        }

        numboundchgs = propagate_and_resolve(ctx, propdomain, status);
        if numboundchgs == -1 {
            return 0;
        } else {
            ncuts += numboundchgs;
        }

        if !is_global_domain(ctx.ms, globaldom, propdomain) {
            let conflictpool = ctx.worker.conflictpool_;
            let pseudocost = ctx.worker.pseudocost_;
            detach!(ctx.ms.mipdata_.pseudocosts[pseudocost], |pc| ctx
                .lp
                .compute_basic_degenerate_duals(
                    ctx.ms,
                    feastol,
                    propdomain,
                    globaldom,
                    conflictpool,
                    pc,
                    true
                ));
        }

        // C++ `HighsTransformedLp transLp(*lp, mipdata.implications, getGlobalDomain())`
        let mut trans_lp = lend_in!(
            ctx.ms.mipdata_.domains,
            propdomain,
            HighsTransformedLp::new(ctx.ms, &*ctx.lp, globaldom)
        );
        if global_domain_infeasible(ctx.ms, globaldom, propdomain) {
            *status = HighsLpRelaxation_Status::kInfeasible;
            return 0;
        }
        let mut lp_aggregator = HighsLpAggregator::new(&*ctx.lp);

        for separator in self.separators.iter_mut() {
            lend_in!(
                ctx.ms.mipdata_.domains,
                propdomain,
                separator.run(ctx.ms, ctx.lp, &mut lp_aggregator, &mut trans_lp, cutpool)
            );
            if global_domain_infeasible(ctx.ms, globaldom, propdomain) {
                *status = HighsLpRelaxation_Status::kInfeasible;
                return 0;
            }
        }

        // user separator hook (mip/hooks.rs); one `is_some` test when no hooks
        if ctx.ms.user_hooks().is_some() {
            crate::mip::hooks::run_user_separator(ctx.ms, &*ctx.lp, &*propdomain, cutpool);
        }

        numboundchgs = propagate_and_resolve(ctx, propdomain, status);
        if numboundchgs == -1 {
            return 0;
        } else {
            ncuts += numboundchgs;
        }

        ncuts + self.separate_cut_pools_and_resolve(ctx, propdomain, status, globaldom, cutpool)
    }

    /// The tail of HiGHS `separationRound`, from `getCutPool().separate(...)` to the end:
    /// collect the violated cuts of the worker's and the global cut pool, add them to the LP and
    /// resolve. Returns the number of cuts added. Split out unchanged so the rusolve pool-only
    /// round can share it.
    fn separate_cut_pools_and_resolve(
        &mut self,
        ctx: &mut SepaCtx<'_>,
        propdomain: &mut HighsDomain,
        status: &mut HighsLpRelaxation_Status,
        globaldom: usize,
        cutpool: usize,
    ) -> HighsInt {
        let feastol = ctx.ms.mipdata_.feastol;
        let mut ncuts: HighsInt = 0;
        let s4 = ctx.ms.options_mip_.mip_cut_selection_hybrid;
        // mip_cut_selection_hybrid: SCIP sepastore.c lines 934-1058 at 32dd425
        if s4 {
            let sol = &ctx.lp.get_solution().col_value;
            let root = is_global_domain(ctx.ms, globaldom, propdomain);
            crate::mip::cut_selection::apply::s4_separate_cutpools(
                ctx.ms,
                sol,
                &*propdomain,
                &mut self.cutset,
                self.s4_randnumgen.as_mut().expect("S4 rng"),
                feastol,
                cutpool,
                root,
            );
        } else {
            let sol = &ctx.lp.get_solution().col_value;
            let user_hooks = ctx.ms.user_hooks(); // rusolve cut-score hook
            detach!(ctx.ms.mipdata_.cutpools[cutpool], |cp| cp.separate_hooked(
                user_hooks.as_ref(),
                sol,
                &*propdomain,
                &mut self.cutset,
                feastol,
                &ctx.ms.mipdata_.cutpools,
                false
            ));
            // Also separate the global cut pool
            if cutpool != 0 {
                detach!(ctx.ms.mipdata_.cutpools[0], |gp| gp.separate_hooked(
                    user_hooks.as_ref(),
                    sol,
                    &*propdomain,
                    &mut self.cutset,
                    feastol,
                    &ctx.ms.mipdata_.cutpools,
                    true
                ));
            }
        }

        if self.cutset.num_cuts() > 0 {
            ncuts += self.cutset.num_cuts();
            ctx.lp.add_cuts(&mut self.cutset);
            *status = lend_in!(
                ctx.ms.mipdata_.workers,
                ctx.worker,
                ctx.lp.resolve_lp(ctx.ms, Some(&mut *propdomain))
            );
            ctx.lp.perform_aging(ctx.ms, true);

            // only for the master domain.
            if is_master_domain(ctx.ms, propdomain)
                && HighsLpRelaxation::unscaled_dual_feasible(*status)
            {
                add_root_redcost(ctx, propdomain);
                if ctx.ms.mipdata_.upper_limit != kHighsInf {
                    propagate_root_redcost(ctx, propdomain);
                }
            }
        }

        ncuts
    }

    /// rusolve mip_node_sepa_pool_only: the round run at a tree node the gate excludes. It keeps
    /// the propagation and the cut-pool part of `separation_round` and skips implied bounds,
    /// cliques, degenerate-dual conflicts, the separators and the user hook, which is where
    /// node separation spends its time (cMIR aggregation and clique extraction from new cuts).
    /// The pool cuts are globally valid, so the node still gets the bound they give, without
    /// the cost of finding new ones.
    fn pool_only_round(
        &mut self,
        ctx: &mut SepaCtx<'_>,
        propdomain: &mut HighsDomain,
        status: &mut HighsLpRelaxation_Status,
    ) -> HighsInt {
        let globaldom = ctx.worker.globaldom_;
        let cutpool = ctx.worker.cutpool_;
        let numboundchgs = propagate_and_resolve(ctx, propdomain, status);
        if numboundchgs == -1 {
            return 0;
        }
        numboundchgs
            + self.separate_cut_pools_and_resolve(ctx, propdomain, status, globaldom, cutpool)
    }

    /// HiGHS: HighsSeparation.cpp::HighsSeparation::separate (lines 181-228)
    /// C++: `void separate(HighsDomain & propdomain)`
    pub fn separate(&mut self, ctx: &mut SepaCtx<'_>, propdomain: &mut HighsDomain) {
        let mut status = ctx.lp.get_status();

        // rusolve mip_node_bigm: cuts built on a node-tightened LP matrix are not globally valid
        if ctx.lp.has_local_mods() {
            return;
        }
        if HighsLpRelaxation::scaled_optimal(status) && !ctx.lp.get_fractional_integers().is_empty()
        {
            // double firstobj = lp->getObjective();
            let firstobj = ctx.ms.mipdata_.rootlpsolobj;
            // rusolve mip_node_sepa_max_rounds: SCIP's separating/maxrounds for local nodes
            // (src/scip/set.c SCIP_DEFAULT_SEPA_MAXROUNDS, read in priceAndCutLoop of
            // src/scip/solve.c, Apache-2.0): the round loop stops after this many rounds. 0 (the
            // default) leaves only the HiGHS stopping rules below.
            let max_rounds = ctx.ms.options_mip_.mip_node_sepa_max_rounds;
            let mut nrounds: HighsInt = 0;

            while ctx.lp.get_objective() < ctx.worker.optimality_limit {
                let lastobj = ctx.lp.get_objective();

                let mut nlpiters = -ctx.lp.get_num_lp_iterations();
                let lp_src_guard = crate::work::lp_kind(crate::work::LP_TREE_SEPA);
                let ncuts = self.separation_round(ctx, propdomain, &mut status);
                drop(lp_src_guard);
                nlpiters += ctx.lp.get_num_lp_iterations();

                if ctx.ms.parallel_lock_active() {
                    ctx.worker.sepa_stats.sepa_lp_iterations += nlpiters;
                } else {
                    ctx.ms.mipdata_.sepa_lp_iterations += nlpiters;
                    ctx.ms.mipdata_.total_lp_iterations += nlpiters;
                }

                if ncuts == 0
                    || !HighsLpRelaxation::scaled_optimal(status)
                    || ctx.lp.get_fractional_integers().is_empty()
                {
                    break;
                }

                nrounds += 1;
                if max_rounds > 0 && nrounds >= max_rounds {
                    break;
                }

                // if the objective improved considerably we continue
                if (ctx.lp.get_objective() - firstobj)
                    <= cpp_max(lastobj - firstobj, ctx.ms.mipdata_.feastol) * 1.01
                {
                    break;
                }
            }
        } else {
            ctx.lp.perform_aging(ctx.ms, true);

            let cutpool = ctx.worker.cutpool_;
            detach!(ctx.ms.mipdata_.cutpools[cutpool], |cp| cp.perform_aging());
        }
    }

    /// HiGHS: HighsSeparation.h::HighsSeparation::setLpRelaxation (lines 30-30)
    /// C++: `void setLpRelaxation(HighsLpRelaxation * lp)`
    // PORT NOTE: `lp` is an index into `ms.mipdata_.lps` (see the struct doc).
    pub fn set_lp_relaxation(&mut self, lp: Option<usize>) {
        self.lp = lp;
    }
}

/// C++ `&propdomain == &mipdata.getDomain()` (`ms.mipdata_.domains[0]`, identity only).
fn is_master_domain(ms: &HighsMipSolver, propdomain: &HighsDomain) -> bool {
    ms.mipdata_.domains[0].holds(propdomain as *const HighsDomain)
}

/// C++ `&propdomain == &mipworker_.getGlobalDomain()` (identity only).
fn is_global_domain(ms: &HighsMipSolver, globaldom: usize, propdomain: &HighsDomain) -> bool {
    ms.mipdata_.domains[globaldom].holds(propdomain as *const HighsDomain)
}

/// C++ `mipworker_.getGlobalDomain().infeasible()`. The global domain may be the detached
/// `propdomain` itself (root), where reading the slot would panic: use `propdomain` then.
fn global_domain_infeasible(
    ms: &HighsMipSolver,
    globaldom: usize,
    propdomain: &HighsDomain,
) -> bool {
    if is_global_domain(ms, globaldom, propdomain) {
        propdomain.infeasible()
    } else {
        ms.mipdata_.domains[globaldom].infeasible()
    }
}

/// C++ `mipdata.redcostfixing.addRootRedcost(mipdata.mipsolver, lp->getSolution().col_dual,
/// lp->getObjective())`. The callee reads `getLp()` / `getDomain()`, so both are lent back.
/// `col_dual` / the objective are cloned before the macro (they alias `ctx.lp`).
fn add_root_redcost(ctx: &mut SepaCtx<'_>, propdomain: &mut HighsDomain) {
    let col_dual = ctx.lp.get_solution().col_dual.clone();
    let objective = ctx.lp.get_objective();
    lend_in!(
        ctx.ms.mipdata_.lps,
        ctx.lp,
        lend_in!(
            ctx.ms.mipdata_.domains,
            propdomain,
            detach!(ctx.ms.mipdata_.redcostfixing, |rc| rc
                .add_root_redcost(ctx.ms, &col_dual, objective))
        )
    );
}

/// C++ `mipdata.redcostfixing.propagateRootRedcost(mipdata.mipsolver)`.
fn propagate_root_redcost(ctx: &mut SepaCtx<'_>, propdomain: &mut HighsDomain) {
    lend_in!(
        ctx.ms.mipdata_.lps,
        ctx.lp,
        lend_in!(
            ctx.ms.mipdata_.domains,
            propdomain,
            detach!(ctx.ms.mipdata_.redcostfixing, |rc| rc
                .propagate_root_redcost(ctx.ms))
        )
    );
}

/// HiGHS: HighsSeparation.cpp::HighsSeparation::separationRound, lambda `propagateAndResolve`
/// (lines 53-101). Returns the number of bound changes, or -1 when the node is infeasible.
fn propagate_and_resolve(
    ctx: &mut SepaCtx<'_>,
    propdomain: &mut HighsDomain,
    status: &mut HighsLpRelaxation_Status,
) -> HighsInt {
    let globaldom = ctx.worker.globaldom_;

    if propdomain.infeasible() || global_domain_infeasible(ctx.ms, globaldom, propdomain) {
        *status = HighsLpRelaxation_Status::kInfeasible;
        propdomain.clear_changed_cols();
        return -1;
    }

    let _ = propdomain.propagate(ctx.ms);
    if propdomain.infeasible() {
        *status = HighsLpRelaxation_Status::kInfeasible;
        propdomain.clear_changed_cols();
        return -1;
    }

    // only modify cliquetable for master worker.
    if is_master_domain(ctx.ms, propdomain) {
        detach!(ctx.ms.mipdata_.cliquetable, |ct| ct
            .cleanup_fixed(ctx.ms, propdomain));
    }

    if global_domain_infeasible(ctx.ms, globaldom, propdomain) {
        *status = HighsLpRelaxation_Status::kInfeasible;
        propdomain.clear_changed_cols();
        return -1;
    }

    let num_bound_chgs = propdomain.get_changed_cols().len() as HighsInt;

    while !propdomain.get_changed_cols().is_empty() {
        let upper_limit = ctx.worker.upper_limit;
        ctx.lp.set_objective_limit(ctx.ms, upper_limit);
        *status = lend_in!(
            ctx.ms.mipdata_.workers,
            ctx.worker,
            ctx.lp.resolve_lp(ctx.ms, Some(&mut *propdomain))
        );
        if !HighsLpRelaxation::scaled_optimal(*status) {
            return -1;
        }

        if is_master_domain(ctx.ms, propdomain)
            && HighsLpRelaxation::unscaled_dual_feasible(*status)
        {
            add_root_redcost(ctx, propdomain);
            if ctx.worker.upper_limit != kHighsInf {
                propagate_root_redcost(ctx, propdomain);
            }
        }
    }

    num_bound_chgs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x1_mss_b_separator_registered_only_with_flag() {
        let mut ms = HighsMipSolver::default();
        ms.options_mip_.mip_ti_precedence_sepa = false;
        ms.options_mip_.mip_ti_list_sched = false;
        ms.options_mip_.mip_netdesign_sepa = false;
        ms.options_mip_.mip_netdesign_heur = false;
        ms.options_mip_.mip_netdesign_cutset = false;
        // Off: the 1.15.1 separator vector, three entries.
        let sepa = HighsSeparation::new(&ms);
        assert_eq!(sepa.separators.len(), 3);
        assert!(
            !sepa
                .separators
                .iter()
                .any(|s| matches!(s, HighsSeparatorBox::MachineSched(_)))
        );

        // On: the machine-scheduling separator is pushed last, clock 999.
        ms.options_mip_.mip_machine_sched_sepa = true;
        let sepa = HighsSeparation::new(&ms);
        assert_eq!(sepa.separators.len(), 4);
        assert!(matches!(
            sepa.separators[3],
            HighsSeparatorBox::MachineSched(_)
        ));
        assert_eq!(sepa.separators[3].base().get_clock_index(), 999);
    }

    /// Small binary knapsack, min -value.x s.t. weight.x <= 200 (12 items).
    fn s4_gate_knapsack_lp() -> HighsLp {
        let weights: [f64; 12] = [
            10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0,
        ];
        let values: [f64; 12] = [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ];
        let n = weights.len();
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = 1;
        lp.col_cost_ = values.iter().map(|v| -v).collect();
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![200.0];
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = (0..=n).map(|j| j as HighsInt).collect();
        lp.a_matrix_.index_ = vec![0; n];
        lp.a_matrix_.value_ = weights.to_vec();
        lp
    }

    /// Brute-force optimum of [`s4_gate_knapsack_lp`] (min -value.x).
    fn s4_gate_brute_force() -> f64 {
        let weights = [
            10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0,
        ];
        let values = [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ];
        let mut best = 0.0f64;
        for mask in 0u32..(1 << 12) {
            let (mut w, mut v) = (0.0f64, 0.0f64);
            for i in 0..12 {
                if mask & (1 << i) != 0 {
                    w += weights[i];
                    v += values[i];
                }
            }
            if w <= 200.0 && v > best {
                best = v;
            }
        }
        -best
    }

    /// Solves the knapsack with the given explicit flag (None: default options).
    /// Returns (objective, num_nodes, total_lp_iterations).
    fn s4_gate_run(hybrid: Option<bool>) -> (f64, i64, i64) {
        let handle = std::thread::spawn(move || {
            crate::parallel::highs_parallel::initialize_scheduler(1);
            let lp = s4_gate_knapsack_lp();
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.presolve = kHighsOffString.to_string();
            options.mip_heuristic_effort = 0.0;
            options.mip_heuristic_run_feasibility_jump = false;
            options.mip_heuristic_run_rens = false;
            options.mip_heuristic_run_rins = false;
            if let Some(flag) = hybrid {
                options.mip_cut_selection_hybrid = flag;
            }
            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let out = (
                ms.solution_objective_,
                ms.mipdata_.num_nodes,
                ms.mipdata_.total_lp_iterations,
            );
            drop(ms);
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
            out
        });
        handle.join().unwrap()
    }

    #[test]
    fn s4_gate_off_matches_1_15_1() {
        // Default options and an explicit `false` must be bit-identical (rule F).
        let (obj_default, nodes_default, iters_default) = s4_gate_run(None);
        let (obj_off, nodes_off, iters_off) = s4_gate_run(Some(false));
        assert_eq!(
            obj_default.to_bits(),
            obj_off.to_bits(),
            "objective bits differ"
        );
        assert_eq!(nodes_default, nodes_off, "num_nodes differ");
        assert_eq!(iters_default, iters_off, "total_lp_iterations differ");
    }

    #[test]
    fn s4_gate_on_same_optimum() {
        let expected = s4_gate_brute_force();
        let (obj, _, _) = s4_gate_run(Some(true));
        assert!(
            (obj - expected).abs() < 1e-6,
            "flag on: objective {} != brute force {}",
            obj,
            expected
        );
    }

    #[test]
    fn s4_sepafield_off_is_none() {
        let mut ms = HighsMipSolver::default();
        ms.options_mip_.mip_ti_precedence_sepa = false;
        ms.options_mip_.mip_ti_list_sched = false;
        ms.options_mip_.mip_netdesign_sepa = false;
        ms.options_mip_.mip_netdesign_heur = false;
        ms.options_mip_.mip_netdesign_cutset = false;
        let sepa = HighsSeparation::new(&ms);
        assert!(sepa.s4_randnumgen.is_none());
        assert_eq!(sepa.separators.len(), 3);
    }

    #[test]
    fn s4_sepafield_on_is_seeded() {
        let mut ms = HighsMipSolver::default();
        ms.options_mip_.mip_cut_selection_hybrid = true;
        // random_seed default is 0, so the shift is 0 and the seed is 0x5EED
        let sepa = HighsSeparation::new(&ms);
        let rng = sepa.s4_randnumgen.as_ref().unwrap();
        let expected = crate::mip::cut_selection::S4RandNumGen::new(0x5EED);
        assert_eq!(rng.seed, expected.seed);
        assert_eq!(rng.xor_seed, expected.xor_seed);
        assert_eq!(rng.mwc_seed, expected.mwc_seed);
        assert_eq!(rng.cst_seed, expected.cst_seed);
    }

    #[test]
    fn s4_sepafield_on_random_seed_shift() {
        let mut ms = HighsMipSolver::default();
        ms.options_mip_.mip_cut_selection_hybrid = true;
        ms.options_mip_.random_seed = 1;
        // shift = 1, multiplier = 10: offset (1 + 6*10) * 1 = 61
        let sepa = HighsSeparation::new(&ms);
        let rng = sepa.s4_randnumgen.as_ref().unwrap();
        let expected = crate::mip::cut_selection::S4RandNumGen::new(0x5EED + 61);
        assert_eq!(rng.seed, expected.seed);
        assert_eq!(rng.xor_seed, expected.xor_seed);
        assert_eq!(rng.mwc_seed, expected.mwc_seed);
        assert_eq!(rng.cst_seed, expected.cst_seed);
    }
}
