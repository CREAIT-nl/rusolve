//! The deterministic work counter.
//!
//! HiGHS has no "work" or "effort" option for its solve, but its dual and primal simplex keep a
//! deterministic cost clock, the *synthetic tick* (`HEkk::total_synthetic_tick_`,
//! `HFactor::build_synthetic_tick`, `HVector::synthetic_tick`; HFactor.cpp:118, 691, 816, 875,
//! 1368, 1610, 1832, 1931, 1991, 2027; HEkkDual.cpp:2198-2199; HEkkPrimal.cpp:1471, 1545-1546;
//! HEkkDualMulti.cpp:666-667). HiGHS uses it to decide when to reinvert; it counts the operations
//! of each FTRAN, BTRAN, update and INVERT weighted by the nonzeros touched, so it does not depend
//! on the machine, the compiler or the number of threads.
//!
//! The port keeps that clock per `HEkk` and resets it at every INVERT. This module adds one
//! monotone total of the same ticks, fed from the same lines that feed `total_synthetic_tick_`.
//! The search never reads it (reporting only), so it cannot move a trace.
//!
//! The totals are per solve, not per process: `HighsMipSolver::run` installs a [`WorkCounters`] in a thread-local slot, the
//! scheduler wrapper in `parallel::highs_parallel::erase_task_lifetime` re-installs it around
//! every task the solve spawns, and the increment sites add to whatever is installed. Work with
//! no solve installed goes to a process-wide total.
//!
//! Ticks are integer-valued sums of counts times small integer weights, added as `u64`, so the
//! total does not depend on the order in which parallel workers add to it.

use std::cell::Cell;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

/// The four work counters of one solve (or of the process, for [`WorkCounters::global`]).
#[derive(Default, Debug)]
pub struct WorkCounters {
    simplex_ticks: AtomicU64,
    propagation_units: AtomicU64,
    fj_units: AtomicU64,
    presolve_units: AtomicU64,
    /// Reporting-only event tallies (see [`Event`]); not part of `work`.
    events: [AtomicU64; Event::COUNT],
    /// Reporting-only LP tallies by source (see [`LpSrc`]): simplex iterations and LP solves.
    lp_iters: [AtomicU64; LP_SRC_COUNT],
    lp_solves: [AtomicU64; LP_SRC_COUNT],
}

/// Reporting-only attribution of every LP relaxation solve: which part
/// of the search asked for it. Set by [`LpSrcGuard`]s around the callers; read by nobody in the
/// search. A heuristic guard (outermost wins) overrides the kind guards inside it, so a dive or
/// sub-MIP started by a heuristic is charged to that heuristic.
pub type LpSrc = u8;
pub const LP_UNTAGGED: LpSrc = 0;
pub const LP_ROOT: LpSrc = 1;
pub const LP_ROOT_SEPA: LpSrc = 2;
pub const LP_TREE_SEPA: LpSrc = 3;
pub const LP_STRONG_BRANCH: LpSrc = 4;
/// First node-resolve source; `LP_NODE0 + origin * 4 + depth_bucket`.
pub const LP_NODE0: LpSrc = 5;
pub const LP_ORIGIN_CHILD: u8 = 0;
pub const LP_ORIGIN_SIBLING: u8 = 1;
pub const LP_ORIGIN_JUMP: u8 = 2;
pub const LP_H_RENS: LpSrc = 17;
pub const LP_H_RINS: LpSrc = 18;
pub const LP_H_ROOT_REDCOST: LpSrc = 19;
pub const LP_H_FEAS_PUMP: LpSrc = 20;
pub const LP_H_CENTRAL: LpSrc = 21;
pub const LP_H_ROUNDING: LpSrc = 22;
pub const LP_H_GUIDED_DIVE: LpSrc = 23;
pub const LP_H_SCIP_SUBMIP: LpSrc = 24;
pub const LP_H_LS_REPAIR: LpSrc = 25;
pub const LP_SRC_COUNT: usize = 26;
pub const LP_SRC_NAMES: [&str; LP_SRC_COUNT] = [
    "untagged",
    "root_lp",
    "root_separation",
    "tree_separation",
    "strong_branching",
    "node_child_depth_0_5",
    "node_child_depth_6_15",
    "node_child_depth_16_30",
    "node_child_depth_31_plus",
    "node_sibling_depth_0_5",
    "node_sibling_depth_6_15",
    "node_sibling_depth_16_30",
    "node_sibling_depth_31_plus",
    "node_jump_depth_0_5",
    "node_jump_depth_6_15",
    "node_jump_depth_16_30",
    "node_jump_depth_31_plus",
    "heur_rens",
    "heur_rins",
    "heur_root_reduced_cost",
    "heur_feasibility_pump",
    "heur_central_rounding",
    "heur_rounding_shifting",
    "heur_guided_dive",
    "heur_scip_submip",
    "heur_ls_int_lp_repair",
];

/// The node-resolve source for a node of `depth` reached by `origin` (`LP_ORIGIN_*`).
pub fn lp_node_src(origin: u8, depth: usize) -> LpSrc {
    let b = match depth {
        0..=5 => 0,
        6..=15 => 1,
        16..=30 => 2,
        _ => 3,
    };
    LP_NODE0 + origin * 4 + b
}

thread_local! {
    static LP_HEUR: Cell<LpSrc> = const { Cell::new(0) };
    static LP_KIND: Cell<LpSrc> = const { Cell::new(0) };
}

/// Sets the LP source of the current thread until dropped.
pub struct LpSrcGuard {
    heur: bool,
    prev: LpSrc,
}

/// A kind guard (innermost wins; a heuristic guard outside it still takes precedence).
pub fn lp_kind(src: LpSrc) -> LpSrcGuard {
    LpSrcGuard {
        heur: false,
        prev: LP_KIND.with(|c| c.replace(src)),
    }
}

/// A heuristic guard (outermost wins: a no-op inside another heuristic).
pub fn lp_heur(src: LpSrc) -> LpSrcGuard {
    LpSrcGuard {
        heur: true,
        prev: LP_HEUR.with(|c| {
            let p = c.get();
            if p == 0 {
                c.set(src);
            }
            p
        }),
    }
}

impl Drop for LpSrcGuard {
    fn drop(&mut self) {
        if self.heur {
            LP_HEUR.with(|c| c.set(self.prev));
        } else {
            LP_KIND.with(|c| c.set(self.prev));
        }
    }
}

/// Records one LP solve of `iters` simplex iterations against the current source.
#[inline]
pub fn add_lp_solve(iters: u64) {
    let h = LP_HEUR.with(Cell::get);
    let src = if h != 0 { h } else { LP_KIND.with(Cell::get) } as usize;
    with_current(|c| {
        c.lp_iters[src].fetch_add(iters, Ordering::Relaxed);
        c.lp_solves[src].fetch_add(1, Ordering::Relaxed);
    });
}

/// Reporting-only event tallies of one solve. They never feed `work()` and
/// the search never reads them.
#[derive(Clone, Copy, Debug)]
#[repr(usize)]
pub enum Event {
    /// `obj_prop_propagate` calls that passed `obj_prop_should_be_propagated`.
    ObjPropRuns,
    /// Bound changes pushed by `obj_prop_propagate` itself.
    ObjPropBoundChanges,
    /// `obj_prop_propagate` calls that proved the domain infeasible.
    ObjPropCutoffs,
    /// `HighsRedcostFixing::propagate_red_cost` calls (local, at a node).
    RedcostLocalRuns,
    /// Bound changes the local reduced-cost fixing pushed on the node domain.
    RedcostLocalBoundChanges,
    /// Local reduced-cost fixing calls that left the node domain infeasible.
    RedcostLocalCutoffs,
    /// Root reduced-cost fixing calls with lurking bounds.
    RedcostRootRuns,
    /// Domain-change stack growth of the root domain across root reduced-cost fixing
    /// (includes the propagation it triggers).
    RedcostRootBoundChanges,
    /// Conflict analysis, source Tree: `conflict_analysis*` calls that passed the entry guards.
    ConfCallsTree,
    /// Conflict analysis, source Heur: `conflict_analysis*` calls that passed the entry guards.
    ConfCallsHeur,
    /// Conflict analysis, source InfLp: `conflict_analysis*` calls that passed the entry guards.
    ConfCallsInfLp,
    /// Conflict analysis, source BoundLp: `conflict_analysis*` calls that passed the entry guards.
    ConfCallsBoundLp,
    /// Conflict analysis, source Degen: `conflict_analysis*` calls that passed the entry guards.
    ConfCallsDegen,
    /// Conflict analysis, source Tree: calls whose infeasibility explanation succeeded (a conflict set was built).
    ConfExplainedTree,
    /// Conflict analysis, source Heur: calls whose infeasibility explanation succeeded (a conflict set was built).
    ConfExplainedHeur,
    /// Conflict analysis, source InfLp: calls whose infeasibility explanation succeeded (a conflict set was built).
    ConfExplainedInfLp,
    /// Conflict analysis, source BoundLp: calls whose infeasibility explanation succeeded (a conflict set was built).
    ConfExplainedBoundLp,
    /// Conflict analysis, source Degen: calls whose infeasibility explanation succeeded (a conflict set was built).
    ConfExplainedDegen,
    /// Conflict analysis, source Tree: conflict constraints added to the pool.
    ConfAddedTree,
    /// Conflict analysis, source Heur: conflict constraints added to the pool.
    ConfAddedHeur,
    /// Conflict analysis, source InfLp: conflict constraints added to the pool.
    ConfAddedInfLp,
    /// Conflict analysis, source BoundLp: conflict constraints added to the pool.
    ConfAddedBoundLp,
    /// Conflict analysis, source Degen: conflict constraints added to the pool.
    ConfAddedDegen,
    /// Conflict analysis, source Tree: literals in the conflict constraints added to the pool.
    ConfLenTree,
    /// Conflict analysis, source Heur: literals in the conflict constraints added to the pool.
    ConfLenHeur,
    /// Conflict analysis, source InfLp: literals in the conflict constraints added to the pool.
    ConfLenInfLp,
    /// Conflict analysis, source BoundLp: literals in the conflict constraints added to the pool.
    ConfLenBoundLp,
    /// Conflict analysis, source Degen: literals in the conflict constraints added to the pool.
    ConfLenDegen,
    /// Bound changes the conflict pool's propagation pushed on a domain.
    ConfPropBoundChanges,
    /// Domains the conflict pool's propagation proved infeasible.
    ConfPropCutoffs,
    /// Linear conflict cuts handed to the cut pool by `generate_conflict`.
    ConfCutsAdded,
    /// Conflicts deleted from a pool (aging).
    ConfRemoved,
}

impl Event {
    pub const COUNT: usize = 32;
    pub const NAMES: [&'static str; Event::COUNT] = [
        "objprop_runs",
        "objprop_bound_changes",
        "objprop_cutoffs",
        "redcost_local_runs",
        "redcost_local_bound_changes",
        "redcost_local_cutoffs",
        "redcost_root_runs",
        "redcost_root_bound_changes",
        "conf_calls_tree",
        "conf_calls_heur",
        "conf_calls_inflp",
        "conf_calls_boundlp",
        "conf_calls_degen",
        "conf_explained_tree",
        "conf_explained_heur",
        "conf_explained_inflp",
        "conf_explained_boundlp",
        "conf_explained_degen",
        "conf_added_tree",
        "conf_added_heur",
        "conf_added_inflp",
        "conf_added_boundlp",
        "conf_added_degen",
        "conf_len_tree",
        "conf_len_heur",
        "conf_len_inflp",
        "conf_len_boundlp",
        "conf_len_degen",
        "conf_prop_bound_changes",
        "conf_prop_cutoffs",
        "conf_cuts_added",
        "conf_removed",
    ];
}

static GLOBAL: WorkCounters = WorkCounters {
    simplex_ticks: AtomicU64::new(0),
    propagation_units: AtomicU64::new(0),
    fj_units: AtomicU64::new(0),
    presolve_units: AtomicU64::new(0),
    events: [const { AtomicU64::new(0) }; Event::COUNT],
    lp_iters: [const { AtomicU64::new(0) }; LP_SRC_COUNT],
    lp_solves: [const { AtomicU64::new(0) }; LP_SRC_COUNT],
};

thread_local! {
    /// The counters of the solve this thread is working for (null: none). Const-initialised, so
    /// an increment is one TLS read, a null test and one relaxed add.
    static CURRENT: Cell<*const WorkCounters> = const { Cell::new(std::ptr::null()) };
}

impl WorkCounters {
    /// The process-wide counters: they receive the work done outside any solve scope and the
    /// totals of every finished solve scope (the four work counters; the reporting-only event
    /// tallies stay per solve and are not folded in).
    pub fn global() -> &'static WorkCounters {
        &GLOBAL
    }

    /// The event tallies now (see [`Event`]).
    pub fn events(&self) -> [u64; Event::COUNT] {
        std::array::from_fn(|i| self.events[i].load(Ordering::Relaxed))
    }

    /// LP simplex iterations and solves by source now (see [`LpSrc`]).
    pub fn lp_by_source(&self) -> ([u64; LP_SRC_COUNT], [u64; LP_SRC_COUNT]) {
        (
            std::array::from_fn(|i| self.lp_iters[i].load(Ordering::Relaxed)),
            std::array::from_fn(|i| self.lp_solves[i].load(Ordering::Relaxed)),
        )
    }

    /// The counters now.
    pub fn snapshot(&self) -> WorkSnapshot {
        WorkSnapshot {
            simplex_ticks: self.simplex_ticks.load(Ordering::Relaxed),
            propagation_units: self.propagation_units.load(Ordering::Relaxed),
            fj_units: self.fj_units.load(Ordering::Relaxed),
            presolve_units: self.presolve_units.load(Ordering::Relaxed),
        }
    }
}

/// Runs `f` on the counters of the current thread's solve, or the process-wide ones.
#[inline(always)]
fn with_current<R>(f: impl FnOnce(&WorkCounters) -> R) -> R {
    let p = CURRENT.with(Cell::get);
    if p.is_null() {
        f(&GLOBAL)
    } else {
        // SAFETY: a non-null CURRENT is set only by `SolveScope` / `TaskScope`, which hold an
        // `Arc` to the pointee for as long as the pointer is installed, and reset it before
        // releasing that `Arc`.
        f(unsafe { &*p })
    }
}

/// The current thread's solve counters, if a solve is installed (a new strong reference).
pub fn current() -> Option<Arc<WorkCounters>> {
    let p = CURRENT.with(Cell::get);
    if p.is_null() {
        None
    } else {
        // SAFETY: `p` is `Arc::as_ptr` of an `Arc` that the installing scope still holds (see
        // `with_current`), so the strong count is at least 1 and `increment_strong_count` +
        // `from_raw` yields one additional owned reference.
        unsafe {
            Arc::increment_strong_count(p);
            Some(Arc::from_raw(p))
        }
    }
}

/// Installs `counters` as the current thread's solve counters until dropped, then restores the
/// previous ones. Used by the scheduler wrapper around each spawned task, so work done by a
/// worker thread is charged to the solve that spawned the task (a worker that steals another
/// solve's task while waiting switches for the duration of that task).
pub struct TaskScope {
    prev: *const WorkCounters,
    _hold: Option<Arc<WorkCounters>>,
}

impl TaskScope {
    /// Installs `counters` (`None`: the process-wide counters).
    pub fn enter(counters: Option<Arc<WorkCounters>>) -> TaskScope {
        let p = counters.as_ref().map_or(std::ptr::null(), Arc::as_ptr);
        let prev = CURRENT.with(|c| c.replace(p));
        TaskScope {
            prev,
            _hold: counters,
        }
    }
}

impl Drop for TaskScope {
    fn drop(&mut self) {
        CURRENT.with(|c| c.set(self.prev));
    }
}

/// One solve's counting scope: [`SolveScope::enter`] at the top of `HighsMipSolver::run`.
/// A nested solve (sub-MIP, heuristic) joins the enclosing scope, so its work is part of the
/// parent's. The outermost scope's totals are added to the process-wide total when it ends.
pub struct SolveScope {
    counters: Arc<WorkCounters>,
    owner: Option<TaskScope>,
}

impl SolveScope {
    pub fn enter() -> SolveScope {
        match current() {
            Some(counters) => SolveScope {
                counters,
                owner: None,
            },
            None => {
                let counters = Arc::new(WorkCounters::default());
                SolveScope {
                    owner: Some(TaskScope::enter(Some(counters.clone()))),
                    counters,
                }
            }
        }
    }

    /// The counters this scope reports.
    pub fn counters(&self) -> &Arc<WorkCounters> {
        &self.counters
    }
}

impl Drop for SolveScope {
    fn drop(&mut self) {
        if self.owner.take().is_some() {
            let s = self.counters.snapshot();
            GLOBAL
                .simplex_ticks
                .fetch_add(s.simplex_ticks, Ordering::Relaxed);
            GLOBAL
                .propagation_units
                .fetch_add(s.propagation_units, Ordering::Relaxed);
            GLOBAL.fj_units.fetch_add(s.fj_units, Ordering::Relaxed);
            GLOBAL
                .presolve_units
                .fetch_add(s.presolve_units, Ordering::Relaxed);
        }
    }
}

/// Source of the conflict analysis running on this thread (reporting only).
#[derive(Clone, Copy, Debug)]
#[repr(u8)]
pub enum ConflictSrc {
    /// Propagation infeasibility in the tree (node propagation, strong-branching propagation).
    Tree = 0,
    /// Propagation infeasibility inside a primal heuristic.
    Heur = 1,
    /// Infeasible node LP, dual ray proof.
    InfLp = 2,
    /// Node LP cut off by the incumbent, dual solution proof.
    BoundLp = 3,
    /// Reconvergence cuts from degenerate duals and reduced costs.
    Degen = 4,
}

thread_local! {
    static CONFLICT_SRC: Cell<u8> = const { Cell::new(0) };
}

/// Restores the previous conflict source when dropped.
pub struct ConflictSrcGuard(u8);

impl Drop for ConflictSrcGuard {
    fn drop(&mut self) {
        CONFLICT_SRC.with(|c| c.set(self.0));
    }
}

/// Tags the conflict analysis run by this thread until the guard is dropped.
pub fn conflict_src(src: ConflictSrc) -> ConflictSrcGuard {
    ConflictSrcGuard(CONFLICT_SRC.with(|c| c.replace(src as u8)))
}

/// True when the ablation lever `mip_conflict_off_mask` switches off the source this thread is
/// tagged with.
#[inline]
pub fn conflict_src_masked(mask: i32) -> bool {
    mask != 0 && (mask >> CONFLICT_SRC.with(Cell::get)) & 1 != 0
}

/// Which tally of the conflict analysis to bump (see [`conflict_event`]).
#[derive(Clone, Copy, Debug)]
#[repr(usize)]
pub enum ConflictKind {
    Calls = 0,
    Explained = 1,
    Added = 2,
    Len = 3,
}

/// Adds `n` to the tally `kind` of the conflict source this thread is tagged with.
#[inline]
pub fn conflict_event(kind: ConflictKind, n: u64) {
    let src = CONFLICT_SRC.with(Cell::get) as usize;
    let idx = Event::ConfCallsTree as usize + (kind as usize) * 5 + src;
    if n != 0 {
        with_current(|c| c.events[idx].fetch_add(n, Ordering::Relaxed));
    }
}

/// Adds `n` to the event tally `e` of the current solve (reporting only).
#[inline]
pub fn add_event(e: Event, n: u64) {
    if n != 0 {
        with_current(|c| c.events[e as usize].fetch_add(n, Ordering::Relaxed));
    }
}

/// Adds `ticks` simplex synthetic ticks to the current solve's total.
#[inline]
pub fn add_simplex_ticks(ticks: f64) {
    with_current(|c| c.simplex_ticks.fetch_add(ticks as u64, Ordering::Relaxed));
}

/// The simplex synthetic ticks so far of the current thread's solve (the process-wide total
/// when no solve is installed). Take a difference around a region of one solve.
#[inline]
pub fn simplex_ticks() -> u64 {
    with_current(|c| c.simplex_ticks.load(Ordering::Relaxed))
}

/// Adds the nonzeros of a row walked by `HPresolve::getRowPositions` (HPresolve.cpp), the
/// presolve loops' row traversal.
#[inline]
pub fn add_presolve_units(nnz: u64) {
    with_current(|c| c.presolve_units.fetch_add(nnz, Ordering::Relaxed));
}

/// Presolve row-walk units so far (see [`add_presolve_units`]).
#[inline]
pub fn presolve_units() -> u64 {
    with_current(|c| c.presolve_units.load(Ordering::Relaxed))
}

/// Adds `nnz` column nonzeros walked by a `HighsDomain` activity update (HighsDomain.cpp
/// `updateActivityLbChange` / `updateActivityUbChange`: one pass over the column per bound change,
/// the inner loop of propagation, probing and dual-fix probing).
#[inline]
pub fn add_propagation_units(nnz: u64) {
    with_current(|c| c.propagation_units.fetch_add(nnz, Ordering::Relaxed));
}

/// Adds the FeasibilityJump `totalEffort` of one call (HighsFeasibilityJump.cpp:107-108; feasibilityjump.hh:459,
/// 590, 604, 704, 726, 761): HiGHS's own deterministic effort count for that heuristic, in
/// nonzeros touched.
#[inline]
pub fn add_fj_units(effort: u64) {
    with_current(|c| c.fj_units.fetch_add(effort, Ordering::Relaxed));
}

/// Propagation units so far (see [`add_propagation_units`]).
#[inline]
pub fn propagation_units() -> u64 {
    with_current(|c| c.propagation_units.load(Ordering::Relaxed))
}

/// FeasibilityJump effort units so far (see [`add_fj_units`]).
#[inline]
pub fn fj_units() -> u64 {
    with_current(|c| c.fj_units.load(Ordering::Relaxed))
}

/// Ticks per propagation, presolve or FeasibilityJump unit in `work`. Empirical: the weight that
/// maximises the log-log correlation with CPU time on our benchmark set.
pub const NONSIMPLEX_UNIT_TICKS: u64 = 16;

/// The three counters at one moment; take a difference around a solve.
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct WorkSnapshot {
    pub simplex_ticks: u64,
    pub propagation_units: u64,
    pub fj_units: u64,
    pub presolve_units: u64,
}

impl WorkSnapshot {
    /// The counters of the current thread's solve now (process-wide ones outside a solve).
    pub fn now() -> Self {
        WorkSnapshot {
            simplex_ticks: simplex_ticks(),
            propagation_units: propagation_units(),
            fj_units: fj_units(),
            presolve_units: presolve_units(),
        }
    }

    /// Counters since `start`.
    pub fn since(&self, start: &WorkSnapshot) -> WorkSnapshot {
        WorkSnapshot {
            simplex_ticks: self.simplex_ticks.saturating_sub(start.simplex_ticks),
            propagation_units: self
                .propagation_units
                .saturating_sub(start.propagation_units),
            fj_units: self.fj_units.saturating_sub(start.fj_units),
            presolve_units: self.presolve_units.saturating_sub(start.presolve_units),
        }
    }

    /// `work`: simplex ticks plus the other units at `NONSIMPLEX_UNIT_TICKS` ticks each.
    pub fn work(&self) -> u64 {
        self.simplex_ticks
            + NONSIMPLEX_UNIT_TICKS * (self.propagation_units + self.fj_units + self.presolve_units)
    }
}

#[cfg(test)]
mod conflict_counter_tests {
    use super::*;

    /// `conflict_event` writes slot `ConfCallsTree + kind * 5 + src`; the slot's name must be
    /// `conf_<kind>_<src>` so the reported counters match the source that was tagged.
    #[test]
    fn conflict_event_index_matches_event_names() {
        let kinds = [
            (ConflictKind::Calls, "calls"),
            (ConflictKind::Explained, "explained"),
            (ConflictKind::Added, "added"),
            (ConflictKind::Len, "len"),
        ];
        let srcs = [
            (ConflictSrc::Tree, "tree"),
            (ConflictSrc::Heur, "heur"),
            (ConflictSrc::InfLp, "inflp"),
            (ConflictSrc::BoundLp, "boundlp"),
            (ConflictSrc::Degen, "degen"),
        ];
        assert_eq!(Event::ConfCallsTree as usize, 8);
        let scope = SolveScope::enter();
        for (ki, (kind, kname)) in kinds.iter().enumerate() {
            for (si, (src, sname)) in srcs.iter().enumerate() {
                let before = scope.counters().events();
                {
                    let _g = conflict_src(*src);
                    conflict_event(*kind, 3);
                }
                let after = scope.counters().events();
                let idx = 8 + ki * 5 + si;
                assert_eq!(Event::NAMES[idx], format!("conf_{kname}_{sname}"));
                for i in 0..Event::COUNT {
                    let want = if i == idx { 3 } else { 0 };
                    assert_eq!(after[i] - before[i], want, "slot {i} for {kname}/{sname}");
                }
            }
        }
    }

    #[test]
    fn conflict_src_masked_reads_the_tagged_bit() {
        assert!(!conflict_src_masked(0));
        for (bit, src) in [
            ConflictSrc::Tree,
            ConflictSrc::Heur,
            ConflictSrc::InfLp,
            ConflictSrc::BoundLp,
            ConflictSrc::Degen,
        ]
        .into_iter()
        .enumerate()
        {
            let _g = conflict_src(src);
            assert!(conflict_src_masked(1 << bit));
            assert!(!conflict_src_masked(31 ^ (1 << bit)));
        }
    }
}
