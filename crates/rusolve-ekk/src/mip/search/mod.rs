//! HiGHS source: mip/HighsSearch.h (class definitions and header-inline methods)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::mip::symmetry_dynamic::{S3PathNode, S3Stats};
use crate::prelude::*;
use std::sync::Arc;

mod part1;
pub use part1::*;
mod part2;
pub use part2::*;

/// HiGHS: HighsSearch.h::HighsSearch (lines 32-273)
///
/// Removed C++ members:
/// * C++: `HighsMipWorker & mipworker`: removed.
/// * C++: `const HighsMipSolver & mipsolver`: removed.
/// * C++ `HighsLpRelaxation* lp`: removed; `SearchCtx::lp`.
/// * C++ `HighsPseudocost& pseudocost`: removed; `SearchCtx::pseudocost`.
#[derive(Clone, Default, Debug)]
pub struct HighsSearch {
    /// C++ `HighsDomain localdom`.
    pub localdom: HighsDomain,
    /// C++: `HighsRandom random`
    pub random: HighsRandom,
    /// C++: `int64_t nnodes`
    pub nnodes: i64,
    /// C++: `int64_t nleaves`
    pub nleaves: i64,
    /// C++: `int64_t lpiterations`
    pub lpiterations: i64,
    /// C++: `int64_t heurlpiterations`
    pub heurlpiterations: i64,
    /// C++: `int64_t sblpiterations`
    pub sblpiterations: i64,
    /// rusolve, reporting only: how the node about to be evaluated was
    /// reached (`work::LP_ORIGIN_*`); read once by `evaluate_node`, then reset to child (the default: a node made by branching).
    pub eval_origin: u8,
    /// C++: `double upper_limit`
    pub upper_limit: f64,
    /// C++: `HighsCDouble treeweight`
    pub treeweight: HighsCDouble,
    /// C++: `std::vector<HighsInt> inds`
    pub inds: Vec<HighsInt>,
    /// C++: `std::vector<double> vals`
    pub vals: Vec<f64>,
    /// C++: `HighsInt depthoffset`
    pub depthoffset: HighsInt,
    /// C++: `bool inbranching`
    pub inbranching: bool,
    /// C++: `bool inheuristic`
    pub inheuristic: bool,
    /// C++: `bool countTreeWeight`
    pub count_tree_weight: bool,
    /// C++: `ChildSelectionRule childselrule`
    pub childselrule: HighsSearch_ChildSelectionRule,
    /// C++: `std::vector<double> subrootsol`
    pub subrootsol: Vec<f64>,
    /// C++: `std::vector<NodeData> nodestack`
    pub nodestack: Vec<HighsSearch_NodeData>,
    /// C++: `StabilizerOrbitWorkspace stabilizerOrbitWorkspace`
    pub stabilizer_orbit_workspace: StabilizerOrbitWorkspace,
    /// C++: `HighsHashTable<HighsInt, int> reliableatnode`
    pub reliableatnode: HighsHashTable<HighsInt, i32>,
    /// rusolve mip_symmetry_dynamic: dynamic-symmetry statistics for this search; no C++ counterpart.
    pub(crate) s3_stats: S3Stats,
    /// rusolve `mip_plunge_quot`: set by `dive` when it stops because the new child's lower
    /// bound is above the plunge bound; the main loop then ends the plunge. No C++ counterpart.
    pub(crate) plunge_bound_stop: bool,
}

/// HiGHS: C++ plain `enum ReliableFlags` nested in `HighsSearch` (values are associated consts).
pub type HighsSearch_ReliableFlags = HighsInt;
impl HighsSearch {
    pub const kUpReliable: HighsInt = 1;
    pub const kDownReliable: HighsInt = 2;
    pub const kReliable: HighsInt = 3;
}

/// HiGHS: HighsSearch.h::HighsSearch::ChildSelectionRule (enum class)
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsSearch_ChildSelectionRule {
    #[default]
    kUp = 0,
    kDown = 1,
    kRootSol = 2,
    kObj = 3,
    kRandom = 4,
    kBestCost = 5,
    kWorstCost = 6,
    kDisjunction = 7,
    kHybridInferenceCost = 8,
}

/// HiGHS: HighsSearch.h::HighsSearch::NodeResult (enum class)
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsSearch_NodeResult {
    #[default]
    kBoundExceeding = 0,
    kDomainInfeasible = 1,
    kLpInfeasible = 2,
    kBranched = 3,
    kSubOptimal = 4,
    kOpen = 5,
}

/// The C++ `HighsNodeQueue& nodequeue` argument of `currentNodeToQueue` /
/// `openNodesToQueue` / `backtrackPlunge`. Every C++ call site passes
/// `parallelLockActive() ? worker.nodequeue : mipdata_->nodequeue`; both stay
/// attached in `ms` (conflict analysis inside these methods reads
/// `mipdata_->nodequeue`), so the method reaches the queue through
/// [`SearchCtx::node_queue`] at each use instead of holding `&mut`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum HighsSearch_NodeQueueRef {
    /// `ms.mipdata_.nodequeue`
    Global,
    /// `ms.mipdata_.workers[ctx.worker].nodequeue`
    Worker,
}

impl SearchCtx<'_> {
    /// The queue named by `q` (borrow it only for the `emplace_node` call).
    pub fn node_queue(&mut self, q: HighsSearch_NodeQueueRef) -> &mut HighsNodeQueue {
        match q {
            HighsSearch_NodeQueueRef::Global => &mut self.ms.mipdata_.nodequeue,
            HighsSearch_NodeQueueRef::Worker => {
                &mut self.ms.mipdata_.workers[self.worker].nodequeue
            }
        }
    }
}

/// HiGHS: HighsSearch.h::HighsSearch::NodeData (lines 80-112)
#[derive(Clone, Default, Debug)]
pub struct HighsSearch_NodeData {
    /// C++: `double lower_bound`
    pub lower_bound: f64,
    /// C++: `double estimate`
    pub estimate: f64,
    /// C++: `double branching_point`
    pub branching_point: f64,
    /// C++: `double lp_objective`
    pub lp_objective: f64,
    /// C++: `double other_child_lb`
    pub other_child_lb: f64,
    /// C++: `std::shared_ptr<const HighsBasis> nodeBasis`
    pub node_basis: Option<Arc<HighsBasis>>,
    /// C++: `std::shared_ptr<const StabilizerOrbits> stabilizerOrbits`
    pub stabilizer_orbits: Option<Arc<StabilizerOrbits>>,
    /// C++: `HighsDomainChange branchingdecision`
    pub branchingdecision: HighsDomainChange,
    /// C++: `HighsInt domgchgStackPos`
    pub domgchg_stack_pos: HighsInt,
    /// C++: `uint8_t skipDepthCount`
    pub skip_depth_count: u8,
    /// C++: `uint8_t opensubtrees`
    pub opensubtrees: u8,
    /// rusolve mip_symmetry_dynamic: branching path of the node; None when mip_symmetry_dynamic is off.
    pub(crate) s3_path: Option<Arc<S3PathNode>>,
    /// rusolve mip_symmetry_dynamic: path child shared by both children of this node's branching; None when off.
    pub(crate) s3_child_path: Option<Arc<S3PathNode>>,
}

impl HighsSearch_NodeData {
    /// HiGHS: HighsSearch.h::HighsSearch::NodeData::NodeData (lines 98-111)
    /// C++: `NodeData(double parentlb = - kHighsInf, double parentestimate = - kHighsInf, std::shared_ptr<const HighsBasis> parentBasis = nullptr, std::shared_ptr<const StabilizerOrbits> stabilizerOrbits = nullptr)`
    /// C++ default arguments: `parentlb = - kHighsInf`, `parentestimate = - kHighsInf`, `parent_basis = nullptr`, `stabilizer_orbits = nullptr`
    pub fn new(
        parentlb: f64,
        parentestimate: f64,
        parent_basis: Option<Arc<HighsBasis>>,
        stabilizer_orbits: Option<Arc<StabilizerOrbits>>,
    ) -> Self {
        HighsSearch_NodeData {
            lower_bound: parentlb,
            estimate: parentestimate,
            branching_point: 0.0,
            lp_objective: -kHighsInf,
            other_child_lb: parentlb,
            node_basis: parent_basis,
            stabilizer_orbits,
            branchingdecision: HighsDomainChange {
                boundval: 0.0,
                column: -1,
                boundtype: HighsBoundType::kLower,
            },
            domgchg_stack_pos: -1,
            skip_depth_count: 0,
            opensubtrees: 2,
            s3_path: None,
            s3_child_path: None,
        }
    }
}

impl HighsSearch {
    /// HiGHS: HighsSearch.h::HighsSearch::branchingVarReliableAtNodeFlags (lines 125-129)
    /// C++: `int branchingVarReliableAtNodeFlags(HighsInt col) const`
    pub fn branching_var_reliable_at_node_flags(&self, col: HighsInt) -> i32 {
        match self.reliableatnode.find(&col) {
            None => 0,
            Some(&it) => it,
        }
    }

    /// HiGHS: HighsSearch.h::HighsSearch::branchingVarReliableAtNode (lines 131-136)
    /// C++: `bool branchingVarReliableAtNode(HighsInt col) const`
    pub fn branching_var_reliable_at_node(&self, col: HighsInt) -> bool {
        match self.reliableatnode.find(&col) {
            Some(&it) => it == Self::kReliable,
            None => false,
        }
    }

    /// HiGHS: HighsSearch.h::HighsSearch::markBranchingVarUpReliableAtNode (lines 138-140)
    /// C++: `void markBranchingVarUpReliableAtNode(HighsInt col)`
    pub fn mark_branching_var_up_reliable_at_node(&mut self, col: HighsInt) {
        *self.reliableatnode.index_or_insert(&col) |= Self::kUpReliable;
    }

    /// HiGHS: HighsSearch.h::HighsSearch::markBranchingVarDownReliableAtNode (lines 142-144)
    /// C++: `void markBranchingVarDownReliableAtNode(HighsInt col)`
    pub fn mark_branching_var_down_reliable_at_node(&mut self, col: HighsInt) {
        *self.reliableatnode.index_or_insert(&col) |= Self::kDownReliable;
    }

    /// HiGHS: HighsSearch.h::HighsSearch::setLpRelaxation (lines 158-158)
    /// C++: `void setLpRelaxation(HighsLpRelaxation * lp)`
    // PORT NOTE: the LP is not stored: it is `SearchCtx::lp` at every call, so the caller
    // builds the ctx with the LP the C++ would have set (the worker's `lps[lp_]` for the
    // worker's own search, `&heurlp` in the heuristics). Kept as a no-op so call sites map 1:1.
    pub fn set_lp_relaxation(&mut self) {}

    /// HiGHS: HighsSearch.h::HighsSearch::setHeuristic (lines 172-175)
    /// C++: `void setHeuristic(bool inheuristic)`
    pub fn set_heuristic(&mut self, inheuristic: bool) {
        self.inheuristic = inheuristic;
        if inheuristic {
            self.childselrule = HighsSearch_ChildSelectionRule::kHybridInferenceCost;
        }
    }

    /// HiGHS: HighsSearch.h::HighsSearch::hasNode (lines 193-193)
    /// C++: `bool hasNode() const`
    pub fn has_node(&self) -> bool {
        !self.nodestack.is_empty()
    }

    /// HiGHS: HighsSearch.h::HighsSearch::currentNodePruned (lines 195-195)
    /// C++: `bool currentNodePruned() const`
    pub fn current_node_pruned(&self) -> bool {
        self.nodestack.last().unwrap().opensubtrees == 0
    }

    /// HiGHS: HighsSearch.h::HighsSearch::getCurrentEstimate (lines 197-197)
    /// C++: `double getCurrentEstimate() const`
    pub fn get_current_estimate(&self) -> f64 {
        self.nodestack.last().unwrap().estimate
    }

    /// HiGHS: HighsSearch.h::HighsSearch::getCurrentLowerBound (lines 199-199)
    /// C++: `double getCurrentLowerBound() const`
    pub fn get_current_lower_bound(&self) -> f64 {
        self.nodestack.last().unwrap().lower_bound
    }

    /// HiGHS: HighsSearch.h::HighsSearch::getCurrentDepth (lines 201-201)
    /// C++: `HighsInt getCurrentDepth() const`
    pub fn get_current_depth(&self) -> HighsInt {
        self.nodestack.len() as HighsInt + self.depthoffset
    }

    /// HiGHS: HighsSearch.h::HighsSearch::getLocalDomain (lines 240-240)
    /// C++: `HighsDomain & getLocalDomain()`
    pub fn get_local_domain(&mut self) -> &mut HighsDomain {
        &mut self.localdom
    }

    /// HiGHS: HighsSearch.h::HighsSearch::getLocalDomain (lines 242-242)
    /// C++: `const HighsDomain & getLocalDomain() const`
    pub fn get_local_domain_2(&self) -> &HighsDomain {
        &self.localdom
    }

    /// HiGHS: HighsSearch.h::HighsSearch::getPseudoCost (lines 244-244)
    /// C++: `HighsPseudocost & getPseudoCost()`
    // PORT NOTE: the search's pseudocost is `SearchCtx::pseudocost`; associated fn.
    pub fn get_pseudo_cost<'c>(ctx: &'c mut SearchCtx<'_>) -> &'c mut HighsPseudocost {
        ctx.pseudocost
    }

    /// HiGHS: HighsSearch.h::HighsSearch::getPseudoCost (lines 246-246)
    /// C++: `const HighsPseudocost & getPseudoCost() const`
    pub fn get_pseudo_cost_2<'c>(ctx: &'c SearchCtx<'_>) -> &'c HighsPseudocost {
        ctx.pseudocost
    }

    /// Parallel glue (no C++ counterpart): run `f` on worker `w`'s own search
    /// (C++ `*mipdata_->workers[w].search_ptr_`) with the C++ objects it points at:
    /// `search_ptr_` is taken out of the (still attached) worker, the worker's LP
    /// `lps[lp_]` and pseudocost `pseudocosts[pseudocost_]` are detached into the
    /// [`SearchCtx`], and everything is put back afterwards. Callers
    /// (`HighsMipSolver::run`, `evaluateRootNode`) wrap each C++ `search.x(..)` call.
    pub fn with_worker_search<R>(
        ms: &mut HighsMipSolver,
        w: usize,
        f: impl FnOnce(&mut HighsSearch, &mut SearchCtx<'_>) -> R,
    ) -> R {
        let mut search = ms.mipdata_.workers[w]
            .search_ptr_
            .take()
            .expect("HighsSearch::with_worker_search: worker search is running (re-entrant call)");
        let lp_i = ms.mipdata_.workers[w].lp_;
        let pc_i = ms.mipdata_.workers[w].pseudocost_;
        let mut lp = ms.mipdata_.lps[lp_i].take();
        let mut pc = ms.mipdata_.pseudocosts[pc_i].take();
        let result = {
            let mut ctx = SearchCtx::new(ms, w, &mut lp, &mut pc);
            f(&mut search, &mut ctx)
        };
        ms.mipdata_.pseudocosts[pc_i].put(pc);
        ms.mipdata_.lps[lp_i].put(lp);
        ms.mipdata_.workers[w].search_ptr_ = Some(search);
        result
    }

    /// Glue: C++ `getSymmetries().propagateOrbitopes(localdom)`.
    /// `ms.mipdata_.symmetries` is a plain field and the callee takes `ms`, so it is
    /// reached through [`HighsMipSolver::with_symmetries_ref`].
    pub fn propagate_orbitopes(&mut self, ms: &mut HighsMipSolver) -> HighsInt {
        let localdom = &mut self.localdom;
        ms.with_symmetries_ref(|ms, sym| sym.propagate_orbitopes(ms, localdom))
    }

    /// Glue: C++ `orbits->orbitalFixing(localdom)` for a node's stabilizer orbits.
    pub fn orbital_fixing(
        &mut self,
        ms: &mut HighsMipSolver,
        orbits: &StabilizerOrbits,
    ) -> HighsInt {
        let localdom = &mut self.localdom;
        ms.with_symmetries_ref(|ms, sym| orbits.orbital_fixing(ms, sym, localdom))
    }

    /// Glue: C++ `localdom` destructor (`removePropagationDomain` on every pool).
    /// Call where the C++ `HighsSearch` goes out of scope (a heuristic's local
    /// `heur`, `resetSearch`).
    pub fn release(&mut self, ms: &mut HighsMipSolver) {
        self.localdom.release(ms);
    }
}

impl HighsMipSolver {
    /// Glue: C++ `mipdata_->symmetries` for a callee that also takes `ms`
    /// (`propagateOrbitopes`, `orbitalFixing`, `computeStabilizerOrbits`). The field
    /// is moved out for the call and put back; no callee reaches
    /// `mipdata_->symmetries` through `ms` (the domain never reads it).
    pub fn with_symmetries<R>(
        &mut self,
        f: impl FnOnce(&mut HighsMipSolver, &mut HighsSymmetries) -> R,
    ) -> R {
        let mut sym = std::mem::take(&mut self.mipdata_.symmetries);
        let result = f(self, &mut sym);
        self.mipdata_.symmetries = sym;
        result
    }

    /// Glue: read-only C++ `mipdata_->symmetries` (`const` methods
    /// `propagateOrbitopes`, `orbitalFixing`, `computeStabilizerOrbits`) for a callee
    /// that also takes `ms`. Unlike [`Self::with_symmetries`] nothing is moved: the
    /// node-search callers run on parallel workers, and moving the shared field out
    /// would race with the other workers' reads (`getSymmetries()`).
    pub fn with_symmetries_ref<R>(
        &mut self,
        f: impl FnOnce(&mut HighsMipSolver, &HighsSymmetries) -> R,
    ) -> R {
        let sym: *const HighsSymmetries = &self.mipdata_.symmetries;
        // SAFETY: C++ workers read `mipdata_->symmetries` in place through `const`
        // methods; it is written only in serial setup (`runSetup`), never while a
        // search or the callee runs, and no callee writes it through `ms`. The pointer
        // targets a field of `*self`, which outlives `f`. Same SharedMut aliasing
        // discipline as `crate::parallel::shared_mut`.
        f(self, unsafe { &*sym })
    }
}

#[cfg(test)]
mod s3_node_search_tests {
    use super::*;

    #[test]
    fn s3_node_none_default() {
        let node = HighsSearch_NodeData::new(-kHighsInf, -kHighsInf, None, None);
        assert!(node.s3_path.is_none());
        assert!(node.s3_child_path.is_none());

        let stats = S3Stats::default();
        assert_eq!(stats.orbitope_nred, 0);
        assert_eq!(stats.orbitope_ncutoff, 0);
        assert_eq!(stats.lexred_nred, 0);
        assert_eq!(stats.lexred_ncutoff, 0);
    }
}
