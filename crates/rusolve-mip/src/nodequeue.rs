//! Open-node storage and node selection.
//!
//! Two layers, as in HiGHS:
//!
//! * [`NodeQueue`] — the global queue, indexed twice: by lower bound and by a
//!   hybrid key `½·bound + ½·estimate` (HiGHS `NodeLowerRbTree` /
//!   `NodeHybridEstimRbTree`, reference/HiGHS/highs/mip/HighsNodeQueue.cpp:
//!   44–91; `popBestNode` :388 vs `popBestBoundNode` :396). [`NodeQueue::pop`]
//!   takes the best-estimate node, except that after 10 consecutive picks that
//!   were not also best-bound it takes the best-bound node (serial branch of
//!   `installNodes`, HighsMipSolver.cpp:614–628) — estimate-driven search for
//!   good incumbents with a guaranteed bound-improving pick.
//! * [`NodeSelector`] — the plunge: after branching the search dives into the
//!   preferred child; the other child waits on a stack. When a dive ends the
//!   search backtracks to the deepest waiting sibling unless switching
//!   branches higher up looks better by pseudocost score (then the sibling is
//!   parked in the queue), and after 100 nodes the whole stack is flushed to
//!   the queue (HiGHS `backtrackPlunge`, HighsSearch.cpp:1623–1796; plunge
//!   limit 100 from `processNodes`, HighsMipSolver.cpp:910–916).
//!
//! MIT attribution: HiGHS (Copyright (c) HiGHS developers), read as a
//! reference text and reimplemented; red-black trees become `BTreeSet`s, the
//! domain-change-stack length tie-break becomes node depth.

use crate::pseudocost::Pseudocost;
use std::cmp::{Ordering, Reverse};
use std::collections::BTreeSet;

/// Weight of the lower bound in the hybrid key (HighsNodeQueue.cpp:84).
const LB_WEIGHT: f64 = 0.5;
/// Weight of the best-estimate in the hybrid key (HighsNodeQueue.cpp:85).
const ESTIMATE_WEIGHT: f64 = 0.5;
/// Best-estimate picks in a row before one best-bound pick is forced.
const BOUND_PICK_PERIOD: u64 = 10;
/// Nodes processed per queue pick before the plunge stack is flushed.
pub const PLUNGE_LIMIT: u64 = 100;

/// The bound change that created a node: `col ≥ ⌈point⌉` (`up`) or
/// `col ≤ ⌊point⌋`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BranchSide {
    pub col: usize,
    pub point: f64,
    pub up: bool,
}

/// What the selector needs to know about a stored node.
pub trait OpenNode {
    fn lower_bound(&self) -> f64;
    fn estimate(&self) -> f64;
    fn depth(&self) -> u32;
    /// `None` for the root.
    fn side(&self) -> Option<BranchSide>;
}

/// `f64` with the IEEE total order, for use in ordered keys.
#[derive(Debug, Clone, Copy)]
struct Total(f64);

impl PartialEq for Total {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Total {}
impl PartialOrd for Total {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Total {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

/// Lower bound first; shallower, then lower estimate on ties.
type LowerKey = (Total, u32, Total, usize);
/// Hybrid value first; deeper on ties.
type HybridKey = (Total, Reverse<u32>, usize);

pub struct NodeQueue<T> {
    slots: Vec<Option<T>>,
    free: Vec<usize>,
    by_lower: BTreeSet<LowerKey>,
    by_hybrid: BTreeSet<HybridKey>,
    picks: u64,
    last_bound_pick: u64,
}

impl<T: OpenNode> Default for NodeQueue<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: OpenNode> NodeQueue<T> {
    pub fn new() -> Self {
        Self {
            slots: Vec::new(),
            free: Vec::new(),
            by_lower: BTreeSet::new(),
            by_hybrid: BTreeSet::new(),
            picks: 0,
            last_bound_pick: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.by_lower.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_lower.is_empty()
    }

    fn keys(n: &T, id: usize) -> (LowerKey, HybridKey) {
        let (lb, est, d) = (n.lower_bound(), n.estimate(), n.depth());
        (
            (Total(lb), d, Total(est), id),
            (
                Total(LB_WEIGHT * lb + ESTIMATE_WEIGHT * est),
                Reverse(d),
                id,
            ),
        )
    }

    pub fn push(&mut self, n: T) {
        let id = match self.free.pop() {
            Some(id) => id,
            None => {
                self.slots.push(None);
                self.slots.len() - 1
            }
        };
        let (lk, hk) = Self::keys(&n, id);
        self.by_lower.insert(lk);
        self.by_hybrid.insert(hk);
        self.slots[id] = Some(n);
    }

    fn take(&mut self, id: usize) -> T {
        let n = self.slots[id].take().expect("queued slot is occupied");
        let (lk, hk) = Self::keys(&n, id);
        self.by_lower.remove(&lk);
        self.by_hybrid.remove(&hk);
        self.free.push(id);
        n
    }

    /// Smallest lower bound in the queue (`+∞` when empty).
    pub fn min_lower_bound(&self) -> f64 {
        self.by_lower.first().map_or(f64::INFINITY, |k| k.0.0)
    }

    pub fn pop_best_bound(&mut self) -> Option<T> {
        let id = self.by_lower.first()?.3;
        Some(self.take(id))
    }

    pub fn pop_best_estimate(&mut self) -> Option<T> {
        let id = self.by_hybrid.first()?.2;
        Some(self.take(id))
    }

    /// Best-estimate pick, with a best-bound pick forced once
    /// [`BOUND_PICK_PERIOD`] picks in a row were not best-bound.
    pub fn pop(&mut self) -> Option<T> {
        let &(lb, depth, _, _) = self.by_lower.first()?;
        let n = if self.picks - self.last_bound_pick >= BOUND_PICK_PERIOD {
            self.last_bound_pick = self.picks;
            self.pop_best_bound()
        } else {
            let n = self.pop_best_estimate()?;
            if n.lower_bound() == lb.0 && n.depth() == depth {
                self.last_bound_pick = self.picks;
            }
            Some(n)
        };
        self.picks += 1;
        n
    }

    /// Drops every node whose lower bound exceeds `limit`.
    pub fn prune(&mut self, limit: f64) {
        while let Some(&(lb, _, _, id)) = self.by_lower.last() {
            if lb.0 <= limit {
                break;
            }
            self.take(id);
        }
    }
}

/// How a new plunge picks its start node from the queue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeSelection {
    /// HiGHS: hybrid bound/estimate picks with a periodic best-bound pick,
    /// plunge backtracking to siblings ([`NodeQueue::pop`]).
    Hybrid,
    /// Best bound, and no backtracking: when a dive ends its waiting
    /// siblings go to the queue.
    BestBound,
}

/// Plunge stack in front of a [`NodeQueue`].
pub struct NodeSelector<T> {
    selection: NodeSelection,
    queue: NodeQueue<T>,
    /// Siblings still to explore in the current plunge, shallowest first.
    stack: Vec<T>,
    /// Node count when the current plunge started.
    plunge_start: u64,
}

/// Where [`NodeSelector::next`] took a node from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    /// Backtracked to a sibling within the current plunge.
    Plunge,
    /// A new plunge starts from a queued node.
    Queue,
}

impl<T: OpenNode> NodeSelector<T> {
    pub fn new(selection: NodeSelection) -> Self {
        Self {
            selection,
            queue: NodeQueue::new(),
            stack: Vec::new(),
            plunge_start: 0,
        }
    }

    /// Parks the child the plunge did not take.
    pub fn push_sibling(&mut self, n: T) {
        self.stack.push(n);
    }

    pub fn num_open(&self) -> usize {
        self.queue.len() + self.stack.len()
    }

    /// Smallest lower bound among waiting nodes (`+∞` when none).
    pub fn open_bound(&self) -> f64 {
        self.stack
            .iter()
            .map(OpenNode::lower_bound)
            .fold(self.queue.min_lower_bound(), f64::min)
    }

    /// Drops waiting nodes whose bound exceeds `limit`.
    pub fn prune(&mut self, limit: f64) {
        self.queue.prune(limit);
        self.stack.retain(|n| n.lower_bound() <= limit);
    }

    /// Next node once the current dive has ended; `nodes` is the number of
    /// nodes processed so far. `None` when the tree is exhausted.
    pub fn next(
        &mut self,
        ps: &Pseudocost,
        nodes: u64,
        upper_limit: f64,
        feastol: f64,
    ) -> Option<(T, Source)> {
        if self.selection == NodeSelection::BestBound {
            for n in self.stack.drain(..) {
                self.queue.push(n);
            }
            self.queue.prune(upper_limit);
            return self.queue.pop_best_bound().map(|n| (n, Source::Queue));
        }
        if nodes.saturating_sub(self.plunge_start) >= PLUNGE_LIMIT {
            for n in self.stack.drain(..) {
                self.queue.push(n);
            }
        }
        while let Some(n) = self.stack.pop() {
            if n.lower_bound() > upper_limit {
                continue;
            }
            if self.better_to_switch_higher(&n, ps, feastol) {
                self.queue.push(n);
                continue;
            }
            return Some((n, Source::Plunge));
        }
        self.queue.prune(upper_limit);
        let n = self.queue.pop()?;
        self.plunge_start = nodes;
        Some((n, Source::Queue))
    }

    /// HighsSearch.cpp:1738–1770: compare the score of the branch `n` would
    /// open with the score gained by flipping the nearest ancestor that still
    /// has an open branch (the next stack entry).
    fn better_to_switch_higher(&self, n: &T, ps: &Pseudocost, feastol: f64) -> bool {
        let (Some(s), Some(a)) = (n.side(), self.stack.last().and_then(OpenNode::side)) else {
            return false;
        };
        let score = |side: BranchSide, up: bool| {
            if up {
                ps.score_up(side.col, side.point)
            } else {
                ps.score_down(side.col, side.point)
            }
        };
        let node_score = score(s, s.up);
        let inactive = score(a, a.up);
        let active = score(a, !a.up);
        inactive - active > node_score + feastol
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Clone, PartialEq)]
    struct N {
        lb: f64,
        est: f64,
        depth: u32,
        side: Option<BranchSide>,
    }

    impl OpenNode for N {
        fn lower_bound(&self) -> f64 {
            self.lb
        }
        fn estimate(&self) -> f64 {
            self.est
        }
        fn depth(&self) -> u32 {
            self.depth
        }
        fn side(&self) -> Option<BranchSide> {
            self.side
        }
    }

    fn n(lb: f64, est: f64, depth: u32) -> N {
        N {
            lb,
            est,
            depth,
            side: None,
        }
    }

    #[test]
    fn hybrid_pick_with_periodic_best_bound() {
        let mut q = NodeQueue::new();
        // Node 0 has the best bound but a poor estimate; the others rank
        // better by ½lb + ½est.
        q.push(n(0.0, 100.0, 1));
        for k in 0..20 {
            q.push(n(10.0 + k as f64, 10.0 + k as f64, 2));
        }
        assert_eq!(q.min_lower_bound(), 0.0);
        let picks: Vec<f64> = (0..11).map(|_| q.pop().unwrap().lb).collect();
        // Ten estimate picks, then the forced best-bound pick.
        assert_eq!(
            &picks[..10],
            &[10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0]
        );
        assert_eq!(picks[10], 0.0);
        q.prune(25.0);
        assert_eq!(q.len(), 6);
        assert_eq!(q.pop_best_bound().unwrap().lb, 20.0);
    }

    #[test]
    fn estimate_pick_that_is_also_best_bound_resets_the_period() {
        let mut q = NodeQueue::new();
        for k in 0..30 {
            q.push(n(k as f64, k as f64, 1));
        }
        for k in 0..30 {
            assert_eq!(q.pop().unwrap().lb, k as f64);
        }
        assert!(q.pop().is_none());
    }

    #[test]
    fn plunge_backtracks_to_deepest_sibling_then_queue() {
        let ps = Pseudocost::new(4);
        let mut s: NodeSelector<N> = NodeSelector::new(NodeSelection::Hybrid);
        s.queue.push(n(5.0, 5.0, 1));
        s.push_sibling(n(1.0, 1.0, 1));
        s.push_sibling(n(2.0, 2.0, 2));
        let (a, src) = s.next(&ps, 3, f64::INFINITY, 1e-6).unwrap();
        assert_eq!((a.depth, src), (2, Source::Plunge));
        // Bound above the limit: pruned on the way.
        let (b, src) = s
            .next(&ps, 4, 0.5, 1e-6)
            .map_or((None, None), |(n, s)| (Some(n), Some(s)));
        assert!(b.is_none() && src.is_none());
        assert_eq!(s.num_open(), 0);
    }

    #[test]
    fn best_bound_selection_never_backtracks() {
        let ps = Pseudocost::new(4);
        let mut s: NodeSelector<N> = NodeSelector::new(NodeSelection::BestBound);
        s.queue.push(n(1.0, 9.0, 1));
        s.push_sibling(n(2.0, 0.0, 3));
        let (a, src) = s.next(&ps, 1, f64::INFINITY, 1e-6).unwrap();
        assert_eq!((a.lower_bound(), src), (1.0, Source::Queue));
        let (b, _) = s.next(&ps, 2, f64::INFINITY, 1e-6).unwrap();
        assert_eq!(b.lower_bound(), 2.0);
        assert!(s.next(&ps, 3, f64::INFINITY, 1e-6).is_none());
    }

    #[test]
    fn plunge_limit_flushes_stack_to_queue() {
        let ps = Pseudocost::new(1);
        let mut s: NodeSelector<N> = NodeSelector::new(NodeSelection::Hybrid);
        s.push_sibling(n(3.0, 3.0, 1));
        s.push_sibling(n(4.0, 4.0, 2));
        let (_, src) = s.next(&ps, PLUNGE_LIMIT, f64::INFINITY, 1e-6).unwrap();
        assert_eq!(src, Source::Queue);
        assert_eq!(s.num_open(), 1);
        assert_eq!(s.open_bound(), 4.0);
    }

    #[test]
    fn sibling_parked_when_flipping_an_ancestor_scores_higher() {
        let mut ps = Pseudocost::new(2);
        // Column 0: up branches cost a lot, down branches nothing.
        for _ in 0..4 {
            ps.add_observation(0, 1.0, 10.0);
            ps.add_observation(0, -1.0, 0.0);
            ps.add_observation(1, 1.0, 0.0);
            ps.add_observation(1, -1.0, 0.0);
        }
        let side = |col, up| {
            Some(BranchSide {
                col,
                point: 0.5,
                up,
            })
        };
        let mut s: NodeSelector<N> = NodeSelector::new(NodeSelection::Hybrid);
        // Ancestor branched on col 0 and took the down side; its up sibling
        // waits. The deeper sibling (col 1) has a near-zero score.
        s.push_sibling(N {
            side: side(0, true),
            ..n(0.0, 0.0, 1)
        });
        s.push_sibling(N {
            side: side(1, true),
            ..n(0.0, 0.0, 2)
        });
        let (got, src) = s.next(&ps, 2, f64::INFINITY, 1e-6).unwrap();
        // The col-1 sibling went to the queue; plunge continues at col 0.
        assert_eq!((got.side.unwrap().col, src), (0, Source::Plunge));
        assert_eq!(s.queue.len(), 1);
    }
}
