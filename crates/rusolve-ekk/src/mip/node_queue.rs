//! HiGHS source: mip/HighsNodeQueue.h, mip/HighsNodeQueue.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::mip::symmetry_dynamic::S3PathNode;
use crate::prelude::*;
use std::sync::Arc;

/// HiGHS: HighsNodeQueue.h::HighsNodeQueue (lines 28-298)
///
/// Removed C++ members:
/// * C++ `std::unique_ptr<AllocatorState> allocatorState`: the node-set pool allocator is dropped (memory only).
#[derive(Clone, Debug)]
pub struct HighsNodeQueue {
    /// C++: `std::vector<OpenNode> nodes`
    pub nodes: Vec<HighsNodeQueue_OpenNode>,
    /// C++: `std::priority_queue<int64_t, std::vector<int64_t>, std::greater<int64_t>> freeslots`
    pub freeslots: BinaryHeap<Reverse<i64>>,
    /// C++: `NodeSetArray colLowerNodesPtr`
    pub col_lower_nodes_ptr: Vec<NodeSet>,
    /// C++: `NodeSetArray colUpperNodesPtr`
    pub col_upper_nodes_ptr: Vec<NodeSet>,
    /// C++: `int64_t lowerRoot = - 1`
    pub lower_root: i64,
    /// C++: `int64_t lowerMin = - 1`
    pub lower_min: i64,
    /// C++: `int64_t hybridEstimRoot = - 1`
    pub hybrid_estim_root: i64,
    /// C++: `int64_t hybridEstimMin = - 1`
    pub hybrid_estim_min: i64,
    /// C++: `int64_t suboptimalRoot = - 1`
    pub suboptimal_root: i64,
    /// C++: `int64_t suboptimalMin = - 1`
    pub suboptimal_min: i64,
    /// C++: `int64_t numSuboptimal = 0`
    pub num_suboptimal: i64,
    /// C++: `double optimality_limit = kHighsInf`
    pub optimality_limit: f64,
    /// C++: `HighsInt numCol = 0`
    pub num_col: HighsInt,
}

impl Default for HighsNodeQueue {
    fn default() -> Self {
        Self {
            nodes: Default::default(),
            freeslots: Default::default(),
            col_lower_nodes_ptr: Default::default(),
            col_upper_nodes_ptr: Default::default(),
            lower_root: -1,
            lower_min: -1,
            hybrid_estim_root: -1,
            hybrid_estim_min: -1,
            suboptimal_root: -1,
            suboptimal_min: -1,
            num_suboptimal: 0,
            optimality_limit: kHighsInf,
            num_col: 0,
        }
    }
}

/// HiGHS: HighsNodeQueue.h::HighsNodeQueue::OpenNode (lines 142-178)
#[derive(Clone, Default, Debug)]
pub struct HighsNodeQueue_OpenNode {
    /// C++: `std::vector<HighsDomainChange> domchgstack`
    pub domchgstack: Vec<HighsDomainChange>,
    /// C++: `std::vector<HighsInt> branchings`
    pub branchings: Vec<HighsInt>,
    /// C++: `std::vector<NodeSet::iterator> domchglinks`
    pub domchglinks: Vec<NodeSetKey>,
    /// C++: `double lower_bound`
    pub lower_bound: f64,
    /// C++: `double estimate`
    pub estimate: f64,
    /// C++: `HighsInt depth`
    pub depth: HighsInt,
    /// C++: `highs::RbTreeLinks<int64_t> lowerLinks`
    pub lower_links: RbTreeLinks<i64>,
    /// C++: `highs::RbTreeLinks<int64_t> hybridEstimLinks`
    pub hybrid_estim_links: RbTreeLinks<i64>,
    /// rusolve mip_symmetry_dynamic: branching path of this open node; None when mip_symmetry_dynamic is off.
    pub(crate) s3_path: Option<Arc<S3PathNode>>,
    /// rusolve mip_node_install_basis: the basis the node's parent stored (shared with its
    /// siblings), installed when the node is taken from the queue; None when the flag is off.
    pub(crate) basis: Option<Arc<HighsBasis>>,
}

impl HighsNodeQueue {
    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::link_estim (lines 115-119)
    /// C++: `void link_estim(int64_t node)`
    pub fn link_estim(&mut self, node: i64) {
        debug_assert!(node != -1);
        let mut rb_tree = HighsNodeQueue_NodeHybridEstimRbTree::new(self);
        rb_tree.cached_link(node);
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::unlink_estim (lines 121-125)
    /// C++: `void unlink_estim(int64_t node)`
    pub fn unlink_estim(&mut self, node: i64) {
        debug_assert!(node != -1);
        let mut rb_tree = HighsNodeQueue_NodeHybridEstimRbTree::new(self);
        rb_tree.cached_unlink(node);
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::link_lower (lines 127-131)
    /// C++: `void link_lower(int64_t node)`
    pub fn link_lower(&mut self, node: i64) {
        debug_assert!(node != -1);
        let mut rb_tree = HighsNodeQueue_NodeLowerRbTree::new(self);
        rb_tree.cached_link(node);
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::unlink_lower (lines 133-137)
    /// C++: `void unlink_lower(int64_t node)`
    pub fn unlink_lower(&mut self, node: i64) {
        debug_assert!(node != -1);
        let mut rb_tree = HighsNodeQueue_NodeLowerRbTree::new(self);
        rb_tree.cached_unlink(node);
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::link_suboptimal (lines 139-144)
    /// C++: `void link_suboptimal(int64_t node)`
    pub fn link_suboptimal(&mut self, node: i64) {
        debug_assert!(node != -1);
        let mut rb_tree = HighsNodeQueue_SuboptimalNodeRbTree::new(self);
        rb_tree.cached_link(node);
        self.num_suboptimal += 1;
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::unlink_suboptimal (lines 146-151)
    /// C++: `void unlink_suboptimal(int64_t node)`
    pub fn unlink_suboptimal(&mut self, node: i64) {
        debug_assert!(node != -1);
        let mut rb_tree = HighsNodeQueue_SuboptimalNodeRbTree::new(self);
        rb_tree.cached_unlink(node);
        self.num_suboptimal -= 1;
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::link_domchgs (lines 153-171)
    /// C++: `void link_domchgs(int64_t node)`
    pub fn link_domchgs(&mut self, node: i64) {
        debug_assert!(node != -1);
        let n = node as usize;
        let numchgs = self.nodes[n].domchgstack.len() as HighsInt;
        self.nodes[n]
            .domchglinks
            .resize(numchgs as usize, NodeSetKey::default());

        for i in 0..numchgs as usize {
            let val = self.nodes[n].domchgstack[i].boundval;
            let col = self.nodes[n].domchgstack[i].column as usize;
            // PORT NOTE: `emplace(val, node).first` becomes an insert of the key, which is the
            // iterator's identity here (an equal existing key is the same element).
            let key: NodeSetKey = (OrderedF64(val), node);
            match self.nodes[n].domchgstack[i].boundtype {
                HighsBoundType::kLower => {
                    self.col_lower_nodes_ptr[col].insert(key);
                }
                HighsBoundType::kUpper => {
                    self.col_upper_nodes_ptr[col].insert(key);
                }
            }
            self.nodes[n].domchglinks[i] = key;
        }
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::unlink_domchgs (lines 173-190)
    /// C++: `void unlink_domchgs(int64_t node)`
    pub fn unlink_domchgs(&mut self, node: i64) {
        debug_assert!(node != -1);
        let n = node as usize;
        let numchgs = self.nodes[n].domchgstack.len() as HighsInt;

        for i in 0..numchgs as usize {
            let col = self.nodes[n].domchgstack[i].column as usize;
            // PORT NOTE: `erase(iterator)` becomes `remove(&key)` on the stored key.
            let key = self.nodes[n].domchglinks[i];
            match self.nodes[n].domchgstack[i].boundtype {
                HighsBoundType::kLower => {
                    self.col_lower_nodes_ptr[col].remove(&key);
                }
                HighsBoundType::kUpper => {
                    self.col_upper_nodes_ptr[col].remove(&key);
                }
            }
        }

        self.nodes[n].domchglinks.clear();
        self.nodes[n].domchglinks.shrink_to_fit();
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::link (lines 192-205)
    /// C++: `double link(int64_t node)`
    pub fn link(&mut self, node: i64) -> f64 {
        if self.nodes[node as usize].lower_bound > self.optimality_limit {
            debug_assert!(self.nodes[node as usize].estimate != kHighsInf);
            self.nodes[node as usize].estimate = kHighsInf;
            self.link_suboptimal(node);
            self.link_domchgs(node);
            return ldexp_f64(1.0, 1 - self.nodes[node as usize].depth);
        }

        self.link_estim(node);
        self.link_lower(node);
        self.link_domchgs(node);
        0.0
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::unlink (lines 207-216)
    /// C++: `void unlink(int64_t node)`
    pub fn unlink(&mut self, node: i64) {
        if self.nodes[node as usize].estimate == kHighsInf {
            self.unlink_suboptimal(node);
        } else {
            self.unlink_estim(node);
            self.unlink_lower(node);
        }
        self.unlink_domchgs(node);
        self.freeslots.push(Reverse(node));
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::setNumCol (lines 218-234)
    /// C++: `void setNumCol(HighsInt numCol)`
    pub fn set_num_col(&mut self, num_col: HighsInt) {
        if self.num_col == num_col {
            return;
        }
        self.num_col = num_col;
        // PORT NOTE: C++ resets `allocatorState` here; the node-set pool allocator is dropped, so
        // there is nothing to reset.
        if num_col == 0 {
            return;
        }
        self.col_lower_nodes_ptr = (0..num_col as usize).map(|_| NodeSet::new()).collect();
        self.col_upper_nodes_ptr = (0..num_col as usize).map(|_| NodeSet::new()).collect();
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::checkGlobalBounds (lines 236-259)
    /// C++: `void checkGlobalBounds(HighsInt col, double lb, double ub, double feastol, HighsCDouble & treeweight)`
    pub fn check_global_bounds(
        &mut self,
        col: HighsInt,
        lb: f64,
        ub: f64,
        feastol: f64,
        treeweight: &mut HighsCDouble,
    ) {
        // C++ `std::set<int64_t> delnodes`
        let mut delnodes: BTreeSet<i64> = BTreeSet::new();

        // C++ `colLowerNodes[col].lower_bound(std::make_pair(ub + feastol, -1))`: first element >= key
        let prune_start: NodeSetKey = (OrderedF64(ub + feastol), -1);
        for &(_, node) in self.col_lower_nodes_ptr[col as usize].range(prune_start..) {
            delnodes.insert(node);
        }

        // C++ `colUpperNodes[col].upper_bound(std::make_pair(lb - feastol, kHighsIInf))`: elements <= key
        let prune_end: NodeSetKey = (OrderedF64(lb - feastol), kHighsIInf as i64);
        for &(_, node) in self.col_upper_nodes_ptr[col as usize].range(..=prune_end) {
            delnodes.insert(node);
        }

        for delnode in delnodes {
            if self.nodes[delnode as usize].estimate != kHighsInf {
                *treeweight += ldexp_f64(1.0, 1 - self.nodes[delnode as usize].depth);
            }
            self.unlink(delnode);
        }
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::pruneInfeasibleNodes (lines 261-309)
    /// C++: `double pruneInfeasibleNodes(HighsDomain & globaldomain, double feastol)`
    // PORT NOTE: signature changed because propagate needs ms
    pub fn prune_infeasible_nodes(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldomain: &mut HighsDomain,
        feastol: f64,
    ) -> f64 {
        let mut treeweight = HighsCDouble::from(0.0);

        loop {
            if globaldomain.infeasible() {
                break;
            }

            let numchgs = globaldomain.get_domain_change_stack().len();

            debug_assert!(self.num_col == globaldomain.col_lower_.len() as HighsInt);

            for i in 0..self.num_col {
                let lb = globaldomain.col_lower_[i as usize];
                let ub = globaldomain.col_upper_[i as usize];
                self.check_global_bounds(i, lb, ub, feastol, &mut treeweight);
            }

            let numopennodes = self.num_nodes() as usize;
            if numopennodes == 0 {
                break;
            }

            for i in 0..self.num_col as usize {
                if self.col_lower_nodes_ptr[i].len() == numopennodes {
                    let globallb = self.col_lower_nodes_ptr[i].iter().next().unwrap().0.0;
                    if globallb > globaldomain.col_lower_[i] {
                        globaldomain.change_bound_2(
                            ms,
                            HighsBoundType::kLower,
                            i as HighsInt,
                            globallb,
                            HighsDomain_Reason::unspecified(),
                        );
                        if globaldomain.infeasible() {
                            break;
                        }
                    }
                }

                if self.col_upper_nodes_ptr[i].len() == numopennodes {
                    let globalub = self.col_upper_nodes_ptr[i].iter().next_back().unwrap().0.0;
                    if globalub < globaldomain.col_upper_[i] {
                        globaldomain.change_bound_2(
                            ms,
                            HighsBoundType::kUpper,
                            i as HighsInt,
                            globalub,
                            HighsDomain_Reason::unspecified(),
                        );
                        if globaldomain.infeasible() {
                            break;
                        }
                    }
                }
            }

            globaldomain.propagate(ms);

            if numchgs == globaldomain.get_domain_change_stack().len() {
                break;
            }
        }

        f64::from(treeweight)
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::pruneNode (lines 311-317)
    /// C++: `double pruneNode(int64_t nodeId)`
    pub fn prune_node(&mut self, node_id: i64) -> f64 {
        let treeweight = if self.nodes[node_id as usize].estimate != kHighsInf {
            ldexp_f64(1.0, 1 - self.nodes[node_id as usize].depth)
        } else {
            0.0
        };
        self.unlink(node_id);
        treeweight
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::performBounding (lines 319-360)
    /// C++: `double performBounding(double upper_limit)`
    pub fn perform_bounding(&mut self, upper_limit: f64) -> f64 {
        // PORT NOTE: C++ keeps one `lowerTree` view alive across `pruneNode` / `unlink_*`. Those
        // calls mutate the same root/min fields, but `predecessor` reads only the child and parent
        // links, so a fresh view per `predecessor` call gives the same result. `last()` and
        // `empty()` are called before any mutation, as in C++.
        if HighsNodeQueue_NodeLowerRbTree::new(self).empty() {
            return 0.0;
        }

        let mut treeweight = HighsCDouble::from(0.0);

        let mut max_lb_node = HighsNodeQueue_NodeLowerRbTree::new(self).last();
        while max_lb_node != -1 {
            if self.nodes[max_lb_node as usize].lower_bound < upper_limit {
                break;
            }
            let next = HighsNodeQueue_NodeLowerRbTree::new(self).predecessor(max_lb_node);
            treeweight += self.prune_node(max_lb_node);
            max_lb_node = next;
        }

        if self.optimality_limit < upper_limit {
            while max_lb_node != -1 {
                if self.nodes[max_lb_node as usize].lower_bound < self.optimality_limit {
                    break;
                }
                let next = HighsNodeQueue_NodeLowerRbTree::new(self).predecessor(max_lb_node);
                debug_assert!(self.nodes[max_lb_node as usize].estimate != kHighsInf);
                self.unlink_estim(max_lb_node);
                self.unlink_lower(max_lb_node);
                treeweight += ldexp_f64(1.0, 1 - self.nodes[max_lb_node as usize].depth);
                self.nodes[max_lb_node as usize].estimate = kHighsInf;
                self.link_suboptimal(max_lb_node);
                max_lb_node = next;
            }
        }

        if self.num_suboptimal != 0 {
            max_lb_node = HighsNodeQueue_SuboptimalNodeRbTree::new(self).last();
            while max_lb_node != -1 {
                if self.nodes[max_lb_node as usize].lower_bound < upper_limit {
                    break;
                }
                let next = HighsNodeQueue_SuboptimalNodeRbTree::new(self).predecessor(max_lb_node);
                self.unlink(max_lb_node);
                max_lb_node = next;
            }
        }

        f64::from(treeweight)
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::emplaceNode (lines 362-386)
    /// C++: `double emplaceNode(std::vector<HighsDomainChange> && domchgs, std::vector<HighsInt> && branchPositions, double lower_bound, double estimate, HighsInt depth)`
    pub fn emplace_node(
        &mut self,
        domchgs: Vec<HighsDomainChange>,
        branch_positions: Vec<HighsInt>,
        lower_bound: f64,
        estimate: f64,
        depth: HighsInt,
    ) -> f64 {
        let pos: i64;

        debug_assert!(estimate != kHighsInf);

        if self.freeslots.is_empty() {
            pos = self.nodes.len() as i64;
            self.nodes.push(HighsNodeQueue_OpenNode::new_2(
                domchgs,
                branch_positions,
                lower_bound,
                estimate,
                depth,
            ));
        } else {
            pos = self.freeslots.peek().expect("freeslots is non-empty").0;
            self.freeslots.pop();
            self.nodes[pos as usize] = HighsNodeQueue_OpenNode::new_2(
                domchgs,
                branch_positions,
                lower_bound,
                estimate,
                depth,
            );
        }

        debug_assert!(self.nodes[pos as usize].lower_bound == lower_bound);
        debug_assert!(self.nodes[pos as usize].estimate == estimate);
        debug_assert!(self.nodes[pos as usize].depth == depth);

        self.link(pos)
    }

    /// rusolve mip_symmetry_dynamic: copy of emplace_node that also stores the branching path; emplace_node is unchanged.
    pub(crate) fn emplace_node_s3(
        &mut self,
        domchgs: Vec<HighsDomainChange>,
        branch_positions: Vec<HighsInt>,
        lower_bound: f64,
        estimate: f64,
        depth: HighsInt,
        s3_path: Option<Arc<S3PathNode>>,
    ) -> f64 {
        self.emplace_node_ext(
            domchgs,
            branch_positions,
            lower_bound,
            estimate,
            depth,
            s3_path,
            None,
        )
    }

    /// rusolve mip_symmetry_dynamic / mip_node_install_basis: copy of emplace_node that also
    /// stores the branching path and the parent basis; emplace_node is unchanged.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn emplace_node_ext(
        &mut self,
        domchgs: Vec<HighsDomainChange>,
        branch_positions: Vec<HighsInt>,
        lower_bound: f64,
        estimate: f64,
        depth: HighsInt,
        s3_path: Option<Arc<S3PathNode>>,
        basis: Option<Arc<HighsBasis>>,
    ) -> f64 {
        let pos: i64;

        debug_assert!(estimate != kHighsInf);

        let mut node =
            HighsNodeQueue_OpenNode::new_2(domchgs, branch_positions, lower_bound, estimate, depth);
        node.s3_path = s3_path;
        node.basis = basis;
        if self.freeslots.is_empty() {
            pos = self.nodes.len() as i64;
            self.nodes.push(node);
        } else {
            pos = self.freeslots.peek().expect("freeslots is non-empty").0;
            self.freeslots.pop();
            self.nodes[pos as usize] = node;
        }

        debug_assert!(self.nodes[pos as usize].lower_bound == lower_bound);
        debug_assert!(self.nodes[pos as usize].estimate == estimate);
        debug_assert!(self.nodes[pos as usize].depth == depth);

        self.link(pos)
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::popBestNode (lines 388-394)
    /// C++: `OpenNode && popBestNode()`
    pub fn pop_best_node(&mut self) -> HighsNodeQueue_OpenNode {
        let best_node = self.hybrid_estim_min;

        self.unlink(best_node);

        // PORT NOTE: C++ returns `std::move(nodes[bestNode])`, a reference into the vector. The
        // slot is moved-from and only read again by `emplace_node`, which assigns it. Returning
        // the value and leaving `Default` in the slot is equivalent.
        std::mem::take(&mut self.nodes[best_node as usize])
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::popBestBoundNode (lines 396-402)
    /// C++: `OpenNode && popBestBoundNode()`
    pub fn pop_best_bound_node(&mut self) -> HighsNodeQueue_OpenNode {
        let best_bound_node = self.lower_min;

        self.unlink(best_bound_node);

        // PORT NOTE: as in `pop_best_node`, the moved-from slot is left at `Default`.
        std::mem::take(&mut self.nodes[best_bound_node as usize])
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::getBestLowerBound (lines 404-410)
    /// C++: `double getBestLowerBound() const`
    pub fn get_best_lower_bound(&self) -> f64 {
        let lb = if self.lower_min == -1 {
            kHighsInf
        } else {
            self.nodes[self.lower_min as usize].lower_bound
        };

        if self.suboptimal_min == -1 {
            return lb;
        }

        cpp_min(self.nodes[self.suboptimal_min as usize].lower_bound, lb)
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::getBestBoundDomchgStackSize (lines 412-420)
    /// C++: `HighsInt getBestBoundDomchgStackSize() const`
    pub fn get_best_bound_domchg_stack_size(&self) -> HighsInt {
        let domchg_stack_size: HighsInt = if self.lower_min == -1 {
            kHighsIInf
        } else {
            self.nodes[self.lower_min as usize].domchgstack.len() as HighsInt
        };

        if self.suboptimal_min == -1 {
            return domchg_stack_size;
        }

        cpp_min(
            self.nodes[self.suboptimal_min as usize].domchgstack.len() as HighsInt,
            domchg_stack_size,
        )
    }

    /// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::clear (lines 422-444)
    /// C++: `void clear()`
    pub fn clear(&mut self) {
        // PORT NOTE: C++ `clear()` takes its `original == false` branch: it copies every member
        // from a fresh queue with the same `numCol`. That resets `optimality_limit` to kHighsInf.
        self.nodes = Vec::new();
        self.col_lower_nodes_ptr = (0..self.num_col as usize).map(|_| NodeSet::new()).collect();
        self.col_upper_nodes_ptr = (0..self.num_col as usize).map(|_| NodeSet::new()).collect();
        self.freeslots = BinaryHeap::new();
        self.lower_root = -1;
        self.lower_min = -1;
        self.hybrid_estim_root = -1;
        self.hybrid_estim_min = -1;
        self.suboptimal_root = -1;
        self.suboptimal_min = -1;
        self.num_suboptimal = 0;
        self.optimality_limit = kHighsInf;
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::setOptimalityLimit (lines 233-235)
    /// C++: `void setOptimalityLimit(double optimality_limit)`
    pub fn set_optimality_limit(&mut self, optimality_limit: f64) {
        self.optimality_limit = optimality_limit;
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::numNodesUp (lines 249-251)
    /// C++: `int64_t numNodesUp(HighsInt col) const`
    pub fn num_nodes_up(&self, col: HighsInt) -> i64 {
        self.col_lower_nodes_ptr[col as usize].len() as i64
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::numNodesDown (lines 253-255)
    /// C++: `int64_t numNodesDown(HighsInt col) const`
    pub fn num_nodes_down(&self, col: HighsInt) -> i64 {
        self.col_upper_nodes_ptr[col as usize].len() as i64
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::numNodesUp (lines 257-263)
    /// C++: `int64_t numNodesUp(HighsInt col, double val) const`
    pub fn num_nodes_up_2(&self, col: HighsInt, val: f64) -> i64 {
        debug_assert!(self.num_col > col);
        let set = &self.col_lower_nodes_ptr[col as usize];
        // C++ `it = upper_bound(std::make_pair(val, kHighsIInf))`: first element > key
        let key: NodeSetKey = (OrderedF64(val), kHighsIInf as i64);
        // `it == begin` means no element is <= key
        if set.range(..=key).next().is_none() {
            return set.len() as i64;
        }
        // std::distance(it, end): elements > key
        set.range((std::ops::Bound::Excluded(key), std::ops::Bound::Unbounded))
            .count() as i64
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::numNodesDown (lines 265-271)
    /// C++: `int64_t numNodesDown(HighsInt col, double val) const`
    pub fn num_nodes_down_2(&self, col: HighsInt, val: f64) -> i64 {
        debug_assert!(self.num_col > col);
        let set = &self.col_upper_nodes_ptr[col as usize];
        // C++ `it = lower_bound(std::make_pair(val, -1))`: first element >= key
        let key: NodeSetKey = (OrderedF64(val), -1);
        // `it == end` means no element is >= key
        if set.range(key..).next().is_none() {
            return set.len() as i64;
        }
        // std::distance(begin, it): elements < key
        set.range(..key).count() as i64
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::getUpNodes (lines 273-275)
    /// C++: `const NodeSet & getUpNodes(HighsInt col) const`
    pub fn get_up_nodes(&self, col: HighsInt) -> &NodeSet {
        &self.col_lower_nodes_ptr[col as usize]
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::getDownNodes (lines 277-279)
    /// C++: `const NodeSet & getDownNodes(HighsInt col) const`
    pub fn get_down_nodes(&self, col: HighsInt) -> &NodeSet {
        &self.col_upper_nodes_ptr[col as usize]
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::numNodes (lines 291-291)
    /// C++: `int64_t numNodes() const`
    pub fn num_nodes(&self) -> i64 {
        self.nodes.len() as i64 - self.freeslots.len() as i64
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::numActiveNodes (lines 293-295)
    /// C++: `int64_t numActiveNodes() const`
    pub fn num_active_nodes(&self) -> i64 {
        self.nodes.len() as i64 - self.freeslots.len() as i64 - self.num_suboptimal
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::empty (lines 297-297)
    /// C++: `bool empty() const`
    pub fn empty(&self) -> bool {
        self.num_active_nodes() == 0
    }
}

impl HighsNodeQueue_OpenNode {
    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::OpenNode::OpenNode (lines 152-160)
    /// C++: `OpenNode()`
    pub fn new() -> Self {
        Self {
            domchgstack: Vec::new(),
            branchings: Vec::new(),
            domchglinks: Vec::new(),
            lower_bound: -kHighsInf,
            estimate: -kHighsInf,
            depth: 0,
            lower_links: RbTreeLinks::default(),
            hybrid_estim_links: RbTreeLinks::default(),
            s3_path: None,
            basis: None,
        }
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::OpenNode::OpenNode (lines 162-171)
    /// C++: `OpenNode(std::vector<HighsDomainChange> && domchgstack, std::vector<HighsInt> && branchings, double lower_bound, double estimate, HighsInt depth)`
    pub fn new_2(
        domchgstack: Vec<HighsDomainChange>,
        branchings: Vec<HighsInt>,
        lower_bound: f64,
        estimate: f64,
        depth: HighsInt,
    ) -> Self {
        Self {
            domchgstack,
            branchings,
            domchglinks: Vec::new(),
            lower_bound,
            estimate,
            depth,
            lower_links: RbTreeLinks::default(),
            hybrid_estim_links: RbTreeLinks::default(),
            s3_path: None,
            basis: None,
        }
    }

    /// HiGHS: HighsNodeQueue.h::HighsNodeQueue::OpenNode::operator= (lines 173-173)
    /// C++: `OpenNode & operator=(OpenNode && other)`
    pub fn op_assign(&mut self, other: HighsNodeQueue_OpenNode) -> &mut HighsNodeQueue_OpenNode {
        *self = other;
        self
    }

    /// rusolve mip_symmetry_dynamic: clone of the branching path of this node (None when the flag is off).
    pub(crate) fn s3_path_handle(&self) -> Option<Arc<S3PathNode>> {
        self.s3_path.clone()
    }
}

/// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::NodeLowerRbTree (lines 44-65)
///
/// CRTP tree object as a view: C++ holds `HighsNodeQueue* nodeQueue` and
/// the base holds `int64_t& rootNode = lowerRoot`, `int64_t& first_ =
/// lowerMin`. Construct where C++ does (`NodeLowerRbTree rbTree(this);`):
/// `HighsNodeQueue_NodeLowerRbTree::new(self)`. Key:
/// `std::tuple<double, HighsInt, double, int64_t>`. (Accessors translated
/// in the skeleton; the tree algorithms are `util::highs_rb_tree`.)
pub struct HighsNodeQueue_NodeLowerRbTree<'a> {
    pub nodes: &'a mut Vec<HighsNodeQueue_OpenNode>,
    pub root: &'a mut i64,
    pub first: &'a mut i64,
}

impl<'a> HighsNodeQueue_NodeLowerRbTree<'a> {
    /// HiGHS: HighsNodeQueue.cpp::NodeLowerRbTree::NodeLowerRbTree (lines 48-51)
    pub fn new(node_queue: &'a mut HighsNodeQueue) -> Self {
        Self {
            nodes: &mut node_queue.nodes,
            root: &mut node_queue.lower_root,
            first: &mut node_queue.lower_min,
        }
    }
}

impl RbTree for HighsNodeQueue_NodeLowerRbTree<'_> {
    type KeyType = (f64, HighsInt, f64, i64);
    type LinkType = i64;
    fn root_node(&self) -> i64 {
        *self.root
    }
    fn set_root_node(&mut self, root: i64) {
        *self.root = root;
    }
    /// HiGHS: HighsNodeQueue.cpp::NodeLowerRbTree::getRbTreeLinks (lines 56-58)
    fn get_rb_tree_links(&self, node: i64) -> &RbTreeLinks<i64> {
        &self.nodes[node as usize].lower_links
    }
    /// HiGHS: HighsNodeQueue.cpp::NodeLowerRbTree::getRbTreeLinks (lines 53-55)
    fn get_rb_tree_links_mut(&mut self, node: i64) -> &mut RbTreeLinks<i64> {
        &mut self.nodes[node as usize].lower_links
    }
    /// HiGHS: HighsNodeQueue.cpp::NodeLowerRbTree::getKey (lines 59-64)
    fn get_key(&self, node: i64) -> (f64, HighsInt, f64, i64) {
        let n = &self.nodes[node as usize];
        (
            n.lower_bound,
            n.domchgstack.len() as HighsInt,
            n.estimate,
            node,
        )
    }
}

impl CacheMinRbTree for HighsNodeQueue_NodeLowerRbTree<'_> {
    fn first_cached_ref(&self) -> i64 {
        *self.first
    }
    fn set_first_cached(&mut self, first: i64) {
        *self.first = first;
    }
}

/// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::NodeHybridEstimRbTree (lines 67-91)
/// (view over `nodes`, `hybridEstimRoot`, `hybridEstimMin`; key
/// `std::tuple<double, HighsInt, int64_t>`).
pub struct HighsNodeQueue_NodeHybridEstimRbTree<'a> {
    pub nodes: &'a mut Vec<HighsNodeQueue_OpenNode>,
    pub root: &'a mut i64,
    pub first: &'a mut i64,
}

impl<'a> HighsNodeQueue_NodeHybridEstimRbTree<'a> {
    /// HiGHS: HighsNodeQueue.cpp::NodeHybridEstimRbTree::NodeHybridEstimRbTree (lines 72-75)
    pub fn new(node_queue: &'a mut HighsNodeQueue) -> Self {
        Self {
            nodes: &mut node_queue.nodes,
            root: &mut node_queue.hybrid_estim_root,
            first: &mut node_queue.hybrid_estim_min,
        }
    }
}

impl RbTree for HighsNodeQueue_NodeHybridEstimRbTree<'_> {
    type KeyType = (f64, HighsInt, i64);
    type LinkType = i64;
    fn root_node(&self) -> i64 {
        *self.root
    }
    fn set_root_node(&mut self, root: i64) {
        *self.root = root;
    }
    fn get_rb_tree_links(&self, node: i64) -> &RbTreeLinks<i64> {
        &self.nodes[node as usize].hybrid_estim_links
    }
    fn get_rb_tree_links_mut(&mut self, node: i64) -> &mut RbTreeLinks<i64> {
        &mut self.nodes[node as usize].hybrid_estim_links
    }
    /// HiGHS: HighsNodeQueue.cpp::NodeHybridEstimRbTree::getKey (lines 83-90)
    fn get_key(&self, node: i64) -> (f64, HighsInt, i64) {
        const K_LB_WEIGHT: f64 = 0.5;
        const K_ESTIM_WEIGHT: f64 = 0.5;
        let n = &self.nodes[node as usize];
        (
            K_LB_WEIGHT * n.lower_bound + K_ESTIM_WEIGHT * n.estimate,
            -(n.domchgstack.len() as HighsInt),
            node,
        )
    }
}

impl CacheMinRbTree for HighsNodeQueue_NodeHybridEstimRbTree<'_> {
    fn first_cached_ref(&self) -> i64 {
        *self.first
    }
    fn set_first_cached(&mut self, first: i64) {
        *self.first = first;
    }
}

/// HiGHS: HighsNodeQueue.cpp::HighsNodeQueue::SuboptimalNodeRbTree (lines 93-114)
/// (view over `nodes` (`lowerLinks`!), `suboptimalRoot`, `suboptimalMin`;
/// key `std::pair<double, int64_t>`).
pub struct HighsNodeQueue_SuboptimalNodeRbTree<'a> {
    pub nodes: &'a mut Vec<HighsNodeQueue_OpenNode>,
    pub root: &'a mut i64,
    pub first: &'a mut i64,
}

impl<'a> HighsNodeQueue_SuboptimalNodeRbTree<'a> {
    /// HiGHS: HighsNodeQueue.cpp::SuboptimalNodeRbTree::SuboptimalNodeRbTree (lines 98-101)
    pub fn new(node_queue: &'a mut HighsNodeQueue) -> Self {
        Self {
            nodes: &mut node_queue.nodes,
            root: &mut node_queue.suboptimal_root,
            first: &mut node_queue.suboptimal_min,
        }
    }
}

impl RbTree for HighsNodeQueue_SuboptimalNodeRbTree<'_> {
    type KeyType = (f64, i64);
    type LinkType = i64;
    fn root_node(&self) -> i64 {
        *self.root
    }
    fn set_root_node(&mut self, root: i64) {
        *self.root = root;
    }
    fn get_rb_tree_links(&self, node: i64) -> &RbTreeLinks<i64> {
        &self.nodes[node as usize].lower_links
    }
    fn get_rb_tree_links_mut(&mut self, node: i64) -> &mut RbTreeLinks<i64> {
        &mut self.nodes[node as usize].lower_links
    }
    /// HiGHS: HighsNodeQueue.cpp::SuboptimalNodeRbTree::getKey (lines 110-112)
    fn get_key(&self, node: i64) -> (f64, i64) {
        (self.nodes[node as usize].lower_bound, node)
    }
}

impl CacheMinRbTree for HighsNodeQueue_SuboptimalNodeRbTree<'_> {
    fn first_cached_ref(&self) -> i64 {
        *self.first
    }
    fn set_first_cached(&mut self, first: i64) {
        *self.first = first;
    }
}

#[cfg(test)]
mod s3_node_queue_tests {
    use super::*;

    #[test]
    fn s3_node_queue_path_round_trip() {
        let mut q = HighsNodeQueue::default();
        let path = Arc::new(S3PathNode {
            parent: None,
            branch_col: 3,
            depth: 1,
            orbi: Vec::new(),
        });

        q.emplace_node_s3(Vec::new(), Vec::new(), 1.0, 2.0, 1, Some(path.clone()));
        let popped = q.pop_best_node();
        assert!(popped.s3_path_handle().is_some());
        assert!(Arc::ptr_eq(&popped.s3_path_handle().unwrap(), &path));

        q.emplace_node(Vec::new(), Vec::new(), 1.0, 2.0, 1);
        let popped2 = q.pop_best_node();
        assert!(popped2.s3_path_handle().is_none());
    }
}
