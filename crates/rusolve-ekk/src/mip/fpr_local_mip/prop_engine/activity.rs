//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! prop_engine.cpp: `IndexedMinHeap` out-of-line methods (lines 101-191),
//! contribution/activity helpers (lines 700-759), and the `PropEngine`
//! activity and domain-PQ methods (lines 761-872).

use super::{ActUndo, HeapEntry, IndexedMinHeap, PqUndo, PropEngine};
use crate::mip::fpr_local_mip::{Deadline, VarState};
use crate::prelude::*;

impl IndexedMinHeap {
    /// mip-heuristics (36f7d9d): prop_engine.cpp::IndexedMinHeap::reserve (lines 101-104)
    pub(crate) fn reserve(&mut self, ncol: HighsInt) {
        self.pos.clear();
        self.pos.resize(ncol as usize, Self::kNotPresent);
        self.heap.reserve(ncol as usize);
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::IndexedMinHeap::clear (lines 106-111)
    pub(crate) fn clear(&mut self) {
        for entry in &self.heap {
            self.pos[entry.var as usize] = Self::kNotPresent;
        }
        self.heap.clear();
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::IndexedMinHeap::insert (lines 113-118)
    pub(crate) fn insert(&mut self, key: f64, var: HighsInt) {
        let idx = self.heap.len() as HighsInt;
        self.heap.push(HeapEntry { key, var });
        self.pos[var as usize] = idx;
        self.sift_up(idx);
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::IndexedMinHeap::erase (lines 120-139)
    pub(crate) fn erase(&mut self, var: HighsInt) {
        let idx = self.pos[var as usize];
        let last_idx = self.heap.len() as HighsInt - 1;
        self.pos[var as usize] = Self::kNotPresent;
        if idx == last_idx {
            self.heap.pop();
            return;
        }
        // Move last entry into the hole, then re-heapify from there.
        self.heap[idx as usize] = self.heap[last_idx as usize];
        self.pos[self.heap[idx as usize].var as usize] = idx;
        self.heap.pop();
        // The moved entry might be smaller than its new parent or larger than a
        // child — try both directions (at most one will do work).
        if idx > 0
            && Self::entry_less(
                &self.heap[idx as usize],
                &self.heap[((idx - 1) / 2) as usize],
            )
        {
            self.sift_up(idx);
        } else {
            self.sift_down(idx);
        }
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::IndexedMinHeap::update (lines 141-152)
    pub(crate) fn update(&mut self, var: HighsInt, new_key: f64) {
        let idx = self.pos[var as usize];
        let old_key = self.heap[idx as usize].key;
        self.heap[idx as usize].key = new_key;
        if new_key < old_key {
            self.sift_up(idx);
        } else if new_key > old_key {
            self.sift_down(idx);
        }
        // Exact tie on key: position still satisfies heap property (tiebreak by
        // var is unchanged), so no sift needed.
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::IndexedMinHeap::sift_up (lines 154-167)
    fn sift_up(&mut self, idx: HighsInt) {
        let entry = self.heap[idx as usize];
        let mut idx = idx;
        while idx > 0 {
            let parent = (idx - 1) / 2;
            if !Self::entry_less(&entry, &self.heap[parent as usize]) {
                break;
            }
            self.heap[idx as usize] = self.heap[parent as usize];
            self.pos[self.heap[idx as usize].var as usize] = idx;
            idx = parent;
        }
        self.heap[idx as usize] = entry;
        self.pos[entry.var as usize] = idx;
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::IndexedMinHeap::sift_down (lines 169-191)
    fn sift_down(&mut self, idx: HighsInt) {
        let entry = self.heap[idx as usize];
        let size = self.heap.len() as HighsInt;
        let mut idx = idx;
        loop {
            let left = (2 * idx) + 1;
            if left >= size {
                break;
            }
            let right = left + 1;
            let mut smallest = left;
            if right < size
                && Self::entry_less(&self.heap[right as usize], &self.heap[left as usize])
            {
                smallest = right;
            }
            if !Self::entry_less(&self.heap[smallest as usize], &entry) {
                break;
            }
            self.heap[idx as usize] = self.heap[smallest as usize];
            self.pos[self.heap[idx as usize].var as usize] = idx;
            idx = smallest;
        }
        self.heap[idx as usize] = entry;
        self.pos[entry.var as usize] = idx;
    }
}

// --- anonymous namespace (prop_engine.cpp:700-759) --------------------------

/// mip-heuristics (36f7d9d): prop_engine.cpp::Contribution (lines 703-706)
#[derive(Clone, Copy, Debug)]
struct Contribution {
    lo: f64,
    hi: f64,
}

/// mip-heuristics (36f7d9d): prop_engine.cpp::contribution (lines 710-721)
fn contribution(vs: &VarState, a: f64) -> Contribution {
    if a == 0.0 {
        return Contribution { lo: 0.0, hi: 0.0 };
    }
    if vs.fixed {
        return Contribution {
            lo: a * vs.val,
            hi: a * vs.val,
        };
    }
    if a > 0.0 {
        return Contribution {
            lo: a * vs.lb,
            hi: a * vs.ub,
        };
    }
    Contribution {
        lo: a * vs.ub,
        hi: a * vs.lb,
    }
}

/// mip-heuristics (36f7d9d): prop_engine.cpp::add_contribution (lines 723-729)
fn add_contribution(contrib: f64, finite: &mut f64, ninf: &mut HighsInt) {
    if contrib.is_infinite() {
        *ninf += 1;
    } else {
        *finite += contrib;
    }
}

/// mip-heuristics (36f7d9d): prop_engine.cpp::published_min (lines 733-735)
fn published_min(ninf: HighsInt, finite: f64) -> f64 {
    if ninf > 0 { -kHighsInf } else { finite }
}

/// mip-heuristics (36f7d9d): prop_engine.cpp::published_max (lines 737-739)
fn published_max(ninf: HighsInt, finite: f64) -> f64 {
    if ninf > 0 { kHighsInf } else { finite }
}

/// mip-heuristics (36f7d9d): prop_engine.cpp::replace_contribution (lines 744-757)
fn replace_contribution(old_c: f64, new_c: f64, finite: &mut f64, ninf: &mut HighsInt) {
    let old_inf = old_c.is_infinite();
    let new_inf = new_c.is_infinite();
    if !old_inf && !new_inf {
        *finite += new_c - old_c;
        return;
    }
    if old_inf {
        *ninf -= 1;
    } else {
        *finite -= old_c;
    }
    add_contribution(new_c, finite, ninf);
}

// --- PropEngine activity and domain-PQ methods ------------------------------

impl<'a> PropEngine<'a> {
    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::publish_activity (lines 761-764)
    pub(crate) fn publish_activity(&mut self, i: HighsInt) {
        self.min_activity_[i as usize] =
            published_min(self.min_ninf_[i as usize], self.min_finite_[i as usize]);
        self.max_activity_[i as usize] =
            published_max(self.max_ninf_[i as usize], self.max_finite_[i as usize]);
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::init_activities (lines 766-785)
    pub(crate) fn init_activities(&mut self) {
        self.min_activity_.clear();
        self.min_activity_.resize(self.nrow_ as usize, 0.0);
        self.max_activity_.clear();
        self.max_activity_.resize(self.nrow_ as usize, 0.0);
        self.min_finite_.clear();
        self.min_finite_.resize(self.nrow_ as usize, 0.0);
        self.max_finite_.clear();
        self.max_finite_.resize(self.nrow_ as usize, 0.0);
        self.min_ninf_.clear();
        self.min_ninf_.resize(self.nrow_ as usize, 0);
        self.max_ninf_.clear();
        self.max_ninf_.resize(self.nrow_ as usize, 0);
        self.act_undo_.clear();
        self.act_undo_.reserve(4usize * (self.ncol_ as usize)); // grows with nnz touched during DFS

        for i in 0..self.nrow_ {
            for k in self.ar_start_[i as usize]..self.ar_start_[(i + 1) as usize] {
                let cc = contribution(
                    &self.vs_[self.ar_index_[k as usize] as usize],
                    self.ar_value_[k as usize],
                );
                add_contribution(
                    cc.lo,
                    &mut self.min_finite_[i as usize],
                    &mut self.min_ninf_[i as usize],
                );
                add_contribution(
                    cc.hi,
                    &mut self.max_finite_[i as usize],
                    &mut self.max_ninf_[i as usize],
                );
            }
            self.publish_activity(i);
        }
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::update_activities (lines 787-819)
    pub(crate) fn update_activities(&mut self, j: HighsInt, old_vs: &VarState) {
        if self.min_activity_.is_empty() {
            return;
        }

        let kbeg = self.col_start_[j as usize];
        let kend = self.col_start_[(j + 1) as usize];

        for p in kbeg..kend {
            let i = self.col_row_[p as usize];
            let a = self.col_val_[p as usize];
            let old_c = contribution(old_vs, a);
            let new_c = contribution(&self.vs_[j as usize], a);
            // Component comparisons, not a delta: the same infinite contribution
            // before and after compares equal and skips the row.
            if old_c.lo != new_c.lo || old_c.hi != new_c.hi {
                self.act_undo_.push(ActUndo {
                    row: i,
                    old_min_ninf: self.min_ninf_[i as usize],
                    old_max_ninf: self.max_ninf_[i as usize],
                    old_min_finite: self.min_finite_[i as usize],
                    old_max_finite: self.max_finite_[i as usize],
                });
                replace_contribution(
                    old_c.lo,
                    new_c.lo,
                    &mut self.min_finite_[i as usize],
                    &mut self.min_ninf_[i as usize],
                );
                replace_contribution(
                    old_c.hi,
                    new_c.hi,
                    &mut self.max_finite_[i as usize],
                    &mut self.max_ninf_[i as usize],
                );
                self.min_activity_[i as usize] =
                    published_min(self.min_ninf_[i as usize], self.min_finite_[i as usize]);
                self.max_activity_[i as usize] =
                    published_max(self.max_ninf_[i as usize], self.max_finite_[i as usize]);
            }
        }
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::init_domain_pq (lines 821-831)
    pub(crate) fn init_domain_pq(&mut self) {
        self.domain_pq_.clear();
        self.pq_undo_.clear();
        self.pq_undo_.reserve(4 * (self.ncol_ as usize));
        for j in 0..self.ncol_ {
            if self.is_int(j) && !self.vs_[j as usize].fixed {
                let dom = self.vs_[j as usize].ub - self.vs_[j as usize].lb;
                self.domain_pq_.insert(dom, j);
            }
        }
        self.pq_active_ = true;
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::pq_top (lines 833-838)
    pub(crate) fn pq_top(&self) -> HighsInt {
        if self.domain_pq_.empty() {
            return -1;
        }
        self.domain_pq_.top_var()
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::pq_mark (lines 840-842)
    pub(crate) fn pq_mark(&self) -> HighsInt {
        self.pq_undo_.len() as HighsInt
    }

    /// mip-heuristics (36f7d9d): prop_engine.cpp::PropEngine::pq_notify (lines 844-872)
    pub(crate) fn pq_notify(&mut self, j: HighsInt, old_vs: &VarState) {
        if !self.pq_active_ {
            return;
        }
        if !self.is_int(j) {
            return;
        }

        let old_dom = old_vs.ub - old_vs.lb;
        let was_present = !old_vs.fixed;
        let is_present = !self.vs_[j as usize].fixed;
        let new_dom = self.vs_[j as usize].ub - self.vs_[j as usize].lb;

        // Log undo with both old and new state
        self.pq_undo_.push(PqUndo {
            var: j,
            old_dom,
            new_dom,
            was_present,
            is_present,
        });

        // Apply the transition.  `update(var, new_key)` handles the common
        // "still in PQ but key changed" case without the erase+insert pair
        // the std::set version required.
        if was_present && is_present {
            if old_dom != new_dom {
                self.domain_pq_.update(j, new_dom);
            }
        } else if was_present {
            self.domain_pq_.erase(j);
        } else if is_present {
            self.domain_pq_.insert(new_dom, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::fpr_local_mip::{Deadline, VarState};

    #[test]
    fn s5_propa_heap_empty_and_single() {
        let mut heap = IndexedMinHeap::default();
        heap.reserve(4);
        assert!(heap.empty());
        assert_eq!(heap.size(), 0);
        assert!(!heap.contains(0));

        heap.insert(3.0, 1);
        assert!(!heap.empty());
        assert_eq!(heap.size(), 1);
        assert!(heap.contains(1));
        assert_eq!(heap.top_var(), 1);

        heap.erase(1);
        assert!(heap.empty());
        assert!(!heap.contains(1));
    }

    #[test]
    fn s5_propa_heap_ordering_tiebreak() {
        let mut heap = IndexedMinHeap::default();
        heap.reserve(8);
        // Insert in non-sorted order; top_var must reflect min(key, var).
        heap.insert(5.0, 2);
        heap.insert(3.0, 0);
        heap.insert(3.0, 1); // same key as var 0; var 0 wins the tiebreak
        heap.insert(7.0, 3);
        assert_eq!(heap.top_var(), 0);

        heap.erase(0);
        assert_eq!(heap.top_var(), 1); // tiebreak now selects next-smallest var

        heap.erase(1);
        assert_eq!(heap.top_var(), 2);

        heap.erase(2);
        assert_eq!(heap.top_var(), 3);

        heap.erase(3);
        assert!(heap.empty());
    }

    #[test]
    fn s5_propa_heap_update_both_directions() {
        let mut heap = IndexedMinHeap::default();
        heap.reserve(4);
        heap.insert(10.0, 0);
        heap.insert(20.0, 1);
        heap.insert(30.0, 2);
        assert_eq!(heap.top_var(), 0);

        // Decrease-key: var 2 beats var 0 after update.
        heap.update(2, 5.0);
        assert_eq!(heap.top_var(), 2);
        assert!(heap.contains(0));
        assert!(heap.contains(1));

        // Increase-key: var 2 drops past both peers.
        heap.update(2, 100.0);
        assert_eq!(heap.top_var(), 0);

        // No-op update (same key) must be a true no-op.
        heap.update(0, 10.0);
        assert_eq!(heap.top_var(), 0);
    }

    #[test]
    fn s5_propa_heap_erase_root_reheapify() {
        let mut heap = IndexedMinHeap::default();
        heap.reserve(8);
        // Distinct vars with distinct keys (keys chosen to force non-sorted
        // insertion order so sift_up is exercised on the build-up).
        let entries: [(f64, HighsInt); 8] = [
            (7.0, 0),
            (3.0, 1),
            (5.0, 2),
            (1.0, 3),
            (6.0, 4),
            (4.0, 5),
            (8.0, 6),
            (2.0, 7),
        ];
        for &(key, var) in &entries {
            heap.insert(key, var);
        }
        // Pop in ascending-key order: 1.0/3, 2.0/7, 3.0/1, 4.0/5, 5.0/2, 6.0/4,
        // 7.0/0, 8.0/6.  Each pop forces erase-of-root followed by sift_down of
        // the moved-last entry.
        assert_eq!(heap.top_var(), 3);
        heap.erase(3);
        assert_eq!(heap.top_var(), 7);
        heap.erase(7);
        assert_eq!(heap.top_var(), 1);
        heap.erase(1);
        assert_eq!(heap.top_var(), 5);
        heap.erase(5);
        assert_eq!(heap.top_var(), 2);
    }

    #[test]
    fn s5_propa_heap_clear_retains_capacity() {
        let mut heap = IndexedMinHeap::default();
        heap.reserve(4);
        heap.insert(1.0, 0);
        heap.insert(2.0, 1);
        heap.insert(3.0, 2);

        heap.clear();
        assert!(heap.empty());
        assert!(!heap.contains(0));
        assert!(!heap.contains(1));
        assert!(!heap.contains(2));

        // Re-populating after clear must behave like a fresh heap.
        heap.insert(10.0, 2);
        heap.insert(5.0, 0);
        assert_eq!(heap.top_var(), 0);
        assert!(heap.contains(2));
    }

    #[test]
    fn s5_propa_contribution_inf_ninf() {
        // Zero coefficient gives {0, 0} regardless of bounds.
        let vs = VarState {
            lb: -kHighsInf,
            ub: kHighsInf,
            val: 0.0,
            fixed: false,
        };
        let cc = contribution(&vs, 0.0);
        assert_eq!(cc.lo, 0.0);
        assert_eq!(cc.hi, 0.0);

        // replace_contribution with +inf toggling ninf.
        let mut finite: f64 = 5.0;
        let mut ninf: HighsInt = 0;
        replace_contribution(1.0, kHighsInf, &mut finite, &mut ninf);
        assert_eq!(ninf, 1);
        assert_eq!(finite, 4.0); // 5.0 - 1.0 = 4.0

        // published_min/max with ninf > 0.
        assert_eq!(published_min(ninf, finite), -kHighsInf);
        assert_eq!(published_max(ninf, finite), kHighsInf);

        // replace_contribution removing the inf.
        replace_contribution(kHighsInf, 2.0, &mut finite, &mut ninf);
        assert_eq!(ninf, 0);
        assert_eq!(finite, 6.0); // 4.0 + 2.0 = 6.0

        assert_eq!(published_min(ninf, finite), 6.0);
        assert_eq!(published_max(ninf, finite), 6.0);

        // replace_contribution with -inf.
        replace_contribution(2.0, -kHighsInf, &mut finite, &mut ninf);
        assert_eq!(ninf, 1);
        assert_eq!(finite, 4.0); // 6.0 - 2.0 = 4.0

        assert_eq!(published_min(ninf, finite), -kHighsInf);
        assert_eq!(published_max(ninf, finite), kHighsInf);
    }

    #[test]
    fn s5_propa_activity_bound_change() {
        // 2 rows, 2 columns.
        // Row 0: 1*x0 + 2*x1, row 1: 3*x0 + 4*x1
        // CSC: col 0 -> rows 0,1; col 1 -> rows 0,1
        let ar_start: Vec<HighsInt> = vec![0, 2, 4];
        let ar_index: Vec<HighsInt> = vec![0, 1, 0, 1];
        let ar_value: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
        let col_start: Vec<HighsInt> = vec![0, 2, 4];
        let col_row: Vec<HighsInt> = vec![0, 1, 0, 1];
        let col_val: Vec<f64> = vec![1.0, 3.0, 2.0, 4.0];
        let col_lb: Vec<f64> = vec![0.0, 0.0];
        let col_ub: Vec<f64> = vec![10.0, 10.0];
        let row_lo: Vec<f64> = vec![0.0, 0.0];
        let row_hi: Vec<f64> = vec![kHighsInf, kHighsInf];
        let integrality: Vec<HighsVarType> =
            vec![HighsVarType::kContinuous, HighsVarType::kContinuous];

        let vs: Vec<VarState> = vec![
            VarState {
                lb: 0.0,
                ub: 10.0,
                val: 0.0,
                fixed: false,
            },
            VarState {
                lb: 0.0,
                ub: 10.0,
                val: 0.0,
                fixed: false,
            },
        ];

        let mut pe = PropEngine {
            ncol_: 2,
            nrow_: 2,
            ar_start_: &ar_start,
            ar_index_: &ar_index,
            ar_value_: &ar_value,
            col_start_: &col_start,
            col_row_: &col_row,
            col_val_: &col_val,
            col_lb_: &col_lb,
            col_ub_: &col_ub,
            row_lo_: &row_lo,
            row_hi_: &row_hi,
            integrality_: &integrality,
            feastol_: 1e-7,
            nnz_: 4,
            vs_: vs.clone(),
            solution_: vec![0.0, 0.0],
            vs_undo_: Vec::new(),
            sol_undo_: Vec::new(),
            prop_worklist_: Vec::new(),
            prop_in_wl_: Vec::new(),
            prop_work_: 0,
            deadline_: Deadline {
                timer: None,
                limit: 0.0,
            },
            min_activity_: Vec::new(),
            max_activity_: Vec::new(),
            min_finite_: Vec::new(),
            max_finite_: Vec::new(),
            min_ninf_: Vec::new(),
            max_ninf_: Vec::new(),
            act_undo_: Vec::new(),
            pq_active_: false,
            domain_pq_: IndexedMinHeap::default(),
            pq_undo_: Vec::new(),
        };

        pe.init_activities();

        // Row 0: min = 1*0 + 2*0 = 0, max = 1*10 + 2*10 = 30
        assert_eq!(pe.min_activity_[0], 0.0);
        assert_eq!(pe.max_activity_[0], 30.0);
        // Row 1: min = 3*0 + 4*0 = 0, max = 3*10 + 4*10 = 70
        assert_eq!(pe.min_activity_[1], 0.0);
        assert_eq!(pe.max_activity_[1], 70.0);

        // Change x0 upper bound to 5.
        let old_vs = pe.vs_[0];
        pe.vs_[0].ub = 5.0;
        pe.update_activities(0, &old_vs);

        // Row 0: max = 1*5 + 2*10 = 25
        assert_eq!(pe.max_activity_[0], 25.0);
        // Row 1: max = 3*5 + 4*10 = 55
        assert_eq!(pe.max_activity_[1], 55.0);
        // min unchanged
        assert_eq!(pe.min_activity_[0], 0.0);
        assert_eq!(pe.min_activity_[1], 0.0);

        // Now test infinite bound toggling ninf.
        // Reset: set x0 ub back to 10, then set x0 lb to -inf.
        let old_vs2 = pe.vs_[0];
        pe.vs_[0].ub = 10.0;
        pe.vs_[0].lb = -kHighsInf;
        pe.update_activities(0, &old_vs2);

        // Row 0: min contribution of x0 = 1*(-inf) = -inf => ninf=1 => min=-inf
        assert!(pe.min_activity_[0].is_infinite());
        assert!(pe.min_activity_[0] < 0.0);
        // Row 1: min contribution of x0 = 3*(-inf) = -inf => ninf=1 => min=-inf
        assert!(pe.min_activity_[1].is_infinite());
        assert!(pe.min_activity_[1] < 0.0);
        // max unchanged
        assert_eq!(pe.max_activity_[0], 30.0);
        assert_eq!(pe.max_activity_[1], 70.0);

        // Now fix x0 to 3.0 (removes the -inf).
        let old_vs3 = pe.vs_[0];
        pe.vs_[0].lb = 3.0;
        pe.vs_[0].ub = 3.0;
        pe.vs_[0].val = 3.0;
        pe.vs_[0].fixed = true;
        pe.update_activities(0, &old_vs3);

        // Row 0: min = 1*3 + 2*0 = 3, max = 1*3 + 2*10 = 23
        assert_eq!(pe.min_activity_[0], 3.0);
        assert_eq!(pe.max_activity_[0], 23.0);
        // Row 1: min = 3*3 + 4*0 = 9, max = 3*3 + 4*10 = 49
        assert_eq!(pe.min_activity_[1], 9.0);
        assert_eq!(pe.max_activity_[1], 49.0);
    }

    /// Builds an engine by struct literal (see the constructor), with
    /// every column at its structural bounds, unfixed.
    #[allow(clippy::too_many_arguments)]
    fn literal_engine<'a>(
        ncol: HighsInt,
        nrow: HighsInt,
        ar_start: &'a [HighsInt],
        ar_index: &'a [HighsInt],
        ar_value: &'a [f64],
        col_start: &'a [HighsInt],
        col_row: &'a [HighsInt],
        col_val: &'a [f64],
        col_lb: &'a [f64],
        col_ub: &'a [f64],
        row_lo: &'a [f64],
        row_hi: &'a [f64],
        integrality: &'a [HighsVarType],
    ) -> PropEngine<'a> {
        let vs: Vec<VarState> = (0..ncol as usize)
            .map(|j| VarState {
                lb: col_lb[j],
                ub: col_ub[j],
                val: 0.0,
                fixed: false,
            })
            .collect();
        PropEngine {
            ncol_: ncol,
            nrow_: nrow,
            ar_start_: ar_start,
            ar_index_: ar_index,
            ar_value_: ar_value,
            col_start_: col_start,
            col_row_: col_row,
            col_val_: col_val,
            col_lb_: col_lb,
            col_ub_: col_ub,
            row_lo_: row_lo,
            row_hi_: row_hi,
            integrality_: integrality,
            feastol_: 1e-6,
            nnz_: ar_value.len(),
            vs_: vs,
            solution_: vec![0.0; ncol as usize],
            vs_undo_: Vec::new(),
            sol_undo_: Vec::new(),
            prop_worklist_: Vec::new(),
            prop_in_wl_: Vec::new(),
            prop_work_: 0,
            deadline_: Deadline {
                timer: None,
                limit: 0.0,
            },
            min_activity_: Vec::new(),
            max_activity_: Vec::new(),
            min_finite_: Vec::new(),
            max_finite_: Vec::new(),
            min_ninf_: Vec::new(),
            max_ninf_: Vec::new(),
            act_undo_: Vec::new(),
            pq_active_: false,
            domain_pq_: {
                let mut h = IndexedMinHeap::default();
                h.reserve(ncol);
                h
            },
            pq_undo_: Vec::new(),
        }
    }

    /// The activity half of fix() (prop_engine.cpp:245-251): mark fixed at
    /// `value`, then update_activities and pq_notify with the old state.
    fn fix_state(pe: &mut PropEngine<'_>, j: HighsInt, value: f64) {
        let old_vs = pe.vs_[j as usize];
        pe.vs_[j as usize].fixed = true;
        pe.vs_[j as usize].val = value;
        pe.update_activities(j, &old_vs);
        pe.pq_notify(j, &old_vs);
    }

    // Port of test_prop_engine.cpp:594-672 ("fixing a column with an infinite
    // bound keeps row activity finite"), the parts reachable without fix()'s
    // bound checks or backtrack_to.
    #[test]
    fn s5_propa_activity_infinite_bound_fix() {
        // x0 integer [0, +inf), x1 continuous (-inf, +inf), x2 integer [0, 4].
        //   row 0: 2 x0 - x2 <= 10
        //   row 1: x0 + x1 >= 1
        let ar_start: Vec<HighsInt> = vec![0, 2, 4];
        let ar_index: Vec<HighsInt> = vec![0, 2, 0, 1];
        let ar_value: Vec<f64> = vec![2.0, -1.0, 1.0, 1.0];
        let col_start: Vec<HighsInt> = vec![0, 2, 3, 4];
        let col_row: Vec<HighsInt> = vec![0, 1, 1, 0];
        let col_val: Vec<f64> = vec![2.0, 1.0, 1.0, -1.0];
        let col_lb: Vec<f64> = vec![0.0, -kHighsInf, 0.0];
        let col_ub: Vec<f64> = vec![kHighsInf, kHighsInf, 4.0];
        let row_lo: Vec<f64> = vec![-kHighsInf, 1.0];
        let row_hi: Vec<f64> = vec![10.0, kHighsInf];
        let integrality: Vec<HighsVarType> = vec![
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
        ];
        let mut pe = literal_engine(
            3,
            2,
            &ar_start,
            &ar_index,
            &ar_value,
            &col_start,
            &col_row,
            &col_val,
            &col_lb,
            &col_ub,
            &row_lo,
            &row_hi,
            &integrality,
        );
        pe.init_activities();
        assert_eq!(pe.row_min_activity(0), -4.0);
        assert_eq!(pe.row_max_activity(0), kHighsInf);
        assert_eq!(pe.row_min_activity(1), -kHighsInf);
        assert_eq!(pe.row_max_activity(1), kHighsInf);

        fix_state(&mut pe, 0, 1e5);
        assert_eq!(pe.row_min_activity(0), 2e5 - 4.0);
        assert_eq!(pe.row_max_activity(0), 2e5);
        // Row 1 keeps x1's two infinite contributions: still unbounded, not NaN.
        assert_eq!(pe.row_min_activity(1), -kHighsInf);
        assert_eq!(pe.row_max_activity(1), kHighsInf);

        fix_state(&mut pe, 1, -3.0);
        assert_eq!(pe.row_min_activity(1), 1e5 - 3.0);
        assert_eq!(pe.row_max_activity(1), 1e5 - 3.0);
        // One undo record per touched row: x0 touched rows 0 and 1, x1 row 1.
        assert_eq!(pe.act_undo_.len(), 3);
    }

    // Port of test_prop_engine.cpp:674-697 ("a zero coefficient on an
    // unbounded column leaves activity finite").
    #[test]
    fn s5_propa_activity_zero_coefficient() {
        // x0 continuous (-inf, +inf) with a stored 0, x1 integer [0, 3].
        //   row 0: 0 x0 + x1 <= 2
        let ar_start: Vec<HighsInt> = vec![0, 2];
        let ar_index: Vec<HighsInt> = vec![0, 1];
        let ar_value: Vec<f64> = vec![0.0, 1.0];
        let col_start: Vec<HighsInt> = vec![0, 1, 2];
        let col_row: Vec<HighsInt> = vec![0, 0];
        let col_val: Vec<f64> = vec![0.0, 1.0];
        let col_lb: Vec<f64> = vec![-kHighsInf, 0.0];
        let col_ub: Vec<f64> = vec![kHighsInf, 3.0];
        let row_lo: Vec<f64> = vec![-kHighsInf];
        let row_hi: Vec<f64> = vec![2.0];
        let integrality: Vec<HighsVarType> =
            vec![HighsVarType::kContinuous, HighsVarType::kInteger];
        let mut pe = literal_engine(
            2,
            1,
            &ar_start,
            &ar_index,
            &ar_value,
            &col_start,
            &col_row,
            &col_val,
            &col_lb,
            &col_ub,
            &row_lo,
            &row_hi,
            &integrality,
        );
        pe.init_activities();
        assert_eq!(pe.row_min_activity(0), 0.0);
        assert_eq!(pe.row_max_activity(0), 3.0);
        fix_state(&mut pe, 0, 4.0);
        fix_state(&mut pe, 1, 3.0);
        assert_eq!(pe.row_min_activity(0), 3.0);
        assert_eq!(pe.row_max_activity(0), 3.0);
    }

    #[test]
    fn s5_propa_domain_pq_transitions() {
        // Three integer columns with domains 4, 2, 2 and one continuous one,
        // in a single empty-coefficient row (activities are not exercised).
        let ar_start: Vec<HighsInt> = vec![0, 0];
        let ar_index: Vec<HighsInt> = vec![];
        let ar_value: Vec<f64> = vec![];
        let col_start: Vec<HighsInt> = vec![0, 0, 0, 0, 0];
        let col_row: Vec<HighsInt> = vec![];
        let col_val: Vec<f64> = vec![];
        let col_lb: Vec<f64> = vec![0.0, 1.0, 0.0, 0.0];
        let col_ub: Vec<f64> = vec![4.0, 3.0, 2.0, 1.0];
        let row_lo: Vec<f64> = vec![-kHighsInf];
        let row_hi: Vec<f64> = vec![kHighsInf];
        let integrality: Vec<HighsVarType> = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
        ];
        let mut pe = literal_engine(
            4,
            1,
            &ar_start,
            &ar_index,
            &ar_value,
            &col_start,
            &col_row,
            &col_val,
            &col_lb,
            &col_ub,
            &row_lo,
            &row_hi,
            &integrality,
        );
        // Inactive PQ: notify is a no-op and logs nothing.
        fix_state(&mut pe, 0, 0.0);
        assert_eq!(pe.pq_mark(), 0);
        pe.vs_[0].fixed = false;

        pe.init_domain_pq();
        assert!(!pe.domain_pq_.contains(3));
        // Domains 2 and 2 tie; var 1 wins the tiebreak.
        assert_eq!(pe.pq_top(), 1);

        // Narrow x0 to [0, 1]: still present, key 4 -> 1, now the top.
        let old_vs = pe.vs_[0];
        pe.vs_[0].ub = 1.0;
        pe.pq_notify(0, &old_vs);
        assert_eq!(pe.pq_top(), 0);
        assert_eq!(pe.pq_mark(), 1);

        // Fix x0: erased.
        fix_state(&mut pe, 0, 1.0);
        assert!(!pe.domain_pq_.contains(0));
        assert_eq!(pe.pq_top(), 1);
        assert_eq!(pe.pq_mark(), 2);
        let u = pe.pq_undo_[1];
        assert_eq!(u.var, 0);
        assert_eq!(u.old_dom, 1.0);
        assert!(u.was_present && !u.is_present);

        // Unfix x0: reinserted with its current domain.
        let old_vs = pe.vs_[0];
        pe.vs_[0].fixed = false;
        pe.pq_notify(0, &old_vs);
        assert!(pe.domain_pq_.contains(0));
        assert_eq!(pe.pq_top(), 0);

        // A continuous column never logs.
        let old_vs = pe.vs_[3];
        pe.vs_[3].ub = 0.5;
        pe.pq_notify(3, &old_vs);
        assert_eq!(pe.pq_mark(), 3);

        // Fix every integer column: empty PQ reports -1.
        fix_state(&mut pe, 0, 0.0);
        fix_state(&mut pe, 1, 1.0);
        fix_state(&mut pe, 2, 0.0);
        assert_eq!(pe.pq_top(), -1);
    }
}
