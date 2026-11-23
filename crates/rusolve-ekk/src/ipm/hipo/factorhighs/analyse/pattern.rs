//! HiGHS source: ipm/hipo/factorhighs/Analyse.cpp (lines 710-946)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Supernodal pattern, relative indices, storage and critical path.

use super::super::super::{Int, Int64};
use super::super::*;
use super::Analyse;
use crate::ipm::hipo::auxiliary::{childrenLinkedList, counts2Ptr, getDiagStart, snFlops, snSpops};
use crate::util::highs_utils::{cpp_max, cpp_min};

impl<'a> Analyse<'a> {
    /// HiGHS: Analyse.cpp::Analyse::snPattern (lines 710-755)
    pub fn snPattern(&mut self) {
        // number of total indices needed
        let mut indices: Int64 = 0;

        for &i in &self.sn_indices_ {
            indices += i as Int64;
        }

        // allocate space for sn pattern
        self.rows_sn_.resize(indices as usize, 0);
        self.ptr_sn_.resize(self.sn_count_ as usize + 1, 0);

        // keep track of visited supernodes
        let mut mark = vec![-1 as Int; self.sn_count_ as usize];

        // compute column pointers of L
        let mut work = vec![0 as Int64; self.sn_indices_.len()];
        for i in 0..self.sn_indices_.len() as Int {
            work[i as usize] = self.sn_indices_[i as usize] as Int64;
        }
        counts2Ptr(&mut self.ptr_sn_, &mut work);

        // consider each row
        for i in 0..self.n_ {
            // for all entries in the row of lower triangle
            for el in self.ptr_upper_[i as usize]..self.ptr_upper_[(i + 1) as usize] {
                // there is nonzero (i,j)
                let j = self.rows_upper_[el as usize];

                // supernode to which column j belongs to
                let mut snj = self.sn_belong_[j as usize];

                // while supernodes are not yet considered
                while snj != -1 && mark[snj as usize] != i {
                    // we may end up too far
                    if self.sn_start_[snj as usize] > i {
                        break;
                    }

                    // supernode snj is now considered for row i
                    mark[snj as usize] = i;

                    // there is a nonzero entry in supernode snj at row i
                    let idx = work[snj as usize] as usize;
                    self.rows_sn_[idx] = i;
                    work[snj as usize] += 1;

                    // go up the elimination tree
                    snj = self.sn_parent_[snj as usize];
                }
            }
        }
    }

    /// HiGHS: Analyse.cpp::Analyse::relativeIndCols (lines 757-801)
    pub fn relativeIndCols(&mut self) {
        // Find the relative indices of the original column wrt the frontal matrix of
        // the corresponding supernode

        self.relind_cols_.resize(self.nz_ as usize, 0);

        // go through the supernodes
        for sn in 0..self.sn_count_ {
            let ptL_start = self.ptr_sn_[sn as usize];
            let ptL_end = self.ptr_sn_[(sn + 1) as usize];

            // go through the columns of the supernode
            for col in self.sn_start_[sn as usize]..self.sn_start_[(sn + 1) as usize] {
                // go through original column and supernodal column
                let mut ptA = self.ptr_lower_[col as usize];
                let mut ptL = ptL_start;

                // offset wrt ptrLower[col]
                let mut index: Int = 0;

                // size of the column of the original matrix
                let col_size = self.ptr_lower_[(col + 1) as usize] - self.ptr_lower_[col as usize];

                while ptL < ptL_end {
                    // if found all the relative indices that are needed, stop
                    if index == col_size {
                        break;
                    }

                    // check if indices coincide
                    if self.rows_sn_[ptL as usize] == self.rows_lower_[ptA as usize] {
                        // yes: save relative index and move pointers forward
                        // NB: difference fits into Int
                        self.relind_cols_[(self.ptr_lower_[col as usize] + index) as usize] =
                            (ptL - ptL_start) as Int;
                        index += 1;
                        ptL += 1;
                        ptA += 1;
                    } else {
                        // no: move pointer of L forward
                        ptL += 1;
                    }
                }
            }
        }
    }

    /// HiGHS: Analyse.cpp::Analyse::relativeIndClique (lines 803-885)
    pub fn relativeIndClique(&mut self) {
        // Find the relative indices of the child clique wrt the frontal matrix of the
        // parent supernode

        self.relind_clique_
            .resize(self.sn_count_ as usize, Vec::new());
        self.consecutive_sums_
            .resize(self.sn_count_ as usize, Vec::new());

        for sn in 0..self.sn_count_ {
            // if there is no parent, skip supernode
            if self.sn_parent_[sn as usize] == -1 {
                continue;
            }

            // number of nodes in the supernode
            let sn_size = self.sn_start_[(sn + 1) as usize] - self.sn_start_[sn as usize];

            // column of the first node in the supernode
            let _j = self.sn_start_[sn as usize];

            // size of the first column of the supernode
            // C++ narrows the Int64 difference to Int
            let sn_column_size =
                (self.ptr_sn_[(sn + 1) as usize] - self.ptr_sn_[sn as usize]) as Int;

            // size of the clique of the supernode
            let sn_clique_size: Int = sn_column_size - sn_size;

            // count number of assembly operations during factorise
            self.sparse_ops_ += snSpops(sn_clique_size as f64);

            self.relind_clique_[sn as usize].resize(sn_clique_size as usize, 0);

            // iterate through the clique of sn
            let mut ptr_current = self.ptr_sn_[sn as usize] + sn_size as Int64;

            // iterate through the full column of parent sn
            let mut ptr_parent = self.ptr_sn_[self.sn_parent_[sn as usize] as usize];

            // keep track of start and end of parent sn column
            let ptr_parent_start = ptr_parent;
            let ptr_parent_end = self.ptr_sn_[(self.sn_parent_[sn as usize] + 1) as usize];

            // where to write into relind
            let mut index: Int = 0;

            // iterate though the column of the parent sn
            while ptr_parent < ptr_parent_end {
                // if found all the relative indices that are needed, stop
                if index == sn_clique_size {
                    break;
                }

                // check if indices coincide
                if self.rows_sn_[ptr_current as usize] == self.rows_sn_[ptr_parent as usize] {
                    // yes: save relative index and move pointers forward
                    // NB: difference fits into Int
                    self.relind_clique_[sn as usize][index as usize] =
                        (ptr_parent - ptr_parent_start) as Int;
                    index += 1;
                    ptr_parent += 1;
                    ptr_current += 1;
                } else {
                    // no: move pointer of parent forward
                    ptr_parent += 1;
                }
            }

            // Difference between consecutive relative indices.
            // Useful to detect chains of consecutive indices.
            self.consecutive_sums_[sn as usize].resize(sn_clique_size as usize, 0);
            for i in 0..sn_clique_size - 1 {
                self.consecutive_sums_[sn as usize][i as usize] = self.relind_clique_[sn as usize]
                    [(i + 1) as usize]
                    - self.relind_clique_[sn as usize][i as usize];
            }

            // Number of consecutive sums that can be done in one blas call.
            self.consecutive_sums_[sn as usize][sn_clique_size as usize - 1] = 1;
            let mut i = sn_clique_size - 2;
            while i >= 0 {
                if self.consecutive_sums_[sn as usize][i as usize] > 1 {
                    self.consecutive_sums_[sn as usize][i as usize] = 1;
                } else if self.consecutive_sums_[sn as usize][i as usize] == 1 {
                    self.consecutive_sums_[sn as usize][i as usize] =
                        self.consecutive_sums_[sn as usize][(i + 1) as usize] + 1;
                } else {
                    if let Some(logger) = self.logger_ {
                        logger.printInfo("Error in consecutiveSums\n");
                    }
                }
                i -= 1;
            }
        }
    }

    /// HiGHS: Analyse.cpp::Analyse::computeStorage (lines 887-905)
    pub fn computeStorage(&self, fr: Int, sz: Int, fr_entries: &mut Int64, cl_entries: &mut Int64) {
        // compute storage required by frontal and clique, based on the format used

        let cl = fr - sz;

        let mut n_blocks = (sz - 1) / self.nb_ + 1;
        let mut temp: Vec<Int64> = Vec::new();
        *fr_entries = getDiagStart(fr, sz, self.nb_, n_blocks, &mut temp, false);

        // clique is stored as a collection of rectangles
        n_blocks = (cl - 1) / self.nb_ + 1;
        let mut schur_size: Int64 = 0;
        for j in 0..n_blocks {
            let jb = cpp_min(self.nb_, cl - j * self.nb_);
            schur_size += (cl - j * self.nb_) as Int64 * jb as Int64;
        }
        *cl_entries = schur_size;
    }

    /// HiGHS: Analyse.cpp::Analyse::computeCriticalPath (lines 907-946)
    pub fn computeCriticalPath(&mut self) {
        // Compute the critical path within the supernodal elimination tree, and the
        // number of operations along the path. This is the number of operations that
        // need to be done sequentially while doing tree parallelism.

        let mut critical_ops = vec![0.0f64; self.sn_count_ as usize];

        // linked lists of children
        let mut head: Vec<Int> = Vec::new();
        let mut next: Vec<Int> = Vec::new();
        childrenLinkedList(&self.sn_parent_, &mut head, &mut next);

        for sn in 0..self.sn_count_ {
            // supernode size
            let sz = self.sn_start_[(sn + 1) as usize] - self.sn_start_[sn as usize];

            // frontal size (C++ narrows the Int64 difference to Int)
            let fr = (self.ptr_sn_[(sn + 1) as usize] - self.ptr_sn_[sn as usize]) as Int;

            critical_ops[sn as usize] = snFlops(sz as f64, (fr - sz) as f64);
        }

        for sn in 0..self.sn_count_ {
            // leaf nodes
            if head[sn as usize] == -1 {
                continue;
            }

            let mut max_ops: f64 = 0.0;
            let mut child = head[sn as usize];
            while child != -1 {
                // critical_ops of this supernode is max over children of
                // (ops_of_this_sn + critical_ops_of_child)
                max_ops = cpp_max(
                    max_ops,
                    critical_ops[sn as usize] + critical_ops[child as usize],
                );
                child = next[child as usize];
            }
            critical_ops[sn as usize] = max_ops;
        }

        for sn in 0..self.sn_count_ {
            self.critical_ops_ = cpp_max(self.critical_ops_, critical_ops[sn as usize]);
        }
    }
}
