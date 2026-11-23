//! HiGHS source: ipm/hipo/factorhighs/Analyse.h, ipm/hipo/factorhighs/Analyse.cpp (lines 1-367)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The analyse phase: elimination tree, postorder, column counts and
//! fundamental supernodes (lines 20-367); relaxation, pattern and the driver
//! are in the sibling files.

pub mod pattern;
pub mod relax;
pub mod run;

use super::super::auxiliary::logger::Logger;
use super::super::auxiliary::{
    childrenLinkedList, dfsPostorder, inversePerm, permuteSym, permuteVector, processEdge, snFlops,
    snSpops, subtreeSize, transpose,
};
use super::super::{Int, Int64};
use super::data_collector::DataCollector;
use super::symbolic::Symbolic;
use super::*;

/// HiGHS: Analyse.h `class Analyse`
pub struct Analyse<'a> {
    pub ready_: bool,
    pub rows_upper_: Vec<Int>,
    pub ptr_upper_: Vec<Int>,
    pub rows_lower_: Vec<Int>,
    pub ptr_lower_: Vec<Int>,
    pub n_: Int,
    pub nz_: Int,
    pub nz_factor_: Int64,
    pub dense_ops_: f64,
    pub dense_ops_norelax_: f64,
    pub sparse_ops_: f64,
    pub critical_ops_: f64,
    pub signs_: Vec<Int>,
    pub perm_: Vec<Int>,
    pub iperm_: Vec<Int>,
    pub parent_: Vec<Int>,
    pub postorder_: Vec<Int>,
    pub col_count_: Vec<Int>,
    pub rows_sn_: Vec<Int>,
    pub ptr_sn_: Vec<Int64>,
    pub sn_indices_: Vec<Int>,
    pub sn_count_: Int,
    pub artificial_nz_: Int64,
    pub sn_belong_: Vec<Int>,
    pub sn_start_: Vec<Int>,
    pub sn_parent_: Vec<Int>,
    pub fake_nz_: Vec<Int64>,
    pub merged_into_: Vec<Int>,
    pub merged_sn_: Int,
    pub relind_cols_: Vec<Int>,
    pub relind_clique_: Vec<Vec<Int>>,
    pub consecutive_sums_: Vec<Vec<Int>>,
    pub serial_storage_: f64,
    pub max_stack_size_: Int64,
    pub clique_block_start_: Vec<Vec<Int64>>,
    pub nb_: Int,
    pub logger_: Option<&'a Logger>,
    pub data_: &'a DataCollector,
}

impl<'a> Analyse<'a> {
    /// HiGHS: Analyse.cpp::Analyse::Analyse (lines 20-63)
    pub fn new(
        rows: &[Int],
        ptr: &[Int],
        signs: &[Int],
        nb: Int,
        logger: Option<&'a Logger>,
        data: &'a DataCollector,
        perm: &[Int],
    ) -> Self {
        let mut obj = Analyse {
            ready_: false,
            rows_upper_: Vec::new(),
            ptr_upper_: Vec::new(),
            rows_lower_: Vec::new(),
            ptr_lower_: Vec::new(),
            n_: 0,
            nz_: 0,
            nz_factor_: 0,
            dense_ops_: 0.0,
            dense_ops_norelax_: 0.0,
            sparse_ops_: 0.0,
            critical_ops_: 0.0,
            signs_: Vec::new(),
            perm_: Vec::new(),
            iperm_: Vec::new(),
            parent_: Vec::new(),
            postorder_: Vec::new(),
            col_count_: Vec::new(),
            rows_sn_: Vec::new(),
            ptr_sn_: Vec::new(),
            sn_indices_: Vec::new(),
            sn_count_: 0,
            artificial_nz_: 0,
            sn_belong_: Vec::new(),
            sn_start_: Vec::new(),
            sn_parent_: Vec::new(),
            fake_nz_: Vec::new(),
            merged_into_: Vec::new(),
            merged_sn_: 0,
            relind_cols_: Vec::new(),
            relind_clique_: Vec::new(),
            consecutive_sums_: Vec::new(),
            serial_storage_: 0.0,
            max_stack_size_: 0,
            clique_block_start_: Vec::new(),
            nb_: 0,
            logger_: logger,
            data_: data,
        };

        // Input the symmetric matrix to be analysed in CSC format.
        // rows contains the row indices.
        // ptr contains the starting points of each column.
        // Only the lower triangular part is used.
        // signs contains the sign that each pivot should have.

        obj.n_ = ptr.len() as Int - 1;
        obj.nz_ = rows.len() as Int;
        obj.signs_ = signs.to_vec();
        obj.nb_ = nb;

        // Create upper triangular part
        obj.rows_upper_.resize(obj.nz_ as usize, 0);
        obj.ptr_upper_.resize(obj.n_ as usize + 1, 0);
        transpose(ptr, rows, &mut obj.ptr_upper_, &mut obj.rows_upper_);

        // Permute the matrix with identical permutation, to extract upper triangular
        // part, if the input is not lower triangular.
        let mut id_perm = vec![0 as Int; obj.n_ as usize];
        for i in 0..obj.n_ {
            id_perm[i as usize] = i;
        }
        obj.permute(&id_perm);

        // actual number of nonzeros of only upper triangular part
        obj.nz_ = obj.ptr_upper_[obj.ptr_upper_.len() - 1];

        // number of nonzeros potentially changed after Permute.
        obj.rows_upper_.resize(obj.nz_ as usize, 0);

        // double transpose to sort columns
        obj.ptr_lower_.resize(obj.n_ as usize + 1, 0);
        obj.rows_lower_.resize(obj.nz_ as usize, 0);
        transpose(
            &obj.ptr_upper_,
            &obj.rows_upper_,
            &mut obj.ptr_lower_,
            &mut obj.rows_lower_,
        );
        transpose(
            &obj.ptr_lower_,
            &obj.rows_lower_,
            &mut obj.ptr_upper_,
            &mut obj.rows_upper_,
        );

        obj.perm_ = perm.to_vec();
        obj.iperm_.resize(obj.n_ as usize, 0);
        inversePerm(&obj.perm_, &mut obj.iperm_);

        obj.ready_ = true;
        obj
    }

    /// HiGHS: Analyse.cpp::Analyse::permute (lines 65-68)
    pub fn permute(&mut self, iperm: &[Int]) {
        let mut empty_vals: Vec<f64> = Vec::new();
        permuteSym(
            iperm,
            &mut self.ptr_upper_,
            &mut self.rows_upper_,
            &mut empty_vals,
            false,
        );
    }

    /// HiGHS: Analyse.cpp::Analyse::eTree (lines 70-100)
    pub fn eTree(&mut self) {
        // Find elimination tree.
        // It works only for upper triangular matrices.
        // The tree is stored in the vector parent:
        //  parent[i] = j
        // means that j is the parent of i in the tree.
        // For the root(s) of the tree, parent[root] = -1.

        self.parent_.resize(self.n_ as usize, 0);
        let mut ancestor = vec![0 as Int; self.n_ as usize];
        let mut next: Int;

        for j in 0..self.n_ {
            // initialise parent and ancestor, which are still unknown
            self.parent_[j as usize] = -1;
            ancestor[j as usize] = -1;

            let mut el = self.ptr_upper_[j as usize];
            while el < self.ptr_upper_[(j + 1) as usize] {
                let mut i = self.rows_upper_[el as usize];
                while i != -1 && i < j {
                    // next is used to move up the tree
                    next = ancestor[i as usize];

                    // ancestor keeps track of the known part of the tree, to avoid
                    // repeating (aka path compression): from j there is a known path to i
                    ancestor[i as usize] = j;

                    if next == -1 {
                        self.parent_[i as usize] = j;
                    }
                    i = next;
                }
                el += 1;
            }
        }
    }

    /// HiGHS: Analyse.cpp::Analyse::postorder (lines 102-142)
    pub fn postorder(&mut self) {
        // Find a postordering of the elimination tree using depth first search

        self.postorder_.resize(self.n_ as usize, 0);

        // create linked list of children
        let mut head: Vec<Int> = Vec::new();
        let mut next: Vec<Int> = Vec::new();
        childrenLinkedList(&self.parent_, &mut head, &mut next);

        // Execute depth first search only for root node(s)
        let mut start: Int = 0;
        for node in 0..self.n_ {
            if self.parent_[node as usize] == -1 {
                dfsPostorder(node, &mut start, &mut head, &next, &mut self.postorder_);
            }
        }

        // Permute elimination tree based on postorder
        let mut ipost = vec![0 as Int; self.n_ as usize];
        inversePerm(&self.postorder_, &mut ipost);
        let mut new_parent = vec![0 as Int; self.n_ as usize];
        for i in 0..self.n_ {
            if self.parent_[i as usize] != -1 {
                new_parent[ipost[i as usize] as usize] = ipost[self.parent_[i as usize] as usize];
            } else {
                new_parent[ipost[i as usize] as usize] = -1;
            }
        }
        self.parent_ = new_parent;

        // Permute matrix based on postorder
        self.permute(&ipost);

        // double transpose to sort columns and compute lower part
        transpose(
            &self.ptr_upper_,
            &self.rows_upper_,
            &mut self.ptr_lower_,
            &mut self.rows_lower_,
        );
        transpose(
            &self.ptr_lower_,
            &self.rows_lower_,
            &mut self.ptr_upper_,
            &mut self.rows_upper_,
        );

        // Update perm and iperm
        permuteVector(&mut self.perm_, &self.postorder_);
        inversePerm(&self.perm_, &mut self.iperm_);
    }

    /// HiGHS: Analyse.cpp::Analyse::colCount (lines 144-197)
    pub fn colCount(&mut self) {
        // Columns count using skeleton matrix.
        // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

        let mut first = vec![-1 as Int; self.n_ as usize];
        let mut ancestor = vec![-1 as Int; self.n_ as usize];
        let mut max_first = vec![-1 as Int; self.n_ as usize];
        let mut prev_leaf = vec![-1 as Int; self.n_ as usize];

        self.col_count_.resize(self.n_ as usize, 0);

        // find first descendant
        for k in 0..self.n_ {
            let mut j = k;
            self.col_count_[j as usize] = if first[j as usize] == -1 { 1 } else { 0 };
            while j != -1 && first[j as usize] == -1 {
                first[j as usize] = k;
                j = self.parent_[j as usize];
            }
        }

        // each node belongs to a separate set
        for j in 0..self.n_ {
            ancestor[j as usize] = j;
        }

        for k in 0..self.n_ {
            let j = k;

            // if not a root, decrement
            if self.parent_[j as usize] != -1 {
                self.col_count_[self.parent_[j as usize] as usize] -= 1;
            }

            // process edges of matrix
            let mut el = self.ptr_lower_[j as usize];
            while el < self.ptr_lower_[(j + 1) as usize] {
                processEdge(
                    j,
                    self.rows_lower_[el as usize],
                    &first,
                    &mut max_first,
                    &mut self.col_count_,
                    &mut prev_leaf,
                    &mut ancestor,
                );
                el += 1;
            }

            if self.parent_[j as usize] != -1 {
                ancestor[j as usize] = self.parent_[j as usize];
            }
        }

        // sum contributions from each child
        for j in 0..self.n_ {
            if self.parent_[j as usize] != -1 {
                self.col_count_[self.parent_[j as usize] as usize] += self.col_count_[j as usize];
            }
        }

        // compute nonzeros of L
        self.dense_ops_norelax_ = 0.0;
        self.nz_factor_ = 0;
        for j in 0..self.n_ {
            self.nz_factor_ += self.col_count_[j as usize] as Int64;
            self.dense_ops_norelax_ +=
                (self.col_count_[j as usize] - 1) as f64 * (self.col_count_[j as usize] - 1) as f64;
        }
    }

    /// HiGHS: Analyse.cpp::Analyse::fundamentalSupernodes (lines 199-267)
    pub fn fundamentalSupernodes(&mut self) {
        // Find fundamental supernodes.

        // isSN[i] is true if node i is the start of a fundamental supernode
        let mut is_sn = vec![false; self.n_ as usize];

        let mut prev_nonz = vec![-1 as Int; self.n_ as usize];

        // compute sizes of subtrees
        let mut subtree_sizes = vec![0 as Int; self.n_ as usize];
        subtreeSize(&self.parent_, &mut subtree_sizes);

        for j in 0..self.n_ {
            let mut el = self.ptr_lower_[j as usize];
            while el < self.ptr_lower_[(j + 1) as usize] {
                let i = self.rows_lower_[el as usize];
                let k = prev_nonz[i as usize];

                // mark as fundamental sn, nodes which are leaf of subtrees
                if k < j - subtree_sizes[j as usize] + 1 {
                    is_sn[j as usize] = true;
                }

                // mark as fundamental sn, nodes which have more than one child
                if self.parent_[i as usize] != -1
                    && subtree_sizes[i as usize] + 1
                        != subtree_sizes[self.parent_[i as usize] as usize]
                {
                    is_sn[self.parent_[i as usize] as usize] = true;
                }

                prev_nonz[i as usize] = j;
                el += 1;
            }
        }

        // create information about fundamental supernodes
        self.sn_belong_.resize(self.n_ as usize, 0);
        let mut sn_number: Int = -1;
        for i in 0..self.n_ {
            // if isSN[i] is true, then node i is the start of a new supernode
            if is_sn[i as usize] {
                sn_number += 1;
            }

            // mark node i as belonging to the current supernode
            self.sn_belong_[i as usize] = sn_number;
        }

        // number of supernodes found
        self.sn_count_ = self.sn_belong_[self.sn_belong_.len() - 1] + 1;

        // sn_start_ contains pointers to the starting node of each supernode
        self.sn_start_.resize(self.sn_count_ as usize + 1, 0);
        let mut next: Int = 0;
        for i in 0..self.n_ {
            if is_sn[i as usize] {
                self.sn_start_[next as usize] = i;
                next += 1;
            }
        }
        self.sn_start_[next as usize] = self.n_;

        // build supernodal elimination tree
        self.sn_parent_.resize(self.sn_count_ as usize, 0);
        for i in 0..self.sn_count_ - 1 {
            let j = self.parent_[(self.sn_start_[(i + 1) as usize] - 1) as usize];
            if j != -1 {
                self.sn_parent_[i as usize] = self.sn_belong_[j as usize];
            } else {
                self.sn_parent_[i as usize] = -1;
            }
        }
        let last = self.sn_parent_.len() - 1;
        self.sn_parent_[last] = -1;
    }

    /// HiGHS: Analyse.cpp::Analyse::relaxSupernodes (lines 269-292)
    pub fn relaxSupernodes(&mut self) {
        // Choose which strategy to use for supernode amalgamation.
        // Relaxing supernodes is fast to run, so just run again the one that is
        // chosen to set the correct information.

        let mut flops_1: f64 = 0.0;
        let mut spops_1: f64 = 0.0;
        self.relaxSnMaxNz(&mut flops_1, &mut spops_1);
        let ops_1 = flops_1 + kSpopsWeightSn * spops_1;

        let mut flops_2: f64 = 0.0;
        let mut spops_2: f64 = 0.0;
        self.relaxSnNetOps(&mut flops_2, &mut spops_2);
        let ops_2 = flops_2 + kSpopsWeightSn * spops_2;

        let total_ops_reduced = ops_2 < ops_1;
        let spops_reduced_enough = spops_2 < spops_1 * 0.6;

        if total_ops_reduced && spops_reduced_enough {
            self.relaxSnNetOps(&mut flops_2, &mut spops_2);
        } else {
            self.relaxSnMaxNz(&mut flops_1, &mut spops_1);
        }
    }

    /// HiGHS: Analyse.cpp::Analyse::relaxSnNetOps (lines 294-367)
    pub fn relaxSnNetOps(&mut self, flops: &mut f64, spops: &mut f64) {
        // =================================================
        // Build information about supernodes
        // =================================================
        let mut sn_size = vec![0 as Int; self.sn_count_ as usize];
        let mut clique_size = vec![0 as Int; self.sn_count_ as usize];
        // C++ fake_nz_.assign(sn_count_, 0): overwrite, not resize
        self.fake_nz_ = vec![0 as Int64; self.sn_count_ as usize];
        for i in 0..self.sn_count_ {
            sn_size[i as usize] = self.sn_start_[(i + 1) as usize] - self.sn_start_[i as usize];
            clique_size[i as usize] =
                self.col_count_[self.sn_start_[i as usize] as usize] - sn_size[i as usize];
        }

        // build linked lists of children
        let mut first_child: Vec<Int> = Vec::new();
        let mut next_child: Vec<Int> = Vec::new();
        childrenLinkedList(&self.sn_parent_, &mut first_child, &mut next_child);

        // =================================================
        // Merge supernodes
        // =================================================
        // C++ merged_into_.assign(sn_count_, -1): overwrite, not resize
        self.merged_into_ = vec![-1 as Int; self.sn_count_ as usize];
        self.merged_sn_ = 0;

        for sn in 0..self.sn_count_ {
            // sort children by increasing size
            let mut children: Vec<Int> = Vec::new();
            let mut child = first_child[sn as usize];
            while child != -1 {
                children.push(child);
                child = next_child[child as usize];
            }
            // PORT NOTE: std::sort on sn_size, keys can tie (equal-size siblings).
            // The oracle's libc++ std::sort uses sort networks / insertion sort
            // (both stable) below 24 elements, so a stable sort_by matches it
            // there; with 24+ children and tied sizes the C++ order is
            // implementation-defined: a known, reported divergence.
            children.sort_by(|a, b| sn_size[*a as usize].cmp(&sn_size[*b as usize]));

            for &child in &children {
                // how many zero rows would become nonzero
                let rows_filled =
                    sn_size[sn as usize] + clique_size[sn as usize] - clique_size[child as usize];

                // how many zero entries would become nonzero
                let nz_added = rows_filled as Int64 * sn_size[child as usize] as Int64;

                // how many artificial nonzeros would the merged supernode have
                let total_art_nz =
                    nz_added + self.fake_nz_[sn as usize] + self.fake_nz_[child as usize];

                // C++ int arithmetic: the whole sum is evaluated in Int, then converted
                let flops_added = (rows_filled * rows_filled * sn_size[child as usize]
                    + 2 * rows_filled
                        * (sn_size[child as usize]
                            * (sn_size[child as usize] + clique_size[child as usize])
                            - sn_size[child as usize] * (sn_size[child as usize] + 1) / 2))
                    as f64;

                // C++ int arithmetic: clique_size * (clique_size + 1) / 2 is an integer division
                let spops_removed =
                    (clique_size[child as usize] * (clique_size[child as usize] + 1) / 2) as f64;

                let net_ops = flops_added - kSpopsWeightSn * spops_removed;

                if net_ops < 0.0 {
                    // merge child with parent
                    sn_size[sn as usize] += sn_size[child as usize];
                    self.fake_nz_[sn as usize] = total_art_nz;
                    self.merged_sn_ += 1;
                    self.merged_into_[child as usize] = sn;
                }
            }
        }

        *flops = 0.0;
        *spops = 0.0;
        for i in 0..self.sn_count_ {
            if self.merged_into_[i as usize] == -1 {
                *spops += snSpops(clique_size[i as usize] as f64);
                *flops += snFlops(sn_size[i as usize] as f64, clique_size[i as usize] as f64);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn etree_of_arrow() {
        // 3x3 arrow, lower pattern: col 0 {0,2}, col 1 {1,2}, col 2 {2}
        let d = DataCollector::new();
        let mut a = Analyse::new(
            &[0, 2, 1, 2, 2],
            &[0, 2, 4, 5],
            &[1, 1, 1],
            128,
            None,
            &d,
            &[0, 1, 2],
        );
        a.eTree();
        assert_eq!(a.parent_, vec![2, 2, -1]);
    }
}
