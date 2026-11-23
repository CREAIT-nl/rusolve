//! HiGHS source: ipm/hipo/factorhighs/Analyse.cpp (lines 369-708)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Supernode relaxation (amalgamation) by artificial nonzeros and the
//! bookkeeping after it.

use super::super::super::{Int, Int64};
use super::super::*;
use super::Analyse;
use crate::ipm::hipo::auxiliary::{
    childrenLinkedList, inversePerm, permuteVector, snFlops, snSpops, transpose,
};
use crate::lp_data::hconst::kHighsInf;
use crate::util::highs_utils::cpp_min;

impl<'a> Analyse<'a> {
    /// HiGHS: Analyse.cpp::Analyse::doRelaxSnMaxNz (lines 369-486)
    pub fn doRelaxSnMaxNz(
        &mut self,
        flops: &mut f64,
        spops: &mut f64,
        max_artificial_nz: Int64,
    ) -> f64 {
        // HiGHS: Analyse.cpp:18
        const int64_limit: Int64 = i64::MAX;

        // =================================================
        // Build information about supernodes
        // =================================================
        let sn_count_ = self.sn_count_;
        let mut sn_size = vec![0 as Int; sn_count_ as usize];
        let mut clique_size = vec![0 as Int; sn_count_ as usize];
        self.fake_nz_ = vec![0 as Int64; sn_count_ as usize];
        for i in 0..sn_count_ {
            sn_size[i as usize] = self.sn_start_[(i + 1) as usize] - self.sn_start_[i as usize];
            clique_size[i as usize] =
                self.col_count_[self.sn_start_[i as usize] as usize] - sn_size[i as usize];
            self.fake_nz_[i as usize] = 0;
        }

        // build linked lists of children
        let mut first_child: Vec<Int> = Vec::new();
        let mut next_child: Vec<Int> = Vec::new();
        childrenLinkedList(&self.sn_parent_, &mut first_child, &mut next_child);

        // =================================================
        // Merge supernodes
        // =================================================
        self.merged_into_ = vec![-1 as Int; sn_count_ as usize];
        self.merged_sn_ = 0;

        for sn in 0..sn_count_ {
            // keep iterating through the children of the supernode, until there's no
            // more child to merge with

            loop {
                let mut child = first_child[sn as usize];

                // info for first criterion
                let mut nz_fakenz: Int64 = int64_limit;
                let mut size_fakenz: Int = 0;
                let mut child_fakenz: Int = -1;

                while child != -1 {
                    // how many zero rows would become nonzero
                    let rows_filled: Int = sn_size[sn as usize] + clique_size[sn as usize]
                        - clique_size[child as usize];

                    // how many zero entries would become nonzero
                    let nz_added: Int64 =
                        (rows_filled as Int64) * (sn_size[child as usize] as Int64);

                    // how many artificial nonzeros would the merged supernode have
                    let total_art_nz: Int64 =
                        nz_added + self.fake_nz_[sn as usize] + self.fake_nz_[child as usize];

                    // Save child with smallest number of artificial zeros created.
                    // Ties are broken based on size of child.
                    if total_art_nz < nz_fakenz
                        || (total_art_nz == nz_fakenz && size_fakenz < sn_size[child as usize])
                    {
                        nz_fakenz = total_art_nz;
                        size_fakenz = sn_size[child as usize];
                        child_fakenz = child;
                    }

                    child = next_child[child as usize];
                }

                if nz_fakenz <= max_artificial_nz {
                    // merging creates fewer nonzeros than the maximum allowed

                    // update information of parent
                    sn_size[sn as usize] += size_fakenz;
                    self.fake_nz_[sn as usize] = nz_fakenz;

                    // count number of merged supernodes
                    self.merged_sn_ += 1;

                    // save information about merging of supernodes
                    self.merged_into_[child_fakenz as usize] = sn;

                    // remove child from linked list of children
                    child = first_child[sn as usize];
                    if child == child_fakenz {
                        // child_smallest is the first child
                        first_child[sn as usize] = next_child[child_fakenz as usize];
                    } else {
                        while next_child[child as usize] != child_fakenz {
                            child = next_child[child as usize];
                        }
                        // now child is the previous child of child_smallest
                        next_child[child as usize] = next_child[child_fakenz as usize];
                    }
                } else {
                    // no more children can be merged with parent
                    break;
                }
            }
        }

        *flops = 0.0;
        *spops = 0.0;
        for i in 0..sn_count_ {
            if self.merged_into_[i as usize] == -1 {
                *spops += snSpops(clique_size[i as usize] as f64);
                *flops += snFlops(sn_size[i as usize] as f64, clique_size[i as usize] as f64);
            }
        }

        // compute total number of artificial nonzeros and artificial ops for this
        // value of max_artificial_nz
        let mut temp_art_nz: f64 = 0.0;
        let mut temp_art_ops: f64 = 0.0;
        for sn in 0..sn_count_ {
            if self.merged_into_[sn as usize] == -1 {
                temp_art_nz += self.fake_nz_[sn as usize] as f64;
                temp_art_ops +=
                    snFlops(sn_size[sn as usize] as f64, clique_size[sn as usize] as f64);
            }
        }
        // temp_art_nz is computed but unused in HiGHS as well
        let _ = temp_art_nz;
        temp_art_ops -= self.dense_ops_norelax_;

        // if enough fake nz or ops have been added, stop.
        let ratio_fake: f64 = temp_art_ops / (temp_art_ops + self.dense_ops_norelax_);

        ratio_fake
    }

    /// HiGHS: Analyse.cpp::Analyse::relaxSnMaxNz (lines 488-542)
    pub fn relaxSnMaxNz(&mut self, flops: &mut f64, spops: &mut f64) {
        // Child which produces smallest number of fake nonzeros is merged if
        // resulting sn has fewer than max_artificial_nz fake nonzeros.
        // Multiple values of max_artificial_nz are tried, chosen with bisection
        // method, until the percentage of artificial nonzeros is in the range [1,2]%.

        let mut max_artificial_nz: Int64 = kStartThreshRelax as Int64;
        let mut largest_below: Int64 = -1;
        let mut smallest_above: Int64 = -1;

        let mut best_dist_ratio: f64 = kHighsInf;
        let mut best_max_art_nz: Int64 = -1;

        for iter in 0..kMaxIterRelax {
            // relax the supernodes and obtain the ratio of how many new ops have been
            // added with the current value of max_artificial_nz
            let ratio_fake: f64 = self.doRelaxSnMaxNz(flops, spops, max_artificial_nz);

            // store the best ratio, in case a good ratio is never found
            let dist_ratio_fake: f64 = cpp_min(
                (ratio_fake - kLowerRatioRelax).abs(),
                (ratio_fake - kUpperRatioRelax).abs(),
            );
            if dist_ratio_fake < best_dist_ratio {
                best_dist_ratio = dist_ratio_fake;
                best_max_art_nz = max_artificial_nz;
            }

            // try to find ratio in interval [0.01,0.02] using bisection
            if ratio_fake < kLowerRatioRelax {
                // ratio too small
                largest_below = max_artificial_nz;
                if smallest_above == -1 {
                    max_artificial_nz *= 2;
                } else {
                    max_artificial_nz = (largest_below + smallest_above) / 2;
                }
            } else if ratio_fake > kUpperRatioRelax {
                // ratio too large
                smallest_above = max_artificial_nz;
                if largest_below == -1 {
                    max_artificial_nz /= 2;
                } else {
                    max_artificial_nz = (largest_below + smallest_above) / 2;
                }
            } else {
                // good ratio
                return;
            }
        }

        // If reach here, no good ratio was found within kMaxIterRelax
        // To avoid having a catastrophically bad ratio in pathological problems,
        // choose the best ratio found

        self.doRelaxSnMaxNz(flops, spops, best_max_art_nz);
    }

    /// HiGHS: Analyse.cpp::Analyse::afterRelaxSn (lines 544-708)
    pub fn afterRelaxSn(&mut self) {
        // number of new supernodes
        let new_snCount: Int = self.sn_count_ - self.merged_sn_;

        // keep track of number of row indices needed for each supernode
        self.sn_indices_ = vec![0 as Int; new_snCount as usize];

        // =================================================
        // Create supernodal permutation
        // =================================================

        // permutation of supernodes needed after merging
        let mut sn_perm = vec![0 as Int; self.sn_count_ as usize];

        // number of new sn that includes the old sn
        let mut new_id = vec![0 as Int; self.sn_count_ as usize];

        // new sn pointer vector
        let mut new_snStart = vec![0 as Int; (new_snCount + 1) as usize];

        // keep track of the children merged into a given supernode
        let mut received_from: Vec<Vec<Int>> = vec![Vec::new(); self.sn_count_ as usize];

        // index to write into sn_perm
        let mut start_perm: Int = 0;

        // index to write into new_snStart
        let mut snStart_ind: Int = 0;

        // next available number for new sn numbering
        let mut next_id: Int = 0;

        for sn in 0..self.sn_count_ {
            if self.merged_into_[sn as usize] > -1 {
                // Current sn was merged into its parent.
                // Save information about which supernode sn was merged into
                received_from[self.merged_into_[sn as usize] as usize].push(sn);
            } else {
                // Current sn was not merged into its parent.
                // It is one of the new sn.

                // Add merged supernodes to the permutation, recursively.

                snStart_ind += 1;

                let mut toadd: Vec<Int> = Vec::new();
                toadd.push(sn);

                while !toadd.is_empty() {
                    let current: Int = *toadd.last().unwrap();

                    if !received_from[current as usize].is_empty() {
                        for i in 0..received_from[current as usize].len() {
                            let val = received_from[current as usize][i];
                            toadd.push(val);
                        }
                        received_from[current as usize].clear();
                    } else {
                        toadd.pop();
                        sn_perm[start_perm as usize] = current;
                        start_perm += 1;
                        new_id[current as usize] = next_id;

                        // count number of nodes in each new supernode
                        new_snStart[snStart_ind as usize] += self.sn_start_[(current + 1) as usize]
                            - self.sn_start_[current as usize];
                    }
                }

                // keep track of total number of artificial nonzeros
                self.artificial_nz_ += self.fake_nz_[sn as usize];

                // Compute number of indices for new sn.
                // This is equal to the number of columns in the new sn plus the clique
                // size of the original supernode where the children where merged.
                self.sn_indices_[next_id as usize] = new_snStart[snStart_ind as usize]
                    + self.col_count_[self.sn_start_[sn as usize] as usize]
                    - self.sn_start_[(sn + 1) as usize]
                    + self.sn_start_[sn as usize];

                next_id += 1;
            }
        }

        // new_snStart contain the number of cols in each new sn.
        // sum them to obtain the sn pointers.
        for i in 0..new_snCount {
            new_snStart[(i + 1) as usize] += new_snStart[i as usize];
        }

        // include artificial nonzeros in the nonzeros of the factor
        self.nz_factor_ += self.artificial_nz_;

        // compute number of flops needed for the factorisation
        self.dense_ops_ = 0.0;
        for sn in 0..new_snCount {
            let colcount_sn: f64 = self.sn_indices_[sn as usize] as f64;
            for i in 0..(new_snStart[(sn + 1) as usize] - new_snStart[sn as usize]) {
                self.dense_ops_ += (colcount_sn - i as f64 - 1.0) * (colcount_sn - i as f64 - 1.0);
            }
        }

        // =================================================
        // Create nodal permutation
        // =================================================
        // Given the supernodal permutation, find the nodal permutation needed after
        // sn merging.

        // permutation to apply to the existing one
        let mut new_perm = vec![0 as Int; self.n_ as usize];

        // index to write into new_perm
        let mut start: Int = 0;

        for i in 0..self.sn_count_ {
            let sn: Int = sn_perm[i as usize];
            for j in self.sn_start_[sn as usize]..self.sn_start_[(sn + 1) as usize] {
                new_perm[start as usize] = j;
                start += 1;
            }
        }

        // obtain inverse permutation
        let mut new_iperm = vec![0 as Int; self.n_ as usize];
        inversePerm(&new_perm, &mut new_iperm);

        // =================================================
        // Create new sn elimination tree
        // =================================================
        let mut new_snParent = vec![-1 as Int; new_snCount as usize];
        for i in 0..self.sn_count_ {
            if self.sn_parent_[i as usize] == -1 {
                continue;
            }

            let ii: Int = new_id[i as usize];
            let pp: Int = new_id[self.sn_parent_[i as usize] as usize];

            if ii == pp {
                continue;
            }

            new_snParent[ii as usize] = pp;
        }

        // =================================================
        // Save new information
        // =================================================

        // build new snBelong, i.e., the sn to which each column belongs
        for sn in 0..self.sn_count_ {
            for i in self.sn_start_[sn as usize]..self.sn_start_[(sn + 1) as usize] {
                self.sn_belong_[i as usize] = new_id[sn as usize];
            }
        }
        permuteVector(&mut self.sn_belong_, &new_perm);

        permuteVector(&mut self.col_count_, &new_perm);

        // Overwrite previous data
        self.sn_parent_ = new_snParent;
        self.sn_start_ = new_snStart;
        self.sn_count_ = new_snCount;

        // Permute matrix based on new permutation
        self.permute(&new_iperm);

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
        permuteVector(&mut self.perm_, &new_perm);
        inversePerm(&self.perm_, &mut self.iperm_);
    }
}
