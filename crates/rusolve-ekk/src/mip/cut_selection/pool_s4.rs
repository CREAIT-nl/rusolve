//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! candidate collection and commit for the hybrid cut selector.

use crate::mip::PoolEvent;
use crate::mip::cut_pool::HighsCutPool;
use crate::mip::cut_pool::HighsCutSet;
use crate::mip::cut_pool::compute_cut_hash;
use crate::mip::domain::HighsDomain;
use crate::prelude::*;
use crate::util::highs_cdouble::HighsCDouble;
use crate::util::highs_hash::hash_u64_bytes;
use crate::util::highs_pdqsort::pdqsort_by;
use crate::util::highs_utils::cpp_max;

use std::sync::atomic::Ordering;

impl HighsCutPool {
    /// HiGHS: HighsCutPool.cpp::HighsCutPool::separate (lines 229-392)
    /// SCIP 11.0.0 (32dd425) candidate collection: augment of HighsCutPool.cpp::separate, returns the truncated list (lines 229-392)
    #[allow(dead_code)] // not used yet
    pub(crate) fn s4_collect_candidates(
        &mut self,
        sol: &[f64],
        domain: &HighsDomain,
        feastol: f64,
        thread_safe: bool,
    ) -> Vec<(f64, HighsInt)> {
        let nrows = self.matrix_.get_num_rows();

        let mut efficacious_cuts: Vec<(f64, HighsInt)> = Vec::new();

        let mut agelim = self.agelim_;

        let mut num_cuts = self.get_num_cuts() - self.num_lp_cuts;
        while agelim > 1 && num_cuts > self.softlimit_ {
            num_cuts -= self.age_distribution[agelim as usize];
            agelim -= 1;
        }

        for i in 0..nrows {
            let iu = i as usize;
            // cuts with an age of -1 are already in the LP and are therefore skipped
            // Warning: Parallel case tries to add cuts already in current LP.
            // Inefficient. Not sure what happens if added twice.
            // The cut shouldn't have enough violation to be added though.
            if self.ages_[iu] < 0 {
                continue;
            }

            let start = self.matrix_.get_row_start(i);
            let end = self.matrix_.get_row_end(i);

            let mut viol = -self.rhs_[iu];

            {
                let rindex = &self.matrix_.get_a_rindex()[start as usize..end as usize];
                let rvalue = &self.matrix_.get_a_rvalue()[start as usize..end as usize];
                for (&col, &v) in rindex.iter().zip(rvalue) {
                    let solval = sol[col as usize];

                    viol += v * solval;
                }
            }

            // if the cut is not violated more than feasibility tolerance
            // we skip it and increase its age, otherwise we reset its age
            let is_propagated = self.matrix_.columns_linked(i);
            if !thread_safe {
                self.age_distribution[self.ages_[iu] as usize] -= 1;
                if is_propagated {
                    self.prop_rows.remove(&(self.ages_[iu] as HighsInt, i));
                }
            }
            if viol <= feastol {
                if thread_safe {
                    continue;
                }
                self.ages_[iu] += 1;
                if (self.ages_[iu] as HighsInt) >= agelim {
                    let h = compute_cut_hash(
                        &self.matrix_.get_a_rindex()[start as usize..end as usize],
                        &self.matrix_.get_a_rvalue()[start as usize..end as usize],
                        self.maxabscoef_[iu],
                        end - start,
                    );

                    // C++: for (p : propagationDomains) p->cutDeleted(i);
                    // PORT NOTE: see perform_aging (no observers -> no-op push skipped).
                    if !self.events.observers.is_empty() {
                        self.events.push(PoolEvent::CutDeleted {
                            cut: i,
                            only_for_prop: false,
                        });
                    }

                    if is_propagated {
                        self.num_prop_rows -= 1;
                        self.num_prop_nzs -= self.get_row_length(i);
                    }

                    self.matrix_.remove_row(i);
                    self.ages_[iu] = -1;
                    self.rhs_[iu] = kHighsInf;
                    self.age_reset_while_locked_[iu].store(0, Ordering::Relaxed);
                    self.has_synced_[iu] = false;
                    let range_count = self.hash_to_cut_map.count(h);

                    for it in 0..range_count {
                        if self.hash_to_cut_map.get(h, it) == i {
                            self.hash_to_cut_map.erase(h, it);
                            break;
                        }
                    }
                } else {
                    if is_propagated {
                        self.prop_rows.insert((self.ages_[iu] as HighsInt, i));
                    }
                    self.age_distribution[self.ages_[iu] as usize] += 1;
                }
                continue;
            }

            // compute the norm only for those entries that do not sit at their minimal
            // activity in the current solution this avoids the phenomenon that the
            // traditional efficacy gets weaker for stronger cuts E.g. when considering
            // a clique cut which has additional entries whose value in the current
            // solution is 0 then the efficacy gets lower for each such entry even
            // though the cut dominates the clique cut where all those entries are
            // relaxed out.
            let mut rownorm = HighsCDouble::from(0.0);
            let mut num_active_nzs: HighsInt = 0;
            let rindex = &self.matrix_.get_a_rindex()[start as usize..end as usize];
            let rvalue = &self.matrix_.get_a_rvalue()[start as usize..end as usize];
            for (&col, &v) in rindex.iter().zip(rvalue) {
                let col = col as usize;
                let solval = sol[col];
                if v > 0.0 {
                    if solval > domain.col_lower_[col] + feastol {
                        rownorm += v * v;
                        num_active_nzs += 1;
                    }
                } else if solval < domain.col_upper_[col] - feastol {
                    rownorm += v * v;
                    num_active_nzs += 1;
                }
            }

            if !thread_safe {
                self.ages_[iu] = 0;
                self.age_distribution[0] += 1;
                if is_propagated {
                    self.prop_rows.insert((self.ages_[iu] as HighsInt, i));
                }
            }
            let score = viol / ((num_active_nzs as f64) * f64::from(rownorm).sqrt());

            efficacious_cuts.push((score, i));
        }
        debug_assert!(self.prop_rows.len() as HighsInt == self.num_prop_rows);
        if efficacious_cuts.is_empty() {
            return Vec::new();
        }

        // C++ pdqsort comparator; `efficacious_cuts.size()` is captured before the sort.
        let size = efficacious_cuts.len() as u64;
        let hash_key = |cut_idx: HighsInt| -> u64 {
            let x = ((cut_idx as u64) << 32).wrapping_add(size);
            hash_u64_bytes(x as u32, (x >> 32) as u32)
        };
        pdqsort_by(
            &mut efficacious_cuts,
            |a: &(f64, HighsInt), b: &(f64, HighsInt)| {
                if a.0 > b.0 {
                    return true;
                }
                if a.0 < b.0 {
                    return false;
                }
                (hash_key(a.1), a.1) > (hash_key(b.1), b.1)
            },
        );

        // PORT NOTE: C++ `bestObservedScore_` / `minScoreFactor_` alias either the field
        // or a copy under `thread_safe`. The rest of this function does not read the
        // fields, so the local is written back to the field when `!thread_safe`.
        let best_observed_score_ = cpp_max(efficacious_cuts[0].0, self.best_observed_score);
        if !thread_safe {
            self.best_observed_score = best_observed_score_;
        }
        let mut min_score_factor_ = self.min_score_factor;
        let min_score = min_score_factor_ * best_observed_score_;

        // C++ std::upper_bound(begin, end, minScore, comp(mscore, c) = mscore > c.first)
        let mut numefficacious: HighsInt = {
            let mut first = 0usize;
            let mut count = efficacious_cuts.len();
            while count > 0 {
                let step = count / 2;
                let it = first + step;
                if !(min_score > efficacious_cuts[it].0) {
                    first = it + 1;
                    count -= step + 1;
                } else {
                    count = step;
                }
            }
            first as HighsInt
        };

        let lower_threshold = (efficacious_cuts.len() / 20) as HighsInt;
        let upper_threshold = (efficacious_cuts.len() - 1) as HighsInt;

        if numefficacious <= lower_threshold {
            numefficacious = cpp_max(efficacious_cuts.len() / 2, 1usize) as HighsInt;
            min_score_factor_ =
                efficacious_cuts[(numefficacious - 1) as usize].0 / best_observed_score_;
        } else if numefficacious > upper_threshold {
            min_score_factor_ = efficacious_cuts[upper_threshold as usize].0 / best_observed_score_;
        }
        if !thread_safe {
            self.min_score_factor = min_score_factor_;
        }

        efficacious_cuts.truncate(numefficacious as usize);

        efficacious_cuts
    }

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::separate (lines 422-434, 394-396, 437-459)
    /// SCIP 11.0.0 (32dd425) commit of the selected cuts (augment of separate)
    #[allow(dead_code)] // not used yet
    pub(crate) fn s4_commit_selected(
        &mut self,
        cutset: &mut HighsCutSet,
        rows: &[HighsInt],
        thread_safe: bool,
    ) {
        let orignumcuts = cutset.num_cuts();
        let origselectednnz = cutset.a_rindex_.len() as HighsInt;
        let mut selectednnz = origselectednnz;

        for &p in rows.iter() {
            let pu = p as usize;
            self.num_lps_[pu].fetch_add(1, Ordering::Relaxed);
            if !thread_safe {
                self.age_distribution[self.ages_[pu] as usize] -= 1;
                self.num_lp_cuts += 1;
                if self.matrix_.columns_linked(p) {
                    self.prop_rows.remove(&(self.ages_[pu] as HighsInt, p));
                    self.prop_rows.insert((-1, p));
                }
                self.ages_[pu] = -1;
            }
            cutset.cutindices.push(p);
            cutset.cutpools.push(self.index_);
            selectednnz += self.matrix_.get_row_end(p) - self.matrix_.get_row_start(p);
        }

        cutset.resize(selectednnz);

        debug_assert!(cutset.a_rvalue_.len() as HighsInt == selectednnz);
        debug_assert!(cutset.a_rindex_.len() as HighsInt == selectednnz);

        let mut offset = origselectednnz;
        for i in orignumcuts..cutset.num_cuts() {
            let iu = i as usize;
            cutset.a_rstart_[iu] = offset;
            let cut = cutset.cutindices[iu];
            let start = self.matrix_.get_row_start(cut);
            let end = self.matrix_.get_row_end(cut);
            cutset.upper_[iu] = self.rhs_[cut as usize];

            for j in start..end {
                debug_assert!(offset < selectednnz);
                cutset.a_rvalue_[offset as usize] = self.matrix_.get_a_rvalue()[j as usize];
                cutset.a_rindex_[offset as usize] = self.matrix_.get_a_rindex()[j as usize];
                offset += 1;
            }
        }

        let ncuts = cutset.num_cuts() as usize;
        cutset.a_rstart_[ncuts] = offset;

        debug_assert!(self.prop_rows.len() as HighsInt == self.num_prop_rows);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::mip_solver::HighsMipSolver;

    fn build_fixture() -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 4;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0; 4];
        lp.col_lower_ = vec![0.0; 4];
        lp.col_upper_ = vec![1.0; 4];
        lp.integrality_ = vec![HighsVarType::kInteger; 4];
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0; 5];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms
    }

    fn add_three_cuts(ms: &mut HighsMipSolver, pool: &mut HighsCutPool) {
        // cut1: -x0 <= -0.8
        let mut rindex1 = vec![0 as HighsInt];
        let mut rvalue1 = vec![-1.0];
        pool.add_cut(
            ms,
            &mut rindex1,
            &mut rvalue1,
            1,
            -0.8,
            false,
            false,
            false,
            false,
        );

        // cut2: -x1 <= -0.9
        let mut rindex2 = vec![1 as HighsInt];
        let mut rvalue2 = vec![-1.0];
        pool.add_cut(
            ms,
            &mut rindex2,
            &mut rvalue2,
            1,
            -0.9,
            false,
            false,
            false,
            false,
        );

        // cut3: x2 + x3 <= 0.2
        let mut rindex3 = vec![2 as HighsInt, 3 as HighsInt];
        let mut rvalue3 = vec![1.0, 1.0];
        pool.add_cut(
            ms,
            &mut rindex3,
            &mut rvalue3,
            2,
            0.2,
            false,
            false,
            false,
            false,
        );
    }

    #[test]
    fn s4_pool_collect_commit_matches_separate() {
        let feastol = 1e-6;
        let sol = vec![0.5, 0.5, 0.5, 0.5];

        // Pool A: use separate
        let mut ms_a = build_fixture();
        let mut pool_a = HighsCutPool::new(4, 10, 100, 0);
        add_three_cuts(&mut ms_a, &mut pool_a);
        let dom_a = HighsDomain::new(&ms_a);
        let mut cutset_a = HighsCutSet::default();
        pool_a.separate(&sol, &dom_a, &mut cutset_a, feastol, &[], false);
        assert!(
            !cutset_a.cutindices.is_empty(),
            "equivalence test must select cuts"
        );

        // Pool B: use s4_collect_candidates + s4_commit_selected
        let mut ms_b = build_fixture();
        let mut pool_b = HighsCutPool::new(4, 10, 100, 0);
        add_three_cuts(&mut ms_b, &mut pool_b);
        let dom_b = HighsDomain::new(&ms_b);
        let rows: Vec<HighsInt> = pool_b
            .s4_collect_candidates(&sol, &dom_b, feastol, false)
            .iter()
            .map(|&(_, p)| p)
            .collect();
        let mut cutset_b = HighsCutSet::default();
        pool_b.s4_commit_selected(&mut cutset_b, &rows, false);

        assert_eq!(cutset_a.cutindices, cutset_b.cutindices);
        assert_eq!(cutset_a.cutpools, cutset_b.cutpools);
        assert_eq!(cutset_a.a_rstart_, cutset_b.a_rstart_);
        assert_eq!(cutset_a.a_rindex_, cutset_b.a_rindex_);
        assert_eq!(cutset_a.a_rvalue_, cutset_b.a_rvalue_);
        assert_eq!(cutset_a.upper_, cutset_b.upper_);
        assert_eq!(cutset_a.lower_, cutset_b.lower_);

        assert_eq!(pool_a.ages_, pool_b.ages_);
        assert_eq!(pool_a.age_distribution, pool_b.age_distribution);
        assert_eq!(pool_a.num_lp_cuts, pool_b.num_lp_cuts);
        assert_eq!(pool_a.min_score_factor, pool_b.min_score_factor);
        assert_eq!(pool_a.best_observed_score, pool_b.best_observed_score);
        for i in 0..pool_a.num_lps_.len() {
            assert_eq!(
                pool_a.num_lps_[i].load(Ordering::Relaxed),
                pool_b.num_lps_[i].load(Ordering::Relaxed)
            );
        }
    }

    #[test]
    fn s4_pool_thread_safe_leaves_ages() {
        let feastol = 1e-6;
        let sol = vec![0.5, 0.5, 0.5, 0.5];

        let mut ms = build_fixture();
        let mut pool = HighsCutPool::new(4, 10, 100, 0);
        add_three_cuts(&mut ms, &mut pool);
        let dom = HighsDomain::new(&ms);

        // Snapshot ages and age_distribution before collect
        let ages_before: Vec<i16> = pool.ages_.clone();
        let age_dist_before: Vec<HighsInt> = pool.age_distribution.clone();
        let num_lps_before: Vec<i16> = pool
            .num_lps_
            .iter()
            .map(|x| x.load(Ordering::Relaxed))
            .collect();

        let num_lp_cuts_before = pool.num_lp_cuts;
        let best_before = pool.best_observed_score;
        let msf_before = pool.min_score_factor;

        let rows: Vec<HighsInt> = pool
            .s4_collect_candidates(&sol, &dom, feastol, true)
            .iter()
            .map(|&(_, p)| p)
            .collect();
        assert!(!rows.is_empty(), "thread_safe test must select cuts");
        let mut cutset = HighsCutSet::default();
        pool.s4_commit_selected(&mut cutset, &rows, true);
        assert_eq!(cutset.cutindices, rows);

        // With thread_safe = true, ages_ and age_distribution must be unchanged
        assert_eq!(pool.ages_, ages_before);
        assert_eq!(pool.age_distribution, age_dist_before);
        assert_eq!(pool.num_lp_cuts, num_lp_cuts_before);
        assert_eq!(pool.best_observed_score, best_before);
        assert_eq!(pool.min_score_factor, msf_before);

        // num_lps_ is incremented exactly once for each committed cut
        for i in 0..pool.num_lps_.len() {
            let expected = num_lps_before[i] + rows.contains(&(i as HighsInt)) as i16;
            assert_eq!(pool.num_lps_[i].load(Ordering::Relaxed), expected);
        }
    }

    #[test]
    fn s4_pool_no_violated_cut_empty() {
        let feastol = 1e-6;
        // sol = [1.0, 1.0, 0.0, 0.0]:
        //   cut1: -x0 <= -0.8 => -1.0 <= -0.8 => satisfied (viol = -(-0.8) + (-1)*1 = 0.8 - 1 = -0.2 <= feastol)
        //   cut2: -x1 <= -0.9 => -1.0 <= -0.9 => satisfied (viol = -(-0.9) + (-1)*1 = 0.9 - 1 = -0.1 <= feastol)
        //   cut3: x2+x3 <= 0.2 => 0+0 <= 0.2 => satisfied (viol = -(0.2) + 1*0 + 1*0 = -0.2 <= feastol)
        // All cuts have viol <= feastol, so none are efficacious.
        let sol = vec![1.0, 1.0, 0.0, 0.0];

        let mut ms = build_fixture();
        let mut pool = HighsCutPool::new(4, 10, 100, 0);
        add_three_cuts(&mut ms, &mut pool);
        let dom = HighsDomain::new(&ms);

        let rows: Vec<HighsInt> = pool
            .s4_collect_candidates(&sol, &dom, feastol, false)
            .iter()
            .map(|&(_, p)| p)
            .collect();
        assert!(rows.is_empty(), "expected no efficacious cuts");

        let mut cutset = HighsCutSet::default();
        pool.s4_commit_selected(&mut cutset, &rows, false);
        assert!(cutset.cutindices.is_empty());
        assert!(cutset.cutpools.is_empty());
    }
}
