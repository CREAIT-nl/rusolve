//! HiGHS source: mip/HighsCutPool.h, mip/HighsCutPool.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//!
//! All bodies are ported. `compute_cut_hash` uses `vector_hash` from `util/highs_hash.rs`.

use crate::mip::{DomainKey, PoolEvent, Slot};
use crate::prelude::*;
use std::sync::atomic::{AtomicI16, AtomicU8, Ordering};

/// HiGHS: HighsCutPool.h::HighsCutSet (lines 22-52)
#[derive(Clone, Default, Debug)]
pub struct HighsCutSet {
    /// C++: `std::vector<HighsInt> cutindices`
    pub cutindices: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> cutpools`
    pub cutpools: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ARstart_`
    pub a_rstart_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ARindex_`
    pub a_rindex_: Vec<HighsInt>,
    /// C++: `std::vector<double> ARvalue_`
    pub a_rvalue_: Vec<f64>,
    /// C++: `std::vector<double> lower_`
    pub lower_: Vec<f64>,
    /// C++: `std::vector<double> upper_`
    pub upper_: Vec<f64>,
}

/// HiGHS: HighsCutPool.h::HighsCutPool (lines 54-194)
#[derive(Default, Debug)]
pub struct HighsCutPool {
    /// C++: `HighsDynamicRowMatrix matrix_`
    pub matrix_: HighsDynamicRowMatrix,
    /// C++: `std::vector<double> rhs_`
    pub rhs_: Vec<f64>,
    /// C++: `std::vector<int16_t> ages_`
    pub ages_: Vec<i16>,
    /// C++: `std::deque<std::atomic<int16_t>> numLps_`
    pub num_lps_: Vec<AtomicI16>,
    /// C++: `std::deque<std::atomic<uint8_t>> ageResetWhileLocked_`
    pub age_reset_while_locked_: Vec<AtomicU8>,
    /// C++: `std::vector<bool> hasSynced_`
    pub has_synced_: Vec<bool>,
    /// C++: `std::vector<double> rownormalization_`
    pub rownormalization_: Vec<f64>,
    /// C++: `std::vector<double> maxabscoef_`
    pub maxabscoef_: Vec<f64>,
    /// C++: `std::vector<uint8_t> rowintegral`
    pub rowintegral: Vec<u8>,
    /// C++: `std::unordered_multimap<uint64_t, HighsInt> hashToCutMap`
    pub hash_to_cut_map: HighsMultiBuckets,
    /// C++ `std::vector<HighsDomain::CutpoolPropagation*> propagationDomains`: replaced by the pool event log (`PoolEventLog`).
    pub events: PoolEventLog,
    /// C++: `std::set<std::pair<HighsInt, HighsInt>> propRows`
    pub prop_rows: BTreeSet<(HighsInt, HighsInt)>,
    /// C++: `double bestObservedScore`
    pub best_observed_score: f64,
    /// C++: `double minScoreFactor`
    pub min_score_factor: f64,
    /// C++: `double minDensityLim`
    pub min_density_lim: f64,
    /// C++: `HighsInt agelim_`
    pub agelim_: HighsInt,
    /// C++: `HighsInt softlimit_`
    pub softlimit_: HighsInt,
    /// C++: `HighsInt numLpCuts`
    pub num_lp_cuts: HighsInt,
    /// C++: `HighsInt numPropNzs`
    pub num_prop_nzs: HighsInt,
    /// C++: `HighsInt numPropRows`
    pub num_prop_rows: HighsInt,
    /// C++: `std::vector<HighsInt> ageDistribution`
    pub age_distribution: Vec<HighsInt>,
    /// C++: `std::vector<std::pair<HighsInt, double>> sortBuffer`
    pub sort_buffer: Vec<(HighsInt, f64)>,
    /// C++: `HighsInt index_`
    pub index_: HighsInt,
}

impl Clone for HighsCutPool {
    fn clone(&self) -> Self {
        HighsCutPool {
            matrix_: self.matrix_.clone(),
            rhs_: self.rhs_.clone(),
            ages_: self.ages_.clone(),
            num_lps_: self
                .num_lps_
                .iter()
                .map(|a| AtomicI16::new(a.load(Ordering::Relaxed)))
                .collect(),
            age_reset_while_locked_: self
                .age_reset_while_locked_
                .iter()
                .map(|a| AtomicU8::new(a.load(Ordering::Relaxed)))
                .collect(),
            has_synced_: self.has_synced_.clone(),
            rownormalization_: self.rownormalization_.clone(),
            maxabscoef_: self.maxabscoef_.clone(),
            rowintegral: self.rowintegral.clone(),
            hash_to_cut_map: self.hash_to_cut_map.clone(),
            events: self.events.clone(),
            prop_rows: self.prop_rows.clone(),
            best_observed_score: self.best_observed_score,
            min_score_factor: self.min_score_factor,
            min_density_lim: self.min_density_lim,
            agelim_: self.agelim_,
            softlimit_: self.softlimit_,
            num_lp_cuts: self.num_lp_cuts,
            num_prop_nzs: self.num_prop_nzs,
            num_prop_rows: self.num_prop_rows,
            age_distribution: self.age_distribution.clone(),
            sort_buffer: self.sort_buffer.clone(),
            index_: self.index_,
        }
    }
}

/// HiGHS: HighsCutPool.cpp::compute_cut_hash (lines 21-31)
/// C++: `static uint64_t compute_cut_hash(const HighsInt * Rindex, const double * Rvalue, double maxabscoef, const HighsInt Rlen)`
///
/// PORT NOTE: `Rindex`/`Rvalue`/`Rlen` become slices; `rlen` is the length of `rindex`
/// and `rvalue` (the C++ callers always pass the same length).
pub fn compute_cut_hash(
    rindex: &[HighsInt],
    rvalue: &[f64],
    maxabscoef: f64,
    rlen: HighsInt,
) -> u64 {
    let rlen = rlen as usize;

    // The C++ fills a `std::vector<uint32_t>` of value hash codes first. Here the codes
    // are streamed straight into the hash, which reads them in the same order.
    let scale = 1.0 / maxabscoef;
    let value_hash = vector_hash_words(
        rlen,
        rvalue[..rlen].iter().map(|&v| double_hash_code(scale * v)),
    );

    vector_hash(&rindex[..rlen]) ^ (value_hash >> 32)
}

impl HighsCutPool {
    /// HiGHS: HighsCutPool.cpp::HighsCutPool::isDuplicate (lines 47-87)
    /// C++: `bool isDuplicate(size_t hash, double norm, const HighsInt * Rindex, const double * Rvalue, HighsInt Rlen, double rhs)`
    pub fn is_duplicate(
        &mut self,
        hash: usize,
        norm: f64,
        rindex: &[HighsInt],
        rvalue: &[f64],
        rlen: HighsInt,
        rhs: f64,
    ) -> bool {
        let key = hash as u64;
        let count = self.hash_to_cut_map.count(key);
        let arvalue = self.matrix_.get_a_rvalue();
        let arindex = self.matrix_.get_a_rindex();

        for it in 0..count {
            let rowindex = self.hash_to_cut_map.get(key, it);
            let start = self.matrix_.get_row_start(rowindex);
            let end = self.matrix_.get_row_end(rowindex);

            if end - start != rlen {
                continue;
            }
            let (s, e) = (start as usize, end as usize);
            if rindex[..rlen as usize] == arindex[s..e] {
                let mut dotprod = 0.0;

                for (&a, &b) in rvalue[..rlen as usize].iter().zip(&arvalue[s..e]) {
                    dotprod += a * b;
                }

                let parallelism = dotprod * self.rownormalization_[rowindex as usize] * norm;

                if parallelism >= 1.0 - 1e-6 {
                    return true;
                }
            }
        }

        false
    }

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::getParallelism (lines 89-116)
    /// C++: `double getParallelism(HighsInt row1, HighsInt row2) const`
    pub fn get_parallelism(&self, row1: HighsInt, row2: HighsInt) -> f64 {
        let arindex = self.matrix_.get_a_rindex();
        let arvalue = self.matrix_.get_a_rvalue();

        // the two rows as slices, so the merge below indexes them relative to their start
        let (s1, e1) = (
            self.matrix_.get_row_start(row1) as usize,
            self.matrix_.get_row_end(row1) as usize,
        );
        let (s2, e2) = (
            self.matrix_.get_row_start(row2) as usize,
            self.matrix_.get_row_end(row2) as usize,
        );
        let (index1, value1) = (&arindex[s1..e1], &arvalue[s1..e1]);
        let (index2, value2) = (&arindex[s2..e2], &arvalue[s2..e2]);

        let mut dotprod = 0.0;
        let (mut i1, mut i2) = (0usize, 0usize);
        while i1 < index1.len() && i2 < index2.len() {
            let col1 = index1[i1];
            let col2 = index2[i2];

            if col1 < col2 {
                i1 += 1;
            } else if col2 < col1 {
                i2 += 1;
            } else {
                dotprod += value1[i1] * value2[i2];
                i1 += 1;
                i2 += 1;
            }
        }

        dotprod * self.rownormalization_[row1 as usize] * self.rownormalization_[row2 as usize]
    }

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::getParallelism (lines 118-148)
    /// C++: `double getParallelism(HighsInt row1, HighsInt row2, const HighsCutPool & pool2) const`
    /// PORT NOTE: signature names the overload `get_parallelism_2` (the Rust name for the 3-argument overload).
    pub fn get_parallelism_2(&self, row1: HighsInt, row2: HighsInt, pool2: &HighsCutPool) -> f64 {
        let mut i1 = self.matrix_.get_row_start(row1);
        let end1 = self.matrix_.get_row_end(row1);

        let mut i2 = pool2.matrix_.get_row_start(row2);
        let end2 = pool2.matrix_.get_row_end(row2);

        let arindex1 = self.matrix_.get_a_rindex();
        let arvalue1 = self.matrix_.get_a_rvalue();
        let arindex2 = pool2.matrix_.get_a_rindex();
        let arvalue2 = pool2.matrix_.get_a_rvalue();

        let mut dotprod = 0.0;
        while i1 != end1 && i2 != end2 {
            let col1 = arindex1[i1 as usize];
            let col2 = arindex2[i2 as usize];

            if col1 < col2 {
                i1 += 1;
            } else if col2 < col1 {
                i2 += 1;
            } else {
                dotprod += arvalue1[i1 as usize] * arvalue2[i2 as usize];
                i1 += 1;
                i2 += 1;
            }
        }

        dotprod * self.rownormalization_[row1 as usize] * pool2.rownormalization_[row2 as usize]
    }

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::lpCutRemoved (lines 150-160)
    /// C++: `void lpCutRemoved(HighsInt cut, bool thread_safe = false)`
    /// C++ default arguments: `thread_safe = false`
    pub fn lp_cut_removed(&mut self, cut: HighsInt, thread_safe: bool) {
        let cu = cut as usize;
        // C++ `numLps_[cut].fetch_add(-1, relaxed)` returns the old value.
        let n = self.num_lps_[cu].fetch_add(-1, Ordering::Relaxed);
        if thread_safe || (n as HighsInt) > 1 {
            return;
        }
        if self.matrix_.columns_linked(cut) {
            self.prop_rows.remove(&(self.ages_[cu] as HighsInt, cut));
            self.prop_rows.insert((1, cut));
        }
        self.ages_[cu] = 1;
        self.num_lp_cuts -= 1;
        self.age_distribution[1] += 1;
    }

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::performAging (lines 162-227)
    /// C++: `void performAging()`
    pub fn perform_aging(&mut self) {
        let cut_index_end = self.matrix_.get_num_rows();

        let mut agelim = self.agelim_;
        let mut num_available_cuts = self.get_num_available_cuts();
        while agelim > 5 && num_available_cuts > self.softlimit_ {
            num_available_cuts -= self.age_distribution[agelim as usize];
            agelim -= 1;
        }

        for i in 0..cut_index_end {
            let iu = i as usize;
            // Catch buffered changes (should only occur in parallel case)
            // C++ 174: `numLps_[i] > 0` is the implicit atomic load (seq_cst).
            if self.num_lps_[iu].load(Ordering::SeqCst) > 0 && self.ages_[iu] >= 0 {
                // Cut has been added to the LP, but age changes haven't been made
                self.age_distribution[self.ages_[iu] as usize] -= 1;
                if self.matrix_.columns_linked(i) {
                    self.prop_rows.remove(&(self.ages_[iu] as HighsInt, i));
                    self.prop_rows.insert((-1, i));
                }
                self.ages_[iu] = -1;
                self.num_lp_cuts += 1;
                self.age_reset_while_locked_[iu].store(0, Ordering::Relaxed);
            } else if self.num_lps_[iu].load(Ordering::SeqCst) == 0
                && self.ages_[iu] == -1
                && self.rhs_[iu] != kHighsInf
            {
                // Cut was removed from the LP, but age changes haven't been made
                if self.matrix_.columns_linked(i) {
                    self.prop_rows.remove(&(self.ages_[iu] as HighsInt, i));
                    self.prop_rows.insert((1, i));
                }
                self.ages_[iu] = 1;
                self.num_lp_cuts -= 1;
                self.age_distribution[1] += 1;
                self.age_reset_while_locked_[iu].store(0, Ordering::Relaxed);
                continue;
            } else if self.age_reset_while_locked_[iu].load(Ordering::Relaxed) == 1 {
                self.reset_age(i, false);
            }
            self.age_reset_while_locked_[iu].store(0, Ordering::Relaxed);
            if self.ages_[iu] < 0 {
                continue;
            }

            let is_propagated = self.matrix_.columns_linked(i);
            if is_propagated {
                self.prop_rows.remove(&(self.ages_[iu] as HighsInt, i));
            }
            self.age_distribution[self.ages_[iu] as usize] -= 1;
            self.ages_[iu] += 1;

            if (self.ages_[iu] as HighsInt) > agelim {
                // C++: for (p : propagationDomains) p->cutDeleted(i);
                // PORT NOTE: the event log has no observers until a domain registers, and
                // without one the notification is a no-op, so the push is skipped.
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
                self.has_synced_[iu] = false;
            } else {
                if is_propagated {
                    self.prop_rows.insert((self.ages_[iu] as HighsInt, i));
                }
                self.age_distribution[self.ages_[iu] as usize] += 1;
            }
        }

        debug_assert!(self.prop_rows.len() as HighsInt == self.num_prop_rows);
    }

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::separate (lines 229-460)
    /// C++: `void separate(const std::vector<double> & sol, const HighsDomain & domain, HighsCutSet & cutset, double feastol, const std::deque<HighsCutPool> & cutpools, bool thread_safe = false)`
    /// C++ default arguments: `thread_safe = false`
    /// PORT NOTE: signature changed because `self` is one of `ms.mipdata_.cutpools` (the worker
    /// pool or `getCutPool()` = `cutpools[0]`) and C++ also reads the whole deque. `cutpools` is
    /// `&[Slot<HighsCutPool>]`. The caller detaches its own pool and passes the deque:
    /// `detach!(ms.mipdata_.cutpools[k], |pool| pool.separate(sol, &dom, cutset, feastol, &ms.mipdata_.cutpools, ts))`.
    /// The body reads `cutpools[j]` only for `j != self.index_` (HighsCutPool.cpp:403-413), so the
    /// detached slot is never dereferenced. The pool needs no `ms`, so nothing can re-enter it.
    pub fn separate(
        &mut self,
        sol: &[f64],
        domain: &HighsDomain,
        cutset: &mut HighsCutSet,
        feastol: f64,
        cutpools: &[Slot<HighsCutPool>],
        thread_safe: bool,
    ) {
        self.separate_hooked(None, sol, domain, cutset, feastol, cutpools, thread_safe)
    }

    /// `separate` with the rusolve cut-score hook (mip/hooks.rs); `hooks == None` is exactly `separate`.
    #[allow(clippy::too_many_arguments)]
    pub fn separate_hooked(
        &mut self,
        hooks: Option<&crate::mip::hooks::MipHooksHandle>,
        sol: &[f64],
        domain: &HighsDomain,
        cutset: &mut HighsCutSet,
        feastol: f64,
        cutpools: &[Slot<HighsCutPool>],
        thread_safe: bool,
    ) {
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
            // rusolve cut-score hook (mip/hooks.rs)
            let score = match hooks {
                Some(h) => h.cut_score(&crate::mip::hooks::CutScoreInfo {
                    pool_row: i,
                    violation: viol,
                    num_nonzeros: end - start,
                    num_active_nonzeros: num_active_nzs,
                    active_norm: f64::from(rownorm).sqrt(),
                    default_score: score,
                }),
                None => score,
            };

            efficacious_cuts.push((score, i));
        }
        debug_assert!(self.prop_rows.len() as HighsInt == self.num_prop_rows);
        if efficacious_cuts.is_empty() {
            return;
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

        let orignumcuts = cutset.num_cuts();
        let origselectednnz = cutset.a_rindex_.len() as HighsInt;
        let mut selectednnz = origselectednnz;

        for &(_score, p) in efficacious_cuts.iter() {
            let mut discard = false;
            let maxpar = 0.1;
            for (i, &cut_pool_idx) in cutset.cutpools.iter().enumerate() {
                if cut_pool_idx == self.index_ {
                    if self.get_parallelism(cutset.cutindices[i], p) > maxpar {
                        discard = true;
                        break;
                    }
                } else {
                    // Warning: This assumes the cuts in the pool are not changing during,
                    // this query, i.e., the worker's pool and the global pool.
                    // Currently safe, but doesn't generalise to all designs.
                    if self.get_parallelism_2(
                        p,
                        cutset.cutindices[i],
                        &cutpools[cut_pool_idx as usize],
                    ) > maxpar
                    {
                        discard = true;
                        break;
                    }
                }
            }

            if discard {
                continue;
            }

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

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::separateLpCutsAfterRestart (lines 462-501)
    /// C++: `void separateLpCutsAfterRestart(HighsCutSet & cutset)`
    pub fn separate_lp_cuts_after_restart(&mut self, cutset: &mut HighsCutSet) {
        // should only be called after a restart with a fresh row matrix right now
        debug_assert!(self.matrix_.get_num_del_rows() == 0);
        let numcuts = self.matrix_.get_num_rows();

        cutset.cutindices.resize(numcuts as usize, 0);
        cutset.cutpools.resize(numcuts as usize, self.index_);
        // C++ std::iota(cutindices.begin(), cutindices.end(), 0)
        for (k, v) in cutset.cutindices.iter_mut().enumerate() {
            *v = k as HighsInt;
        }
        cutset.resize(self.matrix_.nonzero_capacity() as HighsInt);

        let mut offset: HighsInt = 0;
        for i in 0..cutset.num_cuts() {
            let iu = i as usize;
            self.age_distribution[self.ages_[iu] as usize] -= 1;
            self.num_lp_cuts += 1;
            if self.matrix_.columns_linked(i) {
                self.prop_rows.remove(&(self.ages_[iu] as HighsInt, i));
                self.prop_rows.insert((-1, i));
            }
            // C++ 482: `numLps_[i] = 1` is the implicit atomic store (seq_cst).
            self.num_lps_[iu].store(1, Ordering::SeqCst);
            self.ages_[iu] = -1;
            cutset.a_rstart_[iu] = offset;
            let cut = cutset.cutindices[iu];
            let start = self.matrix_.get_row_start(cut);
            let end = self.matrix_.get_row_end(cut);
            cutset.upper_[iu] = self.rhs_[cut as usize];

            for j in start..end {
                debug_assert!(offset < self.matrix_.nonzero_capacity() as HighsInt);
                cutset.a_rvalue_[offset as usize] = self.matrix_.get_a_rvalue()[j as usize];
                cutset.a_rindex_[offset as usize] = self.matrix_.get_a_rindex()[j as usize];
                offset += 1;
            }
        }

        let ncuts = cutset.num_cuts() as usize;
        cutset.a_rstart_[ncuts] = offset;

        debug_assert!(self.prop_rows.len() as HighsInt == self.num_prop_rows);
    }

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::addCut (lines 503-638)
    /// C++: `HighsInt addCut(const HighsMipSolver & mipsolver, HighsInt * Rindex, double * Rvalue, HighsInt Rlen, double rhs, bool integral = false, bool propagate = true, bool extractCliques = true, bool isConflict = false)`
    /// C++ default arguments: `integral = false`, `propagate = true`, `extract_cliques = true`, `is_conflict = false`
    ///
    /// PORT NOTE: signature changed from `mipsolver: &HighsMipSolver` to `ms: &mut HighsMipSolver`
    /// because the body needs mutable access to the global pool (`isDuplicate` is non-const) and to
    /// `mipdata_->cliquetable` (`extractCliquesFromCut` is non-const). Both are reached through `ms`.
    /// `this == &getCutPool()` is written `self.index_ == 0` (the global pool is `cutpools[0]`
    /// with index 0, and every other pool has index >= 1).
    pub fn add_cut(
        &mut self,
        ms: &mut HighsMipSolver,
        rindex: &mut [HighsInt],
        rvalue: &mut [f64],
        rlen: HighsInt,
        rhs: f64,
        integral: bool,
        propagate: bool,
        extract_cliques: bool,
        is_conflict: bool,
    ) -> HighsInt {
        // C++: mipsolver.mipdata_->debugSolution.checkCut(Rindex, Rvalue, Rlen, rhs);
        ms.mipdata_
            .debug_solution
            .check_cut(rindex, rvalue, rlen, rhs);

        self.sort_buffer.resize(rlen as usize, (0, 0.0));

        // compute 1/||a|| for the cut
        // as it is only computed once
        let mut norm = 0.0;
        let mut maxabscoef = 0.0;
        for i in 0..rlen as usize {
            norm += rvalue[i] * rvalue[i];
            maxabscoef = cpp_max(maxabscoef, rvalue[i].abs());
            self.sort_buffer[i].0 = rindex[i];
            self.sort_buffer[i].1 = rvalue[i];
        }
        // PORT NOTE: C++ `pdqsort_branchless`; the keys are distinct column indices, so the
        // sorted order is unique and the branchy `pdqsort_by` gives the same result.
        pdqsort_by(&mut self.sort_buffer, |a, b| a.0 < b.0);
        for i in 0..rlen as usize {
            rindex[i] = self.sort_buffer[i].0;
            rvalue[i] = self.sort_buffer[i].1;
        }
        let h = compute_cut_hash(rindex, rvalue, maxabscoef, rlen);
        let normalization = 1.0 / norm.sqrt();

        // Warning: This global duplicate check assumes the global pool doesn't
        // have cuts added or deleted during time when local pools can add a cut.
        // PORT NOTE: C++ `this != &getCutPool()` is `self.index_ != 0`; the global pool is cutpools[0].
        if self.index_ != 0
            && ms
                .get_cut_pool()
                .is_duplicate(h as usize, normalization, rindex, rvalue, rlen, rhs)
        {
            return -1;
        }
        if self.is_duplicate(h as usize, normalization, rindex, rvalue, rlen, rhs) {
            return -1;
        }

        // if (Rlen > 0.15 * matrix_.numCols())
        //   printf("cut with len %d not propagated\n", Rlen);
        let mut propagate = propagate;
        if propagate {
            let new_prop_nzs = self.num_prop_nzs + rlen;

            let avg_model_nzs = ms.num_nonzero() as f64 / ms.num_row() as f64;

            let new_avg_prop_nzs = new_prop_nzs as f64 / (self.num_prop_rows + 1) as f64;

            let alpha = 2.0;
            if is_conflict {
                // for conflicts we allow an increased average propagation density
                if new_avg_prop_nzs > cpp_max(alpha * avg_model_nzs, self.min_density_lim) {
                    propagate = false;
                } else {
                    self.num_prop_rows += 1;
                    self.num_prop_nzs = new_prop_nzs;
                }
            } else {
                // for cuts we do not want to accept any dense cuts and don't use the
                // average but its actual length
                if (rlen as f64) >= cpp_max(alpha * avg_model_nzs, self.min_density_lim) {
                    propagate = false;
                } else {
                    self.num_prop_rows += 1;
                    self.num_prop_nzs = new_prop_nzs;
                }
            }
        }

        // if we have more than twice the number of nonzeros of the model in use for
        // propagation we stop propagating the rows with the highest age
        let mut prop_row_excess_nzs = self.num_prop_nzs - 2 * ms.num_nonzero();
        if prop_row_excess_nzs > 0 {
            // C++ `rbegin()` walk: the visited elements are the last `visited` of the set.
            let rev_rows: Vec<(HighsInt, HighsInt)> =
                self.prop_rows.iter().rev().copied().collect();
            let mut visited = 0usize;

            while prop_row_excess_nzs > 0 && visited < rev_rows.len() {
                let len = self.get_row_length(rev_rows[visited].1);
                prop_row_excess_nzs -= len;
                self.num_prop_nzs -= len;
                self.num_prop_rows -= 1;
                visited += 1;
            }

            for &(_, row) in &rev_rows[..visited] {
                self.matrix_.unlink_columns(row);
                // C++: for (p : propagationDomains) p->cutDeleted(row, true);
                // PORT NOTE: see perform_aging (no observers -> no-op push skipped).
                if !self.events.observers.is_empty() {
                    self.events.push(PoolEvent::CutDeleted {
                        cut: row,
                        only_for_prop: true,
                    });
                }
            }

            for key in &rev_rows[..visited] {
                self.prop_rows.remove(key);
            }
        }

        // if no such cut exists we append the new cut
        let rowindex = self.matrix_.add_row(rindex, rvalue, rlen, propagate);
        let count = self.hash_to_cut_map.count(h);
        // C++ `emplace` into an unordered_multimap: the new element goes at the end of its equal range.
        self.hash_to_cut_map.emplace_hint(Some(count), h, rowindex);

        let ri = rowindex as usize;
        if rowindex == self.rhs_.len() as HighsInt {
            self.rhs_.resize(ri + 1, 0.0);
            self.ages_.resize(ri + 1, 0);
            self.num_lps_.resize_with(ri + 1, || AtomicI16::new(0));
            self.age_reset_while_locked_
                .resize_with(ri + 1, || AtomicU8::new(0));
            self.has_synced_.resize(ri + 1, false);
            self.rownormalization_.resize(ri + 1, 0.0);
            self.maxabscoef_.resize(ri + 1, 0.0);
            self.rowintegral.resize(ri + 1, 0);
        }

        // set the right hand side and reset the age
        self.rhs_[ri] = rhs;
        self.ages_[ri] = cpp_max(0, self.agelim_ - 5) as i16;
        self.age_distribution[self.ages_[ri] as usize] += 1;
        self.rowintegral[ri] = integral as u8;
        // C++ 617: `numLps_[rowindex] = 0` is the implicit atomic store (seq_cst).
        self.num_lps_[ri].store(0, Ordering::SeqCst);
        self.age_reset_while_locked_[ri].store(0, Ordering::Relaxed);
        self.has_synced_[ri] = false;
        if propagate {
            self.prop_rows
                .insert((self.ages_[ri] as HighsInt, rowindex));
        }
        debug_assert!(self.prop_rows.len() as HighsInt == self.num_prop_rows);

        self.rownormalization_[ri] = normalization;
        self.maxabscoef_[ri] = maxabscoef;

        // C++: for (p : propagationDomains) p->cutAdded(rowindex, propagate);
        // PORT NOTE: the row snapshot (M3) is taken here, after the row is in the matrix.
        if !self.events.observers.is_empty() {
            let start = self.matrix_.get_row_start(rowindex) as usize;
            let end = self.matrix_.get_row_end(rowindex) as usize;
            self.events.push(PoolEvent::CutAdded {
                cut: rowindex,
                propagate,
                inds: self.matrix_.get_a_rindex()[start..end].to_vec(),
                vals: self.matrix_.get_a_rvalue()[start..end].to_vec(),
                rhs: self.rhs_[ri],
            });
        }

        // PORT NOTE: C++ `this == &getCutPool()` is `self.index_ == 0` (see the note above).
        if extract_cliques && self.index_ == 0 {
            // if this is the global cutpool extract cliques from the cut
            if rlen <= 100 {
                // PORT NOTE: detach the clique table so it can take `ms` (M2), as C++
                // `mipdata_->cliquetable.extractCliquesFromCut(mipsolver, ...)`. Its bound
                // changes replay this pool's events, so `self` is lent back to its slot.
                let pool = self.index_ as usize;
                lend!(
                    ms.mipdata_.cutpools[pool],
                    self,
                    detach!(ms.mipdata_.cliquetable, |ct| ct
                        .extract_cliques_from_cut(ms, rindex, rvalue, rlen, rhs))
                );
            }
        }

        rowindex
    }

    /// HiGHS: HighsCutPool.cpp::HighsCutPool::syncCutPool (lines 640-665)
    /// C++: `void syncCutPool(const HighsMipSolver & mipsolver, HighsCutPool & syncpool)`
    /// PORT NOTE: signature changed from `mipsolver: &HighsMipSolver` to `ms: &mut HighsMipSolver`
    /// because it forwards to `add_cut`, which needs `ms` mutably.
    pub fn sync_cut_pool(&mut self, ms: &mut HighsMipSolver, syncpool: &mut HighsCutPool) {
        let cut_index_end = self.matrix_.get_num_rows();

        for i in 0..cut_index_end {
            let iu = i as usize;
            // Only sync cuts in the LP that are not already synced
            // C++ 648: `numLps_[i] > 0` is the implicit atomic load (seq_cst).
            if (self.num_lps_[iu].load(Ordering::SeqCst) > 0
                || self.age_reset_while_locked_[iu].load(Ordering::Relaxed) == 1)
                && !self.has_synced_[iu]
            {
                let mut rlen: HighsInt = 0;
                let mut rindex: &[HighsInt] = &[];
                let mut rvalue: &[f64] = &[];
                self.get_cut(i, &mut rlen, &mut rindex, &mut rvalue);
                // copy cut into something mutable (addCut reorders so can't take const)
                let mut idxs: Vec<HighsInt> = rindex.to_vec();
                let mut vals: Vec<f64> = rvalue.to_vec();
                syncpool.add_cut(
                    ms,
                    &mut idxs,
                    &mut vals,
                    rlen,
                    self.rhs_[iu],
                    self.rowintegral[iu] != 0,
                    true,
                    true,
                    false,
                );
                self.has_synced_[iu] = true;
            }
        }

        debug_assert!(self.prop_rows.len() as HighsInt == self.num_prop_rows);
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::HighsCutPool (lines 85-98)
    /// C++: `HighsCutPool(HighsInt ncols, HighsInt agelim, HighsInt softlimit, HighsInt index)`
    /// The 1.15.1 constructor: `new_with_age_flag` with `h16_cutpool_age = false`.
    pub fn new(ncols: HighsInt, agelim: HighsInt, softlimit: HighsInt, index: HighsInt) -> Self {
        Self::new_with_age_flag(ncols, agelim, softlimit, index, false)
    }

    /// HiGHS 1.16 (#3574): HighsCutPool.h::HighsCutPool::HighsCutPool (lines 85-98)
    /// C++: `HighsCutPool(HighsInt ncols, HighsInt agelim, HighsInt softlimit, HighsInt index)`
    /// `h16_cutpool_age` is `options_mip_.h16_cutpool_age`. Flag off reproduces 1.15.1; on sizes
    /// `ageDistribution` as `std::max(agelim, 1) + 1`, so `agelim == 0` no longer undersizes it.
    /// Extra argument (not in the C++ signature) carries the HiGHS 1.16 flag. The 4-argument `new`
    /// stays for the existing callers and is the flag-off path.
    pub fn new_with_age_flag(
        ncols: HighsInt,
        agelim: HighsInt,
        softlimit: HighsInt,
        index: HighsInt,
        h16_cutpool_age: bool,
    ) -> Self {
        let mut pool = HighsCutPool {
            matrix_: HighsDynamicRowMatrix::new(ncols),
            agelim_: agelim,
            softlimit_: softlimit,
            num_lp_cuts: 0,
            num_prop_nzs: 0,
            num_prop_rows: 0,
            index_: index,
            ..Default::default()
        };
        // HiGHS 1.16 #3574 (h16_cutpool_age): HighsCutPool.h lines 94 (1.16.0)
        if h16_cutpool_age {
            pool.age_distribution
                .resize((std::cmp::max(agelim, 1) + 1) as usize, 0);
        } else {
            pool.age_distribution.resize((agelim + 1) as usize, 0);
        }
        pool.min_score_factor = 0.9;
        pool.best_observed_score = 0.0;
        pool.min_density_lim = 0.1 * ncols as f64;
        pool
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::getMatrix (lines 99-99)
    /// C++: `const HighsDynamicRowMatrix & getMatrix() const`
    pub fn get_matrix(&self) -> &HighsDynamicRowMatrix {
        &self.matrix_
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::getRhs (lines 101-101)
    /// C++: `const std::vector<double> & getRhs() const`
    pub fn get_rhs(&self) -> &[f64] {
        &self.rhs_
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::resetAge (lines 106-121)
    /// C++: `void resetAge(HighsInt cut, bool thread_safe = false)`
    /// C++ default arguments: `thread_safe = false`
    pub fn reset_age(&mut self, cut: HighsInt, thread_safe: bool) {
        let cu = cut as usize;
        if self.ages_[cu] > 0 {
            if thread_safe {
                self.age_reset_while_locked_[cu].store(1, Ordering::Relaxed);
                return;
            }
            if self.matrix_.columns_linked(cut) {
                self.prop_rows.remove(&(self.ages_[cu] as HighsInt, cut));
                self.prop_rows.insert((0, cut));
            }
            self.age_distribution[self.ages_[cu] as usize] -= 1;
            self.age_distribution[0] += 1;
            self.ages_[cu] = 0;
            self.age_reset_while_locked_[cu].store(0, Ordering::Relaxed);
        }
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::addPropagationDomain (lines 132-134)
    /// C++: `void addPropagationDomain(HighsDomain::CutpoolPropagation * domain)`
    /// PORT NOTE: signature changed from `Option<&mut HighsDomain_CutpoolPropagation>` to
    /// `key: DomainKey`. The observer is identified by its M3 key (the coordinator's ruling),
    /// because the C++ `this` pointer has no Rust equivalent. A C++ nullptr is never registered,
    /// so no `Option` is needed. The cursor is the log end, as a fresh observer starts there.
    pub fn add_propagation_domain(&mut self, key: DomainKey) {
        let end = self.events.end();
        self.events.register(key, end);
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::removePropagationDomain (lines 136-143)
    /// C++: `void removePropagationDomain(HighsDomain::CutpoolPropagation * domain)`
    /// PORT NOTE: signature changed to `key: DomainKey` because the observer is identified by its M3 key (see `add_propagation_domain`).
    /// Keys are unique, so unregistering by key matches C++ removing the matching entry.
    pub fn remove_propagation_domain(&mut self, key: DomainKey) {
        self.events.unregister(key);
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::setAgeLimit (lines 145-148)
    /// C++: `void setAgeLimit(HighsInt agelim)`
    /// `h16_cutpool_age` (extra argument): see `new_with_age_flag`. Flag off is the 1.15.1 body.
    pub fn set_age_limit(&mut self, agelim: HighsInt, h16_cutpool_age: bool) {
        self.agelim_ = agelim;
        // HiGHS 1.16 #3574 (h16_cutpool_age): HighsCutPool.h line 147 (1.16.0)
        if h16_cutpool_age {
            self.age_distribution
                .resize((std::cmp::max(agelim, 1) + 1) as usize, 0);
        } else {
            self.age_distribution.resize((agelim + 1) as usize, 0);
        }
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::increaseNumLps (lines 150-153)
    /// C++: `void increaseNumLps(HighsInt cut, HighsInt n)`
    pub fn increase_num_lps(&mut self, cut: HighsInt, n: HighsInt) {
        let cu = cut as usize;
        // C++ HighsCutPool.h:151 `assert(numLps_[cut] >= 1)`: an implicit (seq_cst) load.
        debug_assert!(self.num_lps_[cu].load(Ordering::SeqCst) >= 1);
        // C++ `fetch_add(n)` on an atomic<int16_t>: the sum wraps in 16 bits.
        self.num_lps_[cu].fetch_add(n as i16, Ordering::Relaxed);
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::cutIsIntegral (lines 162-162)
    /// C++: `bool cutIsIntegral(HighsInt cut) const`
    pub fn cut_is_integral(&self, cut: HighsInt) -> bool {
        self.rowintegral[cut as usize] != 0
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::getNumCuts (lines 164-166)
    /// C++: `HighsInt getNumCuts() const`
    pub fn get_num_cuts(&self) -> HighsInt {
        self.matrix_.get_num_rows() - self.matrix_.get_num_del_rows()
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::getNumAvailableCuts (lines 168-168)
    /// C++: `HighsInt getNumAvailableCuts() const`
    pub fn get_num_available_cuts(&self) -> HighsInt {
        self.get_num_cuts() - self.num_lp_cuts
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::getMaxAbsCutCoef (lines 170-170)
    /// C++: `double getMaxAbsCutCoef(HighsInt cut) const`
    pub fn get_max_abs_cut_coef(&self, cut: HighsInt) -> f64 {
        self.maxabscoef_[cut as usize]
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::getRowNormalization (lines 172-174)
    /// C++: `double getRowNormalization(HighsInt cut) const`
    pub fn get_row_normalization(&self, cut: HighsInt) -> f64 {
        self.rownormalization_[cut as usize]
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::getRowLength (lines 181-183)
    /// C++: `HighsInt getRowLength(HighsInt row) const`
    pub fn get_row_length(&self, row: HighsInt) -> HighsInt {
        self.matrix_.get_row_end(row) - self.matrix_.get_row_start(row)
    }

    /// HiGHS: HighsCutPool.h::HighsCutPool::getCut (lines 185-191)
    /// C++: `void getCut(HighsInt cut, HighsInt & cutlen, const HighsInt *& cutinds, const double *& cutvals) const`
    pub fn get_cut<'s>(
        &'s self,
        cut: HighsInt,
        cutlen: &mut HighsInt,
        cutinds: &mut &'s [HighsInt],
        cutvals: &mut &'s [f64],
    ) {
        let start = self.matrix_.get_row_start(cut);
        *cutlen = self.matrix_.get_row_end(cut) - start;
        // PORT NOTE: a deleted row has range (-1, -1), so C++ yields `cutlen == 0` and
        // pointers one before the arrays that are never dereferenced. Conflict analysis
        // can reach this when a reason cut was aged out of a worker's pool while still on
        // the local domain change stack; the port returns empty slices.
        if start < 0 {
            debug_assert_eq!(*cutlen, 0);
            *cutinds = &[];
            *cutvals = &[];
            return;
        }
        let (s, e) = (start as usize, (start + *cutlen) as usize);
        *cutinds = &self.matrix_.get_a_rindex()[s..e];
        *cutvals = &self.matrix_.get_a_rvalue()[s..e];
    }
}

impl HighsCutSet {
    /// HiGHS: HighsCutPool.h::HighsCutSet::numCuts (lines 31-31)
    /// C++: `HighsInt numCuts() const`
    pub fn num_cuts(&self) -> HighsInt {
        self.cutindices.len() as HighsInt
    }

    /// HiGHS: HighsCutPool.h::HighsCutSet::resize (lines 33-40)
    /// C++: `void resize(HighsInt nnz)`
    pub fn resize(&mut self, nnz: HighsInt) {
        let ncuts = self.num_cuts() as usize;
        self.lower_.resize(ncuts, -kHighsInf);
        self.upper_.resize(ncuts, 0.0);
        self.a_rstart_.resize(ncuts + 1, 0);
        self.a_rindex_.resize(nnz as usize, 0);
        self.a_rvalue_.resize(nnz as usize, 0.0);
    }

    /// HiGHS: HighsCutPool.h::HighsCutSet::clear (lines 42-49)
    /// C++: `void clear()`
    /// PORT NOTE: `lower_` is not cleared, as in C++.
    pub fn clear(&mut self) {
        self.cutindices.clear();
        self.cutpools.clear();
        self.upper_.clear();
        self.a_rstart_.clear();
        self.a_rindex_.clear();
        self.a_rvalue_.clear();
    }

    /// HiGHS: HighsCutPool.h::HighsCutSet::empty (lines 51-51)
    /// C++: `bool empty() const`
    pub fn empty(&self) -> bool {
        self.cutindices.is_empty()
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use super::*;

    #[test]
    fn w5_mip_atom_num_lps_and_age_reset_are_atomic() {
        // `numLps_` (`deque<atomic<int16_t>>`) and `ageResetWhileLocked_` (`deque<atomic<uint8_t>>`) are relaxed atomics (HighsCutPool.h:59-60); `lpCutRemoved` / `performAging` give the C++ values under `parallelLockActive()`.

        // Test 1: lpCutRemoved uses fetch_add(-1, Relaxed) and returns the old value.
        // C++: HighsCutPool.cpp:151 `const HighsInt n = numLps_[cut].fetch_add(-1, relaxed);`
        {
            let mut pool = HighsCutPool::new(10, 10, 100, 0);
            // lpCutRemoved reads columnsLinked(cut): the cut must be a row of the pool matrix.
            pool.matrix_.add_row(&mut [0], &mut [1.0], 1, false);
            // Set up: one cut in the LP with numLps = 3
            pool.num_lps_.resize_with(1, || AtomicI16::new(3));
            pool.ages_.resize(1, 0);
            pool.rhs_.resize(1, 1.0);
            pool.age_distribution.resize(11, 0);
            pool.num_lp_cuts = 1;

            // thread_safe = true: early return after fetch_add, no age changes
            pool.lp_cut_removed(0, true);
            assert_eq!(pool.num_lps_[0].load(Ordering::Relaxed), 2);
            assert_eq!(pool.ages_[0], 0); // unchanged because thread_safe

            // thread_safe = false, n = 2 > 1: early return after fetch_add
            pool.lp_cut_removed(0, false);
            assert_eq!(pool.num_lps_[0].load(Ordering::Relaxed), 1);
            assert_eq!(pool.ages_[0], 0); // unchanged because n > 1

            // thread_safe = false, n = 1: not > 1, so age is set to 1
            pool.lp_cut_removed(0, false);
            assert_eq!(pool.num_lps_[0].load(Ordering::Relaxed), 0);
            assert_eq!(pool.ages_[0], 1);
        }

        // Test 2: increaseNumLps uses fetch_add(n, Relaxed) with wrapping.
        // C++: HighsCutPool.h:152 `numLps_[cut].fetch_add(n, relaxed);`
        {
            let mut pool = HighsCutPool::new(10, 10, 100, 0);
            pool.num_lps_.resize_with(1, || AtomicI16::new(1));

            pool.increase_num_lps(0, 3);
            assert_eq!(pool.num_lps_[0].load(Ordering::Relaxed), 4);

            // Wrapping: i16::MAX + 1 wraps to i16::MIN
            pool.num_lps_[0].store(i16::MAX, Ordering::Relaxed);
            pool.increase_num_lps(0, 1);
            assert_eq!(pool.num_lps_[0].load(Ordering::Relaxed), i16::MIN);
        }

        // Test 3: resetAge with thread_safe=true stores 1 to ageResetWhileLocked_.
        // C++: HighsCutPool.h:109 `ageResetWhileLocked_[cut].store(1, relaxed);`
        {
            let mut pool = HighsCutPool::new(10, 10, 100, 0);
            pool.matrix_.add_row(&mut [0], &mut [1.0], 1, false);
            pool.ages_.resize(1, 5);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(0));
            pool.age_distribution.resize(11, 0);

            pool.reset_age(0, true);
            assert_eq!(pool.age_reset_while_locked_[0].load(Ordering::Relaxed), 1);
            assert_eq!(pool.ages_[0], 5); // age unchanged (early return)
        }

        // Test 4: resetAge with thread_safe=false resets age and stores 0.
        // C++: HighsCutPool.h:119 `ageResetWhileLocked_[cut].store(0, relaxed);`
        {
            let mut pool = HighsCutPool::new(10, 10, 100, 0);
            pool.matrix_.add_row(&mut [0], &mut [1.0], 1, false);
            pool.ages_.resize(1, 5);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(1));
            pool.age_distribution.resize(11, 0);
            pool.age_distribution[5] = 1;

            pool.reset_age(0, false);
            assert_eq!(pool.age_reset_while_locked_[0].load(Ordering::Relaxed), 0);
            assert_eq!(pool.ages_[0], 0);
            assert_eq!(pool.age_distribution[5], 0);
            assert_eq!(pool.age_distribution[0], 1);
        }

        // Test 5: performAging reads numLps_ (implicit seq_cst load) and ageResetWhileLocked_ (relaxed load).
        // C++: HighsCutPool.cpp:174 `if (numLps_[i] > 0 && ages_[i] >= 0)`
        // C++: HighsCutPool.cpp:195 `else if (ageResetWhileLocked_[i].load(relaxed) == 1)`
        {
            let mut pool = HighsCutPool::new(10, 10, 100, 0);
            pool.matrix_.add_row(&mut [0], &mut [1.0], 1, false);
            // Set up a cut that is in the LP (numLps > 0) with age >= 0
            pool.num_lps_.resize_with(1, || AtomicI16::new(1));
            pool.ages_.resize(1, 3);
            pool.rhs_.resize(1, 1.0);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(0));
            pool.age_distribution.resize(11, 0);
            pool.age_distribution[3] = 1;
            pool.num_lp_cuts = 0;

            // performAging should catch the buffered change: numLps > 0 && age >= 0
            pool.perform_aging();
            assert_eq!(pool.ages_[0], -1); // age set to -1 (in LP)
            assert_eq!(pool.age_reset_while_locked_[0].load(Ordering::Relaxed), 0);
            assert_eq!(pool.num_lp_cuts, 1);
        }

        // Test 6: separateLpCutsAfterRestart stores 1 to numLps_ (implicit seq_cst store).
        // C++: HighsCutPool.cpp:482 `numLps_[i] = 1;`
        {
            let mut pool = HighsCutPool::new(10, 10, 100, 0);
            // Add a row to the matrix so get_num_rows() returns 1
            let mut rindex = [0i32, 1];
            let mut rvalue = [1.0f64, 1.0];
            pool.matrix_.add_row(&mut rindex, &mut rvalue, 2, false);
            pool.num_lps_.resize_with(1, || AtomicI16::new(0));
            pool.ages_.resize(1, 2);
            pool.rhs_.resize(1, 5.0);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(0));
            pool.age_distribution.resize(11, 0);
            pool.age_distribution[2] = 1;
            pool.num_lp_cuts = 0;

            let mut cutset = HighsCutSet::default();
            pool.separate_lp_cuts_after_restart(&mut cutset);
            assert_eq!(pool.num_lps_[0].load(Ordering::Relaxed), 1);
            assert_eq!(pool.ages_[0], -1);
        }
    }

    // HiGHS 1.16 #3574 (h16_cutpool_age): C++ HighsCutPool.h lines 94 and 147 (1.16.0) size
    // ageDistribution as std::max(agelim, 1) + 1. The 1.15.1 size is agelim + 1.

    /// HiGHS 1.16 #3574: with agelim 0 the 1.16 constructor sizes ageDistribution to 2 (not 1).
    #[test]
    fn h16_cutpool_age_on_constructor_sizes_zero_limit_to_two() {
        let pool = HighsCutPool::new_with_age_flag(10, 0, 100, 0, true);
        assert_eq!(pool.age_distribution.len(), 2);
        let pool = HighsCutPool::new_with_age_flag(10, 5, 100, 0, true);
        assert_eq!(pool.age_distribution.len(), 6);
    }

    /// HiGHS 1.16 #3574: with the flag off, agelim 0 gives the 1.15.1 size of 1, and `new` is the same.
    #[test]
    fn h16_cutpool_age_off_constructor_keeps_1_15_1_size() {
        let pool = HighsCutPool::new_with_age_flag(10, 0, 100, 0, false);
        assert_eq!(pool.age_distribution.len(), 1);
        let pool = HighsCutPool::new(10, 0, 100, 0);
        assert_eq!(pool.age_distribution.len(), 1);
        let pool = HighsCutPool::new(10, 5, 100, 0);
        assert_eq!(pool.age_distribution.len(), 6);
    }

    /// HiGHS 1.16 #3574: set_age_limit with the flag on uses max(agelim, 1) + 1 as well.
    #[test]
    fn h16_cutpool_age_on_set_age_limit_sizes_zero_limit_to_two() {
        let mut pool = HighsCutPool::new(10, 10, 100, 0);
        pool.set_age_limit(0, true);
        assert_eq!(pool.age_distribution.len(), 2);
        assert_eq!(pool.agelim_, 0);
    }

    /// HiGHS 1.16 #3574: set_age_limit with the flag off keeps the 1.15.1 agelim + 1 size.
    #[test]
    fn h16_cutpool_age_off_set_age_limit_keeps_1_15_1_size() {
        let mut pool = HighsCutPool::new(10, 10, 100, 0);
        pool.set_age_limit(0, false);
        assert_eq!(pool.age_distribution.len(), 1);
        pool.set_age_limit(4, false);
        assert_eq!(pool.age_distribution.len(), 5);
    }
}
