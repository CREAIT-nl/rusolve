//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! solution_pool.h / solution_pool.cpp.
//!
//! The HighsSpinMutex is dropped: the path is single-threaded, so
//! every `{ scoped_lock }` block of the C++ is just a block here.

use std::rc::Rc;

use crate::mip::fpr_local_mip::rng::{S5Rng, s5_uniform_int_i32, s5_uniform_real};
use crate::mip::fpr_local_mip::{
    IntegerColumns, PoolAddResult, PoolEntry, PoolEntryRef, PoolSnapshot, RestartPolicy,
    SolutionPool, kDiversityMinHammingFrac, kDiversityObjTolerance, kImprovementObjMargin,
};
use crate::util::highs_utils::{cpp_max, cpp_min};

/// mip-heuristics (36f7d9d): solution_pool.cpp::<anonymous>::round_half_away (lines 39-47)
fn round_half_away(x: f64) -> f64 {
    const kTwoPow52: f64 = 4503599627370496.0;
    if !(x.abs() < kTwoPow52) {
        return x;
    }
    let t = (x as i64) as f64;
    let frac = x - t;
    t + (frac >= 0.5) as i32 as f64 - (frac <= -0.5) as i32 as f64
}

/// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::hamming_distance (lines 51-62)
fn hamming_distance(index: &[i32], a: &[f64], b: &[f64]) -> i32 {
    let n = cpp_min(a.len(), b.len());
    let mut dist: i32 = 0;
    for &j in index {
        if (j as usize) >= n {
            break; // ascending, so no later column is in both vectors
        }
        dist += (round_half_away(a[j as usize]) != round_half_away(b[j as usize])) as i32;
    }
    dist
}

/// mip-heuristics (36f7d9d): solution_pool.cpp::<anonymous>::far_enough (lines 69-72)
fn far_enough(dist: i32, num_integers: i32) -> bool {
    (dist as f64) / (num_integers as f64) >= kDiversityMinHammingFrac
}

/// mip-heuristics (36f7d9d): solution_pool.cpp::<anonymous>::diversity_victim (lines 94-112)
fn diversity_victim(dists: &[i32], num_integers: i32) -> i32 {
    let mut min_dist: i32 = i32::MAX;
    let mut evict_dist: i32 = i32::MAX;
    let mut most_similar_idx: i32 = -1;
    for (idx, &dist) in dists.iter().enumerate() {
        min_dist = cpp_min(min_dist, dist);
        if idx > 0 && dist < evict_dist {
            evict_dist = dist;
            most_similar_idx = idx as i32;
        }
    }
    // `most_similar_idx < 0` means the best entry is the only one there
    // is, so nothing may be evicted.
    if far_enough(min_dist, num_integers) && most_similar_idx > 0 {
        return most_similar_idx;
    }
    -1
}

impl SolutionPool {
    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::SolutionPool (lines 11-14)
    pub(crate) fn new(capacity: i32) -> Self {
        let mut pool = Self {
            entries: Vec::new(),
            capacity,
            integer_mask: Some(Rc::new(IntegerColumns::default())),
            generation: 0,
            scan_hook: None,
            best_seen: 0.0,
            has_best_seen: false,
        };
        pool.entries
            .reserve(if capacity > 0 { capacity as usize } else { 0 });
        pool
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::set_integer_mask (lines 16-29)
    pub(crate) fn set_integer_mask(&mut self, mask: Vec<bool>) {
        let mut index: Vec<i32> = Vec::with_capacity(mask.len());
        for (j, &m) in mask.iter().enumerate() {
            if m {
                index.push(j as i32);
            }
        }
        let next = Rc::new(IntegerColumns { mask, index });
        // The old mask is released after the lock; single-threaded, so the
        // swap is the whole operation.
        self.integer_mask = Some(next);
        self.generation += 1;
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::try_add (lines 116-119);
    /// the `const std::vector<double>&` overload: the vector is copied into
    /// an Rc before the offer, offered or not.
    pub(crate) fn try_add_slice(&mut self, obj: f64, sol: &[f64], source: i32) -> PoolAddResult {
        self.try_add(obj, Rc::new(sol.to_vec()), source)
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::try_add (lines 128-304)
    pub(crate) fn try_add(&mut self, obj: f64, sol: Rc<Vec<f64>>, source: i32) -> PoolAddResult {
        // The stored hook keeps its `FnMut()` type from mod.rs; the seam
        // that needs the pool as an argument goes through try_add_hooked.
        let mut stored = self.scan_hook.take();
        let result = self.try_add_hooked(obj, sol, source, &mut |_p: &mut SolutionPool| {
            if let Some(h) = stored.as_mut() {
                h();
            }
        });
        self.scan_hook = stored;
        result
    }

    /// The whole for(;;) loop of try_add, with the hook seam taking the
    /// pool as an argument (the C++ `scan_hook_` seam, solution_pool.h:197).
    ///
    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::try_add (lines 128-304)
    pub(crate) fn try_add_hooked(
        &mut self,
        obj: f64,
        mut sol: Rc<Vec<f64>>,
        source: i32,
        hook: &mut dyn FnMut(&mut SolutionPool),
    ) -> PoolAddResult {
        // A non-finite objective or point is refused outright (#174).
        if !obj.is_finite() || !sol.iter().all(|v| v.is_finite()) {
            return PoolAddResult::default();
        }
        // C++ holds an evicted entry and a replaced mask in locals so they
        // are freed after the lock; has no lock, so they are simply
        // dropped where they leave the pool.
        // The scan's state across attempts. `seen` keeps every entry a scan
        // has measured alive alongside its distance, so its address cannot
        // be reused by a new entry, and a rescan after the pool moved
        // measures only entries it has not seen.
        let mut scan: Vec<PoolEntryRef> = Vec::new();
        let mut seen: Vec<(Rc<Vec<f64>>, i32)> = Vec::new();
        let mut dists: Vec<i32> = Vec::new();
        let mut mask: Option<Rc<IntegerColumns>> = None;
        let mut mask_changed = false;
        let mut scanned_generation: u64 = 0;
        let mut scanned = false;
        // Sized here so the snapshot taken under the lock never allocates.
        scan.reserve(if self.capacity > 0 {
            self.capacity as usize
        } else {
            0
        });
        dists.reserve(if self.capacity > 0 {
            self.capacity as usize
        } else {
            0
        });

        loop {
            let mut result = PoolAddResult::default();
            let mut need_scan = false;
            {
                // The best objective offered so far, decided here rather
                // than by the caller because this is the one point where
                // "the best before this offer" is race-free (#116).  A
                // watermark and not entries_.front(): the diversity path
                // below can replace the front entry, so the pool's best
                // objective goes backwards while the solve's does not.
                // The margin mirrors `improving_offers`; see
                // kImprovementObjMargin.
                if !self.has_best_seen {
                    result.improved_best = true;
                } else {
                    let margin = kImprovementObjMargin * cpp_max(1.0, self.best_seen.abs());
                    result.improved_best = obj < self.best_seen - margin;
                }

                // Find insertion point (entries kept sorted, best first).
                // std::lower_bound with cmp = |e, val| e.objective < val.
                let mut insert_sorted =
                    |entries: &mut Vec<PoolEntryRef>,
                     sol: Rc<Vec<f64>>,
                     result: &mut PoolAddResult| {
                        let pos = entries.partition_point(|e| e.objective < obj);
                        entries.insert(
                            pos,
                            PoolEntryRef {
                                objective: obj,
                                source,
                                solution: sol,
                            },
                        );
                        result.accepted = true;
                    };

                if (self.entries.len() as i32) < self.capacity {
                    insert_sorted(&mut self.entries, Rc::clone(&sol), &mut result);
                } else if !self.entries.is_empty()
                    && obj < self.entries[self.entries.len() - 1].objective
                {
                    // Standard path: improves on worst — replace worst.
                    self.entries.pop();
                    insert_sorted(&mut self.entries, Rc::clone(&sol), &mut result);
                } else if !self.integer_mask.as_ref().unwrap().index.is_empty()
                    && !self.entries.is_empty()
                {
                    // Diversity-aware path: pool is full and obj doesn't
                    // beat worst. Accept if (a) integer mask is set, (b) obj
                    // is within tolerance of best, and (c) solution is
                    // sufficiently diverse from all entries.
                    let best_obj = self.entries[0].objective;
                    let gap = (obj - best_obj).abs();
                    // Continuous fallback: fraction of |best_obj|, floored
                    // to avoid a discontinuous jump near zero.
                    let threshold = cpp_max(
                        kDiversityObjTolerance * best_obj.abs(),
                        kDiversityObjTolerance * 1e-6,
                    );

                    if gap <= threshold {
                        if scanned && scanned_generation == self.generation {
                            // The scan describes these very entries: decide.
                            let victim = diversity_victim(
                                &dists,
                                self.integer_mask.as_ref().unwrap().index.len() as i32,
                            );
                            if victim > 0 {
                                // Replace the most similar entry other than
                                // the best.
                                self.entries.remove(victim as usize);
                                insert_sorted(&mut self.entries, Rc::clone(&sol), &mut result);
                            }
                        } else {
                            // Measure outside the lock, then decide again.
                            scan.clear();
                            scan.extend(self.entries.iter().cloned());
                            // `mask != integer_mask_`: a null `mask` (no scan yet)
                            // differs from the never-null pool mask.
                            let differs = match mask.as_ref() {
                                Some(m) => !Rc::ptr_eq(m, self.integer_mask.as_ref().unwrap()),
                                None => true,
                            };
                            if differs {
                                // Distances measured under another mask do
                                // not stand: `seen` is emptied after the
                                // unlock.
                                mask = self.integer_mask.clone();
                                mask_changed = true;
                            }
                            scanned_generation = self.generation;
                            need_scan = true;
                        }
                    }
                }

                if result.accepted {
                    self.generation += 1;
                    // Advance the watermark only on an offer the pool kept,
                    // so a refused solution cannot raise the bar the next
                    // one is judged against. It never moves backwards.
                    if result.improved_best {
                        self.best_seen = obj;
                        self.has_best_seen = true;
                    }
                }
            }

            if !need_scan {
                // `improved_best` implies `accepted` in every branch above,
                // but tie them together rather than leaving that to a
                // reader of the admission policy: a staleness gate must
                // never reset on a solution the pool did not keep.
                result.improved_best = result.improved_best && result.accepted;
                return result;
            }

            // The scan, outside the lock. Entries are immutable, so a
            // pointer seen before is the same vector and its distance
            // stands — under the same mask, which is why a mask change
            // empties the cache.
            if mask_changed {
                seen.clear();
                mask_changed = false;
            }
            // One entry the offer is too close to refuses it whatever the
            // others measure, so the scan stops there and `dists` ends with
            // it: the verdict over that prefix is the verdict over all of
            // them.
            let num_integers = mask.as_ref().unwrap().index.len() as i32;
            dists.clear();
            for ref_ in scan.iter() {
                let entry = Rc::clone(&ref_.solution);
                let hit = seen.iter().position(|(p, _)| Rc::ptr_eq(p, &entry));
                match hit {
                    Some(h) => dists.push(seen[h].1),
                    None => {
                        let d = hamming_distance(&mask.as_ref().unwrap().index, &sol, &entry);
                        dists.push(d);
                        seen.push((entry, d));
                    }
                }
                if !far_enough(*dists.last().unwrap(), num_integers) {
                    break;
                }
            }
            scanned = true;
            hook(self);
        }
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::snapshot (lines 306-312)
    pub(crate) fn snapshot(&self) -> PoolSnapshot {
        if self.entries.is_empty() {
            return PoolSnapshot {
                has_solution: false,
                best_objective: f64::INFINITY,
            };
        }
        PoolSnapshot {
            has_solution: true,
            best_objective: self.entries[0].objective,
        }
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::is_integer_col (lines 314-316)
    fn is_integer_col(mask: &[bool], j: i32) -> bool {
        (j as usize) < mask.len() && mask[j as usize]
    }

    /// Guided crossover: keep integer values the parents agree on, coin-flip
    /// everything else.
    ///
    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::guided_crossover (lines 320-333)
    fn guided_crossover(
        mask: &[bool],
        sol_a: &[f64],
        sol_b: &[f64],
        rng: &mut S5Rng,
        out: &mut Vec<f64>,
    ) {
        let ncol = sol_a.len() as i32;
        out.resize(ncol as usize, 0.0);
        for j in 0..ncol {
            if Self::is_integer_col(mask, j)
                && sol_a[j as usize].round() == sol_b[j as usize].round()
            {
                out[j as usize] = sol_a[j as usize];
            } else {
                // Disagree or continuous — coin flip.
                out[j as usize] = if s5_uniform_int_i32(rng, 0, 1) == 0 {
                    sol_a[j as usize]
                } else {
                    sol_b[j as usize]
                };
            }
        }
    }

    /// Neighborhood crossover: the better parent provides the base and only
    /// disagreeing integer variables are coin-flipped.
    ///
    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::neighborhood_crossover (lines 337-352)
    fn neighborhood_crossover(
        mask: &[bool],
        sol_better: &[f64],
        sol_other: &[f64],
        rng: &mut S5Rng,
        out: &mut Vec<f64>,
    ) {
        let ncol = sol_better.len() as i32;
        out.resize(ncol as usize, 0.0);
        for j in 0..ncol {
            if Self::is_integer_col(mask, j)
                && sol_better[j as usize].round() != sol_other[j as usize].round()
            {
                out[j as usize] = if s5_uniform_int_i32(rng, 0, 1) == 0 {
                    sol_better[j as usize]
                } else {
                    sol_other[j as usize]
                };
            } else {
                // Agree or continuous — keep the better parent's value.
                out[j as usize] = sol_better[j as usize];
            }
        }
    }

    /// Biased copy: half the time draw from the better half of the pool,
    /// half the time from anywhere.
    ///
    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::biased_index (lines 356-361)
    fn biased_index(rng: &mut S5Rng, pool_size: i32) -> i32 {
        if pool_size > 1 && s5_uniform_int_i32(rng, 0, 1) == 0 {
            return s5_uniform_int_i32(rng, 0, ((pool_size + 1) / 2) - 1);
        }
        s5_uniform_int_i32(rng, 0, pool_size - 1)
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::get_restart (lines 363-426)
    pub(crate) fn get_restart(
        &mut self,
        rng: &mut S5Rng,
        out: &mut Vec<f64>,
        policy: RestartPolicy,
    ) -> bool {
        if policy == RestartPolicy::kBestOnly {
            return self.copy_best(out);
        }
        // The three strategies `get_restart` dispatches between.
        #[allow(clippy::enum_variant_names)] // C++ k-prefixed enumerator names
        #[derive(Clone, Copy, PartialEq, Eq)]
        enum Strategy {
            kGuided,
            kNeighborhood,
            kCopy,
        }
        let mut strategy = Strategy::kCopy;
        let first: Rc<Vec<f64>>;
        let mut second: Option<Rc<Vec<f64>>> = None;
        let mask: Rc<IntegerColumns>;
        {
            // Every RNG draw that picks entries is made here, in the order
            // the strategies have always made them; the crossover's own
            // coin flips follow outside the lock, on the same stream.
            if self.entries.is_empty() {
                return false;
            }

            let pool_size = self.entries.len() as i32;
            let roll = s5_uniform_real(rng, 0.0, 1.0);

            // Two distinct parents, drawn uniformly.
            let pick_two_parents = |rng: &mut S5Rng| -> (i32, i32) {
                let a = s5_uniform_int_i32(rng, 0, pool_size - 1);
                let mut b = s5_uniform_int_i32(rng, 0, pool_size - 2);
                if b >= a {
                    b += 1;
                }
                (a, b)
            };

            if roll < 0.4 && pool_size >= 2 {
                let (a, b) = pick_two_parents(rng);
                strategy = Strategy::kGuided;
                first = Rc::clone(&self.entries[a as usize].solution);
                second = Some(Rc::clone(&self.entries[b as usize].solution));
            } else if roll < 0.7 && pool_size >= 2 {
                let (a, b) = pick_two_parents(rng);
                // Better parent = lower index (entries sorted best-first).
                strategy = Strategy::kNeighborhood;
                first = Rc::clone(&self.entries[cpp_min(a, b) as usize].solution);
                second = Some(Rc::clone(&self.entries[cpp_max(a, b) as usize].solution));
            } else {
                let idx = Self::biased_index(rng, pool_size);
                first = Rc::clone(&self.entries[idx as usize].solution);
            }
            mask = Rc::clone(self.integer_mask.as_ref().unwrap());
        }

        match strategy {
            Strategy::kGuided => {
                let second = second.unwrap();
                Self::guided_crossover(&mask.mask, &first, &second, rng, out);
            }
            Strategy::kNeighborhood => {
                let second = second.unwrap();
                Self::neighborhood_crossover(&mask.mask, &first, &second, rng, out);
            }
            Strategy::kCopy => {
                out.clone_from(&first);
            }
        }
        true
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::entry_refs (lines 428-435)
    pub(crate) fn entry_refs(&self, out: &mut Vec<PoolEntryRef>) {
        // Released and sized before the lock: dropping a stale reference can
        // free a vector, and growing `out` allocates.
        out.clear();
        out.reserve(if self.capacity > 0 {
            self.capacity as usize
        } else {
            0
        });
        out.extend(self.entries.iter().cloned());
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::sorted_entries (lines 437-446)
    pub(crate) fn sorted_entries(&self) -> Vec<PoolEntry> {
        let mut refs: Vec<PoolEntryRef> = Vec::new();
        self.entry_refs(&mut refs);
        let mut entries: Vec<PoolEntry> = Vec::with_capacity(refs.len());
        for ref_ in refs.iter() {
            entries.push(PoolEntry {
                objective: ref_.objective,
                solution: (*ref_.solution).clone(),
                source: ref_.source,
            });
        }
        entries
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::copy_best (lines 448-459)
    pub(crate) fn copy_best(&self, out: &mut Vec<f64>) -> bool {
        if self.entries.is_empty() {
            return false;
        }
        // best — entries kept sorted by try_add.
        *out = (*self.entries[0].solution).clone();
        true
    }

    /// mip-heuristics (36f7d9d): solution_pool.cpp::SolutionPool::size (lines 461-464)
    pub(crate) fn size(&self) -> i32 {
        self.entries.len() as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::fpr_local_mip::rng::s5_uniform_int_u64;

    // Source tags.  Any ints do: the pool stores them and compares nothing.
    const kSourceA: i32 = 1;
    const kSourceB: i32 = 2;
    const kSourceC: i32 = 3;

    // A random offer stream built to exercise every admission path: a small
    // pool, objectives in a narrow band so the diversity tolerance admits
    // dominated offers, and binary columns so their Hamming distances vary.
    struct Offer {
        obj: f64,
        sol: Vec<f64>,
    }

    fn make_offers(count: i32, ncol: i32, seed: u64) -> Vec<Offer> {
        let mut rng = S5Rng::new(seed);
        let mut offers: Vec<Offer> = Vec::with_capacity(count as usize);
        for _ in 0..count {
            let mut offer = Offer {
                obj: s5_uniform_real(&mut rng, 100.0, 108.0),
                sol: vec![0.0; ncol as usize],
            };
            for v in offer.sol.iter_mut() {
                *v = s5_uniform_int_i32(&mut rng, 0, 1) as f64;
            }
            offers.push(offer);
        }
        offers
    }

    #[test]
    fn s5_pool_the_pointer_try_add_decides_every_offer_as_the_vector_one_does() {
        const kNcol: i32 = 40;
        let offers = make_offers(400, kNcol, 7);

        let mut by_ref = SolutionPool::new(4);
        let mut by_ptr = SolutionPool::new(4);
        by_ref.set_integer_mask(vec![true; kNcol as usize]);
        by_ptr.set_integer_mask(vec![true; kNcol as usize]);

        let mut diversity_admissions = 0;
        for offer in offers.iter() {
            let worst_before = if by_ref.size() == 4 {
                by_ref.sorted_entries().last().unwrap().objective
            } else {
                f64::INFINITY
            };
            let a = by_ref.try_add_slice(offer.obj, &offer.sol, kSourceA);
            let b = by_ptr.try_add(offer.obj, Rc::new(offer.sol.clone()), kSourceA);
            assert_eq!(a.accepted, b.accepted);
            assert_eq!(a.improved_best, b.improved_best);
            if a.accepted && offer.obj >= worst_before {
                diversity_admissions += 1;
            }
        }
        // The stream must actually reach the path the scan outside the lock
        // changed, or the comparison above proves nothing about it.
        assert!(diversity_admissions > 0);

        let ea = by_ref.sorted_entries();
        let eb = by_ptr.sorted_entries();
        assert_eq!(ea.len(), eb.len());
        for i in 0..ea.len() {
            assert_eq!(ea[i].objective, eb[i].objective);
            assert_eq!(ea[i].source, eb[i].source);
            assert_eq!(ea[i].solution, eb[i].solution);
        }

        // Same entries and the same RNG state give the same restart draws.
        let mut ra = S5Rng::new(11);
        let mut rb = S5Rng::new(11);
        let mut out_a: Vec<f64> = Vec::new();
        let mut out_b: Vec<f64> = Vec::new();
        for _ in 0..200 {
            assert!(by_ref.get_restart(&mut ra, &mut out_a, RestartPolicy::kMixed));
            assert!(by_ptr.get_restart(&mut rb, &mut out_b, RestartPolicy::kMixed));
            assert_eq!(out_a, out_b);
        }
    }

    #[test]
    fn s5_pool_an_entry_refs_pointer_outlives_the_entry_s_eviction() {
        let mut pool = SolutionPool::new(2);
        assert!(pool.try_add_slice(10.0, &[1.0, 2.0], kSourceA).accepted);
        assert!(pool.try_add_slice(20.0, &[3.0, 4.0], kSourceB).accepted);

        let mut refs: Vec<PoolEntryRef> = Vec::new();
        pool.entry_refs(&mut refs);
        assert_eq!(refs.len(), 2);
        assert_eq!(refs[0].objective, 10.0);
        assert_eq!(refs[1].objective, 20.0);
        assert_eq!(refs[1].source, kSourceB);
        let held = Rc::clone(&refs[1].solution);

        // Both better than 20: the obj-20 entry is evicted on the first.
        assert!(pool.try_add_slice(5.0, &[5.0, 6.0], kSourceC).accepted);
        assert!(pool.try_add_slice(1.0, &[7.0, 8.0], kSourceC).accepted);
        for entry in pool.sorted_entries() {
            assert_ne!(entry.objective, 20.0);
        }
        assert_eq!(*held, vec![3.0, 4.0]);

        // A refill reuses `refs`' capacity and replaces its contents.
        let storage = refs.as_ptr();
        pool.entry_refs(&mut refs);
        assert_eq!(refs.as_ptr(), storage);
        assert_eq!(refs.len(), 2);
        assert_eq!(refs[0].objective, 1.0);
        assert_eq!(*refs[0].solution, vec![7.0, 8.0]);
    }

    #[test]
    fn s5_pool_a_non_finite_objective_or_point_is_refused() {
        // #174: a NaN objective used to be admitted and ranked like any other.
        const kNan: f64 = f64::NAN;
        const kInf: f64 = f64::INFINITY;
        let mut pool = SolutionPool::new(4);
        for &obj in &[kNan, kInf, -kInf] {
            let r = pool.try_add_slice(obj, &[1.0, 2.0], kSourceA);
            assert!(!r.accepted);
            assert!(!r.improved_best);
        }
        for &v in &[kNan, kInf, -kInf] {
            let r = pool.try_add_slice(1.0, &[v, 2.0], kSourceA);
            assert!(!r.accepted);
            assert!(!r.improved_best);
        }
        assert!(pool.sorted_entries().is_empty());
        assert!(!pool.snapshot().has_solution);

        // A finite offer after them is the pool's first.
        let r = pool.try_add_slice(3.0, &[1.0, 2.0], kSourceA);
        assert!(r.accepted);
        assert!(r.improved_best);
    }

    // The C++ concurrent-offers test (test_solution_pool.cpp:157-215) is
    // dropped: the path is single-threaded and the mutex is not
    // ported, so there is nothing to exercise.

    // ── The scan outside the lock, re-decided when the pool moved under it ──

    const kRedecideNcol: i32 = 20;

    // Three entries and an offer built so the offer is dominated (it ties
    // the worst), sits within the diversity tolerance of the best, and is
    // two columns from the best and far from the rest: under an all-integer
    // mask it is diverse enough to be admitted.
    fn seed_three(pool: &mut SolutionPool) {
        let best = vec![0.0; kRedecideNcol as usize];
        let mid = vec![1.0; kRedecideNcol as usize];
        let mut worst = vec![1.0; kRedecideNcol as usize];
        worst[0] = 0.0;
        assert!(pool.try_add_slice(100.0, &best, kSourceA).accepted);
        assert!(pool.try_add_slice(101.0, &mid, kSourceA).accepted);
        assert!(pool.try_add_slice(102.0, &worst, kSourceA).accepted);
    }

    fn near_best_offer() -> Vec<f64> {
        let mut offer = vec![0.0; kRedecideNcol as usize];
        offer[0] = 1.0;
        offer[1] = 1.0;
        offer
    }

    fn require_same_pool(a: &mut SolutionPool, b: &mut SolutionPool) {
        let ea = a.sorted_entries();
        let eb = b.sorted_entries();
        assert_eq!(ea.len(), eb.len());
        for i in 0..ea.len() {
            assert_eq!(ea[i].objective, eb[i].objective);
            assert_eq!(ea[i].source, eb[i].source);
            assert_eq!(ea[i].solution, eb[i].solution);
        }
    }

    #[test]
    fn s5_pool_a_mask_changed_mid_scan_is_decided_under_the_new_mask() {
        // The mask the offer ends up judged under marks only the columns
        // where it agrees with the best entry, so it is no longer diverse.
        // Distances measured under the old mask must not be reused against
        // the new one.
        let mut narrow = vec![true; kRedecideNcol as usize];
        narrow[0] = false;
        narrow[1] = false;

        let mut pool = SolutionPool::new(3);
        pool.set_integer_mask(vec![true; kRedecideNcol as usize]);
        seed_three(&mut pool);
        let mut hook_calls = 0;
        let moved = pool.try_add_hooked(
            102.0,
            Rc::new(near_best_offer()),
            kSourceB,
            &mut |p: &mut SolutionPool| {
                if hook_calls == 0 {
                    hook_calls += 1;
                    p.set_integer_mask(narrow.clone());
                } else {
                    hook_calls += 1;
                }
            },
        );
        assert_eq!(hook_calls, 2); // scanned, then rescanned under the new mask

        // One hold of the lock in the final state decides it this way.
        let mut reference = SolutionPool::new(3);
        reference.set_integer_mask(narrow);
        seed_three(&mut reference);
        let direct = reference.try_add_slice(102.0, &near_best_offer(), kSourceB);

        assert!(!direct.accepted);
        assert_eq!(moved.accepted, direct.accepted);
        assert_eq!(moved.improved_best, direct.improved_best);
        require_same_pool(&mut pool, &mut reference);
    }

    #[test]
    fn s5_pool_entries_changed_mid_scan_are_rescanned_before_the_decision() {
        // A peer's solution lands during the scan: the offer's own vector,
        // at an objective inside the diversity band. The decision must then
        // go through a rescan against the entries as they now stand, where
        // the offer duplicates the new best and is refused. Applying the
        // first scan's distances to the new entries instead admits it (they
        // say it is two columns from its nearest entry) and evicts by a
        // stale index.
        let mut pool = SolutionPool::new(3);
        pool.set_integer_mask(vec![true; kRedecideNcol as usize]);
        seed_three(&mut pool);
        let mut hook_calls = 0;
        let moved = pool.try_add_hooked(
            102.0,
            Rc::new(near_best_offer()),
            kSourceB,
            &mut |p: &mut SolutionPool| {
                if hook_calls == 0 {
                    hook_calls += 1;
                    assert!(p.try_add_slice(99.5, &near_best_offer(), kSourceC).accepted);
                } else {
                    hook_calls += 1;
                }
            },
        );
        assert_eq!(hook_calls, 2); // scanned, then rescanned against the new entries

        let mut reference = SolutionPool::new(3);
        reference.set_integer_mask(vec![true; kRedecideNcol as usize]);
        seed_three(&mut reference);
        assert!(
            reference
                .try_add_slice(99.5, &near_best_offer(), kSourceC)
                .accepted
        );
        let direct = reference.try_add_slice(102.0, &near_best_offer(), kSourceB);

        assert!(!direct.accepted);
        assert_eq!(moved.accepted, direct.accepted);
        assert_eq!(moved.improved_best, direct.improved_best);
        require_same_pool(&mut pool, &mut reference);

        // Without the move, the same offer is admitted: the first scan's
        // verdict.
        let mut unmoved = SolutionPool::new(3);
        unmoved.set_integer_mask(vec![true; kRedecideNcol as usize]);
        seed_three(&mut unmoved);
        assert!(
            unmoved
                .try_add_slice(102.0, &near_best_offer(), kSourceB)
                .accepted
        );
    }

    #[test]
    fn s5_pool_a_best_only_restart_is_the_best_entry_drawing_nothing() {
        let mut pool = SolutionPool::new(4);
        let mut out: Vec<f64> = Vec::new();
        let mut rng = S5Rng::new(7);
        assert!(!pool.get_restart(&mut rng, &mut out, RestartPolicy::kBestOnly));
        assert!(
            pool.try_add_slice(20.0, &[0.0, 1.0, 0.0], kSourceA)
                .accepted
        );
        assert!(
            pool.try_add_slice(10.0, &[1.0, 1.0, 0.0], kSourceB)
                .accepted
        );
        assert!(
            pool.try_add_slice(30.0, &[0.0, 0.0, 1.0], kSourceC)
                .accepted
        );

        let mut untouched = rng.clone();
        for _ in 0..50 {
            assert!(pool.get_restart(&mut rng, &mut out, RestartPolicy::kBestOnly));
            assert_eq!(out, vec![1.0, 1.0, 0.0]);
        }
        assert_eq!(rng.next_u64(), untouched.next_u64()); // the same stream: nothing was drawn
        // A new best is what the next restart draws.
        assert!(pool.try_add_slice(5.0, &[1.0, 0.0, 1.0], kSourceA).accepted);
        assert!(pool.get_restart(&mut rng, &mut out, RestartPolicy::kBestOnly));
        assert_eq!(out, vec![1.0, 0.0, 1.0]);
    }

    // The diversity rule as it stood before #195 made the scan fast, kept
    // here as the reference the pool is held to: every entry measured,
    // std::round on every column the flag array marks.
    struct ReferencePool {
        capacity: i32,
        mask: Vec<bool>,
        num_integers: i32,
        diversity_offers: i32,
        diversity_accepts: i32,
        entries: Vec<(f64, Vec<f64>)>,
    }

    impl ReferencePool {
        fn new(capacity: i32, mask: Vec<bool>) -> Self {
            let mut num_integers = 0;
            for &b in mask.iter() {
                num_integers += b as i32;
            }
            Self {
                capacity,
                mask,
                num_integers,
                diversity_offers: 0,
                diversity_accepts: 0,
                entries: Vec::new(),
            }
        }

        fn try_add(&mut self, obj: f64, sol: &[f64]) -> bool {
            if (self.entries.len() as i32) < self.capacity {
                self.insert(obj, sol);
                return true;
            }
            if obj < self.entries.last().unwrap().0 {
                self.entries.pop();
                self.insert(obj, sol);
                return true;
            }
            let best = self.entries[0].0;
            if self.num_integers == 0
                || (obj - best).abs()
                    > cpp_max(
                        kDiversityObjTolerance * best.abs(),
                        kDiversityObjTolerance * 1e-6,
                    )
            {
                return false;
            }
            self.diversity_offers += 1;
            let mut min_dist = i32::MAX;
            let mut evict_dist = i32::MAX;
            let mut victim = -1;
            for idx in 0..self.entries.len() {
                let mut dist = 0;
                for j in 0..sol.len() {
                    if self.mask[j] && sol[j].round() != self.entries[idx].1[j].round() {
                        dist += 1;
                    }
                }
                min_dist = cpp_min(min_dist, dist);
                if idx > 0 && dist < evict_dist {
                    evict_dist = dist;
                    victim = idx as i32;
                }
            }
            if (min_dist as f64) / (self.num_integers as f64) < kDiversityMinHammingFrac
                || victim <= 0
            {
                return false;
            }
            self.diversity_accepts += 1;
            self.entries.remove(victim as usize);
            self.insert(obj, sol);
            true
        }

        fn insert(&mut self, obj: f64, sol: &[f64]) {
            let pos = self.entries.partition_point(|e| e.0 < obj);
            self.entries.insert(pos, (obj, sol.to_vec()));
        }
    }

    #[test]
    fn s5_pool_the_diversity_scan_decides_every_offer_as_the_full_std_round_scan_does() {
        // Mixed columns, a third continuous, and values on and off the
        // rounding ties, negative ones included, so the distances depend on
        // how a value rounds and on which columns count; objectives in a
        // band narrow enough that most offers take the diversity path,
        // where the scan stops at the first entry too close to the offer.
        const kNcol: usize = 30;
        const kOffers: i32 = 3000;
        let mut mask = vec![false; kNcol];
        for j in 0..kNcol {
            mask[j] = j % 3 != 1;
        }
        let values: Vec<f64> = vec![
            -2.5,
            -1.5,
            -1.0,
            -0.5,
            -0.49999999999999994,
            -0.0,
            0.0,
            0.3,
            0.5,
            0.7,
            1.0,
            1.5,
            2.5,
            0.49999999999999994,
        ];
        let mut rng = S5Rng::new(11);
        let pick_value = |rng: &mut S5Rng| -> f64 {
            values[s5_uniform_int_u64(rng, 0, values.len() as u64 - 1) as usize]
        };
        // Offers are a few bank vectors with a column or two redrawn, so
        // they sit near the entries and the scan refuses as well as admits.
        let mut bank: Vec<Vec<f64>> = vec![vec![0.0; kNcol]; 4];
        for base in bank.iter_mut() {
            for v in base.iter_mut() {
                *v = pick_value(&mut rng);
            }
        }
        let mut pool = SolutionPool::new(5);
        pool.set_integer_mask(mask.clone());
        let mut reference = ReferencePool::new(5, mask);
        for _ in 0..kOffers {
            let obj = s5_uniform_real(&mut rng, 100.0, 104.0);
            let mut sol =
                bank[s5_uniform_int_u64(&mut rng, 0, bank.len() as u64 - 1) as usize].clone();
            for v in sol.iter_mut() {
                if s5_uniform_int_i32(&mut rng, 0, kNcol as i32 - 1) == 0 {
                    *v = pick_value(&mut rng);
                }
            }
            let expected = reference.try_add(obj, &sol);
            assert_eq!(pool.try_add_slice(obj, &sol, kSourceA).accepted, expected);
        }
        // Both diversity verdicts were reached, many times over.
        assert!(reference.diversity_accepts > 50);
        assert!(reference.diversity_offers - reference.diversity_accepts > 50);
        let entries = pool.sorted_entries();
        assert_eq!(entries.len(), reference.entries.len());
        for i in 0..entries.len() {
            assert_eq!(entries[i].objective, reference.entries[i].0);
            assert_eq!(entries[i].solution, reference.entries[i].1);
        }
    }

    #[test]
    fn s5_pool_a_column_s_hamming_term_is_std_round_s_verdict_at_every_rounding_edge() {
        // Twenty integer columns, so one differing column is exactly the 5%
        // diversity threshold. The best entry and the offer agree everywhere
        // but column 0, where they hold `y` and `x`; the worst entry also
        // differs in column 19, so it is never the offer's nearest. The
        // offer, inside the diversity band and not beating the worst, is
        // admitted exactly when `x` and `y` round apart.
        const kNcol: usize = 21; // the last column continuous
        let mut mask = vec![true; kNcol];
        mask[kNcol - 1] = false;
        let big: f64 = 4503599627370496.0; // 2^52
        let values: Vec<f64> = vec![
            0.0,
            -0.0,
            0.49999999999999994,
            0.5,
            0.5000000000000001,
            -0.49999999999999994,
            -0.5,
            1.5,
            2.5,
            -2.5,
            1e-300,
            big - 0.5,
            big - 1.0,
            big,
            big + 2.0,
            -(big - 0.5),
            1e300,
            -1e300,
        ];
        for &x in values.iter() {
            for &y in values.iter() {
                let mut pool = SolutionPool::new(2);
                pool.set_integer_mask(mask.clone());
                let mut best = vec![0.0; kNcol];
                best[0] = y;
                let mut worst = best.clone();
                worst[19] = 7.0;
                assert!(pool.try_add_slice(1.0, &best, kSourceA).accepted);
                assert!(pool.try_add_slice(1.05, &worst, kSourceB).accepted);
                let mut offer = best.clone();
                offer[0] = x;
                offer[kNcol - 1] = 3.25; // continuous: never counted
                let got = pool.try_add_slice(1.05, &offer, kSourceC).accepted;
                let want = x.round() != y.round();
                assert_eq!(
                    got, want,
                    "x = {x:x?}, y = {y:x?}: accepted {got}, want {want}"
                );
            }
        }
    }

    // The manifest's hand-checked admission and restart cases.

    #[test]
    fn s5_pool_try_add_refuses_a_worse_offer_at_capacity_without_a_mask() {
        let mut pool = SolutionPool::new(2);
        assert!(pool.try_add_slice(10.0, &[0.0, 0.0], kSourceA).accepted);
        assert!(pool.try_add_slice(20.0, &[1.0, 1.0], kSourceB).accepted);
        assert_eq!(pool.size(), 2);
        // Full, worse than the worst, and no integer mask: refused.
        let r = pool.try_add_slice(30.0, &[1.0, 0.0], kSourceC);
        assert!(!r.accepted);
        assert!(!r.improved_best);
        assert_eq!(pool.size(), 2);
        // Better than the worst: the worst is replaced.
        let r = pool.try_add_slice(15.0, &[0.0, 1.0], kSourceC);
        assert!(r.accepted);
        assert!(!r.improved_best);
        let entries = pool.sorted_entries();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].objective, 10.0);
        assert_eq!(entries[1].objective, 15.0);
        assert_eq!(entries[1].source, kSourceC);
    }

    #[test]
    fn s5_pool_try_add_admits_a_diverse_near_best_offer_and_refuses_a_duplicate() {
        const kNcol: usize = 20;
        let mut pool = SolutionPool::new(2);
        pool.set_integer_mask(vec![true; kNcol]);
        let best = vec![0.0; kNcol];
        let mut worst = best.clone();
        worst[0] = 1.0;
        assert!(pool.try_add_slice(1.0, &best, kSourceA).accepted);
        assert!(pool.try_add_slice(1.05, &worst, kSourceB).accepted);

        // Within tolerance of the best (gap 0.06 <= 0.1) but a copy of it:
        // Hamming distance 0, refused.
        assert!(!pool.try_add_slice(1.06, &best, kSourceC).accepted);
        // Outside tolerance (gap 0.2 > 0.1): refused however diverse.
        let mut far = best.clone();
        for v in far.iter_mut().take(10) {
            *v = 1.0;
        }
        assert!(!pool.try_add_slice(1.2, &far, kSourceC).accepted);

        // Distances 5 (to the best) and 4 (to the worst): 4/20 >= 0.05, so
        // the most similar non-best entry, the worst, is replaced.
        let mut diverse = best.clone();
        for v in diverse.iter_mut().take(5) {
            *v = 1.0;
        }
        let r = pool.try_add_slice(1.08, &diverse, kSourceC);
        assert!(r.accepted);
        assert!(!r.improved_best);
        let entries = pool.sorted_entries();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].objective, 1.0);
        assert_eq!(entries[1].objective, 1.08);
        assert_eq!(entries[1].solution, diverse);

        // The same offer with no integer mask is refused: no diversity path.
        let mut unmasked = SolutionPool::new(2);
        assert!(unmasked.try_add_slice(1.0, &best, kSourceA).accepted);
        assert!(unmasked.try_add_slice(1.05, &worst, kSourceB).accepted);
        assert!(!unmasked.try_add_slice(1.08, &diverse, kSourceC).accepted);
    }

    #[test]
    fn s5_pool_improved_best_needs_the_relative_margin() {
        let mut pool = SolutionPool::new(10);
        let r = pool.try_add_slice(1000.0, &[0.0], kSourceA);
        assert!(r.accepted && r.improved_best);
        // margin = 1e-9 * max(1, 1000) = 1e-6: 1e-7 better is kept but is
        // not an improvement, and does not move the watermark.
        let r = pool.try_add_slice(1000.0 - 1e-7, &[1.0], kSourceA);
        assert!(r.accepted);
        assert!(!r.improved_best);
        let r = pool.try_add_slice(1000.0 - 9e-7, &[2.0], kSourceA);
        assert!(r.accepted);
        assert!(!r.improved_best);
        let r = pool.try_add_slice(1000.0 - 2e-6, &[3.0], kSourceA);
        assert!(r.accepted && r.improved_best);
        // Near zero the margin floors at 1e-9 absolute.
        let mut small = SolutionPool::new(10);
        assert!(small.try_add_slice(0.0, &[0.0], kSourceA).improved_best);
        assert!(!small.try_add_slice(-5e-10, &[1.0], kSourceA).improved_best);
        assert!(small.try_add_slice(-2e-9, &[2.0], kSourceA).improved_best);
        // A refused offer never reports an improvement.
        let mut full = SolutionPool::new(1);
        assert!(full.try_add_slice(5.0, &[0.0], kSourceA).improved_best);
        let r = full.try_add_slice(6.0, &[1.0], kSourceA);
        assert!(!r.accepted && !r.improved_best);
    }

    // rng.txt, raw seed=42: uniform_int(0,1) = 1 0 1 1 1 1 0 1 ...,
    // uniform_int(0,4) and friends follow the same libstdc++ Lemire draw.
    #[test]
    fn s5_pool_biased_index_follows_the_seed_42_draws() {
        let mut rng = S5Rng::new(42);
        let got: Vec<i32> = (0..8)
            .map(|_| SolutionPool::biased_index(&mut rng, 5))
            .collect();
        assert_eq!(got, vec![1, 3, 2, 1, 2, 4, 0, 1]);
        // Sixteen draws: every coin was 0 or 1 followed by one index draw.
        let mut replay = S5Rng::new(42);
        for _ in 0..16 {
            replay.next_u64();
        }
        assert_eq!(rng.next_u64(), replay.next_u64());

        // A pool of one draws no coin, only the index uniform_int(0, 0).
        let mut rng = S5Rng::new(42);
        let got: Vec<i32> = (0..3)
            .map(|_| SolutionPool::biased_index(&mut rng, 1))
            .collect();
        assert_eq!(got, vec![0, 0, 0]);
        let mut replay = S5Rng::new(42);
        for _ in 0..3 {
            replay.next_u64();
        }
        assert_eq!(rng.next_u64(), replay.next_u64());
    }

    // Hand trace on rng.txt raw seed=42 (canonical = raw / 2^64):
    //   call 1: roll 0.814 -> copy; coin 0 -> better half; uniform(0,1)=1 -> E1
    //   call 2: roll 0.701 -> copy; coin 1 -> whole pool; uniform(0,2)=1 -> E1
    //   call 3: roll 0.125 -> guided; a=1, b=0 (b<a, kept): E1 x E0;
    //           j0, j1 disagree and j3 is continuous: flips 1 1 1 -> E0's
    //           values; j2 agrees (0) and is kept
    //   call 4: roll 0.680 -> neighborhood; a=0, b=0 -> 1: E0 base, E1 other;
    //           j0, j1 disagree: flips 1 0 -> 1, 0; j2 agrees; j3 E0's
    //   call 5: roll 0.050 -> guided; a=1, b=0: flips 0 1 1 -> 1, 0, 0.25
    #[test]
    fn s5_pool_get_restart_follows_the_seed_42_hand_trace() {
        let mut pool = SolutionPool::new(10);
        pool.set_integer_mask(vec![true, true, true, false]);
        assert!(
            pool.try_add_slice(2.0, &[1.0, 1.0, 0.0, 0.5], kSourceB)
                .accepted
        );
        assert!(
            pool.try_add_slice(3.0, &[0.0, 1.0, 1.0, 0.75], kSourceC)
                .accepted
        );
        assert!(
            pool.try_add_slice(1.0, &[0.0, 0.0, 0.0, 0.25], kSourceA)
                .accepted
        );

        let mut rng = S5Rng::new(42);
        let mut out: Vec<f64> = Vec::new();
        let want: [[f64; 4]; 5] = [
            [1.0, 1.0, 0.0, 0.5],
            [1.0, 1.0, 0.0, 0.5],
            [0.0, 0.0, 0.0, 0.25],
            [1.0, 0.0, 0.0, 0.25],
            [1.0, 0.0, 0.0, 0.25],
        ];
        for (k, w) in want.iter().enumerate() {
            assert!(pool.get_restart(&mut rng, &mut out, RestartPolicy::kMixed));
            assert_eq!(out, w.to_vec(), "call {}", k + 1);
        }
        // 23 raw draws in all: 3 + 3 + 6 + 5 + 6.
        let mut replay = S5Rng::new(42);
        for _ in 0..23 {
            replay.next_u64();
        }
        assert_eq!(rng.next_u64(), replay.next_u64());

        // An empty pool draws nothing and reports no restart.
        let mut empty = SolutionPool::new(10);
        let mut rng = S5Rng::new(42);
        assert!(!empty.get_restart(&mut rng, &mut out, RestartPolicy::kMixed));
        assert_eq!(rng.next_u64(), S5Rng::new(42).next_u64());
    }
}
