//! Cut pool: stores every generated cut, rejects duplicates, and picks the
//! cuts that enter the LP each round by an efficacy score with an adaptive
//! threshold and a parallelism filter.
//!
//! Follows HiGHS `HighsCutPool` (reference/HiGHS/highs/mip/HighsCutPool.cpp,
//! MIT; read as a reference text and reimplemented): duplicate test
//! `isDuplicate` :47–87 (same support, parallelism ≥ 1 − 1e-6), `addCut`
//! :503–638 (new cuts start at age `agelim − 5`), `separate` :229–460
//! (score = violation / (active nonzeros · ‖active part‖), minimum score
//! `0.9 ·` best observed score adapted between 1/20 and all of the candidates,
//! parallelism above 0.1 against an already selected cut discards),
//! `lpCutRemoved` :150–160. Efficacy-plus-parallelism selection is surveyed
//! in Wesselmann & Suhl, "Implementing cutting plane management and selection
//! techniques", Tech. report, U. Paderborn (2012), and Achterberg (2007) §8.9.

use super::CutClass;
use std::collections::HashMap;

/// HiGHS `mip_pool_age_limit` default (HighsOptions.h:1186).
const AGE_LIMIT: i32 = 30;
/// HiGHS `mip_pool_soft_limit` default.
const SOFT_LIMIT: usize = 10_000;
const MAX_PARALLELISM: f64 = 0.1;

pub(crate) struct PoolCut {
    /// Sorted ascending.
    pub inds: Vec<usize>,
    pub vals: Vec<f64>,
    pub rhs: f64,
    pub integral: bool,
    pub class: CutClass,
    /// 1/‖vals‖.
    norm_inv: f64,
    /// −1: in the LP; otherwise rounds since last violated.
    pub age: i32,
    pub deleted: bool,
}

pub(crate) struct CutPool {
    pub cuts: Vec<PoolCut>,
    by_support: HashMap<u64, Vec<usize>>,
    best_observed_score: f64,
    min_score_factor: f64,
}

fn support_hash(inds: &[usize]) -> u64 {
    // FNV-1a over the sorted support.
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &j in inds {
        h ^= j as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

impl CutPool {
    pub fn new() -> Self {
        Self {
            cuts: Vec::new(),
            by_support: HashMap::new(),
            best_observed_score: 0.0,
            min_score_factor: 0.9,
        }
    }

    /// Adds `Σ vals·x[inds] ≤ rhs` unless a (near-)parallel cut with the same
    /// support is stored. Returns the pool index of a new cut.
    pub fn add(
        &mut self,
        inds: &[usize],
        vals: &[f64],
        rhs: f64,
        integral: bool,
        class: CutClass,
    ) -> Option<usize> {
        let mut entries: Vec<(usize, f64)> =
            inds.iter().copied().zip(vals.iter().copied()).collect();
        entries.sort_unstable_by_key(|e| e.0);
        let inds: Vec<usize> = entries.iter().map(|e| e.0).collect();
        let vals: Vec<f64> = entries.iter().map(|e| e.1).collect();
        let norm = vals.iter().map(|v| v * v).sum::<f64>().sqrt();
        if inds.is_empty() || norm.is_nan() || norm <= 0.0 || !rhs.is_finite() {
            return None;
        }
        let norm_inv = 1.0 / norm;
        let h = support_hash(&inds);
        if let Some(list) = self.by_support.get(&h) {
            for &c in list {
                let other = &self.cuts[c];
                if other.deleted || other.inds != inds {
                    continue;
                }
                let dot: f64 = other.vals.iter().zip(&vals).map(|(a, b)| a * b).sum();
                if dot * other.norm_inv * norm_inv >= 1.0 - 1e-6 {
                    return None;
                }
            }
        }
        let idx = self.cuts.len();
        self.by_support.entry(h).or_default().push(idx);
        self.cuts.push(PoolCut {
            inds,
            vals,
            rhs,
            integral,
            class,
            norm_inv,
            age: (AGE_LIMIT - 5).max(0),
            deleted: false,
        });
        Some(idx)
    }

    fn parallelism(&self, a: usize, b: usize) -> f64 {
        let (ca, cb) = (&self.cuts[a], &self.cuts[b]);
        let (mut i, mut k, mut dot) = (0, 0, 0.0);
        while i < ca.inds.len() && k < cb.inds.len() {
            match ca.inds[i].cmp(&cb.inds[k]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => k += 1,
                std::cmp::Ordering::Equal => {
                    dot += ca.vals[i] * cb.vals[k];
                    i += 1;
                    k += 1;
                }
            }
        }
        dot * ca.norm_inv * cb.norm_inv
    }

    fn delete(&mut self, c: usize) {
        let h = support_hash(&self.cuts[c].inds);
        if let Some(list) = self.by_support.get_mut(&h) {
            list.retain(|&k| k != c);
        }
        let cut = &mut self.cuts[c];
        cut.deleted = true;
        cut.inds = Vec::new();
        cut.vals = Vec::new();
    }

    /// Selects the cuts to add to the LP at point `x` (marks them as LP cuts).
    pub fn separate(
        &mut self,
        x: &[f64],
        lower: &[f64],
        upper: &[f64],
        feastol: f64,
    ) -> Vec<usize> {
        let available = self
            .cuts
            .iter()
            .filter(|c| !c.deleted && c.age >= 0)
            .count();
        let mut agelim = AGE_LIMIT;
        if available > SOFT_LIMIT {
            // Shrink the age limit until the pool fits (HighsCutPool.cpp:243).
            let mut hist = vec![0usize; AGE_LIMIT as usize + 1];
            for c in self.cuts.iter().filter(|c| !c.deleted && c.age >= 0) {
                hist[(c.age.min(AGE_LIMIT)) as usize] += 1;
            }
            let mut n = available;
            while agelim > 1 && n > SOFT_LIMIT {
                n -= hist[agelim as usize];
                agelim -= 1;
            }
        }
        let mut scored: Vec<(f64, usize)> = Vec::new();
        for c in 0..self.cuts.len() {
            let cut = &self.cuts[c];
            if cut.deleted || cut.age < 0 {
                continue;
            }
            let viol = cut
                .inds
                .iter()
                .zip(&cut.vals)
                .map(|(&j, &v)| v * x[j])
                .sum::<f64>()
                - cut.rhs;
            if viol <= feastol {
                let cut = &mut self.cuts[c];
                cut.age += 1;
                if cut.age >= agelim {
                    self.delete(c);
                }
                continue;
            }
            // Norm over entries away from the bound that minimizes their term.
            let (mut sq, mut active) = (0.0, 0usize);
            for (&j, &v) in cut.inds.iter().zip(&cut.vals) {
                let away = if v > 0.0 {
                    x[j] > lower[j] + feastol
                } else {
                    x[j] < upper[j] - feastol
                };
                if away {
                    sq += v * v;
                    active += 1;
                }
            }
            self.cuts[c].age = 0;
            let score = viol / (active.max(1) as f64 * sq.sqrt().max(f64::MIN_POSITIVE));
            scored.push((score, c));
        }
        if scored.is_empty() {
            return Vec::new();
        }
        // Best score first; ties by index (HiGHS breaks them by a hash).
        scored.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        self.best_observed_score = self.best_observed_score.max(scored[0].0);
        let min_score = self.min_score_factor * self.best_observed_score;
        let mut num = scored.partition_point(|s| s.0 >= min_score);
        let lower_threshold = scored.len() / 20;
        let upper_threshold = scored.len() - 1;
        if num <= lower_threshold {
            num = (scored.len() / 2).max(1);
            self.min_score_factor = scored[num - 1].0 / self.best_observed_score;
        } else if num > upper_threshold {
            self.min_score_factor = scored[upper_threshold].0 / self.best_observed_score;
        }
        scored.truncate(num);

        let mut selected: Vec<usize> = Vec::new();
        for &(_, c) in &scored {
            if selected
                .iter()
                .any(|&s| self.parallelism(s, c) > MAX_PARALLELISM)
            {
                continue;
            }
            self.cuts[c].age = -1;
            selected.push(c);
        }
        selected
    }

    /// A cut left the LP: it competes again from age 1.
    pub fn lp_cut_removed(&mut self, c: usize) {
        if !self.cuts[c].deleted {
            self.cuts[c].age = 1;
        }
    }
}
