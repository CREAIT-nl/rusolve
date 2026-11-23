//! Pseudocost, inference and cutoff statistics per integer column.
//!
//! A pseudocost is the average objective gain per unit of bound change seen
//! when branching a column up or down (Bénichou et al., "Experiments in
//! mixed-integer linear programming", Math. Prog. 1:76–94, 1971). Reliability
//! branching (Achterberg, Koch & Martin, "Branching rules revisited", Oper.
//! Res. Lett. 33:42–54, 2005) trusts a column's pseudocost once it has
//! `min_reliable` observations per direction.
//!
//! Reimplemented from HiGHS `HighsPseudocost` (reference/HiGHS/highs/mip/
//! HighsPseudocost.h:187–417, MIT licence, Copyright (c) HiGHS developers),
//! read as a reference text: running means for costs (`addObservation` :187),
//! inferences (`addInferenceObservation` :220), cutoff counts
//! (`addCutoffObservation` :175) and conflict scores
//! (`increaseConflictScoreUp/Down` :147–159); the unreliable-cost blend
//! towards the global average (`getPseudocostUp/Down` with offset, :257–288);
//! the product score normalised by the global average and mapped through
//! `1 − 1/(1+s)` with weights 1 / 1e-2 / 1e-4 / 1e-4 for cost / conflict /
//! cutoff / inference (`getScore` :310–349); the one-sided scores used for
//! plunge decisions (`getScoreUp/Down` :359–417). Default
//! `min_reliable = 8` is HiGHS's `mip_pscost_minreliable`
//! (lp_data/HighsOptions.h:1197).
//!
//! The conflict scores are credited by the conflict analysis of
//! [`crate::conflict`]: every explained conflict grows a decaying weight
//! (`increaseConflictWeight` :123–140, ×1.02 per conflict, rescaled above
//! 1000 so recent conflicts dominate — hybrid/conflict branching,
//! Achterberg & Berthold, "Hybrid branching", CPAIOR 2009; Witzig, Berthold
//! & Heinz, "Experiments with conflict analysis in mixed integer
//! programming", CPAIOR 2017) and adds it to every literal's column on the
//! literal's side.
//!
//! Not yet ported: the LP degeneracy factor (needs node-LP duals;
//! HiGHS `computeLPDegneracy`, HighsLpRelaxation.cpp:452), which in HiGHS
//! divides the cost term and multiplies the conflict/cutoff/inference ones.

/// Floor for every factor and denominator in the scores (HighsPseudocost.h:27).
const MIN_THRESHOLD: f64 = 1e-6;
/// Default observations per direction before a pseudocost is trusted.
pub const DEFAULT_MIN_RELIABLE: u32 = 8;

#[derive(Debug, Clone)]
pub struct Pseudocost {
    cost_up: Vec<f64>,
    cost_down: Vec<f64>,
    samples_up: Vec<u32>,
    samples_down: Vec<u32>,
    inferences_up: Vec<f64>,
    inferences_down: Vec<f64>,
    ninferences_up: Vec<u32>,
    ninferences_down: Vec<u32>,
    cutoffs_up: Vec<u32>,
    cutoffs_down: Vec<u32>,
    /// Mean unit gain over every observation of every column.
    cost_total: f64,
    /// Mean inference count over every observation.
    inferences_total: f64,
    nsamples_total: u64,
    ninferences_total: u64,
    ncutoffs_total: u64,
    min_reliable: u32,
    /// Conflict scores per column and side, in units of the decaying
    /// weight (HiGHS `conflictscoreup/down`, HighsPseudocost.h:38–39).
    conflict_up: Vec<f64>,
    conflict_down: Vec<f64>,
    /// Discount of the stored conflict scores: ×1.02 per explained
    /// conflict, everything rescaled above 1000, so recent conflicts
    /// dominate (`increaseConflictWeight` :123–140).
    conflict_weight: f64,
    /// Sum of the weights of all recorded conflict literals
    /// (`conflict_avg_score` :42; the score divides by its mean per column).
    conflict_avg: f64,
}

impl Default for Pseudocost {
    fn default() -> Self {
        Self {
            cost_up: Vec::new(),
            cost_down: Vec::new(),
            samples_up: Vec::new(),
            samples_down: Vec::new(),
            inferences_up: Vec::new(),
            inferences_down: Vec::new(),
            ninferences_up: Vec::new(),
            ninferences_down: Vec::new(),
            cutoffs_up: Vec::new(),
            cutoffs_down: Vec::new(),
            cost_total: 0.0,
            inferences_total: 0.0,
            nsamples_total: 0,
            ninferences_total: 0,
            ncutoffs_total: 0,
            min_reliable: DEFAULT_MIN_RELIABLE,
            conflict_up: Vec::new(),
            conflict_down: Vec::new(),
            // HiGHS starts the decay weight at 1.0 (HighsPseudocost.cpp:26).
            conflict_weight: 1.0,
            conflict_avg: 0.0,
        }
    }
}

/// Running mean update `mean += (v − mean)/n` after incrementing `n`.
fn push_mean(mean: &mut f64, n: &mut u32, v: f64) {
    *n += 1;
    *mean += (v - *mean) / f64::from(*n);
}

fn push_mean_total(mean: &mut f64, n: &mut u64, v: f64) {
    *n += 1;
    *mean += (v - *mean) / *n as f64;
}

/// Maps a nonnegative ratio into `[0, 1)` (HighsPseudocost.h:345).
fn map_score(s: f64) -> f64 {
    1.0 - 1.0 / (1.0 + s)
}

impl Pseudocost {
    pub fn new(num_col: usize) -> Self {
        Self {
            cost_up: vec![0.0; num_col],
            cost_down: vec![0.0; num_col],
            samples_up: vec![0; num_col],
            samples_down: vec![0; num_col],
            inferences_up: vec![0.0; num_col],
            inferences_down: vec![0.0; num_col],
            ninferences_up: vec![0; num_col],
            ninferences_down: vec![0; num_col],
            cutoffs_up: vec![0; num_col],
            cutoffs_down: vec![0; num_col],
            cost_total: 0.0,
            inferences_total: 0.0,
            nsamples_total: 0,
            ninferences_total: 0,
            ncutoffs_total: 0,
            min_reliable: DEFAULT_MIN_RELIABLE,
            conflict_up: vec![0.0; num_col],
            conflict_down: vec![0.0; num_col],
            conflict_weight: 1.0,
            conflict_avg: 0.0,
        }
    }

    /// Discounts the conflict scores once per explained conflict: ×1.02, and
    /// above 1000 everything is rescaled back, so recent conflicts dominate
    /// the branching score (`increaseConflictWeight`, HighsPseudocost.h:123–140).
    pub fn increase_conflict_weight(&mut self) {
        self.conflict_weight *= 1.02;
        if self.conflict_weight > 1000.0 {
            let scale = 1.0 / self.conflict_weight;
            self.conflict_weight = 1.0;
            self.conflict_avg *= scale;
            for (u, d) in self.conflict_up.iter_mut().zip(&mut self.conflict_down) {
                *u *= scale;
                *d *= scale;
            }
        }
    }

    /// Credits the current conflict weight to `col` on the side of the
    /// conflict's literal — a lower-bound literal credits the up side
    /// (HighsDomain.cpp:3833–3841, `increaseConflictScoreUp/Down`
    /// HighsPseudocost.h:147–159).
    pub fn add_conflict_observation(&mut self, col: usize, up: bool) {
        self.conflict_avg += self.conflict_weight;
        if up {
            self.conflict_up[col] += self.conflict_weight;
        } else {
            self.conflict_down[col] += self.conflict_weight;
        }
    }

    pub fn min_reliable(&self) -> u32 {
        self.min_reliable
    }

    pub fn set_min_reliable(&mut self, m: u32) {
        self.min_reliable = m;
    }

    /// Records that moving `col` by `delta` (> 0 up, < 0 down) raised the LP
    /// objective by `objdelta ≥ 0`.
    pub fn add_observation(&mut self, col: usize, delta: f64, objdelta: f64) {
        debug_assert!(delta != 0.0 && objdelta >= 0.0);
        let unit_gain = objdelta / delta.abs();
        if delta > 0.0 {
            push_mean(&mut self.cost_up[col], &mut self.samples_up[col], unit_gain);
        } else {
            push_mean(
                &mut self.cost_down[col],
                &mut self.samples_down[col],
                unit_gain,
            );
        }
        push_mean_total(&mut self.cost_total, &mut self.nsamples_total, unit_gain);
    }

    /// Records that the `up`/down branch on `col` was infeasible or cut off.
    pub fn add_cutoff_observation(&mut self, col: usize, up: bool) {
        self.ncutoffs_total += 1;
        if up {
            self.cutoffs_up[col] += 1;
        } else {
            self.cutoffs_down[col] += 1;
        }
    }

    /// Records how many bound changes propagation derived from the branch.
    pub fn add_inference_observation(&mut self, col: usize, inferences: usize, up: bool) {
        let v = inferences as f64;
        push_mean_total(&mut self.inferences_total, &mut self.ninferences_total, v);
        if up {
            push_mean(
                &mut self.inferences_up[col],
                &mut self.ninferences_up[col],
                v,
            );
        } else {
            push_mean(
                &mut self.inferences_down[col],
                &mut self.ninferences_down[col],
                v,
            );
        }
    }

    pub fn is_reliable(&self, col: usize) -> bool {
        self.samples_up[col].min(self.samples_down[col]) >= self.min_reliable
    }

    pub fn num_observations_up(&self, col: usize) -> u32 {
        self.samples_up[col]
    }

    pub fn num_observations_down(&self, col: usize) -> u32 {
        self.samples_down[col]
    }

    pub fn avg_pseudocost(&self) -> f64 {
        self.cost_total
    }

    pub fn avg_inferences_up(&self, col: usize) -> f64 {
        self.inferences_up[col]
    }

    pub fn avg_inferences_down(&self, col: usize) -> f64 {
        self.inferences_down[col]
    }

    /// Unreliable costs blend towards the global mean with weight
    /// `0.9 + 0.1·n/min_reliable` on the column's own mean (0 with no
    /// samples) — HighsPseudocost.h:257.
    fn blended(&self, own: f64, n: u32) -> f64 {
        if n == 0 || n < self.min_reliable {
            let w = if n == 0 {
                0.0
            } else {
                0.9 + 0.1 * f64::from(n) / f64::from(self.min_reliable)
            };
            w * own + (1.0 - w) * self.cost_total
        } else {
            own
        }
    }

    /// Estimated gain of rounding `frac` up, with `offset` added to the unit
    /// cost so zero pseudocosts still rank by fractionality.
    pub fn cost_up_blended(&self, col: usize, frac: f64, offset: f64) -> f64 {
        (frac.ceil() - frac) * (offset + self.blended(self.cost_up[col], self.samples_up[col]))
    }

    pub fn cost_down_blended(&self, col: usize, frac: f64, offset: f64) -> f64 {
        (frac - frac.floor()) * (offset + self.blended(self.cost_down[col], self.samples_down[col]))
    }

    /// Estimated gain of rounding up: the column's own mean, or the global
    /// mean when it has no samples (HighsPseudocost.h:290).
    pub fn cost_up(&self, col: usize, frac: f64) -> f64 {
        let unit = if self.samples_up[col] == 0 {
            self.cost_total
        } else {
            self.cost_up[col]
        };
        (frac.ceil() - frac) * unit
    }

    pub fn cost_down(&self, col: usize, frac: f64) -> f64 {
        let unit = if self.samples_down[col] == 0 {
            self.cost_total
        } else {
            self.cost_down[col]
        };
        (frac - frac.floor()) * unit
    }

    fn cutoff_rate(cutoffs: u32, samples: u32) -> f64 {
        f64::from(cutoffs) / (f64::from(cutoffs) + f64::from(samples)).max(1.0)
    }

    fn avg_cutoff_rate(&self) -> f64 {
        let c = self.ncutoffs_total as f64;
        c / (c + self.nsamples_total as f64).max(1.0)
    }

    /// The column's stored conflict scores discounted to the current weight,
    /// and the mean conflict weight per column — the normaliser of the
    /// conflict ratios (`conflictscore*/conflict_weight` and
    /// `conflict_avg_score/(conflict_weight·#columns)`, HighsPseudocost.h:335–339).
    fn conflict_sides(&self, col: usize) -> (f64, f64, f64) {
        let w = self.conflict_weight.max(MIN_THRESHOLD);
        (
            self.conflict_up[col] / w,
            self.conflict_down[col] / w,
            self.conflict_avg / (w * self.conflict_up.len() as f64),
        )
    }

    /// Branching score from estimated up/down gains: product rule, each part
    /// normalised by the squared global average (HighsPseudocost.h:310).
    pub fn score(&self, col: usize, up_gain: f64, down_gain: f64) -> f64 {
        let t = MIN_THRESHOLD;
        let cost = up_gain.max(t) * down_gain.max(t) / (self.cost_total * self.cost_total).max(t);
        let inference = self.inferences_up[col].max(t) * self.inferences_down[col].max(t)
            / (self.inferences_total * self.inferences_total).max(t);
        let avg_cut = self.avg_cutoff_rate();
        let cutoff = Self::cutoff_rate(self.cutoffs_up[col], self.samples_up[col]).max(t)
            * Self::cutoff_rate(self.cutoffs_down[col], self.samples_down[col]).max(t)
            / (avg_cut * avg_cut).max(t);
        let (cu, cd, avg) = self.conflict_sides(col);
        let conflict = cu.max(t) * cd.max(t) / (avg * avg).max(t);
        map_score(cost)
            + 1e-2 * map_score(conflict)
            + 1e-4 * (map_score(cutoff) + map_score(inference))
    }

    /// [`Self::score`] from the pseudocost estimates alone.
    pub fn score_frac(&self, col: usize, frac: f64) -> f64 {
        self.score(col, self.cost_up(col, frac), self.cost_down(col, frac))
    }

    /// One-sided score of the up branch (HighsPseudocost.h:359).
    pub fn score_up(&self, col: usize, frac: f64) -> f64 {
        let t = MIN_THRESHOLD;
        let cost = self.cost_up(col, frac) / self.cost_total.max(t);
        let inference = self.inferences_up[col] / self.inferences_total.max(t);
        let cutoff = Self::cutoff_rate(self.cutoffs_up[col], self.samples_up[col])
            / self.avg_cutoff_rate().max(t);
        let (cu, _, avg) = self.conflict_sides(col);
        map_score(cost)
            + 1e-2 * map_score(cu / avg.max(t))
            + 1e-4 * (map_score(cutoff) + map_score(inference))
    }

    /// One-sided score of the down branch (HighsPseudocost.h:389).
    pub fn score_down(&self, col: usize, frac: f64) -> f64 {
        let t = MIN_THRESHOLD;
        let cost = self.cost_down(col, frac) / self.cost_total.max(t);
        let inference = self.inferences_down[col] / self.inferences_total.max(t);
        let cutoff = Self::cutoff_rate(self.cutoffs_down[col], self.samples_down[col])
            / self.avg_cutoff_rate().max(t);
        let (_, cd, avg) = self.conflict_sides(col);
        map_score(cost)
            + 1e-2 * map_score(cd / avg.max(t))
            + 1e-4 * (map_score(cutoff) + map_score(inference))
    }

    /// Best-estimate of the integer optimum below a node: LP objective plus,
    /// per fractional column, the cheaper rounding direction's estimated gain
    /// (HiGHS `computeBestEstimate`, HighsLpRelaxation.cpp:428). The offset
    /// `feastol·max(|z|,1)/#integer columns` keeps fractionality visible when
    /// pseudocosts are zero.
    pub fn best_estimate(
        &self,
        objective: f64,
        fractional: &[(usize, f64)],
        feastol: f64,
        num_integer: usize,
    ) -> f64 {
        if fractional.is_empty() {
            return objective;
        }
        let offset = feastol * objective.abs().max(1.0) / num_integer.max(1) as f64;
        let increase: f64 = fractional
            .iter()
            .map(|&(j, v)| {
                self.cost_up_blended(j, v, offset)
                    .min(self.cost_down_blended(j, v, offset))
            })
            .sum();
        objective + increase
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn running_means_and_reliability() {
        let mut ps = Pseudocost::new(3);
        ps.set_min_reliable(2);
        ps.add_observation(0, 0.5, 1.0); // up, unit gain 2
        ps.add_observation(0, 0.25, 1.0); // up, unit gain 4
        ps.add_observation(0, -0.5, 0.5); // down, unit gain 1
        assert_eq!(ps.num_observations_up(0), 2);
        assert_eq!(ps.num_observations_down(0), 1);
        assert!(!ps.is_reliable(0));
        assert!((ps.cost_up(0, 2.5) - 0.5 * 3.0).abs() < 1e-12);
        assert!((ps.cost_down(0, 2.5) - 0.5 * 1.0).abs() < 1e-12);
        // Global mean of {2, 4, 1}.
        assert!((ps.avg_pseudocost() - 7.0 / 3.0).abs() < 1e-12);
        // Column 1 has no samples: falls back to the global mean.
        assert!((ps.cost_up(1, 0.75) - 0.25 * 7.0 / 3.0).abs() < 1e-12);
        ps.add_observation(0, -0.5, 1.5);
        assert!(ps.is_reliable(0));
    }

    #[test]
    fn blended_cost_weights_own_mean_by_sample_count() {
        let mut ps = Pseudocost::new(2);
        ps.set_min_reliable(4);
        ps.add_observation(0, 1.0, 10.0);
        ps.add_observation(1, 1.0, 2.0);
        // avg = 6; col 0 has 1 sample: w = 0.9 + 0.1/4 = 0.925.
        let w = 0.925;
        let expect = 0.5 * (0.1 + w * 10.0 + (1.0 - w) * 6.0);
        assert!((ps.cost_up_blended(0, 0.5, 0.1) - expect).abs() < 1e-12);
        // No samples: pure average.
        assert!((ps.cost_down_blended(0, 0.5, 0.0) - 0.5 * 6.0).abs() < 1e-12);
    }

    #[test]
    fn score_prefers_balanced_large_gains_and_cutoffs() {
        let mut ps = Pseudocost::new(3);
        for _ in 0..3 {
            ps.add_observation(0, 1.0, 4.0);
            ps.add_observation(0, -1.0, 4.0);
            ps.add_observation(1, 1.0, 8.0);
            ps.add_observation(1, -1.0, 0.0);
        }
        // Product rule: (4,4) beats (8,0).
        assert!(ps.score_frac(0, 0.5) > ps.score_frac(1, 0.5));
        let before = ps.score_frac(2, 0.5);
        ps.add_cutoff_observation(2, true);
        ps.add_cutoff_observation(2, false);
        assert!(ps.score_frac(2, 0.5) > before);
        // Inferences break ties among otherwise equal columns.
        let mut q = Pseudocost::new(2);
        q.add_inference_observation(0, 5, true);
        q.add_inference_observation(0, 5, false);
        q.add_inference_observation(1, 0, true);
        q.add_inference_observation(1, 0, false);
        assert!(q.score_frac(0, 0.5) > q.score_frac(1, 0.5));
    }

    #[test]
    fn best_estimate_adds_cheaper_direction() {
        let mut ps = Pseudocost::new(2);
        ps.set_min_reliable(1);
        ps.add_observation(0, 1.0, 1.0); // up unit 1
        ps.add_observation(0, -1.0, 3.0); // down unit 3
        // x0 = 0.4: up gain 0.6·1, down gain 0.4·3 ⇒ min 0.6 (offset 0).
        let e = ps.best_estimate(10.0, &[(0, 0.4)], 0.0, 2);
        assert!((e - 10.6).abs() < 1e-12, "{e}");
        assert_eq!(ps.best_estimate(10.0, &[], 1e-6, 2), 10.0);
    }

    #[test]
    fn conflict_term_has_higs_weight_and_normalisation() {
        let mut ps = Pseudocost::new(2);
        ps.add_observation(0, 1.0, 2.0);
        ps.add_observation(0, -1.0, 2.0);
        let before = ps.score_frac(0, 0.5);
        ps.increase_conflict_weight(); // weight = 1.02
        ps.add_conflict_observation(0, true);
        ps.add_conflict_observation(0, false);
        let after = ps.score_frac(0, 0.5);
        // Sides: 1.02/1.02 = 1 each; mean per column 2.04/(1.02·2) = 1;
        // ratio 1·1/1 = 1, mapScore = 0.5. With no conflicts the ratio is
        // the 1e-6 floor, so the gain is 1e-2·(0.5 − 1e-6) (getScore,
        // HighsPseudocost.h:335–349).
        assert!(
            (after - before - 1e-2 * (0.5 - 1e-6)).abs() < 1e-9,
            "{after}"
        );
        // One-sided (getScoreUp/Down :359–417): only the credited side sees
        // it, and with no conflicts at all the term is exactly zero there.
        let mut q = Pseudocost::new(2);
        q.add_observation(1, 1.0, 2.0);
        let up0 = q.score_up(1, 0.5);
        let down0 = q.score_down(1, 0.5);
        q.increase_conflict_weight();
        q.add_conflict_observation(1, true);
        // One literal over 2 columns: mean = 1.02/(1.02·2) = 0.5, ratio 2.
        assert!(
            (q.score_up(1, 0.5) - up0 - 1e-2 * (2.0 / 3.0)).abs() < 1e-9,
            "{}",
            q.score_up(1, 0.5)
        );
        assert_eq!(q.score_down(1, 0.5), down0);
    }

    #[test]
    fn conflict_scores_break_ties_between_equal_columns() {
        let mut ps = Pseudocost::new(2);
        for _ in 0..2 {
            ps.add_observation(0, 1.0, 2.0);
            ps.add_observation(0, -1.0, 2.0);
            ps.add_observation(1, 1.0, 2.0);
            ps.add_observation(1, -1.0, 2.0);
        }
        assert_eq!(ps.score_frac(0, 0.5), ps.score_frac(1, 0.5));
        ps.increase_conflict_weight();
        ps.add_conflict_observation(0, true);
        ps.add_conflict_observation(0, false);
        assert!(ps.score_frac(0, 0.5) > ps.score_frac(1, 0.5));
    }

    #[test]
    fn conflict_weight_rescales_above_1000() {
        let mut ps = Pseudocost::new(1);
        for _ in 0..400 {
            ps.increase_conflict_weight();
            ps.add_conflict_observation(0, true);
        }
        // 1.02^400 ≈ 2300 would overflow the cap without rescaling.
        assert!(ps.conflict_weight <= 1000.0);
        // The discount bounds the effective score: within an epoch the k-th
        // of m credits contributes 1.02^k/1.02^m, a geometric sum below
        // 1.02/0.02 = 51, and a rescaling divides past epochs by ~1000
        // (asserted with slack for the epoch in progress).
        let (cu, _, _) = ps.conflict_sides(0);
        assert!(cu < 60.0, "{cu}");
    }
}
