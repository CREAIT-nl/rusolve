//! HiGHS source: mip/HighsPseudocost.h, mip/HighsPseudocost.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsPseudocost.h::minThreshold (line 27)
const MIN_THRESHOLD: f64 = 1e-6;

// PORT NOTE: C++ keeps a `HighsPseudocostDelta&` from markChanged(col) and updates it between
// writes to other members. Rust cannot hold that borrow across `&mut self`, so each body calls
// markChanged, then re-indexes `deltas[changedpos[col]]`, which is the same element.

/// PORT NOTE: the C++ lambda `mapScore` (HighsPseudocost.h lines 345, 382, 412) is shared by
/// getScore / getScoreUp / getScoreDown as one file-private function.
/// HiGHS: HighsPseudocost.h::HighsPseudocost::getScore (lines 345-345)
fn map_score(score: f64) -> f64 {
    1.0 - 1.0 / (1.0 + score)
}

/// HiGHS: HighsPseudocost.h::HighsPseudocostInitialization (lines 29-51)
#[derive(Clone, Default, Debug)]
pub struct HighsPseudocostInitialization {
    /// C++: `std::vector<double> pseudocostup`
    pub pseudocostup: Vec<f64>,
    /// C++: `std::vector<double> pseudocostdown`
    pub pseudocostdown: Vec<f64>,
    /// C++: `std::vector<HighsInt> nsamplesup`
    pub nsamplesup: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> nsamplesdown`
    pub nsamplesdown: Vec<HighsInt>,
    /// C++: `std::vector<double> inferencesup`
    pub inferencesup: Vec<f64>,
    /// C++: `std::vector<double> inferencesdown`
    pub inferencesdown: Vec<f64>,
    /// C++: `std::vector<HighsInt> ninferencesup`
    pub ninferencesup: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ninferencesdown`
    pub ninferencesdown: Vec<HighsInt>,
    /// C++: `std::vector<double> conflictscoreup`
    pub conflictscoreup: Vec<f64>,
    /// C++: `std::vector<double> conflictscoredown`
    pub conflictscoredown: Vec<f64>,
    /// C++: `double cost_total`
    pub cost_total: f64,
    /// C++: `double inferences_total`
    pub inferences_total: f64,
    /// C++: `double conflict_avg_score`
    pub conflict_avg_score: f64,
    /// C++: `int64_t nsamplestotal`
    pub nsamplestotal: i64,
    /// C++: `int64_t ninferencestotal`
    pub ninferencestotal: i64,
}

/// HiGHS: HighsPseudocost.h::HighsPseudocostDelta (lines 53-67)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HighsPseudocostDelta {
    /// C++: `HighsInt col = - 1`
    pub col: HighsInt,
    /// C++: `HighsInt nsamplesup = 0`
    pub nsamplesup: HighsInt,
    /// C++: `HighsInt nsamplesdown = 0`
    pub nsamplesdown: HighsInt,
    /// C++: `HighsInt ninferencesup = 0`
    pub ninferencesup: HighsInt,
    /// C++: `HighsInt ninferencesdown = 0`
    pub ninferencesdown: HighsInt,
    /// C++: `HighsInt ncutoffsup = 0`
    pub ncutoffsup: HighsInt,
    /// C++: `HighsInt ncutoffsdown = 0`
    pub ncutoffsdown: HighsInt,
    /// C++: `double pseudocostup_sum = 0.0`
    pub pseudocostup_sum: f64,
    /// C++: `double pseudocostdown_sum = 0.0`
    pub pseudocostdown_sum: f64,
    /// C++: `double inferencesup_sum = 0.0`
    pub inferencesup_sum: f64,
    /// C++: `double inferencesdown_sum = 0.0`
    pub inferencesdown_sum: f64,
    /// C++: `double conflictscoreup_sum = 0.0`
    pub conflictscoreup_sum: f64,
    /// C++: `double conflictscoredown_sum = 0.0`
    pub conflictscoredown_sum: f64,
}

impl Default for HighsPseudocostDelta {
    fn default() -> Self {
        Self {
            col: -1,
            nsamplesup: 0,
            nsamplesdown: 0,
            ninferencesup: 0,
            ninferencesdown: 0,
            ncutoffsup: 0,
            ncutoffsdown: 0,
            pseudocostup_sum: 0.0,
            pseudocostdown_sum: 0.0,
            inferencesup_sum: 0.0,
            inferencesdown_sum: 0.0,
            conflictscoreup_sum: 0.0,
            conflictscoredown_sum: 0.0,
        }
    }
}

/// HiGHS: HighsPseudocost.h::HighsPseudocost (lines 69-505)
#[derive(Clone, Default, Debug)]
pub struct HighsPseudocost {
    /// C++: `std::vector<double> pseudocostup`
    pub pseudocostup: Vec<f64>,
    /// C++: `std::vector<double> pseudocostdown`
    pub pseudocostdown: Vec<f64>,
    /// C++: `std::vector<HighsInt> nsamplesup`
    pub nsamplesup: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> nsamplesdown`
    pub nsamplesdown: Vec<HighsInt>,
    /// C++: `std::vector<double> inferencesup`
    pub inferencesup: Vec<f64>,
    /// C++: `std::vector<double> inferencesdown`
    pub inferencesdown: Vec<f64>,
    /// C++: `std::vector<HighsInt> ninferencesup`
    pub ninferencesup: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ninferencesdown`
    pub ninferencesdown: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ncutoffsup`
    pub ncutoffsup: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ncutoffsdown`
    pub ncutoffsdown: Vec<HighsInt>,
    /// C++: `std::vector<double> conflictscoreup`
    pub conflictscoreup: Vec<f64>,
    /// C++: `std::vector<double> conflictscoredown`
    pub conflictscoredown: Vec<f64>,
    /// C++: `std::vector<HighsInt> changedpos`
    pub changedpos: Vec<HighsInt>,
    /// C++: `std::vector<HighsPseudocostDelta> deltas`
    pub deltas: Vec<HighsPseudocostDelta>,
    /// C++: `double conflict_weight`
    pub conflict_weight: f64,
    /// C++: `double conflict_avg_score`
    pub conflict_avg_score: f64,
    /// C++: `double cost_total`
    pub cost_total: f64,
    /// C++: `double inferences_total`
    pub inferences_total: f64,
    /// C++: `double delta_cost_sum`
    pub delta_cost_sum: f64,
    /// C++: `double delta_inferences_sum`
    pub delta_inferences_sum: f64,
    /// C++: `int64_t nsamplestotal`
    pub nsamplestotal: i64,
    /// C++: `int64_t ninferencestotal`
    pub ninferencestotal: i64,
    /// C++: `int64_t ncutoffstotal`
    pub ncutoffstotal: i64,
    /// C++: `int64_t delta_nsamplestotal`
    pub delta_nsamplestotal: i64,
    /// C++: `int64_t delta_ninferencestotal`
    pub delta_ninferencestotal: i64,
    /// C++: `HighsInt minreliable`
    pub minreliable: HighsInt,
    /// C++: `double degeneracyFactor`
    pub degeneracy_factor: f64,
    /// rusolve `mip_pscost_variance_tol`: 0 = off (HiGHS count rule). Otherwise a pseudocost
    /// direction is reliable when the 95% half-width of its mean is within this fraction of the
    /// mean (Hendel 2015, sample variance of pseudo costs), or after 4x `minreliable` samples.
    pub variance_tol: f64,
    /// Welford count, mean and sum of squared deviations of the unit gains seen by this
    /// structure (local; not synced between workers or across restarts).
    pub wn_up: Vec<HighsInt>,
    pub wn_down: Vec<HighsInt>,
    pub wmean_up: Vec<f64>,
    pub wmean_down: Vec<f64>,
    pub wm2_up: Vec<f64>,
    pub wm2_down: Vec<f64>,
}

impl HighsPseudocost {
    /// HiGHS: HighsPseudocost.cpp::HighsPseudocost::HighsPseudocost (lines 12-64)
    /// C++: `HighsPseudocost(const HighsMipSolver & mipsolver)`
    pub fn new(mipsolver: &HighsMipSolver) -> Self {
        let num_col = mipsolver.num_col();
        let n = num_col as usize;
        let mut this = Self {
            pseudocostup: vec![0.0; n],
            pseudocostdown: vec![0.0; n],
            nsamplesup: vec![0; n],
            nsamplesdown: vec![0; n],
            inferencesup: vec![0.0; n],
            inferencesdown: vec![0.0; n],
            ninferencesup: vec![0; n],
            ninferencesdown: vec![0; n],
            ncutoffsup: vec![0; n],
            ncutoffsdown: vec![0; n],
            conflictscoreup: vec![0.0; n],
            conflictscoredown: vec![0.0; n],
            changedpos: vec![-1; n],
            deltas: Vec::with_capacity(cpp_min(256, num_col) as usize),
            conflict_weight: 1.0,
            conflict_avg_score: 0.0,
            cost_total: 0.0,
            inferences_total: 0.0,
            delta_cost_sum: 0.0,
            delta_inferences_sum: 0.0,
            nsamplestotal: 0,
            ninferencestotal: 0,
            ncutoffstotal: 0,
            delta_nsamplestotal: 0,
            delta_ninferencestotal: 0,
            minreliable: mipsolver.options_mip_.mip_pscost_minreliable,
            degeneracy_factor: 1.0,
            variance_tol: mipsolver.options_mip_.mip_pscost_variance_tol,
            wn_up: Vec::new(),
            wn_down: Vec::new(),
            wmean_up: Vec::new(),
            wmean_down: Vec::new(),
            wm2_up: Vec::new(),
            wm2_down: Vec::new(),
        };
        if this.variance_tol > 0.0 {
            this.wn_up = vec![0; n];
            this.wn_down = vec![0; n];
            this.wmean_up = vec![0.0; n];
            this.wmean_down = vec![0.0; n];
            this.wm2_up = vec![0.0; n];
            this.wm2_down = vec![0.0; n];
        }
        if let Some(pscostinit) = &mipsolver.pscostinit {
            this.cost_total = pscostinit.cost_total;
            this.inferences_total = pscostinit.inferences_total;
            this.nsamplestotal = pscostinit.nsamplestotal;
            this.ninferencestotal = pscostinit.ninferencestotal;

            this.conflict_avg_score = pscostinit.conflict_avg_score * (num_col as f64);

            // HiGHS #3250 (mip_presolve_row_append): HighsPseudocost.cpp line 49 (1.16.0)
            // loops over getOrigCols(), skipping the FM objective column.
            let cols: Vec<HighsInt> = if mipsolver.options_mip_.mip_presolve_row_append {
                mipsolver.mipdata_.post_solve_stack.get_orig_cols()
            } else {
                (0..num_col).collect()
            };
            for i in cols {
                let orig_col = mipsolver.mipdata_.post_solve_stack.get_orig_col_index(i) as usize;
                let i = i as usize;
                this.pseudocostup[i] = pscostinit.pseudocostup[orig_col];
                this.nsamplesup[i] = pscostinit.nsamplesup[orig_col];
                this.pseudocostdown[i] = pscostinit.pseudocostdown[orig_col];
                this.nsamplesdown[i] = pscostinit.nsamplesdown[orig_col];
                this.inferencesup[i] = pscostinit.inferencesup[orig_col];
                this.ninferencesup[i] = pscostinit.ninferencesup[orig_col];
                this.inferencesdown[i] = pscostinit.inferencesdown[orig_col];
                this.ninferencesdown[i] = pscostinit.ninferencesdown[orig_col];
                this.conflictscoreup[i] = pscostinit.conflictscoreup[orig_col];
                this.conflictscoredown[i] = pscostinit.conflictscoredown[orig_col];
            }
        }
        this
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::markChanged (lines 100-108)
    /// C++: `HighsPseudocostDelta & markChanged(HighsInt col)`
    pub fn mark_changed(&mut self, col: HighsInt) -> &mut HighsPseudocostDelta {
        debug_assert!(col >= 0 && (col as usize) < self.changedpos.len());
        if self.changedpos[col as usize] == -1 {
            self.changedpos[col as usize] = self.deltas.len() as HighsInt;
            let mut delta = HighsPseudocostDelta::default();
            delta.col = col;
            self.deltas.push(delta);
        }
        let pos = self.changedpos[col as usize] as usize;
        &mut self.deltas[pos]
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::addDeltaAverage (lines 110-117)
    /// C++: `static void addDeltaAverage()`
    pub fn add_delta_average_i32(
        avg: &mut f64,
        count: &mut HighsInt,
        delta_sum: f64,
        delta_count: HighsInt,
    ) {
        if delta_count <= 0 {
            return;
        }
        *avg = (*avg * (*count as f64) + delta_sum) / ((*count + delta_count) as f64);
        *count += delta_count;
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::addDeltaAverage (lines 110-117)
    /// C++: `static void addDeltaAverage()`
    pub fn add_delta_average_i64(avg: &mut f64, count: &mut i64, delta_sum: f64, delta_count: i64) {
        if delta_count <= 0 {
            return;
        }
        *avg = (*avg * (*count as f64) + delta_sum) / ((*count + delta_count) as f64);
        *count += delta_count;
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::increaseConflictWeight (lines 123-140)
    /// C++: `void increaseConflictWeight()`
    pub fn increase_conflict_weight(&mut self) {
        self.conflict_weight *= 1.02;

        if self.conflict_weight > 1000.0 {
            let scale = 1.0 / self.conflict_weight;
            self.conflict_weight = 1.0;
            self.conflict_avg_score *= scale;

            for i in 0..self.conflictscoreup.len() {
                self.conflictscoreup[i] *= scale;
                self.conflictscoredown[i] *= scale;
            }
            for delta in self.deltas.iter_mut() {
                delta.conflictscoreup_sum *= scale;
                delta.conflictscoredown_sum *= scale;
            }
        }
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::setDegeneracyFactor (lines 142-145)
    /// C++: `void setDegeneracyFactor(double degeneracyFactor)`
    pub fn set_degeneracy_factor(&mut self, degeneracy_factor: f64) {
        debug_assert!(degeneracy_factor >= 1.0);
        self.degeneracy_factor = degeneracy_factor;
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::increaseConflictScoreUp (lines 147-152)
    /// C++: `void increaseConflictScoreUp(HighsInt col)`
    pub fn increase_conflict_score_up(&mut self, col: HighsInt) {
        // PORT NOTE: C++ holds `delta` (from markChanged) across updates to other members.
        // Re-index `deltas[changedpos[col]]` after the update; it is the same element.
        self.mark_changed(col);
        let cc = col as usize;
        self.conflictscoreup[cc] += self.conflict_weight;
        self.conflict_avg_score += self.conflict_weight;
        let pos = self.changedpos[cc] as usize;
        self.deltas[pos].conflictscoreup_sum += self.conflict_weight;
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::increaseConflictScoreDown (lines 154-159)
    /// C++: `void increaseConflictScoreDown(HighsInt col)`
    pub fn increase_conflict_score_down(&mut self, col: HighsInt) {
        // PORT NOTE: see increase_conflict_score_up (delta is re-indexed after the update).
        self.mark_changed(col);
        let cc = col as usize;
        self.conflictscoredown[cc] += self.conflict_weight;
        self.conflict_avg_score += self.conflict_weight;
        let pos = self.changedpos[cc] as usize;
        self.deltas[pos].conflictscoredown_sum += self.conflict_weight;
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::setMinReliable (lines 161-161)
    /// C++: `void setMinReliable(HighsInt minreliable)`
    pub fn set_min_reliable(&mut self, minreliable: HighsInt) {
        self.minreliable = minreliable;
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getMinReliable (lines 163-163)
    /// C++: `HighsInt getMinReliable() const`
    pub fn get_min_reliable(&self) -> HighsInt {
        self.minreliable
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getNumObservations (lines 165-167)
    /// C++: `HighsInt getNumObservations(HighsInt col) const`
    pub fn get_num_observations(&self, col: HighsInt) -> HighsInt {
        let cc = col as usize;
        self.nsamplesup[cc] + self.nsamplesdown[cc]
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getNumObservationsUp (lines 169-169)
    /// C++: `HighsInt getNumObservationsUp(HighsInt col) const`
    pub fn get_num_observations_up(&self, col: HighsInt) -> HighsInt {
        self.nsamplesup[col as usize]
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getNumObservationsDown (lines 171-173)
    /// C++: `HighsInt getNumObservationsDown(HighsInt col) const`
    pub fn get_num_observations_down(&self, col: HighsInt) -> HighsInt {
        self.nsamplesdown[col as usize]
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::addCutoffObservation (lines 175-185)
    /// C++: `void addCutoffObservation(HighsInt col, bool upbranch)`
    pub fn add_cutoff_observation(&mut self, col: HighsInt, upbranch: bool) {
        self.mark_changed(col);
        let cc = col as usize;
        let pos = self.changedpos[cc] as usize;
        self.ncutoffstotal += 1;
        if upbranch {
            self.ncutoffsup[cc] += 1;
            self.deltas[pos].ncutoffsup += 1;
        } else {
            self.ncutoffsdown[cc] += 1;
            self.deltas[pos].ncutoffsdown += 1;
        }
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::addObservation (lines 187-218)
    /// C++: `void addObservation(HighsInt col, double delta, double objdelta)`
    pub fn add_observation(&mut self, col: HighsInt, delta: f64, objdelta: f64) {
        debug_assert!(delta != 0.0);
        debug_assert!(objdelta >= 0.0);
        self.mark_changed(col);
        let cc = col as usize;
        let pos = self.changedpos[cc] as usize;
        if delta > 0.0 {
            let unit_gain = objdelta / delta;
            if self.variance_tol > 0.0 {
                Self::welford(
                    &mut self.wn_up[cc],
                    &mut self.wmean_up[cc],
                    &mut self.wm2_up[cc],
                    unit_gain,
                );
            }
            let mut d = unit_gain - self.pseudocostup[cc];
            self.nsamplesup[cc] += 1;
            self.pseudocostup[cc] += d / (self.nsamplesup[cc] as f64);
            self.deltas[pos].nsamplesup += 1;
            self.deltas[pos].pseudocostup_sum += unit_gain;

            d = unit_gain - self.cost_total;
            self.nsamplestotal += 1;
            self.cost_total += d / (self.nsamplestotal as f64);
            self.delta_nsamplestotal += 1;
            self.delta_cost_sum += unit_gain;
        } else {
            let unit_gain = -objdelta / delta;
            if self.variance_tol > 0.0 {
                Self::welford(
                    &mut self.wn_down[cc],
                    &mut self.wmean_down[cc],
                    &mut self.wm2_down[cc],
                    unit_gain,
                );
            }
            let mut d = unit_gain - self.pseudocostdown[cc];
            self.nsamplesdown[cc] += 1;
            self.pseudocostdown[cc] += d / (self.nsamplesdown[cc] as f64);
            self.deltas[pos].nsamplesdown += 1;
            self.deltas[pos].pseudocostdown_sum += unit_gain;

            d = unit_gain - self.cost_total;
            self.nsamplestotal += 1;
            self.cost_total += d / (self.nsamplestotal as f64);
            self.delta_nsamplestotal += 1;
            self.delta_cost_sum += unit_gain;
        }
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::addInferenceObservation (lines 220-241)
    /// C++: `void addInferenceObservation(HighsInt col, HighsInt ninferences, bool upbranch)`
    pub fn add_inference_observation(
        &mut self,
        col: HighsInt,
        ninferences: HighsInt,
        upbranch: bool,
    ) {
        self.mark_changed(col);
        let cc = col as usize;
        let pos = self.changedpos[cc] as usize;
        let mut d = (ninferences as f64) - self.inferences_total;
        self.ninferencestotal += 1;
        self.inferences_total += d / (self.ninferencestotal as f64);
        self.delta_ninferencestotal += 1;
        self.delta_inferences_sum += ninferences as f64;
        if upbranch {
            d = (ninferences as f64) - self.inferencesup[cc];
            self.ninferencesup[cc] += 1;
            self.inferencesup[cc] += d / (self.ninferencesup[cc] as f64);
            self.deltas[pos].ninferencesup += 1;
            self.deltas[pos].inferencesup_sum += ninferences as f64;
        } else {
            d = (ninferences as f64) - self.inferencesdown[cc];
            self.ninferencesdown[cc] += 1;
            self.inferencesdown[cc] += d / (self.ninferencesdown[cc] as f64);
            self.deltas[pos].ninferencesdown += 1;
            self.deltas[pos].inferencesdown_sum += ninferences as f64;
        }
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::isReliable (lines 243-245)
    /// C++: `bool isReliable(HighsInt col) const`
    pub fn is_reliable(&self, col: HighsInt) -> bool {
        let cc = col as usize;
        if self.variance_tol > 0.0 && self.minreliable > 0 {
            return self.variance_reliable(cc, true) && self.variance_reliable(cc, false);
        }
        cpp_min(self.nsamplesup[cc], self.nsamplesdown[cc]) >= self.minreliable
    }

    /// Welford update of (count, mean, M2).
    fn welford(n: &mut HighsInt, mean: &mut f64, m2: &mut f64, x: f64) {
        *n += 1;
        let d = x - *mean;
        *mean += d / (*n as f64);
        *m2 += d * (x - *mean);
    }

    /// rusolve `mip_pscost_variance_tol` (Hendel 2015): one direction is reliable after
    /// 4x `minreliable` samples, or when it has at least two local samples and the 95% half-width
    /// of the mean of its unit gains is within `variance_tol` of the mean. With fewer than two
    /// local samples (restart, worker sync) the HiGHS count rule decides.
    fn variance_reliable(&self, cc: usize, up: bool) -> bool {
        let (ns, wn, mean, m2) = if up {
            (
                self.nsamplesup[cc],
                self.wn_up[cc],
                self.wmean_up[cc],
                self.wm2_up[cc],
            )
        } else {
            (
                self.nsamplesdown[cc],
                self.wn_down[cc],
                self.wmean_down[cc],
                self.wm2_down[cc],
            )
        };
        if ns >= self.minreliable.saturating_mul(4) {
            return true;
        }
        if wn < 2 {
            return ns >= self.minreliable;
        }
        let half = 1.96 * (m2 / ((wn - 1) as f64) / (wn as f64)).sqrt();
        half <= self.variance_tol * mean.abs()
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::isReliableUp (lines 247-249)
    /// C++: `bool isReliableUp(HighsInt col) const`
    pub fn is_reliable_up(&self, col: HighsInt) -> bool {
        self.nsamplesup[col as usize] >= self.minreliable
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::isReliableDown (lines 251-253)
    /// C++: `bool isReliableDown(HighsInt col) const`
    pub fn is_reliable_down(&self, col: HighsInt) -> bool {
        self.nsamplesdown[col as usize] >= self.minreliable
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getAvgPseudocost (lines 255-255)
    /// C++: `double getAvgPseudocost() const`
    pub fn get_avg_pseudocost(&self) -> f64 {
        self.cost_total
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getPseudocostUp (lines 257-271)
    /// C++: `double getPseudocostUp(HighsInt col, double frac, double offset) const`
    pub fn get_pseudocost_up(&self, col: HighsInt, frac: f64, offset: f64) -> f64 {
        let cc = col as usize;
        let up = frac.ceil() - frac;

        let cost = if self.nsamplesup[cc] == 0 || self.nsamplesup[cc] < self.minreliable {
            let weight_ps = if self.nsamplesup[cc] == 0 {
                0.0
            } else {
                0.9 + 0.1 * (self.nsamplesup[cc] as f64) / (self.minreliable as f64)
            };
            let mut cost = weight_ps * self.pseudocostup[cc];
            cost += (1.0 - weight_ps) * self.get_avg_pseudocost();
            cost
        } else {
            self.pseudocostup[cc]
        };
        up * (offset + cost)
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getPseudocostDown (lines 273-288)
    /// C++: `double getPseudocostDown(HighsInt col, double frac, double offset) const`
    pub fn get_pseudocost_down(&self, col: HighsInt, frac: f64, offset: f64) -> f64 {
        let cc = col as usize;
        let down = frac - frac.floor();

        let cost = if self.nsamplesdown[cc] == 0 || self.nsamplesdown[cc] < self.minreliable {
            let weight_ps = if self.nsamplesdown[cc] == 0 {
                0.0
            } else {
                0.9 + 0.1 * (self.nsamplesdown[cc] as f64) / (self.minreliable as f64)
            };
            let mut cost = weight_ps * self.pseudocostdown[cc];
            cost += (1.0 - weight_ps) * self.get_avg_pseudocost();
            cost
        } else {
            self.pseudocostdown[cc]
        };
        down * (offset + cost)
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getPseudocostUp (lines 290-294)
    /// C++: `double getPseudocostUp(HighsInt col, double frac) const`
    pub fn get_pseudocost_up_2(&self, col: HighsInt, frac: f64) -> f64 {
        let cc = col as usize;
        let up = frac.ceil() - frac;
        if self.nsamplesup[cc] == 0 {
            return up * self.cost_total;
        }
        up * self.pseudocostup[cc]
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getPseudocostDown (lines 296-300)
    /// C++: `double getPseudocostDown(HighsInt col, double frac) const`
    pub fn get_pseudocost_down_2(&self, col: HighsInt, frac: f64) -> f64 {
        let cc = col as usize;
        let down = frac - frac.floor();
        if self.nsamplesdown[cc] == 0 {
            return down * self.cost_total;
        }
        down * self.pseudocostdown[cc]
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getConflictScoreUp (lines 302-304)
    /// C++: `double getConflictScoreUp(HighsInt col) const`
    pub fn get_conflict_score_up(&self, col: HighsInt) -> f64 {
        self.conflictscoreup[col as usize] / self.conflict_weight
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getConflictScoreDown (lines 306-308)
    /// C++: `double getConflictScoreDown(HighsInt col) const`
    pub fn get_conflict_score_down(&self, col: HighsInt) -> f64 {
        self.conflictscoredown[col as usize] / self.conflict_weight
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getScore (lines 310-350)
    /// C++: `double getScore(HighsInt col, double upcost, double downcost) const`
    pub fn get_score(&self, col: HighsInt, upcost: f64, downcost: f64) -> f64 {
        let cc = col as usize;
        let cost_score = cpp_max(upcost, MIN_THRESHOLD) * cpp_max(downcost, MIN_THRESHOLD)
            / cpp_max(MIN_THRESHOLD, self.cost_total * self.cost_total);
        let inference_score = cpp_max(self.inferencesup[cc], MIN_THRESHOLD)
            * cpp_max(self.inferencesdown[cc], MIN_THRESHOLD)
            / cpp_max(MIN_THRESHOLD, self.inferences_total * self.inferences_total);

        let cut_off_score_up = (self.ncutoffsup[cc] as f64)
            / cpp_max(
                1.0,
                (self.ncutoffsup[cc] as f64) + (self.nsamplesup[cc] as f64),
            );
        let cut_off_score_down = (self.ncutoffsdown[cc] as f64)
            / cpp_max(
                1.0,
                (self.ncutoffsdown[cc] as f64) + (self.nsamplesdown[cc] as f64),
            );
        let avg_cutoffs = (self.ncutoffstotal as f64)
            / cpp_max(
                1.0,
                (self.ncutoffstotal as f64) + (self.nsamplestotal as f64),
            );

        let cutoff_score = cpp_max(cut_off_score_up, MIN_THRESHOLD)
            * cpp_max(cut_off_score_down, MIN_THRESHOLD)
            / cpp_max(MIN_THRESHOLD, avg_cutoffs * avg_cutoffs);

        let conflict_score_up = self.conflictscoreup[cc] / self.conflict_weight;
        let conflict_score_down = self.conflictscoredown[cc] / self.conflict_weight;
        let conflict_score_avg =
            self.conflict_avg_score / (self.conflict_weight * (self.conflictscoreup.len() as f64));
        let conflict_score = cpp_max(conflict_score_up, MIN_THRESHOLD)
            * cpp_max(conflict_score_down, MIN_THRESHOLD)
            / cpp_max(MIN_THRESHOLD, conflict_score_avg * conflict_score_avg);

        map_score(cost_score) / self.degeneracy_factor
            + self.degeneracy_factor
                * (1e-2 * map_score(conflict_score)
                    + 1e-4 * (map_score(cutoff_score) + map_score(inference_score)))
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getScore (lines 352-357)
    /// C++: `double getScore(HighsInt col, double frac) const`
    pub fn get_score_2(&self, col: HighsInt, frac: f64) -> f64 {
        let upcost = self.get_pseudocost_up_2(col, frac);
        let downcost = self.get_pseudocost_down_2(col, frac);

        self.get_score(col, upcost, downcost)
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getScoreUp (lines 359-387)
    /// C++: `double getScoreUp(HighsInt col, double frac) const`
    pub fn get_score_up(&self, col: HighsInt, frac: f64) -> f64 {
        let cc = col as usize;
        let cost_score =
            self.get_pseudocost_up_2(col, frac) / cpp_max(MIN_THRESHOLD, self.cost_total);
        let inference_score = self.inferencesup[cc] / cpp_max(MIN_THRESHOLD, self.inferences_total);

        let cut_off_score_up = (self.ncutoffsup[cc] as f64)
            / cpp_max(
                1.0,
                (self.ncutoffsup[cc] as f64) + (self.nsamplesup[cc] as f64),
            );
        let avg_cutoffs = (self.ncutoffstotal as f64)
            / cpp_max(
                1.0,
                (self.ncutoffstotal as f64) + (self.nsamplestotal as f64),
            );

        let cutoff_score = cut_off_score_up / cpp_max(MIN_THRESHOLD, avg_cutoffs);

        let conflict_score_up = self.conflictscoreup[cc] / self.conflict_weight;
        let conflict_score_avg =
            self.conflict_avg_score / (self.conflict_weight * (self.conflictscoreup.len() as f64));
        let conflict_score = conflict_score_up / cpp_max(MIN_THRESHOLD, conflict_score_avg);

        map_score(cost_score)
            + (1e-2 * map_score(conflict_score)
                + 1e-4 * (map_score(cutoff_score) + map_score(inference_score)))
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getScoreDown (lines 389-417)
    /// C++: `double getScoreDown(HighsInt col, double frac) const`
    pub fn get_score_down(&self, col: HighsInt, frac: f64) -> f64 {
        let cc = col as usize;
        let cost_score =
            self.get_pseudocost_down_2(col, frac) / cpp_max(MIN_THRESHOLD, self.cost_total);
        let inference_score =
            self.inferencesdown[cc] / cpp_max(MIN_THRESHOLD, self.inferences_total);

        let cut_off_score_down = (self.ncutoffsdown[cc] as f64)
            / cpp_max(
                1.0,
                (self.ncutoffsdown[cc] as f64) + (self.nsamplesdown[cc] as f64),
            );
        let avg_cutoffs = (self.ncutoffstotal as f64)
            / cpp_max(
                1.0,
                (self.ncutoffstotal as f64) + (self.nsamplestotal as f64),
            );

        let cutoff_score = cut_off_score_down / cpp_max(MIN_THRESHOLD, avg_cutoffs);

        let conflict_score_down = self.conflictscoredown[cc] / self.conflict_weight;
        let conflict_score_avg = self.conflict_avg_score
            / (self.conflict_weight * (self.conflictscoredown.len() as f64));
        let conflict_score = conflict_score_down / cpp_max(MIN_THRESHOLD, conflict_score_avg);

        map_score(cost_score)
            + (1e-2 * map_score(conflict_score)
                + 1e-4 * (map_score(cutoff_score) + map_score(inference_score)))
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getAvgInferencesUp (lines 419-419)
    /// C++: `double getAvgInferencesUp(HighsInt col) const`
    pub fn get_avg_inferences_up(&self, col: HighsInt) -> f64 {
        self.inferencesup[col as usize]
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::getAvgInferencesDown (lines 421-423)
    /// C++: `double getAvgInferencesDown(HighsInt col) const`
    pub fn get_avg_inferences_down(&self, col: HighsInt) -> f64 {
        self.inferencesdown[col as usize]
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::flushPseudoCost (lines 425-458)
    /// C++: `void flushPseudoCost(HighsPseudocost & pseudocost)`
    pub fn flush_pseudo_cost(&mut self, pseudocost: &mut HighsPseudocost) {
        debug_assert_eq!(pseudocost.ncutoffsup.len(), self.ncutoffsup.len());
        for delta in pseudocost.deltas.iter() {
            let col = delta.col;
            debug_assert!(col >= 0 && (col as usize) < pseudocost.ncutoffsup.len());
            let cc = col as usize;
            Self::add_delta_average_i32(
                &mut self.pseudocostup[cc],
                &mut self.nsamplesup[cc],
                delta.pseudocostup_sum,
                delta.nsamplesup,
            );
            Self::add_delta_average_i32(
                &mut self.pseudocostdown[cc],
                &mut self.nsamplesdown[cc],
                delta.pseudocostdown_sum,
                delta.nsamplesdown,
            );
            Self::add_delta_average_i32(
                &mut self.inferencesup[cc],
                &mut self.ninferencesup[cc],
                delta.inferencesup_sum,
                delta.ninferencesup,
            );
            Self::add_delta_average_i32(
                &mut self.inferencesdown[cc],
                &mut self.ninferencesdown[cc],
                delta.inferencesdown_sum,
                delta.ninferencesdown,
            );
            let conflict_scale = self.conflict_weight / pseudocost.conflict_weight;
            let conflict_delta_up = conflict_scale * delta.conflictscoreup_sum;
            let conflict_delta_down = conflict_scale * delta.conflictscoredown_sum;
            self.conflictscoreup[cc] += conflict_delta_up;
            self.conflictscoredown[cc] += conflict_delta_down;
            self.conflict_avg_score += conflict_delta_up + conflict_delta_down;
            self.ncutoffsup[cc] += delta.ncutoffsup;
            self.ncutoffsdown[cc] += delta.ncutoffsdown;
            self.ncutoffstotal += (delta.ncutoffsup + delta.ncutoffsdown) as i64;
        }
        Self::add_delta_average_i64(
            &mut self.cost_total,
            &mut self.nsamplestotal,
            pseudocost.delta_cost_sum,
            pseudocost.delta_nsamplestotal,
        );
        Self::add_delta_average_i64(
            &mut self.inferences_total,
            &mut self.ninferencestotal,
            pseudocost.delta_inferences_sum,
            pseudocost.delta_ninferencestotal,
        );
        pseudocost.remove_changed();
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::syncPseudoCost (lines 460-493)
    /// C++: `void syncPseudoCost(HighsPseudocost & pseudocost)`
    pub fn sync_pseudo_cost(&mut self, pseudocost: &mut HighsPseudocost) {
        // PORT NOTE: std::copy(src.begin(), src.end(), dst.begin()) copies src.len() elements into
        // the front of dst; the slice prefix below is that copy.
        pseudocost.pseudocostup[..self.pseudocostup.len()].copy_from_slice(&self.pseudocostup);
        pseudocost.pseudocostdown[..self.pseudocostdown.len()]
            .copy_from_slice(&self.pseudocostdown);
        pseudocost.nsamplesup[..self.nsamplesup.len()].copy_from_slice(&self.nsamplesup);
        pseudocost.nsamplesdown[..self.nsamplesdown.len()].copy_from_slice(&self.nsamplesdown);
        pseudocost.inferencesup[..self.inferencesup.len()].copy_from_slice(&self.inferencesup);
        pseudocost.inferencesdown[..self.inferencesdown.len()]
            .copy_from_slice(&self.inferencesdown);
        pseudocost.ninferencesup[..self.ninferencesup.len()].copy_from_slice(&self.ninferencesup);
        pseudocost.ninferencesdown[..self.ninferencesdown.len()]
            .copy_from_slice(&self.ninferencesdown);
        pseudocost.ncutoffsup[..self.ncutoffsup.len()].copy_from_slice(&self.ncutoffsup);
        pseudocost.ncutoffsdown[..self.ncutoffsdown.len()].copy_from_slice(&self.ncutoffsdown);
        pseudocost.conflictscoreup[..self.conflictscoreup.len()]
            .copy_from_slice(&self.conflictscoreup);
        pseudocost.conflictscoredown[..self.conflictscoredown.len()]
            .copy_from_slice(&self.conflictscoredown);
        pseudocost.conflict_weight = self.conflict_weight;
        pseudocost.conflict_avg_score = self.conflict_avg_score;
        pseudocost.cost_total = self.cost_total;
        pseudocost.inferences_total = self.inferences_total;
        pseudocost.nsamplestotal = self.nsamplestotal;
        pseudocost.ninferencestotal = self.ninferencestotal;
        pseudocost.ncutoffstotal = self.ncutoffstotal;
        pseudocost.remove_changed();
    }

    /// HiGHS: HighsPseudocost.h::HighsPseudocost::removeChanged (lines 495-504)
    /// C++: `void removeChanged()`
    pub fn remove_changed(&mut self) {
        for delta in self.deltas.iter() {
            self.changedpos[delta.col as usize] = -1;
        }
        self.deltas.clear();
        self.delta_cost_sum = 0.0;
        self.delta_inferences_sum = 0.0;
        self.delta_nsamplestotal = 0;
        self.delta_ninferencestotal = 0;
    }
}

impl HighsPseudocostInitialization {
    /// HiGHS: HighsPseudocost.cpp::HighsPseudocostInitialization::HighsPseudocostInitialization (lines 66-93)
    /// C++: `HighsPseudocostInitialization(const HighsPseudocost & pscost, HighsInt maxCount)`
    pub fn new(pscost: &HighsPseudocost, max_count: HighsInt) -> Self {
        let mut this = Self {
            pseudocostup: pscost.pseudocostup.clone(),
            pseudocostdown: pscost.pseudocostdown.clone(),
            nsamplesup: pscost.nsamplesup.clone(),
            nsamplesdown: pscost.nsamplesdown.clone(),
            inferencesup: pscost.inferencesup.clone(),
            inferencesdown: pscost.inferencesdown.clone(),
            ninferencesup: pscost.ninferencesup.clone(),
            ninferencesdown: pscost.ninferencesdown.clone(),
            conflictscoreup: vec![0.0; pscost.conflictscoreup.len()],
            conflictscoredown: vec![0.0; pscost.conflictscoreup.len()],
            cost_total: pscost.cost_total,
            inferences_total: pscost.inferences_total,
            conflict_avg_score: pscost.conflict_avg_score,
            nsamplestotal: cpp_min(1, pscost.nsamplestotal),
            ninferencestotal: cpp_min(1, pscost.ninferencestotal),
        };
        let ncol = this.pseudocostup.len() as HighsInt;
        this.conflict_avg_score /= (ncol as f64) * pscost.conflict_weight;
        for i in 0..ncol as usize {
            this.nsamplesup[i] = cpp_min(this.nsamplesup[i], max_count);
            this.nsamplesdown[i] = cpp_min(this.nsamplesdown[i], max_count);
            this.ninferencesup[i] = cpp_min(this.ninferencesup[i], 1);
            this.ninferencesdown[i] = cpp_min(this.ninferencesdown[i], 1);
            this.conflictscoreup[i] = pscost.conflictscoreup[i] / pscost.conflict_weight;
            this.conflictscoredown[i] = pscost.conflictscoredown[i] / pscost.conflict_weight;
        }
        this
    }

    /// HiGHS: HighsPseudocost.cpp::HighsPseudocostInitialization::HighsPseudocostInitialization (lines 95-135)
    /// C++: `HighsPseudocostInitialization(const HighsPseudocost & pscost, HighsInt maxCount, const presolve::HighsPostsolveStack & postsolveStack)`
    ///
    /// `x1b_orig_cols` (`mip_presolve_row_append`, #3250): loop over getOrigCols()
    /// (1.16.0 line 117) instead of every column.
    pub fn new_2(
        pscost: &HighsPseudocost,
        max_count: HighsInt,
        postsolve_stack: &HighsPostsolveStack,
        x1b_orig_cols: bool,
    ) -> Self {
        let n = postsolve_stack.get_orig_num_col() as usize;
        let mut this = Self {
            pseudocostup: vec![0.0; n],
            pseudocostdown: vec![0.0; n],
            nsamplesup: vec![0; n],
            nsamplesdown: vec![0; n],
            inferencesup: vec![0.0; n],
            inferencesdown: vec![0.0; n],
            ninferencesup: vec![0; n],
            ninferencesdown: vec![0; n],
            conflictscoreup: vec![0.0; n],
            conflictscoredown: vec![0.0; n],
            cost_total: pscost.cost_total,
            inferences_total: pscost.inferences_total,
            conflict_avg_score: pscost.conflict_avg_score,
            nsamplestotal: cpp_min(1, pscost.nsamplestotal),
            ninferencestotal: cpp_min(1, pscost.ninferencestotal),
        };
        let ncols = pscost.pseudocostup.len() as HighsInt;
        this.conflict_avg_score /= (ncols as f64) * pscost.conflict_weight;
        let cols: Vec<HighsInt> = if x1b_orig_cols {
            postsolve_stack.get_orig_cols()
        } else {
            (0..ncols).collect()
        };
        for i in cols {
            let oc = postsolve_stack.get_orig_col_index(i) as usize;
            let i = i as usize;
            this.pseudocostup[oc] = pscost.pseudocostup[i];
            this.pseudocostdown[oc] = pscost.pseudocostdown[i];
            this.nsamplesup[oc] = cpp_min(max_count, pscost.nsamplesup[i]);
            this.nsamplesdown[oc] = cpp_min(max_count, pscost.nsamplesdown[i]);
            this.inferencesup[oc] = pscost.inferencesup[i];
            this.inferencesdown[oc] = pscost.inferencesdown[i];
            this.ninferencesup[oc] = 1;
            this.ninferencesdown[oc] = 1;
            this.conflictscoreup[oc] = pscost.conflictscoreup[i] / pscost.conflict_weight;
            this.conflictscoredown[oc] = pscost.conflictscoredown[i] / pscost.conflict_weight;
        }
        this
    }
}

#[cfg(test)]
mod x1b_psr_tests {
    use super::*;

    fn pscost(ncols: usize) -> HighsPseudocost {
        let vals: Vec<f64> = (0..ncols).map(|i| (i + 1) as f64).collect();
        HighsPseudocost {
            pseudocostup: vals.clone(),
            pseudocostdown: vals.clone(),
            nsamplesup: vec![1; ncols],
            nsamplesdown: vec![1; ncols],
            inferencesup: vals.clone(),
            inferencesdown: vals.clone(),
            conflictscoreup: vals.clone(),
            conflictscoredown: vals,
            conflict_weight: 1.0,
            ..Default::default()
        }
    }

    /// HiGHS backport-PSR-E: with `mip_presolve_row_append` the restart pseudocost
    /// initialisation reads only the original columns (HighsPseudocost.cpp 1.16.0), so a
    /// column appended by presolve (FM objective column) is skipped; on a stack without
    /// appended columns both paths agree.
    #[test]
    fn x1b_psr_e_pseudocost_skips_appended() {
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(0, 2);
        let plain_on = HighsPseudocostInitialization::new_2(&pscost(2), 10, &stack, true);
        let plain_off = HighsPseudocostInitialization::new_2(&pscost(2), 10, &stack, false);
        assert_eq!(plain_on.pseudocostup, plain_off.pseudocostup);
        assert_eq!(plain_on.conflictscoreup, plain_off.conflictscoreup);

        stack.append_col_to_model();
        let init = HighsPseudocostInitialization::new_2(&pscost(3), 10, &stack, true);
        assert_eq!(init.pseudocostup, vec![1.0, 2.0]);
        assert_eq!(init.inferencesdown, vec![1.0, 2.0]);
        assert_eq!(init.ninferencesup, vec![1, 1]);
    }
}
