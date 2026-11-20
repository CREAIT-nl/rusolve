//! The LP relaxation as the separators see it: columns `0..n`, then one
//! *slack column* per LP row (`n + i` stands for the activity `a_iᵀx` with
//! bounds `[slack_lower_i, slack_upper_i]`), plus the sparse accumulator and
//! the row aggregator built on top of it.
//!
//! Modelled on HiGHS `HighsLpRelaxation::slackLower/slackUpper/isColIntegral`
//! (reference/HiGHS/highs/mip/HighsLpRelaxation.cpp:141–171, .h:168–227),
//! `HighsSparseVectorSum` (util/HighsSparseVectorSum.h:32–90) and
//! `HighsLpAggregator` (mip/HighsLpAggregator.cpp:12–56); HiGHS is MIT
//! licensed and was read as a reference text, not translated. Row
//! integrality (all columns integer, all coefficients integral ⇒ integral
//! slack, bounds rounded inward) follows HighsMipSolverData.cpp:949–974.

use crate::problem::MipProblem;

/// One LP row `lower ≤ Σ vals·x[inds] ≤ upper`, already in slack-bound form.
#[derive(Debug, Clone)]
pub(crate) struct LpRow {
    pub inds: Vec<usize>,
    pub vals: Vec<f64>,
    /// Slack bounds: the row bounds, an infinite side replaced by the
    /// activity bound over the global domain (may stay infinite).
    pub lower: f64,
    pub upper: f64,
    /// The LP row's own bounds (infinite sides kept): what decides whether a
    /// row is an equation or which side is tight.
    pub row_lower: f64,
    pub row_upper: f64,
    /// The slack only takes integral values on integer points.
    pub integral: bool,
}

/// Separation view of the current root LP.
pub(crate) struct Relaxation {
    pub num_col: usize,
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
    pub integer: Vec<bool>,
    /// Model rows first (`num_model_rows` of them), then cuts in the LP.
    pub rows: Vec<LpRow>,
    pub num_model_rows: usize,
    pub x: Vec<f64>,
    pub activity: Vec<f64>,
}

/// Minimum / maximum of `Σ vals·x` over the box (±∞ when unbounded).
pub(crate) fn activity_range(
    inds: &[usize],
    vals: &[f64],
    lower: &[f64],
    upper: &[f64],
) -> (f64, f64) {
    let (mut lo, mut hi) = (0.0, 0.0);
    for (&j, &v) in inds.iter().zip(vals) {
        if v > 0.0 {
            lo += v * lower[j];
            hi += v * upper[j];
        } else {
            lo += v * upper[j];
            hi += v * lower[j];
        }
    }
    // ∞ − ∞ cannot occur: each sum only accumulates one infinite sign.
    (lo, hi)
}

impl Relaxation {
    /// Model rows under the global (root) bounds `lower`/`upper`. Free rows
    /// (both sides infinite) are kept with infinite slack bounds; separators
    /// never use them.
    pub fn new(p: &MipProblem, lower: &[f64], upper: &[f64], epsilon: f64, feastol: f64) -> Self {
        let mut rows: Vec<LpRow> = (0..p.num_row)
            .map(|_| LpRow {
                inds: Vec::new(),
                vals: Vec::new(),
                lower: 0.0,
                upper: 0.0,
                row_lower: 0.0,
                row_upper: 0.0,
                integral: true,
            })
            .collect();
        let (ptr, idx, val) = (p.matrix.col_ptr(), p.matrix.row_idx(), p.matrix.values());
        for j in 0..p.num_col {
            for k in ptr[j]..ptr[j + 1] {
                if val[k] != 0.0 {
                    let r = &mut rows[idx[k] as usize];
                    r.inds.push(j);
                    r.vals.push(val[k]);
                }
            }
        }
        for (i, r) in rows.iter_mut().enumerate() {
            r.integral = r
                .inds
                .iter()
                .zip(&r.vals)
                .all(|(&j, &v)| p.integer[j] && (v - v.round()).abs() <= epsilon);
            let (mut lo, mut hi) = (p.row_lower[i], p.row_upper[i]);
            if r.integral {
                if lo.is_finite() {
                    lo = (lo - feastol).ceil();
                }
                if hi.is_finite() {
                    hi = (hi + feastol).floor();
                }
            }
            r.row_lower = lo;
            r.row_upper = hi;
            let (minact, maxact) = activity_range(&r.inds, &r.vals, lower, upper);
            r.lower = if lo.is_finite() { lo } else { minact };
            r.upper = if hi.is_finite() { hi } else { maxact };
        }
        Self {
            num_col: p.num_col,
            lower: lower.to_vec(),
            upper: upper.to_vec(),
            integer: p.integer.clone(),
            num_model_rows: rows.len(),
            rows,
            x: vec![0.0; p.num_col],
            activity: vec![0.0; p.num_row],
        }
    }

    /// Appends a cut `Σ vals·x ≤ rhs` that entered the LP: its slack lies in
    /// `[min activity, rhs]` (HighsLpRelaxation::slackLower, kCutPool case).
    pub fn push_cut(&mut self, inds: &[usize], vals: &[f64], rhs: f64, integral: bool) {
        let (minact, _) = activity_range(inds, vals, &self.lower, &self.upper);
        self.rows.push(LpRow {
            inds: inds.to_vec(),
            vals: vals.to_vec(),
            lower: minact,
            upper: rhs,
            row_lower: f64::NEG_INFINITY,
            row_upper: rhs,
            integral,
        });
    }

    /// Keeps the LP cut rows flagged in `keep` (indexed over cut rows only).
    pub fn retain_cuts(&mut self, keep: &[bool]) {
        let m0 = self.num_model_rows;
        let mut i = 0;
        self.rows.retain(|_| {
            i += 1;
            i - 1 < m0 || keep[i - 1 - m0]
        });
    }

    pub fn set_solution(&mut self, x: &[f64]) {
        self.x.clear();
        self.x.extend_from_slice(x);
        self.activity = self
            .rows
            .iter()
            .map(|r| r.inds.iter().zip(&r.vals).map(|(&j, &v)| v * x[j]).sum())
            .collect();
    }

    /// Size of the column + slack index space.
    pub fn num_ext(&self) -> usize {
        self.num_col + self.rows.len()
    }

    pub fn is_integral(&self, k: usize) -> bool {
        if k < self.num_col {
            self.integer[k]
        } else {
            self.rows[k - self.num_col].integral
        }
    }

    pub fn lb(&self, k: usize) -> f64 {
        if k < self.num_col {
            self.lower[k]
        } else {
            self.rows[k - self.num_col].lower
        }
    }

    pub fn ub(&self, k: usize) -> f64 {
        if k < self.num_col {
            self.upper[k]
        } else {
            self.rows[k - self.num_col].upper
        }
    }

    pub fn is_binary(&self, j: usize) -> bool {
        j < self.num_col && self.integer[j] && self.lower[j] == 0.0 && self.upper[j] == 1.0
    }
}

/// Dense-backed sparse accumulator (HighsSparseVectorSum): `add` keeps a
/// cancelled entry alive as a tiny value so the nonzero list stays exact.
pub(crate) struct SparseSum {
    values: Vec<f64>,
    nonzeros: Vec<usize>,
}

impl SparseSum {
    pub fn new(dim: usize) -> Self {
        Self {
            values: vec![0.0; dim],
            nonzeros: Vec::new(),
        }
    }

    pub fn ensure_dim(&mut self, dim: usize) {
        if self.values.len() < dim {
            self.values.resize(dim, 0.0);
        }
    }

    pub fn add(&mut self, k: usize, v: f64) {
        if self.values[k] != 0.0 {
            self.values[k] += v;
        } else {
            self.values[k] = v;
            self.nonzeros.push(k);
        }
        if self.values[k] == 0.0 {
            self.values[k] = f64::MIN_POSITIVE;
        }
    }

    /// Drops entries for which `is_zero(k, v)` holds and returns the rest,
    /// clearing the accumulator.
    pub fn take(&mut self, mut is_zero: impl FnMut(usize, f64) -> bool) -> (Vec<usize>, Vec<f64>) {
        let mut inds = Vec::with_capacity(self.nonzeros.len());
        let mut vals = Vec::with_capacity(self.nonzeros.len());
        for &k in &self.nonzeros {
            let v = self.values[k];
            self.values[k] = 0.0;
            if v.abs() <= f64::MIN_POSITIVE || is_zero(k, v) {
                continue;
            }
            inds.push(k);
            vals.push(v);
        }
        self.nonzeros.clear();
        (inds, vals)
    }

    pub fn clear(&mut self) {
        for &k in &self.nonzeros {
            self.values[k] = 0.0;
        }
        self.nonzeros.clear();
    }
}

/// Weighted row sums `Σ w_i (a_iᵀx − s_i) = 0` over columns and slacks
/// (HighsLpAggregator::addRow/getCurrentAggregation).
pub(crate) struct Aggregator {
    sum: SparseSum,
    /// Same entries as `sum`, kept so a snapshot does not clear the state.
    entries: Vec<(usize, f64)>,
}

impl Aggregator {
    pub fn new(dim: usize) -> Self {
        Self {
            sum: SparseSum::new(dim),
            entries: Vec::new(),
        }
    }

    pub fn add_row(&mut self, relax: &Relaxation, row: usize, weight: f64) {
        let r = &relax.rows[row];
        for (&j, &v) in r.inds.iter().zip(&r.vals) {
            self.entries.push((j, weight * v));
        }
        self.entries.push((relax.num_col + row, -weight));
    }

    /// The current aggregation as a base row `Σ vals·y ≤ 0` (or its negation),
    /// column entries below `droptol` dropped (slack entries are exact copies
    /// of the weights and are kept).
    pub fn current(
        &mut self,
        num_col: usize,
        droptol: f64,
        negate: bool,
    ) -> (Vec<usize>, Vec<f64>) {
        for &(k, v) in &self.entries {
            self.sum.add(k, v);
        }
        let (inds, mut vals) = self.sum.take(|k, v| k < num_col && v.abs() <= droptol);
        if negate {
            for v in &mut vals {
                *v = -*v;
            }
        }
        (inds, vals)
    }

    pub fn clear(&mut self) {
        self.entries.clear();
        self.sum.clear();
    }

    pub fn ensure_dim(&mut self, dim: usize) {
        self.sum.ensure_dim(dim);
    }
}
