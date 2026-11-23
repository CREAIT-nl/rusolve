//! Variable bounds `x ≤ c·y + d` / `x ≥ c·y + d` (y binary) read off single
//! model rows, and the *implied bound cuts* they give.
//!
//! Extraction: each side of a row is brought to `Σ v_k z_k ≤ r` with
//! `v_k > 0`, `z_k ∈ [0, u_k − l_k]` (z = x − l, or z = u − x for a negative
//! coefficient). Setting a binary's z to 1 leaves `r − v_b` for a non-binary
//! z_j, which becomes a variable bound on x_j when it is tighter than the
//! column's range — HiGHS `HighsCliqueTable::extractCliques`
//! (reference/HiGHS/highs/mip/HighsCliqueTable.cpp:839–900, 1268–1378, MIT;
//! reimplemented). Bounds whose extreme value overshoots the column bound are
//! tightened as in `HighsImplications::cleanupVub/cleanupVlb`
//! (HighsImplications.cpp:795–875) — the bound substitution in
//! [`super::transform`] relies on `u − l` bounding the substituted slack.
//!
//! Implied bound cuts: a violated `x ≤ c·y + d` is added as the cut
//! `x − c·y ≤ d`. HiGHS separates the same inequality from probing
//! implications (`HighsImplications::separateImpliedBounds`,
//! HighsImplications.cpp:583–700); we take the implications from single rows
//! only — this catches the weak facility-location linking (`Σ_i x_ij ≤ n·y_j`
//! ⇒ `x_ij ≤ y_j`), which closes 97–100% of the root gap on our
//! facility-location instances. Literature: Wolsey, *Integer Programming* (1998) §8.9.
// TODO(perf): probing-based implications (HighsImplications::runProbing)
// subsume these and also cover general integers.

use super::relax::Relaxation;

/// `coef·y + constant` with y binary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct VarBound {
    pub coef: f64,
    pub constant: f64,
}

impl VarBound {
    pub fn min_value(&self) -> f64 {
        self.constant + self.coef.min(0.0)
    }
    pub fn max_value(&self) -> f64 {
        self.constant + self.coef.max(0.0)
    }
    pub fn at(&self, y: f64) -> f64 {
        self.coef * y + self.constant
    }
}

/// Per column: variable upper / lower bounds `(binary column, bound)`.
pub(crate) struct VarBounds {
    pub vubs: Vec<Vec<(usize, VarBound)>>,
    pub vlbs: Vec<Vec<(usize, VarBound)>>,
}

/// Caps that keep extraction linear-ish on dense rows (HiGHS caps the
/// implication store with `tooManyVarBounds`).
const MAX_PER_COL: usize = 16;
const MAX_PAIRS_PER_ROW: usize = 2000;

impl VarBounds {
    pub fn empty(n: usize) -> Self {
        Self {
            vubs: vec![Vec::new(); n],
            vlbs: vec![Vec::new(); n],
        }
    }

    /// Extracts variable bounds from every model row of `relax` (both sides).
    pub fn extract(relax: &Relaxation, feastol: f64, epsilon: f64) -> Self {
        let mut vb = Self::empty(relax.num_col);
        let mut z: Vec<(usize, f64, bool)> = Vec::new();
        for i in 0..relax.num_model_rows {
            let row = &relax.rows[i];
            if !row.inds.iter().any(|&j| relax.is_binary(j)) {
                continue;
            }
            for (upper_side, dir) in [(true, 1.0), (false, -1.0)] {
                let raw = relax_row_bound(relax, i, upper_side);
                if !raw.is_finite() {
                    continue;
                }
                let mut rhs = dir * raw;
                z.clear();
                let mut free = false;
                for (&j, &a) in row.inds.iter().zip(&row.vals) {
                    let v = dir * a;
                    if v.abs() < epsilon {
                        continue;
                    }
                    if v < 0.0 {
                        if !relax.upper[j].is_finite() {
                            free = true;
                            break;
                        }
                        rhs -= v * relax.upper[j];
                        z.push((j, -v, true));
                    } else {
                        if !relax.lower[j].is_finite() {
                            free = true;
                            break;
                        }
                        rhs -= v * relax.lower[j];
                        z.push((j, v, false));
                    }
                }
                if free {
                    continue;
                }
                vb.extract_side(relax, &z, rhs, feastol, epsilon);
            }
        }
        vb
    }

    /// `Σ v z ≤ rhs` with every `v > 0`, `z ≥ 0`; `comp` marks z = u − x.
    fn extract_side(
        &mut self,
        relax: &Relaxation,
        z: &[(usize, f64, bool)],
        rhs: f64,
        feastol: f64,
        epsilon: f64,
    ) {
        // Only non-binary, non-fixed partners yield pairs: without one the
        // scan below is |z|² `continue`s for nothing:
        // nw04's all-binary 17k-entry rows
        // spent 80 s here, past a 60 s limit. Same output, O(|z|).
        if !z
            .iter()
            .any(|&(j, _, _)| !relax.is_binary(j) && relax.lower[j] != relax.upper[j])
        {
            return;
        }
        let mut pairs = 0usize;
        for &(b, vb, bcomp) in z {
            if !relax.is_binary(b) {
                continue;
            }
            let implied = rhs - vb;
            for &(j, vj, jcomp) in z {
                if relax.is_binary(j) || relax.lower[j] == relax.upper[j] {
                    continue;
                }
                pairs += 1;
                if pairs > MAX_PAIRS_PER_ROW {
                    return;
                }
                // Range of z_j, capped by the row itself (z_b = 0 case).
                let colub = (relax.upper[j] - relax.lower[j]).min(rhs / vj);
                if !colub.is_finite() {
                    continue;
                }
                let mut impl_ub = implied / vj;
                if relax.integer[j] {
                    impl_ub = (impl_ub + feastol).floor();
                }
                if impl_ub >= colub - feastol {
                    continue;
                }
                // z_j ≤ coef·y + constant in terms of the original binary y.
                let (coef, constant) = if bcomp {
                    // z_b = 1 − y: y = 0 forces z_j ≤ impl_ub.
                    (colub - impl_ub, impl_ub)
                } else {
                    (impl_ub - colub, colub)
                };
                if jcomp {
                    // z_j = u − x ⇒ x ≥ −coef·y + (u − constant).
                    let bound = VarBound {
                        coef: -coef,
                        constant: relax.upper[j] - constant,
                    };
                    self.add_vlb(relax, j, b, bound, feastol, epsilon);
                } else {
                    let bound = VarBound {
                        coef,
                        constant: constant + relax.lower[j],
                    };
                    self.add_vub(relax, j, b, bound, feastol, epsilon);
                }
            }
        }
    }

    fn add_vub(
        &mut self,
        relax: &Relaxation,
        j: usize,
        y: usize,
        mut b: VarBound,
        feastol: f64,
        eps: f64,
    ) {
        let ub = relax.upper[j];
        if b.min_value() >= ub - feastol {
            return; // redundant
        }
        if b.max_value() > ub + eps {
            let newcoef = ub - b.min_value();
            if b.coef > 0.0 {
                b.coef = newcoef;
            } else {
                b.constant = ub;
                b.coef = -newcoef;
            }
        }
        let list = &mut self.vubs[j];
        if let Some(e) = list.iter_mut().find(|e| e.0 == y) {
            // Keep the pointwise tighter of the two (compare at y = 0 and 1).
            if b.at(0.0) <= e.1.at(0.0) && b.at(1.0) <= e.1.at(1.0) {
                e.1 = b;
            }
        } else if list.len() < MAX_PER_COL {
            list.push((y, b));
        }
    }

    fn add_vlb(
        &mut self,
        relax: &Relaxation,
        j: usize,
        y: usize,
        mut b: VarBound,
        feastol: f64,
        eps: f64,
    ) {
        let lb = relax.lower[j];
        if b.max_value() <= lb + feastol {
            return;
        }
        if b.min_value() < lb - eps {
            let newcoef = lb - b.max_value();
            if b.coef < 0.0 {
                b.coef = newcoef;
            } else {
                b.constant = lb;
                b.coef = -newcoef;
            }
        }
        let list = &mut self.vlbs[j];
        if let Some(e) = list.iter_mut().find(|e| e.0 == y) {
            if b.at(0.0) >= e.1.at(0.0) && b.at(1.0) >= e.1.at(1.0) {
                e.1 = b;
            }
        } else if list.len() < MAX_PER_COL {
            list.push((y, b));
        }
    }

    /// Implied bound cuts violated by more than `feastol` at `relax.x`, as
    /// `(inds, vals, rhs)` with `Σ vals·x ≤ rhs`.
    pub fn separate(&self, relax: &Relaxation, feastol: f64) -> Vec<(Vec<usize>, Vec<f64>, f64)> {
        let x = &relax.x;
        let mut out = Vec::new();
        for j in 0..relax.num_col {
            for &(y, b) in &self.vubs[j] {
                if relax.lower[y] == relax.upper[y] {
                    continue;
                }
                // x − c·y ≤ d
                if x[j] - b.at(x[y]) > feastol {
                    out.push((vec![j, y], vec![1.0, -b.coef], b.constant));
                }
            }
            for &(y, b) in &self.vlbs[j] {
                if relax.lower[y] == relax.upper[y] {
                    continue;
                }
                // −x + c·y ≤ −d
                if b.at(x[y]) - x[j] > feastol {
                    out.push((vec![j, y], vec![-1.0, b.coef], -b.constant));
                }
            }
        }
        out
    }
}

/// The finite model-row bound on one side (`upper = true` ⇒ row upper), or
/// ±∞ when that side was infinite in the model. Slack bounds replaced
/// infinite sides with activity bounds, which carry no implication.
fn relax_row_bound(relax: &Relaxation, i: usize, upper: bool) -> f64 {
    let r = &relax.rows[i];
    let (minact, maxact) =
        super::relax::activity_range(&r.inds, &r.vals, &relax.lower, &relax.upper);
    if upper {
        if r.upper < maxact {
            r.upper
        } else {
            f64::INFINITY
        }
    } else if r.lower > minact {
        r.lower
    } else {
        f64::NEG_INFINITY
    }
}
