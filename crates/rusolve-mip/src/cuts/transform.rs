//! Bound substitution for single-row cut generation: every column and slack
//! of a base row `Σ a_k x_k ≤ β` is shifted to its closest (simple or
//! variable) bound so that the transformed variables live in `[0, u_k]`,
//! and cuts derived in that space are mapped back (`untransform`).
//!
//! Algorithm, selection rules and tolerances follow HiGHS
//! `HighsTransformedLp` (reference/HiGHS/highs/mip/HighsTransformedLp.cpp:
//! constructor :15–136, transform :138–420, untransform :422–570; MIT, read
//! as a reference text and reimplemented). Bound substitution with variable
//! bounds is the c-MIR preprocessing of Marchand & Wolsey, "Aggregation and
//! mixed integer rounding to solve MIPs", Operations Research 49(3):363–371
//! (2001), §3.
//!
//! Simplification vs HiGHS: integer columns only ever use simple bounds
//! (HiGHS also substitutes variable bounds on general integers sitting at
//! an integral value); variable-bound ties are not broken by duals or node
//! counts (we have neither yet).

use super::implied::{VarBound, VarBounds};
use super::relax::{Relaxation, SparseSum};

/// Below this a bound range / coefficient is zero (HiGHS small_matrix_value).
pub(crate) const SMALL: f64 = 1e-9;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BoundType {
    SimpleLb,
    SimpleUb,
    VariableLb,
    VariableUb,
}

pub(crate) struct TransformedLp {
    bound_dist: Vec<f64>,
    simple_lb_dist: Vec<f64>,
    simple_ub_dist: Vec<f64>,
    lb_dist: Vec<f64>,
    ub_dist: Vec<f64>,
    best_vub: Vec<Option<(usize, VarBound)>>,
    best_vlb: Vec<Option<(usize, VarBound)>>,
    bound_type: Vec<BoundType>,
    sum: SparseSum,
    feastol: f64,
}

/// A base row in transformed space: `Σ vals_k·y_k ≤ rhs`, `y_k ∈ [0, upper_k]`,
/// LP point `solval`.
pub(crate) struct Transformed {
    pub inds: Vec<usize>,
    pub vals: Vec<f64>,
    pub rhs: f64,
    pub upper: Vec<f64>,
    pub solval: Vec<f64>,
    /// Every integer column ended up with a positive coefficient.
    pub integers_positive: bool,
}

fn clip(d: f64, feastol: f64) -> f64 {
    if d <= feastol { 0.0 } else { d }
}

impl TransformedLp {
    #[allow(clippy::needless_range_loop)]
    pub fn new(relax: &Relaxation, vb: &VarBounds, feastol: f64) -> Self {
        let dim = relax.num_ext();
        let mut t = Self {
            bound_dist: vec![0.0; dim],
            simple_lb_dist: vec![0.0; dim],
            simple_ub_dist: vec![0.0; dim],
            lb_dist: vec![0.0; dim],
            ub_dist: vec![0.0; dim],
            best_vub: vec![None; dim],
            best_vlb: vec![None; dim],
            bound_type: vec![BoundType::SimpleLb; dim],
            sum: SparseSum::new(dim),
            feastol,
        };
        let x = &relax.x;
        for j in 0..relax.num_col {
            let (lb, ub) = (relax.lower[j], relax.upper[j]);
            t.simple_ub_dist[j] = clip(ub - x[j], feastol);
            t.simple_lb_dist[j] = clip(x[j] - lb, feastol);
            if relax.integer[j] {
                t.lb_dist[j] = t.simple_lb_dist[j];
                t.ub_dist[j] = t.simple_ub_dist[j];
            } else {
                if lb == ub {
                    continue;
                }
                let (vub, best_ub) = best_vub(relax, j, &vb.vubs[j], feastol);
                let (vlb, best_lb) = best_vlb(relax, j, &vb.vlbs[j], feastol);
                t.best_vub[j] = vub;
                t.best_vlb[j] = vlb;
                t.lb_dist[j] = clip(x[j] - best_lb, feastol);
                t.ub_dist[j] = clip(best_ub - x[j], feastol);
            }
            t.bound_dist[j] = t.lb_dist[j].min(t.ub_dist[j]);
        }
        for i in 0..relax.rows.len() {
            let k = relax.num_col + i;
            let (lb, ub) = (relax.rows[i].lower, relax.rows[i].upper);
            if lb == ub {
                continue;
            }
            let a = relax.activity[i];
            t.lb_dist[k] = clip(a - lb, feastol);
            t.simple_lb_dist[k] = t.lb_dist[k];
            t.ub_dist[k] = clip(ub - a, feastol);
            t.simple_ub_dist[k] = t.ub_dist[k];
            t.bound_dist[k] = t.lb_dist[k].min(t.ub_dist[k]);
        }
        t
    }

    pub fn bound_distance(&self, k: usize) -> f64 {
        self.bound_dist[k]
    }

    /// Substitutes bounds in `Σ vals·x[inds] ≤ rhs` (inds over columns and
    /// slacks). `integers_positive` asks for positive integer coefficients
    /// where both bounds allow it. `prefer_vbds` substitutes an available
    /// variable bound even where the simple bound would relax the column
    /// away (HighsTransformedLp.cpp:262–293; the mod-k separator wants every
    /// continuous column at zero). `None` if a free variable remains or the
    /// row is empty and trivially satisfied.
    pub fn transform(
        &mut self,
        relax: &Relaxation,
        mut inds: Vec<usize>,
        mut vals: Vec<f64>,
        rhs: f64,
        mut integers_positive: bool,
        prefer_vbds: bool,
    ) -> Option<Transformed> {
        let feastol = self.feastol;
        let mut rhs = rhs;
        let mut used_vb = false;
        let mut i = 0;
        while i < inds.len() {
            let col = inds[i];
            let (lb, ub) = (relax.lb(col), relax.ub(col));
            if ub - lb < SMALL {
                rhs -= lb.min(ub) * vals[i];
                inds.swap_remove(i);
                vals.swap_remove(i);
                continue;
            }
            if lb == f64::NEG_INFINITY && ub == f64::INFINITY {
                self.sum.clear();
                return None;
            }
            let old = self.bound_type[col];
            if relax.is_integral(col) {
                // Decided after the loop (always a simple bound here).
                i += 1;
                continue;
            }
            let (lbd, ubd) = (self.lb_dist[col], self.ub_dist[col]);
            let (vlb, vub) = (self.best_vlb[col], self.best_vub[col]);
            let v = vals[i];
            let bt = if lbd < ubd - feastol {
                if vlb.is_none() {
                    BoundType::SimpleLb
                } else if prefer_vbds || v > 0.0 || self.simple_lb_dist[col] > lbd + feastol {
                    BoundType::VariableLb
                } else {
                    BoundType::SimpleLb
                }
            } else if ubd < lbd - feastol {
                if vub.is_none() {
                    BoundType::SimpleUb
                } else if prefer_vbds || v < 0.0 || self.simple_ub_dist[col] > ubd + feastol {
                    BoundType::VariableUb
                } else {
                    BoundType::SimpleUb
                }
            } else if v > 0.0 {
                if vlb.is_some() {
                    BoundType::VariableLb
                } else if prefer_vbds && vub.is_some() {
                    BoundType::VariableUb
                } else {
                    BoundType::SimpleLb
                }
            } else if vub.is_some() {
                BoundType::VariableUb
            } else if prefer_vbds && vlb.is_some() {
                BoundType::VariableLb
            } else {
                BoundType::SimpleUb
            };
            self.bound_type[col] = bt;
            match bt {
                BoundType::SimpleLb => {
                    if v > 0.0 {
                        // relax away with the lower bound
                        rhs -= lb * v;
                        self.bound_type[col] = old;
                        inds.swap_remove(i);
                        vals.swap_remove(i);
                        continue;
                    }
                }
                BoundType::SimpleUb => {
                    if v < 0.0 {
                        rhs -= ub * v;
                        self.bound_type[col] = old;
                        inds.swap_remove(i);
                        vals.swap_remove(i);
                        continue;
                    }
                }
                BoundType::VariableLb => {
                    let (y, b) = vlb.expect("chosen only when present");
                    rhs -= b.constant * v;
                    self.sum.add(y, v * b.coef);
                    used_vb = true;
                    if v > 0.0 {
                        self.bound_type[col] = old;
                        inds.swap_remove(i);
                        vals.swap_remove(i);
                        continue;
                    }
                }
                BoundType::VariableUb => {
                    let (y, b) = vub.expect("chosen only when present");
                    rhs -= b.constant * v;
                    self.sum.add(y, v * b.coef);
                    used_vb = true;
                    vals[i] = -v;
                    if vals[i] > 0.0 {
                        self.bound_type[col] = old;
                        inds.swap_remove(i);
                        vals.swap_remove(i);
                        continue;
                    }
                }
            }
            i += 1;
        }

        if used_vb {
            for (&k, &v) in inds.iter().zip(&vals) {
                if v != 0.0 {
                    self.sum.add(k, v);
                }
            }
            let (ni, nv) = self.sum.take(|_, v| v.abs() <= SMALL);
            inds = ni;
            vals = nv;
        } else {
            self.sum.clear();
        }

        for (k, &col) in inds.iter().enumerate() {
            if !relax.is_integral(col) {
                continue;
            }
            let (lb, ub) = (relax.lb(col), relax.ub(col));
            self.bound_type[col] = if integers_positive {
                if (lb != f64::NEG_INFINITY && vals[k] > 0.0) || ub == f64::INFINITY {
                    BoundType::SimpleLb
                } else {
                    BoundType::SimpleUb
                }
            } else if self.lb_dist[col] < self.ub_dist[col] {
                BoundType::SimpleLb
            } else {
                BoundType::SimpleUb
            };
        }

        let n = inds.len();
        let mut upper = vec![0.0; n];
        let mut solval = vec![0.0; n];
        for k in 0..n {
            let col = inds[k];
            let (lb, ub) = (relax.lb(col), relax.ub(col));
            upper[k] = ub - lb;
            match self.bound_type[col] {
                BoundType::SimpleLb => {
                    if lb == f64::NEG_INFINITY {
                        self.sum.clear();
                        return None;
                    }
                    rhs -= lb * vals[k];
                    solval[k] = self.lb_dist[col];
                }
                BoundType::SimpleUb => {
                    if ub == f64::INFINITY {
                        self.sum.clear();
                        return None;
                    }
                    rhs -= ub * vals[k];
                    vals[k] = -vals[k];
                    solval[k] = self.ub_dist[col];
                }
                BoundType::VariableLb => solval[k] = self.lb_dist[col],
                BoundType::VariableUb => solval[k] = self.ub_dist[col],
            }
            if relax.is_integral(col) {
                integers_positive = integers_positive && vals[k] > 0.0;
            }
        }
        if n == 0 && rhs >= -feastol {
            return None;
        }
        Some(Transformed {
            inds,
            vals,
            rhs,
            upper,
            solval,
            integers_positive,
        })
    }

    /// Maps a cut `Σ vals·y ≤ rhs` in transformed space back to the original
    /// columns (slacks resubstituted by their rows). `integral` rounds the
    /// result. `None` if a tiny leftover cannot be relaxed with a finite bound.
    pub fn untransform(
        &mut self,
        relax: &Relaxation,
        inds: &[usize],
        vals: &[f64],
        rhs: f64,
        integral: bool,
    ) -> Option<(Vec<usize>, Vec<f64>, f64)> {
        let n = relax.num_col;
        let mut rhs = rhs;
        for (&col, &v) in inds.iter().zip(vals) {
            if v == 0.0 {
                continue;
            }
            match self.bound_type[col] {
                BoundType::VariableLb => {
                    let (y, b) = self.best_vlb[col].expect("vlb type set");
                    rhs += b.constant * v;
                    self.sum.add(y, -v * b.coef);
                    self.sum.add(col, v);
                }
                BoundType::VariableUb => {
                    let (y, b) = self.best_vub[col].expect("vub type set");
                    rhs -= b.constant * v;
                    self.sum.add(y, v * b.coef);
                    self.sum.add(col, -v);
                }
                BoundType::SimpleLb => {
                    if col < n {
                        rhs += v * relax.lower[col];
                        self.sum.add(col, v);
                    } else {
                        let row = &relax.rows[col - n];
                        rhs += v * row.lower;
                        for (&j, &a) in row.inds.iter().zip(&row.vals) {
                            self.sum.add(j, v * a);
                        }
                    }
                }
                BoundType::SimpleUb => {
                    if col < n {
                        rhs -= v * relax.upper[col];
                        self.sum.add(col, -v);
                    } else {
                        let row = &relax.rows[col - n];
                        rhs -= v * row.upper;
                        for (&j, &a) in row.inds.iter().zip(&row.vals) {
                            self.sum.add(j, -v * a);
                        }
                    }
                }
            }
        }
        if integral {
            let (ri, rv) = self.sum.take(|_, v| v.abs() < 0.5);
            return Some((ri, rv.iter().map(|v| v.round()).collect(), rhs.round()));
        }
        let feastol = self.feastol;
        let mut abort = false;
        let (ri, rv) = self.sum.take(|col, v| {
            let a = v.abs();
            if a <= SMALL {
                return true;
            }
            if a <= feastol {
                if v > 0.0 {
                    if relax.lower[col] == f64::NEG_INFINITY {
                        abort = true;
                    } else {
                        rhs -= v * relax.lower[col];
                    }
                } else if relax.upper[col] == f64::INFINITY {
                    abort = true;
                } else {
                    rhs -= v * relax.upper[col];
                }
                return true;
            }
            false
        });
        if abort {
            return None;
        }
        Some((ri, rv, rhs))
    }
}

/// HighsImplications::getBestVub (HighsImplications.cpp:157–228) without the
/// node-count and dual tie-breakers. Returns the chosen bound and the upper
/// bound value at the LP point (the simple bound if none is chosen).
fn best_vub(
    relax: &Relaxation,
    col: usize,
    vubs: &[(usize, VarBound)],
    feastol: f64,
) -> (Option<(usize, VarBound)>, f64) {
    let x = &relax.x;
    let mut best_ub = relax.upper[col];
    let mut min_best_ub = best_ub;
    let mut best_dist = f64::INFINITY;
    let mut best = None;
    let range = relax.upper[col] - relax.lower[col];
    let scale = if range.is_finite() { 1.0 / range } else { 1.0 };
    for &(y, b) in vubs {
        if !b.coef.is_finite() || relax.lower[y] == relax.upper[y] {
            continue;
        }
        let val = b.at(x[y]);
        let mut dist = (val - x[col]).max(0.0);
        let ydist = feastol + if b.coef > 0.0 { 1.0 - x[y] } else { x[y] };
        if dist * dist > ydist * ydist * (1.0 + b.coef * b.coef) {
            continue;
        }
        dist *= scale;
        if dist <= best_dist + feastol {
            let minval = b.min_value();
            if dist < best_dist - feastol || minval < min_best_ub - feastol {
                best_ub = val;
                min_best_ub = minval;
                best = Some((y, b));
                best_dist = dist;
            }
        }
    }
    (best, best_ub)
}

/// HighsImplications::getBestVlb (HighsImplications.cpp:230–298), mirrored.
fn best_vlb(
    relax: &Relaxation,
    col: usize,
    vlbs: &[(usize, VarBound)],
    feastol: f64,
) -> (Option<(usize, VarBound)>, f64) {
    let x = &relax.x;
    let mut best_lb = relax.lower[col];
    let mut max_best_lb = best_lb;
    let mut best_dist = f64::INFINITY;
    let mut best = None;
    let range = relax.upper[col] - relax.lower[col];
    let scale = if range.is_finite() { 1.0 / range } else { 1.0 };
    for &(y, b) in vlbs {
        if !b.coef.is_finite() || relax.lower[y] == relax.upper[y] {
            continue;
        }
        let val = b.at(x[y]);
        let mut dist = (x[col] - val).max(0.0);
        let ydist = feastol + if b.coef > 0.0 { x[y] } else { 1.0 - x[y] };
        if dist * dist > ydist * ydist * (1.0 + b.coef * b.coef) {
            continue;
        }
        dist *= scale;
        if dist <= best_dist + feastol {
            let maxval = b.max_value();
            if dist < best_dist - feastol || maxval > max_best_lb + feastol {
                best_lb = val;
                max_best_lb = maxval;
                best = Some((y, b));
                best_dist = dist;
            }
        }
    }
    (best, best_lb)
}
