//! rusolve extension: LP-based outer approximation (OA) of a convex quadratic objective, for
//! convex MIQP. Design: `docs/design/miqp-oa.md`.
//!
//! This file is rusolve's own (not a HiGHS translation).
//!
//! Second-order cones (mixed-integer SOCP, `docs/design/misocp-oa.md`) share the machinery: an
//! [`OaCone`] is `||y(x)|| <= t(x)` with `y`, `t` affine. The cone function is convex, so the
//! gradient cut at a point is valid everywhere; see [`OaForm::cone_cut`].
//!
//! The objective `c'x + f(x)` with `f(x) = 1/2 x'Qx` (Q positive semidefinite for a
//! minimisation, negative semidefinite for a maximisation) is replaced by `c'x + sum_k eta_k`,
//! one epigraph column per block of Q, and the rows
//!
//! ```text
//!   eta_k >= f_k(xh) + grad f_k(xh)'(x - xh)        (minimisation; <= for a maximisation)
//! ```
//!
//! which are valid for every `xh` because `f_k` is convex (concave). Gradient cuts at the
//! solution of an LP relaxation give a polyhedral relaxation that is refined where it matters:
//! Duran and Grossmann, "An outer-approximation algorithm for a class of mixed-integer nonlinear
//! programs", Math. Programming 36 (1986) 307-339; the single-tree form that separates the cuts
//! inside one branch-and-bound search is Quesada and Grossmann, "An LP/NLP based branch and
//! bound algorithm for convex MINLP optimization problems", Computers & Chemical Engineering 16
//! (1992) 937-947.
//!
//! How the pieces fit here:
//! * **Root tangent.** The caller solves the continuous QP relaxation and passes its optimum to
//!   [`add_root_tangents`]; the tangent there makes the root LP value equal the QP value (the
//!   KKT conditions of the QP are those of the LP with that one cut).
//! * **In the tree.** [`OaSeparator`] is a [`MipHooks::separate`] hook: it cuts the node LP
//!   point `xh` (projected on the node bound box) for every block whose `eta_k` lies below
//!   `f_k`. Cuts pass the numerics filter of [`OaForm::tangent`] and then the cut-pool filters of
//!   the port. Hooks run in the presolved column space, so the tree separator needs
//!   `presolve=off`.
//! * **Exact incumbents.** A MIP solution whose `eta_k < f_k(x) - tol` is not a solution of the
//!   MIQP. [`OaOracle`] is a lazy-constraint oracle (`mip/lazy.rs`): it returns the tangent at
//!   that solution and the separator's earlier cuts the solution violates, and the loop solves
//!   again. The solution it accepts satisfies `eta_k >= f_k - tol` for every block, so its
//!   objective in the model equals the true one up to `tol` per block.
//! * **Disaggregation.** With one block per connected component of Q's off-diagonal graph (a
//!   diagonal Q gives one block per column) each `eta_k` gets its own sparse tangents: the
//!   extended-formulation effect measured by Vielma, Dunning, Huchette and Lubin, "Extended
//!   formulations in mixed integer conic quadratic programming", Math. Programming Computation 9
//!   (2017) 369-418.

use crate::lp_data::highs_run::Highs;
use crate::mip::hooks::{MipHooks, SeparationInfo, UserCut};
use crate::mip::lazy::{LazyConstraints, LazyReport};
use crate::prelude::*;
use std::sync::{Arc, Mutex, MutexGuard};

/// A coefficient below this fraction of the largest gradient entry is dropped from a cut when
/// the variable has a finite bound on the side that keeps the cut valid.
const TINY_REL: f64 = 1e-9;
/// A cut with a larger coefficient is not sent (numerics).
const MAX_COEF: f64 = 1e10;

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(|e| e.into_inner())
}

/// `Q` as symmetric triplets `(i, j, q_ij)`, both triangles present (the model's `hessian`).
#[derive(Clone, Debug, Default)]
pub struct QuadForm {
    pub n: usize,
    pub entries: Vec<(usize, usize, f64)>,
}

impl QuadForm {
    /// `1/2 x'Qx`.
    pub fn value(&self, x: &[f64]) -> f64 {
        0.5 * self
            .entries
            .iter()
            .map(|&(i, j, v)| v * x[i] * x[j])
            .sum::<f64>()
    }
}

/// One separable part of the quadratic: the columns it touches and `Q` restricted to them.
#[derive(Clone, Debug)]
pub struct OaBlock {
    pub cols: Vec<usize>,
    /// Row `k` of the block: `(local column, q)` pairs.
    rows: Vec<Vec<(usize, f64)>>,
    /// Column of the epigraph variable in the OA model.
    pub eta: usize,
}

/// An affine expression `c + sum_k val[k] * z[idx[k]]` over the OA model's columns.
#[derive(Clone, Debug, Default)]
pub struct Aff {
    pub idx: Vec<usize>,
    pub val: Vec<f64>,
    pub c: f64,
}

impl Aff {
    /// Builds the expression, merging repeated columns and dropping exact zeros.
    pub fn new(pairs: &[(usize, f64)], cst: f64) -> Aff {
        let mut p: Vec<(usize, f64)> = pairs.to_vec();
        p.sort_by_key(|&(j, _)| j);
        let mut merged: Vec<(usize, f64)> = Vec::new();
        for (j, v) in p {
            match merged.last_mut() {
                Some(last) if last.0 == j => last.1 += v,
                _ => merged.push((j, v)),
            }
        }
        merged.retain(|&(_, v)| v != 0.0);
        Aff {
            idx: merged.iter().map(|&(j, _)| j).collect(),
            val: merged.iter().map(|&(_, v)| v).collect(),
            c: cst,
        }
    }

    pub fn eval(&self, z: &[f64]) -> f64 {
        self.c
            + self
                .idx
                .iter()
                .zip(&self.val)
                .map(|(&j, &v)| v * z[j])
                .sum::<f64>()
    }
}

/// A second-order cone constraint `||y(x)||_2 <= t(x)`.
#[derive(Clone, Debug)]
pub struct OaCone {
    pub t: Aff,
    pub y: Vec<Aff>,
}

/// One second-order cone in the disaggregated (extended) formulation: the aggregate cone
/// `||y|| <= t` is replaced by `y_i^2 <= z_i t`, `sum z_i = t`, `z >= 0` (Vielma, Dunning,
/// Huchette and Lubin 2017). `zcols` are the `z_i` columns; each `y_i^2 <= z_i t` is an
/// [`OaCone`] of three columns in `OaForm::cones`.
#[derive(Clone, Debug)]
pub struct OaGroup {
    pub t: Aff,
    pub y: Vec<Aff>,
    pub zcols: Vec<usize>,
}

/// A column the OA model has besides the originals and the epigraph columns.
#[derive(Clone, Debug)]
pub struct ExtraCol {
    pub cost: f64,
    pub lo: f64,
    pub hi: f64,
    pub name: String,
}

/// The blocks of a quadratic, the cones, and the sign of the sense.
#[derive(Clone, Debug)]
pub struct OaForm {
    /// Original columns; the epigraph columns follow, then `extra`.
    pub n: usize,
    pub blocks: Vec<OaBlock>,
    /// `+1` minimisation (`f` convex, `eta >= tangent`), `-1` maximisation.
    pub sigma: f64,
    /// Cone constraints cut by [`OaForm::violated`].
    pub cones: Vec<OaCone>,
    /// Disaggregated cones (their `z` columns and linking rows).
    pub groups: Vec<OaGroup>,
    pub extra: Vec<ExtraCol>,
}

/// Initial cuts per cone: tangents at `+-e_i` for at most this many components `i`.
const INIT_CUT_DIRS: usize = 32;
/// Most cone cuts one separation call returns (the most violated).
const MAX_CONE_CUTS: usize = 300;

fn find(parent: &mut [usize], mut a: usize) -> usize {
    while parent[a] != a {
        parent[a] = parent[parent[a]];
        a = parent[a];
    }
    a
}

impl OaForm {
    /// Blocks of `q`. `disaggregate`: one block per connected component of the off-diagonal
    /// graph; otherwise a single block holds every column with a nonzero entry.
    pub fn new(q: &QuadForm, maximize: bool, disaggregate: bool) -> OaForm {
        let n = q.n;
        let mut parent: Vec<usize> = (0..n).collect();
        let mut touched = vec![false; n];
        for &(i, j, v) in &q.entries {
            if v == 0.0 {
                continue;
            }
            touched[i] = true;
            touched[j] = true;
            if !disaggregate || i != j {
                let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                if a != b {
                    parent[a] = b;
                }
            }
        }
        if !disaggregate {
            // one block: join every touched column
            let first = touched.iter().position(|&t| t);
            if let Some(f) = first {
                for col in 0..n {
                    if touched[col] {
                        let (a, b) = (find(&mut parent, col), find(&mut parent, f));
                        if a != b {
                            parent[a] = b;
                        }
                    }
                }
            }
        }
        let mut block_of = vec![usize::MAX; n];
        let mut blocks: Vec<OaBlock> = Vec::new();
        for col in 0..n {
            if !touched[col] {
                continue;
            }
            let r = find(&mut parent, col);
            if block_of[r] == usize::MAX {
                block_of[r] = blocks.len();
                blocks.push(OaBlock {
                    cols: Vec::new(),
                    rows: Vec::new(),
                    eta: n + blocks.len(),
                });
            }
            blocks[block_of[r]].cols.push(col);
        }
        let mut local = vec![0usize; n];
        for b in &mut blocks {
            for (l, &col) in b.cols.iter().enumerate() {
                local[col] = l;
            }
            b.rows = vec![Vec::new(); b.cols.len()];
        }
        for &(i, j, v) in &q.entries {
            if v == 0.0 {
                continue;
            }
            let k = block_of[find(&mut parent, i)];
            blocks[k].rows[local[i]].push((local[j], v));
        }
        OaForm {
            n,
            blocks,
            sigma: if maximize { -1.0 } else { 1.0 },
            cones: Vec::new(),
            groups: Vec::new(),
            extra: Vec::new(),
        }
    }

    /// Columns of the OA model: the originals, one epigraph column per block, the extra columns.
    pub fn num_cols(&self) -> usize {
        self.n + self.blocks.len() + self.extra.len()
    }

    /// Appends an extra column and returns its index.
    pub fn add_extra_col(&mut self, cost: f64, lo: f64, hi: f64, name: String) -> usize {
        let j = self.num_cols();
        self.extra.push(ExtraCol { cost, lo, hi, name });
        j
    }

    /// Adds the cone `||y|| <= t` as it is.
    pub fn add_cone(&mut self, t: Aff, y: Vec<Aff>) {
        self.cones.push(OaCone { t, y });
    }

    /// Adds `||y|| <= t` in the disaggregated formulation: one `z_i` column per component and the
    /// cones `y_i^2 <= z_i t`, written as `||(y_i, (z_i - t)/2)|| <= (z_i + t)/2`.
    pub fn add_cone_disaggregated(&mut self, name: &str, t: Aff, y: Vec<Aff>) {
        let mut zcols = Vec::with_capacity(y.len());
        for (i, yi) in y.iter().enumerate() {
            let zc = self.add_extra_col(0.0, 0.0, kHighsInf, format!("{name}_z{i}"));
            zcols.push(zc);
            let half = |ts: f64| {
                let mut pairs = vec![(zc, 0.5)];
                pairs.extend(t.idx.iter().zip(&t.val).map(|(&j, &v)| (j, ts * 0.5 * v)));
                Aff::new(&pairs, ts * 0.5 * t.c)
            };
            self.cones.push(OaCone {
                t: half(1.0),
                y: vec![yi.clone(), half(-1.0)],
            });
        }
        self.groups.push(OaGroup { t, y, zcols });
    }

    /// The point of the OA model for the original columns `x`: the `z` columns of the
    /// disaggregated cones are set to `t y_i^2 / ||y||^2` (their sum is `t`); other extra and the
    /// epigraph columns are 0 and the caller sets them.
    pub fn lift(&self, x: &[f64]) -> Vec<f64> {
        let mut z = vec![0.0; self.num_cols()];
        let m = x.len().min(self.n);
        z[..m].copy_from_slice(&x[..m]);
        for g in &self.groups {
            let t = g.t.eval(&z).max(0.0);
            let yv: Vec<f64> = g.y.iter().map(|a| a.eval(&z)).collect();
            let s2: f64 = yv.iter().map(|v| v * v).sum();
            let k = g.zcols.len().max(1) as f64;
            for (i, &zc) in g.zcols.iter().enumerate() {
                z[zc] = if s2 > 1e-300 {
                    t * yv[i] * yv[i] / s2
                } else {
                    t / k
                };
            }
        }
        z
    }

    /// Equality rows `sum z_i - t(x) = 0` of the disaggregated cones, as `(rhs, indices, values)`.
    pub fn link_rows(&self) -> Vec<(f64, Vec<HighsInt>, Vec<f64>)> {
        self.groups
            .iter()
            .map(|g| {
                let mut pairs: Vec<(usize, f64)> = g.zcols.iter().map(|&zc| (zc, 1.0)).collect();
                pairs.extend(g.t.idx.iter().zip(&g.t.val).map(|(&j, &v)| (j, -v)));
                let a = Aff::new(&pairs, 0.0);
                (g.t.c, a.idx.iter().map(|&j| j as HighsInt).collect(), a.val)
            })
            .collect()
    }

    /// `||y(z)|| - t(z)` of cone `k`: positive when the point is outside the cone.
    pub fn cone_violation(&self, k: usize, z: &[f64]) -> f64 {
        let cn = &self.cones[k];
        let n2: f64 = cn.y.iter().map(|a| a.eval(z).powi(2)).sum();
        n2.sqrt() - cn.t.eval(z)
    }

    /// The cut `u'y(x) <= t(x)` for a vector `u` with `||u|| = 1` as a row `a z <= rhs`.
    /// Valid for every `x` of the cone because `u'y <= ||y|| <= t` (Cauchy-Schwarz). `None` when
    /// the row is empty or numerically unsafe. `glo` / `ghi` are the global column bounds, used
    /// to drop negligible coefficients without losing validity.
    pub fn cone_cut(&self, k: usize, u: &[f64], glo: &[f64], ghi: &[f64]) -> Option<UserCut> {
        let cn = &self.cones[k];
        let mut acc: Vec<(usize, f64)> = Vec::new();
        let mut rhs = cn.t.c;
        for (ui, yi) in u.iter().zip(&cn.y) {
            if *ui == 0.0 {
                continue;
            }
            rhs -= ui * yi.c;
            acc.extend(yi.idx.iter().zip(&yi.val).map(|(&j, &v)| (j, ui * v)));
        }
        acc.extend(cn.t.idx.iter().zip(&cn.t.val).map(|(&j, &v)| (j, -v)));
        let a = Aff::new(&acc, 0.0);
        let gmax = a.val.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if a.idx.is_empty() || !gmax.is_finite() {
            return None;
        }
        let mut indices: Vec<HighsInt> = Vec::new();
        let mut values: Vec<f64> = Vec::new();
        for (&j, &v) in a.idx.iter().zip(&a.val) {
            if v.abs() <= TINY_REL * gmax && j < glo.len() {
                let m = if v > 0.0 { v * glo[j] } else { v * ghi[j] };
                if m.is_finite() {
                    rhs -= m;
                    continue;
                }
            }
            indices.push(j as HighsInt);
            values.push(v);
        }
        if indices.is_empty() {
            return None;
        }
        if !rhs.is_finite() || values.iter().any(|v| !v.is_finite() || v.abs() > MAX_COEF) {
            return None;
        }
        Some(UserCut {
            indices,
            values,
            rhs,
        })
    }

    /// The gradient cut of cone `k` at the point `z`, and the cone's violation there. At the
    /// apex (`y = 0`, where the norm has no gradient) the subgradient `u = e_0` is used: the cut
    /// `y_0 <= t` is valid and is violated exactly when `t < 0`.
    pub fn cone_tangent(
        &self,
        k: usize,
        z: &[f64],
        glo: &[f64],
        ghi: &[f64],
    ) -> Option<(f64, UserCut)> {
        let cn = &self.cones[k];
        let yv: Vec<f64> = cn.y.iter().map(|a| a.eval(z)).collect();
        let tv = cn.t.eval(z);
        let nrm = yv.iter().map(|v| v * v).sum::<f64>().sqrt();
        if !nrm.is_finite() || !tv.is_finite() {
            return None;
        }
        let u: Vec<f64> = if nrm > 1e-14 {
            yv.iter().map(|v| v / nrm).collect()
        } else {
            let mut e = vec![0.0; yv.len()];
            if let Some(f) = e.first_mut() {
                *f = 1.0;
            }
            e
        };
        self.cone_cut(k, &u, glo, ghi).map(|cut| (nrm - tv, cut))
    }

    /// Cuts every cone has from the start: tangents at `+-e_i` (the cone contains `|y_i| <= t`),
    /// which bound `t` below and keep the first LP from being unbounded.
    pub fn initial_cone_cuts(&self, glo: &[f64], ghi: &[f64]) -> Vec<UserCut> {
        let mut out = Vec::new();
        for (k, cn) in self.cones.iter().enumerate() {
            let m = cn.y.len();
            if m == 0 {
                out.extend(self.cone_cut(k, &[], glo, ghi));
                continue;
            }
            for i in 0..m.min(INIT_CUT_DIRS) {
                for s in [1.0, -1.0] {
                    let mut u = vec![0.0; m];
                    u[i] = s;
                    out.extend(self.cone_cut(k, &u, glo, ghi));
                }
            }
        }
        out
    }

    /// `(f_k(x), gradient over the block's columns)`.
    fn eval(&self, k: usize, x: &[f64]) -> (f64, Vec<f64>) {
        let b = &self.blocks[k];
        let xl: Vec<f64> = b.cols.iter().map(|&col| x[col]).collect();
        let g: Vec<f64> = b
            .rows
            .iter()
            .map(|r| r.iter().map(|&(l, v)| v * xl[l]).sum())
            .collect();
        let f = 0.5 * g.iter().zip(&xl).map(|(a, b)| a * b).sum::<f64>();
        (f, g)
    }

    /// The tangent of block `k` at `x` projected on `[lo, hi]` as a row `a z <= rhs`, or `None`
    /// when it is trivial or numerically unsafe. `glo` / `ghi` are the global bounds, used to
    /// drop negligible coefficients without losing validity.
    ///
    /// With `xp` the projected point and `g = Q xp` on the block, `g'xp = 2 f(xp)`, so the cut
    /// is `sigma (g'x - eta) <= sigma f(xp)`.
    pub fn tangent(
        &self,
        k: usize,
        x: &[f64],
        lo: &[f64],
        hi: &[f64],
        glo: &[f64],
        ghi: &[f64],
    ) -> Option<UserCut> {
        let b = &self.blocks[k];
        let xp: Vec<f64> = (0..self.n)
            .map(|col| {
                let mut v = x[col];
                if v < lo[col] {
                    v = lo[col];
                }
                if v > hi[col] {
                    v = hi[col];
                }
                v
            })
            .collect();
        let (f, g) = self.eval(k, &xp);
        let gmax = g.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if !f.is_finite() || !gmax.is_finite() || gmax == 0.0 {
            return None;
        }
        let s = self.sigma;
        let mut rhs = s * f;
        let mut indices: Vec<HighsInt> = Vec::new();
        let mut values: Vec<f64> = Vec::new();
        for (l, &col) in b.cols.iter().enumerate() {
            let a = s * g[l];
            if a.abs() <= TINY_REL * gmax {
                // a z_c >= m on the global box: the cut without the term, with rhs - m, is valid
                let m = if a > 0.0 { a * glo[col] } else { a * ghi[col] };
                if m.is_finite() {
                    rhs -= m;
                    continue;
                }
            }
            indices.push(col as HighsInt);
            values.push(a);
        }
        if indices.is_empty() {
            return None;
        }
        indices.push(b.eta as HighsInt);
        values.push(-s);
        if !rhs.is_finite() || values.iter().any(|v| !v.is_finite() || v.abs() > MAX_COEF) {
            return None;
        }
        Some(UserCut {
            indices,
            values,
            rhs,
        })
    }

    /// Violation of `cut` at the OA point `z` (original columns then epigraph columns).
    pub fn violation(cut: &UserCut, z: &[f64]) -> f64 {
        cut.indices
            .iter()
            .zip(&cut.values)
            .map(|(&col, &v)| v * z[col as usize])
            .sum::<f64>()
            - cut.rhs
    }

    /// Tangents of every block that the point `z` violates by more than `tol`, most violated
    /// first.
    pub fn violated(
        &self,
        z: &[f64],
        lo: &[f64],
        hi: &[f64],
        glo: &[f64],
        ghi: &[f64],
        tol: f64,
    ) -> Vec<(f64, UserCut)> {
        let mut out = Vec::new();
        if z.len() < self.num_cols() {
            return out;
        }
        for k in 0..self.blocks.len() {
            if let Some(cut) = self.tangent(k, &z[..self.n], lo, hi, glo, ghi) {
                let v = Self::violation(&cut, z);
                if v > tol {
                    out.push((v, cut));
                }
            }
        }
        let mut cone_cuts: Vec<(f64, UserCut)> = Vec::new();
        for k in 0..self.cones.len() {
            if self.cone_violation(k, z) > tol
                && let Some((v, cut)) = self.cone_tangent(k, z, glo, ghi)
                && v > tol
            {
                cone_cuts.push((v, cut));
            }
        }
        cone_cuts.sort_by(|a, b| b.0.total_cmp(&a.0));
        cone_cuts.truncate(MAX_CONE_CUTS);
        out.extend(cone_cuts);
        out.sort_by(|a, b| b.0.total_cmp(&a.0));
        out
    }

    /// The OA model: `lp` with one epigraph column per block (cost 1, `eta >= 0` for a
    /// minimisation and `eta <= 0` for a maximisation: `f` has the sign of the sense).
    pub fn extend_lp(&self, lp: &mut HighsLp) {
        let n = lp.num_col_ as usize;
        debug_assert_eq!(n, self.n);
        for (k, _) in self.blocks.iter().enumerate() {
            lp.col_cost_.push(1.0);
            if self.sigma > 0.0 {
                lp.col_lower_.push(0.0);
                lp.col_upper_.push(kHighsInf);
            } else {
                lp.col_lower_.push(-kHighsInf);
                lp.col_upper_.push(0.0);
            }
            if !lp.integrality_.is_empty() {
                lp.integrality_.push(HighsVarType::kContinuous);
            }
            if lp.col_names_.len() == n + k {
                lp.col_names_.push(format!("oa_eta_{k}"));
            }
            let last = *lp.a_matrix_.start_.last().unwrap_or(&0);
            lp.a_matrix_.start_.push(last);
        }
        for (j, e) in self.extra.iter().enumerate() {
            lp.col_cost_.push(e.cost);
            lp.col_lower_.push(e.lo);
            lp.col_upper_.push(e.hi);
            if !lp.integrality_.is_empty() {
                lp.integrality_.push(HighsVarType::kContinuous);
            }
            if lp.col_names_.len() == n + self.blocks.len() + j {
                lp.col_names_.push(e.name.clone());
            }
            let last = *lp.a_matrix_.start_.last().unwrap_or(&0);
            lp.a_matrix_.start_.push(last);
        }
        let n2 = self.num_cols();
        lp.a_matrix_.num_col_ = n2 as HighsInt;
        lp.num_col_ = n2 as HighsInt;
    }
}

/// Counters of what the engine did.
#[derive(Clone, Debug, Default)]
pub struct OaStats {
    pub root_tangents: u64,
    pub tree_calls: u64,
    pub tree_cuts: u64,
    pub lazy_cuts_new: u64,
    pub lazy_cuts_from_pool: u64,
    pub cuts_unsafe: u64,
}

pub type OaStatsHandle = Arc<Mutex<OaStats>>;

/// Adds the tangent of every block at `x` (the QP relaxation's optimum) as model rows.
pub fn add_root_tangents(
    highs: &mut Highs,
    form: &OaForm,
    x: &[f64],
    glo: &[f64],
    ghi: &[f64],
    stats: &OaStatsHandle,
) -> usize {
    let mut added = 0;
    for k in 0..form.blocks.len() {
        match form.tangent(k, x, glo, ghi, glo, ghi) {
            Some(cut) => {
                highs.add_row(
                    -kHighsInf,
                    cut.rhs,
                    cut.indices.len() as HighsInt,
                    &cut.indices,
                    &cut.values,
                );
                added += 1;
            }
            None => lock(stats).cuts_unsafe += 1,
        }
    }
    lock(stats).root_tangents += added as u64;
    added
}

/// Adds the cone cuts every cone has from the start (see [`OaForm::initial_cone_cuts`]), the
/// tangent of every cone at `z` (the relaxation's optimum, in OA columns; `None` skips it) and the
/// equality rows that tie the disaggregated `z` columns to `t`. Returns the number of rows added.
pub fn add_cone_rows(
    highs: &mut Highs,
    form: &OaForm,
    z: Option<&[f64]>,
    glo: &[f64],
    ghi: &[f64],
    stats: &OaStatsHandle,
) -> usize {
    let mut rows = 0;
    for (rhs, idx, val) in form.link_rows() {
        highs.add_row(rhs, rhs, idx.len() as HighsInt, &idx, &val);
        rows += 1;
    }
    let mut cuts = form.initial_cone_cuts(glo, ghi);
    if let Some(z) = z {
        for k in 0..form.cones.len() {
            match form.cone_tangent(k, z, glo, ghi) {
                Some((_, cut)) => {
                    cuts.push(cut);
                    lock(stats).root_tangents += 1;
                }
                None => lock(stats).cuts_unsafe += 1,
            }
        }
    }
    for cut in &cuts {
        highs.add_row(
            -kHighsInf,
            cut.rhs,
            cut.indices.len() as HighsInt,
            &cut.indices,
            &cut.values,
        );
        rows += 1;
    }
    rows
}

/// Cuts the separator found, kept for the oracle (the tree's pool is lost at each re-solve).
pub type CutMemory = Arc<Mutex<Vec<UserCut>>>;
const MEMORY_CAP: usize = 100_000;

/// The in-tree separator (`MipHooks::separate`).
pub struct OaSeparator {
    pub form: Arc<OaForm>,
    pub glo: Arc<Vec<f64>>,
    pub ghi: Arc<Vec<f64>>,
    pub memory: CutMemory,
    pub stats: OaStatsHandle,
}

impl MipHooks for OaSeparator {
    fn name(&self) -> &str {
        "miqp_oa"
    }

    fn separate(&mut self, info: &SeparationInfo<'_>) -> Vec<UserCut> {
        let nc = self.form.num_cols();
        if info.solution.len() != nc || info.col_lower.len() != nc {
            return Vec::new();
        }
        let tol = 1.5 * info.feastol.max(1e-9);
        let v = self.form.violated(
            info.solution,
            &info.col_lower[..self.form.n],
            &info.col_upper[..self.form.n],
            &self.glo,
            &self.ghi,
            tol,
        );
        let mut st = lock(&self.stats);
        st.tree_calls += 1;
        st.tree_cuts += v.len() as u64;
        drop(st);
        let cuts: Vec<UserCut> = v.into_iter().map(|(_, cut)| cut).collect();
        let mut mem = lock(&self.memory);
        if mem.len() + cuts.len() <= MEMORY_CAP {
            mem.extend(cuts.iter().cloned());
        }
        cuts
    }
}

/// The lazy-constraint oracle that makes incumbents exact.
pub struct OaOracle {
    pub form: Arc<OaForm>,
    pub glo: Arc<Vec<f64>>,
    pub ghi: Arc<Vec<f64>>,
    pub memory: CutMemory,
    pub stats: OaStatsHandle,
    /// Above the loop's own feasibility tolerance, so a returned row is always added.
    pub tol: f64,
    /// Rows from the separator's memory added per round, most violated first.
    pub memory_rows: usize,
    /// Primal polish: maps a candidate (OA columns) to a feasible point of the model (OA
    /// columns) or `None`. The point is handed to the next solve as a start.
    pub polish: Option<Polisher>,
    pending_start: Option<Vec<f64>>,
}

/// A primal heuristic for the oracle: candidate in, feasible point (OA columns) out.
pub type Polisher = Box<dyn FnMut(&[f64]) -> Option<Vec<f64>> + Send>;

impl LazyConstraints for OaOracle {
    fn before_resolve(&mut self, highs: &mut Highs) {
        if let Some(col_value) = self.pending_start.take() {
            let sol = HighsSolution {
                col_value,
                ..Default::default()
            };
            let _ = highs.set_solution(&sol);
        }
    }

    fn lazy_constraints(&mut self, solution: &[f64]) -> Vec<UserCut> {
        if let Some(p) = self.polish.as_mut()
            && let Some(point) = p(solution)
        {
            self.pending_start = Some(point);
        }
        let mut out: Vec<UserCut> = self
            .form
            .violated(
                solution, &self.glo, &self.ghi, &self.glo, &self.ghi, self.tol,
            )
            .into_iter()
            .map(|(_, cut)| cut)
            .collect();
        let new_rows = out.len();
        if !out.is_empty() && self.memory_rows > 0 {
            let mem = lock(&self.memory);
            let mut hits: Vec<(f64, usize)> = mem
                .iter()
                .enumerate()
                .map(|(i, cut)| (OaForm::violation(cut, solution), i))
                .filter(|(v, _)| *v > self.tol)
                .collect();
            hits.sort_by(|a, b| b.0.total_cmp(&a.0));
            hits.truncate(self.memory_rows);
            let n_hits = hits.len();
            out.extend(hits.into_iter().map(|(_, i)| mem[i].clone()));
            drop(mem);
            lock(&self.stats).lazy_cuts_from_pool += n_hits as u64;
        }
        lock(&self.stats).lazy_cuts_new += new_rows as u64;
        out
    }
}

/// Engine settings.
#[derive(Clone, Copy, Debug)]
pub struct OaConfig {
    /// Separate in the tree (needs `presolve=off`).
    pub tree: bool,
    pub max_rounds: usize,
    pub memory_rows: usize,
}

impl Default for OaConfig {
    fn default() -> Self {
        OaConfig {
            tree: true,
            max_rounds: 500,
            memory_rows: 50,
        }
    }
}

/// Result of [`run_oa`].
pub struct OaOutcome {
    pub status: HighsStatus,
    pub lazy: LazyReport,
    pub stats: OaStats,
}

/// Solve the loaded OA model: install the tree separator (when `cfg.tree`), then run the
/// lazy-constraint loop with the exactness oracle within `deadline`.
pub fn run_oa(
    highs: &mut Highs,
    form: Arc<OaForm>,
    glo: Arc<Vec<f64>>,
    ghi: Arc<Vec<f64>>,
    stats: OaStatsHandle,
    cfg: OaConfig,
    deadline: Option<std::time::Instant>,
) -> OaOutcome {
    run_oa_with(highs, form, glo, ghi, stats, cfg, deadline, None)
}

/// [`run_oa`] with a primal polish applied to every candidate of the lazy loop.
#[allow(clippy::too_many_arguments)]
pub fn run_oa_with(
    highs: &mut Highs,
    form: Arc<OaForm>,
    glo: Arc<Vec<f64>>,
    ghi: Arc<Vec<f64>>,
    stats: OaStatsHandle,
    cfg: OaConfig,
    deadline: Option<std::time::Instant>,
    polish: Option<Polisher>,
) -> OaOutcome {
    let memory: CutMemory = Arc::new(Mutex::new(Vec::new()));
    if cfg.tree {
        highs.set_mip_hooks(Box::new(OaSeparator {
            form: Arc::clone(&form),
            glo: Arc::clone(&glo),
            ghi: Arc::clone(&ghi),
            memory: Arc::clone(&memory),
            stats: Arc::clone(&stats),
        }));
    }
    let feastol = highs.options_.mip_feasibility_tolerance;
    let mut oracle = OaOracle {
        form,
        glo,
        ghi,
        memory,
        stats: Arc::clone(&stats),
        tol: 2.0 * feastol.max(1e-9),
        memory_rows: cfg.memory_rows,
        polish,
        pending_start: None,
    };
    let (status, lazy) =
        highs.run_with_lazy_constraints_until(&mut oracle, cfg.max_rounds, deadline);
    if cfg.tree {
        highs.clear_mip_hooks();
    }
    let stats = lock(&stats).clone();
    OaOutcome {
        status,
        lazy,
        stats,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q2() -> QuadForm {
        // f = 1/2 (2 x0^2 + 2 x1^2 + 2 x0 x1) -> Q = [[2,1],[1,2]] with a separate x2: 4 x2^2/2
        QuadForm {
            n: 3,
            entries: vec![
                (0, 0, 2.0),
                (0, 1, 1.0),
                (1, 0, 1.0),
                (1, 1, 2.0),
                (2, 2, 4.0),
            ],
        }
    }

    #[test]
    fn blocks_follow_components() {
        let f = OaForm::new(&q2(), false, true);
        assert_eq!(f.blocks.len(), 2);
        assert_eq!(f.blocks[0].cols, vec![0, 1]);
        assert_eq!(f.blocks[1].cols, vec![2]);
        let g = OaForm::new(&q2(), false, false);
        assert_eq!(g.blocks.len(), 1);
        assert_eq!(g.blocks[0].cols, vec![0, 1, 2]);
    }

    #[test]
    fn tangent_is_valid_and_tight() {
        let q = q2();
        let f = OaForm::new(&q, false, true);
        let inf = f64::INFINITY;
        let (lo, hi) = (vec![-inf; 3], vec![inf; 3]);
        let x = [1.0, 2.0, 3.0];
        let cut = f.tangent(0, &x, &lo, &hi, &lo, &hi).unwrap();
        // the cut is tight at x with eta = f_0(x) = 1/2 (2 + 8 + 4) = 7
        let z = [1.0, 2.0, 3.0, 7.0, 18.0];
        assert!(OaForm::violation(&cut, &z).abs() < 1e-12);
        // valid at other points: eta = f_0(y) satisfies the row
        for y in [[0.0, 0.0], [-3.0, 1.0], [5.0, -2.0]] {
            let fy = 0.5 * (2.0 * y[0] * y[0] + 2.0 * y[1] * y[1] + 2.0 * y[0] * y[1]);
            let z = [y[0], y[1], 0.0, fy, 0.0];
            assert!(OaForm::violation(&cut, &z) <= 1e-12);
        }
        // below f -> violated
        let z = [1.0, 2.0, 3.0, 6.0, 18.0];
        assert!(OaForm::violation(&cut, &z) > 0.9);
    }

    #[test]
    fn maximisation_flips_the_row() {
        // Q = -2 on x0: f = -x0^2, concave; eta <= tangent
        let q = QuadForm {
            n: 1,
            entries: vec![(0, 0, -2.0)],
        };
        let f = OaForm::new(&q, true, true);
        let inf = f64::INFINITY;
        let cut = f
            .tangent(0, &[2.0], &[-inf], &[inf], &[-inf], &[inf])
            .unwrap();
        // eta = -4 at x = 2 is on the tangent; eta = -3 is above f(2) -> violated
        assert!(OaForm::violation(&cut, &[2.0, -4.0]).abs() < 1e-12);
        assert!(OaForm::violation(&cut, &[2.0, -3.0]) > 0.9);
        assert!(OaForm::violation(&cut, &[0.0, 0.0]) <= 1e-12);
    }

    fn unit_cone_form(disaggregate: bool) -> OaForm {
        // ||(x1, x2)|| <= x0 over columns 0..3
        let mut f = OaForm::new(
            &QuadForm {
                n: 3,
                entries: vec![],
            },
            false,
            false,
        );
        let t = Aff::new(&[(0, 1.0)], 0.0);
        let y = vec![Aff::new(&[(1, 1.0)], 0.0), Aff::new(&[(2, 1.0)], 0.0)];
        if disaggregate {
            f.add_cone_disaggregated("c", t, y);
        } else {
            f.add_cone(t, y);
        }
        f
    }

    #[test]
    fn cone_cut_is_valid_inside_and_tight_at_the_point() {
        let f = unit_cone_form(false);
        let inf = f64::INFINITY;
        let (lo, hi) = (vec![-inf; 3], vec![inf; 3]);
        let z = [1.0, 2.0, 3.0]; // outside: ||(2, 3)|| = 3.61 > 1
        let (viol, cut) = f.cone_tangent(0, &z, &lo, &hi).unwrap();
        assert!((viol - (13.0f64.sqrt() - 1.0)).abs() < 1e-12);
        assert!(OaForm::violation(&cut, &z) > 1.0); // cuts the point off
        // every point of the cone satisfies the cut
        for k in 0..200 {
            let a = k as f64 * 0.37;
            let r = 0.1 + (k % 17) as f64 * 0.3;
            let (y1, y2) = (r * a.cos(), r * a.sin());
            for slack in [0.0, 0.5] {
                let p = [r + slack, y1, y2];
                assert!(OaForm::violation(&cut, &p) <= 1e-9, "k {k}");
            }
        }
    }

    #[test]
    fn cone_apex_uses_a_subgradient() {
        let f = unit_cone_form(false);
        let inf = f64::INFINITY;
        let (lo, hi) = (vec![-inf; 3], vec![inf; 3]);
        // at the apex y = 0, t = -1 (outside): the cut y_0 <= t is valid and cuts the point off
        let z = [-1.0, 0.0, 0.0];
        let (viol, cut) = f.cone_tangent(0, &z, &lo, &hi).unwrap();
        assert!((viol - 1.0).abs() < 1e-12);
        assert!(OaForm::violation(&cut, &z) > 0.5);
        for p in [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, -1.0, 1.0]] {
            assert!(OaForm::violation(&cut, &p) <= 1e-9);
        }
    }

    #[test]
    fn disaggregated_lift_satisfies_the_link_row_and_the_sub_cones() {
        let f = unit_cone_form(true);
        assert_eq!(f.cones.len(), 2);
        assert_eq!(f.num_cols(), 5);
        let z = f.lift(&[2.0, 1.0, 1.0]); // inside: ||(1, 1)|| = 1.41 <= 2
        for (rhs, idx, val) in f.link_rows() {
            let act: f64 = idx.iter().zip(&val).map(|(&j, &v)| v * z[j as usize]).sum();
            assert!((act - rhs).abs() < 1e-12);
        }
        for k in 0..f.cones.len() {
            assert!(f.cone_violation(k, &z) <= 1e-12, "sub-cone {k}");
        }
        // a point outside the cone cannot satisfy both sub-cones and the link row
        let out = f.lift(&[1.0, 2.0, 2.0]);
        let worst = (0..f.cones.len())
            .map(|k| f.cone_violation(k, &out))
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(worst > 1e-3);
    }
}
