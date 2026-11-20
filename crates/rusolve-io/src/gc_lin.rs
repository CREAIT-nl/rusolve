//! Linearization of the [`GcModel`] objects.
//!
//! `linearize` derives rows and columns for a solve from the user's objects;
//! the objects themselves stay in the model, and the checker evaluates them,
//! never the derived rows. New columns are appended after the model's `n`
//! columns and new rows after its rows, so cutting a solution back to the
//! user's `n` columns is exact.
//!
//! Techniques and where they come from:
//! - Indicator `z = v -> a'x <= b` by big-M with `M` the largest value of
//!   `a'x - b` over the (propagated) column box: Bonami, Lodi, Tramontani,
//!   Wiese, "On mathematical programming with indicator constraints",
//!   Math. Prog. 151 (2015). `M` above `gc_bigm_max` is a typed error, never
//!   a silent huge coefficient.
//! - SOS1 / SOS2 with finite member bounds: the standard binary formulation
//!   (members' bounds linked to selector binaries; for SOS2 one binary per
//!   adjacent pair), Beale and Tomlin, "Special facilities in a general
//!   mathematical programming system for non-convex problems using ordered
//!   sets of variables" (1970); formulation as in Williams, "Model Building in
//!   Mathematical Programming", 5th ed., ch. 9.
//! - Piecewise-linear `y = f(x)`: convex `f` with `y` pushed down (or concave
//!   with `y` pushed up) is the pure LP epigraph, no integers (Williams, ch. 9;
//!   Vielma, Ahmed, Nemhauser, "Mixed-integer models for nonseparable
//!   piecewise-linear optimization", Oper. Res. 58 (2010), Sec. 1). Any other
//!   `f` uses the incremental formulation of the same paper (Sec. 3, "Inc").
//! - `abs`, `max`, `min`: the convex direction is a pure LP (`r >= a_i`);
//!   otherwise one selector binary per argument and big-M rows (Williams,
//!   ch. 9). The direction test is "decreasing `r` never hurts" read off the
//!   objective and rows, and `r` must not appear in another object.
//! - `and`, `or`, `not` over binaries: exact rows (Williams, ch. 9).

use crate::gc::{FuncKind, GcModel, SosKind};
use crate::model::RowSense;

/// Options of the linearization.
#[derive(Debug, Clone, Copy)]
pub struct LinOptions {
    /// Largest big-M a derived row may carry (`gc_bigm_max`).
    pub bigm_max: f64,
    /// Domain-propagation passes over the rows used to tighten the column
    /// box before computing big-Ms (`0` = the model's own bounds).
    pub propagation_passes: usize,
}

impl Default for LinOptions {
    fn default() -> Self {
        Self {
            bigm_max: 1e6,
            propagation_passes: 3,
        }
    }
}

/// What the linearization reads of the user's model (the LP part).
pub struct LinInput<'a> {
    pub col_lower: &'a [f64],
    pub col_upper: &'a [f64],
    /// Objective coefficients as the user wrote them.
    pub col_cost: &'a [f64],
    pub maximize: bool,
    /// Rows in CSR: `row_start.len() == nrow + 1`.
    pub row_start: &'a [usize],
    pub row_index: &'a [usize],
    pub row_value: &'a [f64],
    pub row_lower: &'a [f64],
    pub row_upper: &'a [f64],
    pub col_names: &'a [String],
}

/// A column the linearization adds (index `n + position`).
#[derive(Debug, Clone, PartialEq)]
pub struct NewCol {
    pub name: String,
    pub lower: f64,
    pub upper: f64,
    pub integer: bool,
}

/// A row the linearization adds: `lower <= terms <= upper`.
#[derive(Debug, Clone, PartialEq)]
pub struct NewRow {
    pub name: String,
    pub terms: Vec<(usize, f64)>,
    pub lower: f64,
    pub upper: f64,
}

/// How one object was linearized, for the envelope.
#[derive(Debug, Clone, PartialEq)]
pub struct Route {
    pub object: String,
    /// `big_m`, `redundant`, `sos1_binary`, `sos2_binary`, `epigraph_lp`,
    /// `incremental`, `abs_lp`, `abs_binary`, `max_lp`, `max_binary`,
    /// `min_lp`, `min_binary`, `logic_rows`.
    pub route: &'static str,
    /// Largest big-M used by the object (0 when none).
    pub max_big_m: f64,
}

/// The derived model part.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Extension {
    pub cols: Vec<NewCol>,
    pub rows: Vec<NewRow>,
    /// Existing columns whose box is intersected with `[lower, upper]`
    /// (the domain of a piecewise-linear `x`).
    pub bound_changes: Vec<(usize, f64, f64)>,
    pub routes: Vec<Route>,
}

/// A typed refusal naming the object.
#[derive(Debug, Clone, PartialEq)]
pub struct GcError {
    pub object: String,
    pub message: String,
}

impl std::fmt::Display for GcError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.object, self.message)
    }
}

impl std::error::Error for GcError {}

/// CSR rows of a CSC matrix (a counting transpose).
pub fn csr_from_csc(
    nrow: usize,
    col_start: &[usize],
    col_index: &[usize],
    col_value: &[f64],
) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
    let mut start = vec![0usize; nrow + 1];
    for &i in col_index {
        start[i + 1] += 1;
    }
    for i in 0..nrow {
        start[i + 1] += start[i];
    }
    let mut next = start.clone();
    let mut index = vec![0usize; col_index.len()];
    let mut value = vec![0.0f64; col_index.len()];
    for j in 0..col_start.len().saturating_sub(1) {
        for k in col_start[j]..col_start[j + 1] {
            let i = col_index[k];
            index[next[i]] = j;
            value[next[i]] = col_value[k];
            next[i] += 1;
        }
    }
    (start, index, value)
}

const INF: f64 = f64::INFINITY;
/// A value this far from a bound counts as "the row is redundant".
const TINY: f64 = 1e-9;

struct Lin<'a> {
    inp: &'a LinInput<'a>,
    opt: LinOptions,
    n: usize,
    /// Effective box (propagated; semi columns include 0).
    lo: Vec<f64>,
    hi: Vec<f64>,
    /// How many gc objects mention each column.
    uses: Vec<u32>,
    /// Column-wise view of the rows, built on first use.
    by_col: Option<Vec<Vec<(usize, f64)>>>,
    ext: Extension,
}

fn err<T>(object: impl Into<String>, message: impl Into<String>) -> Result<T, GcError> {
    Err(GcError {
        object: object.into(),
        message: message.into(),
    })
}

/// Linearizes every general constraint of `gc` (semi-continuous columns are
/// left to the solver, which handles them natively).
pub fn linearize(gc: &GcModel, inp: &LinInput<'_>, opt: LinOptions) -> Result<Extension, GcError> {
    let n = inp.col_lower.len();
    let mut uses = vec![0u32; n];
    for o in &gc.indicators {
        uses[o.bin] += 1;
        for &(j, _) in &o.terms {
            uses[j] += 1;
        }
    }
    for s in &gc.sos {
        for &(j, _) in &s.members {
            uses[j] += 1;
        }
    }
    for p in &gc.pwl {
        uses[p.x] += 1;
        uses[p.y] += 1;
    }
    for f in &gc.funcs {
        uses[f.res] += 1;
        for &a in &f.args {
            uses[a] += 1;
        }
    }
    let mut lin = Lin {
        inp,
        opt,
        n,
        lo: inp.col_lower.to_vec(),
        hi: inp.col_upper.to_vec(),
        uses,
        by_col: None,
        ext: Extension::default(),
    };
    for s in &gc.semi {
        lin.lo[s.col] = lin.lo[s.col].min(0.0);
        lin.hi[s.col] = lin.hi[s.col].max(0.0);
    }
    lin.propagate(gc);
    for o in &gc.indicators {
        lin.indicator(o)?;
    }
    for s in &gc.sos {
        lin.sos(s)?;
    }
    for p in &gc.pwl {
        lin.pwl(p)?;
    }
    for f in &gc.funcs {
        lin.func(f)?;
    }
    Ok(lin.ext)
}

impl Lin<'_> {
    fn name(&self, j: usize) -> &str {
        self.inp.col_names.get(j).map_or("", String::as_str)
    }

    fn add_col(&mut self, name: String, lower: f64, upper: f64, integer: bool) -> usize {
        self.ext.cols.push(NewCol {
            name,
            lower,
            upper,
            integer,
        });
        self.n + self.ext.cols.len() - 1
    }

    fn add_row(&mut self, name: String, terms: Vec<(usize, f64)>, lower: f64, upper: f64) {
        self.ext.rows.push(NewRow {
            name,
            terms,
            lower,
            upper,
        });
    }

    fn route(&mut self, object: &str, route: &'static str, max_big_m: f64) {
        self.ext.routes.push(Route {
            object: object.to_string(),
            route,
            max_big_m,
        });
    }

    /// Bound propagation of the regular rows over the effective box. The
    /// tightened box holds for every feasible point, so it is a valid source
    /// for big-Ms. Semi-continuous columns and integer rounding are left
    /// alone (a semi box already contains 0).
    fn propagate(&mut self, gc: &GcModel) {
        let inp = self.inp;
        let nrow = inp.row_lower.len();
        let semi: Vec<bool> = {
            let mut v = vec![false; self.n];
            for s in &gc.semi {
                v[s.col] = true;
            }
            v
        };
        for _ in 0..self.opt.propagation_passes {
            let mut changed = false;
            for i in 0..nrow {
                let (a, b) = (inp.row_start[i], inp.row_start[i + 1]);
                // activity bounds with the count of infinite contributions
                let (mut min_f, mut max_f) = (0.0f64, 0.0f64);
                let (mut min_inf, mut max_inf) = (0usize, 0usize);
                for k in a..b {
                    let (j, c) = (inp.row_index[k], inp.row_value[k]);
                    let (lo, hi) = (self.lo[j], self.hi[j]);
                    let (cl, ch) = if c > 0.0 { (lo, hi) } else { (hi, lo) };
                    if cl.is_infinite() {
                        min_inf += 1;
                    } else {
                        min_f += c * cl;
                    }
                    if ch.is_infinite() {
                        max_inf += 1;
                    } else {
                        max_f += c * ch;
                    }
                }
                let (rl, ru) = (inp.row_lower[i], inp.row_upper[i]);
                for k in a..b {
                    let (j, c) = (inp.row_index[k], inp.row_value[k]);
                    if semi[j] || c == 0.0 {
                        continue;
                    }
                    let (lo, hi) = (self.lo[j], self.hi[j]);
                    let (cl, ch) = if c > 0.0 { (lo, hi) } else { (hi, lo) };
                    // residual min/max activity without column j
                    let rest_min = if cl.is_infinite() {
                        if min_inf == 1 { Some(min_f) } else { None }
                    } else if min_inf == 0 {
                        Some(min_f - c * cl)
                    } else {
                        None
                    };
                    let rest_max = if ch.is_infinite() {
                        if max_inf == 1 { Some(max_f) } else { None }
                    } else if max_inf == 0 {
                        Some(max_f - c * ch)
                    } else {
                        None
                    };
                    // c x <= ru - rest_min ; c x >= rl - rest_max
                    let mut nl = lo;
                    let mut nh = hi;
                    if ru.is_finite()
                        && let Some(rm) = rest_min
                    {
                        let bound = (ru - rm) / c;
                        if c > 0.0 {
                            nh = nh.min(bound);
                        } else {
                            nl = nl.max(bound);
                        }
                    }
                    if rl.is_finite()
                        && let Some(rm) = rest_max
                    {
                        let bound = (rl - rm) / c;
                        if c > 0.0 {
                            nl = nl.max(bound);
                        } else {
                            nh = nh.min(bound);
                        }
                    }
                    // a tightening below the noise floor is not worth a pass
                    if nl > lo + 1e-6 && nl <= nh + 1e-9 {
                        self.lo[j] = nl;
                        changed = true;
                    }
                    if nh < hi - 1e-6 && nh >= nl - 1e-9 {
                        self.hi[j] = nh;
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }
    }

    /// Range of `sum(terms)` over the box.
    fn activity(&self, terms: &[(usize, f64)]) -> (f64, f64) {
        let (mut mn, mut mx) = (0.0f64, 0.0f64);
        for &(j, c) in terms {
            if c == 0.0 {
                continue;
            }
            let (lo, hi) = (self.lo[j], self.hi[j]);
            let (cl, ch) = if c > 0.0 { (lo, hi) } else { (hi, lo) };
            mn += if cl.is_infinite() {
                -c.abs() * INF
            } else {
                c * cl
            };
            mx += if ch.is_infinite() {
                c.abs() * INF
            } else {
                c * ch
            };
        }
        (mn, mx)
    }

    fn merged(terms: &[(usize, f64)]) -> Vec<(usize, f64)> {
        let mut t = terms.to_vec();
        t.sort_by_key(|&(j, _)| j);
        let mut out: Vec<(usize, f64)> = Vec::with_capacity(t.len());
        for (j, c) in t {
            match out.last_mut() {
                Some(l) if l.0 == j => l.1 += c,
                _ => out.push((j, c)),
            }
        }
        out.retain(|&(_, c)| c != 0.0);
        out
    }

    fn big_m(&self, object: &str, m: f64, what: &str) -> Result<f64, GcError> {
        if !m.is_finite() {
            return err(
                object,
                format!("{what} is unbounded, so no big-M exists (give the columns finite bounds)"),
            );
        }
        if m > self.opt.bigm_max {
            return err(
                object,
                format!(
                    "big-M {m:.6e} for {what} exceeds gc_bigm_max {:.6e} (tighten the bounds or raise gc_bigm_max)",
                    self.opt.bigm_max
                ),
            );
        }
        Ok(m)
    }

    // ----------------------------------------------------------- indicator

    fn indicator(&mut self, o: &crate::gc::Indicator) -> Result<(), GcError> {
        let object = format!("indicator '{}'", o.name);
        let terms = Self::merged(&o.terms);
        let (mn, mx) = self.activity(&terms);
        let mut max_m = 0.0f64;
        let mut parts = 0;
        // z has coefficient `s * M`: active when z == bin_val.
        let sign = if o.bin_val { 1.0 } else { -1.0 };
        let shift = |m: f64| if o.bin_val { m } else { 0.0 };
        if matches!(o.sense, RowSense::LessEqual | RowSense::Equal) {
            let m = mx - o.rhs;
            if m > TINY {
                let m = self.big_m(&object, m, "the row activity")?;
                max_m = max_m.max(m);
                let mut t = terms.clone();
                t.push((o.bin, sign * m));
                self.add_row(format!("{}_le", o.name), t, -INF, o.rhs + shift(m));
                parts += 1;
            }
        }
        if matches!(o.sense, RowSense::GreaterEqual | RowSense::Equal) {
            let m = o.rhs - mn;
            if m > TINY {
                let m = self.big_m(&object, m, "the row activity")?;
                max_m = max_m.max(m);
                let mut t = terms.clone();
                t.push((o.bin, -sign * m));
                self.add_row(format!("{}_ge", o.name), t, o.rhs - shift(m), INF);
                parts += 1;
            }
        }
        self.route(
            &object,
            if parts == 0 { "redundant" } else { "big_m" },
            max_m,
        );
        Ok(())
    }

    // ----------------------------------------------------------------- SOS

    fn sos(&mut self, s: &crate::gc::SosSet) -> Result<(), GcError> {
        let object = format!("sos '{}'", s.name);
        let mut members = s.members.clone();
        members.sort_by(|a, b| a.1.total_cmp(&b.1));
        for &(j, _) in &members {
            if !self.lo[j].is_finite() || !self.hi[j].is_finite() {
                return err(
                    &object,
                    format!(
                        "member '{}' has an infinite bound, so the set cannot be linearized (give it finite bounds)",
                        self.name(j)
                    ),
                );
            }
        }
        let m = members.len();
        if m <= 1 {
            self.route(&object, "redundant", 0.0);
            return Ok(());
        }
        // selectors: one per member (SOS1), one per adjacent pair (SOS2)
        let nsel = if s.kind == SosKind::Sos1 { m } else { m - 1 };
        let sel: Vec<usize> = (0..nsel)
            .map(|k| self.add_col(format!("{}_sel{k}", s.name), 0.0, 1.0, true))
            .collect();
        self.add_row(
            format!("{}_one", s.name),
            sel.iter().map(|&c| (c, 1.0)).collect(),
            -INF,
            1.0,
        );
        let mut max_m = 0.0f64;
        for (k, &(j, _)) in members.iter().enumerate() {
            // selectors that let member k be nonzero
            let owners: Vec<usize> = match s.kind {
                SosKind::Sos1 => vec![sel[k]],
                SosKind::Sos2 => {
                    let mut v = Vec::new();
                    if k > 0 {
                        v.push(sel[k - 1]);
                    }
                    if k < m - 1 {
                        v.push(sel[k]);
                    }
                    v
                }
            };
            let (lo, hi) = (self.lo[j].min(0.0), self.hi[j].max(0.0));
            let (blo, bhi) = (self.lo[j], self.hi[j]);
            // x <= hi * (sum owners); x >= lo * (sum owners); when the bound
            // is zero the row is just x <= 0 / x >= 0.
            let bm = self.big_m(&object, hi.max(-lo), "a member's bound")?;
            max_m = max_m.max(bm);
            if bhi > 0.0 || blo > 0.0 {
                let mut t = vec![(j, 1.0)];
                t.extend(owners.iter().map(|&c| (c, -hi)));
                self.add_row(format!("{}_u{k}", s.name), t, -INF, 0.0);
            } else {
                self.add_row(format!("{}_u{k}", s.name), vec![(j, 1.0)], -INF, 0.0);
            }
            if blo < 0.0 {
                let mut t = vec![(j, 1.0)];
                t.extend(owners.iter().map(|&c| (c, -lo)));
                self.add_row(format!("{}_l{k}", s.name), t, 0.0, INF);
            } else {
                self.add_row(format!("{}_l{k}", s.name), vec![(j, 1.0)], 0.0, INF);
            }
        }
        self.route(
            &object,
            if s.kind == SosKind::Sos1 {
                "sos1_binary"
            } else {
                "sos2_binary"
            },
            max_m,
        );
        Ok(())
    }

    // ------------------------------------------------------------ pressure

    /// Rows of every column, built once.
    fn rows_of(&mut self, j: usize) -> &[(usize, f64)] {
        if self.by_col.is_none() {
            let inp = self.inp;
            let mut v: Vec<Vec<(usize, f64)>> = vec![Vec::new(); self.n];
            for i in 0..inp.row_lower.len() {
                for k in inp.row_start[i]..inp.row_start[i + 1] {
                    v[inp.row_index[k]].push((i, inp.row_value[k]));
                }
            }
            self.by_col = Some(v);
        }
        &self.by_col.as_ref().expect("built")[j]
    }

    /// `down = true`: lowering column `r` never costs objective or feasibility
    /// (and strictly gains objective), so a relaxation `r >= f(..)` is tight.
    /// `down = false`: the mirror image. `r` must not appear in another object.
    fn pushed(&mut self, r: usize, down: bool) -> bool {
        if self.uses[r] > 1 {
            return false;
        }
        let c = self.inp.col_cost[r];
        // minimize: lowering r gains when c > 0; maximize: when c < 0
        let gain_down = if self.inp.maximize { -c } else { c };
        let gain = if down { gain_down } else { -gain_down };
        if gain <= 1e-12 {
            return false;
        }
        let rows: Vec<(usize, f64)> = self.rows_of(r).to_vec();
        for (i, a) in rows {
            let (lo, hi) = (self.inp.row_lower[i], self.inp.row_upper[i]);
            // moving r by `+d`; `down` means d < 0
            let ok = if down == (a > 0.0) {
                // activity falls: only a lower row bound can resist
                lo == -INF
            } else {
                hi == INF
            };
            if !ok {
                return false;
            }
        }
        true
    }

    // ----------------------------------------------------------------- PWL

    fn pwl(&mut self, p: &crate::gc::Pwl) -> Result<(), GcError> {
        let object = format!("pwl '{}'", p.name);
        let (xs, ys) = (&p.xs, &p.ys);
        let np = xs.len();
        let seg = np - 1;
        self.ext.bound_changes.push((p.x, xs[0], xs[np - 1]));
        self.lo[p.x] = self.lo[p.x].max(xs[0]);
        self.hi[p.x] = self.hi[p.x].min(xs[np - 1]);
        let ymin = ys.iter().copied().fold(INF, f64::min);
        let ymax = ys.iter().copied().fold(-INF, f64::max);
        let jump = xs.windows(2).any(|w| w[0] == w[1]);
        // slopes (no jumps)
        let slopes: Vec<f64> = if jump {
            Vec::new()
        } else {
            (0..seg)
                .map(|k| (ys[k + 1] - ys[k]) / (xs[k + 1] - xs[k]))
                .collect()
        };
        let convex = !jump && slopes.windows(2).all(|w| w[1] >= w[0] - 1e-12);
        let concave = !jump && slopes.windows(2).all(|w| w[1] <= w[0] + 1e-12);
        let lines = |me: &mut Self, ge: bool| {
            // y >= (or <=) y_k + s_k (x - x_k) for each segment
            for k in 0..seg {
                let s = slopes[k];
                let c = ys[k] - s * xs[k];
                // y - s x >= c
                me.add_row(
                    format!("{}_seg{k}", p.name),
                    vec![(p.y, 1.0), (p.x, -s)],
                    if ge { c } else { -INF },
                    if ge { INF } else { c },
                );
            }
        };
        if convex && self.lo[p.y] <= ymin + 1e-12 && self.pushed(p.y, true) {
            lines(self, true);
            self.route(&object, "epigraph_lp", 0.0);
            return Ok(());
        }
        if concave && self.hi[p.y] >= ymax - 1e-12 && self.pushed(p.y, false) {
            lines(self, false);
            self.route(&object, "epigraph_lp", 0.0);
            return Ok(());
        }
        // incremental formulation: x = x_0 + sum dx_k d_k, y likewise,
        // d_{k+1} <= z_k <= d_k with z binary
        let d: Vec<usize> = (0..seg)
            .map(|k| self.add_col(format!("{}_d{k}", p.name), 0.0, 1.0, false))
            .collect();
        let z: Vec<usize> = (0..seg.saturating_sub(1))
            .map(|k| self.add_col(format!("{}_z{k}", p.name), 0.0, 1.0, true))
            .collect();
        let mut tx = vec![(p.x, 1.0)];
        let mut ty = vec![(p.y, 1.0)];
        for k in 0..seg {
            tx.push((d[k], -(xs[k + 1] - xs[k])));
            ty.push((d[k], -(ys[k + 1] - ys[k])));
        }
        self.add_row(format!("{}_x", p.name), tx, xs[0], xs[0]);
        self.add_row(format!("{}_y", p.name), ty, ys[0], ys[0]);
        for k in 0..seg.saturating_sub(1) {
            // d_{k+1} - z_k <= 0 ; z_k - d_k <= 0
            self.add_row(
                format!("{}_o{k}", p.name),
                vec![(d[k + 1], 1.0), (z[k], -1.0)],
                -INF,
                0.0,
            );
            self.add_row(
                format!("{}_p{k}", p.name),
                vec![(z[k], 1.0), (d[k], -1.0)],
                -INF,
                0.0,
            );
        }
        self.route(&object, "incremental", 0.0);
        Ok(())
    }

    // ------------------------------------------------------------ functions

    fn func(&mut self, f: &crate::gc::FuncCon) -> Result<(), GcError> {
        let object = format!("{} '{}'", f.kind.word(), f.name);
        let r = f.res;
        match f.kind {
            FuncKind::And | FuncKind::Or => {
                let k = f.args.len() as f64;
                let and = f.kind == FuncKind::And;
                for (i, &a) in f.args.iter().enumerate() {
                    // and: r <= a_i ; or: r >= a_i
                    self.add_row(
                        format!("{}_{i}", f.name),
                        vec![(r, 1.0), (a, -1.0)],
                        if and { -INF } else { 0.0 },
                        if and { 0.0 } else { INF },
                    );
                }
                // and: r >= sum a - (k - 1) ; or: r <= sum a
                let mut t: Vec<(usize, f64)> = f.args.iter().map(|&a| (a, -1.0)).collect();
                t.push((r, 1.0));
                if and {
                    self.add_row(format!("{}_all", f.name), t, 1.0 - k, INF);
                } else {
                    self.add_row(format!("{}_any", f.name), t, -INF, 0.0);
                }
                self.route(&object, "logic_rows", 0.0);
            }
            FuncKind::Not => {
                // r + a = 1
                self.add_row(f.name.clone(), vec![(r, 1.0), (f.args[0], 1.0)], 1.0, 1.0);
                self.route(&object, "logic_rows", 0.0);
            }
            FuncKind::Abs => self.abs(f, &object)?,
            FuncKind::Max | FuncKind::Min => self.minmax(f, &object)?,
        }
        Ok(())
    }

    fn abs(&mut self, f: &crate::gc::FuncCon, object: &str) -> Result<(), GcError> {
        let (r, x) = (f.res, f.args[0]);
        let (l, u) = (self.lo[x], self.hi[x]);
        if l >= 0.0 {
            self.add_row(f.name.clone(), vec![(r, 1.0), (x, -1.0)], 0.0, 0.0);
            self.route(object, "abs_lp", 0.0);
            return Ok(());
        }
        if u <= 0.0 {
            self.add_row(f.name.clone(), vec![(r, 1.0), (x, 1.0)], 0.0, 0.0);
            self.route(object, "abs_lp", 0.0);
            return Ok(());
        }
        let rlo = 0.0;
        if self.lo[r] <= rlo && self.pushed(r, true) {
            // r >= x and r >= -x
            self.add_row(format!("{}_p", f.name), vec![(r, 1.0), (x, -1.0)], 0.0, INF);
            self.add_row(format!("{}_n", f.name), vec![(r, 1.0), (x, 1.0)], 0.0, INF);
            self.route(object, "abs_lp", 0.0);
            return Ok(());
        }
        let m = self.big_m(object, (-l).max(u) * 2.0, "the argument range")?;
        // s = 1: x >= 0 and r = x ; s = 0: x <= 0 and r = -x
        let s = self.add_col(format!("{}_s", f.name), 0.0, 1.0, true);
        self.add_row(format!("{}_p", f.name), vec![(r, 1.0), (x, -1.0)], 0.0, INF);
        self.add_row(format!("{}_n", f.name), vec![(r, 1.0), (x, 1.0)], 0.0, INF);
        // x >= l (1 - s)  ->  x + l s >= l
        self.add_row(format!("{}_xl", f.name), vec![(x, 1.0), (s, l)], l, INF);
        // x <= u s  ->  x - u s <= 0
        self.add_row(format!("{}_xu", f.name), vec![(x, 1.0), (s, -u)], -INF, 0.0);
        // r <= x - 2 l (1 - s)  ->  r - x - 2 l s <= -2 l
        self.add_row(
            format!("{}_rp", f.name),
            vec![(r, 1.0), (x, -1.0), (s, -2.0 * l)],
            -INF,
            -2.0 * l,
        );
        // r <= -x + 2 u s  ->  r + x - 2 u s <= 0
        self.add_row(
            format!("{}_rn", f.name),
            vec![(r, 1.0), (x, 1.0), (s, -2.0 * u)],
            -INF,
            0.0,
        );
        self.route(object, "abs_binary", m);
        Ok(())
    }

    fn minmax(&mut self, f: &crate::gc::FuncCon, object: &str) -> Result<(), GcError> {
        let is_max = f.kind == FuncKind::Max;
        let r = f.res;
        // arguments as (column or constant, lo, hi)
        let mut args: Vec<(Option<usize>, f64, f64)> = f
            .args
            .iter()
            .map(|&a| (Some(a), self.lo[a], self.hi[a]))
            .collect();
        if let Some(c) = f.constant {
            args.push((None, c, c));
        }
        // value range of the function over the box
        let (vlo, vhi) = if is_max {
            (
                args.iter().map(|a| a.1).fold(-INF, f64::max),
                args.iter().map(|a| a.2).fold(-INF, f64::max),
            )
        } else {
            (
                args.iter().map(|a| a.1).fold(INF, f64::min),
                args.iter().map(|a| a.2).fold(INF, f64::min),
            )
        };
        // max is convex (relax with r >= a_i, r pushed down); min is concave
        let tag = if is_max { "max" } else { "min" };
        let lp_ok = if is_max {
            self.lo[r] <= vlo && self.pushed(r, true)
        } else {
            self.hi[r] >= vhi && self.pushed(r, false)
        };
        let ge = is_max;
        let one_side = |me: &mut Self, i: usize, a: &(Option<usize>, f64, f64)| match a.0 {
            Some(c) => me.add_row(
                format!("{}_{i}", f.name),
                vec![(r, 1.0), (c, -1.0)],
                if ge { 0.0 } else { -INF },
                if ge { INF } else { 0.0 },
            ),
            None => me.add_row(
                format!("{}_{i}", f.name),
                vec![(r, 1.0)],
                if ge { a.1 } else { -INF },
                if ge { INF } else { a.1 },
            ),
        };
        if lp_ok {
            for (i, a) in args.iter().enumerate() {
                one_side(self, i, a);
            }
            self.route(object, if is_max { "max_lp" } else { "min_lp" }, 0.0);
            return Ok(());
        }
        // binary route: r >= a_i (max) for all i, and r <= a_i + M_i (1 - z_i)
        // for the selected one; mirrored for min.
        let mut max_m = 0.0f64;
        let mut zs = Vec::new();
        for (i, a) in args.iter().enumerate() {
            one_side(self, i, a);
            // M_i bounds how far the function's value can be from a_i:
            // max: r - a_i <= vhi - a.1 ; min: a_i - r <= a.2 - vlo
            let m = if is_max { vhi - a.1 } else { a.2 - vlo };
            let m = self.big_m(object, m, &format!("the {tag} gap of argument {i}"))?;
            if m <= TINY {
                // the argument always equals the extreme: no selector needed
                zs.push(None);
                continue;
            }
            max_m = max_m.max(m);
            let z = self.add_col(format!("{}_z{i}", f.name), 0.0, 1.0, true);
            zs.push(Some(z));
            // max: r <= a_i + M (1 - z)  ->  r - a_i + M z <= M (+ c)
            // min: r >= a_i - M (1 - z)  ->  r - a_i - M z >= -M (+ c)
            let mut t = vec![(r, 1.0)];
            let mut c_arg = 0.0;
            match a.0 {
                Some(c) => t.push((c, -1.0)),
                None => c_arg = a.1,
            }
            if is_max {
                t.push((z, m));
                self.add_row(format!("{}_s{i}", f.name), t, -INF, m + c_arg);
            } else {
                t.push((z, -m));
                self.add_row(format!("{}_s{i}", f.name), t, c_arg - m, INF);
            }
        }
        // some argument is selected (an always-extreme argument needs none)
        if zs.iter().all(Option::is_some) {
            self.add_row(
                format!("{}_pick", f.name),
                zs.iter().flatten().map(|&z| (z, 1.0)).collect(),
                1.0,
                INF,
            );
        }
        self.route(
            object,
            if is_max { "max_binary" } else { "min_binary" },
            max_m,
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gc::{FuncCon, Indicator, Pwl, SosSet};

    fn input<'a>(
        lo: &'a [f64],
        hi: &'a [f64],
        cost: &'a [f64],
        names: &'a [String],
    ) -> LinInput<'a> {
        LinInput {
            col_lower: lo,
            col_upper: hi,
            col_cost: cost,
            maximize: false,
            row_start: &[0],
            row_index: &[],
            row_value: &[],
            row_lower: &[],
            row_upper: &[],
            col_names: names,
        }
    }

    #[test]
    fn indicator_big_m_from_the_box() {
        let names: Vec<String> = ["x", "y", "b"].iter().map(|s| s.to_string()).collect();
        let (lo, hi, cost) = ([0.0, 0.0, 0.0], [10.0, 5.0, 1.0], [0.0; 3]);
        let mut gc = GcModel::default();
        gc.indicators.push(Indicator {
            name: "i".into(),
            bin: 2,
            bin_val: true,
            terms: vec![(0, 1.0), (1, 2.0)],
            sense: RowSense::LessEqual,
            rhs: 4.0,
        });
        let e = linearize(&gc, &input(&lo, &hi, &cost, &names), LinOptions::default()).unwrap();
        // max activity 20, rhs 4: M = 16, row x + 2y + 16 b <= 20
        assert_eq!(e.rows.len(), 1);
        assert_eq!(e.rows[0].terms, vec![(0, 1.0), (1, 2.0), (2, 16.0)]);
        assert_eq!(e.rows[0].upper, 20.0);
        // unbounded activity is a typed error naming the object
        let hi_inf = [INF, 5.0, 1.0];
        let r = linearize(
            &gc,
            &input(&lo, &hi_inf, &cost, &names),
            LinOptions::default(),
        );
        let er = r.unwrap_err();
        assert!(er.object.contains("indicator 'i'") && er.message.contains("unbounded"));
        // cap
        let o = LinOptions {
            bigm_max: 10.0,
            ..Default::default()
        };
        let er = linearize(&gc, &input(&lo, &hi, &cost, &names), o).unwrap_err();
        assert!(er.message.contains("gc_bigm_max"));
    }

    #[test]
    fn pwl_routes() {
        let names: Vec<String> = ["x", "y"].iter().map(|s| s.to_string()).collect();
        let (lo, hi) = ([0.0, 0.0], [4.0, 9.0]);
        let mut gc = GcModel::default();
        gc.pwl.push(Pwl {
            name: "p".into(),
            x: 0,
            y: 1,
            xs: vec![0.0, 2.0, 4.0],
            ys: vec![0.0, 1.0, 4.0],
        });
        // convex, y pushed down by a positive cost: epigraph
        let e = linearize(
            &gc,
            &input(&lo, &hi, &[0.0, 1.0], &names),
            LinOptions::default(),
        )
        .unwrap();
        assert_eq!(e.routes[0].route, "epigraph_lp");
        assert!(e.cols.is_empty());
        // no pressure: incremental
        let e = linearize(
            &gc,
            &input(&lo, &hi, &[0.0, 0.0], &names),
            LinOptions::default(),
        )
        .unwrap();
        assert_eq!(e.routes[0].route, "incremental");
        assert_eq!(e.cols.len(), 3);
    }

    #[test]
    fn sos_and_functions_name_their_failures() {
        let names: Vec<String> = ["x", "y", "r"].iter().map(|s| s.to_string()).collect();
        let (lo, hi, cost) = ([0.0; 3], [INF, 1.0, 5.0], [0.0; 3]);
        let mut gc = GcModel::default();
        gc.sos.push(SosSet {
            name: "s".into(),
            kind: SosKind::Sos1,
            priority: 0.0,
            members: vec![(0, 1.0), (1, 2.0)],
        });
        let er =
            linearize(&gc, &input(&lo, &hi, &cost, &names), LinOptions::default()).unwrap_err();
        assert!(
            er.object == "sos 's'" && er.message.contains("infinite bound"),
            "{er}"
        );
        let mut gc = GcModel::default();
        gc.funcs.push(FuncCon {
            name: "m".into(),
            kind: FuncKind::Max,
            res: 2,
            args: vec![0, 1],
            constant: None,
        });
        let hi2 = [3.0, 1.0, 5.0];
        let e = linearize(&gc, &input(&lo, &hi2, &cost, &names), LinOptions::default()).unwrap();
        assert_eq!(e.routes[0].route, "max_binary");
    }
}
