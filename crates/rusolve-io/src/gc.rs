//! Semi-continuous columns and general constraints of the shared [`crate::Model`].
//!
//! The model stores the *user's* objects. Nothing here is a linearization: the
//! envelope, `rusolve check` and `rusolve explain` report on these objects,
//! and the linearization path (`rusolve-io/src/gc_lin.rs`) derives rows from
//! them only when a solve asks for it.
//!
//! Object semantics (public documentation of the CPLEX and Gurobi model
//! formats and APIs; the grammars are cited in `lp.rs` and `mps.rs`):
//! - semi-continuous `x`: `x = 0` or `lb <= x <= ub` (semi-integer: also integral);
//! - indicator `z = v -> a'x (sense) b`: the row holds when binary `z` equals `v`;
//! - SOS1: at most one member is nonzero; SOS2: at most two, and adjacent in
//!   weight order;
//! - piecewise-linear `y = f(x)`: `(x, y)` lies on the polyline through the
//!   points; `x` is restricted to `[x_first, x_last]` (no extrapolation: a
//!   documented difference from Gurobi's extension of the end segments, kept
//!   because the domain then comes from the object and needs no extra bound);
//!   repeated x values are a jump (discontinuity) and are allowed in the
//!   object and the checker, but not by the linearization;
//! - `r = abs(a)`, `r = max(a_1, .., a_k, [c])`, `r = min(..)`, `r = and(..)`,
//!   `r = or(..)` (binaries) and `r = not(a)`.

use crate::model::{RowSense, VariableBound};

/// A semi-continuous or semi-integer column. Its bounds are the model's
/// `bounds[col]`: `x = 0` or `lower <= x <= upper`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SemiVar {
    pub col: usize,
    /// Semi-integer: the nonzero branch is integral. `bounds[col].integer`
    /// stays `false` for a semi-integer column (the kind lives here).
    pub integer: bool,
}

/// `bin = bin_val -> sum(terms) sense rhs`.
#[derive(Debug, Clone, PartialEq)]
pub struct Indicator {
    pub name: String,
    pub bin: usize,
    /// The value of the binary that activates the row.
    pub bin_val: bool,
    pub terms: Vec<(usize, f64)>,
    pub sense: RowSense,
    pub rhs: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SosKind {
    Sos1,
    Sos2,
}

/// A special ordered set: members with their weights (the order key).
#[derive(Debug, Clone, PartialEq)]
pub struct SosSet {
    pub name: String,
    pub kind: SosKind,
    pub priority: f64,
    pub members: Vec<(usize, f64)>,
}

/// `y = f(x)` through the points `(xs[i], ys[i])`, `xs` non-decreasing.
#[derive(Debug, Clone, PartialEq)]
pub struct Pwl {
    pub name: String,
    pub x: usize,
    pub y: usize,
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FuncKind {
    Abs,
    Min,
    Max,
    And,
    Or,
    Not,
}

impl FuncKind {
    pub fn word(self) -> &'static str {
        match self {
            FuncKind::Abs => "abs",
            FuncKind::Min => "min",
            FuncKind::Max => "max",
            FuncKind::And => "and",
            FuncKind::Or => "or",
            FuncKind::Not => "not",
        }
    }
}

/// `res = kind(args.., [constant])`.
#[derive(Debug, Clone, PartialEq)]
pub struct FuncCon {
    pub name: String,
    pub kind: FuncKind,
    pub res: usize,
    pub args: Vec<usize>,
    /// The optional constant argument of `min` / `max`.
    pub constant: Option<f64>,
}

/// Everything beyond a plain MIP that a [`crate::Model`] can carry.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GcModel {
    pub semi: Vec<SemiVar>,
    pub indicators: Vec<Indicator>,
    pub sos: Vec<SosSet>,
    pub pwl: Vec<Pwl>,
    pub funcs: Vec<FuncCon>,
}

/// One violated object found by [`GcModel::check`].
#[derive(Debug, Clone, PartialEq)]
pub struct GcViolation {
    /// `semi_continuous`, `indicator`, `sos`, `pwl` or `function`.
    pub kind: &'static str,
    /// The object's name (the column name for a semi-continuous column).
    pub name: String,
    /// How far from satisfying it the point is.
    pub amount: f64,
    pub detail: String,
}

/// Tolerances of [`GcModel::check`].
#[derive(Debug, Clone, Copy)]
pub struct GcTol {
    /// Rounds a binary (and the integer part of a semi-integer column).
    pub integrality: f64,
    /// Rows, function values and `x = 0` of a semi-continuous column.
    pub feasibility: f64,
}

impl Default for GcTol {
    fn default() -> Self {
        Self {
            integrality: 1e-6,
            feasibility: 1e-6,
        }
    }
}

impl GcModel {
    pub fn is_empty(&self) -> bool {
        self.semi.is_empty() && !self.has_general()
    }

    /// Any object other than a semi-continuous column.
    pub fn has_general(&self) -> bool {
        !(self.indicators.is_empty()
            && self.sos.is_empty()
            && self.pwl.is_empty()
            && self.funcs.is_empty())
    }

    /// Per-kind object counts, for messages and the envelope.
    pub fn counts(&self) -> [(&'static str, usize); 5] {
        [
            ("semi_continuous", self.semi.len()),
            ("indicator", self.indicators.len()),
            ("sos", self.sos.len()),
            ("pwl", self.pwl.len()),
            ("function", self.funcs.len()),
        ]
    }

    /// `true` when column `j` is semi-continuous or semi-integer.
    pub fn is_semi(&self, j: usize) -> bool {
        self.semi.iter().any(|s| s.col == j)
    }

    /// Typed refusal for commands that cannot honour semi-continuous columns or
    /// general constraints (they would otherwise drop them and report a wrong
    /// answer). `who` is the command, e.g. `rusolve solve`.
    pub fn reject(&self, who: &str) -> Result<(), String> {
        if let Some(what) = self.first_general() {
            return Err(format!(
                "{who}: the model has {what}; general constraints are only supported by `rusolve mip` and `rusolve check`"
            ));
        }
        if !self.semi.is_empty() {
            return Err(format!(
                "{who}: the model has semi-continuous columns; they are only supported by `rusolve mip` and `rusolve check`"
            ));
        }
        Ok(())
    }

    /// The first object of the model, named, for typed "not supported here"
    /// errors: `indicator 'ind1'`, `sos 's1'`, ...
    pub fn first_general(&self) -> Option<String> {
        if let Some(o) = self.indicators.first() {
            return Some(format!("indicator '{}'", o.name));
        }
        if let Some(o) = self.sos.first() {
            return Some(format!("sos '{}'", o.name));
        }
        if let Some(o) = self.pwl.first() {
            return Some(format!("pwl '{}'", o.name));
        }
        if let Some(o) = self.funcs.first() {
            return Some(format!("{} '{}'", o.kind.word(), o.name));
        }
        None
    }

    /// Structural checks every reader runs once the model is complete: an
    /// indicator variable must be a binary, an SOS member must exist, a PWL
    /// has at least two points with non-decreasing `x`, a function has the
    /// arguments its kind needs. The error names the object.
    pub fn validate(&self, bounds: &[VariableBound], col_names: &[String]) -> Result<(), String> {
        let n = bounds.len();
        let col_ok = |what: &str, j: usize| -> Result<(), String> {
            if j >= n {
                Err(format!("{what}: column index {j} out of range"))
            } else {
                Ok(())
            }
        };
        let is_binary = |j: usize| {
            let b = &bounds[j];
            b.integer && b.lower >= 0.0 && b.upper <= 1.0
        };
        let nm = |j: usize| col_names.get(j).cloned().unwrap_or_else(|| format!("c{j}"));
        for s in &self.semi {
            col_ok("semi-continuous column", s.col)?;
            let b = &bounds[s.col];
            if b.lower > b.upper {
                return Err(format!(
                    "semi-continuous column '{}': lower bound {} exceeds upper bound {}",
                    nm(s.col),
                    b.lower,
                    b.upper
                ));
            }
        }
        for o in &self.indicators {
            col_ok(&format!("indicator '{}'", o.name), o.bin)?;
            if !is_binary(o.bin) {
                return Err(format!(
                    "indicator '{}': variable '{}' is not binary (needs integer bounds within [0, 1])",
                    o.name,
                    nm(o.bin)
                ));
            }
            for &(j, _) in &o.terms {
                col_ok(&format!("indicator '{}'", o.name), j)?;
            }
        }
        for s in &self.sos {
            if s.members.is_empty() {
                return Err(format!("sos '{}' has no members", s.name));
            }
            for &(j, _) in &s.members {
                col_ok(&format!("sos '{}'", s.name), j)?;
            }
        }
        for p in &self.pwl {
            col_ok(&format!("pwl '{}'", p.name), p.x)?;
            col_ok(&format!("pwl '{}'", p.name), p.y)?;
            if p.xs.len() != p.ys.len() || p.xs.len() < 2 {
                return Err(format!(
                    "pwl '{}' needs at least two points with matching x and y counts",
                    p.name
                ));
            }
            if p.xs.windows(2).any(|w| w[1] < w[0]) {
                return Err(format!("pwl '{}': x values must be non-decreasing", p.name));
            }
        }
        for f in &self.funcs {
            col_ok(&format!("{} '{}'", f.kind.word(), f.name), f.res)?;
            for &j in &f.args {
                col_ok(&format!("{} '{}'", f.kind.word(), f.name), j)?;
            }
            let (min_args, max_args) = match f.kind {
                FuncKind::Abs | FuncKind::Not => (1, 1),
                _ => (1, usize::MAX),
            };
            if f.args.len() < min_args || f.args.len() > max_args {
                return Err(format!(
                    "{} '{}' has {} arguments",
                    f.kind.word(),
                    f.name,
                    f.args.len()
                ));
            }
            if matches!(f.kind, FuncKind::And | FuncKind::Or | FuncKind::Not) {
                for &j in f.args.iter().chain(std::iter::once(&f.res)) {
                    if !is_binary(j) {
                        return Err(format!(
                            "{} '{}': variable '{}' is not binary",
                            f.kind.word(),
                            f.name,
                            nm(j)
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// Evaluates every object on the point `x` (all model columns).
    ///
    /// Indicator logic is evaluated on the *rounded* binary (a fractional
    /// binary is itself a violation), so a big-M row that let `z = 1e-7` buy a
    /// violation is caught here ("trickle flow"). `lower`/`upper` are the
    /// model's column bounds (they define the semi-continuous branch).
    pub fn check(
        &self,
        x: &[f64],
        lower: &[f64],
        upper: &[f64],
        col_names: &[String],
        tol: GcTol,
    ) -> Vec<GcViolation> {
        let name = |j: usize| -> String {
            match col_names.get(j) {
                Some(n) if !n.is_empty() => n.clone(),
                _ => format!("c{j}"),
            }
        };
        let mut out = Vec::new();
        for s in &self.semi {
            let v = x[s.col];
            let (lo, hi) = (lower[s.col], upper[s.col]);
            let on_zero = v.abs() <= tol.feasibility;
            let in_range = v >= lo - tol.feasibility && v <= hi + tol.feasibility;
            if !on_zero && !in_range {
                let amount = if v < lo { lo - v } else { v - hi }.min(v.abs());
                out.push(GcViolation {
                    kind: "semi_continuous",
                    name: name(s.col),
                    amount,
                    detail: format!(
                        "value {v} is neither 0 nor in [{lo}, {hi}] (not zero, not in range)"
                    ),
                });
            } else if s.integer && !on_zero {
                let frac = (v - v.round()).abs();
                if frac > tol.integrality {
                    out.push(GcViolation {
                        kind: "semi_continuous",
                        name: name(s.col),
                        amount: frac,
                        detail: format!("semi-integer value {v} is not integral"),
                    });
                }
            }
        }
        for o in &self.indicators {
            let z = x[o.bin];
            let zr = z.round();
            if (z - zr).abs() > tol.integrality {
                out.push(GcViolation {
                    kind: "indicator",
                    name: o.name.clone(),
                    amount: (z - zr).abs(),
                    detail: format!("binary {} = {z} is not integral", name(o.bin)),
                });
                continue;
            }
            if (zr == 1.0) != o.bin_val {
                continue;
            }
            let act: f64 = o.terms.iter().map(|&(j, a)| a * x[j]).sum();
            let viol = match o.sense {
                RowSense::LessEqual => act - o.rhs,
                RowSense::GreaterEqual => o.rhs - act,
                RowSense::Equal => (act - o.rhs).abs(),
            };
            if viol > tol.feasibility {
                out.push(GcViolation {
                    kind: "indicator",
                    name: o.name.clone(),
                    amount: viol,
                    detail: format!(
                        "{} = {} holds (rounded) but the row activity is {act}, rhs {} ({})",
                        name(o.bin),
                        o.bin_val as u8,
                        o.rhs,
                        sense_word(o.sense)
                    ),
                });
            }
        }
        for s in &self.sos {
            if let Some(amount) = sos_violation(s, x, tol.feasibility) {
                out.push(GcViolation {
                    kind: "sos",
                    name: s.name.clone(),
                    amount,
                    detail: format!(
                        "{} needs {}",
                        match s.kind {
                            SosKind::Sos1 => "SOS1: at most one nonzero member",
                            SosKind::Sos2 => "SOS2: at most two nonzero members, adjacent",
                        },
                        "fewer or adjacent nonzeros"
                    ),
                });
            }
        }
        for p in &self.pwl {
            let (xv, yv) = (x[p.x], x[p.y]);
            let d = pwl_distance(p, xv, yv);
            if d > tol.feasibility {
                out.push(GcViolation {
                    kind: "pwl",
                    name: p.name.clone(),
                    amount: d,
                    detail: format!(
                        "({}, {}) = ({xv}, {yv}) is {d} off the polyline",
                        name(p.x),
                        name(p.y)
                    ),
                });
            }
        }
        for f in &self.funcs {
            let (expected, bin_gap) = func_value(f, x, tol.integrality);
            let got = x[f.res];
            let d = (got - expected).abs().max(bin_gap);
            if d > tol.feasibility {
                out.push(GcViolation {
                    kind: "function",
                    name: f.name.clone(),
                    amount: d,
                    detail: format!(
                        "{} = {}(..) should be {expected}, is {got}",
                        name(f.res),
                        f.kind.word()
                    ),
                });
            }
        }
        out
    }
}

fn sense_word(s: RowSense) -> &'static str {
    match s {
        RowSense::LessEqual => "<=",
        RowSense::GreaterEqual => ">=",
        RowSense::Equal => "=",
    }
}

/// How much mass must go to make the set valid; `None` when it already is.
fn sos_violation(s: &SosSet, x: &[f64], tol: f64) -> Option<f64> {
    let mut m: Vec<(f64, f64)> = s.members.iter().map(|&(j, w)| (w, x[j])).collect();
    m.sort_by(|a, b| a.0.total_cmp(&b.0));
    let nz: Vec<usize> = (0..m.len()).filter(|&i| m[i].1.abs() > tol).collect();
    let window = match s.kind {
        SosKind::Sos1 => 1,
        SosKind::Sos2 => 2,
    };
    if nz.len() <= 1 || (nz.len() == window && nz[nz.len() - 1] - nz[0] == window - 1) {
        return None;
    }
    // Cheapest window of `window` consecutive members to keep.
    let total: f64 = nz.iter().map(|&i| m[i].1.abs()).sum();
    let mut best = f64::INFINITY;
    for start in 0..m.len().saturating_sub(window - 1).max(1) {
        let kept: f64 = (start..(start + window).min(m.len()))
            .map(|i| m[i].1.abs())
            .sum();
        best = best.min(total - kept);
    }
    Some(best.max(tol))
}

/// Distance of `(x, y)` from the polyline (0 on it; allows jumps).
fn pwl_distance(p: &Pwl, x: f64, y: f64) -> f64 {
    let n = p.xs.len();
    if n == 0 {
        return 0.0;
    }
    if x < p.xs[0] {
        return p.xs[0] - x;
    }
    if x > p.xs[n - 1] {
        return x - p.xs[n - 1];
    }
    let mut best = f64::INFINITY;
    for i in 0..n.saturating_sub(1).max(1) {
        let (x0, x1) = (p.xs[i], p.xs[(i + 1).min(n - 1)]);
        let (y0, y1) = (p.ys[i], p.ys[(i + 1).min(n - 1)]);
        if x < x0 || x > x1 {
            continue;
        }
        let f = if x1 > x0 {
            y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        } else {
            // a vertical piece: any y between the two values is on it
            if y >= y0.min(y1) && y <= y0.max(y1) {
                y
            } else if y < y0.min(y1) {
                y0.min(y1)
            } else {
                y0.max(y1)
            }
        };
        best = best.min((y - f).abs());
    }
    best
}

/// The value `res` must take and, for and/or/not, how far an argument is from
/// a binary.
fn func_value(f: &FuncCon, x: &[f64], int_tol: f64) -> (f64, f64) {
    let mut gap = 0.0f64;
    let mut bin = |j: usize| -> f64 {
        let v = x[j];
        let r = v.round();
        if (v - r).abs() > int_tol {
            gap = gap.max((v - r).abs());
        }
        r
    };
    let v = match f.kind {
        FuncKind::Abs => x[f.args[0]].abs(),
        FuncKind::Max => f
            .args
            .iter()
            .map(|&j| x[j])
            .chain(f.constant)
            .fold(f64::NEG_INFINITY, f64::max),
        FuncKind::Min => f
            .args
            .iter()
            .map(|&j| x[j])
            .chain(f.constant)
            .fold(f64::INFINITY, f64::min),
        FuncKind::And => {
            let a: Vec<f64> = f.args.iter().map(|&j| bin(j)).collect();
            f64::from(a.iter().all(|&v| v == 1.0))
        }
        FuncKind::Or => {
            let a: Vec<f64> = f.args.iter().map(|&j| bin(j)).collect();
            f64::from(a.contains(&1.0))
        }
        FuncKind::Not => 1.0 - bin(f.args[0]),
    };
    (v, gap)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(n: usize) -> Vec<String> {
        (0..n).map(|j| format!("x{j}")).collect()
    }

    #[test]
    fn semi_membership() {
        let g = GcModel {
            semi: vec![SemiVar {
                col: 0,
                integer: false,
            }],
            ..Default::default()
        };
        let (lo, hi) = ([2.0], [5.0]);
        let t = GcTol::default();
        let chk = |v: f64| g.check(&[v], &lo, &hi, &names(1), t).len();
        assert_eq!(chk(0.0), 0);
        assert_eq!(chk(2.0), 0);
        assert_eq!(chk(5.0), 0);
        assert_eq!(chk(1.0), 1);
        assert_eq!(chk(6.0), 1);
    }

    #[test]
    fn pwl_and_funcs() {
        let p = Pwl {
            name: "p".into(),
            x: 0,
            y: 1,
            xs: vec![0.0, 1.0, 3.0],
            ys: vec![0.0, 1.0, 5.0],
        };
        assert!(pwl_distance(&p, 2.0, 3.0) < 1e-12);
        assert!(pwl_distance(&p, 2.0, 4.0) > 0.9);
        assert!(pwl_distance(&p, 4.0, 5.0) > 0.9);
        let f = FuncCon {
            name: "m".into(),
            kind: FuncKind::Max,
            res: 2,
            args: vec![0, 1],
            constant: Some(3.0),
        };
        assert_eq!(func_value(&f, &[1.0, 2.0, 0.0], 1e-6).0, 3.0);
    }
}
