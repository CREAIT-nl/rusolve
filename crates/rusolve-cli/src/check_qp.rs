//! Solution checking for a convex QP: what
//! `rusolve check` adds when the model has a Hessian.
//!
//! * The objective is recomputed with Q: `offset + c'x + x'Qx/2` (`Model::hessian` holds both
//!   triangles, so `x'Qx/2` is the full sum over the triplets, halved).
//! * With the dual section of a HiGHS raw solution file (`# Dual solution values`) the KKT
//!   conditions of `min c'x + x'Qx/2, rl <= Ax <= ru, l <= x <= u` are measured:
//!   stationarity `c + Qx - A'y - d = 0` (HiGHS convention: `d` is the column dual, `y` the row
//!   dual), dual sign feasibility (`d`, `y` at most 0 above the active lower bound and so on),
//!   and complementarity. A maximisation is checked as the minimisation of the negated
//!   objective (`c`, `Q`, `y`, `d` all negated).
//! * The duality gap is primal minus the Wolfe dual value
//!   `offset - x'Qx/2 + sum(y_i * active row bound) + sum(d_j * active column bound)`
//!   (Nocedal & Wright, Numerical Optimization 2nd ed., sec. 16.7 for the QP dual; Boyd &
//!   Vandenberghe, Convex Optimization, sec. 5.2 for the Lagrange dual). With stationarity
//!   met the gap equals the summed complementarity.
//!
//! `optimal` needs stationarity, dual infeasibility and complementarity within `KKT_TOL`, and a gap
//! within `KKT_TOL * (1 + |p| + |d|)` plus the part the stationarity residual explains
//! (`sum |r_j x_j|`: gap = complementarity + r'x).
//!
//! Without duals only the objective is reported and the KKT verdict is `unknown`; the primal
//! feasibility verdict of `rusolve check` is untouched.

use rusolve_core::Triplet;
use rusolve_ekk::prelude::{HighsLp, ObjSense};
use serde::Serialize;

/// KKT residual tolerance, relative to `1 + max |term|` of each condition. 1e-5, not the 1e-7 of
/// HiGHS's feasibility tolerances: the raw solution file prints values to about seven decimals,
/// so duals read back from it carry rounding of that size (HS118: residual 2e-6 from rounding).
pub const KKT_TOL: f64 = 1e-5;

/// Dual values read from a raw solution file, in the model's column and row order.
pub struct Duals {
    pub col: Vec<f64>,
    pub row: Vec<f64>,
}

#[derive(Serialize, Debug)]
pub struct QpCheck {
    /// `offset + c'x`.
    pub objective_linear: f64,
    /// `x'Qx/2`.
    pub objective_quadratic: f64,
    pub duals_present: bool,
    /// `optimal`, `not_optimal` or `unknown` (no duals).
    pub kkt_verdict: &'static str,
    pub tolerance: f64,
    /// Max |c + Qx - A'y - d|, scaled by 1 + the largest term.
    pub stationarity: Option<f64>,
    /// Largest wrong-sign dual (as a minimisation).
    pub dual_infeasibility: Option<f64>,
    /// Largest |dual * slack to its bound|, scaled.
    pub complementarity: Option<f64>,
    pub dual_objective: Option<f64>,
    /// `|primal - dual|` and its value relative to `1 + |primal| + |dual|`.
    pub duality_gap: Option<f64>,
    pub duality_gap_rel: Option<f64>,
}

/// Parses the `# Dual solution values` section of a raw solution file, if it has one.
pub fn parse_duals(text: &str, lp: &HighsLp) -> Option<Duals> {
    let mut lines = text.lines();
    lines
        .by_ref()
        .find(|l| l.trim() == "# Dual solution values")?;
    let status = lines.by_ref().map(str::trim).find(|t| !t.is_empty())?;
    if status == "None" {
        return None;
    }
    let mut section = |label: &str, n: usize, names: &[String], prefix: &str| -> Option<Vec<f64>> {
        let head = lines.by_ref().find(|l| l.trim_start().starts_with("# "))?;
        if !head.contains(label) {
            return None;
        }
        let mut by_name = std::collections::HashMap::new();
        for l in lines.by_ref().take(n) {
            let mut it = l.split_whitespace();
            by_name.insert(it.next()?.to_string(), it.next()?.parse::<f64>().ok()?);
        }
        (0..n)
            .map(|i| {
                let name = match names.get(i) {
                    Some(s) if !s.is_empty() => s.clone(),
                    _ => format!("{prefix}{i}"),
                };
                by_name.get(&name).copied()
            })
            .collect()
    };
    let col = section("Columns", lp.num_col_ as usize, &lp.col_names_, "c")?;
    let row = section("Rows", lp.num_row_ as usize, &lp.row_names_, "r")?;
    Some(Duals { col, row })
}

/// `Q x` from the symmetric triplets.
fn q_times(q: &[Triplet], x: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; x.len()];
    for &(i, j, v) in q {
        out[i] += v * x[j];
    }
    out
}

pub fn check_qp(lp: &HighsLp, q: &[Triplet], x: &[f64], duals: Option<&Duals>) -> QpCheck {
    let n = x.len();
    let s = if lp.sense_ == ObjSense::kMaximize {
        -1.0
    } else {
        1.0
    };
    let qx = q_times(q, x);
    let quad = 0.5 * x.iter().zip(&qx).map(|(a, b)| a * b).sum::<f64>();
    let lin = lp.offset_ + lp.col_cost_.iter().zip(x).map(|(c, v)| c * v).sum::<f64>();
    let mut rep = QpCheck {
        objective_linear: lin,
        objective_quadratic: quad,
        duals_present: duals.is_some(),
        kkt_verdict: "unknown",
        tolerance: KKT_TOL,
        stationarity: None,
        dual_infeasibility: None,
        complementarity: None,
        dual_objective: None,
        duality_gap: None,
        duality_gap_rel: None,
    };
    let Some(du) = duals else { return rep };
    if du.col.len() != n || du.row.len() != lp.num_row_ as usize {
        return rep;
    }
    // The minimisation form: c' = s c, Q' = s Q, y' = s y, d' = s d.
    let a = &lp.a_matrix_;
    let mut aty = vec![0.0; n];
    let mut ax = vec![0.0; lp.num_row_ as usize];
    for j in 0..n {
        for k in a.start_[j] as usize..a.start_[j + 1] as usize {
            let i = a.index_[k] as usize;
            aty[j] += a.value_[k] * s * du.row[i];
            ax[i] += a.value_[k] * x[j];
        }
    }
    let (mut stat, mut stat_scale, mut stat_x) = (0.0f64, 0.0f64, 0.0f64);
    for j in 0..n {
        let terms = [s * lp.col_cost_[j], s * qx[j], aty[j], s * du.col[j]];
        let r = terms[0] + terms[1] - terms[2] - terms[3];
        stat = stat.max(r.abs());
        stat_x += r.abs() * x[j].abs();
        stat_scale = stat_scale.max(terms.iter().fold(0.0f64, |m, t| m.max(t.abs())));
    }
    let (mut infeas, mut comp) = (0.0f64, 0.0f64);
    let mut bound_terms = 0.0;
    // Columns and rows alike: `dual` (minimisation sign) against the interval [lo, hi] at `v`.
    let mut one = |dual: f64, v: f64, lo: f64, hi: f64| {
        let (slack_lo, slack_hi) = (v - lo, hi - v);
        // d >= 0 needs a finite lower bound, d <= 0 a finite upper one.
        if dual > 0.0 {
            if lo.is_finite() {
                comp = comp.max(dual * slack_lo.abs());
                bound_terms += dual * lo;
            } else {
                infeas = infeas.max(dual);
            }
        } else if dual < 0.0 {
            if hi.is_finite() {
                comp = comp.max(-dual * slack_hi.abs());
                bound_terms += dual * hi;
            } else {
                infeas = infeas.max(-dual);
            }
        }
    };
    for (j, &xj) in x.iter().enumerate() {
        one(s * du.col[j], xj, lp.col_lower_[j], lp.col_upper_[j]);
    }
    for (i, &axi) in ax.iter().enumerate() {
        one(s * du.row[i], axi, lp.row_lower_[i], lp.row_upper_[i]);
    }
    let primal = s * (lin + quad);
    // Wolfe dual in the minimisation form: off' - x'Q'x/2 + bound terms, off' = s*offset.
    let dual_val = s * lp.offset_ - s * quad + bound_terms;
    let gap = (primal - dual_val).abs();
    let gap_rel = gap / (1.0 + primal.abs() + dual_val.abs());
    let stat_rel = stat / (1.0 + stat_scale);
    let comp_rel = comp / (1.0 + primal.abs());
    rep.stationarity = Some(stat_rel);
    rep.dual_infeasibility = Some(infeas);
    rep.complementarity = Some(comp_rel);
    rep.dual_objective = Some(s * dual_val);
    rep.duality_gap = Some(gap);
    rep.duality_gap_rel = Some(gap_rel);
    // The gap is the summed complementarity plus r'x (r the stationarity residual), so a residual
    // inside tolerance still moves the gap by up to sum |r_j x_j|; that part is not a defect of
    // the answer beyond what `stationarity` already says.
    let gap_ok = gap <= KKT_TOL * (1.0 + primal.abs() + dual_val.abs()) + stat_x;
    rep.kkt_verdict = if stat_rel <= KKT_TOL && infeas <= KKT_TOL && comp_rel <= KKT_TOL && gap_ok {
        "optimal"
    } else {
        "not_optimal"
    };
    rep
}
