//! Continuous SOCP and convex QCQP through the
//! Clarabel interior-point solver (Goulart and Chen 2024, arXiv:2405.12762; Apache-2.0).
//!
//! `rusolve solve` and `rusolve mip` land here when the LP/MIP reader refused a model for its
//! quadratic rows or cones ([`rusolve_io::conic::conic_rejected`]); a model with neither never
//! reaches this module, so its code path is exactly the pre-quadratic-rows one.
//!
//! Conic form (Clarabel): `min ½xᵀPx + qᵀx  s.t.  Ax + s = b, s in K`, dual `z in K*`,
//! stationarity `Px + q + Aᵀz = 0`. The model's columns are Clarabel's `x`. Rows, in the order
//! Clarabel wants them (zero cone, nonnegative cone, then second-order cones):
//!
//! * a constraint row `lo <= aᵀx <= hi`: equality when `lo == hi`, else one nonnegative row per
//!   finite side (`aᵀx + s = hi`; `-aᵀx + s = -lo`);
//! * a variable bound: the same, with `a = e_j`;
//! * a user cone: one row `-e_v` per member (`s = x_v`), the rotated cone through the identity
//!   `2 x0 x1 >= ‖w‖²  <=>  ‖(x0 - x1, √2 w)‖ <= x0 + x1`;
//! * a convex quadratic row `aᵀx + ½xᵀQx <= b` with `Q = LLᵀ` ([`crate::quad_convexity`]):
//!   `½‖Lᵀx‖² <= b - aᵀx`, which with `u = b - aᵀx` and `2·u·1 >= ‖w‖²`, `w = Lᵀx` is the
//!   cone `s = (u + 1, u - 1, √2 Lᵀx)` (Lobo, Vandenberghe, Boyd, Lebret 1998, section 2.3).
//!   A `>=` side is the same with `-a`, `-Q`, `-lo`.
//!
//! The user's rows and cones are kept in the model; the checker ([`crate::check_conic`])
//! evaluates those, never this reformulation. Certificates are verified on the conic form.
//!
//! Anything the route cannot take is a typed [`ConicError`] naming the row, column or cone:
//! a nonconvex row or objective, an integer column, an equality with a quadratic term, a
//! reserved cone kind.

use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSettings, DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT, SolverStatus,
    SupportedConeT, ZeroConeT,
};
use rusolve_io::Model;
use rusolve_io::conic::ConeKind;
use serde_json::{Value, json};
use std::fmt;
use std::process::ExitCode;

use crate::check_conic::{DEFAULT_TOL, check_conic_solution, row_interval};
use crate::envelope::{NativeSolved, native_envelope};
use crate::quad_convexity::{Refusal, factor_psd};

/// The `clarabel` crate version in `Cargo.toml` (reported in the envelope).
pub const CLARABEL_VERSION: &str = "0.11.1";

/// Why a model cannot take the conic route.
#[derive(Debug, Clone, PartialEq)]
pub enum ConicError {
    /// A quadratic row whose `Q` has negative curvature on the side it is used on.
    NonconvexRow {
        row: String,
        side: &'static str,
        var: String,
        detail: String,
    },
    NonconvexObjective {
        var: String,
        sense: &'static str,
    },
    IntegerVariable {
        var: String,
    },
    /// `aᵀx + ½xᵀQx = rhs` with `Q != 0` is never convex.
    EqualityQuadRow {
        row: String,
    },
    TooLarge {
        row: String,
        support: usize,
    },
    UnsupportedCone {
        cone: String,
        kind: &'static str,
    },
    Solver(String),
}

impl fmt::Display for ConicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonconvexRow {
                row,
                side,
                var,
                detail,
            } => write!(
                f,
                "quadratic row {row} is not convex on its {side} side: {detail} at column {var}"
            ),
            Self::NonconvexObjective { var, sense } => write!(
                f,
                "the objective is not convex (a {sense} needs a {} Hessian): negative curvature at column {var}",
                if *sense == "minimization" {
                    "PSD"
                } else {
                    "NSD"
                }
            ),
            Self::IntegerVariable { var } => write!(
                f,
                "column {var} is integer in a model with quadratic constraints or cones: only continuous SOCP / convex QCQP is supported"
            ),
            Self::EqualityQuadRow { row } => write!(
                f,
                "quadratic row {row} is an equality with a quadratic term, which is not convex"
            ),
            Self::TooLarge { row, support } => write!(
                f,
                "quadratic row {row} has {support} columns with a quadratic term, more than the dense factor accepts"
            ),
            Self::UnsupportedCone { cone, kind } => {
                write!(
                    f,
                    "cone {cone} is a {kind} cone, which the Clarabel route does not take yet"
                )
            }
            Self::Solver(m) => write!(f, "clarabel: {m}"),
        }
    }
}

fn col_name(m: &Model, j: usize) -> String {
    match m.col_names.get(j) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("c{j}"),
    }
}

fn row_name(m: &Model, i: usize) -> String {
    match m.row_names.get(i) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("r{i}"),
    }
}

/// One conic row before the cone blocks are concatenated.
struct Row {
    entries: Vec<(usize, f64)>,
    b: f64,
}

/// `y_i` (the row's rhs sensitivity) gets `sign * z[pending row]` from each listed row.
struct DualLink {
    row: usize,
    sign: f64,
    group: Group,
    index: usize,
}

#[derive(Clone, Copy)]
enum Group {
    Zero,
    Nonneg,
    Soc(usize),
}

/// The conic form plus what is needed to read the answer back.
pub struct ConicForm {
    p: CscMatrix<f64>,
    q: Vec<f64>,
    a: CscMatrix<f64>,
    b: Vec<f64>,
    cones: Vec<SupportedConeT<f64>>,
    /// `true` when the objective was negated (maximize).
    negated: bool,
    duals: Vec<(usize, f64, usize)>,
    /// (cone name, first conic row, size) of each user cone.
    user_cones: Vec<(String, usize, usize)>,
}

fn csc(m: usize, n: usize, mut t: Vec<(usize, usize, f64)>) -> CscMatrix<f64> {
    t.sort_by_key(|&(i, j, _)| (j, i));
    let mut colptr = vec![0usize; n + 1];
    let (mut rowval, mut nzval): (Vec<usize>, Vec<f64>) = (Vec::new(), Vec::new());
    let mut last: Option<(usize, usize)> = None;
    for (i, j, v) in t {
        if last == Some((i, j)) {
            *nzval.last_mut().unwrap() += v;
            continue;
        }
        last = Some((i, j));
        rowval.push(i);
        nzval.push(v);
        colptr[j + 1] += 1;
    }
    for j in 0..n {
        colptr[j + 1] += colptr[j];
    }
    CscMatrix::new(m, n, colptr, rowval, nzval)
}

fn refusal_error(m: &Model, row: &str, side: &'static str, r: Refusal) -> ConicError {
    match r {
        Refusal::Indefinite { var, other, pivot } => ConicError::NonconvexRow {
            row: row.to_string(),
            side,
            var: col_name(m, var),
            detail: match other {
                Some(o) => format!("indefinite with {} (pivot {pivot:.3e})", col_name(m, o)),
                None => format!("negative pivot {pivot:.3e}"),
            },
        },
        Refusal::TooLarge { support } => ConicError::TooLarge {
            row: row.to_string(),
            support,
        },
    }
}

/// Builds the conic form of `m`, or says why it cannot.
pub fn build_form(m: &Model) -> Result<ConicForm, ConicError> {
    let n = m.num_variables;
    if let Some(j) = m.bounds.iter().position(|b| b.integer) {
        return Err(ConicError::IntegerVariable {
            var: col_name(m, j),
        });
    }
    for c in &m.conic.cones {
        if !matches!(c.kind, ConeKind::Soc | ConeKind::RotatedSoc) {
            return Err(ConicError::UnsupportedCone {
                cone: c.name.clone(),
                kind: c.kind.label(),
            });
        }
    }
    let negated = m.obj_sense == rusolve_io::ObjSense::Maximize;
    let sgn = if negated { -1.0 } else { 1.0 };

    // objective: the Hessian is the full square matrix; Clarabel wants the upper triangle
    let lower: Vec<_> = m
        .hessian
        .iter()
        .filter(|&&(i, j, _)| i >= j)
        .copied()
        .collect();
    if let Err(r) = factor_psd(&lower, sgn) {
        if let Refusal::Indefinite { var, .. } = r {
            return Err(ConicError::NonconvexObjective {
                var: col_name(m, var),
                sense: if negated {
                    "maximization"
                } else {
                    "minimization"
                },
            });
        }
        return Err(ConicError::Solver(
            "objective Hessian is too large to test".into(),
        ));
    }
    let p_trip: Vec<_> = m
        .hessian
        .iter()
        .filter(|&&(i, j, _)| i <= j)
        .map(|&(i, j, v)| (i, j, sgn * v))
        .collect();
    let q: Vec<f64> = m.objective.iter().map(|c| sgn * c).collect();

    let (mut zero, mut nonneg): (Vec<Row>, Vec<Row>) = (Vec::new(), Vec::new());
    let mut socs: Vec<Vec<Row>> = Vec::new();
    let mut links: Vec<DualLink> = Vec::new();
    let mut user_soc: Vec<(String, usize)> = Vec::new();

    // constraint rows: linear part from the CSC, by row
    let (cp, ri, va) = (m.matrix.col_ptr(), m.matrix.row_idx(), m.matrix.values());
    let mut lin: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m.num_constraints];
    for j in 0..n {
        for k in cp[j]..cp[j + 1] {
            lin[ri[k] as usize].push((j, va[k]));
        }
    }
    let mut quad_of = vec![usize::MAX; m.num_constraints];
    for (k, qr) in m.conic.quad_rows.iter().enumerate() {
        quad_of[qr.row] = k;
    }
    for i in 0..m.num_constraints {
        let (lo, hi) = row_interval(m, i);
        let name = row_name(m, i);
        let factor = |sign: f64, side: &'static str| match quad_of[i] {
            usize::MAX => Ok(None),
            k => match factor_psd(&m.conic.quad_rows[k].q, sign) {
                Ok(f) if f.l.is_empty() => Ok(None),
                Ok(f) => Ok(Some(f)),
                Err(r) => Err(refusal_error(m, &name, side, r)),
            },
        };
        if lo == hi {
            if factor(1.0, "upper")?.is_some() || factor(-1.0, "lower")?.is_some() {
                return Err(ConicError::EqualityQuadRow { row: name });
            }
            links.push(DualLink {
                row: i,
                sign: -1.0,
                group: Group::Zero,
                index: zero.len(),
            });
            zero.push(Row {
                entries: lin[i].clone(),
                b: hi,
            });
            continue;
        }
        for (side, sign, bound) in [("upper", 1.0, hi), ("lower", -1.0, lo)] {
            if !bound.is_finite() {
                continue;
            }
            let lin_part: Vec<(usize, f64)> = lin[i].iter().map(|&(j, v)| (j, sign * v)).collect();
            let b = sign * bound;
            // dual: d obj / d rhs_i = -z (upper), +z (lower)
            let link_sign = -sign;
            match factor(sign, side)? {
                None => {
                    links.push(DualLink {
                        row: i,
                        sign: link_sign,
                        group: Group::Nonneg,
                        index: nonneg.len(),
                    });
                    nonneg.push(Row {
                        entries: lin_part,
                        b,
                    });
                }
                Some(f) => {
                    let g = socs.len();
                    let mut rows = vec![
                        Row {
                            entries: lin_part.clone(),
                            b: b + 1.0,
                        },
                        Row {
                            entries: lin_part,
                            b: b - 1.0,
                        },
                    ];
                    for c in &f.l {
                        let e: Vec<(usize, f64)> = f
                            .vars
                            .iter()
                            .zip(c)
                            .filter(|&(_, &v)| v != 0.0)
                            .map(|(&j, &v)| (j, -std::f64::consts::SQRT_2 * v))
                            .collect();
                        rows.push(Row { entries: e, b: 0.0 });
                    }
                    for index in 0..2 {
                        links.push(DualLink {
                            row: i,
                            sign: link_sign,
                            group: Group::Soc(g),
                            index,
                        });
                    }
                    socs.push(rows);
                }
            }
        }
    }

    // variable bounds
    for (j, bd) in m.bounds.iter().enumerate() {
        if bd.lower == bd.upper {
            zero.push(Row {
                entries: vec![(j, 1.0)],
                b: bd.lower,
            });
            continue;
        }
        if bd.upper.is_finite() {
            nonneg.push(Row {
                entries: vec![(j, 1.0)],
                b: bd.upper,
            });
        }
        if bd.lower.is_finite() {
            nonneg.push(Row {
                entries: vec![(j, -1.0)],
                b: -bd.lower,
            });
        }
    }

    // user cones
    for c in &m.conic.cones {
        let v = &c.vars;
        let rows: Vec<Row> = match c.kind {
            ConeKind::Soc => v
                .iter()
                .map(|&j| Row {
                    entries: vec![(j, -1.0)],
                    b: 0.0,
                })
                .collect(),
            _ => {
                let mut r = vec![
                    Row {
                        entries: vec![(v[0], -1.0), (v[1], -1.0)],
                        b: 0.0,
                    },
                    Row {
                        entries: vec![(v[0], -1.0), (v[1], 1.0)],
                        b: 0.0,
                    },
                ];
                r.extend(v[2..].iter().map(|&j| Row {
                    entries: vec![(j, -std::f64::consts::SQRT_2)],
                    b: 0.0,
                }));
                r
            }
        };
        user_soc.push((c.name.clone(), socs.len()));
        socs.push(rows);
    }

    // concatenate: zero, nonneg, soc blocks
    let mut a_trip = Vec::new();
    let mut b = Vec::new();
    let mut cones: Vec<SupportedConeT<f64>> = Vec::new();
    fn push(r: &Row, a_trip: &mut Vec<(usize, usize, f64)>, b: &mut Vec<f64>) {
        let at = b.len();
        for &(j, v) in &r.entries {
            a_trip.push((at, j, v));
        }
        b.push(r.b);
    }
    for r in &zero {
        push(r, &mut a_trip, &mut b);
    }
    if !zero.is_empty() {
        cones.push(ZeroConeT(zero.len()));
    }
    let nonneg_at = zero.len();
    for r in &nonneg {
        push(r, &mut a_trip, &mut b);
    }
    if !nonneg.is_empty() {
        cones.push(NonnegativeConeT(nonneg.len()));
    }
    let mut soc_at = Vec::new();
    for rows in &socs {
        soc_at.push(b.len());
        for r in rows {
            push(r, &mut a_trip, &mut b);
        }
        cones.push(SecondOrderConeT(rows.len()));
    }
    let duals = links
        .iter()
        .map(|l| {
            let at = match l.group {
                Group::Zero => l.index,
                Group::Nonneg => nonneg_at + l.index,
                Group::Soc(g) => soc_at[g] + l.index,
            };
            (l.row, l.sign, at)
        })
        .collect();
    let user_cones = user_soc
        .into_iter()
        .map(|(name, g)| (name, soc_at[g], socs[g].len()))
        .collect();
    let rows = b.len();
    Ok(ConicForm {
        p: csc(n, n, p_trip),
        q,
        a: csc(rows, n, a_trip),
        b,
        cones,
        negated,
        duals,
        user_cones,
    })
}

// ---- certificate verification on the conic form -------------------------------------------

fn mul_a(a: &CscMatrix<f64>, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; a.m];
    for (j, &xj) in x.iter().enumerate().take(a.n) {
        for k in a.colptr[j]..a.colptr[j + 1] {
            y[a.rowval[k]] += a.nzval[k] * xj;
        }
    }
    y
}

fn mul_at(a: &CscMatrix<f64>, z: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; a.n];
    for (j, yj) in y.iter_mut().enumerate() {
        for k in a.colptr[j]..a.colptr[j + 1] {
            *yj += a.nzval[k] * z[a.rowval[k]];
        }
    }
    y
}

/// `P x` for the upper-triangular `P`.
fn mul_p_sym(p: &CscMatrix<f64>, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; p.n];
    for j in 0..p.n {
        for k in p.colptr[j]..p.colptr[j + 1] {
            let i = p.rowval[k];
            y[i] += p.nzval[k] * x[j];
            if i != j {
                y[j] += p.nzval[k] * x[i];
            }
        }
    }
    y
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |m, a| m.max(a.abs()))
}

/// Largest violation of `v in K` (or of `v in K*`, which has the same cones: they are
/// self-dual except the zero cone, whose dual is free).
fn cone_violation(cones: &[SupportedConeT<f64>], v: &[f64], dual: bool) -> f64 {
    let mut at = 0;
    let mut worst = 0.0f64;
    for c in cones {
        match c {
            SupportedConeT::ZeroConeT(d) => {
                if !dual {
                    worst = worst.max(inf_norm(&v[at..at + d]));
                }
                at += d;
            }
            SupportedConeT::NonnegativeConeT(d) => {
                for &s in &v[at..at + d] {
                    worst = worst.max(-s);
                }
                at += d;
            }
            SupportedConeT::SecondOrderConeT(d) => {
                let t = &v[at..at + d];
                let norm = t[1..].iter().map(|a| a * a).sum::<f64>().sqrt();
                worst = worst.max(norm - t[0]);
                at += d;
            }
            _ => worst = f64::INFINITY,
        }
    }
    worst
}

/// A certificate of infeasibility, verified: Farkas `z in K*, Aᵀz = 0, bᵀz < 0` for primal
/// infeasibility; `Px = 0, qᵀx < 0, -Ax in K` for dual infeasibility (Clarabel docs, "Algorithm";
/// Boyd and Vandenberghe 2004, section 5.8 for the alternatives).
fn verify_certificate(f: &ConicForm, primal: bool, x: &[f64], z: &[f64]) -> Value {
    let tol = 1e-5;
    if primal {
        let atz = inf_norm(&mul_at(&f.a, z));
        let btz: f64 = f.b.iter().zip(z).map(|(b, z)| b * z).sum();
        let scale = 1.0 + inf_norm(z);
        let cone = cone_violation(&f.cones, z, true);
        json!({
            "kind": "primal_infeasible_farkas",
            "verified": atz <= tol * scale && cone <= tol * scale && btz < -1e-8 * scale,
            "residuals": {"at_z_inf": atz, "b_t_z": btz, "dual_cone_violation": cone},
        })
    } else {
        let px = inf_norm(&mul_p_sym(&f.p, x));
        let qtx: f64 = f.q.iter().zip(x).map(|(q, x)| q * x).sum();
        let s: Vec<f64> = mul_a(&f.a, x).iter().map(|v| -v).collect();
        let scale = 1.0 + inf_norm(x);
        let cone = cone_violation(&f.cones, &s, false);
        json!({
            "kind": "dual_infeasible_ray",
            "verified": px <= tol * scale && cone <= tol * scale && qtx < -1e-8 * scale,
            "residuals": {"p_x_inf": px, "q_t_x": qtx, "cone_violation": cone},
        })
    }
}

fn status_name(s: SolverStatus) -> &'static str {
    match s {
        SolverStatus::Unsolved => "Unsolved",
        SolverStatus::Solved => "Solved",
        SolverStatus::PrimalInfeasible => "PrimalInfeasible",
        SolverStatus::DualInfeasible => "DualInfeasible",
        SolverStatus::AlmostSolved => "AlmostSolved",
        SolverStatus::AlmostPrimalInfeasible => "AlmostPrimalInfeasible",
        SolverStatus::AlmostDualInfeasible => "AlmostDualInfeasible",
        SolverStatus::MaxIterations => "MaxIterations",
        SolverStatus::MaxTime => "MaxTime",
        SolverStatus::NumericalError => "NumericalError",
        SolverStatus::InsufficientProgress => "InsufficientProgress",
        SolverStatus::CallbackTerminated => "CallbackTerminated",
    }
}

/// What a solve produced, in the model's own sense.
pub struct ConicSolved {
    /// One of the result schema's `status` values.
    pub status: &'static str,
    pub clarabel_status: &'static str,
    pub x: Vec<f64>,
    /// `NaN` without a solution.
    pub objective: f64,
    pub iterations: u32,
    pub time_s: f64,
    /// The `conic` block of the envelope.
    pub detail: Value,
}

/// Solves `m` (continuous SOCP / convex QCQP) with Clarabel.
pub fn solve_model(m: &Model, time_limit: Option<f64>) -> Result<ConicSolved, ConicError> {
    let f = build_form(m)?;
    let settings = DefaultSettings {
        verbose: false,
        time_limit: time_limit.unwrap_or(f64::INFINITY),
        ..DefaultSettings::default()
    };
    let mut solver = DefaultSolver::new(&f.p, &f.q, &f.a, &f.b, &f.cones, settings)
        .map_err(|e| ConicError::Solver(format!("{e:?}")))?;
    solver.solve();
    let sol = &solver.solution;
    let cs = status_name(sol.status);
    let (status, has_x) = match sol.status {
        SolverStatus::Solved | SolverStatus::AlmostSolved => ("Optimal", true),
        SolverStatus::PrimalInfeasible | SolverStatus::AlmostPrimalInfeasible => {
            ("Infeasible", false)
        }
        SolverStatus::DualInfeasible | SolverStatus::AlmostDualInfeasible => ("Unbounded", false),
        SolverStatus::MaxIterations => ("IterationLimit", false),
        SolverStatus::MaxTime => ("TimeLimit", false),
        SolverStatus::CallbackTerminated => ("Interrupted", false),
        SolverStatus::NumericalError => ("Error", false),
        SolverStatus::InsufficientProgress | SolverStatus::Unsolved => ("Unknown", false),
    };
    let mut detail = json!({
        "engine": "clarabel",
        "clarabel_version": CLARABEL_VERSION,
        "clarabel_status": cs,
        "reduced_accuracy": matches!(
            sol.status,
            SolverStatus::AlmostSolved
                | SolverStatus::AlmostPrimalInfeasible
                | SolverStatus::AlmostDualInfeasible
        ),
        "iterations": sol.iterations,
        "conic_rows": f.b.len(),
        "r_prim": sol.r_prim,
        "r_dual": sol.r_dual,
    });
    let mut status = status;
    let mut x = Vec::new();
    let mut objective = f64::NAN;
    if has_x {
        x = sol.x.clone();
        let chk = check_conic_solution(m, &x, DEFAULT_TOL, 0, "", "");
        objective = chk.objective.value;
        detail["objective_primal"] =
            json!(if f.negated { -sol.obj_val } else { sol.obj_val } + m.obj_offset);
        detail["objective_dual"] = json!(
            if f.negated {
                -sol.obj_val_dual
            } else {
                sol.obj_val_dual
            } + m.obj_offset
        );
        let max_violation = [
            chk.summary.column_bounds.max,
            chk.summary.rows.max,
            chk.summary.quad_rows.max,
            chk.summary.cones.max,
        ]
        .into_iter()
        .fold(0.0f64, f64::max);
        detail["checker"] = json!({
            "verdict": chk.verdict,
            "tolerance": chk.tolerance,
            "violations": chk.summary.total_violations,
            "max_violation": max_violation,
        });
        // dual of each constraint row: d objective / d rhs, in the model's sense
        let mut y = vec![0.0; m.num_constraints];
        for &(row, sign, at) in &f.duals {
            y[row] += sign * sol.z[at];
        }
        let flip = if f.negated { -1.0 } else { 1.0 };
        detail["row_duals"] = Value::Array(
            y.iter()
                .enumerate()
                .map(|(i, v)| json!({"name": row_name(m, i), "value": flip * v}))
                .collect(),
        );
        detail["cone_duals"] = Value::Array(
            f.user_cones
                .iter()
                .map(|(name, at, d)| json!({"name": name, "z": sol.z[*at..*at + d]}))
                .collect(),
        );
    } else if status == "Infeasible" || status == "Unbounded" {
        let cert = verify_certificate(&f, status == "Infeasible", &sol.x, &sol.z);
        if cert["verified"] != json!(true) {
            // a certificate that does not check is not an answer
            status = "Unknown";
        }
        detail["certificate"] = cert;
    }
    Ok(ConicSolved {
        status,
        clarabel_status: cs,
        x,
        objective,
        iterations: sol.iterations,
        time_s: sol.solve_time,
        detail,
    })
}

struct Flags {
    json: bool,
    include_solution: bool,
    time_limit: Option<f64>,
}

fn parse_flags(cmd: &str, flags: &[String]) -> Result<Flags, String> {
    let mut f = Flags {
        json: false,
        include_solution: false,
        time_limit: None,
    };
    let mut it = flags.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--format" => match it.next().map(String::as_str) {
                Some("json") => f.json = true,
                Some("text") => {}
                _ => return Err("--format needs json or text for a conic model".into()),
            },
            "--solution" | "--include-solution" => f.include_solution = true,
            "--time-limit" => match it.next().and_then(|v| v.parse::<f64>().ok()) {
                Some(t) if t >= 0.0 && t.is_finite() => f.time_limit = Some(t),
                _ => return Err("--time-limit: want a non-negative number of seconds".into()),
            },
            // accepted and ignored: the conic route is deterministic and single-threaded
            // `--no-warnings` is a no-op here: the conic route has no lint pass
            "--no-warnings" => {}
            "--threads" | "--seed" => {
                if it.next().is_none() {
                    return Err(format!("{a} needs a value"));
                }
            }
            other => {
                return Err(format!(
                    "{cmd}: argument '{other}' is not supported for a model with quadratic constraints or cones (supported: --format json|text, --time-limit, --include-solution)"
                ));
            }
        }
    }
    if f.include_solution && !f.json {
        return Err("--solution needs --format json".into());
    }
    Ok(f)
}

/// `rusolve solve|mip <model> [flags]` for a model the LP/MIP reader refused for its conic
/// content. `args` is what followed the command name.
pub fn run(cmd: &str, args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let Some(path) = args.first() else {
        return crate::diag::fail_msg(
            wj,
            "RS0006",
            2,
            format!("rusolve {cmd}: missing model path"),
        );
    };
    let flags = match parse_flags(cmd, &args[1..]) {
        Ok(f) => f,
        Err(e) => {
            return crate::diag::fail_msg(wj, "RS0003", 2, format!("rusolve {cmd}: {e}"));
        }
    };
    let model = match rusolve_io::read_model_file_conic(path) {
        Ok(m) => m,
        Err(e) => return crate::diag::fail_read(flags.json, path, &e.to_string()),
    };
    let solved = match solve_model(&model, flags.time_limit) {
        Ok(s) => s,
        Err(e) => {
            return crate::diag::fail_msg(flags.json, "RS0012", 3, format!("rusolve {cmd}: {e}"));
        }
    };
    if flags.json {
        let n = NativeSolved {
            path,
            model: &model,
            status: solved.status,
            objective: solved.objective,
            dual_bound: solved.objective,
            nodes: 0,
            lp_iterations: i64::from(solved.iterations),
            time_s: solved.time_s,
            seed: 0,
            threads: 1,
            time_limit: flags.time_limit,
            node_limit: None,
            x: &solved.x,
            include_solution: flags.include_solution,
        };
        let mut v = serde_json::to_value(native_envelope(&n)).expect("envelope serialises");
        v["solver"]["engine"] = json!("clarabel");
        v["conic"] = solved.detail;
        println!("{v}");
    } else {
        println!("status: {}", solved.status);
        if !solved.x.is_empty() {
            println!("objective: {}", solved.objective);
        }
        println!(
            "clarabel: {} in {} iterations ({:.3}s)",
            solved.clarabel_status, solved.iterations, solved.time_s
        );
        if let Some(c) = solved.detail.get("checker") {
            println!("checker: {} ({} violations)", c["verdict"], c["violations"]);
        }
        if let Some(c) = solved.detail.get("certificate") {
            println!("certificate: {} verified={}", c["kind"], c["verified"]);
        }
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_io::read_mps_conic;

    fn solve(src: &str) -> ConicSolved {
        solve_model(&read_mps_conic(src).unwrap(), None).unwrap()
    }

    fn err(src: &str) -> ConicError {
        solve_model(&read_mps_conic(src).unwrap(), None)
            .err()
            .unwrap()
    }

    /// min t  s.t. x + y = 2, t >= ‖(x, y)‖  ->  x = y = 1, t = sqrt 2.
    const SOC: &str = "\
NAME soc
ROWS
 N  obj
 E  e1
COLUMNS
    t  obj  1
    x  e1  1
    y  e1  1
RHS
    rhs  e1  2
BOUNDS
 FR bnd t
 FR bnd x
 FR bnd y
CSECTION k1 0.0 QUAD
    t
    x
    y
ENDATA
";

    #[test]
    fn a_second_order_cone_program_is_solved_and_checked() {
        let s = solve(SOC);
        assert_eq!(s.status, "Optimal");
        assert!((s.objective - 2f64.sqrt()).abs() < 1e-6, "{}", s.objective);
        assert_eq!(s.detail["checker"]["verdict"], "feasible");
        // the cone dual is reported by name
        assert_eq!(s.detail["cone_duals"][0]["name"], "k1");
        // d obj / d rhs of x + y = 2: obj = sqrt(2)/2 * rhs
        let y = s.detail["row_duals"][0]["value"].as_f64().unwrap();
        assert!((y - 2f64.sqrt() / 2.0).abs() < 1e-5, "{y}");
    }

    /// min x + y s.t. x² + y² <= 1 (as a quadratic row, Q = 2 I), x, y free:
    /// optimum -sqrt 2 at x = y = -1/sqrt 2.
    const QCP: &str = "\
NAME qcp
ROWS
 N  obj
 L  q1
COLUMNS
    x  obj  1
    y  obj  1
RHS
    rhs  q1  1
BOUNDS
 FR bnd x
 FR bnd y
QCMATRIX q1
    x  x  1
    y  y  1
ENDATA
";

    #[test]
    fn a_convex_quadratic_row_goes_through_the_cholesky_cone() {
        let s = solve(QCP);
        assert_eq!(s.status, "Optimal");
        assert!((s.objective + 2f64.sqrt()).abs() < 1e-6, "{}", s.objective);
        assert!((s.x[0] + s.x[1] + 2f64.sqrt()).abs() < 1e-6);
        assert_eq!(s.detail["checker"]["verdict"], "feasible");
        // d obj / d rhs: tightening the disc raises the optimum by sqrt 2 / 2 per unit of
        // radius², so the dual is -sqrt(2)/2
        let y = s.detail["row_duals"][0]["value"].as_f64().unwrap();
        assert!((y + 2f64.sqrt() / 2.0).abs() < 1e-5, "{y}");
    }

    #[test]
    fn infeasible_and_unbounded_models_come_with_a_verified_certificate() {
        // x² <= 1 together with x >= 2
        let inf = QCP.replace("BOUNDS\n FR bnd x", "BOUNDS\n LO bnd x 2\n FR bnd z");
        let s = solve(&inf.replace(" FR bnd z\n", ""));
        assert_eq!(s.status, "Infeasible", "{:?}", s.detail);
        assert_eq!(s.detail["certificate"]["verified"], true, "{:?}", s.detail);
        // min -x with x free in a cone-free-direction: t >= ‖x‖ only bounds t
        let unb = SOC
            .replace("    t  obj  1", "    t  obj  0\n    x  obj  -1")
            .replace("    x  e1  1\n", "    x  e1  0\n");
        let s = solve(&unb);
        assert_eq!(s.status, "Unbounded", "{:?}", s.detail);
        assert_eq!(s.detail["certificate"]["verified"], true, "{:?}", s.detail);
    }

    #[test]
    fn a_nonconvex_row_names_the_row_and_column() {
        let e = err(&QCP.replace("    y  y  1", "    y  y  -1"));
        match e {
            ConicError::NonconvexRow { row, var, .. } => {
                assert_eq!(row, "q1");
                assert_eq!(var, "y");
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn a_ge_row_needs_negative_semidefinite_q() {
        // x² + y² >= 1 is nonconvex
        let e = err(&QCP.replace(" L  q1", " G  q1"));
        assert!(
            matches!(e, ConicError::NonconvexRow { side: "lower", .. }),
            "{e:?}"
        );
        // -x² - y² >= -1 is the disc again
        let s = solve(
            &QCP.replace(" L  q1", " G  q1")
                .replace("x  x  1", "x  x  -1")
                .replace("y  y  1", "y  y  -1")
                .replace("rhs  q1  1", "rhs  q1  -1"),
        );
        assert_eq!(s.status, "Optimal");
    }

    #[test]
    fn an_integer_column_is_a_typed_error() {
        let src = QCP.replace(" FR bnd x\n", " BV bnd x\n");
        assert_eq!(err(&src), ConicError::IntegerVariable { var: "x".into() });
    }

    #[test]
    fn a_rotated_cone_is_solved() {
        // min x + y s.t. 2 x y >= 1 as RQUAD (x, y, w) with w fixed to 1: x = y = 1/sqrt 2... 2xy >= 1
        let src = "\
NAME rot
ROWS
 N  obj
COLUMNS
    x  obj  1
    y  obj  1
    w  obj  0
BOUNDS
 FX bnd w 1
CSECTION k1 0.0 RQUAD
    x
    y
    w
ENDATA
";
        let s = solve(src);
        assert_eq!(s.status, "Optimal");
        // 2xy >= 1, minimise x + y: x = y = 1/sqrt 2, objective sqrt 2
        assert!((s.objective - 2f64.sqrt()).abs() < 1e-6, "{}", s.objective);
        assert_eq!(s.detail["checker"]["verdict"], "feasible");
    }
}
