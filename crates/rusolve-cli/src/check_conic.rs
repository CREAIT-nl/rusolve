//! Solution checker for models with quadratic rows and cones: `rusolve check <model> <solution>` takes
//! this path when the LP/MIP reader rejected the model for its quadratic rows
//! or cones ([`rusolve_io::conic::conic_rejected`]).
//!
//! The checker evaluates the user's objects, never a reformulation:
//!
//! * column bounds and integrality, by the rules of [`crate::check`];
//! * linear rows, and quadratic rows `aᵀx + ½ xᵀQx` against the row's
//!   `[lower, upper]` (the model's sense, rhs and range);
//! * cone membership: SOC `x0 >= ‖x1..‖`, rotated SOC `2 x0 x1 >= ‖x2..‖²`
//!   with `x0, x1 >= 0`, both measured as the violation of the equivalent SOC
//!   `‖(x0 - x1, √2 w)‖ <= x0 + x1`;
//! * the objective `cᵀx + ½ xᵀHx + offset`.
//!
//! Tolerances: an interior-point answer is accurate to a relative amount, not an absolute one, so
//! every violation is compared with `tol · (1 + scale)` of its own item, the convention of
//! the QP checker's KKT test ([`crate::check_qp`]: "relative to `1 + max |term|` of each
//! condition"; HiGHS's `getKktFailures` and Clarabel's termination test are relative in the
//! same way). The scale is the largest term the item sums: for a column bound the solution-wide `max |x|` (Clarabel scales a bound row by the
//! global norm); for a row `max(|a_ij x_j|, Σ|quadratic terms|, |bounds|)`; for a cone `‖x_cone‖₂`.
//! One large unrelated variable does not loosen any other item. Integrality stays absolute.
//! The default tolerance is [`DEFAULT_TOL`], looser than HiGHS's 1e-7 for the same reason.
//!
//! Output: one JSON object `rusolve.check-conic/1`
//! (`docs/schema/rusolve-check-conic-1.schema.json`), or text. Exit codes as
//! [`crate::check`]: 0 the check ran, 3 input error, 4 internal.

use rusolve_io::Model;
use rusolve_io::conic::ConeKind;
use serde::Serialize;
use std::collections::HashMap;
use std::process::ExitCode;

use crate::check::{Counts, InputError, parse_solution};

/// Default primal tolerance of the conic checker (interior-point accuracy).
pub const DEFAULT_TOL: f64 = 1e-6;

pub const SCHEMA_VERSION: &str = "rusolve.check-conic/1";

#[derive(Serialize, Debug, Clone)]
pub struct Finding {
    /// `column_bound`, `row`, `quad_row`, `cone` or `integrality`.
    pub kind: &'static str,
    pub name: String,
    /// Activity of the row, value of the column, or the cone's slack
    /// `x0 - ‖x1..‖` (negative when violated).
    pub value: f64,
    pub lower: Option<f64>,
    pub upper: Option<f64>,
    /// Violation amount in the units of the item.
    pub amount: f64,
    /// Cone type, for `cone`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cone_type: Option<&'static str>,
}

#[derive(Serialize, Debug, Default)]
pub struct ConicSummary {
    pub column_bounds: Counts,
    pub rows: Counts,
    pub quad_rows: Counts,
    pub cones: Counts,
    pub integrality: Counts,
    pub total_violations: i64,
}

#[derive(Serialize, Debug)]
pub struct ConicObjective {
    /// `cᵀx + ½ xᵀHx + offset`, in the model's own sense.
    pub value: f64,
    pub offset: f64,
    /// The `½ xᵀHx` part alone.
    pub quadratic_part: f64,
    pub sense: &'static str,
}

#[derive(Serialize, Debug)]
pub struct ConicModelInfo {
    pub path: String,
    pub rows: i64,
    pub cols: i64,
    pub integer_cols: i64,
    pub quad_rows: i64,
    pub cones: i64,
}

#[derive(Serialize, Debug)]
pub struct ConicCheckReport {
    pub schema_version: &'static str,
    pub model: ConicModelInfo,
    pub solution_path: String,
    pub verdict: &'static str,
    pub tolerance: f64,
    pub objective: ConicObjective,
    pub summary: ConicSummary,
    pub findings_total: i64,
    pub findings_shown: i64,
    pub findings: Vec<Finding>,
}

fn finite(v: f64) -> Option<f64> {
    v.is_finite().then_some(v)
}

/// `[lower, upper]` of constraint row `i` from its sense, rhs and range (the
/// same reading as `Model::to_sparse_lp`).
pub fn row_interval(m: &Model, i: usize) -> (f64, f64) {
    use rusolve_io::RowSense::*;
    let (rhs, r) = (m.rhs[i], m.ranges[i]);
    match (m.senses[i], r != 0.0) {
        (LessEqual, false) => (f64::NEG_INFINITY, rhs),
        (LessEqual, true) => (rhs - r.abs(), rhs),
        (GreaterEqual, false) => (rhs, f64::INFINITY),
        (GreaterEqual, true) => (rhs, rhs + r.abs()),
        (Equal, false) => (rhs, rhs),
        (Equal, true) if r > 0.0 => (rhs, rhs + r),
        (Equal, true) => (rhs + r, rhs),
    }
}

fn name_of(names: &[String], ix: usize, prefix: &str) -> String {
    match names.get(ix) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("{prefix}{ix}"),
    }
}

/// The violation of a cone at `x` (0 when inside) and its slack
/// `x0 - ‖x1..‖` (SOC) or the equivalent SOC's slack (rotated).
pub fn cone_violation(kind: ConeKind, x: &[f64], vars: &[usize]) -> (f64, f64) {
    let v: Vec<f64> = vars.iter().map(|&j| x[j]).collect();
    match kind {
        ConeKind::Soc => {
            let norm = v[1..].iter().map(|a| a * a).sum::<f64>().sqrt();
            let slack = v[0] - norm;
            ((-slack).max(0.0), slack)
        }
        ConeKind::RotatedSoc => {
            // 2 x0 x1 >= ‖w‖², x0, x1 >= 0  <=>  ‖(x0 - x1, √2 w)‖ <= x0 + x1
            let w2: f64 = v[2..].iter().map(|a| a * a).sum();
            let norm = ((v[0] - v[1]).powi(2) + 2.0 * w2).sqrt();
            let slack = v[0] + v[1] - norm;
            ((-slack).max(0.0), slack)
        }
        // Reserved kinds are rejected by the readers; reaching here is a bug.
        ConeKind::Exp | ConeKind::Pow(_) => (f64::INFINITY, f64::NEG_INFINITY),
    }
}

/// Checks `x` (model column order) against the model. Pure: no I/O.
pub fn check_conic_solution(
    m: &Model,
    x: &[f64],
    tol: f64,
    max_findings: usize,
    model_path: &str,
    solution_path: &str,
) -> ConicCheckReport {
    let mut findings: Vec<Finding> = Vec::new();
    let mut s = ConicSummary::default();
    // Column bounds only keep the solution-wide scale: Clarabel treats a bound as the conic row
    // `-x + s = 0` and scales its primal residual by the global `max(|b|, |s|)`, so a bound
    // violation of an interior-point answer is relative to the whole solution (the sched_* CBLIB
    // answers carry 5e-6 on a near-zero column next to values of 2.7e4).
    let xscale = 1.0 + x.iter().fold(0.0f64, |m, v| m.max(v.abs()));

    for (j, &v) in x.iter().enumerate() {
        let b = &m.bounds[j];
        let amount = (b.lower - v).max(v - b.upper).max(0.0);
        if amount > tol * xscale {
            s.column_bounds.add(amount);
            findings.push(Finding {
                kind: "column_bound",
                name: name_of(&m.col_names, j, "c"),
                value: v,
                lower: finite(b.lower),
                upper: finite(b.upper),
                amount,
                cone_type: None,
            });
        }
        if b.integer {
            let d = (v - v.round()).abs();
            if d > tol {
                s.integrality.add(d);
                findings.push(Finding {
                    kind: "integrality",
                    name: name_of(&m.col_names, j, "c"),
                    value: v,
                    lower: None,
                    upper: None,
                    amount: d,
                    cone_type: None,
                });
            }
        }
    }

    // Row activities: the linear part from the CSC, then the quadratic part.
    let mut act = vec![0.0f64; m.num_constraints];
    // largest term of each row's activity, for the per-row tolerance scale
    let mut term = vec![0.0f64; m.num_constraints];
    let (cp, ri, va) = (m.matrix.col_ptr(), m.matrix.row_idx(), m.matrix.values());
    for (j, &xj) in x.iter().enumerate() {
        for k in cp[j]..cp[j + 1] {
            let t = va[k] * xj;
            act[ri[k] as usize] += t;
            term[ri[k] as usize] = term[ri[k] as usize].max(t.abs());
        }
    }
    let mut quad_idx = vec![usize::MAX; m.num_constraints];
    for (k, q) in m.conic.quad_rows.iter().enumerate() {
        act[q.row] += q.value(x);
        let abs_sum: f64 = q.monomials().map(|(i, j, c)| (c * x[i] * x[j]).abs()).sum();
        term[q.row] = term[q.row].max(abs_sum);
        quad_idx[q.row] = k;
    }
    for (i, &a) in act.iter().enumerate() {
        let (lo, hi) = row_interval(m, i);
        let amount = (lo - a).max(a - hi).max(0.0);
        let is_quad = quad_idx[i] != usize::MAX;
        let scale = [lo, hi]
            .into_iter()
            .filter(|v| v.is_finite())
            .fold(term[i], |s, v| s.max(v.abs()));
        let limit = tol * (1.0 + scale);
        let counts = if is_quad {
            &mut s.quad_rows
        } else {
            &mut s.rows
        };
        if amount > limit {
            counts.add(amount);
            findings.push(Finding {
                kind: if is_quad { "quad_row" } else { "row" },
                name: name_of(&m.row_names, i, "r"),
                value: a,
                lower: finite(lo),
                upper: finite(hi),
                amount,
                cone_type: None,
            });
        }
    }

    for c in &m.conic.cones {
        let (viol, slack) = cone_violation(c.kind, x, &c.vars);
        let norm = c.vars.iter().map(|&j| x[j] * x[j]).sum::<f64>().sqrt();
        if viol > tol * (1.0 + norm) || viol.is_nan() {
            s.cones.add(viol);
            findings.push(Finding {
                kind: "cone",
                name: c.name.clone(),
                value: slack,
                lower: None,
                upper: None,
                amount: viol,
                cone_type: Some(c.kind.label()),
            });
        }
    }
    s.total_violations = s.column_bounds.count
        + s.rows.count
        + s.quad_rows.count
        + s.cones.count
        + s.integrality.count;

    let mut linear = m.obj_offset;
    for (j, &xj) in x.iter().enumerate() {
        linear += m.objective[j] * xj;
    }
    let quadratic_part: f64 = 0.5
        * m.hessian
            .iter()
            .map(|&(i, j, v)| v * x[i] * x[j])
            .sum::<f64>();
    findings.sort_by(|a, b| {
        b.amount
            .total_cmp(&a.amount)
            .then_with(|| a.kind.cmp(b.kind))
            .then_with(|| a.name.cmp(&b.name))
    });
    let findings_total = findings.len() as i64;
    findings.truncate(max_findings);
    ConicCheckReport {
        schema_version: SCHEMA_VERSION,
        model: ConicModelInfo {
            path: model_path.to_string(),
            rows: m.num_constraints as i64,
            cols: m.num_variables as i64,
            integer_cols: m.bounds.iter().filter(|b| b.integer).count() as i64,
            quad_rows: m.conic.quad_rows.len() as i64,
            cones: m.conic.cones.len() as i64,
        },
        solution_path: solution_path.to_string(),
        verdict: if s.total_violations == 0 {
            "feasible"
        } else {
            "infeasible"
        },
        tolerance: tol,
        objective: ConicObjective {
            value: linear + quadratic_part,
            offset: m.obj_offset,
            quadratic_part,
            sense: match m.obj_sense {
                rusolve_io::ObjSense::Minimize => "minimize",
                rusolve_io::ObjSense::Maximize => "maximize",
            },
        },
        summary: s,
        findings_total,
        findings_shown: findings.len() as i64,
        findings,
    }
}

/// Matches `pairs` to the model's columns by name; a missing, unknown,
/// repeated or NaN column is an input error naming it.
pub fn values_by_name(m: &Model, pairs: &[(String, f64)]) -> Result<Vec<f64>, InputError> {
    let n = m.num_variables;
    let names: Vec<String> = (0..n).map(|i| name_of(&m.col_names, i, "c")).collect();
    let index: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut values: Vec<Option<f64>> = vec![None; n];
    for (name, value) in pairs {
        let Some(&ix) = index.get(name.as_str()) else {
            return Err(InputError(format!(
                "column {name} in the solution file is not a column of the model"
            )));
        };
        if values[ix].is_some() {
            return Err(InputError(format!(
                "column {name} appears more than once in the solution file"
            )));
        }
        if value.is_nan() {
            return Err(InputError(format!("column {name} has the value NaN")));
        }
        values[ix] = Some(*value);
    }
    values
        .into_iter()
        .enumerate()
        .map(|(ix, v)| {
            v.ok_or_else(|| {
                InputError(format!(
                    "column {} of the model is missing from the solution file",
                    names[ix]
                ))
            })
        })
        .collect()
}

/// The human-readable report.
pub fn render_text(r: &ConicCheckReport) -> String {
    let mut o = format!(
        "verdict: {}\nobjective: {} ({}, quadratic part {})\nviolations: {} (bounds {}, rows {}, quadratic rows {}, cones {}, integrality {})\n",
        r.verdict,
        r.objective.value,
        r.objective.sense,
        r.objective.quadratic_part,
        r.summary.total_violations,
        r.summary.column_bounds.count,
        r.summary.rows.count,
        r.summary.quad_rows.count,
        r.summary.cones.count,
        r.summary.integrality.count,
    );
    for f in &r.findings {
        o.push_str(&format!(
            "  {} {}: value {} violated by {}\n",
            f.kind, f.name, f.value, f.amount
        ));
    }
    if r.findings_total > r.findings_shown {
        o.push_str(&format!(
            "  ... {} more\n",
            r.findings_total - r.findings_shown
        ));
    }
    o
}

/// `rusolve check` for a conic model: the model has already been rejected by
/// the LP/MIP reader. `tol` is the primal feasibility tolerance.
pub fn run(
    model_path: &str,
    solution_path: &str,
    json: bool,
    max_findings: usize,
    tol: Option<f64>,
) -> ExitCode {
    let model = match rusolve_io::read_model_file_conic(model_path) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("rusolve check: read {model_path}: {e}");
            return ExitCode::from(3);
        }
    };
    let text = match std::fs::read_to_string(solution_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("rusolve check: read {solution_path}: {e}");
            return ExitCode::from(3);
        }
    };
    let values = parse_solution(&text, solution_path.ends_with(".json"))
        .and_then(|pairs| values_by_name(&model, &pairs));
    let x = match values {
        Ok(x) => x,
        Err(InputError(e)) => {
            eprintln!("rusolve check: {solution_path}: {e}");
            return ExitCode::from(3);
        }
    };
    let tol = tol.unwrap_or(DEFAULT_TOL);
    let report = check_conic_solution(&model, &x, tol, max_findings, model_path, solution_path);
    if json {
        match serde_json::to_string(&report) {
            Ok(s) => println!("{s}"),
            Err(e) => {
                eprintln!("rusolve check: cannot serialise: {e}");
                return ExitCode::from(4);
            }
        }
    } else {
        print!("{}", render_text(&report));
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_io::read_mps_conic;

    /// min x + y  s.t. q1: x + [2x² + 3xy + y²] <= 5, l1: x + y >= 1,
    /// cone k1: t >= ‖(x, y)‖, rotated cone k2: 2 x y >= t², x, y in [0, 10].
    const M: &str = "\
NAME planted
ROWS
 N  obj
 L  q1
 G  l1
COLUMNS
    x  obj  1  q1  1
    x  l1  1
    y  obj  1  l1  1
    t  obj  0
RHS
    rhs  q1  5  l1  1
BOUNDS
 UP bnd x 10
 UP bnd y 10
 UP bnd t 100
QCMATRIX q1
    x  x  2
    x  y  1.5
    y  x  1.5
    y  y  1
CSECTION k1 0.0 QUAD
    t
    x
    y
CSECTION k2 0.0 RQUAD
    x
    y
    t
ENDATA
";

    fn check(x: [f64; 3]) -> ConicCheckReport {
        let m = read_mps_conic(M).unwrap();
        check_conic_solution(&m, &x, 1e-7, 50, "m", "s")
    }

    fn kinds(r: &ConicCheckReport) -> Vec<(&'static str, String)> {
        let mut k: Vec<_> = r
            .findings
            .iter()
            .map(|f| (f.kind, f.name.clone()))
            .collect();
        k.sort();
        k
    }

    #[test]
    fn a_feasible_point_passes_and_the_objective_is_recomputed() {
        // x = y = 1: quad row 1 + 2 + 3 + 1 = 7 > 5, so shrink: x = y = 0.5.
        // row: 0.5 + 0.5 + 0.75 + 0.25 = 2 <= 5; l1: 1 >= 1; t = 0.5:
        // SOC 0.5 >= ‖(.5,.5)‖ = 0.707 violated, so make t big enough for k1
        // and small enough for k2: t = 0.7 (k1 needs t >= 0.7071).
        let r = check([0.5, 0.5, 0.71]);
        assert_eq!(kinds(&r), vec![("cone", "k2".to_string())], "{r:?}");
        // k2: 2 x y = 0.5 >= t² = 0.5041 fails; t = 0.7071 is the boundary.
        let r = check([0.5, 0.5, std::f64::consts::FRAC_1_SQRT_2]);
        assert_eq!(r.verdict, "feasible", "{r:?}");
        assert!((r.objective.value - 1.0).abs() < 1e-12);
    }

    #[test]
    fn planted_quadratic_row_violation_is_named_with_its_amount() {
        let r = check([1.0, 1.0, 0.0]);
        let q = r.findings.iter().find(|f| f.kind == "quad_row").unwrap();
        assert_eq!(q.name, "q1");
        assert!((q.value - 7.0).abs() < 1e-12);
        assert!((q.amount - 2.0).abs() < 1e-12);
    }

    #[test]
    fn planted_cone_violations_are_named() {
        // t = 0: k1 fails (needs t >= ‖(x,y)‖), k2 holds (2xy >= 0).
        let r = check([0.5, 0.5, 0.0]);
        assert_eq!(kinds(&r), vec![("cone", "k1".to_string())]);
        let c = &r.findings[0];
        assert_eq!(c.cone_type, Some("SOC"));
        assert!((c.amount - 0.5f64.hypot(0.5)).abs() < 1e-12);
        // x negative: k2 needs x >= 0.
        let r = check([-0.1, 0.5, 0.0]);
        assert!(
            r.findings
                .iter()
                .any(|f| f.kind == "cone" && f.name == "k2")
        );
        assert!(
            r.findings
                .iter()
                .any(|f| f.kind == "column_bound" && f.name == "x")
        );
    }

    #[test]
    fn a_linear_row_violation_is_still_a_row() {
        let r = check([0.1, 0.1, 0.2]);
        assert!(kinds(&r).contains(&("row", "l1".to_string())), "{r:?}");
    }

    #[test]
    fn objective_includes_the_hessian_part() {
        let mut m = read_mps_conic(M).unwrap();
        m.hessian = vec![(0, 0, 4.0), (1, 1, 2.0)];
        let r = check_conic_solution(&m, &[1.0, 2.0, 0.0], 1e-7, 5, "m", "s");
        // c'x = 3, ½(4·1 + 2·4) = 6
        assert!((r.objective.quadratic_part - 6.0).abs() < 1e-12);
        assert!((r.objective.value - 9.0).abs() < 1e-12);
    }

    /// The old rule scaled every tolerance by `1 + max|x|` over the whole solution, so one huge
    /// unrelated column accepted a violated row. A big column `h` is added to `M`; the point
    /// `x = y = 1` violates q1 by 2 (activity 7 vs 5) and `l1`-style rows must still be judged
    /// on their own terms.
    #[test]
    fn one_huge_unrelated_variable_does_not_loosen_other_rows() {
        let src = M
            .replace("    t  obj  0\n", "    t  obj  0\n    h  obj  0\n")
            .replace(" UP bnd t 100\n", " UP bnd t 100\n UP bnd h 1e9\n");
        let m = read_mps_conic(&src).unwrap();
        // columns: x, y, t, h. h = 1e9 sits in no row and no cone.
        let x = [1.0, 1.0, 0.0, 1.0e9];
        let r = check_conic_solution(&m, &x, 1e-7, 50, "m", "s");
        assert!(
            r.findings
                .iter()
                .any(|f| f.kind == "quad_row" && f.name == "q1"),
            "q1 violated by 2 must be reported: {r:?}"
        );
        // the old global scale: tol * (1 + 1e9) = 100 would have accepted amount 2
        let old_limit = 1e-7 * (1.0 + x.iter().fold(0.0f64, |a, v| a.max(v.abs())));
        assert!(2.0 < std::hint::black_box(old_limit));
        // a violation just above the item's own tolerance is caught, just below is not
        let tiny = |d: f64| {
            // t = ‖(x, y)‖ + d on k1, so only the SOC k1 is off, by d
            let x = [0.5, 0.5, 0.5f64.hypot(0.5) - d, 1.0e9];
            check_conic_solution(&m, &x, 1e-7, 50, "m", "s")
        };
        assert!(
            tiny(1e-5).findings.iter().any(|f| f.name == "k1"),
            "1e-5 off a cone of norm 0.7"
        );
        assert!(!tiny(1e-9).findings.iter().any(|f| f.name == "k1"));
    }

    #[test]
    fn solution_columns_are_matched_by_name() {
        let m = read_mps_conic(M).unwrap();
        let ok = values_by_name(
            &m,
            &[("t".into(), 3.0), ("x".into(), 1.0), ("y".into(), 2.0)],
        )
        .unwrap();
        assert_eq!(ok, vec![1.0, 2.0, 3.0]);
        assert!(values_by_name(&m, &[("x".into(), 1.0)]).is_err());
    }
}
