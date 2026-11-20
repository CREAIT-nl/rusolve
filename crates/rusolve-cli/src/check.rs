//! `rusolve check <model.mps> <solution> [--format json] [--max-findings N]
//! [--primal-feasibility-tolerance X] [--mip-feasibility-tolerance X]`.
//!
//! Tells an agent, by name, everything wrong with a candidate solution: every
//! violated column bound, row and integrality, each with its amount, the
//! objective recomputed with the model offset, and a verdict.
//!
//! ## Where the numbers come from
//!
//! The verdict and the tolerances are the ported HiGHS checks, so the checker
//! agrees with HiGHS by construction; only the per-item listing is ours:
//!
//! * `assess_lp_primal_solution` (HiGHS 1.16.0 lp_data/HighsLpUtils.cpp:2623-2769,
//!   `rusolve_ekk::lp_data::highs_complete_solution`, the function behind
//!   `Highs::assessPrimalSolution`, Highs.cpp:3681) gives the three booleans
//!   `valid`/`integral`/`feasible`. Its tolerance is `mip_feasibility_tolerance`
//!   for a model with integer columns and `primal_feasibility_tolerance` for an
//!   LP (HighsLpUtils.cpp:2641-2643): [`Tolerances::applied`].
//! * `assess_col_primal_solution` (HighsLpUtils.cpp:2591-2619, same file of the
//!   port) gives each column's bound and integrality infeasibility; the row
//!   infeasibility is the `row_infeasibility` rule of
//!   `assess_lp_primal_solution` (primal < lower - tol, primal > upper + tol).
//!   An item is listed exactly when `assess_lp_primal_solution` counts it
//!   (infeasibility above the applied tolerance), so the listing and the
//!   verdict cannot disagree: [`CheckReport::agrees_with_highs`] records that.
//! * `get_kkt_failures` (HighsSolution.cpp:73-561, `rusolve_ekk::lp_data::
//!   presolve_support`) is run as `Highs::callSolveMip` runs it
//!   (`Highs::call_solve_mip`, `highs_run.rs:1138-1150`:
//!   `primal_feasibility_tolerance` temporarily set to
//!   `mip_feasibility_tolerance`) and its `num/max/sum_primal_infeasibilities`
//!   are reported next to ours as `highs.kkt_*`.
//!
//! ## The solution file
//!
//! HiGHS has `readSolutionFile` (HighsLpUtils.cpp:2087-2472) but the port does
//! not, and that reader leaves a column missing from a sparse file at 0 and
//! reports an unknown name only through `getIndexFromName`. The reader here
//! follows its raw-style layout (`Model status` / status / blank /
//! `# Primal solution values` / status / `Objective v` / optional
//! `Objective_bound`, `MIP_rel_gap` / `# Columns N` then `name value` lines,
//! HighsLpUtils.cpp:2160-2290) and makes a missing,
//! unknown or repeated column an input error naming it. The row section and
//! everything after it are ignored: rows are recomputed from the columns. The
//! JSON form is `{"columns":[{"name","value"}]}`, also nested as
//! `{"solution":{"columns":[...]}}` (the `rusolve.result/1` envelope).
//!
//! ## Exit codes (as `mip_highs`)
//!
//! 0 the check ran (an `infeasible` verdict is a result, not an error); 2 usage
//! error; 3 input error (model or solution unreadable, a missing, unknown or
//! duplicate column, a value that is not a number); 4 internal error.

use rusolve_core::diagnostic::Diagnostic;
use rusolve_ekk::lp_data::highs_complete_solution::{
    assess_col_primal_solution, assess_lp_primal_solution,
};
use rusolve_ekk::lp_data::presolve_support::{calculate_row_values_quad_vec, get_kkt_failures};
use rusolve_ekk::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::process::ExitCode;

/// The `schema_version` of the output (`docs/schema/rusolve-check-1.schema.json`).
pub const SCHEMA_VERSION: &str = "rusolve.check/1";

/// Findings listed when `--max-findings` is not given.
pub(crate) const DEFAULT_MAX_FINDINGS: usize = 50;

/// An input error: exit code 3.
#[derive(Debug)]
pub struct InputError(pub String);

/// The tolerances in force (HiGHS defaults unless overridden on the command line).
#[derive(Serialize, Clone, Copy, Debug)]
pub struct Tolerances {
    pub primal_feasibility_tolerance: f64,
    pub mip_feasibility_tolerance: f64,
    /// The one `assess_lp_primal_solution` compares infeasibilities with:
    /// `mip_feasibility_tolerance` if the model has an integer column, else
    /// `primal_feasibility_tolerance` (HighsLpUtils.cpp:2641-2643).
    pub applied: f64,
}

#[derive(Serialize, Debug, Clone)]
pub struct Finding {
    /// `column_bound`, `row`, `integrality`, or `semi_continuous`,
    /// `indicator`, `sos`, `pwl`, `function`.
    pub kind: &'static str,
    pub name: String,
    /// Column value (`column_bound`, `integrality`) or row activity (`row`).
    pub value: f64,
    /// `null` when the bound is infinite or the finding is an integrality one.
    pub lower: Option<f64>,
    pub upper: Option<f64>,
    /// The bound that is violated (`column_bound`, `row`).
    pub violated_bound: Option<f64>,
    /// `integrality` only.
    pub nearest_integer: Option<f64>,
    /// Violation amount: distance past the violated bound, or distance to the
    /// nearest integer.
    pub amount: f64,
    /// General-constraint kinds only: what is violated, in words.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

#[derive(Serialize, Default, Debug, Clone, Copy)]
pub struct Counts {
    pub count: i64,
    pub max: f64,
    pub sum: f64,
}

impl Counts {
    pub(crate) fn add(&mut self, amount: f64) {
        self.count += 1;
        self.max = self.max.max(amount);
        self.sum += amount;
    }
}

#[derive(Serialize, Debug)]
pub struct Summary {
    pub column_bounds: Counts,
    pub rows: Counts,
    pub integrality: Counts,
    /// Semi-continuous membership and general constraints; absent for a
    /// model with none.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub general_constraints: Option<Counts>,
    pub total_violations: i64,
}

/// The ported HiGHS verdicts, for cross-checking.
#[derive(Serialize, Debug)]
pub struct HighsAgreement {
    /// `assess_lp_primal_solution`'s three booleans.
    pub valid: bool,
    pub integral: bool,
    pub feasible: bool,
    /// `get_kkt_failures`' primal infeasibility count, maximum and sum
    /// (columns and rows, at the applied tolerance).
    pub kkt_num_primal_infeasibilities: i64,
    pub kkt_max_primal_infeasibility: f64,
    pub kkt_sum_primal_infeasibilities: f64,
}

#[derive(Serialize, Debug)]
pub struct Objective {
    /// `sum(cost * value) + offset`, the number HiGHS reports (the model's own
    /// sense, not flipped).
    pub value: f64,
    pub offset: f64,
    pub sense: &'static str,
}

#[derive(Serialize, Debug)]
pub struct ModelInfo {
    pub path: String,
    pub rows: i64,
    pub cols: i64,
    pub integer_cols: i64,
}

#[derive(Serialize, Debug)]
pub struct CheckReport {
    pub schema_version: &'static str,
    pub model: ModelInfo,
    pub solution_path: String,
    pub verdict: &'static str,
    pub tolerances: Tolerances,
    pub objective: Objective,
    pub summary: Summary,
    pub highs: HighsAgreement,
    /// `true` when our per-item counts equal what `assess_lp_primal_solution`
    /// and `get_kkt_failures` counted. `false` would be a checker bug.
    pub agrees_with_highs: bool,
    /// All findings before the `--max-findings` cut.
    pub findings_total: i64,
    pub findings_shown: i64,
    /// Largest amount first, cut at `--max-findings`.
    pub findings: Vec<Finding>,
    /// Convex QP: present only when the model has a Hessian. Then `objective.value`
    /// includes `x'Qx/2` and `qp` carries the KKT and duality-gap measures.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qp: Option<crate::check_qp::QpCheck>,
}

fn name_of(names: &[String], ix: usize, prefix: &str) -> String {
    match names.get(ix) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("{prefix}{ix}"),
    }
}

fn finite(v: f64) -> Option<f64> {
    (v.abs() < kHighsInf).then_some(v)
}

/// Checks `col_value` (original column order) against `lp`. Pure: no I/O.
pub fn check_solution(
    lp: &HighsLp,
    col_value: &[f64],
    tol: Tolerances,
    max_findings: usize,
    model_path: &str,
    solution_path: &str,
) -> Result<CheckReport, InputError> {
    check_solution_gc(
        lp,
        None,
        col_value,
        tol,
        max_findings,
        model_path,
        solution_path,
    )
}

/// [`check_solution`] plus the general-constraint objects: semi-continuous membership and
/// indicator, SOS, PWL and function constraints, evaluated on the user's
/// objects (rounded binaries), never on a linearization.
pub fn check_solution_gc(
    lp: &HighsLp,
    gc: Option<&rusolve_io::GcModel>,
    col_value: &[f64],
    tol: Tolerances,
    max_findings: usize,
    model_path: &str,
    solution_path: &str,
) -> Result<CheckReport, InputError> {
    let gc = gc.filter(|g| !g.is_empty());
    // The ported assessors log through `log_options`; silence them (stdout is
    // the report).
    let mut options = HighsOptions {
        output_flag: false,
        primal_feasibility_tolerance: tol.primal_feasibility_tolerance,
        mip_feasibility_tolerance: tol.mip_feasibility_tolerance,
        ..Default::default()
    };
    options.log_options.output_flag = false;
    options.log_options.log_to_console = false;

    let mut solution = HighsSolution {
        col_value: col_value.to_vec(),
        ..Default::default()
    };
    if calculate_row_values_quad_vec(lp, col_value, &mut solution.row_value, -1) != HighsStatus::kOk
    {
        return Err(InputError("could not compute the row activities".into()));
    }
    solution.value_valid = true;

    // The ported booleans.
    let (mut valid, mut integral, mut feasible) = (false, false, false);
    assess_lp_primal_solution(
        "",
        &options,
        lp,
        &solution,
        &mut valid,
        &mut integral,
        &mut feasible,
    );

    // The per-item listing, by the same rules and tolerance.
    let have_integrality = !lp.integrality_.is_empty();
    let mut findings: Vec<Finding> = Vec::new();
    let mut summary = Summary {
        column_bounds: Counts::default(),
        rows: Counts::default(),
        integrality: Counts::default(),
        general_constraints: None,
        total_violations: 0,
    };
    for (ix, &value) in col_value.iter().enumerate() {
        let (lower, upper) = (lp.col_lower_[ix], lp.col_upper_[ix]);
        let ty = if have_integrality {
            lp.integrality_[ix]
        } else {
            HighsVarType::kContinuous
        };
        let (mut col_inf, mut int_inf) = (0.0, 0.0);
        assess_col_primal_solution(
            &options,
            value,
            lower,
            upper,
            ty,
            &mut col_inf,
            &mut int_inf,
        );
        if col_inf > tol.applied {
            summary.column_bounds.add(col_inf);
            findings.push(Finding {
                kind: "column_bound",
                name: name_of(&lp.col_names_, ix, "c"),
                value,
                lower: finite(lower),
                upper: finite(upper),
                violated_bound: Some(if value < lower { lower } else { upper }),
                nearest_integer: None,
                amount: col_inf,
                detail: None,
            });
        }
        // A semi column's membership is reported by the gc checker below
        let semi = gc.is_some_and(|g| g.is_semi(ix));
        if int_inf > options.mip_feasibility_tolerance && !semi {
            summary.integrality.add(int_inf);
            findings.push(Finding {
                kind: "integrality",
                name: name_of(&lp.col_names_, ix, "c"),
                value,
                lower: None,
                upper: None,
                violated_bound: None,
                nearest_integer: Some(value.round()),
                amount: int_inf,
                detail: None,
            });
        }
    }
    for ix in 0..lp.num_row_ as usize {
        let activity = solution.row_value[ix];
        let (lower, upper) = (lp.row_lower_[ix], lp.row_upper_[ix]);
        // HighsLpUtils.cpp assessLpPrimalSolution, row loop
        let mut row_inf = 0.0;
        if activity < lower - tol.applied {
            row_inf = lower - activity;
        } else if activity > upper + tol.applied {
            row_inf = activity - upper;
        }
        if row_inf > tol.applied {
            summary.rows.add(row_inf);
            findings.push(Finding {
                kind: "row",
                name: name_of(&lp.row_names_, ix, "r"),
                value: activity,
                lower: finite(lower),
                upper: finite(upper),
                violated_bound: Some(if activity < lower { lower } else { upper }),
                nearest_integer: None,
                amount: row_inf,
                detail: None,
            });
        }
    }
    let port_violations =
        summary.column_bounds.count + summary.rows.count + summary.integrality.count;
    summary.total_violations = port_violations;
    if let Some(g) = gc {
        let gtol = rusolve_io::GcTol {
            integrality: tol.mip_feasibility_tolerance,
            feasibility: tol.applied,
        };
        let mut counts = Counts::default();
        for v in g.check(
            col_value,
            &lp.col_lower_,
            &lp.col_upper_,
            &lp.col_names_,
            gtol,
        ) {
            counts.add(v.amount);
            findings.push(Finding {
                kind: v.kind,
                name: v.name,
                value: v.amount,
                lower: None,
                upper: None,
                violated_bound: None,
                nearest_integer: None,
                amount: v.amount,
                detail: Some(v.detail),
            });
        }
        summary.total_violations += counts.count;
        summary.general_constraints = Some(counts);
    }

    // get_kkt_failures as call_solve_mip runs it (highs_run.rs:1138-1150).
    let mut kkt_options = options.clone();
    kkt_options.primal_feasibility_tolerance = tol.applied;
    let mut info = HighsInfo::default();
    get_kkt_failures(
        &kkt_options,
        false,
        lp,
        &lp.col_cost_,
        &solution,
        &mut info,
        false,
    );

    // The objective: HiGHS writeLpObjective, `offset + sum cost * value`.
    let mut obj = HighsCDouble::from(lp.offset_);
    for (ix, &value) in col_value.iter().enumerate() {
        obj += lp.col_cost_[ix] * value;
    }

    findings.sort_by(|a, b| {
        b.amount
            .total_cmp(&a.amount)
            .then_with(|| a.kind.cmp(b.kind))
            .then_with(|| a.name.cmp(&b.name))
    });
    let findings_total = findings.len() as i64;
    findings.truncate(max_findings);

    let agrees = (valid && feasible && integral) == (port_violations == 0)
        && info.num_primal_infeasibilities as i64
            == summary.column_bounds.count + summary.rows.count;
    Ok(CheckReport {
        schema_version: SCHEMA_VERSION,
        model: ModelInfo {
            path: model_path.to_string(),
            rows: lp.num_row_ as i64,
            cols: lp.num_col_ as i64,
            integer_cols: lp
                .integrality_
                .iter()
                .filter(|t| **t == HighsVarType::kInteger)
                .count() as i64,
        },
        solution_path: solution_path.to_string(),
        verdict: if summary.total_violations == 0 {
            "feasible"
        } else {
            "infeasible"
        },
        tolerances: tol,
        objective: Objective {
            value: f64::from(obj),
            offset: lp.offset_,
            sense: if lp.sense_ == ObjSense::kMaximize {
                "maximize"
            } else {
                "minimize"
            },
        },
        summary,
        highs: HighsAgreement {
            valid,
            integral,
            feasible,
            kkt_num_primal_infeasibilities: info.num_primal_infeasibilities as i64,
            kkt_max_primal_infeasibility: info.max_primal_infeasibility,
            kkt_sum_primal_infeasibilities: info.sum_primal_infeasibilities,
        },
        agrees_with_highs: agrees,
        findings_total,
        findings_shown: findings.len() as i64,
        findings,
        qp: None,
    })
}

/// Parses a solution file's text into `(name, value)` pairs in file order: the
/// HiGHS raw layout (see the module doc) or, when `json` is set, the JSON form.
pub fn parse_solution(text: &str, json: bool) -> Result<Vec<(String, f64)>, InputError> {
    if json {
        let v: serde_json::Value = serde_json::from_str(text)
            .map_err(|e| InputError(format!("solution file is not valid JSON: {e}")))?;
        let cols = v
            .get("columns")
            .or_else(|| v.get("solution").and_then(|s| s.get("columns")))
            .and_then(|col| col.as_array())
            .ok_or_else(|| InputError("solution JSON has no \"columns\" array".into()))?;
        let mut out = Vec::with_capacity(cols.len());
        for (i, col) in cols.iter().enumerate() {
            let name = col.get("name").and_then(|n| n.as_str()).ok_or_else(|| {
                InputError(format!(
                    "solution JSON: columns[{i}] has no string \"name\""
                ))
            })?;
            let value = col.get("value").and_then(|n| n.as_f64()).ok_or_else(|| {
                InputError(format!(
                    "solution JSON: column {name}: \"value\" is not a number"
                ))
            })?;
            out.push((name.to_string(), value));
        }
        return Ok(out);
    }
    // HighsLpUtils.cpp:2130-2175: "Model status", its value, blank line,
    // "# Primal solution values", then the status keyword ("None" = no primal
    // solution, a warning return in HiGHS, an input error here).
    let mut lines = text.lines();
    if !lines
        .by_ref()
        .any(|l| l.trim() == "# Primal solution values")
    {
        return Err(InputError(
            "solution file has no \"# Primal solution values\" section".into(),
        ));
    }
    let status = lines.by_ref().map(str::trim).find(|t| !t.is_empty());
    if status == Some("None") || status.is_none() {
        return Err(InputError(
            "solution file has no primal solution (\"None\")".into(),
        ));
    }
    // Objective, Objective_bound, MIP_rel_gap lines are skipped until the
    // `# Columns N` header (HighsLpUtils.cpp:2183-2205).
    let mut declared: Option<i64> = None;
    for line in lines.by_ref() {
        let mut it = line.split_whitespace();
        if it.next() == Some("#") && it.next() == Some("Columns") {
            declared = it.next().and_then(|n| n.parse::<i64>().ok());
            if declared.is_none() {
                return Err(InputError(format!(
                    "solution file: bad column header \"{line}\""
                )));
            }
            break;
        }
    }
    let Some(declared) = declared else {
        return Err(InputError(
            "solution file has no \"# Columns N\" line".into(),
        ));
    };
    // A negative count is HiGHS's sparse/partial style (:2216-2221).
    let n = declared.unsigned_abs() as usize;
    let mut out = Vec::with_capacity(n);
    for line in lines.by_ref().take(n) {
        let mut it = line.split_whitespace();
        let (Some(name), Some(v)) = (it.next(), it.next()) else {
            return Err(InputError(format!(
                "solution file: bad column line \"{line}\""
            )));
        };
        let value: f64 = v.parse().map_err(|_| {
            InputError(format!(
                "solution file: column {name}: \"{v}\" is not a number"
            ))
        })?;
        out.push((name.to_string(), value));
    }
    if out.len() != n {
        return Err(InputError(format!(
            "solution file declares {n} columns but lists {}",
            out.len()
        )));
    }
    Ok(out)
}

/// Matches `pairs` to the model's columns by name (the C++ `getIndexFromName`
/// of HighsLpUtils.cpp:2243). A column missing from the file, a name the model
/// does not have, a repeated name and a NaN value are input errors naming the
/// column.
pub fn columns_by_name(lp: &HighsLp, pairs: &[(String, f64)]) -> Result<Vec<f64>, InputError> {
    let n = lp.num_col_ as usize;
    let names: Vec<String> = (0..n).map(|i| name_of(&lp.col_names_, i, "c")).collect();
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
    let mut out = Vec::with_capacity(n);
    for (ix, v) in values.into_iter().enumerate() {
        match v {
            Some(v) => out.push(v),
            None => {
                return Err(InputError(format!(
                    "column {} of the model is missing from the solution file",
                    names[ix]
                )));
            }
        }
    }
    Ok(out)
}

fn bound_text(v: Option<f64>, if_none: &str) -> String {
    v.map_or(if_none.to_string(), |v| format!("{v}"))
}

/// The human-readable report.
pub fn render_text(rep: &CheckReport) -> String {
    let mut s = String::new();
    s.push_str(&format!(
        "verdict: {} (model {}: {} rows, {} cols, {} integer; tolerance {})\n",
        rep.verdict,
        rep.model.path,
        rep.model.rows,
        rep.model.cols,
        rep.model.integer_cols,
        rep.tolerances.applied
    ));
    s.push_str(&format!(
        "objective ({}): {} (offset {})\n",
        rep.objective.sense, rep.objective.value, rep.objective.offset
    ));
    let line = |label: &str, k: &Counts| {
        format!(
            "{label:<14} count {:>6}  max {}  sum {}\n",
            k.count, k.max, k.sum
        )
    };
    s.push_str(&line("column bounds", &rep.summary.column_bounds));
    s.push_str(&line("rows", &rep.summary.rows));
    s.push_str(&line("integrality", &rep.summary.integrality));
    if let Some(k) = &rep.summary.general_constraints {
        s.push_str(&line("general", k));
    }
    if let Some(q) = &rep.qp {
        s.push_str(&format!(
            "qp: objective {} + quadratic {}; kkt {}",
            q.objective_linear, q.objective_quadratic, q.kkt_verdict
        ));
        if let (Some(st), Some(gap)) = (q.stationarity, q.duality_gap_rel) {
            s.push_str(&format!(
                " (stationarity {st:.3e}, dual infeasibility {:.3e}, complementarity {:.3e}, relative gap {gap:.3e})",
                q.dual_infeasibility.unwrap_or(0.0),
                q.complementarity.unwrap_or(0.0)
            ));
        }
        s.push('\n');
    }
    for f in &rep.findings {
        match f.kind {
            "semi_continuous" | "indicator" | "sos" | "pwl" | "function" => s.push_str(&format!(
                "  {:<12} {}: {} (by {})\n",
                f.kind.replace('_', " "),
                f.name,
                f.detail.as_deref().unwrap_or(""),
                f.amount
            )),
            "integrality" => s.push_str(&format!(
                "  integrality  {}: value {} is {} from {}\n",
                f.name,
                f.value,
                f.amount,
                f.nearest_integer.unwrap_or(f.value.round())
            )),
            "row" => s.push_str(&format!(
                "  row          {}: activity {} violates [{}, {}] by {}\n",
                f.name,
                f.value,
                bound_text(f.lower, "-inf"),
                bound_text(f.upper, "inf"),
                f.amount
            )),
            _ => s.push_str(&format!(
                "  column bound {}: value {} violates [{}, {}] by {}\n",
                f.name,
                f.value,
                bound_text(f.lower, "-inf"),
                bound_text(f.upper, "inf"),
                f.amount
            )),
        }
    }
    if rep.findings_shown < rep.findings_total {
        s.push_str(&format!(
            "  ... {} more not shown (--max-findings)\n",
            rep.findings_total - rep.findings_shown
        ));
    }
    s
}

fn usage(json: bool, msg: &str) -> ExitCode {
    crate::diag::fail(
        json,
        Diagnostic::error("RS0006", format!("rusolve check: {msg}")).with_help(
            "usage: rusolve check <model.mps> <solution.sol|solution.json> [--format text|json] \
         [--max-findings N] [--primal-feasibility-tolerance X] [--mip-feasibility-tolerance X]",
        ),
        2,
    )
}

fn input_error(json: bool, msg: &str) -> ExitCode {
    crate::diag::fail_msg(json, "RS0010", 3, format!("rusolve check: {msg}"))
}

/// A solution file that was read but does not parse or does not match the model (RS0011).
fn parse_error(json: bool, msg: &str) -> ExitCode {
    crate::diag::fail_msg(json, "RS0011", 3, format!("rusolve check: {msg}"))
}

/// For a QP: re-read the Hessian, recompute the objective with it and measure the KKT
/// conditions when the raw solution file carries duals ([`crate::check_qp`]).
fn add_qp(
    rep: &mut CheckReport,
    lp: &HighsLp,
    model_path: &str,
    text: &str,
    is_json: bool,
) -> Result<(), String> {
    let (model, _) = rusolve_io::read_model_file_order(model_path)
        .map_err(|e| format!("read {model_path}: {e}"))?;
    let values = parse_solution(text, is_json)
        .and_then(|p| columns_by_name(lp, &p))
        .map_err(|InputError(e)| e)?;
    let duals = if is_json {
        None
    } else {
        crate::check_qp::parse_duals(text, lp)
    };
    let q = crate::check_qp::check_qp(lp, &model.hessian, &values, duals.as_ref());
    rep.objective.value = q.objective_linear + q.objective_quadratic;
    rep.qp = Some(q);
    Ok(())
}

/// `rusolve check`: `args` are the words after the subcommand.
pub fn run_check(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let mut positional: Vec<&String> = Vec::new();
    let mut json = false;
    let mut max_findings = DEFAULT_MAX_FINDINGS;
    let defaults = HighsOptions::default();
    let mut primal_tol = defaults.primal_feasibility_tolerance;
    let mut primal_tol_given = false;
    let mut mip_tol = defaults.mip_feasibility_tolerance;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--format" => match it.next().map(String::as_str) {
                Some("json") => json = true,
                Some("text") => json = false,
                other => return usage(wj, &format!("--format wants text or json, got {other:?}")),
            },
            "--max-findings" => match it.next().and_then(|v| v.parse::<usize>().ok()) {
                Some(n) => max_findings = n,
                None => return usage(wj, "--max-findings wants a non-negative integer"),
            },
            "--primal-feasibility-tolerance" | "--mip-feasibility-tolerance" => {
                match it.next().and_then(|v| v.parse::<f64>().ok()) {
                    Some(x) if x.is_finite() && x >= 0.0 => {
                        if a.starts_with("--primal") {
                            primal_tol = x;
                            primal_tol_given = true;
                        } else {
                            mip_tol = x;
                        }
                    }
                    _ => return usage(wj, &format!("{a} wants a non-negative number")),
                }
            }
            flag if flag.starts_with("--") => return usage(wj, &format!("unknown flag {flag}")),
            _ => positional.push(a),
        }
    }
    let [model_path, solution_path] = positional[..] else {
        return usage(wj, "missing or extra argument: want <model> and <solution>");
    };
    let (lp, has_hessian, gc) = match crate::mip_highs::read_mip_full(model_path) {
        Ok(r) => r,
        // Quadratic rows or cones; the LP/MIP reader refused (and said so).
        Err(_) if rusolve_io::conic::conic_rejected() => {
            return crate::check_conic::run(
                model_path,
                solution_path,
                json,
                max_findings,
                // the conic checker's own default unless the user gave one
                primal_tol_given.then_some(primal_tol),
            );
        }
        Err(e) => {
            return crate::diag::fail_read(
                json,
                model_path,
                e.strip_prefix(&format!("read {model_path}: "))
                    .unwrap_or(&e),
            );
        }
    };
    let text = match std::fs::read_to_string(solution_path) {
        Ok(t) => t,
        Err(e) => return input_error(json, &format!("read {solution_path}: {e}")),
    };
    let is_json = solution_path.ends_with(".json");
    let result = parse_solution(&text, is_json)
        .and_then(|pairs| columns_by_name(&lp, &pairs))
        .and_then(|values| {
            let tol = Tolerances {
                primal_feasibility_tolerance: primal_tol,
                mip_feasibility_tolerance: mip_tol,
                applied: if lp.is_mip() { mip_tol } else { primal_tol },
            };
            check_solution_gc(
                &lp,
                Some(&gc),
                &values,
                tol,
                max_findings,
                model_path,
                solution_path,
            )
        });
    let mut checked = match result {
        Ok(report) => report,
        Err(InputError(e)) => return parse_error(json, &format!("{solution_path}: {e}")),
    };
    if has_hessian && let Err(e) = add_qp(&mut checked, &lp, model_path, &text, is_json) {
        return parse_error(json, &e);
    }
    if json {
        match serde_json::to_string(&checked) {
            Ok(s) => println!("{s}"),
            Err(e) => {
                return crate::diag::fail_msg(
                    json,
                    "RS0090",
                    4,
                    format!("rusolve check: cannot serialise: {e}"),
                );
            }
        }
    } else {
        print!("{}", render_text(&checked));
    }
    ExitCode::SUCCESS
}
