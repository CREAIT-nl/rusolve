//! `rusolve ranging <model.mps> [--rows PAT[,PAT]...] [--cols PAT[,PAT]...] [--format json|text]`
//! Sensitivity by name for an LP, and for a MIP on its
//! fixed-integer LP.
//!
//! ## What HiGHS requires
//!
//! `Highs::getRanging` (HiGHS 1.16.0 lp_data/Highs.cpp:2398) calls
//! `getRangingInterface` (HighsInterface.cpp:1856) which calls `getRangingData`
//! (HighsRanging.cpp:70). That function refuses unless the model status is
//! `kOptimal` ("Cannot get ranging without an optimal solution") and unless the
//! `HEkk` instance is initialised for the solve ("Cannot get ranging without a
//! valid Simplex instance"). The simplex instance is the one of the *original*
//! model, so the LP is solved with `presolve = off` and the simplex solver: a
//! presolve-reduced solve would leave `ekk_instance_` describing the reduced
//! LP, and a crossover/IPM solve leaves none. This command sets both options.
//!
//! ## MIPs
//!
//! HiGHS has no dual solution for a MIP: `Highs::run` ends a MIP solve with
//! only primal values (the MIP solver never produces `solution_.dual_valid`),
//! and `getRanging` on it fails the model-status/simplex checks above. The
//! standard practice (the "fixed model" of common MIP solver APIs) is to fix every integer variable at its incumbent value,
//! drop integrality and read the LP's duals and ranging. This command does
//! that: solve the MIP with the ported `Highs::run`, fix the integer columns at
//! the rounded incumbent, re-solve the LP and report its ranging with
//! `"basis_of": "fixed_integer_lp"`. The ranging is therefore valid only for
//! the incumbent's integer assignment, not for the MIP.
//!
//! ## Indices
//!
//! `in_var_`/`ou_var_` of a `HighsRangingRecord` are simplex variable indices:
//! `0..num_col` are columns and `num_col + i` is row `i`
//! (HighsRanging.cpp: variables are columns then rows). They are reported as
//! `{"kind": "col"|"row", "name": ...}`; `null` for -1 (no such variable).
//!
//! Exit codes follow the exit-code table (see `mip_highs`): 0 the solve ran (a non-optimal
//! status is reported, with no ranging), 2 usage, 3 input, 4 internal (the
//! solve or `get_ranging` returned `kError` on an optimal model).

use crate::mip_highs::read_mip;
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_ranging::{HighsRanging, HighsRangingRecord};
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::{Value, json};
use std::process::ExitCode;

pub const SCHEMA_VERSION: &str = "rusolve.ranging/1";

/// JSON has no inf/nan: non-finite is `null` (same convention as
/// `rusolve.result/1`, envelope.rs `finite`).
fn num(v: f64) -> Value {
    if v.is_finite() { json!(v) } else { Value::Null }
}

/// `*` matches any run of characters; everything else is literal.
fn glob(pat: &str, s: &str) -> bool {
    let parts: Vec<&str> = pat.split('*').collect();
    if parts.len() == 1 {
        return pat == s;
    }
    let (first, last) = (parts[0], parts[parts.len() - 1]);
    if !s.starts_with(first) || s.len() < first.len() + last.len() || !s.ends_with(last) {
        return false;
    }
    let mut rest = &s[first.len()..s.len() - last.len()];
    for mid in &parts[1..parts.len() - 1] {
        match rest.find(mid) {
            Some(p) => rest = &rest[p + mid.len()..],
            None => return false,
        }
    }
    true
}

/// The indices whose name matches some pattern (all when `pats` is empty), or
/// the first pattern that matches nothing.
fn select(names: &[String], pats: &[String]) -> Result<Vec<usize>, String> {
    if pats.is_empty() {
        return Ok((0..names.len()).collect());
    }
    let mut keep = vec![false; names.len()];
    for p in pats {
        let mut any = false;
        for (i, n) in names.iter().enumerate() {
            if glob(p, n) {
                keep[i] = true;
                any = true;
            }
        }
        if !any {
            return Err(p.clone());
        }
    }
    Ok((0..names.len()).filter(|&i| keep[i]).collect())
}

fn status_name(s: HighsBasisStatus) -> &'static str {
    match s {
        HighsBasisStatus::kLower => "lower",
        HighsBasisStatus::kBasic => "basic",
        HighsBasisStatus::kUpper => "upper",
        HighsBasisStatus::kZero => "zero",
        HighsBasisStatus::kNonbasic => "nonbasic",
    }
}

struct Names<'a> {
    cols: &'a [String],
    rows: &'a [String],
}

impl Names<'_> {
    fn var(&self, v: HighsInt) -> Value {
        if v < 0 {
            return Value::Null;
        }
        let v = v as usize;
        if v < self.cols.len() {
            json!({"kind": "col", "name": self.cols[v]})
        } else if v - self.cols.len() < self.rows.len() {
            json!({"kind": "row", "name": self.rows[v - self.cols.len()]})
        } else {
            Value::Null
        }
    }
}

/// One direction of one ranging: the limit of the parameter, the objective
/// there, and the variables entering and leaving the basis.
fn side(r: &HighsRangingRecord, i: usize, n: &Names) -> Value {
    json!({
        "value": num(r.value_[i]),
        "objective": num(r.objective_[i]),
        "in_var": n.var(r.in_var_[i]),
        "out_var": n.var(r.ou_var_[i]),
    })
}

/// Fixes the integer columns of `lp` at their rounded incumbent values and
/// drops integrality: the fixed-integer LP.
fn fixed_integer_lp(lp: &HighsLp, x: &[f64]) -> HighsLp {
    let mut fixed = lp.clone();
    for (j, ty) in lp.integrality_.iter().enumerate() {
        if *ty == HighsVarType::kInteger {
            let v = x[j].round().clamp(lp.col_lower_[j], lp.col_upper_[j]);
            fixed.col_lower_[j] = v;
            fixed.col_upper_[j] = v;
        }
    }
    fixed.integrality_.clear();
    fixed
}

fn usage(json: bool, msg: &str) -> ExitCode {
    crate::diag::fail_usage_text(json, "ranging", msg)
}

pub fn run_ranging(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let Some(path) = args.first() else {
        return usage(wj, "missing MPS path");
    };
    let mut json_out = false;
    let mut row_pats: Vec<String> = Vec::new();
    let mut col_pats: Vec<String> = Vec::new();
    let mut it = args[1..].iter();
    while let Some(flag) = it.next() {
        match flag.as_str() {
            "--format" => match it.next().map(String::as_str) {
                Some("json") => json_out = true,
                Some("text") => json_out = false,
                _ => return usage(wj, "--format needs json or text"),
            },
            "--rows" | "--cols" => match it.next() {
                Some(v) => {
                    let dst = if flag == "--rows" {
                        &mut row_pats
                    } else {
                        &mut col_pats
                    };
                    dst.extend(v.split(',').filter(|s| !s.is_empty()).map(String::from));
                }
                None => return usage(wj, &format!("{flag} needs a name or pattern list")),
            },
            other if !other.starts_with("--") => {
                return usage(wj, &format!("unexpected extra argument '{other}'"));
            }
            other => return usage(wj, &format!("unsupported flag '{other}'")),
        }
    }
    // Sensitivity of a semi-continuous model has no fixed-integer LP to read it from.
    let lp = match read_mip(path) {
        Ok(lp) => lp,
        Err(e) => return crate::diag::fail_read_str(json_out, path, &e),
    };
    if lp.integrality_.iter().any(|k| {
        matches!(
            k,
            rusolve_ekk::prelude::HighsVarType::kSemiContinuous
                | rusolve_ekk::prelude::HighsVarType::kSemiInteger
        )
    }) {
        return crate::diag::fail_msg(
            json_out,
            "RS0012",
            3,
            "rusolve ranging: the model has semi-continuous columns; ranging does not support them",
        );
    }
    match ranging_core(&lp, path, &row_pats, &col_pats) {
        Ok((out, why)) => {
            if json_out {
                println!("{out}");
            } else if out["ranging_valid"] == json!(true) {
                print_text(&out);
            } else {
                println!(
                    "status {}: {}",
                    out["status"].as_str().unwrap_or(""),
                    why.unwrap_or("")
                );
            }
            ExitCode::SUCCESS
        }
        Err(e) => crate::diag::fail_core(json_out, "ranging", &e),
    }
}

/// The body of `rusolve ranging` on a loaded model: the `rusolve.ranging/1`
/// object, and for a model with nothing to range (no incumbent, not optimal)
/// the reason the text report gives. Shared by the CLI and the MCP `ranging`
/// tool. `path` is only the name reported in `input.path`.
pub(crate) fn ranging_core(
    lp: &HighsLp,
    path: &str,
    row_pats: &[String],
    col_pats: &[String],
) -> Result<(Value, Option<&'static str>), crate::cli_core::CoreError> {
    use crate::cli_core::CoreError;
    let is_mip = lp.integrality_.contains(&HighsVarType::kInteger);
    let integer_cols = lp
        .integrality_
        .iter()
        .filter(|t| **t == HighsVarType::kInteger)
        .count();
    let (col_names, row_names) = (lp.col_names_.clone(), lp.row_names_.clone());
    let cols_sel = match select(&col_names, col_pats) {
        Ok(v) => v,
        Err(p) => return Err(CoreError::usage(format!("--cols '{p}' matches no column"))),
    };
    let rows_sel = match select(&row_names, row_pats) {
        Ok(v) => v,
        Err(p) => return Err(CoreError::usage(format!("--rows '{p}' matches no row"))),
    };

    // Step 1 (MIP only): solve the MIP, then fix its integers.
    let (lp_to_solve, mip_status) = if is_mip {
        let mut h = Highs::with_options(crate::base_options());
        h.set_option_value_bool("output_flag", false);
        if h.pass_model(lp.clone()) == HighsStatus::kError || h.run() == HighsStatus::kError {
            return Err(CoreError::internal("the MIP solve failed"));
        }
        let st = h.get_model_status();
        if !h.solution_.value_valid {
            return Ok(report_unsolved(
                path,
                lp,
                integer_cols,
                "mip",
                st,
                "the MIP has no incumbent to fix",
            ));
        }
        (fixed_integer_lp(lp, &h.solution_.col_value), Some(st))
    } else {
        (lp.clone(), None)
    };
    let basis_of = if is_mip { "fixed_integer_lp" } else { "lp" };

    // Step 2: solve the LP with the simplex on the original model (see module doc).
    let mut h = Highs::with_options(HighsOptions::default());
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_str("presolve", "off");
    h.set_option_value_str("solver", "simplex");
    if h.pass_model(lp_to_solve) == HighsStatus::kError || h.run() == HighsStatus::kError {
        return Err(CoreError::internal("the LP solve failed"));
    }
    let status = h.get_model_status();
    if status != HighsModelStatus::kOptimal {
        return Ok(report_unsolved(
            path,
            lp,
            integer_cols,
            basis_of,
            status,
            "no ranging without an optimal LP",
        ));
    }
    let mut ranging = HighsRanging::default();
    if h.get_ranging(&mut ranging) == HighsStatus::kError || !ranging.valid {
        return Err(CoreError::internal(
            "Highs::get_ranging failed on an optimal LP",
        ));
    }

    let sol = &h.solution_;
    let basis = &h.basis_;
    let m = &h.model_lp_;
    let n = Names {
        cols: &col_names,
        rows: &row_names,
    };
    let columns: Vec<Value> = cols_sel
        .iter()
        .map(|&j| {
            json!({
                "name": col_names[j],
                "value": num(sol.col_value[j]),
                "reduced_cost": num(sol.col_dual[j]),
                "basis_status": status_name(basis.col_status[j]),
                "cost": num(m.col_cost_[j]),
                "lower": num(m.col_lower_[j]),
                "upper": num(m.col_upper_[j]),
                "cost_ranging": {
                    "down": side(&ranging.col_cost_dn, j, &n),
                    "up": side(&ranging.col_cost_up, j, &n),
                },
                "bound_ranging": {
                    "down": side(&ranging.col_bound_dn, j, &n),
                    "up": side(&ranging.col_bound_up, j, &n),
                },
            })
        })
        .collect();
    let rows: Vec<Value> = rows_sel
        .iter()
        .map(|&i| {
            json!({
                "name": row_names[i],
                "activity": num(sol.row_value[i]),
                "dual": num(sol.row_dual[i]),
                "basis_status": status_name(basis.row_status[i]),
                "lower": num(m.row_lower_[i]),
                "upper": num(m.row_upper_[i]),
                "bound_ranging": {
                    "down": side(&ranging.row_bound_dn, i, &n),
                    "up": side(&ranging.row_bound_up, i, &n),
                },
            })
        })
        .collect();
    let out = json!({
        "schema_version": SCHEMA_VERSION,
        "solver": {"name": "rusolve", "version": env!("CARGO_PKG_VERSION"), "engine": "highs-port"},
        "input": input_json(path, lp, integer_cols),
        "basis_of": basis_of,
        "status": crate::envelope::result_status(status),
        "mip_status": mip_status.map(crate::envelope::result_status),
        "objective": num(h.get_info().objective_function_value),
        "ranging_valid": true,
        "columns": columns,
        "rows": rows,
    });
    Ok((out, None))
}

fn input_json(path: &str, lp: &HighsLp, integer_cols: usize) -> Value {
    json!({"path": path, "rows": lp.num_row_, "cols": lp.num_col_, "integer_cols": integer_cols})
}

/// The solve ran but there is nothing to range: report the status, exit 0.
fn report_unsolved(
    path: &str,
    lp: &HighsLp,
    integer_cols: usize,
    basis_of: &str,
    status: HighsModelStatus,
    why: &'static str,
) -> (Value, Option<&'static str>) {
    let out = json!({
        "schema_version": SCHEMA_VERSION,
        "solver": {"name": "rusolve", "version": env!("CARGO_PKG_VERSION"), "engine": "highs-port"},
        "input": input_json(path, lp, integer_cols),
        "basis_of": basis_of,
        "status": crate::envelope::result_status(status),
        "mip_status": Value::Null,
        "objective": Value::Null,
        "ranging_valid": false,
        "columns": [],
        "rows": [],
    });
    (out, Some(why))
}

fn cell(v: &Value) -> String {
    match v {
        Value::Null => "inf".to_string(),
        other => format!("{:.6e}", other.as_f64().unwrap_or(f64::NAN)),
    }
}

fn var_name(v: &Value) -> String {
    match v["name"].as_str() {
        Some(n) => format!("{}:{n}", v["kind"].as_str().unwrap_or("")),
        None => "-".to_string(),
    }
}

fn print_text(out: &Value) {
    println!(
        "status {}  objective {}  basis_of {}",
        out["status"].as_str().unwrap_or(""),
        cell(&out["objective"]),
        out["basis_of"].as_str().unwrap_or("")
    );
    println!(
        "\ncolumns: name status value reduced_cost | cost down/up (in,out) | bound down/up (in,out)"
    );
    for col in out["columns"].as_array().into_iter().flatten() {
        let (cd, cu) = (&col["cost_ranging"]["down"], &col["cost_ranging"]["up"]);
        let (bd, bu) = (&col["bound_ranging"]["down"], &col["bound_ranging"]["up"]);
        println!(
            "  {} {} {} {} | {} {} ({}/{}) | {} {} ({}/{})",
            col["name"].as_str().unwrap_or(""),
            col["basis_status"].as_str().unwrap_or(""),
            cell(&col["value"]),
            cell(&col["reduced_cost"]),
            cell(&cd["value"]),
            cell(&cu["value"]),
            var_name(&cu["in_var"]),
            var_name(&cu["out_var"]),
            cell(&bd["value"]),
            cell(&bu["value"]),
            var_name(&bu["in_var"]),
            var_name(&bu["out_var"]),
        );
    }
    println!("\nrows: name status activity dual | bound down/up (in,out)");
    for row in out["rows"].as_array().into_iter().flatten() {
        let (bd, bu) = (&row["bound_ranging"]["down"], &row["bound_ranging"]["up"]);
        println!(
            "  {} {} {} {} | {} {} ({}/{})",
            row["name"].as_str().unwrap_or(""),
            row["basis_status"].as_str().unwrap_or(""),
            cell(&row["activity"]),
            cell(&row["dual"]),
            cell(&bd["value"]),
            cell(&bu["value"]),
            var_name(&bu["in_var"]),
            var_name(&bu["out_var"]),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::glob;

    #[test]
    fn glob_matches() {
        assert!(glob("X*", "X01"));
        assert!(glob("*01", "X01"));
        assert!(glob("a*b*c", "a--b--c"));
        assert!(!glob("a*b*c", "a--c"));
        assert!(glob("x", "x"));
        assert!(!glob("x", "xy"));
    }
}
