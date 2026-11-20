//! HiGHS source (`read_basis_text` only): lp_data/HighsLpUtils.cpp `readBasisStream`
//! (lines 2938-3055) and Highs.cpp `Highs::readBasis` (lines 952-975).
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Warm-start inputs for `rusolve mip`: `--start`, `--basis`, `--add-rows`.
//!
//! * `--start <file>`: a MIP start in the HiGHS raw solution format or the JSON that
//!   `--solution-file <x>.json` writes, matched to the model by column NAME. A start that names
//!   every column is checked with [`crate::check::check_solution`] and then installed with the
//!   ported `Highs::set_solution` (Highs.cpp:2492); a start that names only some columns goes
//!   through the ported `Highs::set_solution_sparse` (Highs.cpp:2576), and the solver completes
//!   it (`complete_solution_from_discrete_assignment`, highs_run.rs). The solver then logs and
//!   adopts it in `HighsMipSolverData::checkAddSolution` (HighsMipSolverData.cpp:837), which
//!   installs a feasible start as the first incumbent.
//! * `--basis <file>`: a HiGHS basis file (`HiGHS_basis_file v1|v2`), read as HiGHS 1.16.0
//!   `readBasisStream` (lp_data/HighsLpUtils.cpp:2938-3055) and installed as `Highs::readBasis`
//!   does (Highs.cpp:952-975: `valid = useful = true`, `alien` kept) through the ported
//!   `Highs::set_basis`. Only an LP solve uses it: the ported `call_solve_mip` starts from an
//!   invalid basis (highs_run.rs:1143), so for a model with integer columns the file is read and
//!   checked but not installed, and the report says so.
//! * `--add-rows <file>`: rows added to the model before the solve, from an MPS fragment
//!   (ROWS / COLUMNS / RHS / RANGES, free format, columns by name; the
//!   EvoCut approach: a cut family is just ordinary rows) or a JSON row list. Every row passes the validity
//!   gate (docs/design/hooks.md): malformed rows are always rejected, and when
//!   `mip_debug_solution_file` is set a row that cuts off the debug solution is rejected as
//!   `mip/hooks.rs` rejects a user cut (the same feasibility tolerance).
//!
//! Errors are input errors (exit 3) naming the offending column or row and the closest names
//! the model has. Nothing here runs unless one of the three flags is given.

use crate::check::{Tolerances, check_solution, columns_by_name, parse_solution};
use rusolve_ekk::prelude::*;
use serde_json::{Value, json};
use std::collections::{HashMap, HashSet};

// ---------------------------------------------------------------- names

fn col_name(lp: &HighsLp, ix: usize) -> String {
    match lp.col_names_.get(ix) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("c{ix}"),
    }
}

fn row_name(lp: &HighsLp, ix: usize) -> String {
    match lp.row_names_.get(ix) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("r{ix}"),
    }
}

fn levenshtein(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut cur = vec![i + 1];
        for (j, cb) in b.iter().enumerate() {
            let sub = prev[j] + usize::from(ca != *cb);
            cur.push(sub.min(prev[j + 1] + 1).min(cur[j] + 1));
        }
        prev = cur;
    }
    prev[b.len()]
}

/// ` (closest model names: a, b)` for the names within edit distance 3, or the empty string.
fn did_you_mean(name: &str, names: &[String]) -> String {
    let lower = name.to_lowercase();
    let mut scored: Vec<(usize, &String)> = names
        .iter()
        .map(|n| (levenshtein(&lower, &n.to_lowercase()), n))
        .filter(|(d, _)| *d <= 3)
        .collect();
    scored.sort();
    scored.truncate(3);
    if scored.is_empty() {
        String::new()
    } else {
        let list: Vec<&str> = scored.iter().map(|(_, n)| n.as_str()).collect();
        format!(" (closest model names: {})", list.join(", "))
    }
}

fn all_col_names(lp: &HighsLp) -> Vec<String> {
    (0..lp.num_col_ as usize).map(|i| col_name(lp, i)).collect()
}

fn read_text(path: &str, what: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("cannot read the {what} file {path}: {e}"))
}

// ---------------------------------------------------------------- --start

/// What `--start` did, for the envelope.
pub struct StartReport {
    pub path: String,
    pub format: &'static str,
    pub columns_given: usize,
    pub complete: bool,
    /// Complete starts only: the checker's objective (model sense, offset included).
    pub objective: Option<f64>,
    /// The start, one value per column, for the incumbent comparison (complete starts).
    pub values: Option<Vec<f64>>,
}

/// Reads `path`, matches it to `lp` by name and checks it; returns the report and the
/// `(column, value)` pairs to hand to the ported `setSolution`.
pub fn read_start(
    lp: &HighsLp,
    path: &str,
    tol: Tolerances,
) -> Result<(StartReport, Vec<(usize, f64)>), String> {
    let text = read_text(path, "start")?;
    let is_json = path.ends_with(".json");
    let pairs = parse_solution(&text, is_json).map_err(|e| format!("{path}: {}", e.0))?;
    let names = all_col_names(lp);
    let index: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut given: Vec<(usize, f64)> = Vec::new();
    let mut seen = HashSet::new();
    let mut unknown: Vec<String> = Vec::new();
    for (name, value) in &pairs {
        let Some(&ix) = index.get(name.as_str()) else {
            unknown.push(format!("{name}{}", did_you_mean(name, &names)));
            continue;
        };
        if !seen.insert(ix) {
            return Err(format!("{path}: column {name} appears more than once"));
        }
        if !value.is_finite() {
            return Err(format!("{path}: column {name} has the value {value}"));
        }
        given.push((ix, *value));
    }
    if !unknown.is_empty() {
        let shown: Vec<&str> = unknown.iter().take(5).map(String::as_str).collect();
        return Err(format!(
            "{path}: {} column name(s) are not columns of the model, so the start does not belong \
             to this model: {}{}. The start is matched by column name; write it with \
             `--solution-file` from a solve of the same model",
            unknown.len(),
            shown.join("; "),
            if unknown.len() > shown.len() {
                "; ..."
            } else {
                ""
            }
        ));
    }
    let complete = given.len() == names.len();
    let mut report = StartReport {
        path: path.to_string(),
        format: if is_json { "json" } else { "highs-raw" },
        columns_given: given.len(),
        complete,
        objective: None,
        values: None,
    };
    if complete {
        let values = columns_by_name(lp, &pairs).map_err(|e| format!("{path}: {}", e.0))?;
        let rep = check_solution(lp, &values, tol, 5, "", path).map_err(|e| e.0)?;
        if rep.verdict != "feasible" {
            let first: Vec<String> = rep
                .findings
                .iter()
                .map(|f| match f.kind {
                    "row" => format!(
                        "row {} has activity {} outside [{}, {}]",
                        f.name,
                        f.value,
                        f.lower.map_or("-inf".into(), |v| v.to_string()),
                        f.upper.map_or("inf".into(), |v| v.to_string())
                    ),
                    "integrality" => format!(
                        "integer column {} has the fractional value {}",
                        f.name, f.value
                    ),
                    _ => format!(
                        "column {} = {} violates its bound {}",
                        f.name,
                        f.value,
                        f.violated_bound.unwrap_or(f64::NAN)
                    ),
                })
                .collect();
            return Err(format!(
                "{path}: the start is infeasible for the model ({} violation(s); {}). Fix those \
                 values, or give only the columns you trust (a partial start is completed by the \
                 solver); `rusolve check <model> {path}` lists every violation",
                rep.findings_total,
                first.join("; ")
            ));
        }
        report.objective = Some(rep.objective.value);
        report.values = Some(values);
    } else {
        // A partial start cannot be checked against the rows (HiGHS completes it first), but
        // bounds and integrality of the given values are checked as setSolution checks bounds
        // (Highs.cpp:2576-2625), by name.
        for &(ix, v) in &given {
            let (lo, up) = (lp.col_lower_[ix], lp.col_upper_[ix]);
            if v < lo - tol.primal_feasibility_tolerance
                || v > up + tol.primal_feasibility_tolerance
            {
                return Err(format!(
                    "{path}: column {} = {v} is outside its bounds [{lo}, {up}]",
                    names[ix]
                ));
            }
            let is_int = lp.integrality_.get(ix) == Some(&HighsVarType::kInteger);
            if is_int && (v - v.round()).abs() > tol.mip_feasibility_tolerance {
                return Err(format!(
                    "{path}: integer column {} has the fractional value {v}",
                    names[ix]
                ));
            }
        }
    }
    Ok((report, given))
}

/// The `start` object of the envelope, after the solve.
pub fn start_json(
    rep: &StartReport,
    accepted: bool,
    sense_max: bool,
    final_objective: Option<f64>,
    final_x: &[f64],
) -> Value {
    let better = |a: f64, b: f64| {
        let d = if sense_max { a - b } else { b - a };
        d > 1e-9 * (1.0 + b.abs())
    };
    let is_start = rep.values.as_ref().map(|v| {
        final_objective.is_some()
            && v.len() == final_x.len()
            && v.iter().zip(final_x).all(|(a, b)| (a - b).abs() <= 1e-9)
    });
    let improved = match (rep.objective, final_objective) {
        (Some(s), Some(f)) => Some(better(f, s)),
        _ => None,
    };
    json!({
        "path": rep.path,
        "format": rep.format,
        "columns_given": rep.columns_given,
        "complete": rep.complete,
        "accepted": accepted,
        "objective": rep.objective,
        // A feasible complete start is installed as the first incumbent by checkAddSolution;
        // a partial one is completed by the solver, so the answer is only known from the log.
        "became_incumbent": if rep.complete { Value::Bool(accepted) } else { Value::Null },
        "final_incumbent_is_start": is_start,
        "improved_on_start": improved,
    })
}

// ---------------------------------------------------------------- --basis

/// Port of HiGHS 1.16.0 `readBasisStream` (HighsLpUtils.cpp:2938-3055) for a `HighsBasis`
/// sized to `lp`. C++ `in_file >> x` is whitespace tokenising, as here. The C++ `assert`s on
/// the keywords become errors. Returns the basis and whether the file was the deprecated v1.
pub fn read_basis_text(lp: &HighsLp, text: &str) -> Result<(HighsBasis, bool), String> {
    let mut tok = text.split_whitespace();
    let mut next = |what: &str| {
        tok.next()
            .ok_or_else(|| format!("basis file ends before {what}"))
    };
    let _highs = next("the header")?;
    let version = next("the format version")?;
    let v1 = version == kHighsBasisFileV1;
    let v2 = version == kHighsBasisFileV2;
    if !v1 && !v2 {
        return Err(format!(
            "readBasisFile: cannot read a basis file for HiGHS {version} (expected \
             HiGHS_basis_file {kHighsBasisFileV2})"
        ));
    }
    let mut basis = HighsBasis {
        col_status: vec![HighsBasisStatus::kLower; lp.num_col_ as usize],
        row_status: vec![HighsBasisStatus::kLower; lp.num_row_ as usize],
        ..HighsBasis::default()
    };
    let keyword = next("Valid or None")?;
    if keyword == "None" {
        return Err(
            "the basis file holds no basis (\"None\"), so there is nothing to warm start from"
                .to_string(),
        );
    }
    if keyword != "Valid" {
        return Err(format!(
            "basis file: expected Valid or None, found {keyword}"
        ));
    }
    let status_of = |v: i64, what: &str| -> Result<HighsBasisStatus, String> {
        match v {
            0 => Ok(HighsBasisStatus::kLower),
            1 => Ok(HighsBasisStatus::kBasic),
            2 => Ok(HighsBasisStatus::kUpper),
            3 => Ok(HighsBasisStatus::kZero),
            4 => Ok(HighsBasisStatus::kNonbasic),
            _ => Err(format!("basis file: {what} has the invalid status {v}")),
        }
    };
    for (section, num, is_col) in [
        ("Columns", lp.num_col_ as usize, true),
        ("Rows", lp.num_row_ as usize, false),
    ] {
        let _hash = next("a section header")?;
        let kw = next("a section header")?;
        if kw != section {
            return Err(format!(
                "basis file: expected the {section} section, found {kw}"
            ));
        }
        let declared: usize = next("the section size")?
            .parse()
            .map_err(|_| format!("basis file: bad {section} count"))?;
        if declared != num {
            return Err(format!(
                "readBasisFile: basis file is for {declared} {}, not {num}; it does not belong \
                 to this model (rows added with --add-rows count)",
                section.to_lowercase()
            ));
        }
        let names: Vec<String> = (0..num)
            .map(|i| {
                if is_col {
                    col_name(lp, i)
                } else {
                    row_name(lp, i)
                }
            })
            .collect();
        let have_names = if is_col {
            lp.col_names_.len() == num
        } else {
            lp.row_names_.len() == num
        };
        let index: HashMap<&str, usize> = names
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        for x in 0..num {
            let (ix, raw) = if v1 {
                (x, next("a status")?.to_string())
            } else {
                let name = next("a name")?;
                let ix = if have_names {
                    *index.get(name).ok_or_else(|| {
                        format!(
                            "readBasisStream: {} name \"{name}\" is not found{}",
                            if is_col { "column" } else { "row" },
                            did_you_mean(name, &names)
                        )
                    })?
                } else {
                    x
                };
                (ix, next("a status")?.to_string())
            };
            let v: i64 = raw
                .parse()
                .map_err(|_| format!("basis file: \"{raw}\" is not a status number"))?;
            let st = status_of(v, &names[ix])?;
            if is_col {
                basis.col_status[ix] = st;
            } else {
                basis.row_status[ix] = st;
            }
        }
    }
    // Highs::readBasis (Highs.cpp:962-969): valid and useful; `alien` stays true, so the solver
    // forms and factors the basis itself and copes with a singular one.
    basis.valid = true;
    basis.useful = true;
    Ok((basis, v1))
}

// ---------------------------------------------------------------- --add-rows

#[derive(Debug, Clone)]
pub struct RowSpec {
    pub name: String,
    pub lower: f64,
    pub upper: f64,
    pub entries: Vec<(String, f64)>,
}

fn num(v: &Value, what: &str) -> Result<f64, String> {
    match v {
        Value::Number(n) => n.as_f64().ok_or_else(|| format!("{what} is not a number")),
        Value::String(s) => match s.to_lowercase().as_str() {
            "inf" | "+inf" | "infinity" => Ok(f64::INFINITY),
            "-inf" | "-infinity" => Ok(f64::NEG_INFINITY),
            _ => s
                .parse()
                .map_err(|_| format!("{what} \"{s}\" is not a number")),
        },
        Value::Null => Err(format!("{what} is null")),
        _ => Err(format!("{what} is not a number")),
    }
}

/// JSON: `{"rows": [row, ...]}` or `[row, ...]`; a row is
/// `{"name": str, "lower": num?, "upper": num?, "coefficients": {col: num} | [{"name": col,
/// "value": num}]}`. A missing bound is infinite; a `"rhs"` with `"sense"` (`<=`, `>=`, `=`)
/// is accepted as the other spelling of the two bounds.
fn parse_rows_json(text: &str) -> Result<Vec<RowSpec>, String> {
    let v: Value = serde_json::from_str(text).map_err(|e| format!("not valid JSON: {e}"))?;
    let rows = v
        .get("rows")
        .unwrap_or(&v)
        .as_array()
        .ok_or("expected a JSON array of rows or {\"rows\": [...]}")?;
    let mut out = Vec::new();
    for (i, r) in rows.iter().enumerate() {
        let name = r
            .get("name")
            .and_then(Value::as_str)
            .map_or_else(|| format!("added_row{i}"), str::to_string);
        let ctx = format!("row {name}");
        let mut lower = f64::NEG_INFINITY;
        let mut upper = f64::INFINITY;
        if let Some(sense) = r.get("sense").and_then(Value::as_str) {
            let rhs = num(
                r.get("rhs").ok_or(format!("{ctx}: sense without rhs"))?,
                "rhs",
            )?;
            match sense {
                "<=" | "L" | "l" => upper = rhs,
                ">=" | "G" | "g" => lower = rhs,
                "=" | "==" | "E" | "e" => {
                    lower = rhs;
                    upper = rhs;
                }
                _ => return Err(format!("{ctx}: sense \"{sense}\" is not <=, >= or =")),
            }
        }
        if let Some(x) = r.get("lower") {
            lower = num(x, &format!("{ctx}: lower"))?;
        }
        if let Some(x) = r.get("upper") {
            upper = num(x, &format!("{ctx}: upper"))?;
        }
        let mut entries = Vec::new();
        match r.get("coefficients") {
            Some(Value::Object(m)) => {
                for (k, x) in m {
                    entries.push((k.clone(), num(x, &format!("{ctx}: coefficient of {k}"))?));
                }
            }
            Some(Value::Array(a)) => {
                for e in a {
                    let k = e
                        .get("name")
                        .and_then(Value::as_str)
                        .ok_or(format!("{ctx}: coefficient entry has no \"name\""))?;
                    let x = e
                        .get("value")
                        .ok_or(format!("{ctx}: coefficient {k} has no \"value\""))?;
                    entries.push((
                        k.to_string(),
                        num(x, &format!("{ctx}: coefficient of {k}"))?,
                    ));
                }
            }
            _ => return Err(format!("{ctx}: no \"coefficients\" object or array")),
        }
        out.push(RowSpec {
            name,
            lower,
            upper,
            entries,
        });
    }
    Ok(out)
}

/// Free-format MPS fragment: ROWS, COLUMNS, RHS and RANGES (HMpsFF semantics for those
/// sections, HMpsFF.cpp; the objective `N` row is skipped). BOUNDS and the rest are errors: a
/// fragment adds rows, it does not touch columns.
fn parse_rows_mps(text: &str) -> Result<Vec<RowSpec>, String> {
    let mut section = String::new();
    let mut order: Vec<String> = Vec::new();
    let mut kind: HashMap<String, char> = HashMap::new();
    let mut entries: HashMap<String, Vec<(String, f64)>> = HashMap::new();
    let mut rhs: HashMap<String, f64> = HashMap::new();
    let mut range: HashMap<String, f64> = HashMap::new();
    for (ln, line) in text.lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('*') {
            continue;
        }
        let toks: Vec<&str> = t.split_whitespace().collect();
        if !line.starts_with(char::is_whitespace) {
            let head = toks[0].to_uppercase();
            match head.as_str() {
                "NAME" => {}
                "ROWS" | "COLUMNS" | "RHS" | "RANGES" => section = head,
                "ENDATA" => break,
                _ => {
                    return Err(format!(
                        "line {}: section {head} is not allowed in a row fragment (use ROWS, \
                         COLUMNS, RHS, RANGES)",
                        ln + 1
                    ));
                }
            }
            continue;
        }
        let num_at = |s: &str| -> Result<f64, String> {
            s.parse()
                .map_err(|_| format!("line {}: \"{s}\" is not a number", ln + 1))
        };
        match section.as_str() {
            "ROWS" => {
                let ty = toks[0].to_uppercase();
                let (Some(name), 1) = (toks.get(1), ty.len()) else {
                    return Err(format!("line {}: bad ROWS line", ln + 1));
                };
                let ch = ty.chars().next().unwrap();
                if !"NLGE".contains(ch) {
                    return Err(format!(
                        "line {}: row type {ty} is not N, L, G or E",
                        ln + 1
                    ));
                }
                kind.insert(name.to_string(), ch);
                if ch != 'N' {
                    order.push(name.to_string());
                    entries.insert(name.to_string(), Vec::new());
                }
            }
            "COLUMNS" => {
                if toks.iter().any(|t| t.contains("MARKER")) {
                    return Err(format!(
                        "line {}: integer markers are not allowed in a row fragment",
                        ln + 1
                    ));
                }
                if toks.len() < 3 || toks.len().is_multiple_of(2) {
                    return Err(format!("line {}: bad COLUMNS line", ln + 1));
                }
                for p in toks[1..].chunks(2) {
                    let v = num_at(p[1])?;
                    match kind.get(p[0]) {
                        Some('N') => {}
                        Some(_) => entries
                            .get_mut(p[0])
                            .unwrap()
                            .push((toks[0].to_string(), v)),
                        None => {
                            return Err(format!(
                                "line {}: row {} is not declared in ROWS",
                                ln + 1,
                                p[0]
                            ));
                        }
                    }
                }
            }
            "RHS" | "RANGES" => {
                // optional set name: an odd token count has it first (HMpsFF.cpp parseRhs)
                let body = if !toks.len().is_multiple_of(2) {
                    &toks[1..]
                } else {
                    &toks[..]
                };
                for p in body.chunks(2) {
                    if p.len() != 2 {
                        return Err(format!("line {}: bad {section} line", ln + 1));
                    }
                    let v = num_at(p[1])?;
                    match kind.get(p[0]) {
                        Some('N') => {}
                        Some(_) => {
                            if section == "RHS" {
                                rhs.insert(p[0].to_string(), v);
                            } else {
                                range.insert(p[0].to_string(), v);
                            }
                        }
                        None => {
                            return Err(format!(
                                "line {}: row {} is not declared in ROWS",
                                ln + 1,
                                p[0]
                            ));
                        }
                    }
                }
            }
            _ => return Err(format!("line {}: data before a section header", ln + 1)),
        }
    }
    let mut out = Vec::new();
    for name in order {
        let r = rhs.get(&name).copied().unwrap_or(0.0);
        let (mut lower, mut upper) = match kind[&name] {
            'L' => (f64::NEG_INFINITY, r),
            'G' => (r, f64::INFINITY),
            _ => (r, r),
        };
        if let Some(&rg) = range.get(&name) {
            // HMpsFF RANGES rule: E uses the sign of R, L and G the absolute value
            match kind[&name] {
                'L' => lower = r - rg.abs(),
                'G' => upper = r + rg.abs(),
                _ if rg < 0.0 => lower = r + rg,
                _ => upper = r + rg,
            }
        }
        out.push(RowSpec {
            entries: entries.remove(&name).unwrap_or_default(),
            name,
            lower,
            upper,
        });
    }
    Ok(out)
}

/// Rows accepted by the gate, in model column indices.
pub struct AcceptedRow {
    pub name: String,
    pub lower: f64,
    pub upper: f64,
    pub index: Vec<HighsInt>,
    pub value: Vec<f64>,
}

pub struct AddRowsReport {
    pub path: String,
    pub format: &'static str,
    pub requested: usize,
    pub debug_gate: bool,
    pub accepted: Vec<AcceptedRow>,
    pub rejected: Vec<(String, String)>,
}

/// Parses `path` and runs the validity gate on each row. `debug_solution` is the debug
/// solution in original column order when `mip_debug_solution_file` is set (the same
/// condition as `mip/hooks.rs`: activity beyond a bound by more than `feastol` cuts it off).
pub fn read_and_gate_rows(
    lp: &HighsLp,
    path: &str,
    debug_solution: Option<&[f64]>,
    feastol: f64,
) -> Result<AddRowsReport, String> {
    let text = read_text(path, "add-rows")?;
    let is_json = path.ends_with(".json");
    let specs = if is_json {
        parse_rows_json(&text)
    } else {
        parse_rows_mps(&text)
    }
    .map_err(|e| format!("{path}: {e}"))?;
    let names = all_col_names(lp);
    let index: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut taken: HashSet<String> = (0..lp.num_row_ as usize).map(|i| row_name(lp, i)).collect();
    let mut report = AddRowsReport {
        path: path.to_string(),
        format: if is_json { "json" } else { "mps-fragment" },
        requested: specs.len(),
        debug_gate: debug_solution.is_some(),
        accepted: Vec::new(),
        rejected: Vec::new(),
    };
    'row: for spec in specs {
        let mut reject = |why: String| report.rejected.push((spec.name.clone(), why));
        if spec.name.is_empty() || spec.name.contains(char::is_whitespace) {
            reject("malformed: the row name is empty or contains whitespace".into());
            continue;
        }
        if !taken.insert(spec.name.clone()) {
            reject(
                "malformed: the name is already used by a row of the model or an earlier added row"
                    .into(),
            );
            continue;
        }
        if spec.lower.is_nan() || spec.upper.is_nan() || spec.lower > spec.upper {
            reject(format!(
                "malformed: bounds [{}, {}] are not an interval",
                spec.lower, spec.upper
            ));
            continue;
        }
        if spec.entries.is_empty() {
            reject("malformed: the row has no coefficients".into());
            continue;
        }
        let mut idx = Vec::new();
        let mut val = Vec::new();
        let mut seen = HashSet::new();
        for (cname, v) in &spec.entries {
            let Some(&ix) = index.get(cname.as_str()) else {
                reject(format!(
                    "malformed: column {cname} is not a column of the model{}",
                    did_you_mean(cname, &names)
                ));
                continue 'row;
            };
            if !v.is_finite() {
                reject(format!("malformed: the coefficient of {cname} is {v}"));
                continue 'row;
            }
            if !seen.insert(ix) {
                reject(format!("malformed: column {cname} appears twice"));
                continue 'row;
            }
            idx.push(ix as HighsInt);
            val.push(*v);
        }
        if let Some(dbg) = debug_solution {
            let act: f64 = idx
                .iter()
                .zip(&val)
                .map(|(&col, v)| dbg[col as usize] * v)
                .sum();
            if act < spec.lower - feastol || act > spec.upper + feastol {
                reject(format!(
                    "cuts off the debug solution (activity {act} outside [{}, {}])",
                    spec.lower, spec.upper
                ));
                continue;
            }
        }
        report.accepted.push(AcceptedRow {
            name: spec.name,
            lower: spec.lower,
            upper: spec.upper,
            index: idx,
            value: val,
        });
    }
    Ok(report)
}

/// Adds the accepted rows to the ported `Highs` model and names them. Returns an error text
/// when `Highs::add_rows` refuses.
pub fn add_rows(
    highs: &mut rusolve_ekk::lp_data::highs_run::Highs,
    rows: &[AcceptedRow],
) -> Result<(), String> {
    if rows.is_empty() {
        return Ok(());
    }
    let (mut starts, mut index, mut value) = (Vec::new(), Vec::new(), Vec::new());
    for r in rows {
        starts.push(index.len() as HighsInt);
        index.extend_from_slice(&r.index);
        value.extend_from_slice(&r.value);
    }
    let lower: Vec<f64> = rows.iter().map(|r| r.lower).collect();
    let upper: Vec<f64> = rows.iter().map(|r| r.upper).collect();
    let old = highs.model_lp_.num_row_ as usize;
    if highs.add_rows(
        rows.len() as HighsInt,
        &lower,
        &upper,
        index.len() as HighsInt,
        &starts,
        &index,
        &value,
    ) == HighsStatus::kError
    {
        return Err("the ported Highs::add_rows refused the rows".into());
    }
    let lp = &mut highs.model_lp_;
    lp.row_names_.resize(old, String::new());
    lp.row_names_.extend(rows.iter().map(|r| r.name.clone()));
    lp.row_names_.resize(lp.num_row_ as usize, String::new());
    Ok(())
}

pub fn add_rows_json(rep: &AddRowsReport) -> Value {
    let accepted: Vec<&str> = rep.accepted.iter().map(|r| r.name.as_str()).collect();
    let rejected: Vec<Value> = rep
        .rejected
        .iter()
        .map(|(n, why)| json!({"name": n, "reason": why}))
        .collect();
    json!({
        "path": rep.path,
        "format": rep.format,
        "requested": rep.requested,
        "added": rep.accepted.len(),
        "accepted": rep.accepted.len(),
        "rejected": rep.rejected.len(),
        "debug_solution_gate": rep.debug_gate,
        "accepted_names": accepted,
        "rejected_rows": rejected,
    })
}

/// Reads the debug solution named by `mip_debug_solution_file` into original column order.
pub fn read_debug_solution(lp: &HighsLp, path: &str) -> Result<Vec<f64>, String> {
    let text = read_text(path, "mip_debug_solution_file")?;
    let pairs =
        parse_solution(&text, path.ends_with(".json")).map_err(|e| format!("{path}: {}", e.0))?;
    columns_by_name(lp, &pairs).map_err(|e| format!("{path}: {}", e.0))
}

/// The `basis` object of the envelope.
pub fn basis_json(path: &str, v1: bool, installed: bool, is_mip: bool) -> Value {
    json!({
        "path": path,
        "format": if v1 { "v1" } else { "v2" },
        "read": true,
        "installed": installed,
        "note": if is_mip {
            "the model has integer columns: the port's MIP solve starts from its own root LP, so the basis was read and checked but not used"
        } else {
            "installed through the ported Highs::set_basis as an LP warm start"
        },
    })
}
