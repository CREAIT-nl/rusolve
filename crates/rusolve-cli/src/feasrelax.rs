//! `rusolve feasrelax <model.mps>`
//! tells an agent the cheapest way to make an infeasible model feasible, by name.
//!
//! ```text
//! rusolve feasrelax <model.mps> [--mode weighted|min-count] [--format text|json]
//!     [--global-lower-penalty W] [--global-upper-penalty W] [--global-rhs-penalty W]
//!     [--penalties FILE.json] [--big-m M] [--time-limit S]
//! ```
//!
//! ## Modes
//!
//! * `weighted` (default): the ported `Highs::feasibility_relaxation`
//!   (`rusolve-ekk/src/lp_data/highs_elasticity.rs`; HiGHS 1.16.0
//!   `Highs.cpp:2409` `Highs::feasibilityRelaxation`, which calls
//!   `HighsInterface.cpp` `Highs::elasticityFilter`). It minimises
//!   `sum(weight * violation)`. The weights are HiGHS's arguments: a negative
//!   one means "do not relax" (`elasticityFilter`: "Negative penalty implies
//!   that the bounds cannot be violated"). `--penalties` supplies per-name
//!   weights, which become HiGHS's `local_*_penalty` arrays and override the
//!   globals.
//! * `min-count`: minimise the number of relaxed constraints (cost 1 each; the
//!   weights only decide *which* sides may be relaxed, negative = may not, and
//!   weight the violations of stage two).
//!   For every relaxable row side and column bound an indicator binary `z`
//!   is added with a big-M, the elastic/indicator formulation of Chinneck,
//!   "Feasibility and Infeasibility in Optimization" (Springer 2008), ch. 6
//!   (elastic programming with 0-1 "deletion" variables), and the
//!   "minimise the number of violated constraints" objective that MIP solver
//!   APIs commonly offer as a feasibility relaxation. `a.x >= L` becomes
//!   `a.x + M z >= L` with `M = L - minact(a.x)`, the row's own activity-bound
//!   gap over the (possibly relaxed) variable box, which is the tightest valid
//!   M; `a.x <= U` becomes `a.x - M z <= U` with `M = maxact(a.x) - U`. A
//!   side whose gap is `<= 0` can never be violated and gets no indicator. A
//!   side whose activity bound is infinite, and every column bound (a relaxed
//!   column is unbounded on that side), uses `--big-m` and is reported in
//!   `notes`. The auxiliary MIP is built through the ported modify API
//!   (`add_col`, `add_row`) on a fresh `Highs` and solved by the port with
//!   `mip_feasibility_tolerance = 1e-9` so that `M z` cannot leak.
//!   A MIP only fixes *which* constraints give way, not by how much, so a second
//!   stage re-solves the weighted relaxation with only those constraints
//!   relaxable (everything else weight -1), which makes the reported violations
//!   minimal for that set. If stage two relaxes more constraints than stage one
//!   (a HiGHS row penalty covers both sides of a row), stage one's point is kept.
//!
//! The report is computed from the relaxed point alone: a side is "relaxed" iff
//! its row activity (or column value) violates the original bound by more than
//! `1e-7` (HiGHS's `primal_feasibility_tolerance`). The edit for it moves the
//! bound exactly onto the activity, and `verified` re-solves the original model
//! with all edits applied and confirms it is feasible.
//!
//! Exit codes: 0 ran (a model that cannot be relaxed is `status:
//! NotRelaxable`, not an error), 2 usage, 3 input, 4 internal.

use crate::base_options;
use crate::mip_highs::read_mip;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::{Value, json};
use std::process::ExitCode;

pub const SCHEMA_VERSION: &str = "rusolve.feasrelax/1";

/// HiGHS's `primal_feasibility_tolerance` default.
const VIOL_TOL: f64 = 1e-7;

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Weighted,
    MinCount,
}

#[derive(Clone)]
pub(crate) struct Args {
    path: String,
    mode: Mode,
    json: bool,
    gl: f64,
    gu: f64,
    gr: f64,
    penalties: Option<String>,
    big_m: f64,
    time_limit: f64,
}

/// Resolved weights, one per row and per column side. Negative: not relaxable.
#[derive(Clone)]
struct Weights {
    row: Vec<f64>,
    col_lo: Vec<f64>,
    col_up: Vec<f64>,
}

/// One relaxed row side or column bound.
struct Item {
    is_row: bool,
    index: usize,
    name: String,
    lower_side: bool,
    original: f64,
    relaxed_to: f64,
    violation: f64,
    weight: f64,
}

enum Outcome {
    /// The relaxed point (one value per original column) and the solver objective.
    Point(Vec<f64>, f64, Vec<String>),
    NotRelaxable(Vec<String>),
    Failed(String),
}

fn usage(json: bool, msg: &str) -> ExitCode {
    crate::diag::fail_usage_text(
        json,
        "feasrelax",
        &format!(
            "{msg}\nusage: rusolve feasrelax <model.mps> [--mode weighted|min-count] [--format text|json] \
             [--global-lower-penalty W] [--global-upper-penalty W] [--global-rhs-penalty W] \
             [--penalties FILE.json] [--big-m M] [--time-limit S]"
        ),
    )
}

fn parse_args(args: &[String]) -> Result<Args, String> {
    let mut it = args.iter();
    let path = it.next().ok_or("missing model path")?.clone();
    let mut a = Args {
        path,
        mode: Mode::Weighted,
        json: false,
        gl: 1.0,
        gu: 1.0,
        gr: 1.0,
        penalties: None,
        big_m: 1e4,
        time_limit: f64::INFINITY,
    };
    fn num(flag: &str, v: Option<&String>) -> Result<f64, String> {
        let v = v.ok_or(format!("{flag} needs a value"))?;
        v.parse::<f64>()
            .ok()
            .filter(|x| !x.is_nan())
            .ok_or(format!("{flag}: '{v}' is not a number"))
    }
    while let Some(flag) = it.next() {
        match flag.as_str() {
            "--mode" => match it.next().map(String::as_str) {
                Some("weighted") => a.mode = Mode::Weighted,
                Some("min-count") => a.mode = Mode::MinCount,
                _ => return Err("--mode wants weighted or min-count".into()),
            },
            "--format" => match it.next().map(String::as_str) {
                Some("json") => a.json = true,
                Some("text") => a.json = false,
                _ => return Err("--format wants text or json".into()),
            },
            "--global-lower-penalty" => a.gl = num(flag, it.next())?,
            "--global-upper-penalty" => a.gu = num(flag, it.next())?,
            "--global-rhs-penalty" => a.gr = num(flag, it.next())?,
            "--penalties" => {
                a.penalties = Some(it.next().ok_or("--penalties needs a file")?.clone());
            }
            "--big-m" => {
                a.big_m = num(flag, it.next())?;
                if !(a.big_m > 0.0 && a.big_m.is_finite()) {
                    return Err("--big-m must be positive and finite".into());
                }
            }
            "--time-limit" => a.time_limit = num(flag, it.next())?,
            other if !other.starts_with("--") => {
                return Err(format!("unexpected extra argument '{other}'"));
            }
            other => return Err(format!("unsupported flag '{other}'")),
        }
    }
    Ok(a)
}

/// `{ "rows": {"name": w}, "cols": {"name": {"lower": w, "upper": w}} }` over the
/// global defaults. `Err` is a usage-level problem (unknown name, wrong type).
fn resolve_weights(lp: &HighsLp, a: &Args, file: Option<&Value>) -> Result<Weights, String> {
    let n = lp.num_col_ as usize;
    let m = lp.num_row_ as usize;
    let mut w = Weights {
        row: vec![a.gr; m],
        col_lo: vec![a.gl; n],
        col_up: vec![a.gu; n],
    };
    let Some(doc) = file else { return Ok(w) };
    let find = |names: &[String], what: &str, name: &str| -> Result<usize, String> {
        names
            .iter()
            .position(|x| x == name)
            .ok_or(format!("penalties: no {what} named '{name}'"))
    };
    let weight = |v: &Value, ctx: &str| -> Result<f64, String> {
        v.as_f64()
            .ok_or(format!("penalties: {ctx} must be a number"))
    };
    if let Some(rows) = doc.get("rows") {
        let rows = rows
            .as_object()
            .ok_or("penalties: 'rows' must be an object")?;
        for (name, v) in rows {
            let i = find(&lp.row_names_, "row", name)?;
            w.row[i] = weight(v, &format!("rows.{name}"))?;
        }
    }
    if let Some(cols) = doc.get("cols") {
        let cols = cols
            .as_object()
            .ok_or("penalties: 'cols' must be an object")?;
        for (name, v) in cols {
            let j = find(&lp.col_names_, "column", name)?;
            let o = v
                .as_object()
                .ok_or(format!("penalties: cols.{name} must be an object"))?;
            if let Some(x) = o.get("lower") {
                w.col_lo[j] = weight(x, &format!("cols.{name}.lower"))?;
            }
            if let Some(x) = o.get("upper") {
                w.col_up[j] = weight(x, &format!("cols.{name}.upper"))?;
            }
        }
    }
    Ok(w)
}

/// Row-wise copy of the constraint matrix: `rows[i] = [(col, value)]`.
fn row_entries(lp: &HighsLp) -> Vec<Vec<(usize, f64)>> {
    let m = lp.num_row_ as usize;
    let n = lp.num_col_ as usize;
    let a = &lp.a_matrix_;
    let mut rows = vec![Vec::new(); m];
    if a.format_ == MatrixFormat::kRowwise {
        for (i, row) in rows.iter_mut().enumerate() {
            for k in a.start_[i] as usize..a.start_[i + 1] as usize {
                row.push((a.index_[k] as usize, a.value_[k]));
            }
        }
    } else {
        for j in 0..n {
            for k in a.start_[j] as usize..a.start_[j + 1] as usize {
                rows[a.index_[k] as usize].push((j, a.value_[k]));
            }
        }
    }
    rows
}

fn solver(a: &Args) -> Highs {
    let mut opts = base_options();
    opts.output_flag = false;
    opts.log_to_console = false;
    opts.threads = 1;
    opts.time_limit = a.time_limit;
    opts.sync_log_options();
    Highs::with_options(opts)
}

/// Stage "weighted": HiGHS `Highs::feasibilityRelaxation` on a copy of the model.
fn solve_weighted(lp: &HighsLp, w: &Weights, a: &Args) -> Outcome {
    let n = lp.num_col_ as usize;
    let mut highs = solver(a);
    if highs.pass_model(lp.clone()) == HighsStatus::kError {
        return Outcome::Failed("pass_model failed".into());
    }
    // Local arrays carry every weight, so the globals are only a fallback that is
    // never read; HiGHS requires at least one relaxable kind to be announced,
    // which an array always does (`has_local_*_penalty`).
    let st = highs.feasibility_relaxation(
        a.gl,
        a.gu,
        a.gr,
        Some(&w.col_lo),
        Some(&w.col_up),
        Some(&w.row),
    );
    let sol = highs.get_solution();
    if st == HighsStatus::kError {
        return Outcome::Failed("feasibility_relaxation returned an error".into());
    }
    if !sol.value_valid || sol.col_value.len() != n {
        return Outcome::NotRelaxable(vec![]);
    }
    let x = sol.col_value.clone();
    // An elastic LP that is infeasible (the non-relaxable bounds conflict) leaves
    // a point that violates a bound HiGHS was told not to touch.
    if hard_violated(lp, &x, w) {
        return Outcome::NotRelaxable(vec![]);
    }
    Outcome::Point(x, highs.get_info().objective_function_value, vec![])
}

/// True iff `x` violates a bound whose weight is negative.
fn hard_violated(lp: &HighsLp, x: &[f64], w: &Weights) -> bool {
    items_of(lp, x, &|_, _, _| 0.0)
        .iter()
        .any(|it| weight_of(w, it.is_row, it.index, it.lower_side) < 0.0)
}

fn weight_of(w: &Weights, is_row: bool, i: usize, lower: bool) -> f64 {
    if is_row {
        w.row[i]
    } else if lower {
        w.col_lo[i]
    } else {
        w.col_up[i]
    }
}

/// Every row side and column bound that `x` violates by more than [`VIOL_TOL`].
/// `wt(is_row, index, lower_side)` supplies the reported weight.
#[allow(clippy::needless_range_loop)]
fn items_of(lp: &HighsLp, x: &[f64], wt: &dyn Fn(bool, usize, bool) -> f64) -> Vec<Item> {
    let rows = row_entries(lp);
    let mut out = Vec::new();
    let mut push = |is_row: bool, i: usize, name: &str, lower_side: bool, bound: f64, v: f64| {
        let violation = if lower_side { bound - v } else { v - bound };
        if violation > VIOL_TOL {
            out.push(Item {
                is_row,
                index: i,
                name: name.to_string(),
                lower_side,
                original: bound,
                relaxed_to: v,
                violation,
                weight: wt(is_row, i, lower_side),
            });
        }
    };
    for (i, row) in rows.iter().enumerate() {
        let act: f64 = row.iter().map(|&(j, v)| v * x[j]).sum();
        let name = lp.row_names_.get(i).map(String::as_str).unwrap_or("");
        if lp.row_lower_[i] > -kHighsInf {
            push(true, i, name, true, lp.row_lower_[i], act);
        }
        if lp.row_upper_[i] < kHighsInf {
            push(true, i, name, false, lp.row_upper_[i], act);
        }
    }
    for j in 0..lp.num_col_ as usize {
        let name = lp.col_names_.get(j).map(String::as_str).unwrap_or("");
        if lp.col_lower_[j] > -kHighsInf {
            push(false, j, name, true, lp.col_lower_[j], x[j]);
        }
        if lp.col_upper_[j] < kHighsInf {
            push(false, j, name, false, lp.col_upper_[j], x[j]);
        }
    }
    out
}

/// Stage 1 of `min-count`: the indicator MIP, built on a fresh `Highs` with the
/// ported `add_col`/`add_row`. Returns the relaxed point of the MIP.
#[allow(clippy::needless_range_loop)]
fn solve_min_count_mip(lp: &HighsLp, w: &Weights, a: &Args) -> Outcome {
    let n = lp.num_col_ as usize;
    let m = lp.num_row_ as usize;
    let rows = row_entries(lp);
    let mut notes = Vec::new();
    let mut fallback_cols = 0usize;
    let mut fallback_rows = 0usize;

    // Column domain of the auxiliary model: a relaxable finite bound is opened by
    // `big_m` (a relaxed column has no natural limit; Chinneck 2008 ch. 6 bounds
    // the elastic variable the same way).
    let mut col_lo = lp.col_lower_.clone();
    let mut col_up = lp.col_upper_.clone();
    let mut col_ind: Vec<(usize, bool, f64)> = Vec::new(); // (col, lower side, cost)
    for j in 0..n {
        if lp.col_lower_[j] > -kHighsInf && w.col_lo[j] >= 0.0 {
            col_lo[j] = lp.col_lower_[j] - a.big_m;
            col_ind.push((j, true, 1.0));
        }
        if lp.col_upper_[j] < kHighsInf && w.col_up[j] >= 0.0 {
            col_up[j] = lp.col_upper_[j] + a.big_m;
            col_ind.push((j, false, 1.0));
        }
    }
    fallback_cols += col_ind.len();

    // Row sides: (row, lower side, cost, M).
    let mut row_ind: Vec<(usize, bool, f64, f64)> = Vec::new();
    for i in 0..m {
        if w.row[i] < 0.0 {
            continue;
        }
        let (mut minact, mut maxact) = (0.0f64, 0.0f64);
        for &(j, v) in &rows[i] {
            if v == 0.0 {
                continue;
            }
            let (lo, up) = if v > 0.0 {
                (col_lo[j], col_up[j])
            } else {
                (col_up[j], col_lo[j])
            };
            minact += v * lo;
            maxact += v * up;
        }
        // `v * inf` is +-inf and survives the sum (minact only ever adds -inf
        // terms, maxact only +inf); a stored zero coefficient is skipped above.
        if lp.row_lower_[i] > -kHighsInf {
            let gap = lp.row_lower_[i] - minact;
            if gap > 0.0 {
                let mm = if gap.is_finite() {
                    gap
                } else {
                    fallback_rows += 1;
                    a.big_m
                };
                row_ind.push((i, true, 1.0, mm));
            }
        }
        if lp.row_upper_[i] < kHighsInf {
            let gap = maxact - lp.row_upper_[i];
            if gap > 0.0 {
                let mm = if gap.is_finite() {
                    gap
                } else {
                    fallback_rows += 1;
                    a.big_m
                };
                row_ind.push((i, false, 1.0, mm));
            }
        }
    }
    if fallback_cols > 0 {
        notes.push(format!(
            "{fallback_cols} column bound(s) are relaxable: their relaxed range is limited to big-M {} \
             (--big-m); a relaxation larger than that is not considered",
            a.big_m
        ));
    }
    if fallback_rows > 0 {
        notes.push(format!(
            "{fallback_rows} row side(s) have an infinite activity bound: big-M {} used (--big-m); \
             the minimum count is exact only if no true violation exceeds it",
            a.big_m
        ));
    }

    let mut highs = solver(a);
    // M z must not leak through the integrality tolerance.
    highs.set_option_value_double("mip_feasibility_tolerance", 1e-9);
    for j in 0..n {
        if highs.add_col(0.0, col_lo[j], col_up[j], 0, &[], &[]) == HighsStatus::kError {
            return Outcome::Failed("add_col failed".into());
        }
        if lp.integrality_.get(j) == Some(&HighsVarType::kInteger) {
            highs.change_col_integrality(j as HighsInt, HighsVarType::kInteger);
        }
    }
    let z0 = n;
    let nz = col_ind.len() + row_ind.len();
    for k in 0..nz {
        let cost = if k < col_ind.len() {
            col_ind[k].2
        } else {
            row_ind[k - col_ind.len()].2
        };
        highs.add_col(cost, 0.0, 1.0, 0, &[], &[]);
        highs.change_col_integrality((z0 + k) as HighsInt, HighsVarType::kInteger);
    }
    let add_row = |highs: &mut Highs, lo: f64, up: f64, idx: Vec<HighsInt>, val: Vec<f64>| {
        highs.add_row(lo, up, idx.len() as HighsInt, &idx, &val)
    };
    // Column bound k: x_j + M z >= l   (or x_j - M z <= u).
    for (k, &(j, lower, _)) in col_ind.iter().enumerate() {
        let (idx, val) = (vec![j as HighsInt, (z0 + k) as HighsInt], (vec![1.0]));
        if lower {
            let mut v = val.clone();
            v.push(a.big_m);
            add_row(&mut highs, lp.col_lower_[j], kHighsInf, idx, v);
        } else {
            let mut v = val.clone();
            v.push(-a.big_m);
            add_row(&mut highs, -kHighsInf, lp.col_upper_[j], idx, v);
        }
    }
    // Row sides, one auxiliary row per relaxable side; sides that are not
    // relaxable (negative weight) stay hard in one shared row.
    for i in 0..m {
        let base_idx: Vec<HighsInt> = rows[i].iter().map(|&(j, _)| j as HighsInt).collect();
        let base_val: Vec<f64> = rows[i].iter().map(|&(_, v)| v).collect();
        if w.row[i] < 0.0 {
            add_row(
                &mut highs,
                lp.row_lower_[i],
                lp.row_upper_[i],
                base_idx,
                base_val,
            );
            continue;
        }
        for (k, &(ri, lower, _, mm)) in row_ind.iter().enumerate() {
            if ri != i {
                continue;
            }
            let mut idx = base_idx.clone();
            let mut val = base_val.clone();
            idx.push((z0 + col_ind.len() + k) as HighsInt);
            if lower {
                val.push(mm);
                add_row(&mut highs, lp.row_lower_[i], kHighsInf, idx, val);
            } else {
                val.push(-mm);
                add_row(&mut highs, -kHighsInf, lp.row_upper_[i], idx, val);
            }
        }
    }
    let run = highs.run();
    let status = highs.get_model_status();
    if run == HighsStatus::kError {
        return Outcome::Failed("the auxiliary MIP failed to run".into());
    }
    match status {
        HighsModelStatus::kInfeasible | HighsModelStatus::kUnboundedOrInfeasible => {
            return Outcome::NotRelaxable(notes);
        }
        HighsModelStatus::kOptimal => {}
        other => {
            if highs.get_info().primal_solution_status != 2 {
                return Outcome::Failed(format!(
                    "the auxiliary MIP stopped with status {} and no solution",
                    crate::envelope::result_status(other)
                ));
            }
            notes.push(format!(
                "the auxiliary MIP stopped with status {} before proving optimality: the count may not be minimal",
                crate::envelope::result_status(other)
            ));
        }
    }
    let x = highs.get_solution().col_value[..n].to_vec();
    let obj = highs.get_info().objective_function_value;
    Outcome::Point(x, obj, notes)
}

/// Re-solves the original model with the edits applied. `(verified, status,
/// objective)`.
fn verify(lp: &HighsLp, items: &[Item], a: &Args) -> (bool, String, Option<f64>) {
    let mut lp = lp.clone();
    for it in items {
        let (lo, up) = if it.is_row {
            (&mut lp.row_lower_[it.index], &mut lp.row_upper_[it.index])
        } else {
            (&mut lp.col_lower_[it.index], &mut lp.col_upper_[it.index])
        };
        if it.lower_side {
            *lo = it.relaxed_to;
        } else {
            *up = it.relaxed_to;
        }
    }
    let mut highs = solver(a);
    if highs.pass_model(lp) == HighsStatus::kError || highs.run() == HighsStatus::kError {
        return (false, "Error".into(), None);
    }
    let st = highs.get_model_status();
    let ok = matches!(
        st,
        HighsModelStatus::kOptimal | HighsModelStatus::kUnbounded
    );
    let obj = (st == HighsModelStatus::kOptimal).then(|| highs.get_info().objective_function_value);
    (ok, crate::envelope::result_status(st).to_string(), obj)
}

fn fnum(x: f64) -> String {
    format!("{x}")
}

/// The relaxation of the chosen mode (the two stages of `min-count` included).
fn solve_mode(lp: &HighsLp, w: &Weights, a: &Args) -> Outcome {
    let weight_fn = |w: &Weights| {
        let w = w.clone();
        move |is_row: bool, i: usize, lower: bool| weight_of(&w, is_row, i, lower)
    };
    match a.mode {
        Mode::Weighted => solve_weighted(lp, w, a),
        Mode::MinCount => match solve_min_count_mip(lp, w, a) {
            Outcome::Point(x1, obj1, mut n1) => {
                // Stage 2: minimal violation for the chosen set.
                let chosen = items_of(lp, &x1, &weight_fn(w));
                let mut w2 = Weights {
                    row: vec![-1.0; w.row.len()],
                    col_lo: vec![-1.0; w.col_lo.len()],
                    col_up: vec![-1.0; w.col_up.len()],
                };
                for it in &chosen {
                    let wt = weight_of(w, it.is_row, it.index, it.lower_side);
                    if it.is_row {
                        w2.row[it.index] = wt;
                    } else if it.lower_side {
                        w2.col_lo[it.index] = wt;
                    } else {
                        w2.col_up[it.index] = wt;
                    }
                }
                let mut a2 = a.clone();
                a2.gl = -1.0;
                a2.gu = -1.0;
                a2.gr = -1.0;
                match solve_weighted(lp, &w2, &a2) {
                    Outcome::Point(x2, _, _)
                        if items_of(lp, &x2, &weight_fn(w)).len() <= chosen.len() =>
                    {
                        n1.push("stage 2: violations minimised over the chosen constraints".into());
                        Outcome::Point(x2, obj1, n1)
                    }
                    _ => {
                        n1.push(
                            "stage 2 (violation minimisation) did not apply: stage 1 point kept"
                                .into(),
                        );
                        Outcome::Point(x1, obj1, n1)
                    }
                }
            }
            other => other,
        },
    }
}

/// One constraint side the minimum-count relaxation gives up (`diagnose`).
pub(crate) struct RelaxedCell {
    pub is_row: bool,
    pub index: usize,
    pub lower_side: bool,
    pub original: f64,
    pub relaxed_to: f64,
    pub violation: f64,
}

/// The `--mode min-count` relaxation of `lp` with default weights, run in
/// process: `Ok(Some((cells, verified)))`, `Ok(None)` when the model cannot be
/// relaxed, `Err` on solver failure. Shared with `diagnose`.
pub(crate) fn min_count_relaxation(
    lp: &HighsLp,
    time_limit: f64,
) -> Result<Option<(Vec<RelaxedCell>, bool)>, String> {
    let tl = time_limit.to_string();
    let a = parse_args(&[
        "-".to_string(),
        "--mode".into(),
        "min-count".into(),
        "--time-limit".into(),
        tl,
    ])?;
    let w = resolve_weights(lp, &a, None)?;
    match solve_mode(lp, &w, &a) {
        Outcome::Failed(m) => Err(m),
        Outcome::NotRelaxable(_) => Ok(None),
        Outcome::Point(x, _, _) => {
            let items = items_of(lp, &x, &|r, i, l| weight_of(&w, r, i, l));
            let (ok, _, _) = verify(lp, &items, &a);
            Ok(Some((
                items
                    .iter()
                    .map(|it| RelaxedCell {
                        is_row: it.is_row,
                        index: it.index,
                        lower_side: it.lower_side,
                        original: it.original,
                        relaxed_to: it.relaxed_to,
                        violation: it.violation,
                    })
                    .collect(),
                ok,
            )))
        }
    }
}

pub fn run_feasrelax(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let a = match parse_args(args) {
        Ok(a) => a,
        Err(e) => return usage(wj, &e),
    };
    let lp = match read_mip(&a.path) {
        Ok(lp) => lp,
        Err(e) => return crate::diag::fail_read_str(a.json, &a.path, &e),
    };
    let file: Option<Value> = match &a.penalties {
        None => None,
        Some(p) => match std::fs::read_to_string(p) {
            Err(e) => {
                return crate::diag::fail_msg(
                    a.json,
                    "RS0010",
                    3,
                    format!("rusolve feasrelax: read {p}: {e}"),
                );
            }
            Ok(t) => match serde_json::from_str(&t) {
                Ok(v) => Some(v),
                Err(e) => {
                    return crate::diag::fail_msg(
                        a.json,
                        "RS0011",
                        3,
                        format!("rusolve feasrelax: parse {p}: {e}"),
                    );
                }
            },
        },
    };
    match feasrelax_core(&lp, &a, file.as_ref()) {
        Ok((doc, text)) => {
            if a.json {
                println!("{doc}");
            } else {
                print!("{text}");
            }
            ExitCode::SUCCESS
        }
        Err(e) => crate::diag::fail_core(a.json, "feasrelax", &e),
    }
}

/// `Args` of [`feasrelax_core`] for a model that is already loaded (the MCP
/// server): the same parser as the CLI, over the words the tool arguments
/// map to. `path` is only the name reported in `input.path`.
pub(crate) fn core_args(
    path: &str,
    flags: &[String],
    inline_penalties: bool,
) -> Result<Args, String> {
    let mut argv = vec![path.to_string()];
    argv.extend(flags.iter().cloned());
    let mut a = parse_args(&argv)?;
    if inline_penalties {
        // `weights.penalties_file` reports where the weights came from.
        a.penalties = Some("(inline)".to_string());
    }
    Ok(a)
}

/// The body of `rusolve feasrelax` on a loaded model: the `rusolve.feasrelax/1`
/// object and the text report. `file` is the parsed `--penalties` document.
/// Shared by the CLI and the MCP `feasrelax` tool.
pub(crate) fn feasrelax_core(
    lp: &HighsLp,
    a: &Args,
    file: Option<&Value>,
) -> Result<(Value, String), crate::cli_core::CoreError> {
    use crate::cli_core::CoreError;
    let a = a.clone();
    let lp = lp.clone();
    let w = match resolve_weights(&lp, &a, file) {
        Ok(w) => w,
        Err(e) => return Err(CoreError::usage(e)),
    };
    let started = std::time::Instant::now();
    let mut notes: Vec<String> = Vec::new();
    let weight_fn = |w: &Weights| {
        let w = w.clone();
        move |is_row: bool, i: usize, lower: bool| weight_of(&w, is_row, i, lower)
    };

    let outcome = solve_mode(&lp, &w, &a);
    let mode_name = if a.mode == Mode::Weighted {
        "weighted"
    } else {
        "min-count"
    };

    let (status, x, solver_obj, items) = match outcome {
        Outcome::Failed(msg) => return Err(CoreError::internal(msg)),
        Outcome::NotRelaxable(n) => {
            notes.extend(n);
            ("NotRelaxable", None, None, Vec::new())
        }
        Outcome::Point(x, obj, n) => {
            notes.extend(n);
            let items = items_of(&lp, &x, &weight_fn(&w));
            let status = if items.is_empty() {
                "Feasible"
            } else {
                "Relaxed"
            };
            (status, Some(x), Some(obj), items)
        }
    };
    let (verified, vstatus, vobj) = if status == "NotRelaxable" {
        (false, "NotChecked".to_string(), None)
    } else {
        verify(&lp, &items, &a)
    };
    let total_violation: f64 = items.iter().map(|i| i.violation).sum();
    let total_penalty: f64 = items.iter().map(|i| i.weight * i.violation).sum();
    // The relaxed point under the original objective (cost . x + offset).
    let obj_original = x.as_ref().map(|x| {
        lp.col_cost_
            .iter()
            .zip(x)
            .map(|(cost, v)| cost * v)
            .sum::<f64>()
            + lp.offset_
    });

    let kind = |it: &Item| if it.is_row { "row" } else { "column" };
    let side = |it: &Item| if it.lower_side { "lower" } else { "upper" };
    let edit_text = |it: &Item| {
        format!(
            "change {} bound of {} {} from {} to {}",
            side(it),
            kind(it),
            it.name,
            fnum(it.original),
            fnum(it.relaxed_to)
        )
    };
    let relaxed_json: Vec<Value> = items
        .iter()
        .map(|it| {
            json!({
                "kind": kind(it), "name": it.name, "index": it.index, "bound": side(it),
                "original_bound": it.original, "relaxed_to": it.relaxed_to,
                "violation": it.violation, "weight": it.weight,
                "penalty": it.weight * it.violation,
            })
        })
        .collect();
    let edits_json: Vec<Value> = items
        .iter()
        .map(|it| {
            json!({
                "kind": kind(it), "name": it.name, "bound": side(it),
                "from": it.original, "to": it.relaxed_to, "text": edit_text(it),
            })
        })
        .collect();

    let n_int = lp
        .integrality_
        .iter()
        .filter(|t| **t == HighsVarType::kInteger)
        .count();
    let doc = json!({
        "schema_version": SCHEMA_VERSION,
        "input": {"path": a.path, "rows": lp.num_row_, "cols": lp.num_col_, "integer_cols": n_int},
        "mode": mode_name,
        "weights": {
            "global_lower_penalty": a.gl, "global_upper_penalty": a.gu,
            "global_rhs_penalty": a.gr, "penalties_file": a.penalties,
        },
        "status": status,
        "relaxed": relaxed_json,
        "edits": edits_json,
        "summary": {
            "relaxed_count": items.len(),
            "total_violation": total_violation,
            "total_penalty": total_penalty,
            "solver_objective": solver_obj,
            "objective_original": obj_original,
            "objective_after_edits": vobj,
        },
        "verified": verified,
        "verification": {"status": vstatus},
        "notes": notes,
        "time_s": started.elapsed().as_secs_f64(),
    });
    let mut text = String::new();
    macro_rules! tprintln {
        ($($arg:tt)*) => {{
            text.push_str(&format!($($arg)*));
            text.push('\n');
        }};
    }
    tprintln!("feasrelax {}  mode {mode_name}  status {status}", a.path);
    match status {
        "NotRelaxable" => tprintln!(
            "no relaxation exists: the bounds that may not be relaxed (negative weight) conflict"
        ),
        "Feasible" => tprintln!("the model is feasible as given: nothing to relax"),
        _ => {
            tprintln!("relaxed {} constraint(s):", items.len());
            for it in &items {
                tprintln!(
                    "  {} {} {} bound {} -> {}  violation {}  weight {}  penalty {}",
                    kind(it),
                    it.name,
                    side(it),
                    fnum(it.original),
                    fnum(it.relaxed_to),
                    fnum(it.violation),
                    fnum(it.weight),
                    fnum(it.weight * it.violation)
                );
            }
            tprintln!(
                "total violation {}  total penalty {}",
                fnum(total_violation),
                fnum(total_penalty)
            );
            if let Some(o) = obj_original {
                tprintln!(
                    "objective of the relaxed point under the original objective: {}",
                    fnum(o)
                );
            }
            tprintln!("edits:");
            for it in &items {
                tprintln!("  {}", edit_text(it));
            }
        }
    }
    if status != "NotRelaxable" {
        let o = vobj.map(fnum).unwrap_or_else(|| "n/a".into());
        tprintln!(
            "verified {verified} (original model with the edits applied: {vstatus}, objective {o})"
        );
    }
    for n in &notes {
        tprintln!("note: {n}");
    }
    Ok((doc, text))
}
