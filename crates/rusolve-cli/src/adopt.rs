//! The diagnostics `solve` and `mip` attach to a
//! run, as warnings on stderr and as the envelope's `diagnostics` array (`rusolve.result/1`,
//! additive: present only when there is one).
//!
//! * Before the solve: the cheap lint rules (`lint::warnings`: big-M against the implied bound,
//!   coefficient range, huge or free columns). They read a copy of the model, so the search cannot
//!   move; `--no-warnings` turns them off.
//! * After the solve (computed only after the status is known, with a bounded budget): an
//!   infeasible model carries a small irreducible infeasible set of rows (`HighsIis`, Chinneck and
//!   Dravnieks 1991, "Locating minimal infeasible constraint sets in linear programs", ORSA J.
//!   Computing 3(2)) and an unbounded one a primal ray (the same ray `rusolve explain` verifies);
//!   both name `explain` and `feasrelax`. A returned solution whose row or bound violation exceeds
//!   the tolerance is a warning naming the row.
//!
//! Layout of the diagnostics: rustc (<https://rustc-dev-guide.rust-lang.org/diagnostics.html>).

use rusolve_core::diagnostic::{Diagnostic, Location};
use rusolve_ekk::lp_data::highs_iis::HighsIis;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::Value;

/// `true` when the command line carries `--no-warnings`.
pub fn no_warnings(args: &[String]) -> bool {
    args.iter().any(|a| a == "--no-warnings")
}

/// Largest model (nonzeros) the post-solve IIS or ray is computed for; above it the help line
/// names `rusolve explain` instead.
const POST_MAX_NNZ: usize = 500_000;
/// Wall budget of the post-solve re-solve, seconds.
const POST_BUDGET_S: f64 = 10.0;
/// Members of an infeasible set or ray support named in the notes.
const NAMED: usize = 12;
/// A returned solution violates a row or bound by more than this (absolute).
const VIOLATION_TOL: f64 = 1e-6;

/// The shared post-solve budget (seconds) for a run that took `elapsed` seconds under the user's
/// `--time-limit` (`inf` when none): half the solve time (at least 0.25 s), never above
/// [`POST_BUDGET_S`], and never more than what is left of the user's limit.
pub fn budget_for(elapsed: f64, time_limit: f64) -> f64 {
    let b = (0.5 * elapsed).clamp(0.25, POST_BUDGET_S);
    if time_limit.is_finite() {
        b.min((time_limit - elapsed).max(0.05))
    } else {
        b
    }
}

/// Print `diags` as rustc-style text on stderr.
pub fn print_warnings(diags: &[Diagnostic]) {
    for d in diags {
        eprint!("{}", d.render_text());
    }
}

/// Add the `diagnostics` array to a result envelope given as JSON text; `text` unchanged when
/// there is nothing to add or it is not an object.
pub fn with_diagnostics(text: &str, diags: &[Diagnostic]) -> String {
    if diags.is_empty() {
        return text.to_string();
    }
    match serde_json::from_str::<Value>(text) {
        Ok(Value::Object(mut m)) => {
            let arr: Vec<Value> = diags
                .iter()
                .filter_map(|d| serde_json::from_str(&d.to_json()).ok())
                .collect();
            m.insert("diagnostics".into(), Value::Array(arr));
            Value::Object(m).to_string()
        }
        _ => text.to_string(),
    }
}

/// Diagnostics after a solve: the infeasible set, the ray, or the violated row.
/// `lp` is the model as solved (the user's rows), `solution` what the solver returned. `warn` is
/// false under `--no-warnings`, which drops the violation warning; the infeasible set and the ray
/// are the result of the run, not warnings about the model, and stay.
pub fn post_solve(
    path: &str,
    lp: &HighsLp,
    status: HighsModelStatus,
    solution: &HighsSolution,
    warn: bool,
    budget: f64,
) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    match status {
        HighsModelStatus::kInfeasible => out.extend(infeasible(path, lp, budget)),
        HighsModelStatus::kUnbounded | HighsModelStatus::kUnboundedOrInfeasible => {
            out.extend(unbounded(path, lp, budget));
        }
        // only a status that carries a solution can have a violated row; the phase-1 point of an
        // infeasible LP is not a returned solution
        _ if warn => out.extend(violation(lp, solution)),
        _ => {}
    }
    out
}

fn small(lp: &HighsLp) -> bool {
    lp.a_matrix_.value_.len() <= POST_MAX_NNZ
}

/// Solve the LP relaxation with the irreducible-IIS strategy (the setting `rusolve explain`
/// uses), within the budget.
fn relaxation(lp: &HighsLp, budget: f64) -> Option<(Highs, HighsModelStatus)> {
    let mut lp = lp.clone();
    lp.integrality_.clear();
    let mut h = Highs::with_options(crate::base_options());
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_int("iis_strategy", 6);
    h.set_option_value_double("time_limit", budget);
    // `get_iis` resets time_limit to iis_time_limit (default inf): the caller sets it to what
    // is left of the shared budget before the IIS call
    if h.pass_model(lp) == HighsStatus::kError || h.run() == HighsStatus::kError {
        return None;
    }
    let mut status = h.get_model_status();
    if status == HighsModelStatus::kUnboundedOrInfeasible {
        h.set_option_value_str("presolve", "off");
        if h.run() == HighsStatus::kError {
            return None;
        }
        status = h.get_model_status();
    }
    Some((h, status))
}

fn infeasible(path: &str, lp: &HighsLp, iis_budget: f64) -> Vec<Diagnostic> {
    infeasible_from(path, lp, iis_budget, None)
}

/// [`infeasible`], reusing a relaxation re-solve that already ran (`pre`), as when the solver
/// stopped at "unbounded or infeasible" and the re-solve says infeasible.
fn infeasible_from(
    path: &str,
    lp: &HighsLp,
    iis_budget: f64,
    pre: Option<(Highs, HighsModelStatus, std::time::Instant)>,
) -> Vec<Diagnostic> {
    let mut d = Diagnostic::warning("RS2010", "the model is infeasible");
    let helps = |d: Diagnostic| {
        d.with_help(format!(
            "`rusolve explain {path}` verifies the set and shows the arithmetic"
        ))
        .with_help(format!(
            "`rusolve feasrelax {path}` finds the cheapest relaxation that makes it feasible"
        ))
    };
    if !small(lp) {
        return vec![helps(d.with_note(
            "the model is too large to compute an infeasible set after the solve",
        ))];
    }
    // the shared budget runs from `t0`: the start of the relaxation, also when `pre` ran it
    let (pre, t0) = match pre {
        Some((h, st, t0)) => (Some((h, st)), t0),
        None => (None, std::time::Instant::now()),
    };
    let Some((mut h, st)) = pre.or_else(|| relaxation(lp, iis_budget)) else {
        return vec![helps(d)];
    };
    if st == HighsModelStatus::kOptimal || st == HighsModelStatus::kUnbounded {
        if !lp
            .integrality_
            .iter()
            .any(|t| *t != HighsVarType::kContinuous)
        {
            return vec![helps(d.with_note(
                "the model has no integer columns and its LP is feasible: the conflict involves the general constraints (or their linearization)",
            ))];
        }
        return vec![helps(d.with_note(
            "the LP relaxation is feasible: the conflict involves integrality (`rusolve explain --mip-iis` searches for it)",
        ))];
    }
    if st != HighsModelStatus::kInfeasible {
        return vec![helps(d.with_note(format!(
            "the post-solve budget ({iis_budget:.2} s) ran out before the LP relaxation was classified; `rusolve explain {path}` has no time limit"
        )))];
    }
    // the IIS gets what the relaxation left of the shared budget
    h.set_option_value_double(
        "iis_time_limit",
        (iis_budget - t0.elapsed().as_secs_f64()).max(0.05),
    );
    let mut iis = HighsIis::default();
    let ist = h.get_iis(&mut iis);
    let (rows, cols) = crate::explain::iis_members(&iis, ist);
    if rows.is_empty() && cols.is_empty() {
        return vec![helps(d)];
    }
    let rname = |i: usize| crate::explain::row_name(lp, i);
    let names: Vec<String> = rows.iter().map(|m| rname(m.index)).collect();
    let irreducible = iis.status_ == kIisModelStatusIrreducible;
    d.message = format!(
        "the model is infeasible: {}{} row(s) conflict ({}{})",
        if irreducible { "" } else { "a subset of " },
        rows.len(),
        names
            .iter()
            .take(NAMED)
            .cloned()
            .collect::<Vec<_>>()
            .join(", "),
        if names.len() > NAMED { ", ..." } else { "" }
    );
    if let Some(m) = rows.first() {
        let (lo, up) = (lp.row_lower_[m.index], lp.row_upper_[m.index]);
        d = d.with_location(Location::Row {
            name: rname(m.index),
            rendered: format!("{lo} <= row <= {up}"),
        });
    }
    if irreducible {
        d = d.with_note(format!(
            "irreducible infeasible set: {} row(s), {} column bound(s); dropping any one member makes the rest feasible",
            rows.len(),
            cols.len()
        ));
    } else {
        d = d.with_note(format!(
            "an infeasible subset of {} row(s), {} column bound(s), not reduced: the post-solve budget ran out (`rusolve explain {path}` computes the irreducible set)",
            rows.len(),
            cols.len()
        ));
    }
    for m in rows.iter().take(NAMED) {
        let (lo, up) = (lp.row_lower_[m.index], lp.row_upper_[m.index]);
        d = d.with_note(format!("row {}: {lo} <= activity <= {up}", rname(m.index)));
    }
    if !cols.is_empty() {
        let cn: Vec<String> = cols
            .iter()
            .take(NAMED)
            .map(|m| crate::explain::col_name(lp, m.index))
            .collect();
        d = d.with_note(format!("column bounds in the set: {}", cn.join(", ")));
    }
    vec![helps(d)]
}

fn unbounded(path: &str, lp: &HighsLp, budget: f64) -> Vec<Diagnostic> {
    let mut d = Diagnostic::warning("RS2011", "the model is unbounded");
    let helps = |d: Diagnostic| {
        d.with_help(format!(
            "`rusolve explain {path}` verifies the ray and shows the objective rate"
        ))
        .with_help("bound the columns on the ray or add the missing constraint; `rusolve lint` flags free and huge-bound columns")
    };
    if !small(lp) {
        return vec![helps(d.with_note(
            "the model is too large to compute a ray after the solve",
        ))];
    }
    let t0 = std::time::Instant::now();
    let Some((mut h, st)) = relaxation(lp, budget) else {
        return vec![helps(d)];
    };
    // "unbounded or infeasible" that the re-solve without presolve settles as infeasible is
    // reported as RS2010, not as a bare RS2011
    if st == HighsModelStatus::kInfeasible {
        return infeasible_from(path, lp, budget, Some((h, st, t0)));
    }
    if st != HighsModelStatus::kUnbounded {
        return vec![helps(d)];
    }
    let mut has_ray = false;
    let mut dir = vec![0.0f64; lp.num_col_ as usize];
    let rs = h.get_primal_ray(&mut has_ray, Some(&mut dir));
    if rs == HighsStatus::kError || !has_ray {
        return vec![helps(d)];
    }
    let mut nz: Vec<(usize, f64)> = dir
        .iter()
        .enumerate()
        .filter(|(_, v)| **v != 0.0)
        .map(|(j, v)| (j, *v))
        .collect();
    nz.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
    d.message = format!(
        "the model is unbounded: the objective improves without limit along a ray over {} column(s) ({})",
        nz.len(),
        nz.iter()
            .take(5)
            .map(|&(j, _)| crate::explain::col_name(lp, j))
            .collect::<Vec<_>>()
            .join(", ")
    );
    for &(j, v) in nz.iter().take(NAMED) {
        d = d.with_note(format!(
            "ray: {} moves by {v:.6e}",
            crate::explain::col_name(lp, j)
        ));
    }
    vec![helps(d)]
}

/// A returned solution that violates a row or a column bound by more than [`VIOLATION_TOL`].
fn violation(lp: &HighsLp, sol: &HighsSolution) -> Vec<Diagnostic> {
    if !sol.value_valid
        || sol.row_value.len() != lp.num_row_ as usize
        || sol.col_value.len() != lp.num_col_ as usize
    {
        return Vec::new();
    }
    let mut worst: Option<(f64, String, String)> = None;
    let mut count = 0usize;
    let mut consider = |v: f64, name: String, what: String| {
        if v > VIOLATION_TOL {
            count += 1;
            if worst.as_ref().is_none_or(|w| v > w.0) {
                worst = Some((v, name, what));
            }
        }
    };
    for i in 0..lp.num_row_ as usize {
        let a = sol.row_value[i];
        let v = (lp.row_lower_[i] - a).max(a - lp.row_upper_[i]);
        consider(
            v,
            crate::explain::row_name(lp, i),
            format!(
                "{} <= activity {a} <= {}",
                lp.row_lower_[i], lp.row_upper_[i]
            ),
        );
    }
    for j in 0..lp.num_col_ as usize {
        let x = sol.col_value[j];
        let v = (lp.col_lower_[j] - x).max(x - lp.col_upper_[j]);
        consider(
            v,
            crate::explain::col_name(lp, j),
            format!("{} <= value {x} <= {}", lp.col_lower_[j], lp.col_upper_[j]),
        );
    }
    let Some((v, name, what)) = worst else {
        return Vec::new();
    };
    vec![
        Diagnostic::warning(
            "RS2012",
            format!("the returned solution violates '{name}' by {v:.3e}"),
        )
        .with_location(Location::Row {
            name: name.clone(),
            rendered: what,
        })
        .with_note(format!(
            "{count} row(s) or bound(s) exceed the tolerance {VIOLATION_TOL:e}"
        ))
        .with_help(
            "`rusolve check <model> --solution-file <sol>` re-verifies the solution independently",
        )
        .with_help("tighten `primal_feasibility_tolerance` or rescale the model (`rusolve lint`)"),
    ]
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    /// x in [0,1], one row x >= 2: infeasible.
    fn tiny() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![1.0];
        lp.row_lower_ = vec![2.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        lp
    }

    fn text(d: &[Diagnostic]) -> String {
        d.iter().map(|x| x.render_text()).collect()
    }

    /// n supply rows `<= 10`, n demand rows `>= 10.1` over n*n arcs: infeasible, and the IIS
    /// needs many reduction steps.
    fn transport(n: usize) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = (n * n) as _;
        lp.num_row_ = (2 * n) as _;
        lp.col_cost_ = vec![1.0; n * n];
        lp.col_lower_ = vec![0.0; n * n];
        lp.col_upper_ = vec![kHighsInf; n * n];
        lp.row_lower_ = [vec![-kHighsInf; n], vec![10.1; n]].concat();
        lp.row_upper_ = [vec![10.0; n], vec![kHighsInf; n]].concat();
        lp.a_matrix_.num_col_ = (n * n) as _;
        lp.a_matrix_.num_row_ = (2 * n) as _;
        for s in 0..n {
            for d in 0..n {
                lp.a_matrix_.start_.push(lp.a_matrix_.index_.len() as _);
                lp.a_matrix_.index_.extend([s as i32, (n + d) as i32]);
                lp.a_matrix_.value_.extend([1.0, 1.0]);
            }
        }
        lp.a_matrix_.start_.push(lp.a_matrix_.index_.len() as _);
        lp
    }

    #[test]
    fn time_limited_iis_is_not_called_irreducible() {
        let t = text(&infeasible("t.lp", &transport(40), 1e-9));
        assert!(!t.contains("irreducible infeasible set"), "{t}");
        assert!(!t.contains("dropping any one member"), "{t}");
    }

    #[test]
    fn exhausted_budget_does_not_claim_a_feasible_relaxation() {
        let t = text(&infeasible("t.lp", &transport(60), 0.0));
        assert!(!t.contains("relaxation is feasible"), "{t}");
        assert!(t.contains("budget"), "{t}");
    }

    #[test]
    fn budget_is_capped_by_the_user_limit() {
        assert!(budget_for(0.1, f64::INFINITY) == 0.25);
        assert!(budget_for(100.0, f64::INFINITY) == POST_BUDGET_S);
        assert!(budget_for(2.0, 3.0) <= 1.0 + 1e-9);
    }

    #[test]
    fn infeasible_status_has_no_violation_warning() {
        let lp = tiny();
        let mut sol = HighsSolution::default();
        sol.value_valid = true;
        sol.col_value = vec![1.0];
        sol.row_value = vec![1.0];
        let d = post_solve("t.lp", &lp, HighsModelStatus::kInfeasible, &sol, true, 1.0);
        assert!(d.iter().all(|x| x.code != "RS2012"));
        let d = post_solve("t.lp", &lp, HighsModelStatus::kOptimal, &sol, true, 1.0);
        assert!(d.iter().any(|x| x.code == "RS2012"));
    }

    /// The solver's "unbounded or infeasible" stop, settled as infeasible by the relaxation
    /// re-solve, is reported as RS2010, not as a bare RS2011.
    #[test]
    fn unbounded_or_infeasible_settled_infeasible_is_rs2010() {
        let lp = tiny();
        let sol = HighsSolution::default();
        let d = post_solve(
            "t.lp",
            &lp,
            HighsModelStatus::kUnboundedOrInfeasible,
            &sol,
            true,
            5.0,
        );
        assert!(d.iter().any(|x| x.code == "RS2010"), "{}", text(&d));
        assert!(d.iter().all(|x| x.code != "RS2011"));
    }
}
