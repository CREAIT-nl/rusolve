//! Typed Python exceptions built from the diag-core `Diagnostic`
//! (docs/design/diagnostics.md). The classes are defined in `python/rusolve/errors.py`; this module
//! builds the diagnostic (code, message, help, the Python line that made the call) and hands its
//! `rusolve.diagnostics/1` document to the class.

use pyo3::prelude::*;
use rusolve_core::diagnostic::{Diagnostic, Location, diagnostics_json};

use rusolve_ekk::lp_data::hconst::HighsModelStatus;
use rusolve_ekk::prelude::HighsLp;

use crate::expr::MIXED_MODELS;

/// The line of Python that made the current call (frames inside `rusolve` itself skipped).
pub(crate) fn caller_location(py: Python<'_>) -> Option<Location> {
    let m = py.import("rusolve.errors").ok()?;
    let t: Option<(String, usize, String)> =
        m.getattr("_caller").ok()?.call0().ok()?.extract().ok();
    t.map(|(path, line, source_line)| Location::Python {
        path,
        line,
        source_line,
    })
}

/// Raise `class` (a name in `rusolve.errors`) carrying `d`; its location is the calling line
/// of Python unless `d` already has one.
pub(crate) fn raise(py: Python<'_>, class: &str, mut d: Diagnostic) -> PyErr {
    if d.location.is_none()
        && let Some(l) = caller_location(py)
    {
        d.location = Some(l);
    }
    let built = (|| -> PyResult<PyErr> {
        let cls = py.import("rusolve.errors")?.getattr(class)?;
        let e = cls.call_method1("_from_json", (d.message.clone(), diagnostics_json(&[d], 1)))?;
        Ok(PyErr::from_value(e))
    })();
    built.unwrap_or_else(|e| e)
}

/// A model-building error. RS4003 when the message is the mixed-Models one, else RS4001.
pub(crate) fn model_err(msg: impl Into<String>) -> PyErr {
    let msg = msg.into();
    Python::attach(|py| {
        let d = if msg.starts_with(MIXED_MODELS) {
            Diagnostic::error("RS4003", msg.clone()).with_help(
                "create every Var, LinExpr and Constraint from the Model it is used with",
            )
        } else {
            Diagnostic::error("RS4001", msg.clone()).with_help(
                "the message names the argument; the docstring of the method lists what it accepts",
            )
        };
        raise(py, "ModelError", d)
    })
}

/// An index outside the model (RS4002).
pub(crate) fn index_err(msg: impl Into<String>) -> PyErr {
    let msg = msg.into();
    Python::attach(|py| {
        raise(
            py,
            "ModelIndexError",
            Diagnostic::error("RS4002", msg)
                .with_help("valid indices are 0..n for the n rows or columns added so far"),
        )
    })
}

/// A primal value was read from a Solution without one (RS4007).
pub(crate) fn no_solution(status: &str) -> PyErr {
    let msg = format!("no primal solution is available (status {status})");
    Python::attach(|py| {
        raise(
            py,
            "NoSolutionError",
            Diagnostic::error("RS4007", msg)
                .with_help("check `Solution.ok` first, or call `Solution.raise_for_status()`"),
        )
    })
}

/// The error `Solution.raise_for_status()` raises for a non-optimal `status` (None when optimal).
pub(crate) fn status_err(status: &str, post: Vec<Diagnostic>) -> Option<PyErr> {
    let (class, d) = match status {
        "Optimal" => return None,
        "Infeasible" => (
            "InfeasibleError",
            Diagnostic::error("RS4004", "the model is infeasible")
                .with_help("`Model.explain()` returns an irreducible infeasible subsystem of rows and bounds")
                .with_help("`Model.feasrelax()` finds the cheapest relaxation that makes it feasible"),
        ),
        "Unbounded" | "InfeasibleOrUnbounded" => (
            "UnboundedError",
            Diagnostic::error(
                "RS4005",
                format!("the model is unbounded or infeasible (status {status})"),
            )
            .with_help("`Model.explain()` shows the unbounded ray or the infeasible subsystem")
            .with_help("bound the variables the ray moves along, or add the missing constraint"),
        ),
        s => (
            "LimitError",
            Diagnostic::error("RS4006", format!("the solve ended without an optimal solution (status {s})"))
                .with_help("raise `time_limit`, or read `Solution.x` if `ok` is true"),
        ),
    };
    Some(Python::attach(|py| raise_with(py, class, d, post)))
}

/// [`raise`], followed by more diagnostics (the post-solve ones) in the same document.
fn raise_with(py: Python<'_>, class: &str, mut d: Diagnostic, more: Vec<Diagnostic>) -> PyErr {
    if d.location.is_none()
        && let Some(l) = caller_location(py)
    {
        d.location = Some(l);
    }
    let mut all = vec![d.clone()];
    all.extend(more);
    let built = (|| -> PyResult<PyErr> {
        let cls = py.import("rusolve.errors")?.getattr(class)?;
        let e = cls.call_method1("_from_json", (d.message.clone(), diagnostics_json(&all, 1)))?;
        Ok(PyErr::from_value(e))
    })();
    built.unwrap_or_else(|e| e)
}

/// The CLI's post-solve diagnostics for an infeasible or unbounded `lp` (`adopt::post_solve`: the
/// irreducible infeasible set RS2010, the ray RS2011) within the CLI's bounded budget, with the
/// help lines naming Python calls instead of CLI commands (docs/design/diagnostics.md).
pub(crate) fn post_solve(lp: &HighsLp, status: &str, time: f64, quad: bool) -> Vec<Diagnostic> {
    let ms = match status {
        "Infeasible" => HighsModelStatus::kInfeasible,
        "Unbounded" => HighsModelStatus::kUnbounded,
        "InfeasibleOrUnbounded" => HighsModelStatus::kUnboundedOrInfeasible,
        _ => return Vec::new(),
    };
    let budget = rusolve_cli::adopt::budget_for(time, f64::INFINITY);
    let sol = rusolve_ekk::prelude::HighsSolution::default();
    let mut out = rusolve_cli::adopt::post_solve("the model", lp, ms, &sol, false, budget);
    // infeasibility does not depend on the objective, but a ray computed without the Hessian
    // is not a ray of a quadratic objective: drop RS2011 then
    if quad {
        out.retain(|d| d.code != "RS2011");
    }
    for d in &mut out {
        let code = d.code;
        for h in &mut d.help {
            if h.starts_with("`rusolve explain") {
                *h = if code == "RS2010" {
                    "`Model.iis()` verifies the set and shows the arithmetic".into()
                } else {
                    "`Model.explain()` verifies the ray and shows the objective rate".into()
                };
            } else if h.starts_with("`rusolve feasrelax") {
                *h = "`Model.feasrelax()` finds the cheapest relaxation that makes it feasible"
                    .into();
            } else if h.contains("`rusolve lint`") {
                *h = h.replace("`rusolve lint`", "`Model.lint()`");
            }
        }
        for n in &mut d.notes {
            *n = n
                .replace(
                    "`rusolve explain --mip-iis` searches",
                    "`Model.iis()` can search",
                )
                .replace(
                    "`rusolve explain the model` computes",
                    "`Model.iis()` computes",
                )
                .replace("`rusolve explain the model` has", "`Model.iis()` has");
        }
    }
    out
}
