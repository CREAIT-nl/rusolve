//! `mip_engine=cpsat`: a MIP solved by Google OR-Tools CP-SAT (Apache-2.0) instead of the
//! HiGHS port. Opt-in only; the default engine is `highs`. OR-Tools is loaded at run time
//! (`rusolve_cpsat::runtime`, search order in its module docs); when it is not installed the
//! selection is a typed error (RS0012, exit 3) that says where it looked and how to install it
//! (`docs/design/cpsat-engine.md`).
//!
//! The engine takes what CP-SAT takes exactly: every column integer with finite bounds,
//! integral coefficients and right-hand sides (`rusolve_cpsat::translate` names the first column
//! or row that is not). The answer is never taken on trust:
//!
//! | CP-SAT | result |
//! |---|---|
//! | `OPTIMAL` | `Optimal`, after the point passes [`crate::check`] on the original model and its objective equals CP-SAT's |
//! | `FEASIBLE` | `TimeLimit` with the incumbent (a limit stopped the search), same check |
//! | `INFEASIBLE` | `Infeasible` |
//! | `UNKNOWN` | `TimeLimit`, no incumbent |
//! | `MODEL_INVALID` | error, with CP-SAT's message |
//!
//! A point the checker rejects, or whose objective disagrees with CP-SAT's, is an internal error
//! (exit 4), never an answer. [`run`] fills the `Highs` object (status, solution, info) so the
//! rest of `rusolve mip` and the MCP `solve` tool build their usual result from it, and returns
//! the `cpsat` object of the envelope.

use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::{HighsModelStatus, HighsSolution, ObjSense, kSolutionStatusFeasible};
use serde_json::{Value, json};

/// The `engine` field of an envelope this module produced.
pub const ENGINE: &str = "cp-sat";

/// Why a CP-SAT run did not produce a result.
#[derive(Debug, Clone, PartialEq)]
pub enum Failure {
    /// The OR-Tools library is not installed or not loadable (exit 3); the message lists the
    /// paths searched and the install routes.
    Unavailable(String),
    /// The model is not one CP-SAT can take exactly, or the option value is wrong (exit 3).
    Refused(String),
    /// CP-SAT rejected the model, returned something that does not decode, or returned a point
    /// the checker rejects (exit 4).
    Internal(String),
}

impl Failure {
    pub fn message(&self) -> String {
        match self {
            Failure::Unavailable(m) => format!("mip_engine=cpsat: {m}"),
            Failure::Refused(m) | Failure::Internal(m) => m.clone(),
        }
    }

    /// The exit code: 3 for input the engine cannot honour, 4 for an internal failure.
    pub fn exit(&self) -> u8 {
        match self {
            Failure::Unavailable(_) | Failure::Refused(_) => 3,
            Failure::Internal(_) => 4,
        }
    }
}

/// Whether `mip_engine` selects CP-SAT; an error for a value that is neither `highs` nor `cpsat`.
pub fn selected(highs: &Highs) -> Result<bool, String> {
    match highs.options_.mip_engine.as_str() {
        "highs" => Ok(false),
        "cpsat" => Ok(true),
        other => Err(format!(
            "mip_engine={other}: the engines are \"highs\" (the port, default) and \"cpsat\" (Google OR-Tools CP-SAT)"
        )),
    }
}

/// The CP-SAT parameters of `o`: time limit, workers (`threads`), seed and the two gaps.
pub fn params_of(o: &HighsOptions) -> rusolve_cpsat::CpParams {
    rusolve_cpsat::CpParams {
        time_limit: o.time_limit.is_finite().then_some(o.time_limit.max(0.0)),
        workers: o.threads.max(1),
        seed: o.random_seed,
        rel_gap: o.mip_rel_gap.max(0.0),
        abs_gap: o.mip_abs_gap.max(0.0),
        log: false,
    }
}

/// The `cpsat` object of the envelope for `out`, solved with `params`.
pub fn report_of(
    out: &rusolve_cpsat::CpOutcome,
    params: &rusolve_cpsat::CpParams,
    bound: Option<f64>,
    objective: Option<f64>,
    check: Value,
    note: Option<&str>,
) -> Value {
    json!({
        "status": out.status.name(),
        "bound": bound,
        "objective": objective,
        "wall_time_s": out.wall_time_s,
        "workers": out.workers,
        "ortools_version": out.ortools_version,
        "ortools_library": out.ortools_library.display().to_string(),
        "objective_kind": out.objective_kind.name(),
        "conflicts": out.conflicts,
        "branches": out.branches,
        "aux_variables": out.aux_variables,
        "time_limit_s": params.time_limit,
        "rel_gap_limit": params.rel_gap,
        "abs_gap_limit": params.abs_gap,
        "seed": params.seed,
        "checked": check,
        "note": note,
    })
}

/// The `cpsat` object of `rusolve version --format json` and of the MCP `server_info`: whether the
/// library loads, from where, which version, and the places searched. Loads the library.
pub fn availability_json() -> Value {
    let a = rusolve_cpsat::runtime::availability();
    json!({
        "available": a.available,
        "library": a.library.map(|p| p.display().to_string()),
        "version": a.version,
        "searched": a.searched,
        "env": rusolve_cpsat::runtime::ENV_LIB,
        "install": ["install.sh --with-cpsat", "pip install 'rusolve[cpsat]'"],
        "error": a.error,
    })
}

/// `true` when [`rusolve_cpsat::runtime`] finds the library (loads it).
pub fn available() -> bool {
    rusolve_cpsat::runtime::library().is_ok()
}

/// Solves `highs.model_lp_` with CP-SAT and fills `highs`' model status, solution and info.
/// `model_path` only labels the checker report. Returns the `cpsat` object of the envelope.
pub fn run(highs: &mut Highs, model_path: &str) -> Result<Value, Failure> {
    run_cp(highs, None, model_path)
}

/// [`run`] for a model with CP objects: they go to CP-SAT next to the rows, and the point is also
/// checked against them (`CpModel::check`) before anything is reported.
pub fn run_cp(
    highs: &mut Highs,
    cp: Option<&rusolve_io::CpModel>,
    model_path: &str,
) -> Result<Value, Failure> {
    use rusolve_cpsat::{CpError, CpStatus};

    let o = &highs.options_;
    let params = params_of(o);
    let tol = crate::check::Tolerances {
        primal_feasibility_tolerance: o.primal_feasibility_tolerance,
        mip_feasibility_tolerance: o.mip_feasibility_tolerance,
        applied: o.mip_feasibility_tolerance,
    };
    let lp = &highs.model_lp_;
    let maximize = lp.sense_ == ObjSense::kMaximize;
    let out = match rusolve_cpsat::solve(lp, cp, &params) {
        Ok(out) => out,
        Err(CpError::Unavailable(e)) => return Err(Failure::Unavailable(e.to_string())),
        Err(CpError::Refused(r)) => return Err(Failure::Refused(format!("mip_engine=cpsat: {r}"))),
        Err(e @ CpError::ModelInvalid(_)) => {
            return Err(Failure::Internal(format!("mip_engine=cpsat: {e}")));
        }
        Err(e) => return Err(Failure::Internal(format!("mip_engine=cpsat: {e}"))),
    };

    let mut status = match out.status {
        CpStatus::Optimal => HighsModelStatus::kOptimal,
        CpStatus::Infeasible => HighsModelStatus::kInfeasible,
        CpStatus::Feasible | CpStatus::Unknown => HighsModelStatus::kTimeLimit,
    };
    let mut check = Value::Null;
    let mut point: Option<(Vec<f64>, f64)> = None;
    if let Some(x) = &out.x {
        // our checker, on the original model, before anything is reported
        let rep = crate::check::check_solution(lp, x, tol, 5, model_path, "cp-sat")
            .map_err(|e| Failure::Internal(format!("mip_engine=cpsat: checker: {}", e.0)))?;
        check = json!({
            "verdict": rep.verdict,
            "violations": rep.summary.total_violations,
            "objective": rep.objective.value,
        });
        if rep.summary.total_violations != 0 {
            let first = rep.findings.first().map_or(String::new(), |f| {
                format!(" (first: {} {} by {:e})", f.kind, f.name, f.amount)
            });
            return Err(Failure::Internal(format!(
                "mip_engine=cpsat: CP-SAT returned a point that violates {} item(s) of the model{first}; \
                 the point is discarded, this is a bug",
                rep.summary.total_violations
            )));
        }
        if let Some(cp) = cp {
            let tol = rusolve_io::CpTol {
                integrality: highs.options_.mip_feasibility_tolerance,
            };
            let bad = cp.check(x, &lp.col_names_, tol);
            if let Some(first) = bad.first() {
                return Err(Failure::Internal(format!(
                    "mip_engine=cpsat: CP-SAT returned a point that violates {} scheduling or \
                     combinatorial object(s) of the model (first: {} {}); the point is discarded, \
                     this is a bug",
                    bad.len(),
                    first.kind,
                    first.name
                )));
            }
            check["cp_violations"] = json!(0);
        }
        let ours = rep.objective.value;
        if let Some(theirs) = out.objective
            && (ours - theirs).abs() > 1e-6 * (1.0 + ours.abs())
        {
            return Err(Failure::Internal(format!(
                "mip_engine=cpsat: CP-SAT's objective {theirs} differs from the objective {ours} of its point; \
                 the point is discarded, this is a bug"
            )));
        }
        point = Some((x.clone(), ours));
    } else if matches!(out.status, CpStatus::Optimal | CpStatus::Feasible) {
        return Err(Failure::Internal(
            "mip_engine=cpsat: CP-SAT reported a solution but returned no point".into(),
        ));
    }

    // a bound that passed the point would make the answer inconsistent: say so, never "optimal"
    let mut bound = out.bound;
    let mut bound_note = None;
    if let (Some((_, obj)), Some(b)) = (&point, bound) {
        let past = if maximize { b < *obj } else { b > *obj };
        if past && (b - obj).abs() > 1e-6 * (1.0 + obj.abs()) {
            bound_note = Some("the bound passed the checked point; status downgraded");
            if status == HighsModelStatus::kOptimal {
                status = HighsModelStatus::kUnknown;
            }
        } else if past {
            bound = Some(*obj);
        }
    }

    let inf = if maximize {
        f64::INFINITY
    } else {
        f64::NEG_INFINITY
    };
    let n = lp.num_col_ as usize;
    let mut row_value = vec![0.0; lp.num_row_ as usize];
    if let Some((x, _)) = &point {
        let a = &lp.a_matrix_;
        for j in 0..n {
            for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                row_value[a.index_[p] as usize] += a.value_[p] * x[j];
            }
        }
    }
    let objective = point.as_ref().map(|(_, o)| *o);
    let dual_bound = match (status, objective) {
        (HighsModelStatus::kOptimal, Some(obj)) => bound.unwrap_or(obj),
        _ => bound.unwrap_or(inf),
    };
    fill(
        highs,
        status,
        point.map(|(x, _)| (x, row_value)),
        objective,
        dual_bound,
    );

    Ok(report_of(
        &out, &params, bound, objective, check, bound_note,
    ))
}

/// The state the rest of `rusolve mip` and the MCP `solve` tool read their result from.
pub fn fill(
    highs: &mut Highs,
    status: HighsModelStatus,
    point: Option<(Vec<f64>, Vec<f64>)>,
    objective: Option<f64>,
    dual_bound: f64,
) {
    highs.model_status_ = status;
    let has_point = point.is_some();
    highs.solution_ = match point {
        Some((col_value, row_value)) => HighsSolution {
            value_valid: true,
            col_value,
            row_value,
            ..Default::default()
        },
        None => HighsSolution::default(),
    };
    let info = &mut highs.info_;
    info.valid = true;
    info.objective_function_value = objective.unwrap_or(0.0);
    info.mip_dual_bound = dual_bound;
    // the gap as HiGHS computes it: |ub - lb| / |ub|, infinite without an incumbent
    info.mip_gap = match objective {
        None => f64::INFINITY,
        Some(obj) if obj != 0.0 => (obj - dual_bound).abs() / obj.abs(),
        Some(_) if dual_bound == 0.0 => 0.0,
        Some(_) => f64::INFINITY,
    };
    // CP-SAT has no tree and no simplex: the counters it does have are in the `cpsat` object
    info.mip_node_count = 0;
    info.simplex_iteration_count = 0;
    info.primal_solution_status = if has_point {
        kSolutionStatusFeasible
    } else {
        0
    };
}
