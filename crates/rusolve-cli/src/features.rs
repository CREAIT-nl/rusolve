//! Instance features for an agent or tuner.
//!
//! * `rusolve classify <model> --format json` (schema `rusolve.classify/1`): the 44 static
//!   features of `rusolve_mip::classify` with a name, value and one-line description each.
//! * `rusolve probe <model> [--node-limit N] [--time-limit S] [--seed K] [--format json|text]`
//!   (schema `rusolve.probe/1`): a short run on the port that reports dynamic features.
//!
//! ## Where every probe number comes from
//!
//! Nothing here changes or feeds the search. The probe makes three ordinary `Highs::run` calls
//! on the port, each with the options of `rusolve mip` (serial, seed 0, `base_options()`):
//!
//! 1. the LP relaxation (integrality dropped): root LP objective, simplex iterations, wall time
//!    and the integer columns that are fractional in its solution. The port keeps no pre-cut
//!    root LP solution (`HighsMipSolverData::firstlpsol` / `rootlpsol` are never assigned), so
//!    the relaxation is solved on its own;
//! 2. the MIP with `mip_max_nodes = 1` (0 would stop before the root LP): stops after the root node, so its wall time is the root
//!    phase (MIP presolve, root LP, cut rounds, root heuristics) and its dual bound is the bound
//!    after the root cut loop (HiGHS prints it as BestBound at the end of the root);
//! 3. the MIP with `mip_max_nodes = N`: tree statistics.
//!
//! Runs 2 and 3 are deterministic for the same seed, so run 3 repeats run 2's root and tree-phase
//! quantities are differences (run 3 minus run 2).
//!
//! MIP counters come from `Highs::mip_probe_` (rusolve-ekk `mip/mip_probe.rs`), a read-only
//! copy of what `HighsMipSolverData` holds when the solve returns, taken only when
//! `Highs::capture_mip_probe_` is set. `HighsInfo` supplies the bounds, gap and node count.
//!
//! ## Not recorded by the port (reported as `null`, listed under `unavailable`)
//!
//! * time and source of the first incumbent: `saveReportMipSolution` keeps neither
//!   (a callback would give the time, but registering one changes which code paths run, which
//!   the probe must not do);
//! * the maximum depth reached: the search keeps only the current depth. The deepest *open*
//!   node left in the queue is reported instead, a lower bound on the depth reached;
//! * root cut rounds and cuts added: `HighsMipSolverData::evaluateRootNode` runs the root
//!   separation loop with a local `HighsSeparation` and a local round counter (part3.rs), both
//!   dropped on return. The cuts alive in the global pool and the LP iterations spent in
//!   separation (`sepa_lp_iterations`) are kept and reported; the per-separator counters of
//!   the worker's separation object cover the tree phase only.
//!
//! `features_version` versions the probe feature list (names, order, meaning, normalisation);
//! bump it on any change. Exit codes follow the exit-code table: 0 ran, 2 usage, 3 input, 4 internal.

use crate::envelope::{SOLVER_ENGINE, SOLVER_NAME, Solver, profile_name, sha256_of_file};
use crate::mip_highs::read_mip;
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::{Value, json};
use std::process::ExitCode;
use std::time::Instant;

pub const CLASSIFY_SCHEMA_VERSION: &str = "rusolve.classify/1";
pub const PROBE_SCHEMA_VERSION: &str = "rusolve.probe/1";
/// Version of the probe feature list (`features` in `rusolve.probe/1`).
pub const PROBE_FEATURES_VERSION: u32 = 1;

fn num(v: f64) -> Value {
    if v.is_finite() { json!(v) } else { Value::Null }
}

fn opt(v: Option<f64>) -> Value {
    v.map_or(Value::Null, num)
}

fn solver_json() -> Value {
    serde_json::to_value(Solver {
        name: SOLVER_NAME,
        version: env!("CARGO_PKG_VERSION"),
        git: option_env!("RUSOLVE_GIT"),
        profile: profile_name(),
        engine: SOLVER_ENGINE,
    })
    .expect("solver serialises")
}

fn usage(json: bool, msg: &str) -> ExitCode {
    crate::diag::fail_usage_text(json, "probe", msg.strip_prefix("probe: ").unwrap_or(msg))
}

/// `rusolve classify <model> --format json`: see the module doc. Exit codes: 0, 3 (input).
pub fn classify_json(path: &str, miplib_strict: bool) -> ExitCode {
    let model = match rusolve_io::read_model_file(path) {
        Ok(m) => m,
        Err(e) => return crate::diag::fail_read(true, path, &e.to_string()),
    };
    if let Err(e) = model.gc.reject("rusolve classify") {
        return crate::diag::fail_msg(true, "RS0012", 3, e.to_string());
    }
    let opts = rusolve_mip::classify::ClassifyOptions { miplib_strict };
    let problem = rusolve_mip::MipProblem::from_model(&model);
    let cl = rusolve_mip::classify::classify(&problem, opts);
    let values = cl.features();
    let features: Vec<Value> = rusolve_mip::classify::FEATURE_NAMES
        .iter()
        .zip(rusolve_mip::classify::FEATURE_DESCRIPTIONS.iter())
        .zip(values.iter())
        .map(|((n, d), v)| json!({"name": n, "value": num(*v), "description": d}))
        .collect();
    // The detail block is the existing `--json` object, unchanged.
    let detail: Value = serde_json::from_str(&cl.to_json()).expect("classify json parses");
    let out = json!({
        "schema_version": CLASSIFY_SCHEMA_VERSION,
        "features_version": rusolve_mip::classify::FEATURE_VERSION,
        "solver": solver_json(),
        "input": {
            "path": path,
            "sha256": sha256_of_file(path).unwrap_or_default(),
            "rows": cl.num_rows,
            "columns": cl.num_cols,
            "nonzeros": cl.nnz,
        },
        "class": crate::route::class_of_file(path),
        "options": {"miplib_strict": miplib_strict},
        "features": features,
        "detail": detail,
    });
    println!("{out}");
    ExitCode::SUCCESS
}

/// What one `Highs::run` of the probe left behind.
struct RunOut {
    status: &'static str,
    wall_s: f64,
    objective: Option<f64>,
    has_solution: bool,
    dual_bound: f64,
    gap: f64,
    nodes: i64,
    lp_iterations: i64,
    x: Vec<f64>,
    snap: Option<rusolve_ekk::mip::mip_probe::MipProbeSnapshot>,
}

fn run_once(
    lp: &HighsLp,
    relax: bool,
    node_limit: Option<HighsInt>,
    time_left: f64,
    seed: HighsInt,
) -> Result<RunOut, String> {
    // The options `rusolve mip` uses (mip_highs.rs trace_options).
    let mut options = HighsOptions {
        output_flag: false,
        threads: 1,
        random_seed: seed,
        time_limit: time_left,
        // the relaxation run stands for the MIP root LP, which the lever never touches
        lp_solver_auto: false,
        ..crate::base_options()
    };
    if let Some(n) = node_limit {
        options.mip_max_nodes = n;
    }
    options.sync_log_options();
    let mut h = Highs::with_options(options);
    h.capture_mip_probe_ = !relax;
    let mut model = lp.clone();
    if relax {
        model.integrality_.clear();
    }
    let t = Instant::now();
    if h.pass_model(model) == HighsStatus::kError || h.run() == HighsStatus::kError {
        return Err("the solve failed".to_string());
    }
    let wall_s = t.elapsed().as_secs_f64();
    let info = h.get_info();
    let has_solution = h.solution_.value_valid;
    Ok(RunOut {
        status: crate::envelope::result_status(h.get_model_status()),
        wall_s,
        objective: has_solution.then_some(info.objective_function_value),
        has_solution,
        dual_bound: info.mip_dual_bound,
        gap: info.mip_gap,
        nodes: info.mip_node_count,
        lp_iterations: i64::from(info.simplex_iteration_count),
        x: h.solution_.col_value.clone(),
        snap: h.mip_probe_.take(),
    })
}

fn frac(part: f64, whole: f64) -> Option<f64> {
    (whole > 0.0).then(|| part / whole)
}

/// `rusolve probe`: see the module doc.
pub fn run_probe(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let Some(path) = args.first() else {
        return usage(wj, "probe: missing model path");
    };
    let mut json_out = false;
    let mut node_limit: HighsInt = 50;
    let mut time_limit = 60.0_f64;
    let mut seed: HighsInt = 0;
    let mut it = args[1..].iter();
    while let Some(flag) = it.next() {
        match flag.as_str() {
            "--format" => match it.next().map(String::as_str) {
                Some("json") => json_out = true,
                Some("text") => json_out = false,
                _ => return usage(wj, "probe: --format needs json or text"),
            },
            "--node-limit" => match it.next().and_then(|v| v.parse::<HighsInt>().ok()) {
                Some(n) if n >= 1 => node_limit = n,
                _ => return usage(wj, "probe: --node-limit needs an integer >= 1"),
            },
            "--time-limit" => match it.next().and_then(|v| v.parse::<f64>().ok()) {
                Some(t) if t > 0.0 => time_limit = t,
                _ => return usage(wj, "probe: --time-limit needs a positive number of seconds"),
            },
            "--seed" => match it.next().and_then(|v| v.parse::<HighsInt>().ok()) {
                Some(s) => seed = s,
                None => return usage(wj, "probe: --seed needs an integer"),
            },
            other if !other.starts_with("--") => {
                return usage(wj, &format!("probe: unexpected extra argument '{other}'"));
            }
            other => return usage(wj, &format!("probe: unsupported flag '{other}'")),
        }
    }
    let lp = match read_mip(path) {
        Ok(lp) => lp,
        Err(e) => return crate::diag::fail_read_str(json_out, path, &e),
    };
    let rows = i64::from(lp.num_row_);
    let cols = i64::from(lp.num_col_);
    let nnz = i64::from(lp.a_matrix_.start_[lp.num_col_ as usize]);
    let int_cols: Vec<usize> = lp
        .integrality_
        .iter()
        .enumerate()
        .filter(|(_, t)| **t == HighsVarType::kInteger)
        .map(|(j, _)| j)
        .collect();
    let is_mip = !int_cols.is_empty();
    let start = Instant::now();
    let left = || (time_limit - start.elapsed().as_secs_f64()).max(0.001);

    let fail = |what: &str, e: String| {
        crate::diag::fail_msg(json_out, "RS0090", 4, format!("rusolve probe: {what}: {e}"))
    };
    // Run 1: the LP relaxation (for a pure LP this is the whole solve).
    let relax = match run_once(&lp, true, None, left(), seed) {
        Ok(r) => r,
        Err(e) => return fail("LP relaxation", e),
    };
    let feastol = HighsOptions::default().mip_feasibility_tolerance;
    let fractional = int_cols
        .iter()
        .filter(|&&j| {
            relax
                .x
                .get(j)
                .is_some_and(|v| (v - v.round()).abs() > feastol)
        })
        .count();
    // Runs 2 and 3.
    let (root, tree) = if is_mip {
        let root = match run_once(&lp, false, Some(1), left(), seed) {
            Ok(r) => r,
            Err(e) => return fail("root-only MIP run", e),
        };
        let tree = match run_once(&lp, false, Some(node_limit), left(), seed) {
            Ok(r) => r,
            Err(e) => return fail("MIP probe run", e),
        };
        (Some(root), Some(tree))
    } else {
        (None, None)
    };
    let (rs, ts) = (
        root.as_ref().and_then(|r| r.snap.as_ref()),
        tree.as_ref().and_then(|r| r.snap.as_ref()),
    );

    // Presolve: the model the MIP solver ended on, against the original.
    let presolve = ts.map(|s| {
        (
            rows - s.model_num_row,
            cols - s.model_num_col,
            nnz - s.model_num_nz,
        )
    });
    // Objective values are in the user's sense; `sgn` maps them to minimisation.
    let sgn = if lp.sense_ == ObjSense::kMaximize {
        -1.0
    } else {
        1.0
    };
    let lp_obj = relax.objective.filter(|_| relax.status == "Optimal");
    let root_bound = root
        .as_ref()
        .map(|r| r.dual_bound)
        .filter(|b| b.is_finite());
    let primal = tree
        .as_ref()
        .and_then(|t| t.objective)
        .or_else(|| root.as_ref().and_then(|r| r.objective));
    let gap_closed = match (lp_obj, root_bound, primal) {
        (Some(l), Some(b), Some(p)) => {
            let (l, b, p) = (sgn * l, sgn * b, sgn * p);
            (p - l > 1e-9).then(|| (b - l) / (p - l))
        }
        _ => None,
    };
    let tree_nodes = match (&root, &tree) {
        (Some(r), Some(t)) => Some(t.nodes - r.nodes),
        _ => None,
    };
    let tree_time = match (&root, &tree) {
        (Some(r), Some(t)) => Some(t.wall_s - r.wall_s),
        _ => None,
    };
    let nodes_per_s = match (tree_nodes, tree_time) {
        (Some(n), Some(t)) if n > 0 && t > 1e-3 => Some(n as f64 / t),
        _ => None,
    };
    let lp_iters_per_node = match (rs, ts, tree_nodes) {
        (Some(r), Some(t), Some(n)) if n > 0 => {
            Some((t.total_lp_iterations - r.total_lp_iterations) as f64 / n as f64)
        }
        _ => None,
    };
    let depth_open = ts.and_then(|s| s.open_nodes_max_depth);

    // The compact feature vector for a tuner.
    let frac_removed =
        |removed: Option<i64>, whole: i64| removed.and_then(|r| frac(r as f64, whole as f64));
    let features: Vec<(&str, Option<f64>, &str)> = vec![
        (
            "presolve.rows_removed_frac",
            frac_removed(presolve.map(|p| p.0), rows),
            "share of rows the MIP presolve removed (after the last restart)",
        ),
        (
            "presolve.cols_removed_frac",
            frac_removed(presolve.map(|p| p.1), cols),
            "share of columns the MIP presolve removed",
        ),
        (
            "presolve.nnz_removed_frac",
            frac_removed(presolve.map(|p| p.2), nnz),
            "share of nonzeros the MIP presolve removed",
        ),
        (
            "root_lp.iterations",
            Some(relax.lp_iterations as f64),
            "simplex iterations of the LP relaxation solve",
        ),
        (
            "root_lp.time_s",
            Some(relax.wall_s),
            "wall seconds of the LP relaxation solve",
        ),
        (
            "root_lp.fractional_frac",
            is_mip.then(|| fractional as f64 / int_cols.len() as f64),
            "share of integer columns fractional in the LP relaxation solution",
        ),
        (
            "root.time_s",
            root.as_ref().map(|r| r.wall_s),
            "wall seconds of the root phase (MIP presolve, root LP, cuts, heuristics)",
        ),
        (
            "root.gap_closed_by_cuts",
            gap_closed,
            "(root bound after cuts - LP objective) / (best incumbent - LP objective)",
        ),
        (
            "root.cuts_in_pool",
            rs.map(|s| s.cutpool_cuts as f64),
            "cuts alive in the global cut pool after the root node (root-only run)",
        ),
        (
            "root.sepa_lp_iterations",
            rs.map(|s| s.sepa_lp_iterations as f64),
            "LP iterations spent in root separation rounds (root-only run)",
        ),
        (
            "tree.nodes_per_second",
            nodes_per_s,
            "nodes explored after the root per second of tree-phase wall time",
        ),
        (
            "tree.lp_iterations_per_node",
            lp_iters_per_node,
            "LP iterations after the root per node explored after the root",
        ),
        (
            "tree.open_nodes_max_depth",
            depth_open.map(|d| d as f64),
            "depth of the deepest open node at the end (lower bound on depth reached)",
        ),
        (
            "bounds.gap_final",
            tree.as_ref().map(|t| t.gap).filter(|g| g.is_finite()),
            "relative MIP gap when the probe stopped",
        ),
        (
            "incumbent.found",
            tree.as_ref().map(|t| f64::from(u8::from(t.has_solution))),
            "1 when the probe found an incumbent within its node limit",
        ),
    ];

    let sepa_json = |s: Option<&rusolve_ekk::mip::mip_probe::MipProbeSnapshot>| -> Value {
        s.map_or(Value::Null, |s| {
            Value::Array(
                s.separators
                    .iter()
                    .map(|x| json!({"name": x.name, "calls": x.num_calls, "cuts_found": x.num_cuts_found}))
                    .collect(),
            )
        })
    };
    let out = json!({
        "schema_version": PROBE_SCHEMA_VERSION,
        "features_version": PROBE_FEATURES_VERSION,
        "solver": solver_json(),
        "input": {
            "path": path,
            "sha256": sha256_of_file(path).unwrap_or_default(),
            "rows": rows,
            "columns": cols,
            "nonzeros": nnz,
            "integer_columns": int_cols.len(),
        },
        "budget": {"node_limit": node_limit, "time_limit_s": time_limit, "seed": seed, "threads": 1},
        "is_mip": is_mip,
        "status": tree.as_ref().map_or(relax.status, |t| t.status),
        "wall_s": start.elapsed().as_secs_f64(),
        "presolve": ts.map_or(Value::Null, |s| json!({
            "rows_before": rows, "cols_before": cols, "nonzeros_before": nnz,
            "rows_after": s.model_num_row, "cols_after": s.model_num_col, "nonzeros_after": s.model_num_nz,
            "rows_removed": rows - s.model_num_row, "cols_removed": cols - s.model_num_col,
            "nonzeros_removed": nnz - s.model_num_nz,
            "integer_columns_after": s.model_num_integer,
            "restarts": s.num_restarts,
        })),
        "root_lp": {
            "source": "LP relaxation solved on its own (no pre-cut root LP is kept by the port)",
            "status": relax.status,
            "objective": opt(lp_obj),
            "iterations": relax.lp_iterations,
            "time_s": relax.wall_s,
            "integer_columns": int_cols.len(),
            "fractional_columns": if is_mip { json!(fractional) } else { Value::Null },
            "fractional_fraction": opt(is_mip.then(|| fractional as f64 / int_cols.len() as f64)),
        },
        "root": root.as_ref().map_or(Value::Null, |r| json!({
            "status": r.status,
            "time_s": r.wall_s,
            "lp_iterations_first_root_lp": rs.map(|s| s.first_root_lp_iterations),
            "lp_iterations_total": rs.map(|s| s.total_lp_iterations),
            "dual_bound_after_cuts": num(r.dual_bound),
            "gap_closed_by_cuts": opt(gap_closed),
            "cutpool_cuts_alive": rs.map(|s| s.cutpool_cuts),
            "lp_iterations_separation": rs.map(|s| s.sepa_lp_iterations),
        })),
        "tree": tree.as_ref().map_or(Value::Null, |t| json!({
            "nodes": t.nodes,
            "nodes_after_root": tree_nodes,
            "leaves": ts.map(|s| s.num_leaves),
            "open_nodes": ts.map(|s| s.open_nodes),
            "open_nodes_max_depth": depth_open,
            "tree_phase_time_s": opt(tree_time.filter(|t| *t > 0.0)),
            "nodes_per_second": opt(nodes_per_s),
            "lp_iterations_total": ts.map(|s| s.total_lp_iterations),
            "lp_iterations_heuristics": ts.map(|s| s.heuristic_lp_iterations),
            "lp_iterations_separation": ts.map(|s| s.sepa_lp_iterations),
            "lp_iterations_strong_branching": ts.map(|s| s.sb_lp_iterations),
            "lp_iterations_per_node": opt(lp_iters_per_node),
            "separators": sepa_json(ts),
        })),
        "bounds": tree.as_ref().map_or(Value::Null, |t| json!({
            "primal_bound": opt(t.objective),
            "dual_bound": num(t.dual_bound),
            "gap": num(t.gap),
            "root_lp_objective": opt(lp_obj),
            "root_bound_after_cuts": opt(root_bound),
        })),
        "incumbent": {
            "found": tree.as_ref().map(|t| t.has_solution),
            "improving_solutions": ts.map(|s| s.num_improving_sols),
            "first_time_s": Value::Null,
            "first_source": Value::Null,
        },
        "unavailable": [
            {"feature": "incumbent.first_time_s", "reason": "the port records no time for the first incumbent; a callback would change code paths"},
            {"feature": "incumbent.first_source", "reason": "the port records no source for the first incumbent"},
            {"feature": "tree.depth_max", "reason": "the search keeps only its current depth; tree.open_nodes_max_depth is a lower bound"},
            {"feature": "root.cut_rounds", "reason": "the root separation loop keeps its round count and per-separator counters in locals that are dropped; the tree-phase separator counters are reported under tree.separators"},
            {"feature": "root.cuts_added", "reason": "only the cuts alive in the global pool are kept (root.cutpool_cuts_alive)"},
        ],
        "features": features.iter().map(|(n, v, d)| json!({"name": n, "value": opt(*v), "description": d})).collect::<Vec<_>>(),
    });
    if json_out {
        println!("{out}");
    } else {
        print_text(path, &out, &features);
    }
    ExitCode::SUCCESS
}

fn print_text(path: &str, out: &Value, features: &[(&str, Option<f64>, &str)]) {
    println!("probe {path}");
    println!(
        "  status {}  is_mip {}  wall {:.3}s  node limit {}",
        out["status"].as_str().unwrap_or("?"),
        out["is_mip"],
        out["wall_s"].as_f64().unwrap_or(0.0),
        out["budget"]["node_limit"]
    );
    println!("  features (version {PROBE_FEATURES_VERSION}):");
    for (n, v, _) in features {
        match v {
            Some(v) => println!("    {n:<32} {v}"),
            None => println!("    {n:<32} n/a"),
        }
    }
}
