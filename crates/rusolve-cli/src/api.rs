//! The diagnostic cores as a library surface, for the Python bindings
//! (`crates/rusolve-py`) and any other in-process caller.
//!
//! Every function here is a thin adapter over a core the CLI and the MCP server already share
//! (`check::check_solution`, `explain::explain_core`, `lint::lint_core`, `ranging::ranging_core`,
//! `feasrelax::feasrelax_core`): the model arrives as a loaded `HighsLp` instead of a file path,
//! the answer leaves as the same `rusolve.<command>/1` JSON the CLI prints, so a Python caller and
//! a shell caller see identical objects. Nothing here changes a solve. Errors are the message the
//! CLI would print after its `rusolve <command>:` prefix.
//!
//! The envelope path follows `mcp_tools.rs::solve` (a model without a file digests a canonical
//! dump).

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::mip_stats::MipStats;
use rusolve_ekk::prelude::*;
use serde_json::Value;

use crate::check;
pub use crate::envelope::Envelope;
use crate::envelope::{self, ModelCounts, Solved};

/// Display name of a model that has no file, used for `input.path` / `model.path`.
pub const DEFAULT_NAME: &str = "<python model>";

/// Values by column name, in the model's column order. `pairs` may name a subset; the others
/// are an error (the checker needs a complete assignment).
pub fn columns_by_name(lp: &HighsLp, pairs: &[(String, f64)]) -> Result<Vec<f64>, String> {
    check::columns_by_name(lp, pairs).map_err(|e| e.0)
}

/// `rusolve.check/1` for a complete assignment `values` (original column order).
pub fn check(
    lp: &HighsLp,
    values: &[f64],
    primal_tolerance: Option<f64>,
    mip_tolerance: Option<f64>,
    max_findings: Option<usize>,
    name: &str,
) -> Result<Value, String> {
    let defaults = HighsOptions::default();
    let primal = primal_tolerance.unwrap_or(defaults.primal_feasibility_tolerance);
    let mip = mip_tolerance.unwrap_or(defaults.mip_feasibility_tolerance);
    for (f, v) in [
        ("primal_feasibility_tolerance", primal),
        ("mip_feasibility_tolerance", mip),
    ] {
        if !(v.is_finite() && v >= 0.0) {
            return Err(format!("{f}: expected a non-negative number, got {v}"));
        }
    }
    let tol = check::Tolerances {
        primal_feasibility_tolerance: primal,
        mip_feasibility_tolerance: mip,
        applied: if lp.is_mip() { mip } else { primal },
    };
    let report = check::check_solution(
        lp,
        values,
        tol,
        max_findings.unwrap_or(check::DEFAULT_MAX_FINDINGS),
        name,
        "<inline>",
    )
    .map_err(|e| e.0)?;
    serde_json::to_value(report).map_err(|e| format!("internal error: {e}"))
}

/// `rusolve.explain/1`: the IIS of an infeasible model, the ray of an unbounded one.
pub fn explain(
    lp: &HighsLp,
    name: &str,
    verify: bool,
    max_members: usize,
) -> Result<Value, String> {
    let a = crate::explain::core_args(name, verify, max_members);
    crate::explain::explain_core(lp.clone(), &a)
        .map(|(v, _)| v)
        .map_err(|e| e.msg)
}

/// `rusolve.diagnose/1`: formulation error or data error, for an infeasible model (the
/// `rusolve diagnose` verdict; a MIP is diagnosed through its LP relaxation). `time_limit` bounds
/// the minimum-count relaxation (seconds, the CLI default is 30).
pub fn diagnose(lp: &HighsLp, name: &str, time_limit: f64) -> Result<Value, String> {
    crate::diagnose::diagnose_core(lp.clone(), name, time_limit)
        .map(|(v, _)| v)
        .map_err(|()| format!("the solver failed on {name}"))
}

/// `rusolve.lint/1`.
pub fn lint(
    lp: &HighsLp,
    name: &str,
    rules: &[String],
    min_severity: Option<&str>,
    max_items: Option<usize>,
    with_solve: bool,
) -> Result<Value, String> {
    let mut argv = vec![name.to_string(), "--format".into(), "json".into()];
    for r in rules {
        argv.push("--rule".into());
        argv.push(r.clone());
    }
    if let Some(s) = min_severity {
        argv.push("--min-severity".into());
        argv.push(s.to_string());
    }
    if let Some(n) = max_items {
        argv.push("--max-items".into());
        argv.push(n.to_string());
    }
    if with_solve {
        argv.push("--with-solve".into());
    }
    let a = crate::lint::parse_args(&argv)?;
    crate::lint::lint_core(lp.clone(), &a)
        .map(|(v, _)| v)
        .map_err(|e| e.msg)
}

/// `rusolve.ranging/1`: the sensitivity ranges of the optimal solution (solves the model).
/// `rows` and `cols` are name patterns (empty: all).
pub fn ranging(
    lp: &HighsLp,
    name: &str,
    rows: &[String],
    cols: &[String],
) -> Result<Value, String> {
    crate::ranging::ranging_core(lp, name, rows, cols)
        .map(|(v, _)| v)
        .map_err(|e| e.msg)
}

/// `rusolve.feasrelax/1`. `flags` are the CLI words after the model (`--mode min-count`,
/// `--global-lower-penalty 2`, ...); `penalties` is the `--penalties` document as JSON text.
pub fn feasrelax(
    lp: &HighsLp,
    name: &str,
    flags: &[String],
    penalties: Option<&str>,
) -> Result<Value, String> {
    let file: Option<Value> = match penalties {
        Some(t) => Some(serde_json::from_str(t).map_err(|e| format!("penalties: {e}"))?),
        None => None,
    };
    let a = crate::feasrelax::core_args(name, flags, file.is_some())
        .map_err(|e| format!("invalid arguments: {e}"))?;
    crate::feasrelax::feasrelax_core(lp, &a, file.as_ref())
        .map(|(v, _)| v)
        .map_err(|e| e.msg)
}

/// SHA-256 over the bytes that define the model (sense, offset, dimensions, costs, bounds,
/// integrality, the column-wise matrix), for a model that is no file. Names are left out:
/// they do not change the problem. Raw bytes rather than a `Debug` dump, so hashing a large model
/// costs a memory pass, not a formatting pass.
pub fn model_digest(lp: &HighsLp) -> String {
    digest_hex(&model_digest_raw(lp))
}

/// Lower-case hex of a digest from [`model_digest_raw`].
pub fn digest_hex(d: &[u8; 32]) -> String {
    d.iter().map(|b| format!("{b:02x}")).collect()
}

/// [`model_digest`] as raw bytes (no hex string, no allocation): the Python API keeps these per
/// solve and formats them only when the envelope is asked for.
pub fn model_digest_raw(lp: &HighsLp) -> [u8; 32] {
    use sha2::{Digest, Sha256};
    use std::fmt::Write;
    fn f64s(h: &mut Sha256, v: &[f64]) {
        h.update((v.len() as u64).to_le_bytes());
        for x in v {
            h.update(x.to_le_bytes());
        }
    }
    /// `fmt::Write` straight into the hasher: the header is formatted without a `String`.
    struct Feed<'a>(&'a mut Sha256);
    impl Write for Feed<'_> {
        fn write_str(&mut self, s: &str) -> std::fmt::Result {
            self.0.update(s.as_bytes());
            Ok(())
        }
    }
    let mut h = Sha256::new();
    let _ = write!(
        Feed(&mut h),
        "{:?}|{}|{}|{}|",
        lp.sense_,
        lp.offset_,
        lp.num_col_,
        lp.num_row_
    );
    for v in [
        &lp.col_cost_,
        &lp.col_lower_,
        &lp.col_upper_,
        &lp.row_lower_,
        &lp.row_upper_,
        &lp.a_matrix_.value_,
    ] {
        f64s(&mut h, v);
    }
    h.update(
        lp.integrality_
            .iter()
            .map(|t| *t as u8)
            .collect::<Vec<u8>>(),
    );
    for v in [&lp.a_matrix_.start_, &lp.a_matrix_.index_] {
        h.update((v.len() as u64).to_le_bytes());
        for x in v.iter() {
            h.update(x.to_le_bytes());
        }
    }
    h.finalize().into()
}

/// The cheap facts of a finished solve that [`EnvParts::envelope`] needs to build the
/// `rusolve.result/1` envelope later. Taking them costs a few copies; building the envelope
/// (strings, JSON structs) is left to the caller that asks for it, so a solve that is never
/// serialised pays almost nothing (bench_api).
pub struct EnvParts {
    name: String,
    status: HighsModelStatus,
    run_status: HighsStatus,
    objective: f64,
    dual_bound: f64,
    gap: f64,
    nodes: i64,
    simplex_iterations: i64,
    has_solution: bool,
    time_s: f64,
    mip_stats: Option<MipStats>,
    seed: HighsInt,
    threads: HighsInt,
    time_limit: f64,
    mip_max_nodes: HighsInt,
    options: Vec<(String, String)>,
    counts: ModelCounts,
    digest: [u8; 32],
    convex_qp: bool,
}

/// The pieces of the `rusolve.result/1` envelope of the solve `h` just finished. `options` are
/// the settings the caller changed, in order (`run.options`); `name` is `input.path`.
pub fn env_parts(
    h: &Highs,
    run_status: HighsStatus,
    time_s: f64,
    name: &str,
    options: &[(String, String)],
) -> EnvParts {
    let opts = h.get_options();
    // `threads = 0` means "let the solver choose" (HiGHS: half the hardware threads, at least
    // one); the schema's `run.threads` is the count in force, at least 1.
    let threads = if opts.threads > 0 {
        opts.threads
    } else {
        // cached: `available_parallelism` reads cgroup files on Linux (tens of microseconds)
        static AUTO: std::sync::OnceLock<HighsInt> = std::sync::OnceLock::new();
        *AUTO.get_or_init(|| {
            std::thread::available_parallelism().map_or(1, |n| (n.get() / 2).max(1)) as HighsInt
        })
    };
    let lp = &h.model_lp_;
    let info = h.get_info();
    EnvParts {
        name: name.to_string(),
        status: h.get_model_status(),
        run_status,
        objective: info.objective_function_value,
        dual_bound: info.mip_dual_bound,
        gap: info.mip_gap,
        nodes: info.mip_node_count,
        simplex_iterations: info.simplex_iteration_count as i64,
        has_solution: h.solution_.value_valid,
        time_s,
        mip_stats: h.mip_stats_.clone(),
        seed: opts.random_seed,
        threads,
        time_limit: opts.time_limit,
        mip_max_nodes: opts.mip_max_nodes,
        options: options.to_vec(),
        counts: ModelCounts {
            rows: lp.num_row_ as i64,
            cols: lp.num_col_ as i64,
            integer_cols: lp
                .integrality_
                .iter()
                .filter(|t| **t == HighsVarType::kInteger)
                .count() as i64,
            nonzeros: lp.a_matrix_.value_.len() as i64,
        },
        digest: model_digest_raw(lp),
        convex_qp: false,
    }
}

impl EnvParts {
    /// Marks the solve as a convex QP (the Python API keeps the Hessian outside the LP, so
    /// the class cannot be read from the model).
    pub fn set_convex_qp(&mut self) {
        self.convex_qp = true;
    }

    /// The `rusolve.result/1` envelope.
    pub fn envelope(&self) -> Envelope {
        let info = HighsInfo {
            objective_function_value: self.objective,
            mip_dual_bound: self.dual_bound,
            mip_gap: self.gap,
            mip_node_count: self.nodes,
            simplex_iteration_count: self.simplex_iterations as HighsInt,
            ..HighsInfo::default()
        };
        let solution = HighsSolution {
            value_valid: self.has_solution,
            ..HighsSolution::default()
        };
        let lp = HighsLp::default();
        let solved = Solved {
            path: &self.name,
            lp: &lp,
            status: self.status,
            run_status: self.run_status,
            info: &info,
            solution: &solution,
            time_s: self.time_s,
            mip_stats: self.mip_stats.as_ref(),
            seed: self.seed,
            threads: self.threads,
            time_limit: self.time_limit,
            node_limit: None,
            mip_max_nodes: self.mip_max_nodes,
            options: &self.options,
            include_solution: false,
            wasm: None,
        };
        let mut env =
            envelope::envelope_with_counts(&solved, &self.counts, digest_hex(&self.digest));
        if self.convex_qp {
            env.class = "convex_qp".to_string();
            env.convexity = Some("convex");
        }
        env
    }
}

/// The envelope as a JSON value, with the fix the schema needs (an LP has no tree: HiGHS leaves
/// `mip_node_count` at -1, the schema has minimum 0; `mcp_tools.rs::solve` clamps the same way)
/// and, when `values` is given, `solution.columns` with the model's column names.
pub fn envelope_value(
    env: &Envelope,
    names: &[String],
    values: Option<&[f64]>,
) -> Result<Value, String> {
    let mut v = serde_json::to_value(env).map_err(|e| format!("internal error: {e}"))?;
    if v["stats"]["nodes"].as_i64().is_some_and(|n| n < 0) {
        v["stats"]["nodes"] = Value::from(0);
    }
    if let Some(x) = values {
        let cols: Vec<Value> = x
            .iter()
            .enumerate()
            .map(|(j, val)| {
                let name = match names.get(j) {
                    Some(n) if !n.is_empty() => n.clone(),
                    _ => format!("c{j}"),
                };
                serde_json::json!({"name": name, "value": val})
            })
            .collect();
        v["solution"] = serde_json::json!({ "columns": cols });
    }
    Ok(v)
}

/// Writes `lp` (with the objective Hessian `hessian`, full square `(row, col, value)` entries, or
/// empty) to `path` in the format its extension names (`.mps`, `.lp`, either with `.gz`), through
/// the writer `rusolve convert` uses. Semi-continuous columns are refused: the file formats keep
/// them as general constraints, which this path does not carry.
pub fn write_model(
    lp: &HighsLp,
    hessian: &[(usize, usize, f64)],
    path: &str,
) -> Result<(), String> {
    use rusolve_io::{ObjSense as IoSense, RowSense, VariableBound};
    let lower = path.to_ascii_lowercase();
    let stem = lower.strip_suffix(".gz").unwrap_or(&lower);
    if !(stem.ends_with(".lp") || stem.ends_with(".mps")) {
        return Err(format!(
            "{path}: unknown output format (expected .lp, .mps, .lp.gz or .mps.gz)"
        ));
    }
    if !hessian.is_empty() && stem.ends_with(".mps") {
        // `rusolve_io::write_mps` writes no QUADOBJ section: refuse rather than drop the Hessian
        return Err(
            "a quadratic objective cannot be written as MPS (the writer has no quadratic section); \
             write .lp"
                .into(),
        );
    }
    if lp.integrality_.iter().any(|t| {
        matches!(
            t,
            HighsVarType::kSemiContinuous | HighsVarType::kSemiInteger
        )
    }) {
        return Err("semi-continuous columns cannot be written yet".into());
    }
    let mut lp = lp.clone();
    lp.a_matrix_.ensure_colwise();
    let n = lp.num_col_ as usize;
    let nr = lp.num_row_ as usize;
    let mut triplets = Vec::with_capacity(lp.a_matrix_.value_.len());
    for j in 0..n {
        for k in lp.a_matrix_.start_[j] as usize..lp.a_matrix_.start_[j + 1] as usize {
            triplets.push((lp.a_matrix_.index_[k] as usize, j, lp.a_matrix_.value_[k]));
        }
    }
    let (mut senses, mut rhs, mut ranges) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..nr {
        let (lo, up) = (lp.row_lower_[i], lp.row_upper_[i]);
        let (lo_fin, up_fin) = (lo.is_finite() && lo > -1e20, up.is_finite() && up < 1e20);
        match (lo_fin, up_fin) {
            (true, true) if lo == up => {
                senses.push(RowSense::Equal);
                rhs.push(lo);
                ranges.push(0.0);
            }
            (true, true) => {
                senses.push(RowSense::GreaterEqual);
                rhs.push(lo);
                ranges.push(up - lo);
            }
            (true, false) => {
                senses.push(RowSense::GreaterEqual);
                rhs.push(lo);
                ranges.push(0.0);
            }
            (false, true) => {
                senses.push(RowSense::LessEqual);
                rhs.push(up);
                ranges.push(0.0);
            }
            (false, false) => {
                // a free row constrains nothing: keep it as `>= -1e30`
                senses.push(RowSense::GreaterEqual);
                rhs.push(-1e30);
                ranges.push(0.0);
            }
        }
    }
    let name_of = |names: &[String], i: usize, p: &str| match names.get(i) {
        Some(s) if !s.is_empty() => s.clone(),
        _ => format!("{p}{i}"),
    };
    let model = rusolve_io::Model {
        name: "rusolve".to_string(),
        num_variables: n,
        num_constraints: nr,
        obj_sense: if lp.sense_ == ObjSense::kMaximize {
            IoSense::Maximize
        } else {
            IoSense::Minimize
        },
        obj_name: Some("obj".to_string()),
        obj_offset: lp.offset_,
        objective: lp.col_cost_.clone(),
        senses,
        rhs,
        ranges,
        bounds: (0..n)
            .map(|j| VariableBound {
                lower: lp.col_lower_[j],
                upper: lp.col_upper_[j],
                integer: lp
                    .integrality_
                    .get(j)
                    .is_some_and(|t| *t == HighsVarType::kInteger),
            })
            .collect(),
        matrix: rusolve_core::CscMatrix::from_triplets(nr, n, &triplets),
        row_names: (0..nr).map(|i| name_of(&lp.row_names_, i, "r")).collect(),
        col_names: (0..n).map(|j| name_of(&lp.col_names_, j, "c")).collect(),
        hessian: hessian.to_vec(),
        gc: Default::default(),
        conic: Default::default(),
    };
    rusolve_io::write_model_file_report(path, &model)
        .map(|_| ())
        .map_err(|e| format!("write {path}: {e}"))
}
