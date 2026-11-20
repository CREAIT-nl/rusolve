//! The result envelope, the solution file by name and names that survive, for
//! the port MIP path: `RUSOLVE_MIP_HIGHS=1 rusolve mip <file>`
//! (`crate::mip_highs::run_mip_highs`).
//!
//! Everything here is *reporting*: no value computed here is read by the
//! solver, so the default search cannot move (`tools/bench/identity.sh` stays
//! `IDENTICAL 20 of 20`).
//!
//! ## The envelope
//!
//! `--format json` makes stdout carry exactly one JSON object of schema
//! `rusolve.result/1`. `docs/schema/rusolve-result-1.schema.json` is the
//! normative schema and `crates/rusolve-cli/tests/envelope.rs` checks a real
//! envelope against it. With `--opt output_flag=true` the HiGHS log is routed to
//! **stderr**: `log_to_console` is switched off and the log is re-emitted from
//! the `kCallbackLogging` callback (`rusolve_ekk::lp_data::highs_callback`, the
//! mechanism `mip-trace` already uses), so no `print!` in the port reaches
//! stdout and the port's IO stays untouched.
//!
//! ## Names
//!
//! `mip_highs::read_mip` fills `HighsLp::col_names_`/`row_names_` as HiGHS's MPS
//! reader does, so `solution.columns[].name` here and the names in the solution
//! file are the model's own. The reported solution is the **original** model's:
//! `Highs::call_solve_mip` copies the MIP solver's solution and computes
//! `solution_.row_value` against `model_lp_.a_matrix_`
//! (`rusolve-ekk/src/lp_data/highs_run.rs:1100-1113`), i.e. indexed by the
//! original columns and rows, while presolve shrinks name vectors only in its
//! own copy (`rusolve-ekk/src/presolve/hpresolve/shrink.rs:22-177`). Names are
//! therefore mapped by *original* index.
//!
//! ## The solution file
//!
//! `--solution-file <path>` writes the HiGHS raw solution file, ported below
//! with the real column and row names; `--solution-file <path>.json` writes the
//! `solution` object of the envelope instead.

use rusolve_ekk::mip::mip_stats::MipStats;
use rusolve_ekk::prelude::*;
use serde::Serialize;
use serde_json::{Map, Value};
use std::io::Write;

/// The `schema_version` of the envelope. The normative schema lives in
/// `docs/schema/rusolve-result-1.schema.json`; bump the `-N` suffix of that
/// file name and this string together on any breaking change.
pub const SCHEMA_VERSION: &str = "rusolve.result/1";

/// The `solver.name` / `solver.engine` constants of the envelope.
pub const SOLVER_NAME: &str = "rusolve";
pub const SOLVER_ENGINE: &str = "highs-port";

/// The envelope `status` for a `HighsModelStatus`
/// (`rusolve-ekk/src/lp_data/hconst.rs:254-277`).
///
/// One explicit arm per variant, so a HiGHS status added by a future rebase is
/// a compile error here rather than a silent alias. Mapping (the `status` column
/// is also the JSON-schema enum; the HiGHS text of
/// `utilModelStatusToString`/`HighsModelUtils.cpp:1415-1481` is in brackets):
///
/// | `HighsModelStatus`        | `status`                | why                                                    |
/// |---------------------------|-------------------------|--------------------------------------------------------|
/// | `kOptimal`                | `Optimal`               |                                                        |
/// | `kInfeasible`             | `Infeasible`            |                                                        |
/// | `kUnbounded`              | `Unbounded`             |                                                        |
/// | `kUnboundedOrInfeasible`  | `InfeasibleOrUnbounded` | the status HiGHS allows when it cannot tell (Highs.cpp:4855-4872) |
/// | `kTimeLimit`              | `TimeLimit`             |                                                        |
/// | `kIterationLimit`         | `IterationLimit`        |                                                        |
/// | `kSolutionLimit`          | `SolutionLimit`, or `NodeLimit` | HiGHS has no separate node-limit status: `mip_max_nodes`, `mip_max_leaves`, `mip_max_improving_sols` and `mip_max_stall_nodes` all set `kSolutionLimit` (`HighsMipSolverData.cpp:2647-2676`, `HighsMipSolver.cpp:1124-1129`). [`envelope`] reports `NodeLimit` only when the node count reached `mip_max_nodes`, so an agent is never told a node limit stopped a solve that the improving-solution, leaf or stall limit stopped. The CSV line keeps its `NodeLimit` spelling for every case (`mip_highs::csv_status`). |
/// | `kInterrupt`              | `Interrupted`           |                                                        |
/// | `kHighsInterrupt`         | `Interrupted`           | an interrupt raised by the solver itself (terminator, `HighsMipSolver.cpp:1211`, `HighsMipSolverData.cpp:3085`), not an error |
/// | `kNotset`                 | `NotSolved`             | no status was ever set (e.g. the model was never run)  |
/// | `kMemoryLimit`            | `MemoryLimit`           | (`Memory limit reached`)                               |
/// | `kModelEmpty`             | `ModelEmpty`            | no columns or no rows; HiGHS returns `kOk` (Highs.cpp:4824-4830) |
/// | `kObjectiveBound`         | `ObjectiveBound`        | the `objective_bound` option was reached               |
/// | `kObjectiveTarget`        | `ObjectiveTarget`       | the `objective_target` option was reached              |
/// | `kUnknown`                | `Unknown`               |                                                        |
/// | `kLoadError`              | `Error`                 |                                                        |
/// | `kModelError`             | `Error`                 | the five HiGHS error statuses (`Highs.cpp:4806-4812`)  |
/// | `kPresolveError`          | `Error`                 |                                                        |
/// | `kSolveError`             | `Error`                 |                                                        |
/// | `kPostsolveError`         | `Error`                 |                                                        |
///
/// `run_status == HighsStatus::kError` is reported as `Error` by the caller,
/// whatever the model status says, exactly as the CSV line does today
/// (`mip_highs.rs:213-217`).
pub fn result_status(status: HighsModelStatus) -> &'static str {
    match status {
        HighsModelStatus::kOptimal => "Optimal",
        HighsModelStatus::kInfeasible => "Infeasible",
        HighsModelStatus::kUnbounded => "Unbounded",
        HighsModelStatus::kUnboundedOrInfeasible => "InfeasibleOrUnbounded",
        HighsModelStatus::kTimeLimit => "TimeLimit",
        HighsModelStatus::kIterationLimit => "IterationLimit",
        HighsModelStatus::kSolutionLimit => "SolutionLimit",
        HighsModelStatus::kInterrupt => "Interrupted",
        HighsModelStatus::kHighsInterrupt => "Interrupted",
        HighsModelStatus::kNotset => "NotSolved",
        HighsModelStatus::kMemoryLimit => "MemoryLimit",
        HighsModelStatus::kModelEmpty => "ModelEmpty",
        HighsModelStatus::kObjectiveBound => "ObjectiveBound",
        HighsModelStatus::kObjectiveTarget => "ObjectiveTarget",
        HighsModelStatus::kUnknown => "Unknown",
        HighsModelStatus::kLoadError
        | HighsModelStatus::kModelError
        | HighsModelStatus::kPresolveError
        | HighsModelStatus::kSolveError
        | HighsModelStatus::kPostsolveError => "Error",
    }
}

/// JSON cannot carry `inf`/`nan`. Non-finite solver values (a dual bound of
/// `-inf` before the first incumbent, a gap of `inf` with no incumbent, the
/// `inf` default of `time_limit`) become JSON `null` rather than the `null`
/// that serde_json would otherwise write for them silently.
fn finite(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}

/// `solver` of the envelope (what to check out and which profile to set to
/// run the same code again).
#[derive(Serialize)]
pub struct Solver {
    pub name: &'static str,
    pub version: &'static str,
    /// `build.rs`: short commit hash, `-dirty` appended when tracked files had
    /// uncommitted changes. Omitted when the build had no git, where an absent
    /// field is honest and a wrong one is not.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub git: Option<&'static str>,
    /// The `RUSOLVE_PROFILE` in force: `recommended` or `highs-1.15.1`.
    pub profile: &'static str,
    pub engine: &'static str,
}

/// The profile name `crate::base_options` selects from `RUSOLVE_PROFILE` (main.rs
/// exits with a usage error on any other value before a solve starts).
pub(crate) fn profile_name() -> &'static str {
    match std::env::var("RUSOLVE_PROFILE").as_deref() {
        Ok("highs-1.15.1") => "highs-1.15.1",
        _ => "recommended",
    }
}

/// Lower-case hex SHA-256 of the file at `path` as read from disk: for a
/// `.gz` model these are the compressed bytes, so the hash can be checked with
/// `shasum -a 256` on the file the user has.
pub(crate) fn sha256_of_file(path: &str) -> std::io::Result<String> {
    use sha2::{Digest, Sha256};
    let digest = Sha256::digest(std::fs::read(path)?);
    Ok(digest.iter().map(|b| format!("{b:02x}")).collect())
}

/// `input` of the envelope: the model as read, before presolve.
#[derive(Serialize)]
pub struct Input {
    pub path: String,
    /// SHA-256 (lower-case hex) of the model file's bytes as read from disk, the
    /// compressed bytes for a `.gz`.
    pub sha256: String,
    pub rows: i64,
    pub cols: i64,
    /// Columns with `HighsVarType::kInteger` in `HighsLp::integrality_` (0 when
    /// the vector is empty, i.e. a pure LP, as HiGHS leaves it).
    pub integer_cols: i64,
    /// Structural nonzeros of `HighsLp::a_matrix_`.
    pub nonzeros: i64,
}

/// `stats` of the envelope.
#[derive(Serialize)]
pub struct Stats {
    /// `HighsInfo::mip_node_count`.
    pub nodes: i64,
    /// `HighsInfo::simplex_iteration_count`: the total LP iterations of the
    /// tree, the same number the CSV line reports.
    pub lp_iterations: i64,
    /// Wall-clock seconds of read + solve + report.
    pub time_s: f64,
    /// Deterministic work of the solve, in simplex synthetic ticks (`rusolve_ekk::work`):
    /// the same model, options and seed give the same number on every run and machine. 0 when
    /// the model took no MIP path (an LP: its ticks are not collected).
    // TODO: a `work_limit` option needs the option registry.
    pub work: u64,
    /// Per-phase statistics of the MIP solve (`rusolve_ekk::mip::mip_stats`); absent when the
    /// model took no MIP path.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phases: Option<Value>,
    /// `stats.wasm` — the loaded Wasm hook component and its call counters. Absent when no
    /// component was loaded (`--wasm-hook` / MCP `wasm_hook`), so every other envelope is
    /// unchanged.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wasm: Option<Value>,
}

/// The `stats.wasm` object (docs/design/wasm-hooks.md): what the component was called and
/// what it did. Built once here so the CLI and the MCP server report the same shape.
pub fn wasm_stats_json(name: &str, s: &rusolve_wasm::WasmStats) -> Value {
    serde_json::json!({
        "component": name,
        "calls": s.calls,
        "ok": s.ok,
        "declined": s.declined,
        "traps": s.traps,
        "out_of_fuel": s.out_of_fuel,
        "invalid": s.invalid,
        "fuel_used": s.fuel_used,
    })
}

/// `run` of the envelope: everything needed to attempt the same run again
/// (this keeps growing).
#[derive(Serialize)]
pub struct Run {
    pub seed: i64,
    pub threads: i64,
    /// `null` when unlimited (the HiGHS default, `kHighsInf`).
    pub time_limit: Option<f64>,
    /// `null` when `--node-limit` was not given (`mip_max_nodes` unlimited).
    pub node_limit: Option<i64>,
    /// The `--opt key=value` settings as given, as a JSON object. A repeated key
    /// keeps its last value, since that is what the solver keeps; serde_json's
    /// default map writes keys sorted, not in argument order.
    pub options: Map<String, Value>,
}

/// `limit.improving`: whether the incumbent was still improving when a limit ended the
/// solve. Only a run with the progress tracker (`--log-format jsonl` or `--watchdog`) sees the
/// improving-solution callbacks, so `tracked` is false (and the rest null) otherwise.
#[derive(Serialize)]
pub struct Improving {
    pub tracked: bool,
    /// Improving solutions seen, the first one included.
    pub count: Option<u64>,
    /// Seconds from the solve start to the last improving solution; null when none was seen.
    pub last_improvement_s: Option<f64>,
    /// `time_s` minus `last_improvement_s`.
    pub since_last_improvement_s: Option<f64>,
    /// Reporting convenience (no solver logic reads it): an incumbent exists and the last
    /// improvement came within the last tenth of the run, `since <= 0.1 * time_s`.
    pub still_improving: Option<bool>,
}

/// `limit` of the envelope: present whenever a limit or an interrupt ended the solve.
/// The incumbent, bound and gap are the top-level `objective`, `dual_bound` and `gap`.
#[derive(Serialize)]
pub struct Limit {
    /// `time`, `node`, `solution`, `iteration`, `memory`, `objective_bound`,
    /// `objective_target` or `interrupt`.
    pub kind: &'static str,
    /// `solver` (HiGHS checked the limit itself), `watchdog_interrupt` (the wall-clock
    /// watchdog requested the interrupt callback) or `watchdog_kill` (the solve did not
    /// react within the grace period and the watchdog reported the last state it had seen).
    pub enforced_by: &'static str,
    pub improving: Improving,
}

/// The limit a status stands for, `None` for a status that is not a limit.
pub(crate) fn limit_kind(status: &str) -> Option<&'static str> {
    Some(match status {
        "TimeLimit" => "time",
        "NodeLimit" => "node",
        "SolutionLimit" => "solution",
        "IterationLimit" => "iteration",
        "MemoryLimit" => "memory",
        "ObjectiveBound" => "objective_bound",
        "ObjectiveTarget" => "objective_target",
        "Interrupted" => "interrupt",
        _ => return None,
    })
}

/// One named column of `solution`.
#[derive(Serialize)]
pub struct Column {
    pub name: String,
    pub value: f64,
}

/// `solution` of the envelope: the incumbent in **original-model** column order.
#[derive(Serialize)]
pub struct Solution {
    pub columns: Vec<Column>,
}

/// The `rusolve.result/1` envelope. Field order is the declaration order.
#[derive(Serialize)]
pub struct Envelope {
    pub schema_version: &'static str,
    pub solver: Solver,
    pub input: Input,
    pub status: &'static str,
    /// `null` unless a primal solution exists.
    pub objective: Option<f64>,
    /// `HighsInfo::mip_dual_bound`, in the minimisation sense HiGHS reports it
    /// in; `null` when infinite (no incumbent yet).
    pub dual_bound: Option<f64>,
    /// `HighsInfo::mip_gap`: the relative gap HiGHS computes in
    /// `HighsMipSolverData::limitsToGap` (`HighsMipSolverData.cpp:2762`, ported
    /// at `mip/mip_solver_data/part1.rs:746`), copied to `info_.mip_gap` at
    /// `lp_data/highs_run.rs:1156`. A fraction, not a percentage. `null` when
    /// there is no incumbent (HiGHS's gap is then `inf`).
    pub gap: Option<f64>,
    pub has_solution: bool,
    pub stats: Stats,
    pub run: Run,
    /// Present only when `--solution` / `--include-solution` was given.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub solution: Option<Solution>,
    /// Present when a limit or an interrupt ended the solve.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub limit: Option<Limit>,
    /// The problem class the router assigned (`lp`, `mip`, `convex_qp`, `conic`, `gc`, ...).
    pub class: String,
    /// The engine that solved the class (`highs-port`, `native`, `qp`, `clarabel`,
    /// `gc-linearization`). `solver.engine` stays the stack that ran.
    pub engine: String,
    /// `convex` for a quadratic class that was solved; absent for a linear one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub convexity: Option<&'static str>,
}

/// What the CLI hands over after a solve: everything the envelope and the
/// solution file are built from, all of it in original-model space.
pub struct Solved<'a> {
    /// The path as the user gave it.
    pub path: &'a str,
    /// The incumbent model as passed to `Highs::pass_model` (`Highs::model_lp_`).
    pub lp: &'a HighsLp,
    /// `Highs::get_model_status()`.
    pub status: HighsModelStatus,
    /// The `HighsStatus` of `Highs::run` (`kError` moves `status` to `Error`).
    pub run_status: HighsStatus,
    pub info: &'a HighsInfo,
    pub solution: &'a HighsSolution,
    pub time_s: f64,
    /// `Highs::mip_stats_`: the statistics of the MIP solve.
    pub mip_stats: Option<&'a MipStats>,
    pub seed: HighsInt,
    pub threads: HighsInt,
    pub time_limit: f64,
    pub node_limit: Option<HighsInt>,
    /// `HighsOptions::mip_max_nodes` in effect (from `--node-limit` or `--opt`), to tell
    /// a node limit apart from the other limits behind `kSolutionLimit`.
    pub mip_max_nodes: HighsInt,
    /// The `--opt` settings in argument order.
    pub options: &'a [(String, String)],
    /// `--solution` / `--include-solution`.
    pub include_solution: bool,
    /// The `stats.wasm` object (`wasm_stats_json`), `Some` only when a Wasm hook component
    /// was loaded for this run.
    pub wasm: Option<Value>,
}

/// The name of original column `ix`: `HighsLp::col_names_` when it is there, else
/// HiGHS's minimal generated name (`kHighsMinimalColNamePrefix` = `"c"` + index,
/// `HConst.h:189`, the prefix `normaliseNames` falls back to at
/// `HighsModelUtils.cpp:503-516`).
fn col_name(lp: &HighsLp, ix: usize) -> String {
    match lp.col_names_.get(ix) {
        Some(name) if !name.is_empty() => name.clone(),
        _ => format!("c{ix}"),
    }
}

/// As [`col_name`] for rows (`kHighsMinimalrowNamePrefix` = `"r"`,
/// `HConst.h:190`).
fn row_name(lp: &HighsLp, ix: usize) -> String {
    match lp.row_names_.get(ix) {
        Some(name) if !name.is_empty() => name.clone(),
        _ => format!("r{ix}"),
    }
}

/// The `solution` object: the incumbent, one entry per original column, in
/// column order. Empty when there is no primal solution.
pub fn solution_of(solved: &Solved) -> Solution {
    let lp = solved.lp;
    if !solved.solution.value_valid {
        return Solution {
            columns: Vec::new(),
        };
    }
    let columns = (0..lp.num_col_ as usize)
        .map(|ix| Column {
            name: col_name(lp, ix),
            // C++ indexes `col_value[ix]` after asserting it is long enough
            // (HighsModelUtils.cpp:277-280); rusolve reports NaN rather than
            // panicking if a port bug leaves the vector short.
            value: solved
                .solution
                .col_value
                .get(ix)
                .copied()
                .unwrap_or(f64::NAN),
        })
        .collect();
    Solution { columns }
}

/// The `run.options` object from the `--opt` pairs (last value of a repeated
/// key wins).
pub(crate) fn options_of(options: &[(String, String)]) -> Map<String, Value> {
    let mut map = Map::new();
    for (key, value) in options {
        map.insert(key.clone(), Value::String(value.clone()));
    }
    map
}

/// The `stats.phases` object. Every number is read from what the solve recorded
/// (`rusolve_ekk::mip::mip_stats`); seconds are wall time and vary between runs, everything else
/// is deterministic for threads = 1. Phase times: `presolve` from the start of the solve to the end
/// of MIP presolve, `setup` the rest of the setup, `root.node_time_s` the root node (root LP, root
/// separation, root heuristics), `search.time_s` everything after it. The port's separate clocks
/// for the root LP and for heuristics are stubs, so neither has its own time.
fn phases_of(m: &MipStats) -> Value {
    use serde_json::json;
    let reductions: Map<String, Value> = m
        .reductions
        .iter()
        .map(|(name, n)| (name.clone(), json!(n)))
        .collect();
    let separators: Vec<Value> = m
        .separators
        .iter()
        .map(|s| {
            json!({"name": s.name, "calls": s.calls, "cuts_found": s.cuts_found, "time_s": s.time_s})
        })
        .collect();
    let sources: Vec<Value> = m
        .heuristics
        .iter()
        .map(|h| json!({"source": h.source, "improving_incumbents": h.improving}))
        .collect();
    json!({
        "presolve": {
            "time_s": m.presolve_s,
            "status": m.presolve_status,
            "restarts": m.num_restarts,
            "rows_before": m.rows_before,
            "cols_before": m.cols_before,
            "nonzeros_before": m.nnz_before,
            "rows_after": m.rows_after,
            "cols_after": m.cols_after,
            "nonzeros_after": m.nnz_after,
            "reductions": reductions,
        },
        "setup": {"time_s": m.setup_s},
        "root": {
            "lp_iterations": m.root_lp_iterations,
            "lp_objective": finite(m.root_lp_objective),
            "node_time_s": m.root_node_s,
        },
        "cuts": {
            "rounds": m.sepa_rounds,
            "pool_cuts": m.cutpool_cuts,
            "separators": separators,
        },
        "heuristics": {"improving_incumbents_by_source": sources},
        "search": {"time_s": m.search_s, "nodes": m.num_nodes, "leaves": m.num_leaves},
        "lp_iterations": {
            "total": m.lp_iterations_total,
            "strong_branching": m.lp_iterations_strong_branching,
            "separation": m.lp_iterations_separation,
            "heuristics": m.lp_iterations_heuristics,
            "by_source": rusolve_ekk::work::LP_SRC_NAMES
                .iter()
                .enumerate()
                .filter(|(i, _)| m.lp_solves_by_source[*i] > 0)
                .map(|(i, k)| {
                    (
                        (*k).to_string(),
                        json!({"iterations": m.lp_iters_by_source[i], "solves": m.lp_solves_by_source[i]}),
                    )
                })
                .collect::<Map<String, Value>>(),
        },
        "work": {
            "simplex_ticks": m.work_parts.simplex_ticks,
            "propagation_units": m.work_parts.propagation_units,
            "feasibility_jump_units": m.work_parts.fj_units,
            "presolve_units": m.work_parts.presolve_units,
            "unit_ticks": rusolve_ekk::work::NONSIMPLEX_UNIT_TICKS,
        },
        "propagation_events": rusolve_ekk::work::Event::NAMES
            .iter()
            .zip(m.events.iter())
            .map(|(k, v)| ((*k).to_string(), json!(v)))
            .collect::<Map<String, Value>>(),
        "first_incumbent": m.first_incumbent_s.map(|t| json!({
            "time_s": t,
            "work": m.first_incumbent_work,
        })),
        "primal_dual_integral": finite(m.primal_dual_integral),
        "total_time_s": m.total_s,
    })
}

/// Builds the `rusolve.result/1` envelope of a finished solve.
pub fn envelope(solved: &Solved) -> Envelope {
    // The model was read from this path moments ago, so a failure here
    // is not an input error but a vanished file.
    let sha256 =
        sha256_of_file(solved.path).unwrap_or_else(|e| panic!("hash {}: {e}", solved.path));
    envelope_with_digest(solved, sha256)
}

/// [`envelope`] for a model that is not a file (the Python API): `sha256` is the digest of
/// the model as the caller defines it, and `solved.path` is only a display name.
pub fn envelope_with_digest(solved: &Solved, sha256: String) -> Envelope {
    let lp = solved.lp;
    let counts = ModelCounts {
        rows: lp.num_row_ as i64,
        cols: lp.num_col_ as i64,
        integer_cols: lp
            .integrality_
            .iter()
            .filter(|t| **t == HighsVarType::kInteger)
            .count() as i64,
        // HiGHS: HighsSparseMatrix::countNnz, the nnz of the column-wise
        // matrix. The port has no countNnz, so the value buffer length.
        nonzeros: lp.a_matrix_.value_.len() as i64,
    };
    envelope_with_counts(solved, &counts, sha256)
}

/// The model-size facts `input` reports, taken when the solve finished so that a caller (the
/// Python API) can build the envelope later without keeping the model.
#[derive(Clone, Copy, Debug)]
pub struct ModelCounts {
    pub rows: i64,
    pub cols: i64,
    pub integer_cols: i64,
    pub nonzeros: i64,
}

/// [`envelope_with_digest`] with the model facts given. Reads nothing from `solved.lp` (a
/// caller that built `counts` earlier may pass any model there), except for `include_solution`.
pub fn envelope_with_counts(solved: &Solved, counts: &ModelCounts, sha256: String) -> Envelope {
    let info = solved.info;
    let has_solution = solved.solution.value_valid;
    let status = if solved.run_status == HighsStatus::kError {
        "Error"
    } else if solved.status == HighsModelStatus::kSolutionLimit
        && info.mip_node_count >= solved.mip_max_nodes as i64
    {
        "NodeLimit"
    } else {
        result_status(solved.status)
    };
    let route_fields = crate::route::envelope_fields(counts.integer_cols, SOLVER_ENGINE);
    Envelope {
        schema_version: SCHEMA_VERSION,
        solver: Solver {
            name: SOLVER_NAME,
            version: env!("CARGO_PKG_VERSION"),
            git: option_env!("RUSOLVE_GIT"),
            profile: profile_name(),
            engine: SOLVER_ENGINE,
        },
        input: Input {
            path: solved.path.to_string(),
            sha256,
            rows: counts.rows,
            cols: counts.cols,
            integer_cols: counts.integer_cols,
            nonzeros: counts.nonzeros,
        },
        status,
        objective: if has_solution {
            finite(info.objective_function_value)
        } else {
            None
        },
        dual_bound: finite(info.mip_dual_bound),
        gap: finite(info.mip_gap),
        has_solution,
        stats: Stats {
            nodes: info.mip_node_count,
            lp_iterations: info.simplex_iteration_count as i64,
            time_s: solved.time_s,
            work: solved.mip_stats.map_or(0, |m| m.work),
            phases: solved.mip_stats.map(phases_of),
            wasm: solved.wasm.clone(),
        },
        run: Run {
            seed: solved.seed as i64,
            threads: solved.threads as i64,
            time_limit: finite(solved.time_limit),
            node_limit: solved.node_limit.map(|n| n as i64),
            options: options_of(solved.options),
        },
        solution: solved.include_solution.then(|| solution_of(solved)),
        limit: limit_kind(status).map(|kind| Limit {
            kind,
            enforced_by: "solver",
            improving: Improving {
                tracked: false,
                count: None,
                last_improvement_s: None,
                since_last_improvement_s: None,
                still_improving: None,
            },
        }),
        class: route_fields.0,
        engine: route_fields.1,
        convexity: route_fields.2,
    }
}

/// `solver.engine` of an envelope built by [`native_envelope`].
pub const NATIVE_ENGINE: &str = "native";

/// What the native stack (`rusolve-mip`, the legacy `rusolve-lp` engines) hands over after a
/// solve, so it prints the same `rusolve.result/1` envelope as the port. Everything is in
/// the caller's objective sense, offset included.
#[derive(Clone, Copy)]
pub struct NativeSolved<'a> {
    pub path: &'a str,
    pub model: &'a rusolve_io::Model,
    /// One of the schema's `status` values.
    pub status: &'static str,
    /// `NaN` without an incumbent.
    pub objective: f64,
    pub dual_bound: f64,
    pub nodes: i64,
    pub lp_iterations: i64,
    pub time_s: f64,
    pub seed: i64,
    pub threads: i64,
    pub time_limit: Option<f64>,
    pub node_limit: Option<i64>,
    /// The incumbent, in column order; empty when there is none.
    pub x: &'a [f64],
    pub include_solution: bool,
}

/// The relative gap as HiGHS computes it in `HighsMipSolverData::limitsToGap`
/// (`HighsMipSolverData.cpp:2762`): `|ub - lb| / |ub|`, 0 when both are 0, infinite without an
/// incumbent. Non-finite values become `null` in the envelope.
fn native_gap(primal: f64, dual: f64) -> f64 {
    if !primal.is_finite() {
        f64::INFINITY
    } else if primal == 0.0 {
        if dual == 0.0 { 0.0 } else { f64::INFINITY }
    } else {
        (primal - dual).abs() / primal.abs()
    }
}

/// Builds the `rusolve.result/1` envelope of a native-stack solve.
pub fn native_envelope(n: &NativeSolved) -> Envelope {
    let has_solution = !n.x.is_empty();
    let m = n.model;
    let route_fields = crate::route::envelope_fields(
        m.bounds.iter().filter(|b| b.integer).count() as i64,
        NATIVE_ENGINE,
    );
    Envelope {
        schema_version: SCHEMA_VERSION,
        solver: Solver {
            name: SOLVER_NAME,
            version: env!("CARGO_PKG_VERSION"),
            git: option_env!("RUSOLVE_GIT"),
            profile: profile_name(),
            engine: NATIVE_ENGINE,
        },
        input: Input {
            path: n.path.to_string(),
            sha256: sha256_of_file(n.path).unwrap_or_else(|e| panic!("hash {}: {e}", n.path)),
            rows: m.num_constraints as i64,
            cols: m.num_variables as i64,
            integer_cols: m.bounds.iter().filter(|b| b.integer).count() as i64,
            nonzeros: m.matrix.values().len() as i64,
        },
        status: n.status,
        objective: if has_solution {
            finite(n.objective)
        } else {
            None
        },
        dual_bound: finite(n.dual_bound),
        gap: if has_solution {
            finite(native_gap(n.objective, n.dual_bound))
        } else {
            None
        },
        has_solution,
        stats: Stats {
            nodes: n.nodes,
            lp_iterations: n.lp_iterations,
            time_s: n.time_s,
            // Port-only: the native search has no work counter or phase statistics
            work: 0,
            phases: None,
            // Hooks are port-only too (the native search has no MipHooks)
            wasm: None,
        },
        run: Run {
            seed: n.seed,
            threads: n.threads,
            time_limit: n.time_limit.and_then(finite),
            node_limit: n.node_limit,
            options: Map::new(),
        },
        solution: n.include_solution.then(|| native_solution(n)),
        // The native search has no callbacks, so improvement is never tracked there
        limit: limit_kind(n.status).map(|kind| Limit {
            kind,
            enforced_by: "solver",
            improving: Improving {
                tracked: false,
                count: None,
                last_improvement_s: None,
                since_last_improvement_s: None,
                still_improving: None,
            },
        }),
        class: route_fields.0,
        engine: route_fields.1,
        convexity: route_fields.2,
    }
}

/// The native incumbent by column name (`c<index>` where the model has no name, as HiGHS).
fn native_solution(n: &NativeSolved) -> Solution {
    let m = n.model;
    Solution {
        columns: if n.x.is_empty() {
            Vec::new()
        } else {
            (0..m.num_variables)
                .map(|ix| Column {
                    name: match m.col_names.get(ix) {
                        Some(name) if !name.is_empty() => name.clone(),
                        _ => format!("c{ix}"),
                    },
                    value: n.x.get(ix).copied().unwrap_or(f64::NAN),
                })
                .collect()
        },
    }
}

/// The native envelope as the single JSON object `--format json` puts on stdout.
pub fn native_envelope_json(n: &NativeSolved) -> String {
    serde_json::to_string(&native_envelope(n)).unwrap_or_else(|e| format!("{{\"error\":{e:?}}}"))
}

/// The `solution` object for `--solution-file <path>.json` on the native stack.
pub fn native_solution_file_json(n: &NativeSolved) -> String {
    serde_json::to_string_pretty(&native_solution(n))
        .unwrap_or_else(|e| format!("{{\"error\":{e:?}}}"))
}

/// The envelope as the single JSON object `--format json` puts on stdout.
pub fn envelope_json(solved: &Solved) -> String {
    // Serialising our own plain structs: the only failure mode is a non-finite
    // double, which serde_json writes as `null`.
    serde_json::to_string(&envelope(solved)).unwrap_or_else(|e| format!("{{\"error\":{e:?}}}"))
}

/// The `solution` object of the envelope, for `--solution-file <path>.json`:
/// the same shape, so `{"columns": []}` when there is no incumbent.
pub fn solution_file_json(solved: &Solved) -> String {
    serde_json::to_string_pretty(&solution_of(solved))
        .unwrap_or_else(|e| format!("{{\"error\":{e:?}}}"))
}

/// Writes `--solution-file <path>`: the HiGHS raw solution file of a
/// finished solve, with the model's own column and row names.
pub fn write_solution_file(
    path: &str,
    log_options: HighsLogOptions,
    solved: &Solved,
) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut file = std::io::BufWriter::new(file);
    write_model_solution(
        &mut file,
        &log_options,
        solved.lp,
        solved.status,
        solved.solution,
        solved.info,
    );
    Ok(())
}

/// HiGHS 1.16.0 lp_data/HighsModelUtils.cpp::writeMipBoundGap (lines 205-216):
/// the dual bound and relative gap lines of a solved MIP.
///
/// `FILE*` -> `&mut dyn Write` and `highsFprintfString` -> `write_all`, as in
/// `rusolve_ekk::lp_data::model_utils::write_objective_value`. The `log_options`
/// argument of the C++ writer is dropped: it only feeds `highsFprintfString`'s
/// `stdout` branch, which is never taken when writing to a file.
fn write_mip_bound_gap(file: &mut dyn Write, mip_dual_bound: f64, mip_gap: f64) {
    let s = format!(
        "Objective_bound {}\nMIP_rel_gap {}\n",
        highs_double_to_string(mip_dual_bound, kHighsSolutionValueToStringTolerance),
        highs_double_to_string(mip_gap, kHighsSolutionValueToStringTolerance)
    );
    let _ = file.write_all(s.as_bytes());
}

/// HiGHS 1.16.0 lp_data/HighsModelUtils.cpp::writeModelSolution (lines 268-359)
/// for the default raw style (`sparse = false`, `partial = false`), reached from
/// `writeSolutionFile` (lines 529-601, the final `else` at lines 588-599) via
/// `Highs::writeSolution(filename, kSolutionStyleRaw)` (lp_data/Highs.cpp:3843
/// -3868; default style `kSolutionStyleRaw` at Highs.h:271-272). The objective
/// and primal-solution blocks are the already-ported writers in
/// `rusolve_ekk::lp_data::model_utils` (`write_lp_objective`, lines 189-194, and
/// `write_primal_solution`, lines 204-244 of the same file).
///
/// rusolve deviations, none of them reachable on the port MIP path:
/// * `Highs::writeSolution` calls `normaliseNames` first (Highs.cpp:3853,
///   HighsModelUtils.cpp:470-516) to replace absent, duplicated or
///   space-containing names. The port gives every column and constraint row a
///   non-empty unique name from the MPS file, and the reader rejects duplicate
///   row names outright (`rusolve-io/src/mps.rs:345`) where HiGHS warns and
///   drops both name vectors (`io/HMpsFF.cpp:64-78`).
/// * The `kSolutionStyleSparse` / `kSolutionStylePartial` / `kSolutionStylePretty`
///   / `kSolutionStyleGlpsol*` branches of `writeSolutionFile` are not ported,
///   nor is the `ranging` suffix of `Highs::writeSolution` (Highs.cpp:3871-3884),
///   which HiGHS rejects for a MIP.
/// * The `assert`s on name and solution sizes (lines 276-287) become index
///   guards: a short vector writes `nan` rather than panicking.
pub fn write_model_solution(
    file: &mut dyn Write,
    log_options: &HighsLogOptions,
    lp: &HighsLp,
    model_status: HighsModelStatus,
    solution: &HighsSolution,
    info: &HighsInfo,
) {
    let have_primal = solution.value_valid;
    let have_dual = solution.dual_valid;
    let _ = file.write_all(b"Model status\n");
    let s = format!("{}\n", util_model_status_to_string(model_status));
    let _ = file.write_all(s.as_bytes());
    let _ = file.write_all(b"\n# Primal solution values\n");
    if !have_primal || info.primal_solution_status == kSolutionStatusNone {
        let _ = file.write_all(b"None\n");
    } else {
        let _ = file.write_all(if info.primal_solution_status == kSolutionStatusFeasible {
            b"Feasible\n"
        } else {
            // C++ asserts the third case is infeasible (line 300).
            b"Infeasible\n"
        });
        // C++ calls writeModelObjective (lines 180-187), which adds the Hessian
        // term. `read_mip` never sets a Hessian, so the LP objective is the same
        // number here.
        write_lp_objective(file, log_options, lp, &solution.col_value);
        if lp.is_mip() && model_status == HighsModelStatus::kOptimal && info.mip_gap < kHighsInf {
            write_mip_bound_gap(file, info.mip_dual_bound, info.mip_gap);
        }
        write_primal_solution(file, log_options, lp, &solution.col_value, false);
        let _ = file.write_all(format!("# Rows {}\n", lp.num_row_).as_bytes());
        for ix in 0..lp.num_row_ as usize {
            let value = solution.row_value.get(ix).copied().unwrap_or(f64::NAN);
            let s = format!(
                "{} {}\n",
                row_name(lp, ix),
                highs_double_to_string(value, kHighsSolutionValueToStringTolerance)
            );
            let _ = file.write_all(s.as_bytes());
        }
    }
    let _ = file.write_all(b"\n# Dual solution values\n");
    if !have_dual || info.dual_solution_status == kSolutionStatusNone {
        let _ = file.write_all(b"None\n");
    } else {
        let _ = file.write_all(if info.dual_solution_status == kSolutionStatusFeasible {
            b"Feasible\n"
        } else {
            b"Infeasible\n"
        });
        let _ = file.write_all(format!("# Columns {}\n", lp.num_col_).as_bytes());
        for ix in 0..lp.num_col_ as usize {
            let value = solution.col_dual.get(ix).copied().unwrap_or(f64::NAN);
            let s = format!(
                "{} {}\n",
                col_name(lp, ix),
                highs_double_to_string(value, kHighsSolutionValueToStringTolerance)
            );
            let _ = file.write_all(s.as_bytes());
        }
        let _ = file.write_all(format!("# Rows {}\n", lp.num_row_).as_bytes());
        for ix in 0..lp.num_row_ as usize {
            let value = solution.row_dual.get(ix).copied().unwrap_or(f64::NAN);
            let s = format!(
                "{} {}\n",
                row_name(lp, ix),
                highs_double_to_string(value, kHighsSolutionValueToStringTolerance)
            );
            let _ = file.write_all(s.as_bytes());
        }
    }
    // Highs::writeSolution appends the basis for the raw style
    // (Highs.cpp:3867-3870) through writeBasisFile (HighsLpUtils.cpp:2828-2844):
    // the version line, then `None`, because `basis_.valid` is false - a MIP run
    // leaves the basis invalidated (`Highs::call_solve_mip` writes solution_ and
    // info_ only, after `invalidate_basis`).
    let _ = file.write_all(b"\n# Basis\n");
    let _ = file.write_all(format!("HiGHS_basis_file {kHighsBasisFileV2}\nNone\n").as_bytes());
    let _ = file.flush();
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every `HighsModelStatus` variant, so the test below fails to compile if
    /// the port gains or loses one.
    const ALL_MODEL_STATUSES: [HighsModelStatus; 20] = [
        HighsModelStatus::kNotset,
        HighsModelStatus::kLoadError,
        HighsModelStatus::kModelError,
        HighsModelStatus::kPresolveError,
        HighsModelStatus::kSolveError,
        HighsModelStatus::kPostsolveError,
        HighsModelStatus::kModelEmpty,
        HighsModelStatus::kOptimal,
        HighsModelStatus::kInfeasible,
        HighsModelStatus::kUnboundedOrInfeasible,
        HighsModelStatus::kUnbounded,
        HighsModelStatus::kObjectiveBound,
        HighsModelStatus::kObjectiveTarget,
        HighsModelStatus::kTimeLimit,
        HighsModelStatus::kIterationLimit,
        HighsModelStatus::kUnknown,
        HighsModelStatus::kSolutionLimit,
        HighsModelStatus::kInterrupt,
        HighsModelStatus::kMemoryLimit,
        HighsModelStatus::kHighsInterrupt,
    ];

    /// The envelope's `status` enum, read from the committed schema rather than
    /// repeated here, so the mapping and the schema cannot drift apart. The
    /// matching test on the produced envelope is
    /// `envelope_reports_a_status_the_schema_allows` in `tests/envelope.rs`.
    fn schema_status_enum() -> Vec<String> {
        const SCHEMA: &str = include_str!("../../../docs/schema/rusolve-result-1.schema.json");
        let schema: serde_json::Value = serde_json::from_str(SCHEMA).unwrap();
        schema["properties"]["status"]["enum"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect()
    }

    /// The status mapping names every HiGHS model status explicitly and
    /// only ever returns a value of the envelope's closed enum - and every value
    /// of that enum is reachable, so the schema carries no dead value.
    #[test]
    fn envelope_status_mapping_covers_every_model_status() {
        assert_eq!(
            ALL_MODEL_STATUSES.len(),
            20,
            "the port has a different number of HighsModelStatus variants"
        );
        let allowed = schema_status_enum();
        let mut seen = std::collections::BTreeSet::new();
        for status in ALL_MODEL_STATUSES {
            let name = result_status(status);
            assert!(
                allowed.iter().any(|a| a == name),
                "{status:?} maps to {name:?}, which the schema's status enum does not allow"
            );
            seen.insert(name);
        }
        // `NodeLimit` is not a `result_status` value: `envelope` refines
        // `kSolutionLimit` into it when the node count reached `mip_max_nodes`
        // (`envelope_separates_the_node_limit_from_the_other_solution_limits`).
        seen.insert("NodeLimit");
        for allowed in &allowed {
            assert!(
                seen.contains(allowed.as_str()),
                "{allowed} is in the schema's status enum but no model status maps to it"
            );
        }
        // The distinctions an agent acts on survive the mapping.
        assert_eq!(result_status(HighsModelStatus::kOptimal), "Optimal");
        assert_eq!(
            result_status(HighsModelStatus::kUnboundedOrInfeasible),
            "InfeasibleOrUnbounded"
        );
        assert_ne!(
            result_status(HighsModelStatus::kUnboundedOrInfeasible),
            result_status(HighsModelStatus::kUnbounded)
        );
        assert_ne!(
            result_status(HighsModelStatus::kTimeLimit),
            result_status(HighsModelStatus::kSolutionLimit)
        );
        assert_ne!(
            result_status(HighsModelStatus::kIterationLimit),
            result_status(HighsModelStatus::kTimeLimit)
        );
        assert_ne!(
            result_status(HighsModelStatus::kNotset),
            result_status(HighsModelStatus::kSolveError)
        );
        assert_ne!(
            result_status(HighsModelStatus::kInterrupt),
            result_status(HighsModelStatus::kModelError)
        );
    }
}
