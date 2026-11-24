//! The tools of `rusolve mcp` and the model handles behind them.
//!
//! A handle (`model_id`) owns one live ported `Highs` instance. Its
//! `model_lp_` is the `HighsLp` of the model as edited so far, and the instance
//! keeps its last basis, so after `modify` an LP re-solve starts from that
//! basis (HiGHS `Highs::run` skips presolve and warm-starts when `basis_.valid`;
//! `changeColBounds`, `changeRowBounds`, `changeCoeff`, `addRow` and the rest
//! of `highs_modify_api.rs` keep the basis valid, `Highs.cpp` 2846-3295). There
//! is no session-global "current model": every tool takes a `model_id`.
//!
//! Each tool is a thin wrapper over code the subcommands already use:
//! `solve` -> `envelope::envelope` (`rusolve.result/1`), `check_solution` ->
//! `check::check_solution`, `explain` -> `explain::explain_core`, `lint` ->
//! `lint::lint_core`, `ranging` -> `ranging::ranging_core`, `feasrelax` ->
//! `feasrelax::feasrelax_core`. The CLI paths call the same cores, so the JSON
//! an agent gets is the JSON `--format json` prints.

use crate::envelope;
use crate::mip_highs::read_mip;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::{Map, Value, json};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

type Args = Map<String, Value>;
type ToolResult = Result<Value, String>;

const SCHEMA_RESULT: &str = include_str!("../../../docs/schema/rusolve-result-1.schema.json");
const SCHEMA_CHECK: &str = include_str!("../../../docs/schema/rusolve-check-1.schema.json");
const SCHEMA_EXPLAIN: &str = include_str!("../../../docs/schema/rusolve-explain-1.schema.json");
const SCHEMA_LINT: &str = include_str!("../../../docs/schema/rusolve-lint-1.schema.json");
const SCHEMA_RANGING: &str = include_str!("../../../docs/schema/rusolve-ranging-1.schema.json");
const SCHEMA_FEASRELAX: &str = include_str!("../../../docs/schema/rusolve-feasrelax-1.schema.json");

/// Columns `get_solution` returns when the caller does not say.
const DEFAULT_PAGE: usize = 1000;

// ------------------------------------------------------------------ handles

struct Model {
    /// A file that exists for the lifetime of the handle: the loaded path, or
    /// the temp file `mps_text` was written to. `envelope::envelope` hashes it.
    source_path: String,
    /// `Some` when `source_path` is ours to delete.
    temp: Option<PathBuf>,
    highs: Highs,
    edits_applied: usize,
    /// A solve ran since the last edit, so `get_solution` has an answer.
    solved: bool,
    /// A model the router sends to an engine other than the LP/MIP port (a convex QP, a
    /// conic model, general constraints). `highs` is then an empty placeholder.
    routed: Option<Routed>,
    /// Scheduling and combinatorial (CP) objects from `build_model`'s `cp` section. No file format
    /// carries them, so the handle keeps them; `highs` holds the user's model and `solve` runs
    /// the CP path (`cp_solve.rs`) in this process.
    cp: Option<rusolve_io::CpModel>,
    /// The last CP solve: the solver (its solution cut to the user's columns) and the `cp` report.
    cp_result: Option<CpResult>,
}

struct CpResult {
    highs: Highs,
}

/// A handle whose class `rusolve solve` routes to another engine (`route.rs`).
/// `solve` runs `rusolve-cli solve <file> --format json` as a child process, the same
/// command a user would run, so the answer is the CLI's and a cancel can kill the child.
struct Routed {
    class: crate::route::Class,
    /// The `solution.columns` of the last solve, for `get_solution`.
    columns: Option<Vec<Value>>,
    status: Option<String>,
    objective: Option<f64>,
}

/// What a running tool call can see of its request: the progress token, the cancel flag and the
/// stdout writer (`mcp.rs`). Progress follows
/// <https://modelcontextprotocol.io/specification/2026-07-28/basic/utilities/progress>, cancellation
/// <https://modelcontextprotocol.io/specification/2026-07-28/basic/utilities/cancellation>.
pub struct Ctx {
    pub token: Option<Value>,
    pub cancel: Arc<AtomicBool>,
    pub out: Arc<crate::mcp::Out>,
    /// (counter, time of the last notification); `progress` must increase with each one.
    last: Mutex<(u64, Option<Instant>)>,
}

/// Minimum seconds between two progress notifications of one call.
const PROGRESS_INTERVAL: Duration = Duration::from_millis(200);

impl Ctx {
    pub fn new(token: Option<Value>, cancel: Arc<AtomicBool>, out: Arc<crate::mcp::Out>) -> Self {
        Self {
            token,
            cancel,
            out,
            last: Mutex::new((0, None)),
        }
    }

    pub fn cancelled(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }

    /// `notifications/progress`, only when the request carried a `progressToken`, and at most one
    /// per [`PROGRESS_INTERVAL`] (the first one always goes out).
    pub fn progress(&self, message: impl FnOnce() -> String) {
        let Some(token) = &self.token else { return };
        let mut last = self.last.lock().unwrap_or_else(|e| e.into_inner());
        if last.1.is_some_and(|t| t.elapsed() < PROGRESS_INTERVAL) {
            return;
        }
        last.0 += 1;
        last.1 = Some(Instant::now());
        self.out.send(&json!({
            "jsonrpc": "2.0",
            "method": "notifications/progress",
            "params": {"progressToken": token, "progress": last.0, "message": message()},
        }));
    }
}

impl Model {
    /// The name reported as `path` by the tools that take a loaded model.
    fn display(&self, id: &str) -> String {
        if self.edits_applied > 0 {
            format!("{id} (edited)")
        } else {
            self.source_path.clone()
        }
    }
}

pub struct Registry {
    models: BTreeMap<String, Model>,
    next_id: u32,
    dir: PathBuf,
}

impl Default for Registry {
    fn default() -> Self {
        Self::new()
    }
}

impl Registry {
    pub fn new() -> Self {
        Self {
            models: BTreeMap::new(),
            next_id: 1,
            dir: std::env::temp_dir().join(format!("rusolve-mcp-{}", std::process::id())),
        }
    }

    /// Removes the temp files of `mps_text` models.
    pub fn shutdown(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }

    pub fn has_tool(&self, name: &str) -> bool {
        TOOL_NAMES.contains(&name)
    }

    pub fn tool_list(&self) -> Vec<Value> {
        tool_defs()
    }

    pub fn call(&mut self, name: &str, args: &Args, ctx: &Arc<Ctx>) -> ToolResult {
        match name {
            "load_model" => self.load_model(args),
            "build_model" => self.build_model(args),
            "solve" => {
                let ctx = ctx.clone();
                self.with_model(args, move |id, m, a| solve(id, m, a, &ctx))
            }
            "get_solution" => self.with_model(args, get_solution),
            "modify" => self.with_model(args, linear_only("modify", modify)),
            "check_solution" => {
                self.with_model(args, linear_only("check_solution", check_solution))
            }
            "explain" => self.with_model(args, linear_only("explain", explain)),
            "lint" => self.with_model(args, linear_only("lint", lint)),
            "ranging" => self.with_model(args, linear_only("ranging", ranging)),
            "feasrelax" => self.with_model(args, linear_only("feasrelax", feasrelax)),
            "free_model" => self.free_model(args),
            other => Err(format!("unknown tool {other}")),
        }
    }

    fn with_model(
        &mut self,
        args: &Args,
        f: impl FnOnce(&str, &mut Model, &Args) -> ToolResult,
    ) -> ToolResult {
        let id = req_str(args, "model_id")?;
        let known: Vec<String> = self.models.keys().cloned().collect();
        let Some(m) = self.models.get_mut(&id) else {
            return Err(format!(
                "model_id: no model '{id}' (loaded: {}); call load_model first",
                if known.is_empty() {
                    "none".to_string()
                } else {
                    known.join(", ")
                }
            ));
        };
        f(&id, m, args)
    }

    fn load_model(&mut self, args: &Args) -> ToolResult {
        only_fields(args, &["path", "mps_text", "lp_text", "mof_text"])?;
        let path = opt_str(args, "path")?;
        let texts = [
            ("mps_text", opt_str(args, "mps_text")?, "mps"),
            ("lp_text", opt_str(args, "lp_text")?, "lp"),
            ("mof_text", opt_str(args, "mof_text")?, "mof.json"),
        ];
        let given: Vec<&str> = std::iter::once(("path", path.is_some()))
            .chain(texts.iter().map(|(n, t, _)| (*n, t.is_some())))
            .filter(|(_, g)| *g)
            .map(|(n, _)| n)
            .collect();
        match given.len() {
            0 => {
                return Err(
                    "path / mps_text / lp_text / mof_text: one is required (a model file path, or the model as MPS, LP or MathOptFormat text)".into(),
                );
            }
            1 => {}
            _ => {
                return Err(format!(
                    "{}: give exactly one of path, mps_text, lp_text, mof_text, not several",
                    given.join(" / ")
                ));
            }
        }
        if let Some(p) = path {
            return self.register(p, None, "path");
        }
        let (what, text, ext) = texts
            .into_iter()
            .find_map(|(n, t, e)| t.map(|t| (n, t, e)))
            .expect("one text given");
        self.load_text(&text, ext, what)
    }

    /// Writes `text` to a temp file of the right extension (the reader picks the format from it)
    /// and registers it.
    fn load_text(&mut self, text: &str, ext: &str, what: &str) -> ToolResult {
        std::fs::create_dir_all(&self.dir)
            .map_err(|e| format!("cannot create {}: {e}", self.dir.display()))?;
        let file = self.dir.join(format!("m{}.{ext}", self.next_id));
        std::fs::write(&file, text).map_err(|e| format!("cannot write {}: {e}", file.display()))?;
        self.register(file.to_string_lossy().into_owned(), Some(file), what)
    }

    /// A model an agent writes as JSON, through the same reader as a loaded one
    /// (`mcp_build.rs`).
    fn build_model(&mut self, args: &Args) -> ToolResult {
        only_fields(
            args,
            &["name", "sense", "variables", "rows", "objective", "cp"],
        )?;
        let mps = crate::mcp_build::to_mps(args)?;
        let cp = crate::mcp_build::parse_cp(args)?.filter(|k| !k.is_empty());
        let mut summary = self.load_text(&mps, "mps", "build_model")?;
        let Some(cp) = cp else {
            return Ok(summary);
        };
        // The objects are checked against the model as written: every column integer with finite
        // bounds, presence binary, sizes nonnegative. A bad one frees the handle again.
        let id = summary["model_id"].as_str().unwrap_or_default().to_string();
        let m = self.models.get_mut(&id).expect("just registered");
        if let Err(e) = crate::cp_solve::validate(&m.highs.model_lp_, &cp) {
            let at = crate::mcp_build::cp_error_path(&cp, &e.object).unwrap_or_else(|| "cp".into());
            let msg = format!("{at}: RS4011: {e}");
            if let Some(t) = self.models.remove(&id).and_then(|m| m.temp) {
                let _ = std::fs::remove_file(t);
            }
            return Err(msg);
        }
        let engine = crate::cp_solve::select_engine(&cp, &m.highs.options_).name();
        let counts: serde_json::Map<String, Value> = cp
            .counts()
            .iter()
            .filter(|(_, n)| *n > 0)
            .map(|(k, n)| (k.to_string(), json!(n)))
            .collect();
        m.cp = Some(cp);
        summary["class"] = json!("cp");
        summary["engine"] = json!(engine);
        summary["cp_objects"] = Value::Object(counts);
        Ok(summary)
    }

    /// Reads `source_path`, classifies it (router), and keeps a handle: a live `Highs` for an
    /// LP or MIP, a routed handle for a class another engine solves, a typed error for a class
    /// nothing solves.
    fn register(&mut self, source_path: String, temp: Option<PathBuf>, what: &str) -> ToolResult {
        let id = format!("m{}", self.next_id);
        let fail = |temp: &Option<PathBuf>, msg: String| {
            if let Some(t) = temp {
                let _ = std::fs::remove_file(t);
            }
            msg
        };
        let shown = |e: String| {
            if temp.is_some() {
                format!("{what}: {}", e.replace(&source_path, &format!("<{what}>")))
            } else {
                format!("{what}: {e}")
            }
        };
        if rusolve_io::may_have_special_content(&source_path)
            && let Ok(model) = rusolve_io::read_model_file_conic(&source_path)
        {
            use crate::route::Class;
            let verdict = crate::route::classify_model(&model);
            let class = verdict.class;
            if !matches!(class, Class::Lp | Class::Mip) {
                // A convex MIQP and a mixed-integer conic model load as routed handles: their
                // engine (`oa`) runs only when `solve` is given the lever in `options`.
                let oa_class = matches!(class, Class::MiConic)
                    || (class == Class::Miqp && verdict.convex_miqp);
                let Some(engine) = class.engine().or(oa_class.then_some("oa")) else {
                    return Err(fail(
                        &temp,
                        format!(
                            "{what}: {} models are not supported: {}",
                            class.name(),
                            verdict.reason
                        ),
                    ));
                };
                let summary = json!({
                    "model_id": id,
                    "name": model.name,
                    "rows": model.num_constraints,
                    "cols": model.num_variables,
                    "integer_cols": model.bounds.iter().filter(|b| b.integer).count(),
                    "sense": if model.obj_sense == rusolve_io::ObjSense::Maximize { "maximize" } else { "minimize" },
                    "names": {"rows": names_summary(&model.row_names), "cols": names_summary(&model.col_names)},
                    "class": class.name(),
                    "engine": engine,
                });
                self.next_id += 1;
                self.models.insert(
                    id,
                    Model {
                        source_path,
                        temp,
                        highs: Highs::with_options(solve_options()),
                        edits_applied: 0,
                        solved: false,
                        routed: Some(Routed {
                            class,
                            columns: None,
                            status: None,
                            objective: None,
                        }),
                        cp: None,
                        cp_result: None,
                    },
                );
                return Ok(summary);
            }
        }
        let lp = match read_mip(&source_path) {
            Ok(lp) => lp,
            Err(e) => return Err(fail(&temp, shown(e))),
        };
        let class = if count_integer(&lp) > 0 { "mip" } else { "lp" };
        let summary = json!({
            "model_id": id,
            "name": lp.model_name_,
            "rows": lp.num_row_,
            "cols": lp.num_col_,
            "integer_cols": count_integer(&lp),
            "nonzeros": lp.a_matrix_.value_.len(),
            "sense": if lp.sense_ == ObjSense::kMinimize { "minimize" } else { "maximize" },
            "names": {"rows": names_summary(&lp.row_names_), "cols": names_summary(&lp.col_names_)},
            "class": class,
            "engine": "highs-port",
        });
        let mut highs = Highs::with_options(solve_options());
        if highs.pass_model(lp) == HighsStatus::kError {
            return Err(fail(
                &temp,
                "internal error: Highs::pass_model rejected the model".into(),
            ));
        }
        self.next_id += 1;
        self.models.insert(
            id,
            Model {
                source_path,
                temp,
                highs,
                edits_applied: 0,
                solved: false,
                routed: None,
                cp: None,
                cp_result: None,
            },
        );
        Ok(summary)
    }

    fn free_model(&mut self, args: &Args) -> ToolResult {
        only_fields(args, &["model_id"])?;
        let id = req_str(args, "model_id")?;
        match self.models.remove(&id) {
            Some(m) => {
                if let Some(t) = m.temp {
                    let _ = std::fs::remove_file(t);
                }
                Ok(json!({"model_id": id, "freed": true, "remaining": self.models.len()}))
            }
            None => Err(format!("model_id: no model '{id}' to free")),
        }
    }
}

/// The options a handle's `Highs` starts from: the CLI's `trace_options`
/// (mip_highs.rs): serial, seed 0, no logging, over `base_options`. A log on
/// stdout would corrupt the protocol, so `log_to_console` is also forced off
/// again after every option the caller sets.
fn solve_options() -> HighsOptions {
    let mut o = HighsOptions {
        output_flag: false,
        threads: 1,
        random_seed: 0,
        ..crate::base_options()
    };
    o.sync_log_options();
    o
}

fn count_integer(lp: &HighsLp) -> usize {
    lp.integrality_
        .iter()
        .filter(|t| **t == HighsVarType::kInteger)
        .count()
}

fn names_summary(names: &[String]) -> Value {
    json!({"count": names.len(), "first": names.iter().take(5).collect::<Vec<_>>()})
}

// ------------------------------------------------------------ argument access

fn type_name(v: &Value) -> &'static str {
    match v {
        Value::Null => "null",
        Value::Bool(_) => "a boolean",
        Value::Number(_) => "a number",
        Value::String(_) => "a string",
        Value::Array(_) => "an array",
        Value::Object(_) => "an object",
    }
}

/// Rejects an argument the tool does not have, naming it (the schema says
/// `additionalProperties: false`).
fn only_fields(args: &Args, allowed: &[&str]) -> Result<(), String> {
    for k in args.keys() {
        if !allowed.contains(&k.as_str()) {
            return Err(format!(
                "{k}: unknown argument (this tool takes: {})",
                allowed.join(", ")
            ));
        }
    }
    Ok(())
}

fn opt_str(a: &Args, f: &str) -> Result<Option<String>, String> {
    match a.get(f) {
        None => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.clone())),
        Some(v) => Err(format!("{f}: expected a string, got {}", type_name(v))),
    }
}

fn req_str(a: &Args, f: &str) -> Result<String, String> {
    opt_str(a, f)?.ok_or_else(|| format!("{f}: required"))
}

fn opt_f64(a: &Args, f: &str) -> Result<Option<f64>, String> {
    match a.get(f) {
        None => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(v) => Err(format!("{f}: expected a number, got {}", type_name(v))),
    }
}

fn opt_usize(a: &Args, f: &str) -> Result<Option<usize>, String> {
    match a.get(f) {
        None => Ok(None),
        Some(Value::Number(n)) => n
            .as_u64()
            .map(|x| Some(x as usize))
            .ok_or_else(|| format!("{f}: expected a non-negative integer, got {n}")),
        Some(v) => Err(format!(
            "{f}: expected a non-negative integer, got {}",
            type_name(v)
        )),
    }
}

fn opt_bool(a: &Args, f: &str) -> Result<Option<bool>, String> {
    match a.get(f) {
        None => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(v) => Err(format!("{f}: expected a boolean, got {}", type_name(v))),
    }
}

fn opt_str_list(a: &Args, f: &str) -> Result<Vec<String>, String> {
    match a.get(f) {
        None => Ok(Vec::new()),
        Some(Value::Array(items)) => items
            .iter()
            .enumerate()
            .map(|(i, v)| {
                v.as_str()
                    .map(str::to_string)
                    .ok_or_else(|| format!("{f}[{i}]: expected a string, got {}", type_name(v)))
            })
            .collect(),
        Some(v) => Err(format!(
            "{f}: expected an array of strings, got {}",
            type_name(v)
        )),
    }
}

fn opt_object<'a>(a: &'a Args, f: &str) -> Result<Option<&'a Map<String, Value>>, String> {
    match a.get(f) {
        None => Ok(None),
        Some(Value::Object(m)) => Ok(Some(m)),
        Some(v) => Err(format!("{f}: expected an object, got {}", type_name(v))),
    }
}

// --------------------------------------------------------------------- solve

/// Load a Wasm hook component (`docs/design/wasm-hooks.md`) for one solve: a tool error naming
/// the file on any load failure (unreadable, invalid module, an import, a bad export
/// signature), never a crash. Returns the component name, its counters handle and the hooks.
fn load_wasm_hook(
    path: &str,
    fuel: Option<u64>,
    memory: Option<usize>,
) -> Result<
    (
        String,
        rusolve_wasm::WasmStatsHandle,
        rusolve_wasm::WasmHooks,
    ),
    String,
> {
    let name = std::path::Path::new(path)
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string());
    let bytes = std::fs::read(path).map_err(|e| format!("wasm_hook: cannot read {path}: {e}"))?;
    let mut limits = rusolve_wasm::WasmLimits::default();
    if let Some(f) = fuel {
        limits.fuel_per_call = f;
    }
    if let Some(m) = memory {
        limits.max_memory_bytes = m;
    }
    let hooks = rusolve_wasm::WasmHooks::new(&name, &bytes, limits)
        .map_err(|e| format!("wasm_hook: {path}: {e}"))?;
    let stats = hooks.stats();
    Ok((name, stats, hooks))
}

/// `options` of `solve`: a flat `{name: value}` map through the ported
/// `Highs::set_option_value_*` setters, picked by the JSON type of the value.
// The tool's `options` schema lists the registry names (`commands::option_properties`);
// a bad value is still reported by the setter's refusal.
fn apply_option(h: &mut Highs, key: &str, v: &Value) -> Result<String, String> {
    let bad = || format!("options.{key}: unknown option or invalid value {v}");
    let status = match v {
        Value::Bool(b) => {
            let s = h.set_option_value_bool(key, *b);
            if s == HighsStatus::kError {
                h.set_option_value_str(key, if *b { "true" } else { "false" })
            } else {
                s
            }
        }
        Value::Number(n) => {
            let int = n.as_i64().and_then(|i| HighsInt::try_from(i).ok());
            match (int, n.as_f64()) {
                (Some(i), _) => {
                    let s = h.set_option_value_int(key, i);
                    if s == HighsStatus::kError {
                        h.set_option_value_double(key, i as f64)
                    } else {
                        s
                    }
                }
                (None, Some(x)) => h.set_option_value_double(key, x),
                _ => return Err(bad()),
            }
        }
        Value::String(s) => h.set_option_value_str(key, s),
        _ => {
            return Err(format!(
                "options.{key}: expected a string, number or boolean, got {}",
                type_name(v)
            ));
        }
    };
    if status == HighsStatus::kError {
        return Err(bad());
    }
    Ok(match v {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    })
}

/// Wraps a tool that needs the live LP/MIP of the handle: on a routed handle (a QP, conic or
/// general-constraint model) it fails with the tools that do work there.
fn linear_only(
    tool: &'static str,
    f: fn(&str, &mut Model, &Args) -> ToolResult,
) -> impl FnOnce(&str, &mut Model, &Args) -> ToolResult {
    move |id, m, args| {
        if let Some(r) = &m.routed {
            return Err(format!(
                "model_id: '{id}' is a {} model, solved by the {} engine; {tool} works on LP and MIP models only (solve, get_solution and free_model work on it)",
                r.class.name(),
                r.class.engine().unwrap_or("router"),
            ));
        }
        // a model with CP objects: only `check_solution` evaluates them; the others would work
        // on the plain rows and silently ignore the objects
        if tool != "check_solution"
            && let Some(o) = m.cp.as_ref().and_then(|k| k.first_object())
        {
            return Err(format!(
                "model_id: RS4012: '{id}' has scheduling or combinatorial objects (first: {o}); {tool} would ignore them. solve, get_solution, check_solution and free_model support them"
            ));
        }
        f(id, m, args)
    }
}

/// Progress of the interrupt callbacks of the port while a handle solves: cancel through
/// `user_interrupt` (HiGHS `HighsCallbackInput`, `kCallbackMipInterrupt` and
/// `kCallbackSimplexInterrupt`), and a progress notification of the node or iteration count.
fn install_callbacks(h: &mut Highs, ctx: &Arc<Ctx>) {
    let ctx = ctx.clone();
    h.set_callback(HighsCallbackFunctionType::new(
        move |kind: HighsInt,
              _msg: &str,
              out: &HighsCallbackOutput,
              cb_in: &mut HighsCallbackInput| {
            if ctx.cancelled() {
                cb_in.user_interrupt = true;
                return;
            }
            if kind == kCallbackMipInterrupt {
                ctx.progress(|| {
                    format!(
                        "MIP: {} nodes, {} LP iterations, dual bound {}, gap {}",
                        out.mip_node_count,
                        out.mip_total_lp_iterations,
                        fmt_num(out.mip_dual_bound),
                        fmt_num(out.mip_gap),
                    )
                });
            } else if kind == kCallbackSimplexInterrupt {
                ctx.progress(|| format!("simplex: {} iterations", out.simplex_iteration_count));
            }
        },
    ));
    h.start_callback(kCallbackMipInterrupt);
    h.start_callback(kCallbackSimplexInterrupt);
}

fn fmt_num(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.6e}")
    } else {
        "none".into()
    }
}

/// `solve` on a routed handle: the CLI's own `solve` in a child process.
fn solve_routed(_id: &str, m: &mut Model, args: &Args, ctx: &Ctx) -> ToolResult {
    use std::io::Read;
    use std::process::{Command, Stdio};
    let class = m.routed.as_ref().expect("routed handle").class;
    // `options` on a convex MIQP or mixed-integer conic handle: only the outer-approximation
    // levers (`miqp_oa*`, `misocp_oa*`), handed to the child as RUSOLVE_EKK_OPTS.
    let oa_class = matches!(
        class,
        crate::route::Class::Miqp | crate::route::Class::MiConic
    );
    let mut oa_opts: Vec<String> = Vec::new();
    if oa_class && let Some(o) = args.get("options") {
        let map = o
            .as_object()
            .ok_or_else(|| "options: must be an object {option name: value}".to_string())?;
        for (k, v) in map {
            if !(k.starts_with("miqp_oa") || k.starts_with("misocp_oa")) {
                return Err(format!(
                    "options.{k}: a {} model takes only the miqp_oa* / misocp_oa* options here",
                    class.name()
                ));
            }
            let val = match v {
                Value::Bool(b) => b.to_string(),
                Value::Number(n) => n.to_string(),
                Value::String(s) if !s.contains(',') => s.clone(),
                _ => return Err(format!("options.{k}: must be a boolean, number or string")),
            };
            oa_opts.push(format!("{k}={val}"));
        }
    }
    for k in [
        "options",
        "warm_start",
        "wasm_hook",
        "wasm_fuel",
        "wasm_memory",
    ] {
        if args.contains_key(k) && !(k == "options" && oa_class) {
            return Err(format!(
                "{k}: not supported on a {} model (the {} engine takes time_limit and include_solution)",
                class.name(),
                class.engine().unwrap_or("routed"),
            ));
        }
    }
    let time_limit = opt_f64(args, "time_limit")?;
    if let Some(t) = time_limit
        && t < 0.0
    {
        return Err(format!("time_limit: must be >= 0, got {t}"));
    }
    let include_solution = opt_bool(args, "include_solution")?.unwrap_or(false);
    let exe =
        std::env::current_exe().map_err(|e| format!("internal error: no executable path: {e}"))?;
    let mut cmd = Command::new(exe);
    cmd.arg("solve")
        .arg(&m.source_path)
        .args(["--format", "json", "--include-solution"]);
    if let Some(t) = time_limit {
        cmd.args(["--time-limit", &t.to_string()]);
    }
    if oa_class {
        // an explicit `options` replaces the inherited setting; none leaves the environment alone
        if !oa_opts.is_empty() {
            cmd.env("RUSOLVE_EKK_OPTS", oa_opts.join(","));
        }
    }
    let mut child = cmd
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| {
            format!(
                "internal error: cannot start the {} engine: {e}",
                class.name()
            )
        })?;
    // drain both pipes on threads so a large solution cannot block the child
    let mut so = child.stdout.take().expect("piped");
    let mut se = child.stderr.take().expect("piped");
    let out_t = std::thread::spawn(move || {
        let mut b = String::new();
        let _ = so.read_to_string(&mut b);
        b
    });
    let err_t = std::thread::spawn(move || {
        let mut b = String::new();
        let _ = se.read_to_string(&mut b);
        b
    });
    let start = Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(st)) => break st,
            Ok(None) => {}
            Err(e) => return Err(format!("internal error: waiting for the engine: {e}")),
        }
        if ctx.cancelled() {
            let _ = child.kill();
            let _ = child.wait();
            m.solved = false;
            return Err("cancelled".into());
        }
        ctx.progress(|| {
            format!(
                "{} engine running, {:.1} s",
                class.engine().unwrap_or("routed"),
                start.elapsed().as_secs_f64()
            )
        });
        std::thread::sleep(Duration::from_millis(10));
    };
    let stdout = out_t.join().unwrap_or_default();
    let stderr = err_t.join().unwrap_or_default();
    let mut env: Value = serde_json::from_str(stdout.trim()).map_err(|_| {
        format!(
            "the {} engine exited with {status} and no JSON result: {}",
            class.name(),
            stderr.trim().lines().last().unwrap_or("(no message)")
        )
    })?;
    if env.get("schema_version").and_then(Value::as_str) == Some("rusolve.error/1") {
        return Err(env["error"]["message"]
            .as_str()
            .unwrap_or("the engine refused the model")
            .to_string());
    }
    let columns = env
        .get("solution")
        .and_then(|s| s.get("columns"))
        .and_then(Value::as_array)
        .cloned();
    if let Some(r) = m.routed.as_mut() {
        r.status = env
            .get("status")
            .and_then(Value::as_str)
            .map(str::to_string);
        r.objective = env.get("objective").and_then(Value::as_f64);
        r.columns = columns;
    }
    m.solved = m.routed.as_ref().is_some_and(|r| r.columns.is_some());
    if !include_solution && let Some(o) = env.as_object_mut() {
        o.remove("solution");
    }
    // The QP and conic engines print the raw HiGHS conventions (no tree: nodes -1, no thread
    // count: 0); the result schema (minimum 1 and 0) is clamped as on the LP path.
    if env["stats"]["nodes"].as_i64().is_some_and(|n| n < 0) {
        env["stats"]["nodes"] = json!(0);
    }
    if env["run"]["threads"].as_i64().is_some_and(|n| n < 1) {
        env["run"]["threads"] = json!(1);
    }
    Ok(env)
}

/// `solve` on a handle with CP objects (`cp_solve.rs`): in this process, the answer checked on the
/// user's objects. `options` and `time_limit` go to the HiGHS port that solves the linearization.
fn solve_cp(_id: &str, m: &mut Model, args: &Args, ctx: &Arc<Ctx>) -> ToolResult {
    for k in ["warm_start", "wasm_hook", "wasm_fuel", "wasm_memory"] {
        if args.contains_key(k) {
            return Err(format!(
                "{k}: not supported on a model with CP objects (solve takes options, time_limit and include_solution)"
            ));
        }
    }
    let mut settings: Vec<(String, String)> = Vec::new();
    if let Some(opts) = opt_object(args, "options")? {
        for (k, v) in opts {
            let shown = apply_option(&mut m.highs, k, v)?;
            settings.push((k.clone(), shown));
        }
    }
    if let Some(t) = opt_f64(args, "time_limit")? {
        if t < 0.0 {
            return Err(format!("time_limit: must be >= 0, got {t}"));
        }
        m.highs.set_option_value_double("time_limit", t);
        settings.push(("time_limit".into(), t.to_string()));
    }
    let include_solution = opt_bool(args, "include_solution")?.unwrap_or(false);
    m.highs.set_option_value_bool("log_to_console", false);
    let options = m.highs.get_options().clone();
    let user_lp = m.highs.model_lp_.clone();
    let cp = m.cp.clone().expect("cp handle");
    let start = Instant::now();
    m.solved = false;
    m.cp_result = None;
    let run_ctx = ctx.clone();
    let solved = crate::cp_solve::solve_with(user_lp, cp, options, move |h| {
        install_callbacks(h, &run_ctx);
        let r = h.run();
        h.stop_callback(kCallbackMipInterrupt);
        h.stop_callback(kCallbackSimplexInterrupt);
        r
    })
    .map_err(|e| {
        let d = e.diagnostic();
        format!("{}: {e}", d.code)
    })?;
    let time_s = start.elapsed().as_secs_f64();
    m.solved = solved.run != HighsStatus::kError;
    let h = &solved.highs;
    let opts = h.get_options();
    let env_in = envelope::Solved {
        path: &m.source_path,
        lp: &solved.user_lp,
        status: h.get_model_status(),
        run_status: solved.run,
        info: h.get_info(),
        solution: &h.solution_,
        time_s,
        seed: opts.random_seed,
        threads: opts.threads,
        time_limit: opts.time_limit,
        node_limit: None,
        mip_max_nodes: opts.mip_max_nodes,
        options: &settings,
        include_solution,
        mip_stats: h.mip_stats_.as_ref(),
        wasm: None,
    };
    let mut env = serde_json::to_value(envelope::envelope(&env_in))
        .map_err(|e| format!("internal error: cannot serialise the result: {e}"))?;
    if env["stats"]["nodes"].as_i64().is_some_and(|n| n < 0) {
        env["stats"]["nodes"] = json!(0);
    }
    env["class"] = json!("cp");
    env["engine"] = json!(solved.engine);
    env["cp"] = solved.report.clone();
    if let Some(cps) = solved.report.get("cpsat") {
        env["cpsat"] = cps.clone();
    }
    let input = env["input"].as_object_mut().expect("envelope.input");
    input.insert("edits_applied".into(), json!(m.edits_applied));
    m.cp_result = Some(CpResult {
        highs: solved.highs,
    });
    Ok(env)
}

fn solve(id: &str, m: &mut Model, args: &Args, ctx: &Arc<Ctx>) -> ToolResult {
    only_fields(
        args,
        &[
            "model_id",
            "options",
            "time_limit",
            "warm_start",
            "include_solution",
            "wasm_hook",
            "wasm_fuel",
            "wasm_memory",
        ],
    )?;
    if m.routed.is_some() {
        return solve_routed(id, m, args, ctx);
    }
    if m.cp.is_some() {
        return solve_cp(id, m, args, ctx);
    }
    let mut settings: Vec<(String, String)> = Vec::new();
    if let Some(opts) = opt_object(args, "options")? {
        for (k, v) in opts {
            let shown = apply_option(&mut m.highs, k, v)?;
            settings.push((k.clone(), shown));
        }
    }
    if let Some(t) = opt_f64(args, "time_limit")? {
        if t < 0.0 {
            return Err(format!("time_limit: must be >= 0, got {t}"));
        }
        m.highs.set_option_value_double("time_limit", t);
        settings.push(("time_limit".into(), t.to_string()));
    }
    let include_solution = opt_bool(args, "include_solution")?.unwrap_or(false);
    let warm = opt_bool(args, "warm_start")?.unwrap_or(true);
    // `options: {"mip_engine": "cpsat"}` (docs/design/cpsat-engine.md) replaces the solve below
    // with CP-SAT; a Wasm hook is a tool error, and so is a missing OR-Tools library (RS0012)
    let cpsat = crate::cpsat_solve::selected(&m.highs).map_err(|e| format!("options: {e}"))?;
    if cpsat && args.contains_key("wasm_hook") {
        return Err("wasm_hook: not supported with mip_engine=cpsat".into());
    }
    // Keep stdout for the protocol whatever the options said.
    m.highs.set_option_value_bool("log_to_console", false);
    if !warm {
        m.highs.clear_solver();
    }
    // `mip_stats_` is only written by a MIP solve that finishes; clear it so a failed or
    // LP re-solve of this handle does not report the previous solve's statistics
    m.highs.mip_stats_ = None;
    // A Wasm hook component for this solve only — a load failure is a tool error naming the
    // file, and the hooks are cleared afterwards, so a `solve` without `wasm_hook` stays the
    // unmodified port search (identity 20/20).
    let mut wasm_report = None;
    let wasm_hook = opt_str(args, "wasm_hook")?;
    let wasm_fuel = opt_usize(args, "wasm_fuel")?.map(|f| f as u64);
    let wasm_memory = opt_usize(args, "wasm_memory")?;
    if wasm_hook.is_none() && wasm_fuel.is_some() {
        return Err("wasm_fuel: needs wasm_hook (it configures the component)".into());
    }
    if wasm_hook.is_none() && wasm_memory.is_some() {
        return Err("wasm_memory: needs wasm_hook (it configures the component)".into());
    }
    if let Some(p) = &wasm_hook {
        let (name, stats, hooks) = load_wasm_hook(p, wasm_fuel, wasm_memory)?;
        m.highs.set_mip_hooks(Box::new(hooks));
        wasm_report = Some((name, stats));
    }
    let start = std::time::Instant::now();
    // a `threads` change between two solves on this server's thread
    m.highs.reset_scheduler_for_threads();
    // Cancel and progress ride on the port's interrupt callbacks for this solve only; they
    // read, never steer, the search (no callback means the unmodified run).
    install_callbacks(&mut m.highs, ctx);
    let mut cpsat_report = None;
    let run_status = if cpsat {
        match crate::cpsat_solve::run(&mut m.highs, &m.source_path) {
            Ok(report) => {
                cpsat_report = Some(report);
                HighsStatus::kOk
            }
            Err(f) => {
                m.highs.stop_callback(kCallbackMipInterrupt);
                m.highs.stop_callback(kCallbackSimplexInterrupt);
                return Err(f.message());
            }
        }
    } else {
        m.highs.run()
    };
    m.highs.stop_callback(kCallbackMipInterrupt);
    m.highs.stop_callback(kCallbackSimplexInterrupt);
    m.highs.clear_mip_hooks();
    let time_s = start.elapsed().as_secs_f64();
    m.solved = run_status != HighsStatus::kError;
    let h = &m.highs;
    let opts = h.get_options();
    let solved = envelope::Solved {
        path: &m.source_path,
        lp: &h.model_lp_,
        status: h.get_model_status(),
        run_status,
        info: h.get_info(),
        solution: &h.solution_,
        time_s,
        seed: opts.random_seed,
        threads: opts.threads,
        time_limit: opts.time_limit,
        node_limit: None,
        mip_max_nodes: opts.mip_max_nodes,
        options: &settings,
        include_solution,
        mip_stats: h.mip_stats_.as_ref(),
        // The component's counters, only when wasm_hook loaded one for this solve
        wasm: wasm_report.map(|(name, st)| {
            let s = st.lock().unwrap_or_else(|e| e.into_inner()).clone();
            envelope::wasm_stats_json(&name, &s)
        }),
    };
    let mut env = serde_json::to_value(envelope::envelope(&solved))
        .map_err(|e| format!("internal error: cannot serialise the result: {e}"))?;
    // An LP has no tree: HiGHS leaves mip_node_count at -1, which the schema
    // (minimum 0) rejects. The CLI prints the raw -1; the MCP result is
    // clamped so structuredContent validates against its outputSchema.
    if env["stats"]["nodes"].as_i64().is_some_and(|n| n < 0) {
        env["stats"]["nodes"] = json!(0);
    }
    if let Some(report) = cpsat_report {
        env["engine"] = json!(crate::cpsat_solve::ENGINE);
        env["cpsat"] = report;
    }
    let input = env["input"].as_object_mut().expect("envelope.input");
    input.insert("edits_applied".into(), json!(m.edits_applied));
    if m.edits_applied > 0 {
        input.insert("path".into(), json!(m.display(id)));
        input.insert("sha256".into(), json!(model_digest(&h.model_lp_)));
    }
    Ok(env)
}

/// SHA-256 of a canonical dump of the model, for a handle whose model is no
/// longer any file (`input.sha256` still identifies what was solved).
fn model_digest(lp: &HighsLp) -> String {
    use sha2::{Digest, Sha256};
    let dump = format!(
        "{:?}|{:?}|{}|{}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}",
        lp.sense_,
        lp.offset_,
        lp.num_col_,
        lp.num_row_,
        lp.col_cost_,
        lp.col_lower_,
        lp.col_upper_,
        lp.row_lower_,
        lp.row_upper_,
        lp.integrality_,
        lp.a_matrix_.start_,
        lp.a_matrix_.index_,
        lp.a_matrix_.value_,
        lp.col_names_,
        lp.row_names_,
    );
    Sha256::digest(dump.as_bytes())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn get_solution(id: &str, m: &mut Model, args: &Args) -> ToolResult {
    only_fields(args, &["model_id", "names", "offset", "limit"])?;
    let names = opt_str_list(args, "names")?;
    let offset = opt_usize(args, "offset")?.unwrap_or(0);
    let limit = opt_usize(args, "limit")?.unwrap_or(DEFAULT_PAGE);
    if let Some(r) = &m.routed {
        return get_solution_routed(id, r, m.solved, &names, offset, limit);
    }
    // a CP handle: the solution is the solver's after the CP solve (cut to the user's columns);
    // names and counts come from the user's model
    let h = m.cp_result.as_ref().map_or(&m.highs, |r| &r.highs);
    if !m.solved || !h.solution_.value_valid {
        return Err(format!(
            "model_id: model '{id}' has no solution (never solved, edited since the last solve, or the last solve found no primal solution); call solve first"
        ));
    }
    let lp = &m.highs.model_lp_;
    let col_name = |j: usize| {
        lp.col_names_
            .get(j)
            .cloned()
            .unwrap_or_else(|| format!("c{j}"))
    };
    let value = |j: usize| h.solution_.col_value.get(j).copied();
    let finite = |x: Option<f64>| x.filter(|v| v.is_finite());
    let total = lp.num_col_ as usize;
    let (columns, truncated, count) = if names.is_empty() {
        let end = (offset + limit).min(total);
        let cols: Vec<Value> = (offset.min(total)..end)
            .map(|j| json!({"name": col_name(j), "value": finite(value(j))}))
            .collect();
        (cols, end < total, total)
    } else {
        let index: std::collections::HashMap<String, usize> =
            (0..total).map(|j| (col_name(j), j)).collect();
        let mut cols = Vec::new();
        for (k, n) in names.iter().enumerate() {
            let Some(&j) = index.get(n) else {
                return Err(format!("names[{k}]: no column named '{n}' in model '{id}'"));
            };
            cols.push(json!({"name": n, "value": finite(value(j))}));
        }
        let n = cols.len();
        (cols, false, n)
    };
    let info = h.get_info();
    Ok(json!({
        "model_id": id,
        "status": envelope::result_status(h.get_model_status()),
        "objective": finite(Some(info.objective_function_value)),
        "total_columns": total,
        "returned": count.min(columns.len()),
        "offset": if names.is_empty() { offset } else { 0 },
        "truncated": truncated,
        "columns": columns,
    }))
}

/// `get_solution` for a routed handle: pages the `solution.columns` the child returned.
fn get_solution_routed(
    id: &str,
    r: &Routed,
    solved: bool,
    names: &[String],
    offset: usize,
    limit: usize,
) -> ToolResult {
    let cols = match (&r.columns, solved) {
        (Some(cs), true) => cs,
        _ => {
            return Err(format!(
                "model_id: model '{id}' has no solution (never solved, or the last solve found no primal solution); call solve first"
            ));
        }
    };
    let total = cols.len();
    let pick = |col: &Value| json!({"name": col["name"], "value": col["value"]});
    let (columns, truncated) = if names.is_empty() {
        let end = (offset + limit).min(total);
        let page: Vec<Value> = cols[offset.min(total)..end].iter().map(pick).collect();
        (page, end < total)
    } else {
        let mut page = Vec::new();
        for (k, n) in names.iter().enumerate() {
            let Some(col) = cols.iter().find(|x| x["name"].as_str() == Some(n)) else {
                return Err(format!("names[{k}]: no column named '{n}' in model '{id}'"));
            };
            page.push(pick(col));
        }
        (page, false)
    };
    Ok(json!({
        "model_id": id,
        "status": r.status,
        "objective": r.objective,
        "total_columns": total,
        "returned": columns.len(),
        "offset": if names.is_empty() { offset } else { 0 },
        "truncated": truncated,
        "columns": columns,
    }))
}

// -------------------------------------------------------------------- modify

fn find_index(names: &[String], what: &str, name: &str, at: &str) -> Result<HighsInt, String> {
    names
        .iter()
        .position(|n| n == name)
        .map(|i| i as HighsInt)
        .ok_or_else(|| format!("{at}: no {what} named '{name}'"))
}

/// A bound field of an edit: absent keeps the current value, `null` is
/// infinite on that side, a number sets it.
fn bound_edit(e: &Args, f: &str, current: f64, inf: f64, at: &str) -> Result<f64, String> {
    match e.get(f) {
        None => Ok(current),
        Some(Value::Null) => Ok(inf),
        Some(Value::Number(n)) => Ok(n.as_f64().unwrap_or(current)),
        Some(v) => Err(format!(
            "{at}.{f}: expected a number or null (infinite), got {}",
            type_name(v)
        )),
    }
}

fn edit_num(e: &Args, f: &str, at: &str) -> Result<f64, String> {
    match e.get(f) {
        Some(Value::Number(n)) => n
            .as_f64()
            .ok_or_else(|| format!("{at}.{f}: not a finite number")),
        Some(v) => Err(format!("{at}.{f}: expected a number, got {}", type_name(v))),
        None => Err(format!("{at}.{f}: required")),
    }
}

fn edit_str(e: &Args, f: &str, at: &str) -> Result<String, String> {
    match e.get(f) {
        Some(Value::String(s)) if !s.is_empty() => Ok(s.clone()),
        Some(Value::String(_)) => Err(format!("{at}.{f}: must not be empty")),
        Some(v) => Err(format!("{at}.{f}: expected a string, got {}", type_name(v))),
        None => Err(format!("{at}.{f}: required")),
    }
}

/// `coefs`: `{name: value}` -> sorted parallel index and value vectors over `names`.
fn coef_vectors(
    e: &Args,
    names: &[String],
    what: &str,
    at: &str,
) -> Result<(Vec<HighsInt>, Vec<f64>), String> {
    let mut pairs: Vec<(HighsInt, f64)> = Vec::new();
    match e.get("coefs") {
        None => {}
        Some(Value::Object(m)) => {
            for (n, v) in m {
                let i = find_index(names, what, n, &format!("{at}.coefs"))?;
                let x = v.as_f64().ok_or_else(|| {
                    format!("{at}.coefs.{n}: expected a number, got {}", type_name(v))
                })?;
                pairs.push((i, x));
            }
        }
        Some(v) => {
            return Err(format!(
                "{at}.coefs: expected an object {{{what} name: coefficient}}, got {}",
                type_name(v)
            ));
        }
    }
    pairs.sort_by_key(|p| p.0);
    Ok(pairs.into_iter().unzip())
}

fn rejected(kind: &str, at: &str) -> String {
    format!("{at}: HiGHS rejected the {kind} edit")
}

fn apply_edit(h: &mut Highs, e: &Args, at: &str) -> Result<(), String> {
    let kind = edit_str(e, "kind", at)?;
    let ok = |s: HighsStatus, kind: &str| {
        if s == HighsStatus::kError {
            Err(rejected(kind, at))
        } else {
            Ok(())
        }
    };
    let allowed: &[&str] = match kind.as_str() {
        "col_bounds" | "row_bounds" => &["kind", "name", "lower", "upper"],
        "cost" => &["kind", "name", "cost"],
        "coef" => &["kind", "row", "col", "value"],
        "add_row" => &["kind", "name", "lower", "upper", "coefs"],
        "add_col" => &["kind", "name", "cost", "lower", "upper", "integer", "coefs"],
        "del_row" | "del_col" => &["kind", "name"],
        other => {
            return Err(format!(
                "{at}.kind: unknown edit kind '{other}' (col_bounds, row_bounds, cost, coef, add_row, add_col, del_row, del_col)"
            ));
        }
    };
    for k in e.keys() {
        if !allowed.contains(&k.as_str()) {
            return Err(format!(
                "{at}.{k}: not a field of a {kind} edit (it takes: {})",
                allowed.join(", ")
            ));
        }
    }
    let inf = kHighsInf;
    match kind.as_str() {
        "col_bounds" => {
            let name = edit_str(e, "name", at)?;
            let j = find_index(
                &h.model_lp_.col_names_,
                "column",
                &name,
                &format!("{at}.name"),
            )?;
            let (lo, up) = (
                h.model_lp_.col_lower_[j as usize],
                h.model_lp_.col_upper_[j as usize],
            );
            let lo = bound_edit(e, "lower", lo, -inf, at)?;
            let up = bound_edit(e, "upper", up, inf, at)?;
            ok(h.change_col_bounds(j, lo, up), &kind)
        }
        "row_bounds" => {
            let name = edit_str(e, "name", at)?;
            let i = find_index(&h.model_lp_.row_names_, "row", &name, &format!("{at}.name"))?;
            let (lo, up) = (
                h.model_lp_.row_lower_[i as usize],
                h.model_lp_.row_upper_[i as usize],
            );
            let lo = bound_edit(e, "lower", lo, -inf, at)?;
            let up = bound_edit(e, "upper", up, inf, at)?;
            ok(h.change_row_bounds(i, lo, up), &kind)
        }
        "cost" => {
            let name = edit_str(e, "name", at)?;
            let j = find_index(
                &h.model_lp_.col_names_,
                "column",
                &name,
                &format!("{at}.name"),
            )?;
            let cost = edit_num(e, "cost", at)?;
            ok(h.change_col_cost(j, cost), &kind)
        }
        "coef" => {
            let r = edit_str(e, "row", at)?;
            let cname = edit_str(e, "col", at)?;
            let i = find_index(&h.model_lp_.row_names_, "row", &r, &format!("{at}.row"))?;
            let j = find_index(
                &h.model_lp_.col_names_,
                "column",
                &cname,
                &format!("{at}.col"),
            )?;
            let v = edit_num(e, "value", at)?;
            ok(h.change_coeff(i, j, v), &kind)
        }
        "add_row" => {
            let name = edit_str(e, "name", at)?;
            if h.model_lp_.row_names_.contains(&name) {
                return Err(format!("{at}.name: a row named '{name}' already exists"));
            }
            let (idx, val) = coef_vectors(e, &h.model_lp_.col_names_, "column", at)?;
            let lo = bound_edit(e, "lower", -inf, -inf, at)?;
            let up = bound_edit(e, "upper", inf, inf, at)?;
            ok(h.add_row(lo, up, idx.len() as HighsInt, &idx, &val), &kind)?;
            let new = h.model_lp_.num_row_ - 1;
            ok(h.pass_row_name(new, &name), &kind)
        }
        "add_col" => {
            let name = edit_str(e, "name", at)?;
            if h.model_lp_.col_names_.contains(&name) {
                return Err(format!("{at}.name: a column named '{name}' already exists"));
            }
            let (idx, val) = coef_vectors(e, &h.model_lp_.row_names_, "row", at)?;
            let cost = match e.get("cost") {
                None => 0.0,
                Some(_) => edit_num(e, "cost", at)?,
            };
            let lo = bound_edit(e, "lower", 0.0, -inf, at)?;
            let up = bound_edit(e, "upper", inf, inf, at)?;
            let integer = opt_bool(e, "integer").map_err(|m| format!("{at}.{m}"))?;
            ok(
                h.add_col(cost, lo, up, idx.len() as HighsInt, &idx, &val),
                &kind,
            )?;
            let new = h.model_lp_.num_col_ - 1;
            ok(h.pass_col_name(new, &name), &kind)?;
            if integer == Some(true) {
                ok(h.change_col_integrality(new, HighsVarType::kInteger), &kind)?;
            }
            Ok(())
        }
        "del_row" => {
            let name = edit_str(e, "name", at)?;
            let i = find_index(&h.model_lp_.row_names_, "row", &name, &format!("{at}.name"))?;
            ok(h.delete_rows_set(1, &[i]), &kind)
        }
        "del_col" => {
            let name = edit_str(e, "name", at)?;
            let j = find_index(
                &h.model_lp_.col_names_,
                "column",
                &name,
                &format!("{at}.name"),
            )?;
            ok(h.delete_cols_set(1, &[j]), &kind)
        }
        _ => unreachable!("kind checked above"),
    }
}

/// Edits by name through the ported modify API (`highs_modify_api.rs`),
/// in order. An edit that fails stops the batch; the ones before it stay
/// applied and the message says how many.
fn modify(id: &str, m: &mut Model, args: &Args) -> ToolResult {
    only_fields(args, &["model_id", "edits"])?;
    let edits = match args.get("edits") {
        Some(Value::Array(a)) => a,
        Some(v) => {
            return Err(format!(
                "edits: expected an array of edit objects, got {}",
                type_name(v)
            ));
        }
        None => return Err("edits: required (an array of edit objects)".into()),
    };
    m.highs.set_option_value_bool("log_to_console", false);
    let mut applied = 0usize;
    for (k, e) in edits.iter().enumerate() {
        let at = format!("edits[{k}]");
        let Some(obj) = e.as_object() else {
            return Err(finish_err(
                m,
                applied,
                format!("{at}: expected an object, got {}", type_name(e)),
            ));
        };
        if let Err(msg) = apply_edit(&mut m.highs, obj, &at) {
            return Err(finish_err(m, applied, msg));
        }
        applied += 1;
        m.edits_applied += 1;
        m.solved = false;
    }
    Ok(modify_summary(id, m, applied))
}

fn finish_err(m: &Model, applied: usize, msg: String) -> String {
    if applied == 0 {
        format!("{msg}; no edit was applied")
    } else {
        format!(
            "{msg}; the {applied} edit(s) before it were applied (model now has {} edits)",
            m.edits_applied
        )
    }
}

fn modify_summary(id: &str, m: &Model, applied: usize) -> Value {
    let lp = &m.highs.model_lp_;
    json!({
        "model_id": id,
        "applied": applied,
        "edits_total": m.edits_applied,
        "rows": lp.num_row_,
        "cols": lp.num_col_,
        "integer_cols": count_integer(lp),
        "basis_kept": m.highs.basis_.valid,
    })
}

// ------------------------------------------------------ the subcommand cores

fn check_solution(id: &str, m: &mut Model, args: &Args) -> ToolResult {
    use crate::check;
    only_fields(
        args,
        &[
            "model_id",
            "solution",
            "primal_feasibility_tolerance",
            "mip_feasibility_tolerance",
            "max_findings",
        ],
    )?;
    let defaults = HighsOptions::default();
    let primal = opt_f64(args, "primal_feasibility_tolerance")?
        .unwrap_or(defaults.primal_feasibility_tolerance);
    let mip =
        opt_f64(args, "mip_feasibility_tolerance")?.unwrap_or(defaults.mip_feasibility_tolerance);
    for (f, v) in [
        ("primal_feasibility_tolerance", primal),
        ("mip_feasibility_tolerance", mip),
    ] {
        if !(v.is_finite() && v >= 0.0) {
            return Err(format!("{f}: expected a non-negative number, got {v}"));
        }
    }
    let max_findings = opt_usize(args, "max_findings")?.unwrap_or(check::DEFAULT_MAX_FINDINGS);
    let Some(sol) = args.get("solution") else {
        return Err(
            "solution: required ({\"values\": {column: value}} or {\"columns\": [{name, value}]})"
                .into(),
        );
    };
    let Some(sol) = sol.as_object() else {
        return Err(format!(
            "solution: expected an object, got {}",
            type_name(sol)
        ));
    };
    let mut pairs: Vec<(String, f64)> = Vec::new();
    if let Some(values) = sol.get("values") {
        let Some(map) = values.as_object() else {
            return Err("solution.values: expected an object {column name: number}".into());
        };
        for (n, v) in map {
            let x = v.as_f64().ok_or_else(|| {
                format!(
                    "solution.values.{n}: expected a number, got {}",
                    type_name(v)
                )
            })?;
            pairs.push((n.clone(), x));
        }
    } else if let Some(cols) = sol
        .get("columns")
        .or_else(|| sol.get("solution").and_then(|s| s.get("columns")))
    {
        let Some(cols) = cols.as_array() else {
            return Err("solution.columns: expected an array of {name, value}".into());
        };
        for (i, col) in cols.iter().enumerate() {
            let name = col
                .get("name")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("solution.columns[{i}].name: required string"))?;
            let x = col.get("value").and_then(Value::as_f64).ok_or_else(|| {
                format!("solution.columns[{i}].value: required number (column {name})")
            })?;
            pairs.push((name.to_string(), x));
        }
    } else {
        return Err(
            "solution: needs \"values\" ({column: value}) or \"columns\" ([{name, value}])".into(),
        );
    }
    let lp = &m.highs.model_lp_;
    let values = check::columns_by_name(lp, &pairs).map_err(|e| format!("solution: {}", e.0))?;
    let tol = check::Tolerances {
        primal_feasibility_tolerance: primal,
        mip_feasibility_tolerance: mip,
        applied: if lp.is_mip() { mip } else { primal },
    };
    let report = check::check_solution_objects(
        lp,
        None,
        m.cp.as_ref(),
        &values,
        tol,
        max_findings,
        &m.display(id),
        "<inline>",
    )
    .map_err(|e| format!("solution: {}", e.0))?;
    serde_json::to_value(report).map_err(|e| format!("internal error: {e}"))
}

fn explain(id: &str, m: &mut Model, args: &Args) -> ToolResult {
    only_fields(args, &["model_id", "verify", "max_members"])?;
    let verify = opt_bool(args, "verify")?.unwrap_or(true);
    let max_members = opt_usize(args, "max_members")?.unwrap_or(1000);
    let a = crate::explain::core_args(&m.display(id), verify, max_members);
    crate::explain::explain_core(m.highs.model_lp_.clone(), &a)
        .map(|(v, _)| v)
        .map_err(|e| e.msg)
}

fn lint(id: &str, m: &mut Model, args: &Args) -> ToolResult {
    only_fields(
        args,
        &[
            "model_id",
            "rules",
            "min_severity",
            "max_items",
            "with_solve",
        ],
    )?;
    let mut argv = vec![m.display(id), "--format".into(), "json".into()];
    for r in opt_str_list(args, "rules")? {
        argv.push("--rule".into());
        argv.push(r);
    }
    if let Some(s) = opt_str(args, "min_severity")? {
        argv.push("--min-severity".into());
        argv.push(s);
    }
    if let Some(n) = opt_usize(args, "max_items")? {
        argv.push("--max-items".into());
        argv.push(n.to_string());
    }
    if opt_bool(args, "with_solve")? == Some(true) {
        argv.push("--with-solve".into());
    }
    let a = crate::lint::parse_args(&argv).map_err(|e| {
        if e.starts_with("unknown rule") {
            format!("rules: {e}")
        } else if e.contains("--min-severity") {
            format!("min_severity: {e}")
        } else {
            e
        }
    })?;
    crate::lint::lint_core(m.highs.model_lp_.clone(), &a)
        .map(|(v, _)| v)
        .map_err(|e| e.msg)
}

fn ranging(id: &str, m: &mut Model, args: &Args) -> ToolResult {
    only_fields(args, &["model_id", "rows", "cols"])?;
    let rows = opt_str_list(args, "rows")?;
    let cols = opt_str_list(args, "cols")?;
    crate::ranging::ranging_core(&m.highs.model_lp_, &m.display(id), &rows, &cols)
        .map(|(v, _)| v)
        .map_err(|e| {
            e.msg
                .replace("--cols '", "cols: pattern '")
                .replace("--rows '", "rows: pattern '")
        })
}

fn feasrelax(id: &str, m: &mut Model, args: &Args) -> ToolResult {
    only_fields(
        args,
        &[
            "model_id",
            "mode",
            "global_lower_penalty",
            "global_upper_penalty",
            "global_rhs_penalty",
            "penalties",
            "big_m",
            "time_limit",
        ],
    )?;
    let mut flags: Vec<String> = Vec::new();
    if let Some(mode) = opt_str(args, "mode")? {
        if mode != "weighted" && mode != "min-count" {
            return Err(format!(
                "mode: expected \"weighted\" or \"min-count\", got \"{mode}\""
            ));
        }
        flags.extend(["--mode".into(), mode]);
    }
    for (f, flag) in [
        ("global_lower_penalty", "--global-lower-penalty"),
        ("global_upper_penalty", "--global-upper-penalty"),
        ("global_rhs_penalty", "--global-rhs-penalty"),
        ("big_m", "--big-m"),
        ("time_limit", "--time-limit"),
    ] {
        if let Some(x) = opt_f64(args, f)? {
            if f == "big_m" && !(x > 0.0 && x.is_finite()) {
                return Err(format!("{f}: must be positive and finite, got {x}"));
            }
            flags.extend([flag.into(), x.to_string()]);
        }
    }
    let penalties = opt_object(args, "penalties")?.map(|o| Value::Object(o.clone()));
    let a = crate::feasrelax::core_args(&m.display(id), &flags, penalties.is_some())
        .map_err(|e| format!("invalid arguments: {e}"))?;
    crate::feasrelax::feasrelax_core(&m.highs.model_lp_, &a, penalties.as_ref())
        .map(|(v, _)| v)
        .map_err(|e| {
            if e.code == 2 {
                format!("penalties: {}", e.msg)
            } else {
                e.msg
            }
        })
}

// ------------------------------------------------------------------- schemas

const TOOL_NAMES: [&str; 11] = [
    "load_model",
    "build_model",
    "solve",
    "get_solution",
    "modify",
    "check_solution",
    "explain",
    "lint",
    "ranging",
    "feasrelax",
    "free_model",
];

fn schema_of(text: &str) -> Value {
    serde_json::from_str(text).expect("committed schema is valid JSON")
}

fn model_id_prop() -> Value {
    json!({"type": "string", "description": "A handle returned by load_model."})
}

fn input(props: Value, required: &[&str]) -> Value {
    json!({
        "type": "object",
        "properties": props,
        "required": required,
        "additionalProperties": false,
    })
}

fn names_summary_schema() -> Value {
    json!({
        "type": "object",
        "required": ["count", "first"],
        "properties": {
            "count": {"type": "integer", "minimum": 0},
            "first": {"type": "array", "items": {"type": "string"}},
        },
    })
}

fn tool(
    name: &str,
    title: &str,
    description: &str,
    input: Value,
    output: Value,
    ro: bool,
) -> Value {
    json!({
        "name": name,
        "title": title,
        "description": description,
        "inputSchema": input,
        "outputSchema": output,
        "annotations": {"readOnlyHint": ro, "openWorldHint": false},
    })
}

/// The tool list, in the fixed order of [`TOOL_NAMES`].
fn tool_defs() -> Vec<Value> {
    let id = || json!({"model_id": model_id_prop()});
    let with_id = |extra: Value| {
        let mut p = id();
        p.as_object_mut()
            .unwrap()
            .extend(extra.as_object().unwrap().clone());
        p
    };
    let edit_schema = json!({
        "type": "object",
        "required": ["kind"],
        "additionalProperties": false,
        "description": "One edit. kind col_bounds/row_bounds: name, lower?, upper? (absent keeps the value, null is infinite). cost: name, cost. coef: row, col, value. add_row: name, lower?, upper?, coefs {column: value}. add_col: name, cost?, lower? (default 0), upper? (default infinite), integer?, coefs {row: value}. del_row / del_col: name.",
        "properties": {
            "kind": {"enum": ["col_bounds", "row_bounds", "cost", "coef", "add_row", "add_col", "del_row", "del_col"]},
            "name": {"type": "string", "description": "Row or column name."},
            "row": {"type": "string"},
            "col": {"type": "string"},
            "lower": {"type": ["number", "null"]},
            "upper": {"type": ["number", "null"]},
            "cost": {"type": "number"},
            "value": {"type": "number"},
            "integer": {"type": "boolean"},
            "coefs": {"type": "object", "additionalProperties": {"type": "number"}},
        },
    });
    vec![
        tool(
            "load_model",
            "Load a model",
            "Load a model and return a model_id handle with its size, name summary, class and engine. Every other tool takes the model_id. Give exactly one of path, mps_text, lp_text (CPLEX LP format) or mof_text (MathOptFormat JSON). LP and MIP models get the full tool set. Convex QP, SOCP/QCQP and general-constraint models (indicator, SOS, piecewise-linear, min/max/abs/and/or, semi-continuous) are routed to their engine and support solve, get_solution and free_model. Convex MIQP and mixed-integer conic models load too, but solve runs their experimental outer-approximation engine only with options {\"miqp_oa\": true} or {\"misocp_oa\": true} (off by default: slower), and otherwise refuses with that hint. A class nothing solves (nonconvex Q, nonconvex quadratic rows, general constraints with quadratic terms) is refused with the class and the reason.",
            input(
                json!({
                    "path": {"type": "string", "description": "Path of a model file: MPS, LP, AMPL .nl or MathOptFormat .mof.json, optionally .gz/.bz2/.xz/.zst."},
                    "mps_text": {"type": "string", "description": "The model as fixed- or free-format MPS text."},
                    "lp_text": {"type": "string", "description": "The model as CPLEX LP-format text (Minimize/Maximize, Subject To, Bounds, General/Binary, End). A quadratic objective is written in brackets: obj: x + [ x^2 + 2 x * y ] / 2."},
                    "mof_text": {"type": "string", "description": "The model as MathOptFormat JSON text (a .mof.json file's content)."},
                }),
                &[],
            ),
            json!({
                "type": "object",
                "required": ["model_id", "rows", "cols", "integer_cols", "names"],
                "properties": {
                    "model_id": {"type": "string"},
                    "name": {"type": "string"},
                    "rows": {"type": "integer", "minimum": 0},
                    "cols": {"type": "integer", "minimum": 0},
                    "integer_cols": {"type": "integer", "minimum": 0},
                    "nonzeros": {"type": "integer", "minimum": 0},
                    "sense": {"enum": ["minimize", "maximize"]},
                    "class": {"type": "string", "description": "The router's problem class: lp, mip, convex_qp, conic, gc, gc_mip, miqp, mi_conic."},
                    "engine": {"type": "string", "description": "The engine that solves the class: highs-port, hipo (QP), clarabel (conic), gc-linearization, oa (outer approximation, MIQP and mixed-integer conic, lever off by default)."},
                    "names": {
                        "type": "object",
                        "required": ["rows", "cols"],
                        "properties": {"rows": names_summary_schema(), "cols": names_summary_schema()},
                    },
                },
            }),
            false,
        ),
        tool(
            "build_model",
            "Build a model from JSON",
            "Build an LP/MIP from JSON, no file needed, and return a model_id like load_model. An optional cp section adds scheduling and combinatorial objects (class cp): intervals [{name, start, size, end, presence?}], no_overlap [{name, intervals}], cumulative [{name, intervals, demands, capacity}], all_different [{name, variables}], element [{name, index, values, target}], table [{name, variables, tuples, negated?}]. Every variable a cp object names (start, end, presence, index, target, variables) must be declared in variables, type integer (binary for presence), with both bounds finite: upper defaults to infinity, so write it. An interval does not create its end: declare the end variable and name it; start, size and end are all required (only presence is optional). size, demands, values and capacity take an integer (a JSON number without a fraction, so 3 not 3.0) or a variable name; a size or demand must be nonnegative. no_overlap and cumulative name intervals by their interval name, not by variable. The reply is the handle summary; its model_id is in structuredContent and in the JSON text of content[0].text. solve checks the answer on these objects (envelope cp.verified); modify, explain, lint, ranging and feasrelax refuse such a model (RS4012). variables: [{name, type? continuous|integer|binary, lower?, upper?, objective?}] (bounds default 0 and infinity, null is infinite). rows: [{name, terms: {variable: coefficient}, then either sense (<=, >=, =) with rhs, or lower and/or upper}]. objective: {terms?: {variable: coefficient}, offset?} adds to the per-variable objective. sense is minimize (default) or maximize. A bad entry is an error naming its JSON path, for example rows[2].terms.z. For a quadratic or general-constraint model write LP text and use load_model.",
            input(
                json!({
                    "name": {"type": "string", "description": "Model name (no whitespace)."},
                    "sense": {"enum": ["minimize", "maximize"]},
                    "variables": {
                        "type": "array",
                        "minItems": 1,
                        "items": {
                            "type": "object",
                            "required": ["name"],
                            "additionalProperties": false,
                            "properties": {
                                "name": {"type": "string"},
                                "type": {"enum": ["continuous", "integer", "binary"]},
                                "lower": {"type": ["number", "null"]},
                                "upper": {"type": ["number", "null"]},
                                "objective": {"type": "number", "description": "Objective coefficient."},
                            },
                        },
                    },
                    "rows": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "required": ["name", "terms"],
                            "additionalProperties": false,
                            "properties": {
                                "name": {"type": "string"},
                                "terms": {"type": "object", "additionalProperties": {"type": "number"}, "description": "{variable name: coefficient}."},
                                "sense": {"enum": ["<=", ">=", "="]},
                                "rhs": {"type": "number", "description": "With sense."},
                                "lower": {"type": ["number", "null"]},
                                "upper": {"type": ["number", "null"]},
                            },
                        },
                    },
                    "objective": {
                        "type": "object",
                        "additionalProperties": false,
                        "properties": {
                            "terms": {"type": "object", "additionalProperties": {"type": "number"}},
                            "offset": {"type": "number"},
                        },
                    },
                    "cp": {
                        "type": "object",
                        "description": "Scheduling and combinatorial objects. Every variable they name must be declared in variables as an integer (binary for presence) with finite lower and upper bounds. Variables and intervals are named by string.",
                        "additionalProperties": false,
                        "properties": {
                            "intervals": {"type": "array", "items": {"type": "object", "description": "{name, start, size, end, presence?}: start and end are required names of declared integer variables with finite bounds (the interval does not create them); size is an integer or a variable name (nonnegative); presence is an optional binary variable. start + size = end when present."}},
                            "no_overlap": {"type": "array", "items": {"type": "object", "description": "{name, intervals: [interval names]}"}},
                            "cumulative": {"type": "array", "items": {"type": "object", "description": "{name, intervals, demands (numbers or variables, one per interval), capacity (number or variable)}"}},
                            "all_different": {"type": "array", "items": {"type": "object", "description": "{name, variables}"}},
                            "element": {"type": "array", "items": {"type": "object", "description": "{name, index (variable, 0-based), values (numbers or variables), target (variable)}"}},
                            "table": {"type": "array", "items": {"type": "object", "description": "{name, variables, tuples: [[int,...]], negated?}"}},
                        },
                    },
                }),
                &["variables"],
            ),
            json!({
                "type": "object",
                "required": ["model_id", "rows", "cols", "integer_cols", "names"],
                "properties": {
                    "model_id": {"type": "string"},
                    "name": {"type": "string"},
                    "rows": {"type": "integer", "minimum": 0},
                    "cols": {"type": "integer", "minimum": 0},
                    "integer_cols": {"type": "integer", "minimum": 0},
                    "nonzeros": {"type": "integer", "minimum": 0},
                    "sense": {"enum": ["minimize", "maximize"]},
                    "class": {"type": "string"},
                    "engine": {"type": "string"},
                    "names": {
                        "type": "object",
                        "required": ["rows", "cols"],
                        "properties": {"rows": names_summary_schema(), "cols": names_summary_schema()},
                    },
                },
            }),
            false,
        ),
        tool(
            "solve",
            "Solve a model",
            "Solve the model handle (progress notifications with a progressToken, notifications/cancelled stops it) and return the rusolve.result/1 envelope (status, objective, dual bound, gap, stats.lp_iterations, run settings). After modify an LP re-solves from the kept basis (warm); set warm_start false to force a cold solve. Option values persist on the handle. On a convex MIQP or mixed-integer conic handle, options takes only {\"miqp_oa\": true} or {\"misocp_oa\": true} (experimental, off by default). With wasm_hook, a fuel-metered WebAssembly component scores branching and cuts and picks open nodes for this solve only; its counters appear as stats.wasm. options {\"mip_engine\": \"cpsat\"} solves an all-integer model, with or without CP objects, with OR-Tools CP-SAT when its library is installed (`rusolve version --format json` shows cpsat.available; a missing library is RS0012 naming the install routes); the answer is checked on your model and objects, and the envelope has a cpsat object.",
            input(
                with_id(json!({
                    // The option names and types come from the option registry
                    // (`commands::option_properties`); ranges and descriptions are `rusolve options --json`.
                    "options": {
                        "type": "object",
                        "description": "Flat {option name: value} map applied with the HiGHS option setters (for example {\"presolve\": \"off\", \"mip_rel_gap\": 0.01}). Names are the registry's; `rusolve options --json` lists type, range, default and description of each.",
                        "properties": crate::commands::option_properties(),
                        "additionalProperties": false,
                    },
                    "time_limit": {"type": "number", "minimum": 0, "description": "Seconds."},
                    "warm_start": {"type": "boolean", "description": "false clears the kept basis first (a cold solve). Default true."},
                    "include_solution": {"type": "boolean", "description": "Embed the incumbent by column name in the result. Default false; get_solution pages it instead."},
                    "wasm_hook": {"type": "string", "description": "Path of a .wasm/.wat hook component (docs/design/wasm-hooks.md): branch/cut scores and node choice from sandboxed WebAssembly, applied to this solve only. A component can only reorder, never change the answer; a bad component is a tool error naming the file."},
                    "wasm_fuel": {"type": "integer", "minimum": 0, "description": "Fuel (interpreter steps) per hook call with wasm_hook. Default 1000000; a call out of fuel falls back to the HiGHS rule."},
                    "wasm_memory": {"type": "integer", "minimum": 1, "description": "Linear-memory cap in bytes for the wasm_hook component. Default 1048576 (1 MiB); a component declaring more fails to load, which is a tool error naming the file."},
                })),
                &["model_id"],
            ),
            schema_of(SCHEMA_RESULT),
            false,
        ),
        tool(
            "get_solution",
            "Get the solution by name",
            "Column values of the last solve, by column name. Pass names for specific columns, or page through all of them with offset and limit (default 1000).",
            input(
                with_id(json!({
                    "names": {"type": "array", "items": {"type": "string"}},
                    "offset": {"type": "integer", "minimum": 0},
                    "limit": {"type": "integer", "minimum": 0},
                })),
                &["model_id"],
            ),
            json!({
                "type": "object",
                "required": ["model_id", "status", "total_columns", "returned", "truncated", "columns"],
                "properties": {
                    "model_id": {"type": "string"},
                    "status": {"type": "string"},
                    "objective": {"type": ["number", "null"]},
                    "total_columns": {"type": "integer", "minimum": 0},
                    "returned": {"type": "integer", "minimum": 0},
                    "offset": {"type": "integer", "minimum": 0},
                    "truncated": {"type": "boolean"},
                    "columns": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "required": ["name", "value"],
                            "properties": {"name": {"type": "string"}, "value": {"type": ["number", "null"]}},
                        },
                    },
                },
            }),
            true,
        ),
        tool(
            "modify",
            "Edit a model",
            "Apply edits by row/column name to the model handle: col_bounds, row_bounds, cost, coef, add_row, add_col, del_row, del_col. The basis is kept, so the next solve of an LP is warm. Edits apply in order; one that fails stops the batch (earlier edits stay applied).",
            input(
                with_id(json!({"edits": {"type": "array", "items": edit_schema}})),
                &["model_id", "edits"],
            ),
            json!({
                "type": "object",
                "required": ["model_id", "applied", "edits_total", "rows", "cols", "integer_cols", "basis_kept"],
                "properties": {
                    "model_id": {"type": "string"},
                    "applied": {"type": "integer", "minimum": 0},
                    "edits_total": {"type": "integer", "minimum": 0},
                    "rows": {"type": "integer", "minimum": 0},
                    "cols": {"type": "integer", "minimum": 0},
                    "integer_cols": {"type": "integer", "minimum": 0},
                    "basis_kept": {"type": "boolean"},
                },
            }),
            false,
        ),
        tool(
            "check_solution",
            "Check a solution",
            "Check a candidate solution against the model's bounds, rows and integrality, by name, and list the violations (rusolve.check/1). The solution is {\"values\": {column: value}} or {\"columns\": [{name, value}]} (a solve result with include_solution works).",
            input(
                with_id(json!({
                    "solution": {"type": "object"},
                    "primal_feasibility_tolerance": {"type": "number", "minimum": 0},
                    "mip_feasibility_tolerance": {"type": "number", "minimum": 0},
                    "max_findings": {"type": "integer", "minimum": 0},
                })),
                &["model_id", "solution"],
            ),
            schema_of(SCHEMA_CHECK),
            true,
        ),
        tool(
            "explain",
            "Explain infeasibility or unboundedness",
            "Solve the model's LP (the relaxation for a MIP) and, when it is infeasible, return an irreducible infeasible set or Farkas certificate by row/column name; when unbounded, a primal ray (rusolve.explain/1). The certificate is verified unless verify is false.",
            input(
                with_id(json!({
                    "verify": {"type": "boolean"},
                    "max_members": {"type": "integer", "minimum": 0},
                })),
                &["model_id"],
            ),
            schema_of(SCHEMA_EXPLAIN),
            true,
        ),
        tool(
            "lint",
            "Lint a model",
            "Static modelling checks (big-M, coefficient ranges, scaling, empty rows/columns, near-parallel rows, ...) with a fix hint per finding (rusolve.lint/1).",
            input(
                with_id(json!({
                    "rules": {"type": "array", "items": {"type": "string"}, "description": "Rule ids to run; default all."},
                    "min_severity": {"enum": ["info", "warning", "error"]},
                    "max_items": {"type": "integer", "minimum": 0},
                    "with_solve": {"type": "boolean", "description": "Also run the solve-based condition-estimate rule."},
                })),
                &["model_id"],
            ),
            schema_of(SCHEMA_LINT),
            true,
        ),
        tool(
            "ranging",
            "Sensitivity ranging",
            "Cost and bound ranging, reduced costs, duals and basis status by name for an LP (for a MIP, on the LP with integers fixed at the incumbent) (rusolve.ranging/1). rows and cols are name patterns ('*' wildcard); default all.",
            input(
                with_id(json!({
                    "rows": {"type": "array", "items": {"type": "string"}},
                    "cols": {"type": "array", "items": {"type": "string"}},
                })),
                &["model_id"],
            ),
            schema_of(SCHEMA_RANGING),
            true,
        ),
        tool(
            "feasrelax",
            "Cheapest feasibility relaxation",
            "For an infeasible model, the cheapest set of bound relaxations that makes it feasible, as named edits you can pass to modify (rusolve.feasrelax/1). mode weighted minimises the weighted violation; min-count minimises the number of relaxed constraints.",
            input(
                with_id(json!({
                    "mode": {"enum": ["weighted", "min-count"]},
                    "global_lower_penalty": {"type": "number"},
                    "global_upper_penalty": {"type": "number"},
                    "global_rhs_penalty": {"type": "number"},
                    "penalties": {"type": "object", "description": "{\"rows\": {name: weight}, \"cols\": {name: {\"lower\": w, \"upper\": w}}}; a negative weight forbids relaxing."},
                    "big_m": {"type": "number", "exclusiveMinimum": 0},
                    "time_limit": {"type": "number", "minimum": 0},
                })),
                &["model_id"],
            ),
            schema_of(SCHEMA_FEASRELAX),
            true,
        ),
        tool(
            "free_model",
            "Free a model handle",
            "Release a model handle and its solver state.",
            input(id(), &["model_id"]),
            json!({
                "type": "object",
                "required": ["model_id", "freed", "remaining"],
                "properties": {
                    "model_id": {"type": "string"},
                    "freed": {"type": "boolean"},
                    "remaining": {"type": "integer", "minimum": 0},
                },
            }),
            false,
        ),
    ]
}
