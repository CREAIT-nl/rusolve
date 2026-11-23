//! The CLI front of the rusolve-ekk MIP port.
//!
//! * `rusolve mip <model.(mps|lp|nl|mof.json)[.gz|.bz2|.xz|.zst]> [--csv | --format json] [--solution | --include-solution]
//!   [--solution-file PATH] [--time-limit S] [--node-limit N] [--threads N] [--no-presolve]
//!   [--seed N] [--mip-rel-gap X] [--mip-abs-gap X] [--profile FILE] [--options-file FILE]
//!   [--opt key=value]... [--wasm-hook PATH [--wasm-fuel N] [--wasm-memory BYTES]]`
//!   (the flag table is `commands::MIP_FLAGS`; `rusolve mip --help` prints it. `--profile` is a
//!   `rusolve.profile/1` file, see `profile.rs`; order of precedence: profile, options
//!   file, `--opt`)
//!   by default (`RUSOLVE_MIP_HIGHS=0` selects the native search): the MIP
//!   solved by the ported `Highs::run` -> `callSolveMip`, printing the same
//!   `status,objective,bound,nodes,lp_iters,time_s,x` CSV as the native
//!   `rusolve-mip` path, so a script can A/B the two.
//! * `rusolve mip-trace <model.mps> [time_limit_s]`: the Rust half of the
//!   acceptance harness; prints the lines of tools/highs-oracle/mip_trace.cpp
//!   (`log`, `cb`, `sol`, `final`), which must match it byte for byte.
//!
//! ## Output modes
//!
//! `--csv` (and `--format csv`) is the harness contract and is byte-identical to
//! before the envelope existed. `--format json` replaces it with one
//! `rusolve.result/1` JSON object on stdout
//! ([`envelope`]) and moves the HiGHS log to stderr; `--solution` adds the
//! incumbent by name to that object, and `--solution-file` writes the
//! HiGHS raw solution file, or the envelope's `solution` object when the path
//! ends in `.json`. `--csv` with `--format json` is a usage error.
//!
//! The model is the `HighsLp` HMpsFF would build (`ekk::highs_lp_from_model`,
//! plus the names the reader keeps, and the integrality vector: rusolve-io already
//! applies HMpsFF's marker-binary and `BV`/`LI`/`UI` rules). The driver itself is
//! ours; the options it sets are the ones mip_trace.cpp sets.
//!
//! ## Exit codes
//!
//! | code | meaning |
//! |---|---|
//! | 0 | the solve ran, whatever the model status: `Optimal`, `Infeasible`, `TimeLimit`, `NodeLimit` and the rest are reported in the output, not in the code |
//! | 2 | usage error: an unsupported flag, a flag with no value, a value that will not parse, a bad `--opt`, `--csv` together with `--format json`, `--wasm-fuel`/`--wasm-memory` without `--wasm-hook` |
//! | 3 | input error: the model file is missing, unreadable, or does not parse; the `--wasm-hook` component file is missing, unreadable, or is not a loadable module |
//! | 4 | internal error: `Highs::pass_model` or `Highs::run` returned `HighsStatus::kError`, or a requested `--solution-file` could not be written |
//!
//! The table applies in every output mode, because nothing that consumes this
//! path decides on the code: `tools/bench/abtest.py` records `rc` but builds its row
//! from the CSV fields, `tools/bench/identity.sh` reads `tail -1` of stdout
//! through a pipeline and never looks at the code, and `tools/diff/mip_diff.py`
//! records `exit` while parsing the CSV regardless. The native `rusolve mip` path keeps
//! the same codes (`run_mip` in `main.rs`); a flag only the other stack has is a usage
//! error that names the stack (`commands::flag_stack_error`).
//!
//! ## Who parses the flags
//!
//! [`commands::parse_flags`] does: it walks the arguments, looks every flag up in
//! `commands::MIP_FLAGS`, and refuses an unknown flag, a flag of the other stack, a missing
//! value, a value that is not of the flag's kind or not in its list, and a repeat of a flag
//! that is not repeatable - all exit 2, all naming the flag and the value. What a flag *does*
//! is still here, in the same order as the flags appear, and so are the checks between flags
//! (`--csv` with `--format json`, `--wasm-fuel` without `--wasm-hook`, `--watchdog` without a
//! time limit) and the range checks of the port's own types (`--threads >= 1`, a `--seed`
//! that does not fit a `HighsInt`). A flag added to the registry is therefore parsed,
//! documented and typed here without a line of parser code.

use rusolve_ekk::lp_data::highs_option_records as registry;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use std::process::ExitCode;
use std::sync::{Arc, Mutex};

/// `true` unless `RUSOLVE_MIP_HIGHS=0`: `rusolve mip` runs the port by default.
pub fn enabled() -> bool {
    std::env::var("RUSOLVE_MIP_HIGHS").map_or(true, |v| v != "0")
}

/// `HighsLp` with integrality and the file's names, from any format `rusolve_io` reads
/// (MPS, LP, AMPL `.nl`, MathOptFormat JSON, optionally gzip/bzip2/xz/zstd compressed; the
/// same formats on both stacks): an MPS file in its HMpsFF file order, the others in CSC order
/// (`rusolve_io::read_model_file_order`). Shared by every subcommand that loads a model for the port (check, iis,
/// ranging, lint, ...), so they all see the same names and column order.
pub(crate) fn read_mip(path: &str) -> Result<HighsLp, String> {
    read_mip_q(path).map(|(lp, _)| lp)
}

/// [`read_mip`] plus whether the model carries a quadratic objective (a Hessian), which the
/// MIP stack ignores: `rusolve mip` refuses it.
pub(crate) fn read_mip_q(path: &str) -> Result<(HighsLp, bool), String> {
    let (lp, has_hessian, gc) = read_mip_full(path)?;
    // A model with general constraints is only understood by the commands
    // that take a `GcModel` (solve, check). Every other one names the object.
    if let Some(what) = gc.first_general() {
        return Err(format!(
            "read {path}: {what} (general constraints are supported by `rusolve mip` and `rusolve check` only)"
        ));
    }
    Ok((lp, has_hessian))
}

/// Everything a model file carries for the port: the LP, whether it has a Hessian and
/// its general-constraint objects.
pub(crate) fn read_mip_full(path: &str) -> Result<(HighsLp, bool, rusolve_io::GcModel), String> {
    let (model, order) =
        rusolve_io::read_model_file_order(path).map_err(|e| format!("read {path}: {e}"))?;
    let has_hessian = !model.hessian.is_empty();
    let mut lp = crate::ekk::highs_lp_from_model(&model, &order);
    // The names the reader kept
    // survive into the port, as HiGHS's MPS reader sets them. `HMpsFF::parse`
    // moves the two vectors it built into the LP (1.16.0 io/HMpsFF.cpp:106-107,
    // `lp.row_names_ = std::move(row_names); lp.col_names_ = std::move(col_names)`).
    // What it stores:
    // * `row_names_` is the constraint rows in ROWS order. The objective row is
    //   not a row: `N` rows push nothing onto `row_lower`/`row_upper`
    //   (:632-651) and the first one names the objective instead of a row
    //   (`objective_name`, :683-687); further `N` rows are free rows, mapped to
    //   -2 and skipped (:671-675). rusolve-io does the same (`mps.rs:349-369`
    //   pushes a name only for a non-`N` row, in `con_index` order).
    // * `col_names_` is every column in COLUMNS order (`HMpsFF.cpp:499`,
    //   :863; `mps.rs:126-127`).
    // * A duplicated name makes HiGHS warn and clear both vectors
    //   (:64-78), so it never hands the solver a partial name list; rusolve-io
    //   rejects a duplicate row name outright (`mps.rs:345`) and cannot name two
    //   columns alike (`get_col_idx` looks the name up first), so these vectors
    //   are always the complete ones HiGHS would keep.
    crate::ekk::fill_names(&mut lp, &model);
    // Highs::readModel (Highs.cpp:814) names the model after the file
    lp.model_name_ = extract_model_name(path);
    if model.bounds.iter().any(|b| b.integer) || !model.gc.semi.is_empty() {
        lp.integrality_ = model
            .bounds
            .iter()
            .map(|b| {
                if b.integer {
                    HighsVarType::kInteger
                } else {
                    HighsVarType::kContinuous
                }
            })
            .collect();
        // Semi-continuous columns: the port's `without_semi_variables` handles these.
        for s in &model.gc.semi {
            lp.integrality_[s.col] = if s.integer {
                HighsVarType::kSemiInteger
            } else {
                HighsVarType::kSemiContinuous
            };
        }
    }
    Ok((lp, has_hessian, model.gc))
}

/// HiGHS: io/Filereader.cpp::extractModelName (lines 75-87): the file name without
/// its directory and its extension (and a trailing `.gz`). rusolve also reads bzip2, xz and
/// zstd files, which HiGHS does not, so their extensions are stripped the same way.
fn extract_model_name(filename: &str) -> String {
    let mut name = filename.to_string();
    if let Some(found) = name.rfind(['/', '\\']) {
        name = name[found + 1..].to_string();
    }
    let mut found = name.rfind('.');
    // C++ `name.substr(found + 1)` with `found == npos` is the whole name
    let ext = match found {
        Some(f) => &name[f + 1..],
        None => name.as_str(),
    };
    if matches!(ext, "gz" | "bz2" | "xz" | "zst" | "zstd" | "gzip") {
        let f = found.expect("a compression extension has a dot or is the whole name");
        name.truncate(f);
        found = name.rfind('.');
    }
    if let Some(f) = found {
        name.truncate(f);
    }
    name
}

/// The options mip_trace.cpp sets: serial, seed 0, logging to the callback
/// only.
fn trace_options(time_limit: f64, random_seed: HighsInt) -> HighsOptions {
    HighsOptions {
        output_flag: false,
        threads: 1,
        random_seed,
        time_limit,
        ..crate::base_options()
    }
}

/// The CSV status spelling of `rusolve mip` (`MipStatus` `Debug` names),
/// else HiGHS `modelStatusToString`.
fn csv_status(status: HighsModelStatus) -> String {
    match status {
        HighsModelStatus::kOptimal => "Optimal".to_string(),
        HighsModelStatus::kInfeasible => "Infeasible".to_string(),
        HighsModelStatus::kUnbounded => "Unbounded".to_string(),
        HighsModelStatus::kTimeLimit => "TimeLimit".to_string(),
        HighsModelStatus::kSolutionLimit => "NodeLimit".to_string(),
        other => util_model_status_to_string(other).to_string(),
    }
}

/// `rusolve mip` routed to the port: see the module doc for the flags, the
/// output modes and the exit-code table.
pub fn run_mip_highs(args: &[String]) -> ExitCode {
    let Some(path) = args.first() else {
        return crate::diag::fail_msg(
            crate::diag::wants_json(args),
            "RS0006",
            2,
            "rusolve mip: missing model path (see `rusolve mip --help`)",
        );
    };
    // The command registry parses the command line (`commands::parse_flags`): it looks
    // every flag up in `commands::MIP_FLAGS`, refuses an unknown flag, a flag of the other
    // stack, a missing value, a value that is not of the flag's kind or not in its list, and
    // a repeat of a flag that is not repeatable (each an exit-2 usage error naming the flag).
    // What each flag *does* stays below, and so do the checks between flags.
    let parsed = match crate::commands::parse_flags("mip", crate::commands::Stack::Port, &args[1..])
    {
        Ok(parsed) => parsed,
        Err(err) => return crate::diag::fail_usage(crate::diag::wants_json(args), &err),
    };
    if let Some(extra) = parsed.positionals().first() {
        return crate::diag::fail_msg(
            crate::diag::wants_json(args),
            "RS0006",
            2,
            format!("rusolve mip: unexpected argument '{extra}' (see `rusolve mip --help`)"),
        );
    }
    // `--format json` puts one envelope object on stdout and the log on
    // stderr instead of the CSV line
    let format = parsed.get_str("--format");
    let csv = parsed.has("--csv") || format == Some("csv");
    let json = format == Some("json");
    // Include the named incumbent in the envelope (both spellings mean the
    // same; the envelope is small unless one of them is given)
    let include_solution = parsed.has("--solution") || parsed.has("--include-solution");
    // Write the solution file (HiGHS raw format, or the envelope's
    // `solution` object when the path ends in `.json`)
    let solution_file = parsed.get_str("--solution-file").map(str::to_owned);
    // Warm-start read side (warm.rs): MIP start, LP basis, rows added before the solve
    let start_file = parsed.get_str("--start").map(str::to_owned);
    let basis_file = parsed.get_str("--basis").map(str::to_owned);
    let add_rows_file = parsed.get_str("--add-rows").map(str::to_owned);
    // wasm-hooks.md: Wasm hook component and its budgets (`--wasm-hook PATH`,
    // `--wasm-fuel N`, `--wasm-memory BYTES`; the budgets need the hook)
    let wasm_hook = parsed.get_str("--wasm-hook").map(str::to_owned);
    let wasm_fuel = parsed.get_u64("--wasm-fuel");
    let wasm_memory = parsed.get_u64("--wasm-memory").map(|m| m as usize);
    // `--log-format jsonl` streams progress events as NDJSON on stderr
    let jsonl = parsed.get_str("--log-format") == Some("jsonl");
    // `--watchdog` (needs a time limit) and `--watchdog-grace S`
    let watchdog = parsed.has("--watchdog");
    let mut watchdog_grace = 5.0_f64;
    if let Some(g) = parsed.get_f64("--watchdog-grace") {
        if g.is_nan() || g < 0.0 {
            return crate::diag::fail_msg(
                json,
                "RS0003",
                2,
                "--watchdog-grace needs a non-negative number of seconds",
            );
        }
        watchdog_grace = g;
    }
    let mut time_limit = kHighsInf;
    if let Some(t) = parsed.get_f64("--time-limit") {
        time_limit = t;
    }
    // The port's limits are `HighsInt`: a value too large for it is the same usage error a
    // hand-parsed `--node-limit`/`--seed` gave.
    let mut node_limit: Option<HighsInt> = None;
    if let Some(n) = parsed.get_i64("--node-limit") {
        match HighsInt::try_from(n) {
            Ok(n) => node_limit = Some(n),
            Err(_) => {
                return crate::diag::fail_msg(json, "RS0003", 2, "--node-limit needs an integer");
            }
        }
    }
    let mut seed: HighsInt = 0;
    if let Some(s) = parsed.get_i64("--seed") {
        match HighsInt::try_from(s) {
            Ok(s) => seed = s,
            Err(_) => {
                return crate::diag::fail_msg(json, "RS0003", 2, "--seed needs an integer");
            }
        }
    }
    // `--threads` is accepted on both stacks. Without the flag the port's MIP search
    // is serial (`trace_options` pins `threads: 1`). With it the value reaches the
    // `threads` option, which sizes the task scheduler; under the default profile
    // (`h16_mip_workers`, `parallel=choose`) the tree search then runs N workers.
    let mut threads: Option<HighsInt> = None;
    if let Some(n) = parsed.get_i64("--threads") {
        match HighsInt::try_from(n) {
            Ok(n) if n >= 1 => threads = Some(n),
            _ => {
                return crate::diag::fail_msg(json, "RS0003", 2, "--threads needs an integer >= 1");
            }
        }
    }
    // `key=value` option settings in argument order; --mip-rel-gap and --mip-abs-gap
    // are shorthands for two of them and --no-presolve for a third, so every value is
    // checked the same way. `--opt` is a `KeyValue` flag, so the `=` and the key are
    // already checked by the registry.
    let settings: Vec<(String, String)> = parsed
        .flags()
        .iter()
        .filter_map(|given| match given.name {
            "--opt" => given
                .value
                .as_deref()
                .and_then(|v| v.split_once('='))
                .map(|(k, v)| (k.to_string(), v.to_string())),
            "--mip-rel-gap" | "--mip-abs-gap" => given
                .value
                .clone()
                .map(|v| (given.name.trim_start_matches("--").replace('-', "_"), v)),
            "--no-presolve" => Some(("presolve".to_string(), "off".to_string())),
            _ => None,
        })
        .collect();
    // --options-file: HiGHS' options file, applied after the profile and before --opt,
    // so a command line still wins over the file.
    let options_file = parsed.get_str("--options-file").map(str::to_owned);
    // --profile is a rusolve.profile/1 file (profile.rs); it replaces the starting options
    // and is applied before --options-file and --opt, so the command line still wins.
    let profile_file = parsed.get_str("--profile").map(str::to_owned);
    if csv && json {
        // `--csv` is the byte-for-byte contract of tools/bench/abtest.py and
        // tools/bench/identity.sh; refusing the combination is honest where
        // silently dropping one of them would not be.
        return crate::diag::fail_msg(
            json,
            "RS0003",
            2,
            "--csv and --format json cannot be combined: --format json replaces the CSV line",
        );
    }
    if wasm_hook.is_none() && (wasm_fuel.is_some() || wasm_memory.is_some()) {
        return crate::diag::fail_msg(
            json,
            "RS0003",
            2,
            "--wasm-fuel and --wasm-memory need --wasm-hook PATH (they configure the component)",
        );
    }
    let profile = match &profile_file {
        Some(p) => match crate::profile::Profile::read(p) {
            Ok(profile) => Some(profile),
            Err((code, msg)) => {
                return crate::diag::fail_coded(json, "mip", code, &msg);
            }
        },
        None => None,
    };
    let (mut lp, has_hessian, gc) = match read_mip_full(path) {
        Ok(p) => p,
        Err(e) => {
            return crate::diag::fail_read(
                crate::diag::wants_json(args),
                path,
                e.strip_prefix(&format!("read {path}: ")).unwrap_or(&e),
            );
        }
    };
    if has_hessian {
        return crate::diag::fail_msg(
            json,
            "RS0012",
            3,
            format!("rusolve mip: {}", crate::qp::unsupported("mip")),
        );
    }
    // A model with general constraints is linearized here; one without
    // keeps the port's path untouched.
    let mut gc_plan: Option<crate::gc_solve::GcPlan> = None;
    if gc.has_general() {
        if start_file.is_some() || basis_file.is_some() || add_rows_file.is_some() {
            return crate::diag::fail_msg(
                json,
                "RS0003",
                2,
                format!(
                    "rusolve mip: --start, --basis and --add-rows are not supported with general constraints ({})",
                    gc.first_general().unwrap_or_default()
                ),
            );
        }
        let bigm = parsed.get_f64("--gc-bigm-max").unwrap_or(1e6);
        match crate::gc_solve::prepare(&mut lp, gc, bigm) {
            Ok(plan) => gc_plan = Some(plan),
            Err(e) => {
                return crate::diag::fail_msg(
                    json,
                    "RS0012",
                    3,
                    format!("rusolve mip: cannot linearize {e}"),
                );
            }
        }
    }
    // The cheap lint rules, on a copy, before the solve (never touches the search)
    let warn = !parsed.has("--no-warnings");
    let mut diags = if warn {
        crate::lint::warnings(&lp)
    } else {
        Vec::new()
    };
    crate::adopt::print_warnings(&diags);
    let mut options = match &profile {
        // `Profile::read` already ran `effective()` once, so it cannot fail here
        Some(p) => HighsOptions {
            output_flag: false,
            threads: 1,
            random_seed: seed,
            time_limit,
            ..p.effective().expect("a validated profile applies")
        },
        None => trace_options(time_limit, seed),
    };
    if let Some(p) = &profile
        && let Some(msg) = crate::profile::region_message(p, path)
    {
        eprintln!("{msg}");
    }
    if let Some(n) = node_limit {
        options.mip_max_nodes = n;
    }
    if let Some(n) = threads {
        options.threads = n;
    }
    options.sync_log_options();
    let start = std::time::Instant::now();
    let mut highs = Highs::with_options(options);
    if let Some(path) = &options_file {
        let entries = match registry::read_options_file(path) {
            Ok(entries) => entries,
            Err(err) => {
                return crate::diag::fail_msg(json, "RS0010", 2, format!("rusolve mip: {err}"));
            }
        };
        for entry in &entries {
            if registry::option_type_name(&entry.name).is_none() {
                return crate::diag::unknown_option(json, "mip", &entry.name);
            }
            if highs.set_option_value_str(&entry.name, &entry.value) == HighsStatus::kError {
                return crate::diag::fail_msg(
                    json,
                    "RS0004",
                    2,
                    format!(
                        "rusolve mip: unknown option or invalid value: {}={} ({} line {})",
                        entry.name, entry.value, path, entry.line
                    ),
                );
            }
        }
    }
    for (key, value) in &settings {
        if registry::option_type_name(key).is_none() {
            return crate::diag::unknown_option(json, "mip", key);
        }
        if highs.set_option_value_str(key, value) == HighsStatus::kError {
            return crate::diag::fail_msg(
                json,
                "RS0004",
                2,
                format!("rusolve mip: unknown option or invalid value: {key}={value}"),
            );
        }
    }
    // Wasm hooks (docs/design/wasm-hooks.md, rusolve's own extension surface — not a HiGHS port):
    // a sandboxed Wasm hook component for branching scores, cut scores and node selection.
    // A bad component file (unreadable, invalid module, an import, a bad export signature) is
    // an input error naming the file (exit code 3). Without --wasm-hook nothing is loaded and
    // no hooks are installed, so the search is the unmodified port one (identity 20/20).
    let mut wasm_report: Option<(String, rusolve_wasm::WasmStatsHandle)> = None;
    if let Some(p) = &wasm_hook {
        let name = std::path::Path::new(p)
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| p.clone());
        let bytes = match std::fs::read(p) {
            Ok(bytes) => bytes,
            Err(e) => {
                return crate::diag::fail_msg(
                    json,
                    "RS0010",
                    3,
                    format!("rusolve mip: --wasm-hook: cannot read {p}: {e}"),
                );
            }
        };
        let mut limits = rusolve_wasm::WasmLimits::default();
        if let Some(f) = wasm_fuel {
            limits.fuel_per_call = f;
        }
        if let Some(m) = wasm_memory {
            limits.max_memory_bytes = m;
        }
        match rusolve_wasm::WasmHooks::new(&name, &bytes, limits) {
            Ok(hooks) => {
                wasm_report = Some((name, hooks.stats()));
                highs.set_mip_hooks(Box::new(hooks));
            }
            Err(e) => {
                return crate::diag::fail_msg(
                    json,
                    "RS0011",
                    3,
                    format!("rusolve mip: --wasm-hook: {p}: {e}"),
                );
            }
        }
    }
    // With `--format json`, stdout carries exactly one JSON object. The port
    // writes its log with `print!` (rusolve-ekk/src/io/highs_io.rs:150-191), so
    // the log is moved off stdout by switching `log_to_console` off and
    // re-emitting the same text from the `kCallbackLogging` callback - the
    // mechanism `mip-trace` uses below, and no change to the port's IO. Setting
    // it after the `--opt` loop means `--opt log_to_console=true` cannot put the
    // log back on stdout; it is installed only when logging is on at all, so a
    // plain `--format json` run registers no callback and the search never sees
    // one. The callback only ever writes: it never sets `user_interrupt` or
    // `user_has_solution`, and no other callback kind is started, so
    // `interrupt_from_callback_with_data`
    // (rusolve-ekk/src/mip/mip_solver_data/part3.rs:1203-1215) returns early at
    // every other fire site.
    // The progress tracker is the only callback when it is asked for (the text log is
    // not written then: stderr carries events only); `rusolve_cli::progress` documents it.
    let tracker = (jsonl || watchdog).then(|| crate::progress::Tracker::new(start, jsonl));
    let watchdog_limit = highs.get_options().time_limit;
    if watchdog && !watchdog_limit.is_finite() {
        return crate::diag::fail_msg(
            json,
            "RS0003",
            2,
            "--watchdog needs a time limit (--time-limit S)",
        );
    }
    if let Some(tracker) = &tracker {
        highs.set_option_value_bool("log_to_console", false);
        highs.set_callback(tracker.callback());
        highs.start_callback(kCallbackMipImprovingSolution);
        highs.start_callback(kCallbackMipInterrupt);
        tracker.phase("presolve");
    } else if json && highs.get_options().output_flag {
        highs.set_option_value_bool("log_to_console", false);
        highs.set_callback(HighsCallbackFunctionType::new(
            move |kind: HighsInt,
                  message: &str,
                  _out: &HighsCallbackOutput,
                  _cb_in: &mut HighsCallbackInput| {
                if kind == kCallbackLogging {
                    eprint!("{message}");
                }
            },
        ));
        highs.start_callback(kCallbackLogging);
    }
    let pass_status = highs.pass_model(lp);
    // Warm-start read side (warm.rs), after the model is in and before the solve: rows first (a
    // start and a basis are matched against the model including them), then the start, then
    // the basis. Without the three flags nothing below runs.
    let mut add_rows_report = None;
    let mut start_report = None;
    let mut basis_report = None;
    let is_mip_model = highs.model_lp_.is_mip();
    if pass_status != HighsStatus::kError {
        if let Some(p) = &add_rows_file {
            let dbg_path = highs.get_options().mip_debug_solution_file.clone();
            let feastol = highs.get_options().mip_feasibility_tolerance;
            let gated = (|| {
                let dbg = if dbg_path.is_empty() {
                    None
                } else {
                    Some(crate::warm::read_debug_solution(
                        &highs.model_lp_,
                        &dbg_path,
                    )?)
                };
                crate::warm::read_and_gate_rows(&highs.model_lp_, p, dbg.as_deref(), feastol)
            })();
            match gated {
                Ok(rep) => {
                    if let Err(e) = crate::warm::add_rows(&mut highs, &rep.accepted) {
                        return crate::diag::fail_msg(
                            json,
                            "RS0090",
                            4,
                            format!("rusolve mip: {p}: {e}"),
                        );
                    }
                    for (name, why) in &rep.rejected {
                        eprintln!("rusolve mip: --add-rows: row {name} rejected: {why}");
                    }
                    add_rows_report = Some(rep);
                }
                Err(e) => {
                    return crate::diag::fail_msg(
                        json,
                        "RS0011",
                        3,
                        format!("rusolve mip: --add-rows: {e}"),
                    );
                }
            }
        }
        if let Some(p) = &start_file {
            let o = highs.get_options();
            let tol = crate::check::Tolerances {
                primal_feasibility_tolerance: o.primal_feasibility_tolerance,
                mip_feasibility_tolerance: o.mip_feasibility_tolerance,
                applied: if is_mip_model {
                    o.mip_feasibility_tolerance
                } else {
                    o.primal_feasibility_tolerance
                },
            };
            match crate::warm::read_start(&highs.model_lp_, p, tol) {
                Ok((rep, given)) => {
                    // the ported setSolution (Highs.cpp:2492 full, :2576 partial)
                    let st = if rep.complete {
                        let sol = HighsSolution {
                            col_value: rep.values.clone().unwrap_or_default(),
                            ..Default::default()
                        };
                        highs.set_solution(&sol)
                    } else {
                        let idx: Vec<HighsInt> = given.iter().map(|g| g.0 as HighsInt).collect();
                        let val: Vec<f64> = given.iter().map(|g| g.1).collect();
                        highs.set_solution_sparse(idx.len() as HighsInt, &idx, &val)
                    };
                    if st == HighsStatus::kError {
                        return crate::diag::fail_msg(
                            json,
                            "RS0011",
                            3,
                            format!("rusolve mip: --start: the ported setSolution rejected {p}"),
                        );
                    }
                    start_report = Some(rep);
                }
                Err(e) => {
                    return crate::diag::fail_msg(
                        json,
                        "RS0011",
                        3,
                        format!("rusolve mip: --start: {e}"),
                    );
                }
            }
        }
        if let Some(p) = &basis_file {
            let read = std::fs::read_to_string(p)
                .map_err(|e| format!("cannot read the basis file {p}: {e}"))
                .and_then(|t| crate::warm::read_basis_text(&highs.model_lp_, &t));
            match read {
                Ok((basis, v1)) => {
                    let installed = !is_mip_model
                        && highs.set_basis(&basis, "readBasis") != HighsStatus::kError;
                    if !is_mip_model && !installed {
                        return crate::diag::fail_msg(
                            json,
                            "RS0011",
                            3,
                            format!(
                                "rusolve mip: --basis: {p}: the basis is not consistent with the model"
                            ),
                        );
                    }
                    basis_report = Some(crate::warm::basis_json(p, v1, installed, is_mip_model));
                }
                Err(e) => {
                    return crate::diag::fail_msg(
                        json,
                        "RS0011",
                        3,
                        format!("rusolve mip: --basis: {e}"),
                    );
                }
            }
        }
    }
    if watchdog && let Some(tracker) = &tracker {
        // The envelope of this model without a solve result, completed by the watchdog from
        // the last state the callbacks saw if the solve does not return in time.
        let template = crate::envelope::envelope(&crate::envelope::Solved {
            path,
            lp: &highs.model_lp_,
            status: HighsModelStatus::kNotset,
            run_status: HighsStatus::kOk,
            info: &HighsInfo::default(),
            solution: &HighsSolution::default(),
            time_s: 0.0,
            seed,
            threads: highs.get_options().threads,
            time_limit,
            node_limit,
            mip_max_nodes: highs.get_options().mip_max_nodes,
            options: &settings,
            include_solution: false,
            mip_stats: None,
            // built before the run: the hook counters are not known yet
            wasm: None,
        });
        let (jsonl, json, csv) = (jsonl, json, csv);
        tracker.spawn_watchdog(
            watchdog_limit,
            watchdog_grace,
            Box::new(move |tracker| {
                let env = tracker.hard_stop_envelope(template);
                if jsonl {
                    tracker.emit_final(&env);
                }
                let objective = env.objective.unwrap_or(f64::NAN);
                let bound = env.dual_bound.unwrap_or(f64::NEG_INFINITY);
                if json {
                    println!(
                        "{}",
                        serde_json::to_string(&env)
                            .unwrap_or_else(|e| format!("{{\"error\":{e:?}}}"))
                    );
                } else if csv {
                    println!("status,objective,bound,nodes,lp_iters,time_s,x");
                    println!(
                        "TimeLimit,{objective},{bound},{},{},{},[]",
                        env.stats.nodes, env.stats.lp_iterations, env.stats.time_s
                    );
                } else {
                    println!(
                        "status:    {}",
                        util_model_status_to_string(HighsModelStatus::kTimeLimit)
                    );
                    println!("objective: {objective}");
                    println!("bound:     {bound}");
                    println!(
                        "nodes:     {} ({} LP iterations, {:.3}s)",
                        env.stats.nodes, env.stats.lp_iterations, env.stats.time_s
                    );
                }
            }),
        );
    }
    let run_status = if pass_status == HighsStatus::kError {
        HighsStatus::kError
    } else {
        // `--threads N` starts an N-worker scheduler for this run
        highs.reset_scheduler_for_threads();
        highs.run()
    };
    // Check the user's objects, repair trickle flow once, cut the solution
    let gc_outcome = gc_plan.as_ref().map(|plan| {
        let o = highs.get_options();
        let tol = rusolve_io::GcTol {
            integrality: o.mip_feasibility_tolerance,
            feasibility: o.mip_feasibility_tolerance,
        };
        crate::gc_solve::finish(plan, &mut highs, tol)
    });
    if let Some(o) = &gc_outcome {
        if let Some(t) = o.report.get("trickle_flow") {
            eprintln!("rusolve mip: trickle flow repaired once: {t}");
        }
        if let Some(v) = o.report["violations"].as_array().filter(|v| !v.is_empty()) {
            eprintln!(
                "rusolve mip: WARNING: the incumbent violates {} general constraint(s) of the model (first: {})",
                v.len(),
                v[0]
            );
        }
    }
    let secs = start.elapsed().as_secs_f64();
    let status = highs.get_model_status();
    let info = highs.get_info();
    let has_solution = highs.solution_.value_valid;
    let objective = if has_solution {
        info.objective_function_value
    } else {
        f64::NAN
    };
    // Claim the final report before the post-solve work, so a watchdog hard stop during
    // it cannot replace a proved status with TimeLimit.
    if let Some(tracker) = &tracker
        && !tracker.claim_final()
    {
        loop {
            std::thread::sleep(std::time::Duration::from_secs(1));
        }
    }
    // Infeasible set, ray, violated row (after the solve, bounded budget: never more than
    // what is left of the user's --time-limit)
    let post = crate::adopt::post_solve(
        path,
        // the infeasible set and the ray need the linearized model (it carries the general
        // constraints); the violation check needs the user's rows
        match &gc_plan {
            Some(p)
                if !matches!(
                    status,
                    HighsModelStatus::kInfeasible
                        | HighsModelStatus::kUnbounded
                        | HighsModelStatus::kUnboundedOrInfeasible
                ) =>
            {
                &p.user_lp
            }
            _ => &highs.model_lp_,
        },
        status,
        &highs.solution_,
        warn,
        crate::adopt::budget_for(secs, time_limit),
    );
    crate::adopt::print_warnings(&post);
    diags.extend(post);
    let solved = crate::envelope::Solved {
        path,
        lp: gc_plan.as_ref().map_or(&highs.model_lp_, |p| &p.user_lp),
        status,
        run_status,
        info,
        solution: &highs.solution_,
        time_s: secs,
        mip_stats: highs.mip_stats_.as_ref(),
        seed,
        threads: highs.get_options().threads,
        time_limit,
        node_limit,
        mip_max_nodes: highs.get_options().mip_max_nodes,
        options: &settings,
        include_solution,
        // The component's counters, only when a --wasm-hook component ran
        wasm: wasm_report.map(|(name, st)| {
            let s = st.lock().unwrap_or_else(|e| e.into_inner()).clone();
            crate::envelope::wasm_stats_json(&name, &s)
        }),
    };
    // warm.rs reports (read side): only when one of the flags was given
    let mut extras: Vec<(&str, serde_json::Value)> = Vec::new();
    if let Some(rep) = &start_report {
        let final_obj = has_solution.then_some(info.objective_function_value);
        extras.push((
            "start",
            crate::warm::start_json(
                rep,
                true,
                highs.model_lp_.sense_ == ObjSense::kMaximize,
                final_obj,
                &highs.solution_.col_value,
            ),
        ));
    }
    if let Some(rep) = &add_rows_report {
        extras.push(("added_rows", crate::warm::add_rows_json(rep)));
    }
    if let Some(o) = &gc_outcome {
        extras.push(("general_constraints", o.report.clone()));
    }
    if let Some(b) = basis_report {
        extras.push(("basis", b));
    }
    // Whoever claims the final report first (this thread, or the watchdog at its hard
    // deadline) is the only one that prints it.
    let mut tracked_env = None;
    if let Some(tracker) = &tracker {
        let mut env = crate::envelope::envelope(&solved);
        tracker.annotate(&mut env);
        tracker.emit_final(&env);
        tracked_env = Some(env);
    }
    if json {
        let mut text = match &tracked_env {
            Some(env) => {
                serde_json::to_string(env).unwrap_or_else(|e| format!("{{\"error\":{e:?}}}"))
            }
            None => crate::envelope::envelope_json(&solved),
        };
        if !extras.is_empty()
            && let Ok(serde_json::Value::Object(mut m)) = serde_json::from_str(&text)
        {
            for (k, v) in &extras {
                m.insert(k.to_string(), v.clone());
            }
            text = serde_json::Value::Object(m).to_string();
        }
        println!("{}", crate::adopt::with_diagnostics(&text, &diags));
    } else if csv {
        println!("status,objective,bound,nodes,lp_iters,time_s,x");
        let xs: Vec<String> = if has_solution {
            highs
                .solution_
                .col_value
                .iter()
                .map(|v| format!("{v}"))
                .collect()
        } else {
            Vec::new()
        };
        println!(
            "{},{},{},{},{},{secs},[{}]",
            if run_status == HighsStatus::kError {
                "Error".to_string()
            } else {
                csv_status(status)
            },
            objective,
            info.mip_dual_bound,
            info.mip_node_count,
            info.simplex_iteration_count,
            xs.join(";")
        );
    } else {
        println!("status:    {}", util_model_status_to_string(status));
        println!("objective: {objective}");
        println!("bound:     {}", info.mip_dual_bound);
        println!(
            "nodes:     {} ({} LP iterations, {secs:.3}s)",
            info.mip_node_count, info.simplex_iteration_count
        );
        for (k, v) in &extras {
            println!("{k}: {v}");
        }
    }
    if csv {
        // the CSV contract keeps stdout; the warm-start reports go to stderr there
        for (k, v) in &extras {
            eprintln!("{k}: {v}");
        }
    }
    // After the result lines, so stdout is complete whatever happens to the
    // file. `.json` writes the envelope's `solution` object, anything else the
    // HiGHS raw solution file with the real names.
    if let Some(sol_path) = &solution_file {
        let written = if sol_path.ends_with(".json") {
            std::fs::write(sol_path, crate::envelope::solution_file_json(&solved))
                .map_err(|e| e.to_string())
        } else {
            crate::envelope::write_solution_file(
                sol_path,
                highs.get_options().log_options.clone(),
                &solved,
            )
            .map_err(|e| e.to_string())
        };
        if let Err(e) = written {
            // the envelope is already on stdout: the error stays on stderr, one JSON document
            return crate::diag::fail_msg(
                false,
                "RS0090",
                4,
                format!("rusolve mip: cannot write {sol_path}: {e}"),
            );
        }
    }
    if run_status == HighsStatus::kError {
        ExitCode::from(4)
    } else {
        ExitCode::SUCCESS
    }
}

/// State of the trace callback (mip_trace.cpp `TraceState`).
#[derive(Default)]
struct TraceState {
    have_last: bool,
    last_nodes: i64,
    last_db: f64,
    last_pb: f64,
}

/// C `%.17g`.
fn fmt_g17(v: f64) -> String {
    rusolve_ekk::util::c_printf::c_fmt_g(v, 17)
}

/// C `atof`: the longest numeric prefix of `s` after leading whitespace (0.0 when there is
/// none), with the `inf`, `infinity` and `nan` words as strtod reads them.
fn c_atof(s: &str) -> f64 {
    let t = s.trim_start();
    let (sign, body) = match t.as_bytes().first() {
        Some(b'-') => (-1.0, &t[1..]),
        Some(b'+') => (1.0, &t[1..]),
        _ => (1.0, t),
    };
    let lower = body.to_ascii_lowercase();
    if lower.starts_with("inf") {
        return sign * f64::INFINITY;
    }
    if lower.starts_with("nan") {
        return f64::NAN;
    }
    let b = body.as_bytes();
    let mut i = 0;
    let mut digits = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
        digits += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return 0.0;
    }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let exp_start = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > exp_start {
            i = j;
        }
    }
    sign * body[..i].parse::<f64>().unwrap_or(0.0)
}

/// mip_trace.cpp `stripTimes`: replace every `<digits>[.<digits>]s` token by `<t>s`.
fn strip_times(s: &str) -> String {
    let b = s.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i].is_ascii_digit() && (i == 0 || !b[i - 1].is_ascii_alphanumeric()) {
            let mut j = i;
            while j < b.len() && (b[j].is_ascii_digit() || b[j] == b'.') {
                j += 1;
            }
            if j < b.len()
                && b[j] == b's'
                && (j + 1 == b.len() || !b[j + 1].is_ascii_alphanumeric())
            {
                out.extend_from_slice(b"<t>s");
                i = j + 1;
                continue;
            }
            out.extend_from_slice(&b[i..j]);
            i = j;
            continue;
        }
        out.push(b[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// mip_trace.cpp `keepLogLine`: the selected HiGHS log lines.
fn keep_log_line(line: &str) -> bool {
    const KEYS: [&str; 23] = [
        "Presolv",
        "reductions",
        "rows,",
        "Restart",
        "Solving report",
        "Status",
        "Primal bound",
        "Dual bound",
        "Gap",
        "Solution status",
        "Nodes",
        "LP iterations",
        "Objective",
        "Coefficient ranges",
        "Cost",
        "Bound",
        "RHS",
        "Matrix",
        "MIP has",
        "symmetr",
        "Symmetr",
        "clique",
        "Clique",
    ];
    if line.contains("Timing") {
        return false;
    }
    // the timing block of the solving report
    if line.contains("(Presolve)") || line.contains("(Solve)") || line.contains("(Postsolve)") {
        return false;
    }
    KEYS.iter().any(|k| line.contains(k))
}

/// mip_trace.cpp `callback`: prints the `log`, `cb` and `sol` lines.
fn mip_trace_callback(
    st: &mut TraceState,
    kind: HighsInt,
    message: &str,
    out: &HighsCallbackOutput,
) {
    match kind {
        x if x == kCallbackLogging => {
            // a log message may hold several lines
            let text = message.trim_end_matches(['\n', '\r']);
            for line in text.split('\n') {
                if !line.is_empty() && keep_log_line(line) {
                    println!("log {}", strip_times(line));
                }
            }
        }
        x if x == kCallbackMipInterrupt => {
            let nodes = out.mip_node_count;
            if st.have_last
                && nodes == st.last_nodes
                && out.mip_dual_bound.to_bits() == st.last_db.to_bits()
                && out.mip_primal_bound.to_bits() == st.last_pb.to_bits()
            {
                return;
            }
            st.have_last = true;
            st.last_nodes = nodes;
            st.last_db = out.mip_dual_bound;
            st.last_pb = out.mip_primal_bound;
            println!(
                "cb nodes={} lpit={} db={} pb={}",
                nodes,
                out.mip_total_lp_iterations,
                fmt_g17(out.mip_dual_bound),
                fmt_g17(out.mip_primal_bound)
            );
        }
        x if x == kCallbackMipImprovingSolution => {
            println!(
                "sol nodes={} lpit={} obj={}",
                out.mip_node_count,
                out.mip_total_lp_iterations,
                fmt_g17(out.objective_function_value)
            );
        }
        _ => {}
    }
}

/// `rusolve mip-trace <model.mps> [time_limit_s]`: the trace of
/// tools/highs-oracle/mip_trace.cpp from the port: install a `HighsCallback::user_callback` that prints the
/// `log` (filtered, timings stripped), `cb` (deduplicated) and `sol` lines
/// exactly as mip_trace.cpp `callback`, start the logging / MIP interrupt /
/// improving-solution callbacks, run, then print the `final` line with
/// `%.17g` doubles.
///
/// Argument rule (mip_trace.cpp `main`): `args[1]` is the time limit only when it has no
/// `=`; every later `key=value` goes through `Highs::setOptionValue(key, value)` (string
/// overload, after the defaults, exit 1 on `kError`). Defaults are mip_trace.cpp's:
/// `output_flag = true`, `log_to_console = false` (the log reaches stdout only through
/// the `kCallbackLogging` callback), `threads = 1`, `parallel = off`, `random_seed = 0`,
/// `time_limit`.
pub fn run_mip_trace(args: &[String]) -> ExitCode {
    let Some(path) = args.first() else {
        eprintln!("usage: rusolve mip-trace <model.mps> [time_limit_s]");
        return ExitCode::from(2);
    };
    // argv[2] is the time limit when it has no '='; the rest are key=value.
    let (time_limit, first_option) = match args.get(1) {
        Some(arg) if !arg.contains('=') => (c_atof(arg), 2),
        _ => (3600.0, 1),
    };
    // mip_trace.cpp is the HiGHS 1.15.1 oracle: HiGHS 1.16 options off unless set below
    let mut highs = Highs::with_options(HighsOptions::default());
    // mip_trace.cpp defaults: the log goes to the callback only
    highs.set_option_value_bool("output_flag", true);
    highs.set_option_value_bool("log_to_console", false);
    highs.set_option_value_int("threads", 1);
    // mip_trace.cpp pins parallel=off for HiGHS 1.16 (#3406); under 1.15.1
    // semantics it is the same search as the default choose
    highs.set_option_value_str("parallel", "off");
    highs.set_option_value_int("random_seed", 0);
    // the oracle propagates continuous bounds with the HiGHS 0.3 threshold
    highs.set_option_value_double("mip_cont_prop_improve", 0.3);
    highs.set_option_value_double("mip_pscost_variance_tol", 0.0);
    // structure levers are default-on (they act only when the structure exists); the oracle has neither
    highs.set_option_value_bool("mip_presolve_xor_gf2", false);
    highs.set_option_value_bool("mip_cdcl_uip", false);
    highs.set_option_value_bool("mip_zero_obj_cdcl", false);
    highs.set_option_value_bool("mip_activation_dff", false);
    highs.set_option_value_bool("mip_dff_separator", false);
    highs.set_option_value_bool("mip_presolve_clique_varbound", false);
    highs.set_option_value_bool("mip_orbitope_root_fix", false);
    highs.set_option_value_bool("mip_ti_precedence_sepa", false);
    highs.set_option_value_bool("mip_ti_list_sched", false);
    highs.set_option_value_bool("mip_improve_early", false);
    highs.set_option_value_bool("mip_first_inc", false);
    highs.set_option_value_bool("mip_first_inc_root", false);
    highs.set_option_value_bool("mip_netdesign_sepa", false);
    highs.set_option_value_bool("mip_netdesign_heur", false);
    highs.set_option_value_bool("mip_netdesign_cutset", false);
    highs.set_option_value_bool("mip_root_bound_stall", false);
    highs.set_option_value_bool("mip_lattice", false);
    highs.set_option_value_bool("mip_lattice_kernel", false);
    highs.set_option_value_bool("mip_heur_lp_quota", false);
    highs.set_option_value_double("time_limit", time_limit);
    for arg in &args[first_option..] {
        let Some(eq) = arg.find('=').filter(|&e| e > 0) else {
            eprintln!("mip-trace: expected key=value, got '{arg}'");
            return ExitCode::from(2);
        };
        let (key, value) = (&arg[..eq], &arg[eq + 1..]);
        if highs.set_option_value_str(key, value) == HighsStatus::kError {
            eprintln!("mip-trace: unknown option or invalid value: {key}={value}");
            return ExitCode::from(1);
        }
    }
    let lp = match read_mip(path) {
        Ok(lp) => lp,
        Err(e) => {
            eprintln!("mip-trace: cannot read {path}: {e}");
            return ExitCode::from(1);
        }
    };
    if highs.pass_model(lp) == HighsStatus::kError {
        eprintln!("mip-trace: cannot read {path}");
        return ExitCode::from(1);
    }
    let state = Arc::new(Mutex::new(TraceState::default()));
    let trace = Arc::clone(&state);
    let callback = HighsCallbackFunctionType::new(
        move |kind: HighsInt,
              message: &str,
              out: &HighsCallbackOutput,
              _cb_in: &mut HighsCallbackInput| {
            mip_trace_callback(&mut trace.lock().unwrap(), kind, message, out);
        },
    );
    highs.set_callback(callback);
    // C++ ignores the return values of these calls
    highs.start_callback(kCallbackLogging);
    highs.start_callback(kCallbackMipInterrupt);
    highs.start_callback(kCallbackMipImprovingSolution);
    highs.run();
    let info = highs.get_info();
    println!(
        "final status={} obj={} bound={} nodes={} lp_iters={}",
        util_model_status_to_string(highs.get_model_status()),
        fmt_g17(info.objective_function_value),
        fmt_g17(info.mip_dual_bound),
        info.mip_node_count,
        info.simplex_iteration_count
    );
    ExitCode::SUCCESS
}

/// `rusolve presolve-dump <model.mps> [key=value ...]`: run the port's
/// `Highs::presolve` (the MIP presolve for a MIP) at the mip_trace settings
/// and print the postsolve-stack reduction sequence and the presolved model
/// (all doubles %.17g), in the format of the C++ presolve_dump driver used
/// to find the first presolve divergence against the HiGHS 1.16 oracle
/// (tools/highs-oracle).
pub fn run_presolve_dump(args: &[String]) -> ExitCode {
    let Some(path) = args.first() else {
        eprintln!("usage: rusolve presolve-dump <model.mps> [key=value ...]");
        return ExitCode::from(2);
    };
    let mut highs = Highs::with_options(HighsOptions::default());
    highs.set_option_value_bool("output_flag", false);
    highs.set_option_value_int("threads", 1);
    highs.set_option_value_str("parallel", "off");
    highs.set_option_value_int("random_seed", 0);
    // the oracle propagates continuous bounds with the HiGHS 0.3 threshold
    highs.set_option_value_double("mip_cont_prop_improve", 0.3);
    highs.set_option_value_double("mip_pscost_variance_tol", 0.0);
    // structure levers are default-on (they act only when the structure exists); the oracle has neither
    highs.set_option_value_bool("mip_presolve_xor_gf2", false);
    highs.set_option_value_bool("mip_cdcl_uip", false);
    highs.set_option_value_bool("mip_zero_obj_cdcl", false);
    highs.set_option_value_bool("mip_activation_dff", false);
    highs.set_option_value_bool("mip_dff_separator", false);
    highs.set_option_value_bool("mip_presolve_clique_varbound", false);
    highs.set_option_value_bool("mip_orbitope_root_fix", false);
    highs.set_option_value_bool("mip_ti_precedence_sepa", false);
    highs.set_option_value_bool("mip_ti_list_sched", false);
    highs.set_option_value_bool("mip_improve_early", false);
    highs.set_option_value_bool("mip_first_inc", false);
    highs.set_option_value_bool("mip_first_inc_root", false);
    highs.set_option_value_bool("mip_netdesign_sepa", false);
    highs.set_option_value_bool("mip_netdesign_heur", false);
    highs.set_option_value_bool("mip_netdesign_cutset", false);
    highs.set_option_value_bool("mip_root_bound_stall", false);
    highs.set_option_value_bool("mip_lattice", false);
    highs.set_option_value_bool("mip_lattice_kernel", false);
    highs.set_option_value_bool("mip_heur_lp_quota", false);
    for arg in &args[1..] {
        let Some(eq) = arg.find('=').filter(|&e| e > 0) else {
            eprintln!("presolve-dump: expected key=value, got '{arg}'");
            return ExitCode::from(2);
        };
        let (key, value) = (&arg[..eq], &arg[eq + 1..]);
        if highs.set_option_value_str(key, value) == HighsStatus::kError {
            eprintln!("presolve-dump: unknown option or invalid value: {key}={value}");
            return ExitCode::from(1);
        }
    }
    let lp = match read_mip(path) {
        Ok(lp) => lp,
        Err(e) => {
            eprintln!("presolve-dump: cannot read {path}: {e}");
            return ExitCode::from(1);
        }
    };
    if highs.pass_model(lp) == HighsStatus::kError {
        eprintln!("presolve-dump: cannot read {path}");
        return ExitCode::from(1);
    }
    highs.presolve();
    let ps = &highs.presolve_.data_.post_solve_stack;
    let types: Vec<_> = ps.debug_reduction_types().collect();
    println!("nred {}", types.len());
    // PD_DETAIL: append the rows/cols each DuplicateRow/RedundantRow/FixedCol recorded
    let detail = std::env::var_os("PD_DETAIL").is_some();
    for (i, t) in types.iter().enumerate() {
        let name = format!("{t:?}");
        let extra = if detail {
            ps.debug_reduction_detail(i)
        } else {
            String::new()
        };
        println!("red {i} {}{extra}", name.trim_start_matches('k'));
    }
    let lp = highs.get_presolved_lp();
    let mut a = lp.a_matrix_.clone();
    a.ensure_colwise();
    println!(
        "dims rows={} cols={} nnz={} offset={}",
        lp.num_row_,
        lp.num_col_,
        a.num_nz(),
        fmt_g17(lp.offset_)
    );
    let join = |v: &[HighsInt]| v.iter().map(|x| format!(" {x}")).collect::<String>();
    println!("origcol{}", join(ps.get_orig_cols_index()));
    println!("origrow{}", join(ps.get_orig_rows_index()));
    for j in 0..lp.num_col_ as usize {
        let int = lp.integrality_.get(j).map_or(0, |&t| t as i32);
        let mut line = format!(
            "col {j} c={} l={} u={} i={int} :",
            fmt_g17(lp.col_cost_[j]),
            fmt_g17(lp.col_lower_[j]),
            fmt_g17(lp.col_upper_[j])
        );
        for k in a.start_[j] as usize..a.start_[j + 1] as usize {
            line.push_str(&format!(" {}:{}", a.index_[k], fmt_g17(a.value_[k])));
        }
        println!("{line}");
    }
    for i in 0..lp.num_row_ as usize {
        println!(
            "row {i} l={} u={}",
            fmt_g17(lp.row_lower_[i]),
            fmt_g17(lp.row_upper_[i])
        );
    }
    ExitCode::SUCCESS
}
