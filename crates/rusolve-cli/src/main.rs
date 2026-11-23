//! `rusolve solve model.mps [--csv]` — thin CLI over the LP pipeline;
//! `rusolve mip …` for MIPs and `rusolve classify …` for instance structure.
//!
//! The solver-space LP enters through the `SparseLp` front door
//! (`Model::to_sparse_lp`): the MPS reader's CSC matrix is carried
//! through without a dense intermediate.

use rusolve_io::RowSense;
#[cfg(test)]
use rusolve_io::read_mps;
#[cfg(test)]
use rusolve_lp::Lp;
use rusolve_lp::{PfiConfig, PfiEngine, SimplexConfig, SimplexSolver};
use std::process::ExitCode;

pub mod adopt;
pub mod api;
mod check;
mod check_conic;
mod check_qp;
mod cli_core;
mod commands;
mod conic_solve;
mod convert;
mod diag;
mod diagnose;
mod ekk;
mod ekk_seq;
mod env_registry;
mod envelope;
mod explain;
mod explain_mip;
mod feasrelax;
mod features;
mod gc_solve;
mod lint;
mod mcp;
mod mcp_build;
mod mcp_resources;
mod mcp_tools;
mod mip_highs;
mod options_cmd;
mod profile;
mod progress;
mod qp;
mod quad_convexity;
mod ranging;
mod route;
mod tune;
mod tune_race;
mod warm;

/// The HiGHS options every CLI solve starts from: the best-measured
/// [`HighsOptions::recommended`] profile, or HiGHS 1.15.1 semantics (every HiGHS 1.16
/// option off) under `RUSOLVE_PROFILE=highs-1.15.1`.
pub(crate) fn base_options() -> rusolve_ekk::lp_data::highs_options::HighsOptions {
    use rusolve_ekk::lp_data::highs_options::HighsOptions;
    match std::env::var("RUSOLVE_PROFILE").as_deref() {
        Ok("highs-1.15.1") => HighsOptions::default(),
        Ok("recommended") | Err(_) => HighsOptions::recommended(),
        Ok(other) => {
            eprintln!("RUSOLVE_PROFILE={other}: expected recommended or highs-1.15.1");
            std::process::exit(2);
        }
    }
}

/// Test convenience: the solver-space [`Lp`] of a model (the CLI solve path
/// itself enters through `Model::to_sparse_lp` + `solve_sparse`).
#[cfg(test)]
fn model_to_lp(m: &rusolve_io::Model) -> Result<Lp, String> {
    // SparseLp front door; the row-form `Lp` comes back via the expansion
    // (no fp arithmetic, no dense intermediate).
    Ok(m.to_sparse_lp()?.to_lp())
}

pub fn run(args: &[String]) -> ExitCode {
    // `rusolve --explain RSnnnn` (docs/design/diagnostics.md)
    if let Some(code) = diag::explain(args) {
        return code;
    }
    if let Some(cmd @ ("solve" | "mip")) = args.first().map(String::as_str)
        && let Some(path) = args.get(1)
        && !path.starts_with("--")
        && !args.iter().any(|a| a == "--help" || a == "-h")
    {
        // One classifier and routing table for both commands (route.rs). `None` is a plain
        // LP or MIP, or a class that keeps its legacy path: it runs exactly as before.
        if let Some(code) = route::dispatch(cmd, path, &args[1..]) {
            return code;
        }
    }
    let may_go_conic = matches!(args.first().map(String::as_str), Some("solve" | "mip"));
    diag::defer_conic_json(may_go_conic);
    let code = run_legacy(args);
    diag::defer_conic_json(false);
    // Safety net: a legacy reader refused a conic model the pre-scan did not announce (and said
    // so on stderr); take it to the Clarabel route anyway. A model with neither cones nor
    // quadratic rows never sets the flag, so its path is the pre-quadratic-rows one.
    if code != ExitCode::SUCCESS
        && rusolve_io::conic::conic_rejected()
        && let Some(cmd @ ("solve" | "mip")) = args.first().map(String::as_str)
    {
        route::note(route::Class::Conic);
        return conic_solve::run(cmd, &args[1..]);
    }
    code
}

fn run_legacy(args: &[String]) -> ExitCode {
    // `--help` / `-h` / `<cmd> --help` (and `--json`) print the registry's usage and exit 0,
    // before any command reads its arguments as a path.
    if let Some(code) = commands::handle_help(args) {
        return code;
    }
    let json = diag::wants_json(args);
    let mut args = args.iter();
    let Some(cmd) = args.next() else {
        eprint!("{}", commands::usage_text());
        return ExitCode::from(2);
    };
    if cmd == "env" {
        return run_env(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "dump" {
        let Some(path) = args.next() else {
            eprintln!("missing MPS path");
            return ExitCode::FAILURE;
        };
        let m = match rusolve_io::read_model_file(path) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("read {path}: {e}");
                return ExitCode::FAILURE;
            }
        };
        if let Err(e) = m.gc.reject("rusolve dump") {
            eprintln!("{e}");
            return ExitCode::from(3);
        }
        // machine-readable model dump for differential testing
        println!("n_vars {}", m.num_variables);
        println!("n_cons {}", m.num_constraints);
        println!("obj_sense {:?}", m.obj_sense);
        for (j, c) in m.objective.iter().enumerate() {
            if *c != 0.0 {
                println!("obj {j} {c:?}");
            }
        }
        for i in 0..m.num_constraints {
            println!(
                "sense {i} {:?}",
                match m.senses[i] {
                    RowSense::LessEqual => 'L',
                    RowSense::GreaterEqual => 'G',
                    RowSense::Equal => 'E',
                }
            );
            if m.rhs[i] != 0.0 {
                println!("rhs {i} {:?}", m.rhs[i]);
            }
            if m.ranges[i] != 0.0 {
                println!("range {i} {:?}", m.ranges[i]);
            }
        }
        for j in 0..m.num_variables {
            for k in m.matrix.col_ptr()[j]..m.matrix.col_ptr()[j + 1] {
                println!(
                    "nz {j} {} {:?}",
                    m.matrix.row_idx()[k],
                    m.matrix.values()[k]
                );
            }
        }
        for (j, b) in m.bounds.iter().enumerate() {
            println!("bound {j} {:?} {:?}", b.lower, b.upper);
        }
        return ExitCode::SUCCESS;
    }
    if cmd == "diagnose" {
        // Formulation or data? (crates/rusolve-cli/src/diagnose.rs)
        return diagnose::run(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "explain" {
        return explain::run(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "convert" {
        return convert::run(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "lint" {
        // Static model lint (crates/rusolve-cli/src/lint.rs)
        return lint::run(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "options" {
        // The option registry, so an agent can see what it may change.
        return options_cmd::run_options(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "profile" {
        // Versioned profile files (crates/rusolve-cli/src/profile.rs)
        return profile::run(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "tune" {
        // Tuning harness (crates/rusolve-cli/src/tune.rs)
        return tune::run(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "probe" {
        // Dynamic instance features (crates/rusolve-cli/src/features.rs)
        return features::run_probe(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "mcp" {
        // Model Context Protocol server over stdio (crates/rusolve-cli/src/mcp.rs)
        return mcp::run(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "classify" {
        return run_classify(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "ekk-seq" {
        // Acceptance harness: the ported Highs API replaying
        // tools/highs-oracle/api_seq.cpp.
        let Some(path) = args.next() else {
            eprintln!("usage: rusolve ekk-seq <model.mps>");
            return ExitCode::FAILURE;
        };
        return ekk_seq::run_ekk_seq(path);
    }
    if cmd == "check" {
        // solution checker
        return check::run_check(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "feasrelax" {
        // cheapest feasibility relaxation, by name.
        return feasrelax::run_feasrelax(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "mip-trace" {
        // Acceptance harness: the MIP port replaying
        // tools/highs-oracle/mip_trace.cpp.
        return mip_highs::run_mip_trace(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "presolve-dump" {
        // presolve reduction sequence vs the HiGHS 1.16 C++ oracle (tools/highs-oracle)
        return mip_highs::run_presolve_dump(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "ranging" {
        // sensitivity by name on the port
        return ranging::run_ranging(&args.map(|a| a.to_string()).collect::<Vec<_>>());
    }
    if cmd == "mip" {
        let args = args.map(|a| a.to_string()).collect::<Vec<_>>();
        // the rusolve-ekk HiGHS MIP port is the default (RUSOLVE_MIP_HIGHS=0: native search).
        if mip_highs::enabled() {
            return mip_highs::run_mip_highs(&args);
        }
        return run_mip(&args);
    }
    if cmd != "solve" {
        return diag::unknown_command(json, cmd);
    }
    let Some(path) = args.next() else {
        return diag::fail_msg(
            json,
            "RS0006",
            2,
            "rusolve solve: missing model path (see `rusolve solve --help`)",
        );
    };
    // The LP clock starts before the read, as HiGHS's run clock does.
    let t_start = std::time::Instant::now();
    let flags: Vec<String> = args.map(|a| a.to_string()).collect();
    // The same strictness as `mip`: an unknown flag or a bad value is a usage error (2),
    // not silently ignored. `--format`, `--solution` and `--include-solution` are shared.
    let fmt = match parse_solve_flags(&flags) {
        Ok(f) => f,
        Err(e) => {
            let mut d = rusolve_core::diagnostic::Diagnostic::error(
                if e.starts_with("unknown argument") {
                    "RS0002"
                } else {
                    "RS0003"
                },
                format!("rusolve solve: {e}"),
            );
            if let Some(bad) = e
                .strip_prefix("unknown argument '")
                .and_then(|r| r.split('\'').next())
                && let Some(c) = commands::find("solve")
            {
                // The closest flag of `solve`
                d.help.extend(rusolve_core::diagnostic::did_you_mean(
                    bad,
                    c.flags.iter().map(|f| f.name),
                ));
            }
            return diag::fail(json, d, 2);
        }
    };
    let csv = fmt.csv;
    let json = fmt.json;
    let include_solution = fmt.include_solution;
    // Iteration policy: a stand-alone
    // LP solve has NO pivot cap by default and is bounded by the wall clock
    // instead — HiGHS's own policy: `simplex_iteration_limit` defaults to
    // kHighsIInf "for LPs, but not subproblems in the MIP solver" and
    // `time_limit` governs (lp_data/HighsOptions.h, MIT; policy, not code).
    // The library default (`SimplexConfig::max_iterations` = 100k) stays:
    // it is the MIP search's per-LP backstop. The old 100k cap killed fome13
    // (HiGHS: 139,647 iters) and qap15 (292k) by policy, not by algorithm.
    // `--iter-limit N` restores an explicit cap (100000 = the earlier CLI default).
    let flag_val = |name: &str| -> Option<&str> {
        flags
            .iter()
            .position(|a| a == name)
            .and_then(|k| flags.get(k + 1))
            .map(|v| v.as_str())
    };
    let time_limit: Option<f64> = match flag_val("--time-limit") {
        Some(v) => match v.parse::<f64>() {
            Ok(t) if t >= 0.0 && t.is_finite() => Some(t),
            _ => {
                return diag::fail_msg(
                    json,
                    "RS0003",
                    2,
                    format!("--time-limit: want a non-negative number of seconds, got '{v}'"),
                );
            }
        },
        None => None,
    };
    let max_iterations: usize = match flag_val("--iter-limit") {
        Some(v) => match v.parse::<usize>() {
            Ok(n) => n,
            Err(_) => {
                return diag::fail_msg(
                    json,
                    "RS0003",
                    2,
                    format!("--iter-limit: want a pivot count, got '{v}'"),
                );
            }
        },
        // `usize::MAX / 2`: the same "uncapped" sentinel the MIP root uses
        // under a deadline (no `iters - start` overflow in the loops).
        None => usize::MAX / 2,
    };
    let no_presolve = flags.iter().any(|a| a == "--no-presolve");
    let duals = flags.iter().any(|a| a == "--duals");

    // rusolve-ekk engine (faithful port of the HiGHS serial simplex through
    // the ported Highs::run driver, no presolve; crates/rusolve-ekk). The
    // default since a Netlib + LP-47 benchmark (LP-47 sgm 0.174 s vs
    // 0.193 s legacy). It reads the file itself (HiGHS
    // MPS semantics: file-order columns, first-duplicate-wins). `--ekk` is
    // kept as a no-op for old scripts.
    if use_ekk_engine(&flags) {
        return ekk::solve_ekk_cli(
            path,
            csv,
            duals,
            time_limit,
            flag_val("--iter-limit").map(|_| max_iterations),
            ekk::JsonOut {
                json,
                include_solution,
            },
            flag_val("--solution-file"),
        );
    }
    let model = match rusolve_io::read_model_file(path) {
        Ok(m) => m,
        Err(e) => {
            return diag::fail_read(json, path, &e.to_string());
        }
    };
    if !model.hessian.is_empty() {
        eprintln!(
            "rusolve solve: {}",
            qp::unsupported("solve (legacy engine, RUSOLVE_EKK=0)")
        );
        return ExitCode::from(3);
    }
    let sp = match model.to_sparse_lp() {
        Ok(sp) => sp,
        Err(e) => {
            eprintln!("convert: {e}");
            return ExitCode::from(3);
        }
    };
    // Interior-point path: opt-in only; the
    // default simplex path below is untouched.
    if flags.iter().any(|a| a == "--ipm") || std::env::var("RUSOLVE_IPM").is_ok_and(|v| v != "0") {
        return solve_ipm_cli(
            &sp,
            &model,
            path,
            no_presolve,
            csv,
            duals,
            time_limit,
            t_start,
            (json, include_solution),
        );
    }
    let mut solver = SimplexSolver::new(SimplexConfig {
        max_iterations,
        presolve: !no_presolve,
        // Sparse solve-path lever:
        // RUSOLVE_SPARSE=1 runs the direct CSC-only `Standard` build when
        // presolve/dual are off (both consume the `Lp` frame); otherwise
        // `solve_sparse` delegates via the (bit-identical) `Lp` hop. Lever
        // off is the same delegate path. The days of "build the dense Lp
        // first, then convert" are over: the MPS reader's CSC reaches
        // the solver untouched.
        sparse: std::env::var("RUSOLVE_SPARSE").is_ok(),
        // SH pricing is DEFAULT ON (A/B: -37% iters, wall tie,
        // scsd6 recovered); RUSOLVE_NO_SH=1 is the lever back to Dantzig.
        sh_pricing: std::env::var("RUSOLVE_NO_SH").is_err(),
        // Implied-equality conversion: re-measured with DTE v2 in
        // (the co-landing condition) and still harmful — scrs8 p2 312→1827,
        // e226 +51% iters. RUSOLVE_PS_IMPEQ=1 stays the A/B lever, default off.
        presolve_implied_eq: std::env::var("RUSOLVE_PS_IMPEQ").is_ok(),
        // Power-of-2 equilibration scaling: DEFAULT OFF behind the
        // RUSOLVE_SCALE=1 lever (Netlib A/B: tuff -77%, scrs8 -219 iters; but
        // e226 577→3840 phase-1 stall — see SimplexConfig::scaling docs).
        scaling: std::env::var("RUSOLVE_SCALE").is_ok(),
        // Phase-2 entering tie perturbation: DEFAULT ON (fires only on
        // 4×bland_after Bland stalls ⇒ non-stalling runs bit-identical);
        // RUSOLVE_NO_P2PERT=1 = lever back.
        p2_tie_perturb: std::env::var("RUSOLVE_NO_P2PERT").is_err(),
        // Phase-1 re-perturbation: DEFAULT ON (fires only on an already
        // fired eps whose w probe is bit-frozen across a stall window);
        // RUSOLVE_NO_P1REPERT=1 is the lever back (sticky-eps crawl).
        p1_re_perturb: std::env::var("RUSOLVE_NO_P1REPERT").is_err(),
        // General-scale DupCol merge:
        // DEFAULT ON after an A/B (PRE-bucket sgm it_r
        // 1.317->0.996, net -665 iters on LP-47, Netlib 5226->5179);
        // RUSOLVE_PS_NO_PARCOL=1 is the lever back.
        presolve_parcol: std::env::var("RUSOLVE_PS_NO_PARCOL").is_err(),
        presolve_agg: std::env::var("RUSOLVE_PRESOLVE_AGG").as_deref() != Ok("0"),
        presolve_parcol_dom: std::env::var("RUSOLVE_PS_NO_PARCOLDOM").is_err(),
        presolve_vub_eq: std::env::var("RUSOLVE_PS_NO_VUBEQ").is_err(),
        // Dual simplex engine: DEFAULT ON for `solve` since a 47
        // MIPLIB-relaxation benchmark (47/47
        // Optimal vs primal 44/47, 0 objective mismatches, 0.38x primal
        // time). RUSOLVE_PRIMAL=1 = lever back; the MIP search keeps its own
        // `MipConfig::lp` (primal cold, `--dual-lp` opt-in).
        dual: std::env::var("RUSOLVE_PRIMAL").is_err(),
        // Dual accuracy regime:
        // default ON in the engine; RUSOLVE_DUAL_ACC=0 = lever off.
        dual_acc: std::env::var("RUSOLVE_DUAL_ACC").is_ok_and(|v| v != "0"),
        // Policy arm: with RUSOLVE_PRIMAL=1, route cold primal-
        // infeasible starts to the dual engine instead of paying the
        // artificial-basis phase 1 (measured: phase 1 owns 63-100 % of
        // primal iterations on all 19 LP-47 blowups; HiGHS pays no primal
        // phase 1 anywhere — dual chosen 46/46). Default OFF: a cold dual start
        // lost on the MIP benchmark tier.
        p1_route_dual: std::env::var("RUSOLVE_P1_ROUTE_DUAL").is_ok(),
        // Stall-threshold probe lever (diagnostics: does Bland CONTAIN or
        // CAUSE the tuff phase-1 grind? huge value = Bland never engages).
        bland_after: std::env::var("RUSOLVE_BLAND_AFTER")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(25),
        // Bounded Bland episode budget (iterations; released early on
        // bland_after non-degenerate pivots or the randomized-deviation rung;
        // re-entry needs doubled stall thresholds: 25, 50, 100, ...).
        // RUSOLVE_BLAND_EPISODE=0 = legacy sticky Bland (A/B lever back).
        bland_episode: std::env::var("RUSOLVE_BLAND_EPISODE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(25),
        // Default flipped to Sparse (Netlib A/B: geo-mean 50.8×→31.3×);
        // RUSOLVE_DENSE_LU=1 is the A/B lever back to the old engine.
        pfi: PfiConfig {
            engine: if std::env::var("RUSOLVE_DENSE_LU").is_ok() {
                PfiEngine::Dense
            } else {
                PfiEngine::Sparse
            },
            ..PfiConfig::default()
        },
        ..SimplexConfig::default()
    });
    solver.set_deadline(time_limit.map(|t| t_start + std::time::Duration::from_secs_f64(t)));
    let sol = match solver.solve_sparse(&sp) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("solve: {e}");
            return ExitCode::from(4);
        }
    };
    // An MPS RHS entry on the objective (N) row is an objective
    // constant. HiGHS folds it into the *reported* objective only
    // (`HMpsFF.cpp:1065–1084` → `lp.offset_`, added at reporting by
    // `HighsInterface.cpp:3195`); the simplex itself optimizes raw `c·x`.
    // Shift the print/CSV line the same way — and only when a constant is
    // present, so every offset-free instance stays bit-identical.
    let reported_obj = if model.obj_offset != 0.0 {
        sol.objective + model.obj_offset
    } else {
        sol.objective
    };
    if json {
        let optimal = sol.status == rusolve_lp::SimplexStatus::Optimal;
        let native = envelope::NativeSolved {
            path,
            model: &model,
            status: match sol.status {
                rusolve_lp::SimplexStatus::Optimal => "Optimal",
                rusolve_lp::SimplexStatus::Infeasible => "Infeasible",
                rusolve_lp::SimplexStatus::Unbounded => "Unbounded",
                rusolve_lp::SimplexStatus::IterationLimit => "IterationLimit",
            },
            objective: reported_obj,
            // an LP solved to optimality has its objective as its bound
            dual_bound: if optimal { reported_obj } else { f64::NAN },
            nodes: 0,
            lp_iterations: (sol.phase1_iterations + sol.phase2_iterations) as i64,
            time_s: t_start.elapsed().as_secs_f64(),
            seed: 0,
            threads: 1,
            time_limit,
            node_limit: None,
            x: if optimal { &sol.x } else { &[] },
            include_solution,
        };
        println!("{}", envelope::native_envelope_json(&native));
    } else if csv {
        println!("status,objective,phase1_iters,phase2_iters,x");
        let xs: Vec<String> = sol.x.iter().map(|v| format!("{v}")).collect();
        println!(
            "{:?},{},{},{},[{}]",
            sol.status,
            reported_obj,
            sol.phase1_iterations,
            sol.phase2_iterations,
            xs.join(";")
        );
    } else {
        println!("status:  {:?}", sol.status);
        println!("objective: {}", reported_obj);
        println!(
            "iterations: {} (phase1) + {} (phase2){}",
            sol.phase1_iterations,
            sol.phase2_iterations,
            if sol.bland_triggered { " [bland]" } else { "" }
        );
        let first5: Vec<String> = sol.x.iter().take(5).map(|v| format!("{v:.6}")).collect();
        println!("x[:5]: [{}]", first5.join(", "));
    }
    if sol.status == rusolve_lp::SimplexStatus::IterationLimit
        && let Some(t) = time_limit
        && t_start.elapsed().as_secs_f64() >= t
    {
        eprintln!("rusolve: time limit {t} s reached");
    }
    // Opt-in dual/ray dump (never part of the --csv contract lines).
    if duals {
        if sol.duals_valid {
            println!("row_duals: {:?}", sol.row_duals);
            println!("reduced_costs: {:?}", sol.reduced_costs);
        } else {
            println!("duals: unavailable (presolve-ON path or non-optimal status)");
        }
        if let Some(d) = &sol.primal_ray {
            println!("primal_ray: {d:?}");
        }
        if let Some(p) = &sol.dual_ray {
            println!("dual_ray: {p:?}");
        }
    }
    // The solve ran; the status is in the output, not in the exit code.
    ExitCode::SUCCESS
}

/// What `rusolve solve`'s flags ask of the output.
struct SolveFmt {
    csv: bool,
    json: bool,
    include_solution: bool,
}

/// Validates every flag of `rusolve solve` against `commands::SOLVE_FLAGS` semantics: an unknown
/// flag, a missing or unparsable value and `--csv` with `--format json` are usage errors.
/// The value parsing of `--time-limit` / `--iter-limit` stays where it was in `run`.
fn parse_solve_flags(flags: &[String]) -> Result<SolveFmt, String> {
    let mut fmt = SolveFmt {
        csv: false,
        json: false,
        include_solution: false,
    };
    let mut duals = false;
    let mut it = flags.iter();
    while let Some(f) = it.next() {
        match f.as_str() {
            "--csv" => fmt.csv = true,
            "--format" => match it.next().map(String::as_str) {
                Some("json") => fmt.json = true,
                Some("csv") => fmt.csv = true,
                Some("text") => {}
                _ => return Err("--format needs json, csv or text".to_string()),
            },
            "--solution" | "--include-solution" => fmt.include_solution = true,
            "--duals" => duals = true,
            "--no-presolve" | "--ipm" | "--legacy-lp" | "--ekk" | "--no-warnings" => {}
            "--time-limit" | "--iter-limit" | "--solution-file" => {
                if it.next().is_none() {
                    return Err(format!("{f} needs a value"));
                }
            }
            other => {
                return Err(format!(
                    "unknown argument '{other}' (see `rusolve solve --help`)"
                ));
            }
        }
    }
    if fmt.csv && fmt.json {
        return Err(
            "--csv and --format json cannot be combined: --format json replaces the CSV line"
                .to_string(),
        );
    }
    if fmt.json && duals {
        return Err("--duals cannot be combined with --format json".to_string());
    }
    if fmt.include_solution && !fmt.json {
        return Err("--solution needs --format json".to_string());
    }
    Ok(fmt)
}

/// Engine choice for `rusolve solve`: ekk unless the IPM is asked for, the
/// legacy simplex is asked for (`--legacy-lp` / `RUSOLVE_EKK=0`), or a
/// legacy-only `RUSOLVE_*` lever is set — census and diff scripts
/// set those to study the legacy engine's internals, and ekk would ignore
/// them silently. An explicit `--ekk` / `RUSOLVE_EKK=1` overrides the levers.
///
/// "legacy-only lever" is a registry fact, not an allowlist: the variable has
/// `EnvScope::LegacyLp` (only the legacy engine reads it) and `EnvKind::Lever`
/// (it is not test plumbing) in `env_registry`. Anything the port, the MIP or
/// the test harness reads no longer moves the engine (no silent switches).
fn use_ekk_engine(flags: &[String]) -> bool {
    let has = |f: &str| flags.iter().any(|a| a == f);
    if has("--ipm") || std::env::var("RUSOLVE_IPM").is_ok_and(|v| v != "0") {
        return false;
    }
    let ekk_env = std::env::var("RUSOLVE_EKK").ok();
    if has("--legacy-lp") || ekk_env.as_deref() == Some("0") {
        return false;
    }
    if has("--ekk") || ekk_env.is_some() {
        return true;
    }
    if let Some(v) = env_registry::legacy_engine_lever() {
        eprintln!(
            "rusolve: {} is a legacy-engine lever; solving with the legacy simplex",
            v.name
        );
        return false;
    }
    true
}

/// `rusolve env [--json]`: every environment variable the solver reads,
/// with what honours it and what its values mean.
fn run_env(args: &[String]) -> ExitCode {
    for flag in args {
        if flag != "--json" {
            eprintln!("rusolve env: unknown flag '{flag}' (only --json is accepted)");
            return ExitCode::from(2);
        }
    }
    if args.iter().any(|a| a == "--json") {
        println!("{}", env_registry::table_json());
    } else {
        print!("{}", env_registry::table_text());
    }
    ExitCode::SUCCESS
}

/// `rusolve solve --ipm` (or `RUSOLVE_IPM=1`): interior-point LP, no
/// crossover. Presolve uses the same knobs/levers as the simplex path; the
/// output lines keep the simplex format (`iterations:` reports IPM
/// iterations as phase 2) so the harness parsers stay unchanged.
#[allow(clippy::too_many_arguments)]
fn solve_ipm_cli(
    sp: &rusolve_lp::SparseLp,
    model: &rusolve_io::Model,
    path: &str,
    no_presolve: bool,
    csv: bool,
    duals: bool,
    time_limit: Option<f64>,
    t_start: std::time::Instant,
    (json, include_solution): (bool, bool),
) -> ExitCode {
    let obj_offset = model.obj_offset;
    use rusolve_lp::ipm::{IpmConfig, IpmStatus};
    let cfg = IpmConfig {
        deadline: time_limit.map(|t| t_start + std::time::Duration::from_secs_f64(t)),
        max_iter: std::env::var("RUSOLVE_IPM_MAXIT")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(200),
        ..IpmConfig::default()
    };
    let ps = SimplexConfig {
        presolve_implied_eq: std::env::var("RUSOLVE_PS_IMPEQ").is_ok(),
        presolve_parcol: std::env::var("RUSOLVE_PS_NO_PARCOL").is_err(),
        presolve_agg: std::env::var("RUSOLVE_PRESOLVE_AGG").as_deref() != Ok("0"),
        presolve_parcol_dom: std::env::var("RUSOLVE_PS_NO_PARCOLDOM").is_err(),
        presolve_vub_eq: std::env::var("RUSOLVE_PS_NO_VUBEQ").is_err(),
        ..SimplexConfig::default()
    };
    let res = if no_presolve {
        rusolve_lp::ipm::solve_sparse(sp, &cfg)
    } else {
        rusolve_lp::ipm::solve_presolved(sp, &cfg, &ps)
    };
    let sol = match res {
        Ok(s) => s,
        Err(e) => {
            eprintln!("solve: {e}");
            return ExitCode::from(4);
        }
    };
    let reported_obj = if obj_offset != 0.0 {
        sol.objective + obj_offset
    } else {
        sol.objective
    };
    let st = &sol.stats;
    if json {
        let optimal = sol.status == IpmStatus::Optimal;
        let native = envelope::NativeSolved {
            path,
            model,
            status: match sol.status {
                IpmStatus::Optimal => "Optimal",
                IpmStatus::PrimalInfeasible => "Infeasible",
                IpmStatus::DualInfeasible => "Unbounded",
                IpmStatus::IterationLimit => "IterationLimit",
                IpmStatus::TimeLimit => "TimeLimit",
                _ => "Unknown",
            },
            objective: reported_obj,
            dual_bound: if optimal { reported_obj } else { f64::NAN },
            nodes: 0,
            lp_iterations: sol.iterations as i64,
            time_s: t_start.elapsed().as_secs_f64(),
            seed: 0,
            threads: 1,
            time_limit,
            node_limit: None,
            x: if optimal { &sol.x } else { &[] },
            include_solution,
        };
        println!("{}", envelope::native_envelope_json(&native));
    } else if csv {
        println!("status,objective,phase1_iters,phase2_iters,x");
        let xs: Vec<String> = sol.x.iter().map(|v| format!("{v}")).collect();
        println!(
            "{:?},{},{},{},[{}]",
            sol.status,
            reported_obj,
            0,
            sol.iterations,
            xs.join(";")
        );
    } else {
        println!("status:  {:?}", sol.status);
        println!("objective: {}", reported_obj);
        println!("iterations: 0 (phase1) + {} (phase2) [ipm]", sol.iterations);
        let first5: Vec<String> = sol.x.iter().take(5).map(|v| format!("{v:.6}")).collect();
        println!("x[:5]: [{}]", first5.join(", "));
    }
    eprintln!(
        "rusolve ipm stats: m {} n {} nnz(A) {} nnz(M) {} nnz(L) {} supernodes {} dense_cols {} dense_split {} correctors {} skipped_pivots {} order {:.3}s symbolic {:.3}s factor {:.3}s solve {:.3}s ipm {:.3}s retried {}",
        st.m,
        st.n,
        st.nnz_a,
        st.nnz_m,
        st.nnz_l,
        st.supernodes,
        st.dense_cols,
        st.dense_split,
        st.correctors,
        st.skipped_pivots_max,
        st.t_order,
        st.t_symbolic,
        st.t_factor,
        st.t_solve,
        st.t_ipm,
        st.retried
    );
    if sol.status == IpmStatus::TimeLimit
        && let Some(t) = time_limit
    {
        eprintln!("rusolve: time limit {t} s reached");
    }
    if duals {
        if sol.duals_valid {
            println!("row_duals: {:?}", sol.row_duals);
            println!("reduced_costs: {:?}", sol.reduced_costs);
        } else {
            println!("duals: unavailable (presolve-ON path or non-optimal status)");
        }
    }
    ExitCode::SUCCESS
}

/// `rusolve mip <model.mps> [--csv] [--no-presolve] [--node-limit N] [--time-limit S] [--no-cuts]
/// [--tree-cuts] [--no-tree-cuts] [--tree-cut-effort X] [--lp-age-limit N] [--modk] [--no-modk] [--gomory] [--no-gomory]
/// [--no-warm] [--dual-lp] [--root-dual] [--no-root-dual] [--no-redcost] [--submip] [--no-submip] [--root-dive] [--no-root-dive] [--alns] [--no-alns] [--localmip] [--no-localmip] [--conflict] [--no-conflict] [--restart] [--no-restart] [--conflict-cutoff] [--no-conflict-cutoff] [--conflict-score] [--no-conflict-score] [--sb-invest] [--no-sb-invest] [--warm-dse] [--warm-unfit-last] [--restart-min-gap G] [--cutoff V] [--heur-off LIST] [--heur-on LIST] [--fj-restart-effort F]
/// [--branch SPEC] [--threads N] [--seed N]`, LIST =
/// `trivial,fj,simple,randomized,zi,shifting,diving` or `all`; SPEC as in
/// `BranchConfig::apply` (`RUSOLVE_BRANCH` is applied first, for A/B runs
/// through scripts that pass no flags).
/// CSV contract (tools/diff/mip_diff.py): header
/// `status,objective,bound,nodes,lp_iters,time_s,x`, x as `[v0;v1;…]`.
/// Exit codes are the exit-code table of `mip_highs.rs`: 0 the solve ran (the status is in the
/// output), 2 usage, 3 input, 4 internal. `--format json`, `--solution` and
/// `--solution-file <p>.json` print the same `rusolve.result/1` envelope as the port.
fn run_mip(args: &[String]) -> ExitCode {
    let wj = diag::wants_json(args);
    let usage = |msg: &str| diag::fail_usage_text(wj, "mip", msg);
    let Some(path) = args.first() else {
        return usage("missing MPS path");
    };
    let mut config = rusolve_mip::MipConfig::default();
    if let Ok(spec) = std::env::var("RUSOLVE_BRANCH")
        && let Err(e) = config.branching.apply(&spec)
    {
        return usage(&format!("RUSOLVE_BRANCH: {e}"));
    }
    // Gap audit: `audit ...` stderr lines and
    // LP-solve timing; no effect on the search.
    config.audit = std::env::var("RUSOLVE_AUDIT").is_ok_and(|v| v != "0");
    let mut csv = false;
    // The shared surface of `commands::MIP_FLAGS`, as on the port path
    let mut json = false;
    let mut include_solution = false;
    let mut solution_file: Option<String> = None;
    let mut rest = args[1..].iter();
    while let Some(flag) = rest.next() {
        match flag.as_str() {
            "--csv" => csv = true,
            "--no-presolve" => config.presolve = false,
            "--heur-off" | "--heur-on" => {
                let Some(v) = rest.next() else {
                    return usage(&format!("{flag} needs a value"));
                };
                let set = if flag == "--heur-on" {
                    config.heuristics.enable(v)
                } else {
                    config.heuristics.disable(v)
                };
                if let Err(e) = set {
                    return usage(&e.to_string());
                }
            }
            "--fj-restart-effort" => match rest.next().map(|v| v.parse::<f64>()) {
                Some(Ok(f)) if f >= 0.0 => config.heuristics.fj_restart_effort = f,
                _ => {
                    return usage("--fj-restart-effort needs a value >= 0");
                }
            },
            "--no-cuts" => config.cuts.enabled = false,
            "--tree-cuts" => config.cuts.tree = true,
            "--no-tree-cuts" => config.cuts.tree = false,
            "--modk" => config.cuts.modk = true,
            "--tree-age-ticks" => config.cuts.tree_age_ticks = true,
            "--no-tree-age-ticks" => config.cuts.tree_age_ticks = false,
            "--no-modk" => config.cuts.modk = false,
            "--probe-cliques" => config.cuts.probe_cliques = true,
            "--clique-extend" => config.cuts.clique_extend = true,
            "--no-clique-extend" => config.cuts.clique_extend = false,
            "--no-probe-cliques" => config.cuts.probe_cliques = false,
            "--gomory" => config.cuts.gomory = true,
            "--no-gomory" => config.cuts.gomory = false,
            "--drop-all-cuts-at-end" => config.cuts.drop_all_at_end = true,
            "--drop-gomory-at-end" => config.cuts.drop_gomory_at_end = true,
            "--efficacy-rounds" => config.cuts.efficacy_rounds = true,
            "--efficacy-end" => config.cuts.efficacy_end = true,
            "--no-efficacy-end" => config.cuts.efficacy_end = false,
            "--gomory-park" => config.cuts.gomory_park = true,
            "--no-gomory-park" => config.cuts.gomory_park = false,
            "--root-cutoff-row" => config.cuts.root_cutoff_row = true,
            "--no-root-cutoff-row" => config.cuts.root_cutoff_row = false,
            "--root-cutoff-tighten" => config.cuts.root_cutoff_tighten = true,
            "--no-root-cutoff-tighten" => config.cuts.root_cutoff_tighten = false,
            "--gomory-density" => match rest.next().map(|v| v.parse::<f64>()) {
                Some(Ok(f)) if f >= 0.0 => config.cuts.gomory_density = f,
                _ => {
                    return usage("--gomory-density needs a value >= 0");
                }
            },
            "--no-warm" => config.warm_start = false,
            "--dual-lp" => config.lp.dual = true,
            "--no-redcost" => config.redcost = false,
            "--submip" => config.submip = true,
            "--no-submip" => config.submip = false,
            "--root-dive" => {
                // Lever: the whole dive-shaped root chain in one flag.
                config.root_dive = true;
                config.submip = true;
            }
            "--no-root-dive" => config.root_dive = false,
            "--alns" => config.alns = true,
            "--no-alns" => config.alns = false,
            "--localmip" => config.localmip = true,
            "--no-localmip" => config.localmip = false,
            "--conflict" => config.conflict = true,
            "--no-conflict" => config.conflict = false,
            "--restart" => config.restart = true,
            "--no-restart" => config.restart = false,
            "--clique-fixings" => config.clique_fixings = true,
            "--no-clique-fixings" => config.clique_fixings = false,
            "--root-probing" => config.root_probing = true,
            "--no-root-probing" => config.root_probing = false,
            "--conflict-cutoff" => config.conflict_cutoff = true,
            "--conflict-score" => config.conflict_score = true,
            "--no-conflict-score" => config.conflict_score = false,
            "--sb-invest" => config.branching.sb_invest = true,
            "--no-sb-invest" => config.branching.sb_invest = false,
            "--root-dual" => config.root_lp_dual = true,
            "--no-root-dual" => config.root_lp_dual = false,
            "--warm-dse" => config.warm_dse = true,
            "--warm-unfit-last" => config.warm_unfit_last = true,
            "--no-warm-dse" => config.warm_dse = false,
            "--no-conflict-cutoff" => config.conflict_cutoff = false,
            "--cutoff" => match rest.next().map(|v| v.parse::<f64>()) {
                // Audit lever: known-objective cutoff (HiGHS
                // `objective_bound`), caller sense.
                Some(Ok(c)) if !c.is_nan() => config.objective_cutoff = Some(c),
                _ => {
                    return usage("--cutoff needs a number");
                }
            },
            "--branch" => {
                let Some(v) = rest.next() else {
                    return usage("--branch needs a value");
                };
                if let Err(e) = config.branching.apply(v) {
                    return usage(&e.to_string());
                }
            }
            "--threads" | "--seed" => {
                let Some(v) = rest.next() else {
                    return usage(&format!("{flag} needs a value"));
                };
                let ok = if flag == "--threads" {
                    v.parse()
                        .ok()
                        .filter(|&n: &usize| n >= 1)
                        .map(|n| config.threads = n)
                        .is_some()
                } else {
                    v.parse().map(|n| config.seed = n).is_ok()
                };
                if !ok {
                    return usage(&format!("bad value for {flag}: {v}"));
                }
            }
            "--node-limit" | "--time-limit" | "--tree-cut-effort" | "--lp-age-limit"
            | "--weak-rounds" | "--flat-weak-rounds" | "--max-rounds" | "--restart-min-gap" => {
                let Some(v) = rest.next() else {
                    return usage(&format!("{flag} needs a value"));
                };
                let ok = if flag == "--node-limit" {
                    v.parse().map(|n| config.node_limit = Some(n)).is_ok()
                } else if flag == "--tree-cut-effort" {
                    v.parse().map(|e| config.cuts.tree_effort = e).is_ok()
                } else if flag == "--lp-age-limit" {
                    // 0 = drop a row at the first aging check that finds it
                    // slack; HiGHS bounds mip_lp_age_limit at [0, int16::MAX]
                    // with default 10 (HighsOptions.h:1175–1181).
                    v.parse().map(|n| config.cuts.lp_age_limit = n).is_ok()
                } else if flag == "--weak-rounds" {
                    // root-loop early stop: N weak rounds in a row end it
                    // (cuts/mod.rs CutConfig::weak_rounds, default 4).
                    v.parse().map(|n| config.cuts.weak_rounds = n).is_ok()
                } else if flag == "--flat-weak-rounds" {
                    // weak-round budget before the root bound first moves
                    // (CutConfig::flat_weak_rounds).
                    v.parse().map(|n| config.cuts.flat_weak_rounds = n).is_ok()
                } else if flag == "--restart-min-gap" {
                    // Restart root-gap bar (restart::MIN_ROOT_GAP).
                    v.parse().map(|g| config.restart_min_gap = g).is_ok()
                } else if flag == "--max-rounds" {
                    // root-loop round cap; default = HiGHS 2·√maxTreeSizeLog2
                    v.parse().map(|n| config.cuts.max_rounds = Some(n)).is_ok()
                } else {
                    v.parse().map(|t| config.time_limit = Some(t)).is_ok()
                };
                if !ok {
                    return usage(&format!("bad value for {flag}: {v}"));
                }
            }
            "--format" => match rest.next().map(String::as_str) {
                Some("json") => json = true,
                Some("csv") => csv = true,
                Some("text") => {}
                _ => {
                    return usage("--format needs json, csv or text");
                }
            },
            // A no-op on the native stack (it has no lint pass), accepted so a command line
            // moves between stacks unchanged
            "--no-warnings" => {}
            "--solution" | "--include-solution" => include_solution = true,
            "--solution-file" => match rest.next() {
                Some(p) => solution_file = Some(p.clone()),
                None => {
                    return usage("--solution-file needs a path");
                }
            },
            other => {
                return match commands::flag_stack_error("mip", other, commands::Stack::Native) {
                    Some(msg) => diag::fail_msg(wj, "RS0002", 2, msg),
                    None => diag::fail_msg(
                        wj,
                        "RS0002",
                        2,
                        format!("rusolve mip: unknown flag '{other}' (see `rusolve mip --help`)"),
                    ),
                };
            }
        }
    }
    if csv && json {
        return usage(
            "--csv and --format json cannot be combined: --format json replaces the CSV line",
        );
    }
    if let Some(p) = &solution_file
        && !p.ends_with(".json")
    {
        return usage(
            "--solution-file on the native search writes the json solution object, so PATH must end in .json (the HiGHS raw format is port-only: unset RUSOLVE_MIP_HIGHS)",
        );
    }
    let model = match rusolve_io::read_model_file(path) {
        Ok(m) => m,
        Err(e) => {
            return diag::fail_read(wj, path, &e.to_string());
        }
    };
    // The native engine has no semi-continuous or general-constraint
    // support; the port path (`rusolve mip` default) has.
    if let Some(what) = model.gc.first_general() {
        return diag::fail_msg(
            wj,
            "RS0012",
            3,
            format!(
                "rusolve mip (native engine): {what} is not supported here; use the default engine"
            ),
        );
    }
    if !model.gc.semi.is_empty() {
        return diag::fail_msg(
            wj,
            "RS0012",
            3,
            "rusolve mip (native engine): semi-continuous columns are not supported here; use the default engine",
        );
    }
    if !model.hessian.is_empty() {
        return diag::fail_msg(
            wj,
            "RS0012",
            3,
            format!("rusolve mip: {}", qp::unsupported("mip (native search)")),
        );
    }
    // The MPS objective constant
    // (`Model::obj_offset`) is reporting-only, as in `solve` above and in
    // HiGHS (`lp.offset_` added to the reported objective,
    // HighsInterface.cpp:3195): every objective-space value printed below
    // goes through `reported`, and `--cutoff` (given in that same reported
    // space, e.g. a known optimum) is shifted back into raw `c·x` space for
    // the solver. Only when a constant is present, so offset-free instances
    // stay bit-identical (no `-0.0 + 0.0` sign flips).
    let offset = model.obj_offset;
    let reported = |v: f64| if offset != 0.0 { v + offset } else { v };
    if offset != 0.0
        && let Some(c) = config.objective_cutoff.as_mut()
    {
        *c -= offset;
    }
    if config.audit {
        // The `audit` lines of the solver are raw `c·x` values; this lets
        // the harness lift them to the reported objective.
        eprintln!("audit model obj_offset={offset:e}");
    }
    let problem = rusolve_mip::MipProblem::from_model(&model);
    let start = std::time::Instant::now();
    let result = rusolve_mip::solve(&problem, &config);
    let secs = start.elapsed().as_secs_f64();
    let sol = match result {
        Ok(s) => s,
        Err(e) => {
            if csv {
                println!("status,objective,bound,nodes,lp_iters,time_s,x");
                println!("Error,nan,nan,0,0,{secs},[]");
            }
            return diag::fail_msg(wj, "RS0090", 4, format!("mip: {e}"));
        }
    };
    if config.audit {
        let l = &sol.lp_stats;
        eprintln!(
            "audit end status={:?} objective={:e} bound={:e} nodes={} lp_iters={} \
             lp_s={:.6} probe_s={:.6} probe_iters={} restarts={} root_retries={} time={secs:.6}",
            sol.status,
            reported(sol.objective),
            reported(sol.bound),
            sol.nodes,
            sol.lp_iterations,
            l.lp_nanos as f64 * 1e-9,
            l.probe_nanos as f64 * 1e-9,
            l.probe_iterations,
            sol.restarts,
            l.root_retries
        );
    }
    // The envelope the port prints, from the native result. `--solution-file <p>.json`
    // writes its `solution` object (the raw HiGHS file format is port-only).
    let native = envelope::NativeSolved {
        path,
        model: &model,
        status: match sol.status {
            rusolve_mip::MipStatus::Optimal => "Optimal",
            rusolve_mip::MipStatus::Infeasible => "Infeasible",
            rusolve_mip::MipStatus::Unbounded => "Unbounded",
            rusolve_mip::MipStatus::NodeLimit => "NodeLimit",
            rusolve_mip::MipStatus::TimeLimit => "TimeLimit",
        },
        objective: reported(sol.objective),
        dual_bound: reported(sol.bound),
        nodes: sol.nodes as i64,
        lp_iterations: sol.lp_iterations as i64,
        time_s: secs,
        seed: config.seed as i64,
        threads: config.threads as i64,
        time_limit: config.time_limit,
        node_limit: config.node_limit.map(|n| n as i64),
        x: &sol.x,
        include_solution,
    };
    if let Some(sol_path) = &solution_file
        && let Err(e) = std::fs::write(sol_path, envelope::native_solution_file_json(&native))
    {
        return diag::fail_msg(
            wj,
            "RS0090",
            4,
            format!("rusolve mip: cannot write {sol_path}: {e}"),
        );
    }
    if json {
        println!("{}", envelope::native_envelope_json(&native));
    } else if csv {
        println!("status,objective,bound,nodes,lp_iters,time_s,x");
        let xs: Vec<String> = sol.x.iter().map(|v| format!("{v}")).collect();
        println!(
            "{:?},{},{},{},{},{secs},[{}]",
            sol.status,
            reported(sol.objective),
            reported(sol.bound),
            sol.nodes,
            sol.lp_iterations,
            xs.join(";")
        );
    } else {
        println!("status:    {:?}", sol.status);
        println!("objective: {}", reported(sol.objective));
        println!("bound:     {}", reported(sol.bound));
        println!(
            "nodes:     {} ({} LP iterations, {secs:.3}s)",
            sol.nodes, sol.lp_iterations
        );
        if sol.restarts > 0 {
            println!("restarts:  {}", sol.restarts);
        }
        let pb = &sol.probing;
        if pb.probed > 0 {
            println!(
                "probing:   {} of {} binaries probed, {} infeasible directions, {} fixings, \
                 {} tightenings, {} 2-cliques; {:.3}s",
                pb.probed,
                pb.candidates,
                pb.infeasible_probes,
                pb.fixings,
                pb.tightenings,
                pb.cliques,
                pb.seconds
            );
        }
        let l = &sol.lp_stats;
        println!(
            "LPs:       {} warm ({} cutoff, {} infeasible re-checked, {} overturned, {} ray-trusted), {} cold, \
             {} fallbacks ({} unfit, {} rejected, {} iteration limit, {} error); \
             iterations {} warm / {} cold",
            l.warm,
            l.cutoffs,
            l.infeasible_checks,
            l.infeasible_overturned,
            l.infeasible_trusted,
            l.cold,
            l.fallbacks(),
            l.fallback_unfit,
            l.fallback_rejected,
            l.fallback_iter_limit,
            l.fallback_error,
            l.warm_iterations,
            l.cold_iterations
        );
        if l.unfit_recovered > 0 {
            println!(
                "           {} unfit snapshots warm-started from the latest basis",
                l.unfit_recovered
            );
        }
        if l.no_lp_branches + l.unsolved_pruned > 0 {
            println!(
                "unsolved:  {} node LPs branched on without an LP, {} fixed nodes pruned unsolved",
                l.no_lp_branches, l.unsolved_pruned
            );
        }
        println!(
            "strong branching: {} probes ({} at the iteration cap), {} LP iterations ({:.1}% of all)",
            l.probes,
            l.probe_iter_limits,
            l.probe_iterations,
            100.0 * l.probe_iterations as f64 / sol.lp_iterations.max(1) as f64
        );
        let r = &sol.redcost;
        println!(
            "redcost:   {} local tightenings on {} of {} nodes with duals, {} root (global) \
             tightenings, {} LP re-solves, {} nodes pruned; duals {} from LP / {} from basis / \
             {} unavailable",
            r.local_tightenings,
            r.nodes_tightened,
            r.nodes_with_duals,
            r.root_tightenings,
            r.resolves,
            r.pruned,
            l.duals_lp,
            l.duals_basis,
            l.duals_failed
        );
        let c = &sol.conflict;
        println!(
            "conflict:  {} proofs ({} from cutoffs) -> {} conflicts stored ({} duplicate or too long, {} without \
             usable proof), {} aged out; {} bound tightenings, {} nodes pruned by propagation",
            c.proofs,
            c.cutoff_proofs,
            c.added,
            c.rejected,
            c.no_proof,
            c.aged_out,
            c.tightenings,
            c.pruned
        );
        let s = &sol.submip;
        println!(
            "submip:    RENS {} calls / {} improving, RINS {} calls / {} improving, {} skipped; \
             {} nodes, {} LP iterations, {:.3}s",
            s.rens_calls,
            s.rens_found,
            s.rins_calls,
            s.rins_found,
            s.skipped,
            s.nodes,
            s.lp_iterations,
            s.seconds
        );
        if s.alns_calls.iter().any(|&c| c > 0) {
            let arms: Vec<String> = rusolve_mip::alns::Arm::ALL
                .iter()
                .enumerate()
                .map(|(k, a)| format!("{} {}/{}", a.name(), s.alns_found[k], s.alns_calls[k]))
                .collect();
            println!(
                "alns:      {} (improving/calls), {:.3}s",
                arms.join(", "),
                s.alns_seconds
            );
        }
        if let Some(f) = sol.first_incumbent {
            println!(
                "first incumbent: {:?} after {} nodes, {:.4}s, value {}",
                f.source,
                f.nodes,
                f.seconds,
                reported(f.value)
            );
        }
        for f in sol.incumbents.iter().skip(1) {
            println!(
                "  improving: {} after {} nodes, {:.4}s {:?}",
                reported(f.value),
                f.nodes,
                f.seconds,
                f.source
            );
        }
        let c = &sol.cuts;
        if c.rounds > 0 {
            println!(
                "cuts:      {} rounds, root bound {} -> {} ({} LP iterations, {:.3}s){}",
                c.rounds,
                reported(c.root_bound_before),
                reported(c.root_bound_after),
                c.lp_iterations,
                c.time_s,
                if c.cutoff_rows > 0 {
                    " +cutoff row"
                } else {
                    ""
                }
            );
            for class in rusolve_mip::CutClass::ALL {
                let k = class as usize;
                if c.generated[k] > 0 {
                    println!(
                        "  {:<16} generated {:>6}  applied {:>6}  in LP {:>6}",
                        class.name(),
                        c.generated[k],
                        c.applied[k],
                        c.in_lp[k]
                    );
                }
            }
            if c.aggregated > 0 {
                println!("  ({} generated from multi-row aggregations)", c.aggregated);
            }
            if c.probe_probed > 0 {
                println!(
                    "  probing: {} binaries probed, {} cliques, {} fixings, {} implied-bound cuts ({:.3}s)",
                    c.probe_probed, c.probe_cliques, c.probe_fixings, c.probe_ibcuts, c.probe_s
                );
            }
            if c.tree_rounds > 0 {
                println!(
                    "  tree: {} calls, {} rounds, {} cuts applied ({} LP iterations)",
                    c.tree_calls, c.tree_rounds, c.tree_applied, c.tree_lp_iterations
                );
            }
            if c.lp_failures > 0 {
                println!(
                    "  {} cut-loop LP re-solves failed (loops stopped on the last good LP)",
                    c.lp_failures
                );
            }
        }
    }
    ExitCode::SUCCESS
}

/// `rusolve classify <model.mps> [--json] [--miplib-strict]`: MIPLIB 2017
/// row-type histogram, column / objective statistics, structure detectors
/// and the feature vector (rusolve_mip::classify). No solve.
fn run_classify(args: &[String]) -> ExitCode {
    let wj = diag::wants_json(args);
    let usage = |msg: &str| diag::fail_usage_text(wj, "classify", msg);
    let Some(path) = args.first() else {
        return usage("missing model path");
    };
    let mut json = false;
    let mut format_json = false; // `--format json` -> rusolve.classify/1
    let mut opts = rusolve_mip::classify::ClassifyOptions::default();
    let mut flags = args[1..].iter();
    while let Some(flag) = flags.next() {
        match flag.as_str() {
            "--json" => json = true,
            "--miplib-strict" => opts.miplib_strict = true,
            "--format" => match flags.next().map(String::as_str) {
                Some("json") => format_json = true,
                Some("text") => format_json = false,
                _ => return usage("--format needs json or text"),
            },
            other if !other.starts_with("--") => {
                return usage(&format!("unexpected extra argument '{other}'"));
            }
            other => return usage(&format!("unknown flag '{other}'")),
        }
    }
    if format_json {
        return features::classify_json(path, opts.miplib_strict);
    }
    let model = match rusolve_io::read_model_file(path) {
        Ok(m) => m,
        Err(e) => return diag::fail_read(wj, path, &e.to_string()),
    };
    if let Some(what) = model.gc.first_general() {
        return diag::fail_msg(
            wj,
            "RS0012",
            3,
            format!("classify: {what} is not supported"),
        );
    }
    let problem = rusolve_mip::MipProblem::from_model(&model);
    let c = rusolve_mip::classify::classify(&problem, opts);
    if json {
        println!("{}", c.to_json());
    } else {
        print!("{}", c.to_text());
    }
    ExitCode::SUCCESS
}

pub fn main_entry() -> ExitCode {
    // A `RUSOLVE_*` name no crate reads is a typo, and a typo must not
    // change which engine runs. Fail before any subcommand, with the usage
    // code of the exit-code table in `mip_highs`.
    let args = std::env::args().skip(1).collect::<Vec<String>>();
    if let Some((msg, hint)) = env_registry::unknown_env() {
        let mut d = rusolve_core::diagnostic::Diagnostic::error("RS0005", msg);
        if let Some(h) = hint {
            d = d.with_help(h);
        }
        return diag::fail(diag::wants_json(&args), d, 2);
    }
    run(&args)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_lp::SimplexConfig;

    fn solve_file(name: &str, cfg: SimplexConfig) -> rusolve_lp::SimplexSolution {
        SimplexSolver::new(cfg).solve(&lp_file(name)).unwrap()
    }

    fn lp_file(name: &str) -> Lp {
        let path = format!("../../benches/data/netlib/{name}.mps");
        let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let m = read_mps(&text).unwrap_or_else(|e| panic!("{name}: {e}"));
        model_to_lp(&m).unwrap()
    }

    fn solve_repro(name: &str, cfg: SimplexConfig) -> rusolve_lp::SimplexSolution {
        let path = format!("../../benches/data/repros/{name}.mps");
        let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let m = read_mps(&text).unwrap_or_else(|e| panic!("{name}: {e}"));
        let lp = model_to_lp(&m).unwrap();
        SimplexSolver::new(cfg).solve(&lp).unwrap()
    }

    #[test]
    fn a2_bounded_bland_episode_exits_stall_on_packing_lps() {
        // The three set_packing LPs hit the 100k-iteration cap in phase 2
        // (or ground 10,021 iters on the root) because the engaged-Bland state
        // was terminal; with Bland disabled entirely they solve in 174–313 p2
        // iters. Bounded Bland episodes must reach the correct optimum in
        // <1,000 phase-2 iterations each, both with and without
        // phase-2 perturbation (the perturbation must not be the only way
        // out — the release itself is the fix).
        for pert in [true, false] {
            for (name, want_obj) in [
                ("set_packing_s1_node_stall", 746.0),
                ("set_packing_s2_node_stall", 660.0),
                ("set_packing_s1_root_slow", 758.3333333333333),
            ] {
                let s = solve_repro(
                    name,
                    SimplexConfig {
                        presolve: false,
                        p2_tie_perturb: pert,
                        ..Default::default()
                    },
                );
                assert_eq!(
                    s.status,
                    rusolve_lp::SimplexStatus::Optimal,
                    "{name} pert={pert}"
                );
                assert!(
                    (s.objective - want_obj).abs() <= 1e-9 * want_obj.abs() + 1e-6,
                    "{name} pert={pert} obj {} want {want_obj}",
                    s.objective
                );
                assert!(
                    s.bland_triggered,
                    "{name} pert={pert}: repro must actually engage Bland (else it does not test the Bland budget)"
                );
                assert!(
                    s.phase2_iterations < 1000,
                    "{name} pert={pert} p2={} iters (Bland never exits)",
                    s.phase2_iterations
                );
            }
        }
    }

    #[test]
    fn a2_sticky_bland_lever_reproduces_the_terminal_stall() {
        // The legacy semantics (bland_episode=0 ⇒ Bland never releases) must
        // still hit the iteration limit on the node repros WITH the phase-2
        // perturbation off — pins the diagnosis (the bug is in EXITING Bland,
        // not pricing) and the lever's A/B value.
        let s = solve_repro(
            "set_packing_s1_node_stall",
            SimplexConfig {
                presolve: false,
                p2_tie_perturb: false,
                bland_episode: 0,
                ..Default::default()
            },
        );
        assert_eq!(s.status, rusolve_lp::SimplexStatus::IterationLimit);
    }

    #[test]
    fn a2_episode_release_preserves_optima_and_never_worsens() {
        // Bounded episodes must not change any optimum on the Netlib set, and
        // the default config must not be worse than the sticky-legacy lever in
        // total iterations.
        let mut def_total = 0usize;
        let mut sticky_total = 0usize;
        for name in [
            "adlittle", "afiro", "blend", "boeing2", "bore3d", "brandy", "capri", "e226", "finnis",
            "grow7", "kb2", "lotfi", "sc105", "sc50a", "scagr25", "scagr7", "scorpion", "scrs8",
            "scsd1", "scsd6", "sctap1", "share1b", "share2b", "stair", "standata", "standmps",
            "tuff", "vtp.base",
        ] {
            let def = solve_file(name, SimplexConfig::default());
            let sticky = solve_file(
                name,
                SimplexConfig {
                    bland_episode: 0,
                    ..Default::default()
                },
            );
            assert_eq!(def.status, sticky.status, "{name} status");
            if def.objective.is_finite() {
                assert!(
                    (def.objective - sticky.objective).abs()
                        <= 1e-9 * def.objective.abs().max(1.0) + 1e-6,
                    "{name} obj def={} sticky={}",
                    def.objective,
                    sticky.objective
                );
            }
            def_total += def.phase1_iterations + def.phase2_iterations;
            sticky_total += sticky.phase1_iterations + sticky.phase2_iterations;
        }
        assert!(
            def_total <= sticky_total,
            "episodes worsened Netlib totals: def={def_total} sticky={sticky_total}"
        );
    }

    #[test]
    fn sparse_lever_bit_identical_on_netlib_no_presolve() {
        // Corpus pin (a representative subset with presolve off — the only
        // mode the CSC-only path serves): the
        // sparse front door must reproduce the dense solve bit for bit —
        // identical per-phase iteration counts, objective bits, x bits.
        for name in [
            "afiro", "blend", "boeing2", "e226", "scsd6", "sctap1", "stair", "tuff",
        ] {
            let lp = lp_file(name);
            let off = SimplexSolver::new(SimplexConfig {
                presolve: false,
                ..Default::default()
            })
            .solve(&lp)
            .unwrap();
            let sp = rusolve_lp::SparseLp::from_lp(&lp);
            let on = SimplexSolver::new(SimplexConfig {
                presolve: false,
                sparse: true,
                ..Default::default()
            })
            .solve_sparse(&sp)
            .unwrap();
            assert_eq!(on.status, off.status, "{name}");
            assert_eq!(
                on.objective.to_bits(),
                off.objective.to_bits(),
                "{name} obj"
            );
            assert_eq!(on.phase1_iterations, off.phase1_iterations, "{name} p1");
            assert_eq!(on.phase2_iterations, off.phase2_iterations, "{name} p2");
            for (i, (a, b)) in on.x.iter().zip(&off.x).enumerate() {
                assert_eq!(a.to_bits(), b.to_bits(), "{name} x[{i}]");
            }
        }
    }

    #[test]
    fn p2_perturbation_reduces_stalled_phase2_iters() {
        // On the phase-2 degenerate grinders the perturbed run terminates
        // with FEWER phase-2 iterations and the same objective (measured:
        // scsd6 p2 946→277, tuff p2 308→278).
        // Pinned to bland_episode=0 (legacy sticky Bland): the bounded
        // Bland episodes give the same runs a second escape (release to fast
        // pricing), which would erase the
        // perturbation-only delta this test measures.
        // Aggregator OFF: the presolve aggregator removes the
        // tuff phase-2 stall itself (p2 104 both ways), which would erase
        // the delta. Presolve stays ON: primal on raw tuff with this
        // config misreports Optimal at obj 0 (bounds violated 1.7e2), a
        // separate open defect.
        // DupCol merge OFF (same isolation convention as the aggregator
        // pin): the duplicate-column merge being on by default
        // changes tuff's reduced model and erases the p2 delta (198=198);
        // this test measures the PERTURBATION, so it freezes the presolve
        // state the numbers above were taken under.
        for (name, want) in [("scsd6", 50.5), ("tuff", 0.2921477650936129)] {
            let off = solve_file(
                name,
                SimplexConfig {
                    p2_tie_perturb: false,
                    bland_episode: 0,
                    presolve_agg: false,
                    presolve_parcol: false,
                    ..Default::default()
                },
            );
            let on = solve_file(
                name,
                SimplexConfig {
                    p2_tie_perturb: true,
                    bland_episode: 0,
                    presolve_agg: false,
                    presolve_parcol: false,
                    ..Default::default()
                },
            );
            assert_eq!(on.status, off.status, "{name}");
            assert!(
                (on.objective - off.objective).abs() <= 1e-9 * off.objective.abs().max(1.0),
                "{name} obj on={} off={}",
                on.objective,
                off.objective
            );
            assert!(
                (on.objective - want).abs() <= 1e-9 * want.abs().max(1.0) + 1e-6,
                "{name} obj {} want {want}",
                on.objective
            );
            assert!(
                on.phase2_iterations < off.phase2_iterations,
                "{name} p2 on={} off={}",
                on.phase2_iterations,
                off.phase2_iterations
            );
            assert!(
                on.phase1_iterations == off.phase1_iterations,
                "{name} p1 moved"
            );
        }
    }

    #[test]
    fn p2_perturbation_bit_identical_when_trigger_never_fires() {
        // Instances whose phase 2 never reaches the
        // escalation rung never build eps ⇒ lever on/off is bit-identical
        // (objective, primal vector, per-phase iteration counts).
        for name in ["afiro", "blend", "kb2"] {
            let off = solve_file(
                name,
                SimplexConfig {
                    p2_tie_perturb: false,
                    ..Default::default()
                },
            );
            let on = solve_file(
                name,
                SimplexConfig {
                    p2_tie_perturb: true,
                    ..Default::default()
                },
            );
            assert_eq!(on.status, off.status, "{name}");
            assert_eq!(
                on.objective.to_bits(),
                off.objective.to_bits(),
                "{name} obj"
            );
            assert_eq!(on.phase1_iterations, off.phase1_iterations, "{name} p1");
            assert_eq!(on.phase2_iterations, off.phase2_iterations, "{name} p2");
            for (i, (a, b)) in on.x.iter().zip(&off.x).enumerate() {
                assert_eq!(a.to_bits(), b.to_bits(), "{name} x[{i}]");
            }
        }
    }
}
