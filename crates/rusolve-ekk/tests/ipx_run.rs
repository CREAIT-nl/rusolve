//! Golden test for IPX end to end: `Highs::passModel` +
//! `Highs::run` with `solver=ipx` through `solveLp` -> `solveLpIpx`
//! (ipm/ipx_wrapper.rs) -> `ipx::LpSolver::Solve` (IPM driver, crossover)
//! -> the HiGHS solution and basis copies (lp_data/highs_solution_ipx.rs).
//!
//! Fixtures: `ipx/run/<lp>.jsonl`, written by `highs-oracle ipx`
//! (tools/highs-oracle/src/ipx.cpp runIpx) for every netlib LP and the
//! synthetic `tools/highs-oracle/lps/*.mps`. The header holds the LP HiGHS
//! read and the `useHipo` probe; then, per configuration (`xon`
//! run_crossover=on, `xoff` run_crossover=off, `ipm` solver=ipm, `dual`
//! ipx_dualize_strategy=1), a `run_header`, the oracle (b) per-iteration
//! records of `highs-oracle-ipx-hooks.patch` and a `run_result` (statuses,
//! HighsInfo, HighsSolution, HighsBasis).
//!
//! Per configuration:
//! 1. (b) The Rust trace writer (`ipm/ipx/oracle_trace.rs`) captures the
//!    same records; they must be byte-identical (every `mu`, step size,
//!    residual, objective, KKT/basis counter and iterate hash after every
//!    IPM iteration, `driver_end` and `solve_end` `ipx::Info`).
//! 2. (c) The `run_result` is reproduced bit for bit.
//!
//! The `ipx/trace/<lp>.<cfg>.jsonl` fixtures (`highs-oracle ipx-trace`)
//! are byte-diffed the same way by `ipx_trace_matches_highs`, on the LP of
//! the `ipx/run` fixture of the same name.
//!
//! Committed: afiro, ipxmix, ipxtall, ipxdep; the rest is target-only
//! (`RUSOLVE_REQUIRE_TARGET_FIXTURES=1`). The `dual` configuration trips
//! the `debug_assert!` of model/part1.rs (C++ assert compiled out in the
//! NDEBUG oracle) and replays in release builds only.

mod common;
mod port;

use common::records::LpRec;
use common::*;
use port::*;
use rusolve_ekk::ipm::ipx::oracle_trace;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;

/// The HiGHS options of `oracleBaseOptions` + `solver=ipx` + `overrides`.
fn options_for(overrides: &str) -> HighsOptions {
    let mut o = HighsOptions {
        output_flag: false,
        threads: 1,
        parallel: "off".to_string(),
        presolve: "off".to_string(),
        random_seed: 0,
        solver: "ipx".to_string(),
        ..port::profile_options()
    };
    for kv in overrides.split(',') {
        if kv == "default" {
            continue;
        }
        let (k, v) = kv.split_once('=').expect("key=value override");
        match k {
            "run_crossover" => o.run_crossover = v.to_string(),
            "solver" => o.solver = v.to_string(),
            "ipx_dualize_strategy" => o.ipx_dualize_strategy = v.parse().unwrap(),
            other => panic!("unknown override {other}"),
        }
    }
    o.sync_log_options();
    o
}

/// One configuration of a fixture: its header, hook lines and result.
struct Run<'a> {
    cfg: String,
    overrides: String,
    hook: Vec<&'a str>,
    result: Rec<'a>,
}

fn runs<'a>(fx: &'a Fixture, text: &'a str) -> Vec<Run<'a>> {
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), fx.len(), "{}", fx.file);
    let mut out = Vec::new();
    let mut i = 1;
    while i < fx.len() {
        let h = fx.rec(i);
        assert_eq!(h.t(), "run_header", "{}", h.ctx);
        let mut j = i + 1;
        while fx.rec(j).t() != "run_result" {
            j += 1;
        }
        out.push(Run {
            cfg: h.str("cfg").to_string(),
            overrides: h.str("overrides").to_string(),
            hook: lines[i + 1..j].to_vec(),
            result: fx.rec(j),
        });
        i = j + 1;
    }
    out
}

/// (b): the first differing hook line, field by field.
fn diff_trace(chk: &mut Checker, want: &[&str], got: &str) {
    let got: Vec<&str> = got.lines().collect();
    for (k, (w, g)) in want.iter().zip(got.iter()).enumerate() {
        if w != g {
            let wv: serde_json::Value = serde_json::from_str(w).unwrap();
            let gv: serde_json::Value = serde_json::from_str(g).unwrap();
            for m in diff_values(&format!("trace[{k}]"), &wv, &gv) {
                chk.push(&m.path, m.expected, m.actual, m.note);
            }
            if chk.ok() {
                chk.push(
                    &format!("trace[{k}]"),
                    w.to_string(),
                    g.to_string(),
                    String::new(),
                );
            }
            return;
        }
    }
    chk.int("trace.len", want.len() as i32, got.len() as i32);
}

fn replay(fails: &mut Failures, fx: &Fixture) {
    let text = std::fs::read_to_string(&fx.path).unwrap();
    let header = fx.first("ipx_header");
    assert_eq!(
        header.int("hipo_available"),
        0,
        "{}: HiPO in the oracle",
        fx.file
    );
    let lp = lp_from_rec(&LpRec::from_rec(&header.obj("lp")));
    for run in runs(fx, &text) {
        if cfg!(debug_assertions) && run.overrides.contains("ipx_dualize_strategy=1") {
            eprintln!("SKIP {} {} (dualize=1 in a debug build)", fx.file, run.cfg);
            continue;
        }
        let ctx = format!("{} [{}]", fx.file, run.cfg);
        let lp = lp.clone();
        fails.run(&ctx.clone(), || {
            let mut chk = Checker::new(ctx);
            let mut highs = Highs::with_options(options_for(&run.overrides));
            chk.int("pass_model", HighsStatus::kOk as i32, highs.pass_model(lp) as i32);
            oracle_trace::capture_begin();
            let run_status = highs.run();
            let trace = oracle_trace::capture_end();
            diff_trace(&mut chk, &run.hook, &trace);
            let r = &run.result;
            chk.int("run_status", r.int("run_status"), run_status as i32)
                .int(
                    "model_status",
                    r.int("model_status"),
                    highs.get_model_status() as i32,
                );
            let w = r.obj("info");
            let a = highs.get_info();
            macro_rules! ii {
                ($($f:ident)*) => { $( chk.int(concat!("info.", stringify!($f)), w.int(stringify!($f)), a.$f as i32); )* };
            }
            macro_rules! dd {
                ($($f:ident)*) => { $( chk.f64(concat!("info.", stringify!($f)), w.f64(stringify!($f)), a.$f); )* };
            }
            ii!(valid simplex_iteration_count ipm_iteration_count crossover_iteration_count
                pdlp_iteration_count qp_iteration_count primal_solution_status
                dual_solution_status basis_validity);
            dd!(objective_function_value);
            ii!(num_primal_infeasibilities);
            dd!(max_primal_infeasibility sum_primal_infeasibilities);
            ii!(num_dual_infeasibilities);
            dd!(max_dual_infeasibility sum_dual_infeasibilities);
            ii!(num_relative_primal_infeasibilities);
            dd!(max_relative_primal_infeasibility);
            ii!(num_relative_dual_infeasibilities);
            dd!(max_relative_dual_infeasibility);
            ii!(num_primal_residual_errors);
            dd!(max_primal_residual_error);
            ii!(num_dual_residual_errors);
            dd!(max_dual_residual_error);
            ii!(num_relative_primal_residual_errors);
            dd!(max_relative_primal_residual_error);
            ii!(num_relative_dual_residual_errors);
            dd!(max_relative_dual_residual_error);
            ii!(num_complementarity_violations);
            dd!(max_complementarity_violation primal_dual_objective_error);
            let sol = highs.get_solution();
            let basis = highs.get_basis();
            let col_status: Vec<i32> = basis.col_status.iter().map(|s| *s as i32).collect();
            let row_status: Vec<i32> = basis.row_status.iter().map(|s| *s as i32).collect();
            chk.bool("value_valid", r.bool("value_valid"), sol.value_valid)
                .bool("dual_valid", r.bool("dual_valid"), sol.dual_valid)
                .f64s("col_value", &r.f64s("col_value"), &sol.col_value)
                .f64s("col_dual", &r.f64s("col_dual"), &sol.col_dual)
                .f64s("row_value", &r.f64s("row_value"), &sol.row_value)
                .f64s("row_dual", &r.f64s("row_dual"), &sol.row_dual)
                .bool("basis_valid", r.bool("basis_valid"), basis.valid)
                .bool("basis_alien", r.bool("basis_alien"), basis.alien)
                .ints("col_status", &r.ints("col_status"), &col_status)
                .ints("row_status", &r.ints("row_status"), &row_status);
            Some(chk)
        });
        // A panic inside the run leaves the capture installed.
        oracle_trace::capture_end();
    }
}

#[test]
fn ipx_run_matches_highs() {
    let mut fails = Failures::default();
    let mut files = 0;
    for fx in fixtures_in("ipx/run") {
        replay(&mut fails, &fx);
        files += 1;
    }
    assert!(files > 0, "no ipx/run fixtures");
    fails.finish("ipx_run");
}

/// Oracle (b) on the `ipx/trace` fixtures: hook lines byte-identical and
/// the `result` record bit for bit.
#[test]
fn ipx_trace_matches_highs() {
    let lps: Vec<(String, LpRec)> = fixtures_in("ipx/run")
        .iter()
        .map(|fx| {
            let h = fx.first("ipx_header");
            (h.str("name").to_string(), LpRec::from_rec(&h.obj("lp")))
        })
        .collect();
    let mut fails = Failures::default();
    let mut traces = 0;
    for fx in fixtures_in("ipx/trace") {
        let header = fx.rec(0);
        let (name, cfg) = (header.str("name"), header.str("cfg"));
        if cfg!(debug_assertions) && cfg.contains("ipx_dualize_strategy=1") {
            eprintln!("SKIP {} (dualize=1 in a debug build)", fx.file);
            continue;
        }
        let lp = match lps.iter().find(|(n, _)| n == name) {
            Some((_, lp)) => lp_from_rec(lp),
            None => panic!("{}: no ipx/run fixture for {name}", fx.file),
        };
        let text = std::fs::read_to_string(&fx.path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let hook = lines[1..lines.len() - 1].to_vec();
        let r = fx.rec(fx.len() - 1);
        let ctx = fx.file.clone();
        fails.run(&ctx.clone(), || {
            let mut chk = Checker::new(ctx);
            let mut highs = Highs::with_options(options_for(cfg));
            chk.int(
                "pass_model",
                HighsStatus::kOk as i32,
                highs.pass_model(lp) as i32,
            );
            oracle_trace::capture_begin();
            let run_status = highs.run();
            let trace = oracle_trace::capture_end();
            diff_trace(&mut chk, &hook, &trace);
            let a = highs.get_info();
            let sol = highs.get_solution();
            let basis = highs.get_basis();
            let col_status: Vec<i32> = basis.col_status.iter().map(|s| *s as i32).collect();
            let row_status: Vec<i32> = basis.row_status.iter().map(|s| *s as i32).collect();
            chk.int("run_status", r.int("run_status"), run_status as i32)
                .int(
                    "model_status",
                    r.int("model_status"),
                    highs.get_model_status() as i32,
                )
                .f64(
                    "objective_function_value",
                    r.f64("objective_function_value"),
                    a.objective_function_value,
                )
                .int(
                    "ipm_iteration_count",
                    r.int("ipm_iteration_count"),
                    a.ipm_iteration_count,
                )
                .int(
                    "crossover_iteration_count",
                    r.int("crossover_iteration_count"),
                    a.crossover_iteration_count,
                )
                .int(
                    "simplex_iteration_count",
                    r.int("simplex_iteration_count"),
                    a.simplex_iteration_count,
                )
                .int(
                    "primal_solution_status",
                    r.int("primal_solution_status"),
                    a.primal_solution_status,
                )
                .int(
                    "dual_solution_status",
                    r.int("dual_solution_status"),
                    a.dual_solution_status,
                )
                .bool("basis_valid", r.bool("basis_valid"), basis.valid)
                .ints("col_status", &r.ints("col_status"), &col_status)
                .ints("row_status", &r.ints("row_status"), &row_status)
                .f64s("col_value", &r.f64s("col_value"), &sol.col_value)
                .f64s("col_dual", &r.f64s("col_dual"), &sol.col_dual)
                .f64s("row_value", &r.f64s("row_value"), &sol.row_value)
                .f64s("row_dual", &r.f64s("row_dual"), &sol.row_dual);
            Some(chk)
        });
        oracle_trace::capture_end();
        traces += 1;
    }
    assert!(traces > 0, "no ipx/trace fixtures");
    fails.finish("ipx_trace");
}
