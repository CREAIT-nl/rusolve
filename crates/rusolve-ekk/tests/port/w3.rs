//! PDLP oracle glue (cuPDLP-C and HiPDLP through `Highs::run`): fixture
//! discovery, the base options with the `pdlp_*` overrides, the
//! `w3_header` LP and the `w3_result` check. Used by
//! `tests/t_w3_{cupdlp,hipdlp}.rs`; record layout in
//! tools/highs-oracle/FORMAT.md ("W3 records").

use crate::common::records::LpRec;
use crate::common::*;
use crate::port::w1::{diff_trace, hook_lines, wave_fixtures};
use crate::port::*;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::pdlp::oracle_trace;
use rusolve_ekk::prelude::*;

/// Every `w3/<dir>/*.jsonl` fixture available (`dir` is `cupdlp` or
/// `hipdlp`). Empty after a SKIP banner when none was generated, unless
/// RUSOLVE_REQUIRE_TARGET_FIXTURES=1.
pub fn w3_fixtures(dir: &str) -> Vec<Fixture> {
    wave_fixtures(&format!("w3/{dir}/"), "PDLP oracle commands")
}

/// The oracle's base options (output off, serial, presolve off,
/// `random_seed=0`) plus the `overrides` of the `w3_header` (`k=v,k=v`).
pub fn w3_options(overrides: &str) -> HighsOptions {
    let mut o = HighsOptions {
        output_flag: false,
        threads: 1,
        parallel: "off".to_string(),
        presolve: "off".to_string(),
        random_seed: 0,
        ..crate::port::profile_options()
    };
    for kv in overrides.split(',') {
        if kv == "default" || kv.is_empty() {
            continue;
        }
        let (k, v) = kv.split_once('=').expect("key=value override");
        let i = || v.parse::<HighsInt>().unwrap();
        match k {
            "solver" => o.solver = v.to_string(),
            "pdlp_iteration_limit" => o.pdlp_iteration_limit = i(),
            "pdlp_scaling_mode" => o.pdlp_scaling_mode = i(),
            "pdlp_ruiz_iterations" => o.pdlp_ruiz_iterations = i(),
            "pdlp_restart_strategy" => o.pdlp_restart_strategy = i(),
            "pdlp_cupdlpc_restart_method" => o.pdlp_cupdlpc_restart_method = i(),
            "pdlp_step_size_strategy" => o.pdlp_step_size_strategy = i(),
            "pdlp_features_off" => o.pdlp_features_off = i(),
            "pdlp_optimality_tolerance" => o.pdlp_optimality_tolerance = v.parse().unwrap(),
            other => panic!("unknown PDLP override {other}"),
        }
    }
    o.sync_log_options();
    o
}

/// The `w3_header` record (first line of every PDLP fixture).
pub fn w3_header(fx: &Fixture) -> Rec<'_> {
    let h = fx.rec(0);
    assert_eq!(h.t(), "w3_header", "{}: first record", fx.file);
    h
}

/// The LP the oracle passed to `Highs::passModel` (as read from the MPS).
pub fn w3_lp(fx: &Fixture) -> HighsLp {
    lp_from_rec(&LpRec::from_rec(&w3_header(fx).obj("lp")))
}

/// The `w3_result` fields (FORMAT.md, "W3 records") against the port after
/// `run`.
pub fn check_w3_result(chk: &mut Checker, r: &Rec, run_status: HighsStatus, highs: &Highs) {
    let info = highs.get_info();
    let sol = highs.get_solution();
    chk.int("run_status", r.int("run_status"), run_status as i32)
        .int(
            "model_status",
            r.int("model_status"),
            highs.get_model_status() as i32,
        )
        .f64(
            "objective_function_value",
            r.f64("objective_function_value"),
            info.objective_function_value,
        )
        .int(
            "pdlp_iteration_count",
            r.int("pdlp_iteration_count"),
            info.pdlp_iteration_count,
        )
        .int(
            "primal_solution_status",
            r.int("primal_solution_status"),
            info.primal_solution_status,
        )
        .int(
            "dual_solution_status",
            r.int("dual_solution_status"),
            info.dual_solution_status,
        )
        .int(
            "num_primal_infeasibilities",
            r.int("num_primal_infeasibilities"),
            info.num_primal_infeasibilities,
        )
        .f64(
            "max_primal_infeasibility",
            r.f64("max_primal_infeasibility"),
            info.max_primal_infeasibility,
        )
        .f64(
            "sum_primal_infeasibilities",
            r.f64("sum_primal_infeasibilities"),
            info.sum_primal_infeasibilities,
        )
        .int(
            "num_dual_infeasibilities",
            r.int("num_dual_infeasibilities"),
            info.num_dual_infeasibilities,
        )
        .f64(
            "max_dual_infeasibility",
            r.f64("max_dual_infeasibility"),
            info.max_dual_infeasibility,
        )
        .f64(
            "sum_dual_infeasibilities",
            r.f64("sum_dual_infeasibilities"),
            info.sum_dual_infeasibilities,
        )
        .bool(
            "basis_valid",
            r.bool("basis_valid"),
            highs.get_basis().valid,
        )
        .f64s("col_value", &r.f64s("col_value"), &sol.col_value)
        .f64s("col_dual", &r.f64s("col_dual"), &sol.col_dual)
        .f64s("row_value", &r.f64s("row_value"), &sol.row_value)
        .f64s("row_dual", &r.f64s("row_dual"), &sol.row_dual);
}

/// Replays one PDLP fixture: `passModel(lp)` + `run` with the header
/// overrides, capturing the port's PDLP hook lines (`hooks`); they must
/// equal the fixture's hook lines byte for byte, then `w3_result` bit for
/// bit.
pub fn replay_w3(fx: &Fixture, hooks: &[&str]) -> Option<Checker> {
    let text = std::fs::read_to_string(&fx.path).unwrap();
    let header = w3_header(fx);
    let want = hook_lines(fx, &text, hooks);
    let result = fx.first("w3_result");
    let mut chk = Checker::new(fx.file.clone());
    let mut highs = Highs::with_options(w3_options(header.str("overrides")));
    chk.int(
        "pass_model",
        HighsStatus::kOk as i32,
        highs.pass_model(w3_lp(fx)) as i32,
    );
    oracle_trace::capture_begin();
    let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| highs.run()));
    let got = oracle_trace::capture_end();
    let run_status = match run {
        Ok(s) => s,
        Err(e) => std::panic::resume_unwind(e),
    };
    diff_trace(&mut chk, &want, &got);
    check_w3_result(&mut chk, &result, run_status, &highs);
    Some(chk)
}
