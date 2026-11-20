//! HiPO oracle glue (HiPO through `Highs::run`): fixture discovery, the base
//! options with the `hipo_*` overrides (and `hipo_available = true`), the
//! `w4_header` LP, the `w4_result` check and the standalone replay of one
//! `hipo_ord` record through the ported orderings. Used by
//! `tests/t_w4_{orderings,hipo}.rs`; record layout in
//! tools/highs-oracle/FORMAT.md ("W4 records").

use crate::common::records::LpRec;
use crate::common::*;
use crate::port::w1::{diff_trace, hook_lines, wave_fixtures};
use crate::port::*;
use rusolve_ekk::ipm::hipo::factorhighs::kMetisSeed;
use rusolve_ekk::ipm::hipo::oracle_trace;
use rusolve_ekk::ipm::hipo::orderings::amd::{AMD_CONTROL, AMD_INFO, AMD_OK};
use rusolve_ekk::ipm::hipo::orderings::extras;
use rusolve_ekk::ipm::hipo::orderings::metis::{
    METIS_NOPTIONS, METIS_OK, METIS_OPTION_DBGLVL, METIS_OPTION_NO2HOP, METIS_OPTION_SEED,
};
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;

/// Every `w4/<dir>/*.jsonl` fixture available (`dir` is `orderings` or
/// `hipo`). Empty after a SKIP banner when none was generated, unless
/// RUSOLVE_REQUIRE_TARGET_FIXTURES=1.
pub fn w4_fixtures(dir: &str) -> Vec<Fixture> {
    wave_fixtures(&format!("w4/{dir}/"), "HiPO oracle commands")
}

/// The oracle's base options (output off, serial, presolve off,
/// `random_seed=0`, HiPO available) plus the `overrides` of the `w4_header`
/// (`k=v,k=v`). The oracle sets `threads=1` and `parallel=off`, so HiPO's
/// `setParallel` resolves to serial.
pub fn w4_options(overrides: &str) -> HighsOptions {
    let mut o = HighsOptions {
        output_flag: false,
        threads: 1,
        parallel: "off".to_string(),
        presolve: "off".to_string(),
        random_seed: 0,
        hipo_available: true,
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
            "hipo_system" => o.hipo_system = v.to_string(),
            "hipo_parallel_type" => o.hipo_parallel_type = v.to_string(),
            "hipo_ordering" => o.hipo_ordering = v.to_string(),
            "hipo_block_size" => o.hipo_block_size = i(),
            "run_crossover" => o.run_crossover = v.to_string(),
            "ipm_iteration_limit" => o.ipm_iteration_limit = i(),
            "ipm_optimality_tolerance" => o.ipm_optimality_tolerance = v.parse().unwrap(),
            other => panic!("unknown HiPO override {other}"),
        }
    }
    o.sync_log_options();
    o
}

/// The `w4_header` record (first line of every HiPO fixture).
pub fn w4_header(fx: &Fixture) -> Rec<'_> {
    let h = fx.rec(0);
    assert_eq!(h.t(), "w4_header", "{}: first record", fx.file);
    h
}

/// The LP the oracle passed to `Highs::passModel` (as read from the MPS).
pub fn w4_lp(fx: &Fixture) -> HighsLp {
    lp_from_rec(&LpRec::from_rec(&w4_header(fx).obj("lp")))
}

/// The `w4_result` fields (FORMAT.md, "W4 records") against the port after
/// `run`.
pub fn check_w4_result(chk: &mut Checker, r: &Rec, run_status: HighsStatus, highs: &Highs) {
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
            "ipm_iteration_count",
            r.int("ipm_iteration_count"),
            info.ipm_iteration_count,
        )
        .int(
            "crossover_iteration_count",
            r.int("crossover_iteration_count"),
            info.crossover_iteration_count,
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

/// Replays one HiPO fixture: `passModel(lp)` + `run` with the header
/// overrides, capturing the port's HiPO hook lines (`hooks`); they must
/// equal the fixture's hook lines byte for byte, then `w4_result` bit for
/// bit.
pub fn replay_w4(fx: &Fixture, hooks: &[&str]) -> Option<Checker> {
    let text = std::fs::read_to_string(&fx.path).unwrap();
    let header = w4_header(fx);
    let want = hook_lines(fx, &text, hooks);
    let result = fx.first("w4_result");
    let mut chk = Checker::new(fx.file.clone());
    let mut highs = Highs::with_options(w4_options(header.str("overrides")));
    chk.int(
        "pass_model",
        HighsStatus::kOk as i32,
        highs.pass_model(w4_lp(fx)) as i32,
    );
    oracle_trace::capture_begin();
    let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| highs.run()));
    let got = oracle_trace::capture_end();
    let run_status = match run {
        Ok(s) => s,
        Err(e) => std::panic::resume_unwind(e),
    };
    diff_trace(&mut chk, &want, &got);
    check_w4_result(&mut chk, &result, run_status, &highs);
    Some(chk)
}

/// Re-runs the ordering of one `hipo_ord` record on its recorded pattern
/// (`n`, `ptr`, `rows`: the full symmetric CSC without diagonal that
/// `FactorHiGHSSolver::chooseOrdering` passes), with the same arguments as
/// the C++ lambda `run_ordering_and_analyse` (FactorHiGHSSolver.cpp
/// 397-433), and checks `perm` and `failure`.
pub fn replay_hipo_ord(chk: &mut Checker, r: &Rec) {
    let ordering = r.str("ordering");
    let mut n: HighsInt = r.int("n");
    let ptr = r.ints("ptr");
    let rows = r.ints("rows");
    let mut perm = vec![0 as HighsInt; n as usize];
    let failure = match ordering {
        "metis" => {
            let mut options = [0 as HighsInt; METIS_NOPTIONS];
            extras::metis::set_default_options(&mut options);
            options[METIS_OPTION_SEED] = kMetisSeed;
            options[METIS_OPTION_DBGLVL] = 0;
            options[METIS_OPTION_NO2HOP] = 1;
            let mut iperm = vec![0 as HighsInt; n as usize];
            let status = extras::metis::nodeND(
                &mut n,
                &ptr,
                &rows,
                None,
                Some(&options),
                &mut perm,
                &mut iperm,
            );
            status != METIS_OK
        }
        "amd" => {
            let mut control = [0.0; AMD_CONTROL];
            extras::amd::set_defaults(&mut control);
            let mut info = [0.0; AMD_INFO];
            let status =
                extras::amd::order(n, &ptr, &rows, &mut perm, Some(&control), Some(&mut info));
            status != AMD_OK
        }
        "rcm" => {
            let adj_num = *ptr.last().unwrap();
            extras::rcm::genrcm(n, adj_num, &ptr, &rows, &mut perm) != 0
        }
        other => panic!("unknown hipo_ord ordering {other}"),
    };
    chk.bool("failure", r.int("failure") != 0, failure);
    chk.ints("perm", &r.ints("perm"), &perm);
}
