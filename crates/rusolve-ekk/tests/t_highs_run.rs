//! Golden test for the driver around `HEkk::solve`: `Highs::passModel` +
//! `Highs::run` (src/lp_data/highs_run.rs), `solveLp` (highs_solve.rs),
//! `HApp::solveLpSimplex` (simplex/happ.rs), `considerScaling`/`scaleLp`/
//! `unscaleSolution` (lp_scale.rs) and `assessLp` (lp_assess.rs).
//!
//! Fixtures: every `simplex/<lp>.{dual,primal}.jsonl` trace. The oracle ran
//! `Highs::run` on the MPS relaxation (`tools/highs-oracle/src/simplex.cpp`
//! runSimplex) and dumped the `result` record (unscaled solution, basis,
//! info). The test needs the LP HiGHS was given; the trace holds it only in
//! scaled form (`ekk_solve_enter.lp`), so it is rebuilt with
//! `HighsLp::unapplyScale`. The scale factors are powers of two, so the
//! unscaled LP is bit-exact (the post-`assessLp` LP, which `assessLp` leaves
//! unchanged).
//!
//! Per trace:
//! 1. `considerScaling` on the rebuilt LP reproduces `ekk_solve_enter.lp`
//!    (scale factors, scaled matrix/bounds/costs) bit-exactly.
//! 2. `Highs::passModel` + `Highs::run` reproduces the `result` record:
//!    run/model status, objective, iteration count, solution statuses and
//!    infeasibility counts, basis and every solution vector (bit-exact).

mod common;
mod port;

use common::events::{EkkSolveEnter, RunResult};
use common::*;
use port::*;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::lp_data::lp_scale::consider_scaling;
use rusolve_ekk::prelude::*;

/// The LP the oracle passed to HiGHS, rebuilt from the scaled LP HEkk saw.
fn unscaled_lp(enter: &EkkSolveEnter) -> HighsLp {
    let mut lp = lp_from_rec(&enter.lp);
    lp.is_moved_ = false;
    lp.unapply_scale();
    lp.clear_scale();
    lp
}

fn run_gate(fx: &Fixture) -> Option<Checker> {
    let enter = EkkSolveEnter::from_rec(&recs_of(fx, "ekk_solve_enter")[0]);
    let result = RunResult::from_rec(&recs_of(fx, "result")[0]);
    let options = options_from_json(&enter.options, &fx.file);
    let lp = unscaled_lp(&enter);
    let mut chk = Checker::new(fx.file.clone());

    // 1. Scaling reproduces the LP HEkk received.
    let mut scaled = lp.clone();
    consider_scaling(&options, &mut scaled);
    let (w, a) = (&enter.lp, lp_rec_of(&scaled));
    chk.int("lp.num_col", w.num_col, a.num_col)
        .int("lp.num_row", w.num_row, a.num_row)
        .int("lp.sense", w.sense, a.sense)
        .f64("lp.offset", w.offset, a.offset)
        .f64s("lp.col_cost", &w.col_cost, &a.col_cost)
        .f64s("lp.col_lower", &w.col_lower, &a.col_lower)
        .f64s("lp.col_upper", &w.col_upper, &a.col_upper)
        .f64s("lp.row_lower", &w.row_lower, &a.row_lower)
        .f64s("lp.row_upper", &w.row_upper, &a.row_upper)
        .int("lp.a.format", w.a.format, a.a.format)
        .ints("lp.a.start", &w.a.start, &a.a.start)
        .ints("lp.a.index", &w.a.index, &a.a.index)
        .f64s("lp.a.value", &w.a.value, &a.a.value)
        .bool("lp.is_scaled", w.is_scaled, a.is_scaled)
        .int("lp.scale.strategy", w.scale.strategy, a.scale.strategy)
        .bool(
            "lp.scale.has_scaling",
            w.scale.has_scaling,
            a.scale.has_scaling,
        )
        .int("lp.scale.num_col", w.scale.num_col, a.scale.num_col)
        .int("lp.scale.num_row", w.scale.num_row, a.scale.num_row)
        .f64("lp.scale.cost", w.scale.cost, a.scale.cost)
        .f64s("lp.scale.col", &w.scale.col, &a.scale.col)
        .f64s("lp.scale.row", &w.scale.row, &a.scale.row);

    // 2. Highs::passModel + Highs::run.
    let mut highs = Highs::with_options(options);
    let pass_status = highs.pass_model(lp);
    chk.int(
        "pass_model status",
        HighsStatus::kOk as i32,
        pass_status as i32,
    );
    let run_status = highs.run();
    let info = highs.get_info();
    let sol = highs.get_solution();
    let basis = highs.get_basis();
    chk.int("run_status", result.run_status, run_status as i32)
        .int(
            "model_status",
            result.model_status,
            highs.get_model_status() as i32,
        )
        .f64(
            "objective_function_value",
            result.objective_function_value,
            info.objective_function_value,
        )
        .int(
            "simplex_iteration_count",
            result.simplex_iteration_count,
            info.simplex_iteration_count,
        )
        .int(
            "primal_solution_status",
            result.primal_solution_status,
            info.primal_solution_status,
        )
        .int(
            "dual_solution_status",
            result.dual_solution_status,
            info.dual_solution_status,
        )
        .int(
            "num_primal_infeasibilities",
            result.num_primal_infeasibilities,
            info.num_primal_infeasibilities,
        )
        .int(
            "num_dual_infeasibilities",
            result.num_dual_infeasibilities,
            info.num_dual_infeasibilities,
        )
        .bool("basis_valid", result.basis_valid, basis.valid);
    let col_status: Vec<i32> = basis.col_status.iter().map(|s| *s as i32).collect();
    let row_status: Vec<i32> = basis.row_status.iter().map(|s| *s as i32).collect();
    chk.ints("col_status", &result.col_status, &col_status)
        .ints("row_status", &result.row_status, &row_status)
        .f64s("col_value", &result.col_value, &sol.col_value)
        .f64s("col_dual", &result.col_dual, &sol.col_dual)
        .f64s("row_value", &result.row_value, &sol.row_value)
        .f64s("row_dual", &result.row_dual, &sol.row_dual);
    Some(chk)
}

fn run_all(files: Vec<String>, what: &str) {
    let mut fails = Failures::default();
    for file in files {
        let Some(fx) = load(&file) else { continue };
        fails.run(&file, || run_gate(&fx));
    }
    fails.finish(what);
}

fn committed() -> Vec<String> {
    let mut v = committed_trace_files(TraceKind::Dual);
    v.extend(committed_trace_files(TraceKind::Primal));
    v
}

#[test]
fn highs_run_matches_highs_committed() {
    run_all(committed(), "t_highs_run (committed traces)");
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_highs_run -- --ignored"]
fn highs_run_matches_highs_all() {
    let mut v = trace_files(TraceKind::Dual);
    v.extend(trace_files(TraceKind::Primal));
    run_all(v, "t_highs_run (all traces)");
}
