//! Golden test for iCrash through `Highs::run`: the
//! `w1/icrash/<lp>.jsonl` fixtures (`highs-oracle w1-icrash`).
//!
//! Per fixture: `passModel` + `run` with the header's overrides
//! (`icrash=true,simplex_strategy=4`; the iCrash block of `Highs::run`
//! asserts the primal simplex strategy) reproduces
//! 1. every `icrash_iter` record (`ICrashInfo::details`, without `time`),
//! 2. the `icrash_info` record (`Highs::icrash_info_`, without `total_time`),
//! 3. the `w1_result` record (the crossover solution), bit for bit.

mod common;
mod port;

use common::*;
use port::w1::*;
use port::*;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;

fn replay(fx: &Fixture) -> Option<Checker> {
    let header = w1_header(fx);
    let mut chk = Checker::new(fx.file.clone());
    let mut highs = Highs::with_options(w1_options(header.str("overrides")));
    chk.int(
        "pass_model",
        HighsStatus::kOk as i32,
        highs.pass_model(header_lp(fx)) as i32,
    );
    let run_status = highs.run();
    let info = &highs.icrash_info_;
    let iters: Vec<Rec> = fx.of_type("icrash_iter").collect();
    chk.int("details.len", iters.len() as i32, info.details.len() as i32);
    for (k, (w, a)) in iters.iter().zip(info.details.iter()).enumerate() {
        chk.int(&format!("details[{k}].num"), w.int("num"), a.num)
            .f64(&format!("details[{k}].weight"), w.f64("weight"), a.weight)
            .f64(
                &format!("details[{k}].lambda_norm_2"),
                w.f64("lambda_norm_2"),
                a.lambda_norm_2,
            )
            .f64(
                &format!("details[{k}].lp_objective"),
                w.f64("lp_objective"),
                a.lp_objective,
            )
            .f64(
                &format!("details[{k}].quadratic_objective"),
                w.f64("quadratic_objective"),
                a.quadratic_objective,
            )
            .f64(
                &format!("details[{k}].residual_norm_2"),
                w.f64("residual_norm_2"),
                a.residual_norm_2,
            );
    }
    let w = fx.first("icrash_info");
    chk.int(
        "num_iterations",
        w.int("num_iterations"),
        info.num_iterations,
    )
    .f64(
        "final_lp_objective",
        w.f64("final_lp_objective"),
        info.final_lp_objective,
    )
    .f64(
        "final_quadratic_objective",
        w.f64("final_quadratic_objective"),
        info.final_quadratic_objective,
    )
    .f64(
        "final_residual_norm_2",
        w.f64("final_residual_norm_2"),
        info.final_residual_norm_2,
    )
    .f64(
        "starting_weight",
        w.f64("starting_weight"),
        info.starting_weight,
    )
    .f64("final_weight", w.f64("final_weight"), info.final_weight)
    .f64s("x_values", &w.f64s("x_values"), &info.x_values);
    check_w1_result(&mut chk, &fx.first("w1_result"), run_status, &highs);
    Some(chk)
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w1_icrash -- --ignored"]
fn w1_icrash_matches_highs() {
    let fxs = w1_fixtures("icrash");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || replay(fx));
    }
    fails.finish("t_w1_icrash");
}
