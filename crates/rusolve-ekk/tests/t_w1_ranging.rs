//! Golden test for `Highs::getRanging`: the
//! `w1/ranging/<lp>.jsonl` fixtures (`highs-oracle w1-ranging`).
//!
//! Per fixture: `passModel` + `run` (simplex, base options) reproduces the
//! `w1_result` record, then `getRanging` reproduces `ranging_status` and the
//! six `ranging` records (`value_`, `objective_`, `in_var_`, `ou_var_`) bit
//! for bit.

mod common;
mod port;

use common::*;
use port::w1::*;
use port::*;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;

fn record<'a>(r: &'a HighsRanging, which: &str) -> &'a HighsRangingRecord {
    match which {
        "col_cost_up" => &r.col_cost_up,
        "col_cost_dn" => &r.col_cost_dn,
        "col_bound_up" => &r.col_bound_up,
        "col_bound_dn" => &r.col_bound_dn,
        "row_bound_up" => &r.row_bound_up,
        "row_bound_dn" => &r.row_bound_dn,
        w => panic!("unknown ranging record {w}"),
    }
}

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
    check_w1_result(&mut chk, &fx.first("w1_result"), run_status, &highs);
    let mut ranging = HighsRanging::default();
    let status = highs.get_ranging(&mut ranging);
    let want = fx.first("ranging_status");
    chk.int("ranging_status", want.int("status"), status as i32)
        .bool("ranging.valid", want.bool("valid"), ranging.valid);
    let recs = recs_of(fx, "ranging");
    chk.int("ranging records", 6, recs.len() as i32);
    for r in recs {
        let which = r.str("which");
        let a = record(&ranging, which);
        chk.f64s(&format!("{which}.value"), &r.f64s("value"), &a.value_)
            .f64s(
                &format!("{which}.objective"),
                &r.f64s("objective"),
                &a.objective_,
            )
            .ints(&format!("{which}.in_var"), &r.ints("in_var"), &a.in_var_)
            .ints(&format!("{which}.ou_var"), &r.ints("ou_var"), &a.ou_var_);
    }
    Some(chk)
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w1_ranging -- --ignored"]
fn w1_ranging_matches_highs() {
    let fxs = w1_fixtures("ranging");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || replay(fx));
    }
    fails.finish("t_w1_ranging");
}
