//! Golden test for `Highs::getIis`: the
//! `w1/iis/<lp>.jsonl` fixtures (`highs-oracle w1-iis`).
//!
//! Per `iis_result` record (one per configuration `cfg`, each on a fresh
//! `Highs`): `passModel` + `run` + `getIis` with the record's `overrides`
//! reproduce the run/model status, the `getIis` status and every `HighsIis`
//! field the oracle dumps (`valid_`, `status_`, `strategy_`, the index /
//! bound / status vectors and the `info_` counts; the times are not
//! compared).

mod common;
mod port;

use common::*;
use port::w1::*;
use port::*;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;

fn replay(fx: &Fixture, r: &Rec) -> Option<Checker> {
    let mut chk = Checker::new(format!("{} cfg={}", fx.file, r.str("cfg")));
    let mut highs = Highs::with_options(w1_options(r.str("overrides")));
    chk.int(
        "pass_model",
        HighsStatus::kOk as i32,
        highs.pass_model(header_lp(fx)) as i32,
    );
    let run_status = highs.run();
    chk.int("run_status", r.int("run_status"), run_status as i32)
        .int(
            "model_status",
            r.int("model_status"),
            highs.get_model_status() as i32,
        );
    let mut iis = HighsIis::default();
    let status = highs.get_iis(&mut iis);
    let info = r.obj("info");
    chk.int("get_iis_status", r.int("get_iis_status"), status as i32)
        .bool("valid", r.bool("valid"), iis.valid_)
        .int("status", r.int("status"), iis.status_)
        .int("strategy", r.int("strategy"), iis.strategy_)
        .ints("col_index", &r.ints("col_index"), &iis.col_index_)
        .ints("row_index", &r.ints("row_index"), &iis.row_index_)
        .ints("col_bound", &r.ints("col_bound"), &iis.col_bound_)
        .ints("row_bound", &r.ints("row_bound"), &iis.row_bound_)
        .ints("col_status", &r.ints("col_status"), &iis.col_status_)
        .ints("row_status", &r.ints("row_status"), &iis.row_status_)
        .int(
            "info.num_lp_solved",
            info.int("num_lp_solved"),
            iis.info_.num_lp_solved,
        )
        .int(
            "info.sum_simplex_iteration_counts",
            info.int("sum_simplex_iteration_counts"),
            iis.info_.sum_simplex_iteration_counts,
        )
        .int(
            "info.min_simplex_iteration_count",
            info.int("min_simplex_iteration_count"),
            iis.info_.min_simplex_iteration_count,
        )
        .int(
            "info.max_simplex_iteration_count",
            info.int("max_simplex_iteration_count"),
            iis.info_.max_simplex_iteration_count,
        );
    Some(chk)
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w1_iis -- --ignored"]
fn w1_iis_matches_highs() {
    let fxs = w1_fixtures("iis");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        w1_header(fx);
        for r in recs_of(fx, "iis_result") {
            fails.run(&r.ctx.clone(), || replay(fx, &r));
        }
    }
    fails.finish("t_w1_iis");
}
