//! Golden test for QUASS through `Highs::run`: the
//! `w1/qp/<inst>.jsonl` fixtures (`highs-oracle w1-qp`).
//!
//! Per fixture (FORMAT.md, "W1 records"):
//! 1. `Highs::passModel(lp)` + `passHessian(hessian)` + `run` with the
//!    `w1_header` overrides, capturing the port's QUASS hook lines
//!    (`qpsolver/oracle_trace.rs`); they must equal the fixture's
//!    `qp_iter` / `qp_end` lines **byte for byte** (the first differing line
//!    is reported field by field).
//! 2. The `w1_result` record is reproduced bit for bit.
//!
//! All these fixtures are target-only; with none generated the test prints a
//! SKIP banner (RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a failure).

mod common;
mod port;

use common::*;
use port::w1::*;
use port::*;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use rusolve_ekk::qpsolver::oracle_trace;

fn replay(fx: &Fixture) -> Option<Checker> {
    let text = std::fs::read_to_string(&fx.path).unwrap();
    let header = w1_header(fx);
    let want = hook_lines(fx, &text, &["qp_iter", "qp_end"]);
    let result = fx.first("w1_result");
    let mut chk = Checker::new(fx.file.clone());
    let mut highs = Highs::with_options(w1_options(header.str("overrides")));
    chk.int(
        "pass_model",
        HighsStatus::kOk as i32,
        highs.pass_model(header_lp(fx)) as i32,
    );
    chk.int(
        "pass_hessian",
        HighsStatus::kOk as i32,
        highs.pass_hessian(header_hessian(fx)) as i32,
    );
    oracle_trace::capture_begin();
    let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| highs.run()));
    let got = oracle_trace::capture_end();
    let run_status = match run {
        Ok(s) => s,
        Err(e) => std::panic::resume_unwind(e),
    };
    diff_trace(&mut chk, &want, &got);
    check_w1_result(&mut chk, &result, run_status, &highs);
    Some(chk)
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w1_qp -- --ignored"]
fn w1_qp_matches_highs() {
    let fxs = w1_fixtures("qp");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || replay(fx));
    }
    fails.finish("t_w1_qp");
}
