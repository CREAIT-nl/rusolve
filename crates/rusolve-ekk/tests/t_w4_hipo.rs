//! Golden test for HiPO through `Highs::run` (`solver=hipo`): the `w4/hipo/<lp>.jsonl` and
//! `w4/orderings/<lp>.jsonl` fixtures (`highs-oracle w4-orderings` / `w4-hipo`).
//!
//! Per fixture (FORMAT.md, "W4 records"):
//! 1. `Highs::passModel(lp)` + `run` with the `w4_header` overrides and
//!    `hipo_available = true`, capturing the port's HiPO hook lines
//!    (`ipm/hipo/oracle_trace.rs`); they must equal the fixture's
//!    `hipo_ord` / `hipo_fact` / `hipo_iter` / `hipo_end` lines **byte for
//!    byte** (the first differing line is reported field by field).
//! 2. The `w4_result` record is reproduced bit for bit.
//!
//! All these fixtures are target-only; with none generated the test prints a
//! SKIP banner (RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a failure).

mod common;
mod port;

use port::w4::*;
use port::*;

const HOOKS: &[&str] = &["hipo_ord", "hipo_fact", "hipo_iter", "hipo_end"];

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w4_hipo -- --ignored"]
fn w4_hipo_matches_highs() {
    let fxs = w4_fixtures("hipo");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || replay_w4(fx, HOOKS));
    }
    fails.finish("t_w4_hipo");
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w4_hipo -- --ignored"]
fn w4_hipo_orderings_runs_match_highs() {
    let fxs = w4_fixtures("orderings");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || replay_w4(fx, HOOKS));
    }
    fails.finish("t_w4_hipo_orderings");
}
