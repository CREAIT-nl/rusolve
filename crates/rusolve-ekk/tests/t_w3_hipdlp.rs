//! Golden test for HiPDLP through `Highs::run` (`solver=hipdlp`): the
//! `w3/hipdlp/<lp>.jsonl` fixtures (`highs-oracle w3-hipdlp`).
//!
//! Per fixture (FORMAT.md, "W3 records"):
//! 1. `Highs::passModel(lp)` + `run` with the `w3_header` overrides,
//!    capturing the port's HiPDLP hook lines (`pdlp/oracle_trace.rs`); they
//!    must equal the fixture's `hipdlp_check` / `hipdlp_restart` /
//!    `hipdlp_end` lines **byte for byte** (the first differing line is
//!    reported field by field).
//! 2. The `w3_result` record is reproduced bit for bit.
//!
//! All these fixtures are target-only; with none generated the test prints a
//! SKIP banner (RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a failure).

mod common;
mod port;

use port::w3::*;
use port::*;

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w3_hipdlp -- --ignored"]
fn w3_hipdlp_matches_highs() {
    let fxs = w3_fixtures("hipdlp");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || {
            replay_w3(fx, &["hipdlp_check", "hipdlp_restart", "hipdlp_end"])
        });
    }
    fails.finish("t_w3_hipdlp");
}
