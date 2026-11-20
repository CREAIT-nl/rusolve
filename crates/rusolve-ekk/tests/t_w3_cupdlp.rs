//! Golden test for cuPDLP-C through `Highs::run` (`solver=pdlp`): the
//! `w3/cupdlp/<lp>.jsonl` fixtures (`highs-oracle w3-cupdlp`).
//!
//! Per fixture (FORMAT.md, "W3 records"):
//! 1. `Highs::passModel(lp)` + `run` with the `w3_header` overrides,
//!    capturing the port's cuPDLP hook lines (`pdlp/oracle_trace.rs`); they
//!    must equal the fixture's `pdlp_check` / `pdlp_restart` / `pdlp_end`
//!    lines **byte for byte** (the first differing line is reported field
//!    by field).
//! 2. The `w3_result` record is reproduced bit for bit.
//!
//! All these fixtures are target-only; with none generated the test prints a
//! SKIP banner (RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a failure).

mod common;
mod port;

use port::w3::*;
use port::*;

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w3_cupdlp -- --ignored"]
fn w3_cupdlp_matches_highs() {
    let fxs = w3_fixtures("cupdlp");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || {
            replay_w3(fx, &["pdlp_check", "pdlp_restart", "pdlp_end"])
        });
    }
    fails.finish("t_w3_cupdlp");
}
