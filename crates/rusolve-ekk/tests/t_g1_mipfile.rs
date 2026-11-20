//! Golden test for the MIP improving-solution file: the `g1/mipfile/<case>.jsonl` fixtures
//! (`highs-oracle g1-mipfile`): `run` with `mip_improving_solution_file` set: the file text is compared byte for byte (`g1_file`).
//!
//! Per fixture (FORMAT.md, "G1 records"): `Highs::passModel(lp)` with the
//! `g1_header` overrides, then each `g1_call` in order (its return status
//! and outputs bit for bit), each `g1_result` against the port state.
//!
//! All these fixtures are target-only; with none generated the test prints a
//! SKIP banner (RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a failure).

mod common;
mod port;

use port::g1::*;

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_g1_mipfile -- --ignored"]
fn g1_mipfile_matches_highs() {
    run_g1("mipfile", "t_g1_mipfile");
}
