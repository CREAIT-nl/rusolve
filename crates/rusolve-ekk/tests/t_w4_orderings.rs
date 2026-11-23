//! Golden test for the HiPO fill-reducing orderings: every `hipo_ord` record of the `w4/orderings/<lp>.jsonl` and
//! `w4/hipo/<lp>.jsonl` fixtures (`highs-oracle w4-orderings` / `w4-hipo`).
//!
//! Per record (FORMAT.md, "W4 records"): the recorded pattern (`n`, `ptr`,
//! `rows`) goes through the ported METIS / AMD / RCM with the arguments of
//! `FactorHiGHSSolver::chooseOrdering`; `perm` and `failure` must match
//! exactly. This isolates the orderings from the IPM: a mismatch here is a
//! METIS / AMD / RCM bug, never a HiPO one.
//!
//! All these fixtures are target-only; with none generated the test prints a
//! SKIP banner (RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a failure).

mod common;
mod port;

use common::*;
use port::w4::*;
use port::*;

/// Replays the `hipo_ord` records of one ordering (`metis`, `amd`, `rcm`).
/// One test per library, so that each ordering can be
/// run on its own: `-- --ignored w4_orderings_<lib>`.
fn orderings_match_highs(which: &str) {
    let mut fxs = w4_fixtures("orderings");
    fxs.extend(w4_fixtures("hipo"));
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        for (k, r) in fx.of_type("hipo_ord").enumerate() {
            if r.str("ordering") != which {
                continue;
            }
            let ctx = format!(
                "{} hipo_ord[{k}] {} {}",
                fx.file,
                r.str("nla"),
                r.str("ordering")
            );
            fails.run(&ctx, || {
                let mut chk = Checker::new(ctx.clone());
                replay_hipo_ord(&mut chk, &r);
                Some(chk)
            });
        }
    }
    fails.finish(&format!("t_w4_orderings ({which})"));
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w4_orderings -- --ignored"]
fn w4_orderings_metis_match_highs() {
    orderings_match_highs("metis");
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w4_orderings -- --ignored"]
fn w4_orderings_amd_match_highs() {
    orderings_match_highs("amd");
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --release --test t_w4_orderings -- --ignored"]
fn w4_orderings_rcm_match_highs() {
    orderings_match_highs("rcm");
}
