//! Golden test for manifest row 26, `src/util/hfactor/refactor.rs`
//! (HFactorRefactor.cpp: `HFactor::rebuild`, `RefactorInfo::clear`).
//!
//! Fixtures: the `refactor` records of `factor/<lp>.jsonl` (afiro committed;
//! others target-only), written by the second half of `buildCase` in
//! `tools/highs-oracle/src/factor.cpp`: a fresh FT factor is set up on the
//! same basis, given the pivot sequence recorded by the first build
//! (`refactor_pivot_row/var/type`, `refactor_build_synthetic_tick`, all in
//! the record) with `refactor_info_.use = true`, and `build()` then takes
//! the `rebuild` path instead of INVERT. Compared: the full build result
//! (as `t_hfactor_build`) and the solve batch on the rebuilt factor.
//!
//! Only full-rank bases have a refactor record (the oracle skips the
//! rank-deficient ones).

mod common;
mod port;

use common::*;
use port::*;
use rusolve_ekk::prelude::*;

#[test]
fn refactor_rebuild_matches_highs() {
    let mut fails = Failures::default();
    for fx in fixtures_in("factor") {
        let (a, name) = factor_matrix(&fx);
        for rec in recs_of(&fx, "refactor") {
            let basis = rec.str("basis");
            let mut bi = rec.ints("basic_index");
            let mut g = HFactor::default();
            setup_factor(&mut g, &a, kUpdateMethodFt);
            g.refactor_info_ = RefactorInfo {
                r#use: true,
                pivot_row: rec.ints("refactor_pivot_row"),
                pivot_var: rec.ints("refactor_pivot_var"),
                pivot_type: rec.i8s("refactor_pivot_type"),
                build_synthetic_tick: rec.f64("refactor_build_synthetic_tick"),
            };
            let rd = g.build(&mut bi, None);
            let mut chk = Checker::new(format!(
                "{} {name} basis={basis} (HFactor::build with refactor_info_.use = true -> HFactor::rebuild)",
                rec.ctx
            ));
            check_build_result(&mut chk, &rec.build_result("result"), &g, rd, &bi);
            check_solves(&mut chk, "solves", &rec.solves("solves"), &mut g);
            fails.take(&chk);
        }
    }
    fails.finish("t_hfactor_refactor");
}
