//! Golden test for manifest rows 19-22, `src/util/hfactor/{utils,setup,
//! build,kernel}.rs` (HFactor.cpp: setupGeneral, build, buildSimple,
//! buildKernel, buildHandleRankDeficiency, buildFinish, ...).
//!
//! Fixtures: the `build` records of `factor/<lp>.jsonl` (afiro committed;
//! adlittle, blend, kb2, sc50a, share2b target-only), written by
//! `tools/highs-oracle/src/factor.cpp` (`buildCase`) on the *scaled* LP
//! matrix HEkk factorises. Five bases per LP: `slack`, `mid` (HiGHS dual
//! simplex stopped at half the optimal iteration count), `optimal`, and two
//! rank-deficient ones (`rank_deficient_dup`, `rank_deficient_mixed`).
//!
//! Replay: `HFactor::setupGeneral(num_col, num_row, num_row, a, ...,
//! kDefaultPivotThreshold, kDefaultPivotTolerance, kHighsDebugLevelMin,
//! NULL, true, kUpdateMethodFt)` then `build()`; compare the return value,
//! the permuted `basic_index`, the no-pivot lists, ticks, counts, the
//! recorded `refactor_info_` pivot sequence and every INVERT array
//! (`getInvert()`). FTRAN/BTRAN after the build are in `t_hfactor_solve`.

mod common;
mod port;

use common::*;
use port::*;

#[test]
fn build_matches_highs() {
    let mut fails = Failures::default();
    for fx in fixtures_in("factor") {
        let (a, name) = factor_matrix(&fx);
        for rec in recs_of(&fx, "build") {
            let basis = rec.str("basis");
            let (f, rd, bi) = build_factor(&a, &rec.ints("basic_index"));
            let mut chk = Checker::new(format!(
                "{} {name} basis={basis} (HFactor::build, update_method=FT)",
                rec.ctx
            ));
            check_build_result(&mut chk, &rec.build_result("result"), &f, rd, &bi);
            fails.take(&chk);
        }
    }
    fails.finish("t_hfactor_build");
}
