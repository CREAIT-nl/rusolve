//! Golden test for manifest row 23, `src/util/hfactor/solve.rs`
//! (HFactor.cpp: ftranL/U, btranL/U, solveHyper, solveMatrixT and the
//! FT/PF/MPF/APF solve variants), entered through `ftranCall`/`btranCall`
//! (row 20).
//!
//! Fixtures: the `solves` lists of the `build` records of
//! `factor/<lp>.jsonl` (afiro committed; others target-only). After each
//! build the oracle solves, on the same factor, six RHS (`unit` x3, a
//! structural `column`, `sparse`, `dense`) at expected densities 1.0, 0.1,
//! 0.01, FTRAN and BTRAN, plus the `std::vector<double>` overloads on the
//! dense RHS. The output HVector is compared in full (count, index order,
//! every array bit, synthetic_tick).
//!
//! The factor is rebuilt here exactly as in `t_hfactor_build`; if that test
//! fails for a basis, fix it first: these solves then fail too. The update
//! solves (FT/PF/MPF/APF) are exercised by `t_hfactor_update`.

mod common;
mod port;

use common::*;
use port::*;

#[test]
fn solves_after_build_match_highs() {
    let mut fails = Failures::default();
    for fx in fixtures_in("factor") {
        let (a, name) = factor_matrix(&fx);
        for rec in recs_of(&fx, "build") {
            let basis = rec.str("basis");
            let (mut f, _, _) = build_factor(&a, &rec.ints("basic_index"));
            let mut chk = Checker::new(format!(
                "{} {name} basis={basis} (factor from HFactor::build; see t_hfactor_build)",
                rec.ctx
            ));
            check_solves(&mut chk, "solves", &rec.solves("solves"), &mut f);
            fails.take(&chk);
        }
    }
    fails.finish("t_hfactor_solve");
}
