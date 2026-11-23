//! Golden test for the 1.16.0 clique-row normalisation (HiGHS #3286 +
//! #3318): the `g1/clique_normalise/<case>.jsonl`
//! fixtures (`highs-oracle g1-clique-normalise`, tools/highs-oracle/src/clqn.cpp).
//!
//! Per fixture: the `g1_header` model goes into a HighsMipSolver (new_data,
//! init, setupDomainPropagation), then `HPresolve::ok_set_input_mip` and one
//! `normalise_clique_rows` call with `mip_clique_normalise` on, as the
//! 1.16.0 TestPresolveRules normalisation cases set it up. The
//! `clqn_result` record is checked bit for bit: the call's result, the
//! postsolve reduction count, column bounds and deletions, and every row
//! (bounds, deletion, nonzeros in `getRowVector` order, implied activity
//! bounds). The same run with `mip_coef_change_impl_bounds` on (#3318's
//! implied-bound bookkeeping) must give the same state.
//!
//! The fixtures exist only for the 1.16.0 oracle (normaliseCliqueRows is
//! after 1.15.1): `ORACLE_HIGHS=1.16.0 tools/highs-oracle/gen_fixtures.sh`,
//! then `RUSOLVE_FIXTURES_DIR=target/highs-fixtures-1.16.0`. Without them the
//! test prints a SKIP banner (RUSOLVE_REQUIRE_TARGET_FIXTURES=1 makes that a
//! failure).

mod common;
mod port;

use common::*;
use port::g1::*;
use port::*;
use rusolve_ekk::prelude::*;

/// The presolve state after the call, in the `clqn_result` layout.
struct State {
    result: i32,
    num_reductions: i64,
    col_lower: Vec<f64>,
    col_upper: Vec<f64>,
    col_deleted: Vec<i32>,
    row_lower: Vec<f64>,
    row_upper: Vec<f64>,
    row_deleted: Vec<i32>,
    row_start: Vec<i32>,
    row_index: Vec<i32>,
    row_value: Vec<f64>,
    sum_lower: Vec<f64>,
    sum_upper: Vec<f64>,
    num_inf_sum_lower: Vec<i32>,
    num_inf_sum_upper: Vec<i32>,
}

fn normalise(lp: &HighsLp, coef_change: bool) -> State {
    let mut options = g1_options("default");
    options.mip_clique_normalise = true;
    options.mip_coef_change_impl_bounds = coef_change;
    let mut cb = HighsCallback::default();
    let solution = HighsSolution::default();
    let mut ms = HighsMipSolver::new(&mut cb, &options, lp, &solution, false, 0);
    ms.timer_.start();
    ms.new_data();
    ms.init();
    ms.setup_domain_propagation();
    let mut postsolve_stack = std::mem::take(&mut ms.mipdata_.post_solve_stack);

    let mut presolve = HPresolve::default();
    assert!(presolve.ok_set_input_mip(Box::new(ms), -1));
    let result = presolve.normalise_clique_rows(&mut postsolve_stack);

    let m = &presolve.model;
    let mut s = State {
        result: result as i32,
        num_reductions: postsolve_stack.num_reductions() as i64,
        col_lower: m.col_lower_.clone(),
        col_upper: m.col_upper_.clone(),
        col_deleted: presolve
            .col_deleted
            .iter()
            .map(|&d| (d != 0) as i32)
            .collect(),
        row_lower: m.row_lower_.clone(),
        row_upper: m.row_upper_.clone(),
        row_deleted: presolve
            .row_deleted
            .iter()
            .map(|&d| (d != 0) as i32)
            .collect(),
        row_start: Vec::new(),
        row_index: Vec::new(),
        row_value: Vec::new(),
        sum_lower: Vec::new(),
        sum_upper: Vec::new(),
        num_inf_sum_lower: Vec::new(),
        num_inf_sum_upper: Vec::new(),
    };
    for row in 0..m.num_row_ {
        s.row_start.push(s.row_index.len() as i32);
        if presolve.row_deleted[row as usize] != 0 {
            s.sum_lower.push(0.0);
            s.sum_upper.push(0.0);
            s.num_inf_sum_lower.push(0);
            s.num_inf_sum_upper.push(0);
            continue;
        }
        let mut it = presolve.get_row_vector(row);
        while let Some(nz) = it.next_nz(&presolve) {
            s.row_index.push(nz.index());
            s.row_value.push(nz.value());
        }
        let b = &presolve.implied_row_bounds;
        s.sum_lower.push(b.get_sum_lower(row, 0.0));
        s.sum_upper.push(b.get_sum_upper(row, 0.0));
        s.num_inf_sum_lower.push(b.get_num_inf_sum_lower(row));
        s.num_inf_sum_upper.push(b.get_num_inf_sum_upper(row));
    }
    s.row_start.push(s.row_index.len() as i32);
    s
}

fn check(chk: &mut Checker, tag: &str, r: &Rec, s: &State) {
    let k = |f: &str| format!("{tag}.{f}");
    chk.int(&k("result"), r.int("result"), s.result);
    chk.int(
        &k("num_reductions"),
        r.i64("num_reductions"),
        s.num_reductions,
    );
    chk.f64s(&k("col_lower"), &r.f64s("col_lower"), &s.col_lower);
    chk.f64s(&k("col_upper"), &r.f64s("col_upper"), &s.col_upper);
    chk.ints(&k("col_deleted"), &r.ints("col_deleted"), &s.col_deleted);
    chk.f64s(&k("row_lower"), &r.f64s("row_lower"), &s.row_lower);
    chk.f64s(&k("row_upper"), &r.f64s("row_upper"), &s.row_upper);
    chk.ints(&k("row_deleted"), &r.ints("row_deleted"), &s.row_deleted);
    chk.ints(&k("row_start"), &r.ints("row_start"), &s.row_start);
    chk.ints(&k("row_index"), &r.ints("row_index"), &s.row_index);
    chk.f64s(&k("row_value"), &r.f64s("row_value"), &s.row_value);
    chk.f64s(&k("sum_lower"), &r.f64s("sum_lower"), &s.sum_lower);
    chk.f64s(&k("sum_upper"), &r.f64s("sum_upper"), &s.sum_upper);
    chk.ints(
        &k("num_inf_sum_lower"),
        &r.ints("num_inf_sum_lower"),
        &s.num_inf_sum_lower,
    );
    chk.ints(
        &k("num_inf_sum_upper"),
        &r.ints("num_inf_sum_upper"),
        &s.num_inf_sum_upper,
    );
}

fn replay(fx: &Fixture) -> Option<Checker> {
    let lp = g1_lp(fx);
    let r = fx.first("clqn_result");
    let mut chk = Checker::new(fx.file.clone());
    check(&mut chk, "flag", &r, &normalise(&lp, false));
    check(&mut chk, "flag+coef_change", &r, &normalise(&lp, true));
    Some(chk)
}

#[test]
#[ignore = "1.16.0 target-only fixtures (ORACLE_HIGHS=1.16.0 tools/highs-oracle/gen_fixtures.sh): RUSOLVE_FIXTURES_DIR=target/highs-fixtures-1.16.0 cargo test -p rusolve-ekk --release --test t_g1_clique_normalise -- --ignored"]
fn g1_clique_normalise_matches_highs() {
    let fxs = g1_fixtures("clique_normalise");
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || replay(fx));
    }
    fails.finish("t_g1_clique_normalise");
}
