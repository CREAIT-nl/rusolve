//! Golden test for the IPX model port (`src/ipm/ipx/model/`).
//!
//! Fixtures: `ipx/model/<lp>.<cfg>.jsonl`, written by `highs-oracle
//! ipx-model` (tools/highs-oracle/src/ipx.cpp) from the real `ipx::Model` of
//! HiGHS 1.15.1: `afiro`, `ipxmix` (free, boxed, fixed and flipped
//! variables, a ranged row, power-of-two equilibration) and `ipxtall`
//! (dualized by the default rule), each also with `dualize=1` (`.dual`),
//! committed; the netlib set target-only.
//!
//! Replay, bit for bit:
//! * `model`: `Model::Load` of the `input` arrays (IpxWrapper.cpp
//!   fillInIpxData) under the `params` record: errflag, the dualization
//!   decision, AI / AIt / A, b, c, bounds, norms, the scaled user data,
//!   flipped / boxed variables, colscale / rowscale, `GetInfo`.
//! * `presolve_starting_point`, `presolve_ipm_starting_point`,
//!   `postsolve_interior` (+ `EvaluateInteriorSolution`): the user <->
//!   solver space transforms of the fixed vectors in the record.

mod common;
mod port;

use std::rc::Rc;

use common::*;
use port::ipx::*;
use port::*;
use rusolve_ekk::ipm::ipx::model::Model;
use rusolve_ekk::ipm::ipx::{Info, Vector};

fn check_model(chk: &mut Checker, e: &Rec, model: &Model) {
    let v = model.oracle_view();
    chk.bool("dualized", e.bool("dualized"), v.dualized);
    for (k, a) in [
        ("num_rows", v.num_rows),
        ("num_cols", v.num_cols),
        ("num_dense_cols", v.num_dense_cols),
        ("nz_dense", v.nz_dense),
        ("num_constr", v.num_constr),
        ("num_eqconstr", v.num_eqconstr),
        ("num_var", v.num_var),
        ("num_free_var", v.num_free_var),
        ("num_entries", v.num_entries),
    ] {
        chk.int(k, e.int(k), a);
    }
    check_ipx_matrix(chk, "AI", &e.obj("AI"), v.AI);
    check_ipx_matrix(chk, "AIt", &e.obj("AIt"), v.AIt);
    check_ipx_matrix(chk, "A", &e.obj("A"), v.A);
    for (k, a) in [
        ("b", v.b),
        ("c", v.c),
        ("lb", v.lb),
        ("ub", v.ub),
        ("scaled_obj", v.scaled_obj),
        ("scaled_rhs", v.scaled_rhs),
        ("scaled_lbuser", v.scaled_lbuser),
        ("scaled_ubuser", v.scaled_ubuser),
        ("colscale", v.colscale),
        ("rowscale", v.rowscale),
    ] {
        chk.f64s(k, &e.f64s(k), a);
    }
    for (k, a) in [
        ("norm_bounds", v.norm_bounds),
        ("norm_c", v.norm_c),
        ("norm_obj", v.norm_obj),
        ("norm_rhs", v.norm_rhs),
        ("offset", v.offset),
    ] {
        chk.f64(k, e.f64(k), a);
    }
    chk.ints("boxed_vars", &e.ints("boxed_vars"), v.boxed_vars);
    chk.ints("flipped_vars", &e.ints("flipped_vars"), v.flipped_vars);
    let ct: Vec<i32> = v.constr_type.iter().map(|&c| c as i32).collect();
    chk.ints("constr_type", &e.ints("constr_type"), &ct);
    chk.bool(
        "filippo",
        e.bool("filippo"),
        model.filippo_dualization_test(),
    );
    let mut info = Info::default();
    model.get_info(&mut info);
    check_ipx_info(chk, "info", &e.obj("info"), &info);
}

fn opt(v: &[f64]) -> Option<&[f64]> {
    Some(v)
}

fn replay_transform(chk: &mut Checker, r: &Rec, model: &Model) {
    let (n, m) = (model.cols(), model.rows());
    match r.t() {
        "presolve_starting_point" => {
            let (x, s, y, z) = (
                r.f64s("x_user"),
                r.f64s("slack_user"),
                r.f64s("y_user"),
                r.f64s("z_user"),
            );
            let (mut xs, mut ys, mut zs) = (Vector::new(n + m), Vector::new(m), Vector::new(n + m));
            model.presolve_starting_point(
                opt(&x),
                opt(&s),
                opt(&y),
                opt(&z),
                &mut xs,
                &mut ys,
                &mut zs,
            );
            check_vec(chk, "out", r, "x", &xs);
            check_vec(chk, "out", r, "y", &ys);
            check_vec(chk, "out", r, "z", &zs);
        }
        "presolve_ipm_starting_point" => {
            let u: Vec<Vec<f64>> = [
                "x_user",
                "xl_user",
                "xu_user",
                "slack_user",
                "y_user",
                "zl_user",
                "zu_user",
            ]
            .iter()
            .map(|k| r.f64s(k))
            .collect();
            let mut xs = Vector::new(n + m);
            let mut xls = Vector::new(n + m);
            let mut xus = Vector::new(n + m);
            let mut ys = Vector::new(m);
            let mut zls = Vector::new(n + m);
            let mut zus = Vector::new(n + m);
            let err = model.presolve_ipm_starting_point(
                opt(&u[0]),
                opt(&u[1]),
                opt(&u[2]),
                opt(&u[3]),
                opt(&u[4]),
                opt(&u[5]),
                opt(&u[6]),
                &mut xs,
                &mut xls,
                &mut xus,
                &mut ys,
                &mut zls,
                &mut zus,
            );
            chk.int("errflag", r.int("errflag"), err);
            check_vec(chk, "out", r, "x", &xs);
            check_vec(chk, "out", r, "xl", &xls);
            check_vec(chk, "out", r, "xu", &xus);
            check_vec(chk, "out", r, "y", &ys);
            check_vec(chk, "out", r, "zl", &zls);
            check_vec(chk, "out", r, "zu", &zus);
        }
        "postsolve_interior" => {
            let x = vector_of(r, "x_solver");
            let xl = vector_of(r, "xl_solver");
            let xu = vector_of(r, "xu_solver");
            let y = vector_of(r, "y_solver");
            let zl = vector_of(r, "zl_solver");
            let zu = vector_of(r, "zu_solver");
            let (nu, mu) = (r.f64s("x").len(), r.f64s("slack").len());
            let (mut xu_, mut xlu, mut xuu, mut zlu, mut zuu) = (
                vec![0.0; nu],
                vec![0.0; nu],
                vec![0.0; nu],
                vec![0.0; nu],
                vec![0.0; nu],
            );
            let (mut su, mut yu) = (vec![0.0; mu], vec![0.0; mu]);
            model.postsolve_interior_solution(
                &x,
                &xl,
                &xu,
                &y,
                &zl,
                &zu,
                Some(&mut xu_),
                Some(&mut xlu),
                Some(&mut xuu),
                Some(&mut su),
                Some(&mut yu),
                Some(&mut zlu),
                Some(&mut zuu),
            );
            check_vec(chk, "out", r, "x", &xu_);
            check_vec(chk, "out", r, "xl", &xlu);
            check_vec(chk, "out", r, "xu", &xuu);
            check_vec(chk, "out", r, "slack", &su);
            check_vec(chk, "out", r, "y", &yu);
            check_vec(chk, "out", r, "zl", &zlu);
            check_vec(chk, "out", r, "zu", &zuu);
            let mut info = Info::default();
            model.evaluate_interior_solution(&x, &xl, &xu, &y, &zl, &zu, &mut info);
            check_ipx_info(chk, "info", &r.obj("info"), &info);
        }
        t => chk.push("t", "known record".into(), t.into(), String::new()),
    }
}

/// `dualize=1` trips `assert(dualize == 0 || dualize == -1)`-style debug
/// assertions the NDEBUG oracle compiles out: those fixtures replay in release builds only.
pub fn skip_in_debug(fx: &Fixture) -> bool {
    cfg!(debug_assertions) && fx.first("params").int("dualize") > 0
}

#[test]
fn ipx_model_matches_highs() {
    let mut fails = Failures::default();
    let mut dualized = 0;
    let mut flipped = 0;
    for fx in fixtures_in("ipx/model") {
        if skip_in_debug(&fx) {
            eprintln!("SKIP {} (dualize=1 in a debug build)", fx.file);
            continue;
        }
        let control = ipx_control(&fx);
        let input = fx.first("input");
        let e = fx.first("model");
        let mut model = Model::new();
        let err = ipx_load_model(&control, &input, &mut model);
        let model = Rc::new(model);
        fails.run(&e.ctx, || {
            let mut chk = Checker::new(e.ctx.clone());
            chk.int("errflag", e.int("errflag"), err);
            check_model(&mut chk, &e, &model);
            Some(chk)
        });
        dualized += e.bool("dualized") as usize;
        flipped += !e.ints("flipped_vars").is_empty() as usize;
        if err != 0 {
            continue;
        }
        for r in fx.recs().skip(3) {
            fails.run(&r.ctx.clone(), || {
                let mut chk = Checker::new(r.ctx.clone());
                replay_transform(&mut chk, &r, &model);
                Some(chk)
            });
        }
    }
    assert!(dualized > 0, "no dualized model in the fixtures");
    assert!(
        flipped > 0,
        "no model with flipped variables in the fixtures"
    );
    fails.finish("ipx_model");
}
