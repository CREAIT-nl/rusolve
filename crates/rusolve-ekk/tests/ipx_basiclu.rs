//! Golden test for the BASICLU port (`src/ipm/basiclu/`).
//!
//! Fixtures: `basiclu/<lp>.jsonl` (sc50a committed; afiro, kb2, adlittle,
//! blend, share2b target-only), written by `tools/highs-oracle/src/basiclu.cpp`
//! from the real BASICLU of HiGHS 1.15.1 through its object API, on bases of
//! the scaled LP HEkk factorises (slack, mid-simplex, optimal, rank deficient;
//! the last two also with IPX's strict absolute pivot tolerance 1e-3; and a
//! long m-1 update run from the slack basis that reaches the L/U compression
//! of `lu_update`), plus
//! `basiclu/{dense70,rand120,neardep40}.jsonl` (committed): synthetic
//! matrices that reach `lu_pivot_any` (pivot columns with more than 64
//! off-diagonals), many file expansions / garbage collections, and (near
//! dependent columns under the strict tolerance) `lu_remove_col`.
//!
//! Replay, record by record and bit for bit:
//! * `factorize`: `basiclu_obj_initialize(m)` + `basiclu_obj_factorize(B)`;
//!   status, memory sizes after reallocation, the whole `xstore[0..512]`
//!   scalar block, `basiclu_obj_get_factors` (permutations, L, U) and the
//!   solve batch (dense `'N'`/`'T'`, sparse `'N'`/`'T'` with the pattern
//!   order of `ilhs`).
//! * `update_step` (Forrest-Tomlin): `solve_for_update` `'N'` of the
//!   entering column and `'T'` of the leaving position (both with the
//!   solution), `basiclu_obj_update(xtbl)` with the port's own `lhs[p]`,
//!   then the scalar block and the solve batch.
//! * `refactorize`: `basiclu_obj_factorize` of the updated basis on the
//!   used store, then as `factorize`.

mod common;
mod port;

use common::*;
use port::*;
use rusolve_ekk::ipm::basiclu::LuInt;
use rusolve_ekk::ipm::basiclu::basiclu_object::*;
use rusolve_ekk::ipm::basiclu::{
    BASICLU_ABS_PIVOT_TOLERANCE, BASICLU_LNZ, BASICLU_MEMORYL, BASICLU_MEMORYU, BASICLU_MEMORYW,
};
use rusolve_ekk::ipm::basiclu::{BASICLU_OK, BASICLU_UNZ};

/// `xstore[0..512]` as (non-zero bit positions, values).
fn check_xhead(chk: &mut Checker, field: &str, e: &Rec, xstore: &[f64]) {
    check_nz(chk, field, e, &xstore[..512]);
}

/// A vector stored by the oracle as its non-zero bit positions `i` and
/// values `v`; every other entry must be `+0.0`.
fn check_nz(chk: &mut Checker, field: &str, e: &Rec, x: &[f64]) {
    let (ei, ev) = (e.ints("i"), e.f64s("v"));
    let mut ai = Vec::new();
    let mut av = Vec::new();
    for (k, &v) in x.iter().enumerate() {
        if v.to_bits() != 0 {
            ai.push(k as i32);
            av.push(v);
        }
    }
    chk.ints(&format!("{field}.i"), &ei, &ai);
    if ei == ai {
        chk.f64s(&format!("{field}.v"), &ev, &av);
    }
}

fn check_mem(chk: &mut Checker, rec: &Rec, obj: &BasicluObject) {
    chk.int(
        "lmem",
        rec.int("lmem"),
        obj.xstore[BASICLU_MEMORYL] as LuInt,
    );
    chk.int(
        "umem",
        rec.int("umem"),
        obj.xstore[BASICLU_MEMORYU] as LuInt,
    );
    chk.int(
        "wmem",
        rec.int("wmem"),
        obj.xstore[BASICLU_MEMORYW] as LuInt,
    );
}

fn factorize(obj: &mut BasicluObject, b: &Rec) -> LuInt {
    basiclu_obj_factorize(
        obj,
        &b.ints("bbegin"),
        &b.ints("bend"),
        &b.ints("bi"),
        &b.f64s("bx"),
    )
}

fn check_factors(chk: &mut Checker, e: &Rec, obj: &mut BasicluObject, m: usize) {
    let lnz = obj.xstore[BASICLU_LNZ] as LuInt as usize;
    let unz = obj.xstore[BASICLU_UNZ] as LuInt as usize;
    let mut rowperm = vec![0; m];
    let mut colperm = vec![0; m];
    let mut lcolptr = vec![0; m + 1];
    let mut lrowidx = vec![0; m + lnz];
    let mut lvalue = vec![0.0; m + lnz];
    let mut ucolptr = vec![0; m + 1];
    let mut urowidx = vec![0; m + unz];
    let mut uvalue = vec![0.0; m + unz];
    let status = basiclu_obj_get_factors(
        obj,
        Some(&mut rowperm),
        Some(&mut colperm),
        Some(&mut lcolptr),
        Some(&mut lrowidx),
        Some(&mut lvalue),
        Some(&mut ucolptr),
        Some(&mut urowidx),
        Some(&mut uvalue),
    );
    chk.int("factors.status", e.int("status"), status);
    chk.ints("factors.rowperm", &e.ints("rowperm"), &rowperm);
    chk.ints("factors.colperm", &e.ints("colperm"), &colperm);
    chk.ints("factors.lcolptr", &e.ints("lcolptr"), &lcolptr);
    chk.ints("factors.lrowidx", &e.ints("lrowidx"), &lrowidx);
    chk.f64s("factors.lvalue", &e.f64s("lvalue"), &lvalue);
    chk.ints("factors.ucolptr", &e.ints("ucolptr"), &ucolptr);
    chk.ints("factors.urowidx", &e.ints("urowidx"), &urowidx);
    chk.f64s("factors.uvalue", &e.f64s("uvalue"), &uvalue);
    check_xhead(chk, "factors.xhead", &e.obj("xhead"), &obj.xstore);
}

struct Rhs {
    dense: Vec<f64>,
    sparse_i: Vec<LuInt>,
    sparse_v: Vec<f64>,
}

fn check_solves(chk: &mut Checker, solves: &[Rec], obj: &mut BasicluObject, rhs: &Rhs, m: usize) {
    for (k, s) in solves.iter().enumerate() {
        let op = s.str("op");
        let trans = s.str("trans").as_bytes()[0];
        let f = format!("solves[{k}]({op} {})", trans as char);
        match op {
            "dense" => {
                let mut lhs = vec![0.0; m];
                let status = basiclu_obj_solve_dense(obj, &rhs.dense, &mut lhs, trans);
                chk.int(&format!("{f}.status"), s.int("status"), status);
                chk.f64s(&format!("{f}.lhs"), &s.f64s("lhs"), &lhs);
            }
            "sparse" => {
                let status = basiclu_obj_solve_sparse(
                    obj,
                    rhs.sparse_i.len() as LuInt,
                    &rhs.sparse_i,
                    &rhs.sparse_v,
                    trans,
                );
                chk.int(&format!("{f}.status"), s.int("status"), status);
                chk.int(&format!("{f}.nzlhs"), s.int("nzlhs"), obj.nzlhs);
                let nz = obj.nzlhs.max(0) as usize;
                chk.ints(&format!("{f}.ilhs"), &s.ints("ilhs"), &obj.ilhs[..nz]);
                check_nz(chk, &format!("{f}.lhs"), &s.obj("lhs"), &obj.lhs[..m]);
            }
            _ => panic!("unknown solve op {op:?}"),
        }
        if s.has("xhead") {
            check_xhead(chk, &format!("{f}.xhead"), &s.obj("xhead"), &obj.xstore);
        }
    }
}

/// Replays one basis: its `factorize` record, the `update_step` records
/// that follow it and the closing `refactorize` record.
fn replay(fails: &mut Failures, fx: &Fixture, start: usize) -> usize {
    let rec = fx.rec(start);
    let m = rec.int("m") as usize;
    let lp = rec.str("lp").to_string();
    let basis = rec.str("basis").to_string();
    let rhs = Rhs {
        dense: rec.f64s("rhs_dense"),
        sparse_i: rec.ints("rhs_sparse_i"),
        sparse_v: rec.f64s("rhs_sparse_v"),
    };
    let mut obj = BasicluObject::default();
    assert_eq!(basiclu_obj_initialize(&mut obj, m as LuInt), BASICLU_OK);
    obj.xstore[BASICLU_ABS_PIVOT_TOLERANCE] = rec.f64("abstol");

    let mut chk = Checker::new(format!("{} {lp} basis={basis} factorize", rec.ctx));
    let status = factorize(&mut obj, &rec.obj("b"));
    chk.int("status", rec.int("status"), status);
    check_mem(&mut chk, &rec, &obj);
    check_xhead(&mut chk, "xhead", &rec.obj("xhead"), &obj.xstore);
    check_factors(&mut chk, &rec.obj("factors"), &mut obj, m);
    check_solves(&mut chk, &rec.list("solves"), &mut obj, &rhs, m);
    let ok = fails.take(&chk);
    let mut next = start + 1;
    if !ok {
        // Later records depend on this state; skip them.
        while next < fx.len() && fx.rec(next).t() != "factorize" {
            next += 1;
        }
        return next;
    }

    while next < fx.len() {
        let r = fx.rec(next);
        match r.t() {
            "update_step" => {
                let mut chk = Checker::new(format!(
                    "{} {lp} basis={basis} update_step k={}",
                    r.ctx,
                    r.int("k")
                ));
                let p = r.int("pos_out");
                let col_i = r.ints("col_i");
                let col_v = r.f64s("col_v");
                let st = basiclu_obj_solve_for_update(
                    &mut obj,
                    col_i.len() as LuInt,
                    &col_i,
                    Some(&col_v),
                    b'N',
                    1,
                );
                chk.int("ftran_status", r.int("ftran_status"), st);
                chk.int("ftran_nzlhs", r.int("ftran_nzlhs"), obj.nzlhs);
                let nz = obj.nzlhs.max(0) as usize;
                chk.ints("ftran_ilhs", &r.ints("ftran_ilhs"), &obj.ilhs[..nz]);
                check_nz(&mut chk, "ftran_lhs", &r.obj("ftran_lhs"), &obj.lhs[..m]);
                let xtbl = obj.lhs[p as usize];
                let st = basiclu_obj_solve_for_update(&mut obj, 0, &[p], None, b'T', 1);
                chk.int("btran_status", r.int("btran_status"), st);
                chk.int("btran_nzlhs", r.int("btran_nzlhs"), obj.nzlhs);
                let nz = obj.nzlhs.max(0) as usize;
                chk.ints("btran_ilhs", &r.ints("btran_ilhs"), &obj.ilhs[..nz]);
                check_nz(&mut chk, "btran_lhs", &r.obj("btran_lhs"), &obj.lhs[..m]);
                chk.f64("xtbl", r.f64("xtbl"), xtbl);
                let st = basiclu_obj_update(&mut obj, xtbl);
                chk.int("update_status", r.int("update_status"), st);
                check_mem(&mut chk, &r, &obj);
                check_xhead(&mut chk, "xhead", &r.obj("xhead"), &obj.xstore);
                if r.has("solves") {
                    check_solves(&mut chk, &r.list("solves"), &mut obj, &rhs, m);
                }
                next += 1;
                if !fails.take(&chk) {
                    while next < fx.len() && fx.rec(next).t() != "factorize" {
                        next += 1;
                    }
                    return next;
                }
            }
            "refactorize" => {
                let mut chk = Checker::new(format!("{} {lp} basis={basis} refactorize", r.ctx));
                let status = factorize(&mut obj, &r.obj("b"));
                chk.int("status", r.int("status"), status);
                check_mem(&mut chk, &r, &obj);
                check_xhead(&mut chk, "xhead", &r.obj("xhead"), &obj.xstore);
                check_factors(&mut chk, &r.obj("factors"), &mut obj, m);
                check_solves(&mut chk, &r.list("solves"), &mut obj, &rhs, m);
                fails.take(&chk);
                return next + 1;
            }
            _ => return next,
        }
    }
    next
}

#[test]
fn basiclu_matches_highs() {
    let mut fails = Failures::default();
    let mut num_updates = 0;
    for fx in fixtures_in("basiclu") {
        num_updates += fx.of_type("update_step").count();
        let mut k = 0;
        while k < fx.len() {
            if fx.rec(k).t() == "factorize" {
                let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    replay(&mut fails, &fx, k)
                }));
                match r {
                    Ok(next) => k = next,
                    Err(e) => {
                        let msg = e
                            .downcast_ref::<String>()
                            .cloned()
                            .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                            .unwrap_or_default();
                        let ctx = fx.rec(k).ctx.clone();
                        fails.fail(format!("{ctx}: PANIC inside the port: {msg}"));
                        k += 1;
                        while k < fx.len() && fx.rec(k).t() != "factorize" {
                            k += 1;
                        }
                    }
                }
            } else {
                k += 1;
            }
        }
    }
    assert!(num_updates > 0, "no update_step records");
    fails.finish("ipx_basiclu");
}
