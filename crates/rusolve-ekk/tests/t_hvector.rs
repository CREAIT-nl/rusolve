//! Golden tests for manifest row 4, `src/util/hvector.rs` (HVectorBase.cpp).
//!
//! Fixtures: the `hvector_*` records of `unit/<lp>.jsonl` (afiro, kb2, sc50a
//! committed; adlittle, blend, share2b target-only), written by
//! `tools/highs-oracle/src/unit.cpp`. The oracle's calls, replayed here on
//! an HVector rebuilt from the record's `in` (`setup(size)`, `clear()`, then
//! count / index / array / synthetic_tick):
//!
//! - `hvector_tight`: `t = in; t.tight()` (also with `count = -1`).
//! - `hvector_pack`: `p = in; p.packFlag = true; p.pack()`; `out` carries the
//!   pack fields.
//! - `hvector_norm2`: `in.norm2()`.
//! - `hvector_reIndex`: `r = in; r.reIndex()` (also with `count = -1`).
//! - `hvector_saxpy`: `s = in; s.saxpy(pivotX, &pivot)`.
//!
//! `setup` and `clear` are covered by every other test that builds vectors
//! (`hvector_from_rec`), and by `setup_and_clear_reset_state` below.

mod common;
mod port;

use common::*;
use port::*;
use rusolve_ekk::prelude::*;

fn ctx(rec: &Rec, what: &str) -> String {
    format!("{} (HVector::{what})", rec.ctx)
}

#[test]
fn tight_matches_highs() {
    for fx in fixtures_in("unit") {
        let mut n = 0;
        for rec in recs_of(&fx, "hvector_tight") {
            let input = rec.hv("in");
            let mut v = hvector_from_rec(&input);
            v.tight();
            let mut chk = Checker::new(format!("{} in.count={}", ctx(&rec, "tight"), input.count));
            check_hv(&mut chk, "out", &rec.hv("out"), &v);
            chk.finish();
            n += 1;
        }
        assert_eq!(n, 2, "{}: expected 2 hvector_tight records", fx.file);
    }
}

#[test]
fn pack_matches_highs() {
    for fx in fixtures_in("unit") {
        for rec in recs_of(&fx, "hvector_pack") {
            let mut v = hvector_from_rec(&rec.hv("in"));
            v.pack_flag = true;
            v.pack();
            let mut chk = Checker::new(ctx(&rec, "pack"));
            check_hv(&mut chk, "out", &rec.hv("out"), &v);
            chk.finish();
        }
    }
}

#[test]
fn norm2_matches_highs() {
    for fx in fixtures_in("unit") {
        for rec in recs_of(&fx, "hvector_norm2") {
            let v = hvector_from_rec(&rec.hv("in"));
            let mut chk = Checker::new(ctx(&rec, "norm2"));
            chk.f64("out", rec.f64("out"), v.norm2());
            chk.finish();
        }
    }
}

#[test]
fn re_index_matches_highs() {
    for fx in fixtures_in("unit") {
        let mut n = 0;
        for rec in recs_of(&fx, "hvector_reIndex") {
            let input = rec.hv("in");
            let mut v = hvector_from_rec(&input);
            v.re_index();
            let mut chk =
                Checker::new(format!("{} in.count={}", ctx(&rec, "reIndex"), input.count));
            check_hv(&mut chk, "out", &rec.hv("out"), &v);
            chk.finish();
            n += 1;
        }
        assert_eq!(n, 2, "{}: expected 2 hvector_reIndex records", fx.file);
    }
}

#[test]
fn saxpy_matches_highs() {
    for fx in fixtures_in("unit") {
        for (k, rec) in recs_of(&fx, "hvector_saxpy").into_iter().enumerate() {
            let pivot_x = rec.f64("pivotX");
            let pivot = hvector_from_rec(&rec.hv("pivot"));
            let mut v = hvector_from_rec(&rec.hv("in"));
            v.saxpy(pivot_x, &pivot);
            let mut chk = Checker::new(format!(
                "{} case {k} pivotX={} pivot.count={}",
                ctx(&rec, "saxpy"),
                f64_hex(pivot_x),
                pivot.count
            ));
            check_hv(&mut chk, "out", &rec.hv("out"), &v);
            chk.finish();
        }
    }
}

/// `setup(size)` sizes every array to `size`, zeroes `array`, sets
/// `count = 0`, `packFlag = false`, `synthetic_tick = 0`; `clear()` then
/// leaves a zero vector (HVectorBase.cpp `setup`, `clear`, `clearScalars`).
/// Not a fixture comparison: the observable postconditions of the C++ code.
#[test]
fn setup_and_clear_reset_state() {
    let mut v = HVector::default();
    v.setup(5);
    let ctx = "HVector::setup(5)";
    assert_eq!(v.size, 5, "{ctx}: size");
    assert_eq!(v.count, 0, "{ctx}: count");
    assert_eq!(v.array.len(), 5, "{ctx}: array.len()");
    assert_eq!(v.index.len(), 5, "{ctx}: index.len()");
    assert!(!v.pack_flag, "{ctx}: packFlag");
    assert_eq!(v.synthetic_tick.to_bits(), 0, "{ctx}: synthetic_tick");
    v.index[0] = 3;
    v.array[3] = -2.5;
    v.count = 1;
    v.synthetic_tick = 7.0;
    v.clear();
    let ctx = "HVector::clear() after one entry";
    assert_eq!(v.count, 0, "{ctx}: count");
    assert!(
        v.array.iter().all(|x| x.to_bits() == 0),
        "{ctx}: array must be all +0.0, got {:?}",
        v.array
    );
    assert_eq!(v.synthetic_tick.to_bits(), 0, "{ctx}: synthetic_tick");
    assert!(!v.pack_flag, "{ctx}: packFlag");
}
