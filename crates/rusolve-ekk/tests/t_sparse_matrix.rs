//! Golden tests for manifest rows 13-16, `src/util/sparse_matrix/*`
//! (HighsSparseMatrix.cpp), with HighsCDouble (row 1) and
//! HighsSparseVectorSum (row 8) underneath the `quad = 1` products.
//!
//! Fixtures: `unit/<lp>.jsonl` (afiro, kb2, sc50a committed; adlittle,
//! blend, share2b target-only), written by `tools/highs-oracle/src/unit.cpp`
//! on the *unscaled* LP matrix as read from MPS. Every record carries its
//! inputs, and each test replays the oracle's exact call:
//!
//! - `createRowwise`: `ar.createRowwise(a)`.
//! - `createRowwisePartitioned`: `arp.createRowwisePartitioned(a, NULL)` and
//!   with the mixed partition `part[j] = (j%3==1)?0:1` (columns),
//!   `part[num_col+i] = (i%3==1)?1:0` (rows). The record stores only
//!   `part[0..num_col)`; the row part is rebuilt from the formula and the
//!   column part is checked against the record.
//! - `priceByColumn(quad, res, col)` on `a`; `priceByRow(quad, res, col)` and
//!   `priceByRowWithSwitch(quad, res, col, ed, from, sd)` on the row-wise or
//!   partitioned matrix (`matrix` field); `res` is `setup(num_col)+clear()`.
//! - `computeDot(col, j)` for each `use_col`.
//! - `collectAj(v, col, multiplier)` accumulated into ONE vector over the
//!   `steps` (so step k depends on steps < k: read the first failing step).
//! - `update(var_in, var_out, a)` on the partitioned matrix, from `before`,
//!   over the `steps` (sequential: first failing step only).

mod common;
mod port;

use common::*;
use port::*;
use rusolve_ekk::prelude::*;

fn lp_matrix(fx: &Fixture) -> (HighsSparseMatrix, String) {
    let lp = fx.first("lp");
    (matrix_from_rec(&lp.mat("a")), lp.str("name").to_string())
}

/// The mixed partition of unit.cpp, full length num_col + num_row.
fn mixed_partition(num_col: i32, num_row: i32) -> Vec<i8> {
    let mut part = vec![1i8; (num_col + num_row) as usize];
    for j in 0..num_col {
        part[j as usize] = if j % 3 == 1 { 0 } else { 1 };
    }
    for i in 0..num_row {
        part[(num_col + i) as usize] = if i % 3 == 1 { 1 } else { 0 };
    }
    part
}

fn rowwise(a: &HighsSparseMatrix) -> HighsSparseMatrix {
    let mut ar = HighsSparseMatrix::default();
    ar.create_rowwise(a);
    ar
}

fn partitioned(a: &HighsSparseMatrix, part: Option<&[i8]>) -> HighsSparseMatrix {
    let mut arp = HighsSparseMatrix::default();
    arp.create_rowwise_partitioned(a, part);
    arp
}

fn new_vec(size: i32) -> HVector {
    let mut v = HVector::default();
    v.setup(size);
    v.clear();
    v
}

#[test]
fn create_rowwise_matches_highs() {
    for fx in fixtures_in("unit") {
        let (a, name) = lp_matrix(&fx);
        for rec in recs_of(&fx, "createRowwise") {
            let ar = rowwise(&a);
            let mut chk = Checker::new(format!(
                "{} {name} (HighsSparseMatrix::createRowwise)",
                rec.ctx
            ));
            chk.mat("out", &rec.mat("out"), &matrix_view(&ar));
            chk.finish();
        }
    }
}

#[test]
fn create_rowwise_partitioned_matches_highs() {
    for fx in fixtures_in("unit") {
        let (a, name) = lp_matrix(&fx);
        let mut n = 0;
        for rec in recs_of(&fx, "createRowwisePartitioned") {
            let null_partition = rec.bool("null_partition");
            let part = mixed_partition(a.num_col_, a.num_row_);
            let mut chk = Checker::new(format!(
                "{} {name} null_partition={null_partition} (HighsSparseMatrix::createRowwisePartitioned)",
                rec.ctx
            ));
            let arp = if null_partition {
                partitioned(&a, None)
            } else {
                chk.i8s(
                    "in_partition[0..num_col) (test's partition formula)",
                    &rec.i8s("in_partition"),
                    &part[..a.num_col_ as usize],
                );
                partitioned(&a, Some(&part))
            };
            chk.mat("out", &rec.mat("out"), &matrix_view(&arp));
            chk.finish();
            n += 1;
        }
        assert_eq!(
            n, 2,
            "{}: expected 2 createRowwisePartitioned records",
            fx.file
        );
    }
}

#[test]
fn price_by_column_matches_highs() {
    for fx in fixtures_in("unit") {
        let (a, name) = lp_matrix(&fx);
        for rec in recs_of(&fx, "priceByColumn") {
            let quad = rec.int("quad") != 0;
            let col = hvector_from_rec(&rec.hv("in"));
            let mut res = new_vec(a.num_col_);
            a.price_by_column(quad, &mut res, &col, DEBUG_REPORT_OFF);
            let mut chk = Checker::new(format!(
                "{} {name} in_kind={} quad={} (HighsSparseMatrix::priceByColumn)",
                rec.ctx,
                rec.str("in_kind"),
                quad as i32
            ));
            check_hv(&mut chk, "out", &rec.hv("out"), &res);
            chk.finish();
        }
    }
}

#[test]
fn price_by_row_matches_highs() {
    for fx in fixtures_in("unit") {
        let (a, name) = lp_matrix(&fx);
        let ar = rowwise(&a);
        let arp = partitioned(&a, Some(&mixed_partition(a.num_col_, a.num_row_)));
        for rec in recs_of(&fx, "priceByRow") {
            let quad = rec.int("quad") != 0;
            let which = rec.str("matrix");
            let m = if which == "partitioned" { &arp } else { &ar };
            let col = hvector_from_rec(&rec.hv("in"));
            let mut res = new_vec(a.num_col_);
            m.price_by_row(quad, &mut res, &col, DEBUG_REPORT_OFF);
            let mut chk = Checker::new(format!(
                "{} {name} matrix={which} in_kind={} quad={} (HighsSparseMatrix::priceByRow)",
                rec.ctx,
                rec.str("in_kind"),
                quad as i32
            ));
            check_hv(&mut chk, "out", &rec.hv("out"), &res);
            chk.finish();
        }
    }
}

#[test]
fn price_by_row_with_switch_matches_highs() {
    for fx in fixtures_in("unit") {
        let (a, name) = lp_matrix(&fx);
        let ar = rowwise(&a);
        let arp = partitioned(&a, Some(&mixed_partition(a.num_col_, a.num_row_)));
        for rec in recs_of(&fx, "priceByRowWithSwitch") {
            let quad = rec.int("quad") != 0;
            let which = rec.str("matrix");
            let m = if which == "partitioned" { &arp } else { &ar };
            let ed = rec.f64("expected_density");
            let from = rec.int("from_index");
            let sd = rec.f64("switch_density");
            let col = hvector_from_rec(&rec.hv("in"));
            let mut res = new_vec(a.num_col_);
            m.price_by_row_with_switch(quad, &mut res, &col, ed, from, sd, DEBUG_REPORT_OFF);
            let mut chk = Checker::new(format!(
                "{} {name} matrix={which} in_kind={} quad={} expected_density={ed} from_index={from} switch_density={sd} (HighsSparseMatrix::priceByRowWithSwitch)",
                rec.ctx,
                rec.str("in_kind"),
                quad as i32
            ));
            check_hv(&mut chk, "out", &rec.hv("out"), &res);
            chk.finish();
        }
    }
}

#[test]
fn compute_dot_matches_highs() {
    for fx in fixtures_in("unit") {
        let (a, name) = lp_matrix(&fx);
        for rec in recs_of(&fx, "computeDot") {
            let col = hvector_from_rec(&rec.hv("in"));
            let use_col = rec.ints("use_col");
            let got: Vec<f64> = use_col
                .iter()
                .map(|&j| a.compute_dot_hvector(&col, j))
                .collect();
            let mut chk = Checker::new(format!(
                "{} {name} in_kind={} (HighsSparseMatrix::computeDot; out[k] is use_col[k], use_col={use_col:?})",
                rec.ctx,
                rec.str("in_kind")
            ));
            chk.f64s("out", &rec.f64s("out"), &got);
            chk.finish();
        }
    }
}

#[test]
fn collect_aj_sequence_matches_highs() {
    for fx in fixtures_in("unit") {
        let (a, name) = lp_matrix(&fx);
        for rec in recs_of(&fx, "collectAj") {
            let mut v = new_vec(rec.int("size"));
            for (k, st) in rec.list("steps").iter().enumerate() {
                let col = st.int("col");
                let mul = st.f64("multiplier");
                a.collect_aj(&mut v, col, mul);
                let mut chk = Checker::new(format!(
                    "{} {name} step {k} (col={col}, multiplier={}) after steps 0..{k} accumulated into one HVector (HighsSparseMatrix::collectAj)",
                    rec.ctx,
                    f64_hex(mul)
                ));
                check_hv(&mut chk, "after", &st.hv("after"), &v);
                chk.finish();
            }
        }
    }
}

#[test]
fn update_partitioned_sequence_matches_highs() {
    for fx in fixtures_in("unit") {
        let (a, name) = lp_matrix(&fx);
        for rec in recs_of(&fx, "updatePartitioned") {
            let part = mixed_partition(a.num_col_, a.num_row_);
            let before = rec.mat("before");
            let mut chk = Checker::new(format!(
                "{} {name} (inputs of the HighsSparseMatrix::update sequence)",
                rec.ctx
            ));
            chk.i8s(
                "in_partition[0..num_col) (test's partition formula)",
                &rec.i8s("in_partition"),
                &part[..a.num_col_ as usize],
            );
            chk.finish();
            let mut m = matrix_from_rec(&before);
            for (k, st) in rec.list("steps").iter().enumerate() {
                let var_in = st.int("var_in");
                let var_out = st.int("var_out");
                m.update(var_in, var_out, &a);
                let mut chk = Checker::new(format!(
                    "{} {name} step {k} (var_in={var_in}, var_out={var_out}), applied after steps 0..{k} starting from `before` (HighsSparseMatrix::update)",
                    rec.ctx
                ));
                chk.mat("after", &st.mat("after"), &matrix_view(&m));
                chk.finish();
            }
        }
    }
}
