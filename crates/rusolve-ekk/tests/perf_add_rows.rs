//! Timing probe: adding rows one at a time to a wide
//! model must cost O(nnz of the row), not O(num_col) (HighsMatrixUtils.cpp:153 uses a hash map).
//! Run: cargo test -p rusolve-ekk --release --test perf_add_rows -- --ignored --nocapture

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use std::time::Instant;

fn quiet() -> Highs {
    let mut h = Highs::default();
    h.set_option_value_bool("output_flag", false);
    h
}

#[test]
#[ignore]
fn add_row_one_at_a_time_200k_cols() {
    let ncol: HighsInt = 200_000;
    let nrow: HighsInt = 100_000;
    let mut h = quiet();
    let cost = vec![0.0; ncol as usize];
    let lo = vec![0.0; ncol as usize];
    let up = vec![1.0; ncol as usize];
    h.add_cols(ncol, &cost, &lo, &up, 0, &[], &[], &[]);
    let t = Instant::now();
    for r in 0..nrow {
        let idx = [r, r + ncol / 2];
        h.add_row(0.0, 1.0, 2, &idx, &[1.0, 1.0]);
    }
    let dt = t.elapsed().as_secs_f64();
    println!(
        "ADD_ROW_LOOP {nrow} rows: {dt:.3} s ({:.2} us/row)",
        dt * 1e6 / nrow as f64
    );
}

#[test]
#[ignore]
#[allow(clippy::field_reassign_with_default)]
fn pass_model_3m_nnz() {
    let (nrow, ncol, per_col) = (100_000usize, 100_000usize, 30usize);
    let mut lp = HighsLp::default();
    lp.num_col_ = ncol as HighsInt;
    lp.num_row_ = nrow as HighsInt;
    lp.col_cost_ = vec![1.0; ncol];
    lp.col_lower_ = vec![0.0; ncol];
    lp.col_upper_ = vec![1.0; ncol];
    lp.row_lower_ = vec![0.0; nrow];
    lp.row_upper_ = vec![1e3; nrow];
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = ncol as HighsInt;
    lp.a_matrix_.num_row_ = nrow as HighsInt;
    let mut s = 1u64;
    let (mut start, mut index, mut value) = (vec![0 as HighsInt], Vec::new(), Vec::new());
    for _ in 0..ncol {
        let mut rows: Vec<HighsInt> = Vec::new();
        while rows.len() < per_col {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let r = ((s >> 33) as usize % nrow) as HighsInt;
            if !rows.contains(&r) {
                rows.push(r);
            }
        }
        rows.sort();
        for r in rows {
            index.push(r);
            value.push(1.0);
        }
        start.push(index.len() as HighsInt);
    }
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    let mut h = quiet();
    let t = Instant::now();
    h.pass_model(lp);
    println!("PASS_MODEL 3M nnz: {:.3} s", t.elapsed().as_secs_f64());
}
