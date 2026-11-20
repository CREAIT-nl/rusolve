//! `rusolve ekk-seq <model.mps>`: the Rust half of the rusolve-ekk
//! model-modification acceptance harness. It replays, through the ported `Highs` API, the exact
//! call sequence of tools/highs-oracle/api_seq.cpp (stages 0-3, the column
//! and row choice rules, `presolve = "off"` after stage 0) and prints the
//! same lines; the two outputs must be identical.
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//! (Only the API call sequence; the driver itself is ours.)

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use std::process::ExitCode;

/// Reads `path` as `ekk::solve_ekk_cli` does (`read_mps_file_file_order`,
/// `highs_lp_from_model`; integrality is already dropped there, the
/// `clearIntegrality` of the C++ driver), `output_flag = false`,
/// `pass_model`, then:
/// * stage 0: `run()`; then `options_.presolve = "off"`.
/// * stage 1: column `j` and its new bound exactly as api_seq.cpp; one
///   `change_cols_bounds_set(1, &[j], &[lo], &[up])`; `run()`.
/// * stage 2: the row of api_seq.cpp via `add_rows(1, &[-kHighsInf],
///   &[upper], nnz, &[0], &index, &value)`; `run()` (or print
///   `stage 2 skip` and stop).
/// * stage 3: `delete_rows_mask` of the last row; `run()`.
///
/// Each stage prints `stage <k> status=<s> iters=<n> obj=<o>` with `s` the
/// HiGHS `modelStatusToString` text, `n` = `info_.simplex_iteration_count`
/// and `o` the C `%.17g` text (17 significant digits, trailing zeros
/// stripped, exponent form below 1e-5 / at or above 1e17; Rust `{}` and
/// `{:e}` are NOT that text). Check against api_seq on afiro first.
pub fn run_ekk_seq(path: &str) -> ExitCode {
    // C printf "%.17g".
    fn fmt_g17(v: f64) -> String {
        rusolve_ekk::util::c_printf::c_fmt_g(v, 17)
    }

    // One stage line of api_seq.cpp `report`.
    fn report(h: &Highs, stage: i32) {
        println!(
            "stage {} status={} iters={} obj={}",
            stage,
            util_model_status_to_string(h.get_model_status()),
            h.get_info().simplex_iteration_count,
            fmt_g17(h.get_info().objective_function_value)
        );
    }

    let (model, order) = match rusolve_io::read_mps_file_file_order(path) {
        Ok(x) => x,
        Err(e) => {
            eprintln!("read {path}: {e}");
            return ExitCode::FAILURE;
        }
    };
    if let Err(e) = model.gc.reject("rusolve ekk-seq") {
        eprintln!("{e}");
        return ExitCode::FAILURE;
    }
    let lp = crate::ekk::highs_lp_from_model(&model, &order);
    let mut options = HighsOptions {
        output_flag: false,
        ..crate::base_options()
    };
    options.sync_log_options();
    let mut h = Highs::with_options(options);
    if h.pass_model(lp) == HighsStatus::kError {
        eprintln!("ekk: passModel failed");
        return ExitCode::FAILURE;
    }

    let _ = h.run();
    report(&h, 0);
    h.options_.presolve = "off".to_string();

    let n = h.get_lp().num_col_;
    if n == 0 {
        return ExitCode::SUCCESS;
    }
    let mut x: Vec<f64> = h.get_solution().col_value.clone();
    let mut j: HighsInt = 0;
    for k in 1..n {
        if x[k as usize] > x[j as usize] {
            j = k;
        }
    }
    let mut lo = h.get_lp().col_lower_[j as usize];
    let mut up = h.get_lp().col_upper_[j as usize];
    let xj = x[j as usize];
    if xj > 0.0 {
        // std::max(lo, 0.5 * x[j])
        let v = 0.5 * xj;
        up = if lo < v { v } else { lo };
    } else {
        // std::min(up, x[j] + 1)
        let v = xj + 1.0;
        lo = if v < up { v } else { up };
    }
    h.change_cols_bounds_set(1, &[j], &[lo], &[up]);
    let _ = h.run();
    report(&h, 1);

    x = h.get_solution().col_value.clone();
    let mut index: Vec<HighsInt> = Vec::new();
    let mut value: Vec<f64> = Vec::new();
    let mut sum: f64 = 0.0;
    for (k, &xk) in x.iter().enumerate().take(n as usize) {
        if xk != 0.0 {
            index.push(k as HighsInt);
            value.push(1.0);
            sum += xk;
        }
    }
    if index.is_empty() {
        println!("stage 2 skip");
        return ExitCode::SUCCESS;
    }
    let row_lower = -kHighsInf;
    let row_upper = if sum > 0.0 { 0.9 * sum } else { sum - 1.0 };
    let start: HighsInt = 0;
    h.add_rows(
        1,
        &[row_lower],
        &[row_upper],
        index.len() as HighsInt,
        &[start],
        &index,
        &value,
    );
    let _ = h.run();
    report(&h, 2);

    let num_row = h.get_num_row();
    let mut mask: Vec<HighsInt> = vec![0; num_row as usize];
    if let Some(last) = mask.last_mut() {
        *last = 1;
    }
    h.delete_rows_mask(&mut mask);
    let _ = h.run();
    report(&h, 3);
    ExitCode::SUCCESS
}
