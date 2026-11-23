//! HiGHS source: lp_data/HighsLpUtils.cpp (isLessInfeasibleDSECandidate, appendColsToLpVectors, appendRowsToLpVectors)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use crate::lp_data::highs_lp::HighsLp;
use crate::prelude::*;

/// HiGHS: HighsLpUtils.cpp::appendColsToLpVectors (lines 1588-1614, 27 lines)
pub fn append_cols_to_lp_vectors(
    lp: &mut HighsLp,
    num_new_col: HighsInt,
    col_cost: &[f64],
    col_lower: &[f64],
    col_upper: &[f64],
) {
    debug_assert!(num_new_col >= 0);
    if num_new_col == 0 {
        return;
    }
    let new_num_col: HighsInt = lp.num_col_ + num_new_col;
    lp.col_cost_.resize(new_num_col as usize, 0.0);
    lp.col_lower_.resize(new_num_col as usize, 0.0);
    lp.col_upper_.resize(new_num_col as usize, 0.0);
    let have_integrality: bool = !lp.integrality_.is_empty();
    if have_integrality {
        debug_assert!(lp.integrality_.len() as HighsInt == lp.num_col_);
        lp.integrality_
            .resize(new_num_col as usize, HighsVarType::kContinuous);
    }
    let have_names: bool = !lp.col_names_.is_empty();
    if have_names {
        lp.col_names_.resize(new_num_col as usize, String::new());
    }
    for new_col in 0..num_new_col {
        let i_col: HighsInt = lp.num_col_ + new_col;
        lp.col_cost_[i_col as usize] = col_cost[new_col as usize];
        lp.col_lower_[i_col as usize] = col_lower[new_col as usize];
        lp.col_upper_[i_col as usize] = col_upper[new_col as usize];
        // Cannot guarantee to create unique names, so name is blank
        if have_names {
            lp.col_names_[i_col as usize] = String::new();
        }
        if have_integrality {
            lp.integrality_[i_col as usize] = HighsVarType::kContinuous;
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::appendRowsToLpVectors (lines 1616-1634, 19 lines)
pub fn append_rows_to_lp_vectors(
    lp: &mut HighsLp,
    num_new_row: HighsInt,
    row_lower: &[f64],
    row_upper: &[f64],
) {
    debug_assert!(num_new_row >= 0);
    if num_new_row == 0 {
        return;
    }
    let new_num_row: HighsInt = lp.num_row_ + num_new_row;
    lp.row_lower_.resize(new_num_row as usize, 0.0);
    lp.row_upper_.resize(new_num_row as usize, 0.0);
    let have_names: bool = !lp.row_names_.is_empty();
    if have_names {
        lp.row_names_.resize(new_num_row as usize, String::new());
    }

    for new_row in 0..num_new_row {
        let i_row: HighsInt = lp.num_row_ + new_row;
        lp.row_lower_[i_row as usize] = row_lower[new_row as usize];
        lp.row_upper_[i_row as usize] = row_upper[new_row as usize];
        // Cannot guarantee to create unique names, so name is blank
        if have_names {
            lp.row_names_[i_row as usize] = String::new();
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::isLessInfeasibleDSECandidate (lines 3162-3200, 39 lines)
pub fn is_less_infeasible_dse_candidate(log_options: &HighsLogOptions, lp: &HighsLp) -> bool {
    let mut max_col_num_en: HighsInt = -1;
    let max_allowed_col_num_en: HighsInt = 24;
    let max_assess_col_num_en: HighsInt = cpp_max(9, max_allowed_col_num_en);
    let max_average_col_num_en: HighsInt = 6;
    let mut col_length_k: Vec<HighsInt> = vec![0; (1 + max_assess_col_num_en) as usize];
    let mut lidse_candidate: bool = true;
    for col in 0..lp.num_col_ {
        // Check limit on number of entries in the column has not been breached
        let col_num_en: HighsInt =
            lp.a_matrix_.start_[(col + 1) as usize] - lp.a_matrix_.start_[col as usize];
        max_col_num_en = cpp_max(col_num_en, max_col_num_en);
        if col_num_en > max_assess_col_num_en {
            return false;
        }
        col_length_k[col_num_en as usize] += 1;
        for en in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
            let value: f64 = lp.a_matrix_.value_[en as usize];
            // All nonzeros must be +1 or -1
            if value.abs() != 1.0 {
                return false;
            }
        }
    }
    let mut average_col_num_en: f64 = lp.a_matrix_.start_[lp.num_col_ as usize] as f64;
    average_col_num_en = average_col_num_en / (lp.num_col_ as f64);
    lidse_candidate = lidse_candidate && average_col_num_en <= (max_average_col_num_en as f64);
    crate::highs_log_dev!(
        log_options,
        HighsLogType::kInfo,
        "LP {} has all |entries|=1; max column count = {} (limit {}); average column count = {} (limit {}): LP {} a candidate for LiDSE\n",
        lp.model_name_,
        max_col_num_en,
        max_allowed_col_num_en,
        average_col_num_en,
        max_average_col_num_en,
        if lidse_candidate { "is" } else { "is not" }
    );
    lidse_candidate
}
