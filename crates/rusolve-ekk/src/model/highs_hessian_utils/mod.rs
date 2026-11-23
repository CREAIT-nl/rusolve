//! HiGHS source: model/HighsHessianUtils.cpp (lines 22-318; see `normalise.rs` for 320-627)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Free functions; reached by path (`crate::model::highs_hessian_utils::*`
//! is also in the prelude: no name clashes).

pub mod normalise;
pub use normalise::*;

use crate::model::highs_hessian::HighsHessian;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g;
use crate::util::highs_matrix_utils::{assess_matrix, assess_matrix_dimensions};

/// HiGHS: HighsHessianUtils.cpp::assessHessian (lines 22-90)
pub fn assess_hessian(hessian: &mut HighsHessian, options: &HighsOptions) -> HighsStatus {
    let mut return_status = HighsStatus::kOk;

    return_status = interpret_call_status(
        &options.log_options,
        assess_hessian_dimensions(options, hessian),
        return_status,
        "assessHessianDimensions",
    );
    if return_status == HighsStatus::kError {
        return return_status;
    }

    // If the Hessian has no columns there is nothing left to test
    if hessian.dim_ == 0 {
        hessian.clear();
        return HighsStatus::kOk;
    }

    // Assess the Hessian matrix
    //
    // The start of column 0 must be zero.
    if hessian.start_[0] != 0 {
        highs_log_user(
            &options.log_options,
            HighsLogType::kError,
            &format!(
                "Hessian has nonzero value ({}) for the start of column 0\n",
                hessian.start_[0]
            ),
        );
        return HighsStatus::kError;
    }
    // Assess Q, summing duplicates, but deferring the assessment of
    // values (other than those which are identically zero)
    let sum_duplicates = true;
    let mut matrix_p_end = Vec::new();
    let call_status = assess_matrix(
        &options.log_options,
        "Hessian",
        hessian.dim_,
        hessian.dim_,
        false,
        &mut hessian.start_,
        &mut matrix_p_end,
        &mut hessian.index_,
        &mut hessian.value_,
        0.0,
        kHighsInf,
        sum_duplicates,
        options.h16_assess_small,
    );
    return_status = interpret_call_status(
        &options.log_options,
        call_status,
        return_status,
        "assessMatrix",
    );
    if return_status == HighsStatus::kError {
        return return_status;
    }
    // Transform the Hessian to pure column-wise lower triangle format
    let call_status = normalise_hessian(options, hessian);
    return_status = interpret_call_status(
        &options.log_options,
        call_status,
        return_status,
        "normaliseHessian",
    );
    if return_status == HighsStatus::kError {
        return return_status;
    }
    // Assess values in Q
    let mut matrix_p_end = Vec::new();
    let call_status = assess_matrix(
        &options.log_options,
        "Hessian",
        hessian.dim_,
        hessian.dim_,
        false,
        &mut hessian.start_,
        &mut matrix_p_end,
        &mut hessian.index_,
        &mut hessian.value_,
        options.small_matrix_value,
        options.large_matrix_value,
        false,
        options.h16_assess_small,
    );
    return_status = interpret_call_status(
        &options.log_options,
        call_status,
        return_status,
        "assessMatrix",
    );
    if return_status == HighsStatus::kError {
        return return_status;
    }

    let mut hessian_num_nz = hessian.num_nz();
    // If the Hessian has nonzeros, complete its diagonal with explicit
    // zeros if necessary
    if hessian_num_nz != 0 {
        complete_hessian_diagonal(options, hessian);
        hessian_num_nz = hessian.num_nz();
    }
    // If entries have been removed from the matrix, resize the index
    // and value vectors
    if (hessian.index_.len() as HighsInt) > hessian_num_nz {
        hessian.index_.resize(hessian_num_nz as usize, 0);
    }
    if (hessian.value_.len() as HighsInt) > hessian_num_nz {
        hessian.value_.resize(hessian_num_nz as usize, 0.0);
    }

    if return_status != HighsStatus::kError {
        return_status = HighsStatus::kOk;
    }
    if return_status != HighsStatus::kOk {
        highs_log_dev(
            &options.log_options,
            HighsLogType::kInfo,
            &format!(
                "assessHessian returns HighsStatus = {}\n",
                highs_status_to_string(return_status)
            ),
        );
    }
    return_status
}

/// HiGHS: HighsHessianUtils.cpp::assessHessianDimensions (lines 92-102)
pub fn assess_hessian_dimensions(
    options: &HighsOptions,
    hessian: &mut HighsHessian,
) -> HighsStatus {
    if hessian.dim_ == 0 {
        return HighsStatus::kOk;
    }

    // Assess the Hessian dimensions and vector sizes
    let hessian_p_end: Vec<HighsInt> = Vec::new();
    let partitioned = false;
    assess_matrix_dimensions(
        &options.log_options,
        hessian.dim_,
        partitioned,
        &hessian.start_,
        &hessian_p_end,
        &hessian.index_,
        &hessian.value_,
    )
}

/// HiGHS: HighsHessianUtils.cpp::completeHessianDiagonal (lines 104-169)
pub fn complete_hessian_diagonal(options: &HighsOptions, hessian: &mut HighsHessian) {
    // Count the number of missing diagonal entries
    let mut num_missing_diagonal_entries: HighsInt = 0;
    let dim = hessian.dim_;
    let num_nz = hessian.num_nz();
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        let i_el = hessian.start_[i_col as usize];
        if i_el < num_nz {
            if hessian.index_[i_el as usize] != i_col {
                num_missing_diagonal_entries += 1;
            }
        } else {
            num_missing_diagonal_entries += 1;
        }
        i_col += 1;
    }
    if num_missing_diagonal_entries > 0 {
        highs_log_dev(
            &options.log_options,
            HighsLogType::kInfo,
            &format!(
                "Hessian has dimension {} and {} nonzeros: inserting {} zeros onto the diagonal\n",
                dim, num_nz, num_missing_diagonal_entries
            ),
        );
    }
    debug_assert!(num_missing_diagonal_entries >= dim - num_nz);
    if num_missing_diagonal_entries == 0 {
        return;
    }
    // There are missing diagonal entries to be inserted as explicit zeros
    let new_num_nz = hessian.num_nz() + num_missing_diagonal_entries;
    let mut to_i_el = new_num_nz;
    hessian.index_.resize(new_num_nz as usize, 0);
    hessian.value_.resize(new_num_nz as usize, 0.0);
    let mut next_start = hessian.num_nz();
    hessian.start_[dim as usize] = to_i_el;
    let mut num_missing_diagonal_entries_added: HighsInt = 0;
    let mut i_col = dim - 1;
    while i_col >= 0 {
        // Shift the entries that are sure to be off-diagonal
        let mut i_el = next_start - 1;
        while i_el > hessian.start_[i_col as usize] {
            debug_assert!(hessian.index_[i_el as usize] != i_col);
            to_i_el -= 1;
            hessian.index_[to_i_el as usize] = hessian.index_[i_el as usize];
            hessian.value_[to_i_el as usize] = hessian.value_[i_el as usize];
            i_el -= 1;
        }
        // Now consider any first entry. If there is none, or if it's not
        // the diagonal, then there is no diagonal entry for this column
        let no_diagonal_entry = if hessian.start_[i_col as usize] < next_start {
            let i_el = hessian.start_[i_col as usize];
            // Copy the first entry
            to_i_el -= 1;
            hessian.index_[to_i_el as usize] = hessian.index_[i_el as usize];
            hessian.value_[to_i_el as usize] = hessian.value_[i_el as usize];
            // If the first entry isn't the diagonal, then there is no
            // diagonal entry for this column
            hessian.index_[i_el as usize] != i_col
        } else {
            true
        };
        if no_diagonal_entry {
            // There is no diagonal entry, so have insert an explicit zero
            to_i_el -= 1;
            hessian.index_[to_i_el as usize] = i_col;
            hessian.value_[to_i_el as usize] = 0.0;
            num_missing_diagonal_entries_added += 1;
            debug_assert!(num_missing_diagonal_entries_added <= num_missing_diagonal_entries);
        }
        next_start = hessian.start_[i_col as usize];
        hessian.start_[i_col as usize] = to_i_el;
        i_col -= 1;
    }
    debug_assert!(to_i_el == 0);
}

/// HiGHS: HighsHessianUtils.cpp::okHessianDiagonal (lines 171-209)
///
/// C++ default arguments: sense = ObjSense::kMinimize
pub fn ok_hessian_diagonal(
    options: &HighsOptions,
    hessian: &mut HighsHessian,
    sense: ObjSense,
) -> bool {
    let mut min_diagonal_value = kHighsInf;
    let mut max_diagonal_value = -kHighsInf;
    let dim = hessian.dim_;
    let sense_sign = sense as HighsInt;
    let mut num_illegal_diagonal_value: HighsInt = 0;
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        // Assumes that the diagonal entry is always first, possibly with explicit
        // zero value
        let i_el = hessian.start_[i_col as usize];
        debug_assert!(hessian.index_[i_el as usize] == i_col);
        let diagonal_value = sense_sign as f64 * hessian.value_[i_el as usize];
        min_diagonal_value = cpp_min(diagonal_value, min_diagonal_value);
        max_diagonal_value = cpp_max(diagonal_value, max_diagonal_value);
        // Diagonal entries signed by sense must be non-negative
        if diagonal_value < 0.0 {
            num_illegal_diagonal_value += 1;
        }
        i_col += 1;
    }

    let certainly_not_positive_semidefinite = num_illegal_diagonal_value > 0;
    if certainly_not_positive_semidefinite {
        if sense == ObjSense::kMinimize {
            highs_log_user(
                &options.log_options,
                HighsLogType::kError,
                &format!(
                    "Hessian has {} diagonal entries in [{}, 0) so is not positive semidefinite for minimization\n",
                    num_illegal_diagonal_value,
                    c_fmt_g(min_diagonal_value, 6)
                ),
            );
        } else {
            highs_log_user(
                &options.log_options,
                HighsLogType::kError,
                &format!(
                    "Hessian has {} diagonal entries in (0, {}] so is not negative semidefinite for maximization\n",
                    num_illegal_diagonal_value,
                    c_fmt_g(-min_diagonal_value, 6)
                ),
            );
        }
    }
    !certainly_not_positive_semidefinite
}

/// HiGHS: HighsHessianUtils.cpp::extractTriangularHessian (lines 211-252)
pub fn extract_triangular_hessian(
    options: &HighsOptions,
    hessian: &mut HighsHessian,
) -> HighsStatus {
    // Viewing the Hessian column-wise, remove any entries in the strict
    // upper triangle
    let mut return_status = HighsStatus::kOk;
    let dim = hessian.dim_;
    let mut nnz: HighsInt = 0;
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        let nnz0 = nnz;
        let mut i_el = hessian.start_[i_col as usize];
        while i_el < hessian.start_[(i_col + 1) as usize] {
            let i_row = hessian.index_[i_el as usize];
            if i_row < i_col {
                i_el += 1;
                continue;
            }
            hessian.index_[nnz as usize] = i_row;
            hessian.value_[nnz as usize] = hessian.value_[i_el as usize];
            if i_row == i_col && nnz > nnz0 {
                // Diagonal entry is not first in column so swap it in
                hessian.index_[nnz as usize] = hessian.index_[nnz0 as usize];
                hessian.value_[nnz as usize] = hessian.value_[nnz0 as usize];
                hessian.index_[nnz0 as usize] = i_row;
                hessian.value_[nnz0 as usize] = hessian.value_[i_el as usize];
            }
            nnz += 1;
            i_el += 1;
        }
        hessian.start_[i_col as usize] = nnz0;
        i_col += 1;
    }
    let num_ignored_nz = hessian.start_[dim as usize] - nnz;
    debug_assert!(num_ignored_nz >= 0);
    if num_ignored_nz != 0 {
        if hessian.format_ == HessianFormat::kTriangular {
            highs_log_user(
                &options.log_options,
                HighsLogType::kWarning,
                &format!(
                    "Ignored {} entries of Hessian in opposite triangle\n",
                    num_ignored_nz
                ),
            );
            return_status = HighsStatus::kWarning;
        }
        hessian.start_[dim as usize] = nnz;
    }
    debug_assert!(hessian.start_[dim as usize] == nnz);
    hessian.format_ = HessianFormat::kTriangular;
    return_status
}

/// HiGHS: HighsHessianUtils.cpp::triangularToSquareHessian (lines 254-318)
pub fn triangular_to_square_hessian(
    hessian: &HighsHessian,
    start: &mut Vec<HighsInt>,
    index: &mut Vec<HighsInt>,
    value: &mut Vec<f64>,
) {
    let dim = hessian.dim_;
    if dim <= 0 {
        start.clear();
        start.resize(1, 0);
        return;
    }
    debug_assert!(hessian.format_ == HessianFormat::kTriangular);
    // Cannot use the following, since it puts the diagonal entry first
    // in each column and, strangely, the active set QP solver seems to
    // need the entries column-by-column in row order
    // (C++ lines 266-271: commented-out toSquare() block dropped)
    let nnz = hessian.start_[dim as usize];
    let square_nnz = nnz + (nnz - dim);
    start.resize((dim + 1) as usize, 0);
    index.resize(square_nnz as usize, 0);
    value.resize(square_nnz as usize, 0.0);
    let mut length: Vec<HighsInt> = vec![0; dim as usize];
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        let i_row = hessian.index_[hessian.start_[i_col as usize] as usize];
        debug_assert!(i_row == i_col);
        length[i_col as usize] += 1;
        let mut i_el = hessian.start_[i_col as usize] + 1;
        while i_el < hessian.start_[(i_col + 1) as usize] {
            let i_row = hessian.index_[i_el as usize];
            debug_assert!(i_row > i_col);
            length[i_row as usize] += 1;
            length[i_col as usize] += 1;
            i_el += 1;
        }
        i_col += 1;
    }
    start[0] = 0;
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        start[(i_col + 1) as usize] = start[i_col as usize] + length[i_col as usize];
        i_col += 1;
    }
    debug_assert!(square_nnz == start[dim as usize]);
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        let i_el = hessian.start_[i_col as usize];
        let i_row = hessian.index_[i_el as usize];
        let to_el = start[i_col as usize];
        index[to_el as usize] = i_row;
        value[to_el as usize] = hessian.value_[i_el as usize];
        start[i_col as usize] += 1;
        let mut i_el = hessian.start_[i_col as usize] + 1;
        while i_el < hessian.start_[(i_col + 1) as usize] {
            let i_row = hessian.index_[i_el as usize];
            let to_el = start[i_row as usize];
            index[to_el as usize] = i_col;
            value[to_el as usize] = hessian.value_[i_el as usize];
            start[i_row as usize] += 1;
            let to_el = start[i_col as usize];
            index[to_el as usize] = i_row;
            value[to_el as usize] = hessian.value_[i_el as usize];
            start[i_col as usize] += 1;
            i_el += 1;
        }
        i_col += 1;
    }
    start[0] = 0;
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        start[(i_col + 1) as usize] = start[i_col as usize] + length[i_col as usize];
        i_col += 1;
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    fn hess(start: Vec<HighsInt>, index: Vec<HighsInt>, value: Vec<f64>) -> HighsHessian {
        HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: start,
            index_: index,
            value_: value,
        }
    }

    #[test]
    fn complete_diagonal_inserts_explicit_zero() {
        let options = HighsOptions::default();
        let mut h = hess(vec![0, 2, 2], vec![0, 1], vec![2.0, 1.0]);
        complete_hessian_diagonal(&options, &mut h);
        assert_eq!(h.start_, vec![0, 2, 3]);
        assert_eq!(h.index_, vec![0, 1, 1]);
        assert_eq!(h.value_, vec![2.0, 1.0, 0.0]);
        assert!(ok_hessian_diagonal(&options, &mut h, ObjSense::kMinimize));
    }

    #[test]
    fn diagonal_sign_by_sense() {
        let options = HighsOptions::default();
        let mut neg = hess(vec![0, 1, 2], vec![0, 1], vec![2.0, -1.0]);
        assert!(!ok_hessian_diagonal(
            &options,
            &mut neg,
            ObjSense::kMinimize
        ));
        let mut negdef = hess(vec![0, 1, 2], vec![0, 1], vec![-2.0, -1.0]);
        assert!(ok_hessian_diagonal(
            &options,
            &mut negdef,
            ObjSense::kMaximize
        ));
    }

    #[test]
    fn triangular_to_square_row_order() {
        let h = hess(vec![0, 2, 3], vec![0, 1, 1], vec![2.0, 1.0, 4.0]);
        let (mut s, mut i, mut v) = (Vec::new(), Vec::new(), Vec::new());
        triangular_to_square_hessian(&h, &mut s, &mut i, &mut v);
        assert_eq!(s, vec![0, 2, 4]);
        assert_eq!(i, vec![0, 1, 0, 1]);
        assert_eq!(v, vec![2.0, 1.0, 1.0, 4.0]);
    }

    #[test]
    fn extract_triangular_square_drops_upper() {
        // Square input, col 0 = rows [1, 0], col 1 = rows [0, 1]. The diagonal of
        // col 0 is swapped in with nnz == iEl, so C++ lines 229 then 231 read
        // value_[1] after overwriting it: the diagonal takes the value 3, not 2.
        let options = HighsOptions::default();
        let mut h = hess(vec![0, 2, 4], vec![1, 0, 0, 1], vec![3.0, 2.0, 3.0, 4.0]);
        h.format_ = HessianFormat::kSquare;
        assert_eq!(
            extract_triangular_hessian(&options, &mut h),
            HighsStatus::kOk
        );
        assert_eq!(h.format_, HessianFormat::kTriangular);
        assert_eq!(h.start_, vec![0, 2, 3]);
        assert_eq!(h.index_[..3], [0, 1, 1]);
        assert_eq!(h.value_[..3], [3.0, 3.0, 4.0]);
    }

    #[test]
    fn extract_triangular_warns_on_opposite_triangle() {
        let options = HighsOptions::default();
        let mut h = hess(vec![0, 1, 3], vec![0, 0, 1], vec![1.0, 5.0, 2.0]);
        assert_eq!(
            extract_triangular_hessian(&options, &mut h),
            HighsStatus::kWarning
        );
        assert_eq!(h.start_, vec![0, 1, 2]);
        assert_eq!(h.index_[..2], [0, 1]);
        assert_eq!(h.value_[..2], [1.0, 2.0]);
    }

    #[test]
    fn assess_early_returns() {
        let options = HighsOptions::default();
        let mut empty = hess(vec![0], vec![], vec![]);
        empty.dim_ = 0;
        assert_eq!(assess_hessian(&mut empty, &options), HighsStatus::kOk);
        assert_eq!(empty.start_, vec![0]);
        let mut bad_start = hess(vec![1, 2, 3], vec![0, 1, 1], vec![2.0, 1.0, 4.0]);
        assert_eq!(
            assess_hessian(&mut bad_start, &options),
            HighsStatus::kError
        );
    }

    #[test]
    fn triangular_to_square_empty() {
        let mut h = hess(vec![0], vec![], vec![]);
        h.dim_ = 0;
        let (mut s, mut i, mut v) = (vec![7, 8], Vec::new(), Vec::new());
        triangular_to_square_hessian(&h, &mut s, &mut i, &mut v);
        assert_eq!(s, vec![0]);
    }

    #[test]
    fn assess_accepts_valid_triangular() {
        let options = HighsOptions::default();
        let mut h = hess(vec![0, 2, 3], vec![0, 1, 1], vec![2.0, 1.0, 4.0]);
        assert_eq!(assess_hessian(&mut h, &options), HighsStatus::kOk);
        assert_eq!(h.num_nz(), 3);
    }
}
