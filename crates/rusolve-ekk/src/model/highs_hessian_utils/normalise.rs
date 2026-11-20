//! HiGHS source: model/HighsHessianUtils.cpp (lines 320-627)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//!

use crate::model::highs_hessian::HighsHessian;
use crate::prelude::*;
use crate::util::c_printf::{c_fmt_g, c_fmt_g_width};

/// HiGHS: HighsHessianUtils.cpp::normaliseHessian (lines 320-566)
pub fn normalise_hessian(options: &HighsOptions, hessian: &mut HighsHessian) -> HighsStatus {
    let dim = hessian.dim_;
    let triangular = hessian.format_ == HessianFormat::kTriangular;
    let square = !triangular;
    debug_assert!((hessian.format_ == HessianFormat::kSquare) == square);
    let mut upper = HighsSparseMatrix::default();
    let mut upper_length: Vec<HighsInt> = vec![0; dim as usize];
    let mut num_upper_nz: HighsInt = 0;
    let mut num_lower_nz: HighsInt = 0;
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        let mut i_el = hessian.start_[i_col as usize];
        while i_el < hessian.start_[(i_col + 1) as usize] {
            let i_row = hessian.index_[i_el as usize];
            if i_row < i_col {
                upper_length[i_row as usize] += 1;
                num_upper_nz += 1;
            } else if i_row > i_col {
                num_lower_nz += 1;
            }
            i_el += 1;
        }
        i_col += 1;
    }
    if triangular && num_upper_nz == 0 {
        return HighsStatus::kOk;
    }
    // If square, do quick check for asymmetry
    if square && num_upper_nz != num_lower_nz {
        highs_log_user(
            &options.log_options,
            HighsLogType::kError,
            &format!(
                "Hessian has {} / {} lower / upper triangular entries so is not symmetric\n",
                num_lower_nz, num_upper_nz
            ),
        );
        return HighsStatus::kError;
    }
    if num_upper_nz > 0 {
        // Form the HighsSparseMatrix of strict upper triangular entries
        num_lower_nz = 0;
        // Define the upper starts, and use upper_length to store the
        // element number for the next nonzero in each row
        let mut i_row: HighsInt = 0;
        while i_row < dim {
            upper
                .start_
                .push(upper.start_[i_row as usize] + upper_length[i_row as usize]);
            upper_length[i_row as usize] = upper.start_[i_row as usize];
            i_row += 1;
        }
        upper.index_.resize(num_upper_nz as usize, 0);
        upper.value_.resize(num_upper_nz as usize, 0.0);
        let mut i_col = 0;
        while i_col < dim {
            // Loop through the Hessian column, moving its upper triangular
            // nonzeros to upper, so have to shift the lower nonzeros
            let from_el = hessian.start_[i_col as usize];
            // Define the new start for this lower column
            hessian.start_[i_col as usize] = num_lower_nz;
            let mut i_el = from_el;
            while i_el < hessian.start_[(i_col + 1) as usize] {
                let i_row = hessian.index_[i_el as usize];
                if i_row < i_col {
                    // Store the nonzero in the upper row
                    upper.index_[upper_length[i_row as usize] as usize] = i_col;
                    upper.value_[upper_length[i_row as usize] as usize] =
                        hessian.value_[i_el as usize];
                    upper_length[i_row as usize] += 1;
                } else {
                    // Shift the nonzero in the lower column
                    hessian.index_[num_lower_nz as usize] = i_row;
                    hessian.value_[num_lower_nz as usize] = hessian.value_[i_el as usize];
                    num_lower_nz += 1;
                }
                i_el += 1;
            }
            i_col += 1;
        }
        hessian.start_[dim as usize] = num_lower_nz;
        // Check that upper_length has reached the start of the next row
        let mut i_row = 0;
        while i_row < dim {
            debug_assert!(upper_length[i_row as usize] == upper.start_[(i_row + 1) as usize]);
            i_row += 1;
        }
        debug_assert!(upper_length[(dim - 1) as usize] == num_upper_nz);
    } else {
        upper.start_.resize((dim + 1) as usize, 0);
    }
    // Have to work from a copy of the Hessian if it is triangular and
    // there are upper triangular entries to insert
    let hessian_copy = if triangular {
        hessian.clone()
    } else {
        HighsHessian::default()
    };

    // Determine the lower (upper) off-diagonal entries in each column
    // (row) of the Hessian
    let mut lower_on_below_diagonal: Vec<f64> = vec![0.0; dim as usize];
    let mut upper_off_diagonal: Vec<f64> = vec![0.0; dim as usize];
    // Should be no duplicates
    let mut debug_num_duplicate: HighsInt = 0;
    let mut num_summation: HighsInt = 0;
    let mut num_upper_triangle: HighsInt = 0;
    let mut num_hessian_el: HighsInt = 0;
    let k_square_hessian_asymmetry_tolerance: f64 = 1e-10;
    let mut num_illegal_asymmetry: HighsInt = 0;
    let mut min_illegal_asymmetry: f64 = kHighsInf;
    let mut max_illegal_asymmetry: f64 = 0.0;
    let mut num_ok_asymmetry: HighsInt = 0;
    let mut min_ok_asymmetry: f64 = kHighsInf;
    let mut max_ok_asymmetry: f64 = 0.0;

    // Dropped: lines 417 (const bool expensive_2821_check = true;)
    let mut i_col: HighsInt = 0;
    while i_col < dim {
        let mut i_el = (if triangular { &hessian_copy } else { &*hessian }).start_[i_col as usize];
        while i_el
            < (if triangular { &hessian_copy } else { &*hessian }).start_[(i_col + 1) as usize]
        {
            let i_row = (if triangular { &hessian_copy } else { &*hessian }).index_[i_el as usize];
            debug_assert!(i_row >= i_col);
            if lower_on_below_diagonal[i_row as usize] != 0.0 {
                debug_num_duplicate += 1;
            }
            lower_on_below_diagonal[i_row as usize] =
                (if triangular { &hessian_copy } else { &*hessian }).value_[i_el as usize];
            i_el += 1;
        }
        // Now look at the corresponding row of the upper triangular
        // matrix - deliberately referring to it as iCol and picking out
        // its indices as iRow
        let mut i_el = upper.start_[i_col as usize];
        while i_el < upper.start_[(i_col + 1) as usize] {
            let i_row = upper.index_[i_el as usize];
            debug_assert!(i_row > i_col);
            if upper_off_diagonal[i_row as usize] != 0.0 {
                debug_num_duplicate += 1;
            }
            upper_off_diagonal[i_row as usize] = upper.value_[i_el as usize];
            if square {
                // When square, ensure that the upper triangular value matches
                // the corresponding lower triangular value to within the
                // tolerance used by JuMP
                let asymmetry = (upper_off_diagonal[i_row as usize]
                    - lower_on_below_diagonal[i_row as usize])
                    .abs();
                if asymmetry > k_square_hessian_asymmetry_tolerance {
                    num_illegal_asymmetry += 1;
                    min_illegal_asymmetry = cpp_min(asymmetry, min_illegal_asymmetry);
                    max_illegal_asymmetry = cpp_max(asymmetry, max_illegal_asymmetry);
                } else if asymmetry != 0.0 {
                    num_ok_asymmetry += 1;
                    min_ok_asymmetry = cpp_min(asymmetry, min_ok_asymmetry);
                    max_ok_asymmetry = cpp_max(asymmetry, max_ok_asymmetry);
                    let average = (upper_off_diagonal[i_row as usize]
                        + lower_on_below_diagonal[i_row as usize])
                        * 0.5;
                    upper_off_diagonal[i_row as usize] = average;
                    lower_on_below_diagonal[i_row as usize] = average;
                }
                // Don't zero the upper off diagonal entry, so that it's
                // possible to check that nonzeros in the lower off diagonal
                // match the upper off diagonal entry
            } else {
                // When triangular, add the upper triangular value to the
                // lower triangular value, and zero the upper triangular value
                // so that it doesn't generate a duplicate
                debug_assert!(triangular);
                num_upper_triangle += 1;
                if lower_on_below_diagonal[i_row as usize] != 0.0 {
                    num_summation += 1;
                }
                lower_on_below_diagonal[i_row as usize] += upper_off_diagonal[i_row as usize];
                upper_off_diagonal[i_row as usize] = 0.0;
            }
            i_el += 1;
        }
        // Now gather the nonzeros, zeroing lower_on_below_diagonal and
        // upper_off_diagonal
        let from_el = (if triangular { &hessian_copy } else { &*hessian }).start_[i_col as usize];
        hessian.start_[i_col as usize] = num_hessian_el;
        // Assign any nonzero diagonal entry
        let mut i_row = i_col;
        if lower_on_below_diagonal[i_row as usize] != 0.0 {
            hessian.index_[num_hessian_el as usize] = i_row;
            hessian.value_[num_hessian_el as usize] = lower_on_below_diagonal[i_row as usize];
            num_hessian_el += 1;
            lower_on_below_diagonal[i_row as usize] = 0.0;
        }
        // Assign nonzeros below the diagonal
        let mut i_el = from_el;
        while i_el
            < (if triangular { &hessian_copy } else { &*hessian }).start_[(i_col + 1) as usize]
        {
            let i_row = (if triangular { &hessian_copy } else { &*hessian }).index_[i_el as usize];
            if lower_on_below_diagonal[i_row as usize] != 0.0 {
                hessian.index_[num_hessian_el as usize] = i_row;
                hessian.value_[num_hessian_el as usize] = lower_on_below_diagonal[i_row as usize];
                num_hessian_el += 1;
                lower_on_below_diagonal[i_row as usize] = 0.0;
            }
            i_el += 1;
        }
        // Assign nonzeros above the diagonal
        let mut i_el = upper.start_[i_col as usize];
        while i_el < upper.start_[(i_col + 1) as usize] {
            let i_row = upper.index_[i_el as usize];
            if triangular {
                debug_assert!(upper_off_diagonal[i_row as usize] == 0.0);
            }
            // Look for nonzeros in lower_on_below_diagonal created by
            // adding an entry that was in the upper triangle and added in
            if lower_on_below_diagonal[i_row as usize] != 0.0 {
                hessian.index_[num_hessian_el as usize] = i_row;
                hessian.value_[num_hessian_el as usize] = lower_on_below_diagonal[i_row as usize];
                num_hessian_el += 1;
                lower_on_below_diagonal[i_row as usize] = 0.0;
            }
            // Any upper entry retained for checking symmetry for square
            // Hessians needs to be zeroed
            upper_off_diagonal[i_row as usize] = 0.0;
            i_el += 1;
        }
        // Dropped: lines 509-518 (expensive_2821_check block)
        i_col += 1;
    } // Loop iCol = 0; iCol < dim; iCol++
    hessian.start_[dim as usize] = num_hessian_el;
    hessian.format_ = HessianFormat::kTriangular;
    hessian.index_.resize(num_hessian_el as usize, 0);
    hessian.value_.resize(num_hessian_el as usize, 0.0);

    debug_assert!(debug_num_duplicate == 0);

    let mut warning_found = false;
    let mut error_found = false;
    if num_ok_asymmetry != 0 {
        debug_assert!(square);
        highs_log_user(
            &options.log_options,
            HighsLogType::kInfo,
            &format!(
                "Square Hessian contains {} non-symmetr{} in [{}, {}] within tolerance of {}\n",
                num_ok_asymmetry,
                if num_ok_asymmetry == 1 { "y" } else { "ies" },
                c_fmt_g(min_ok_asymmetry, 2),
                c_fmt_g(max_ok_asymmetry, 2),
                c_fmt_g(k_square_hessian_asymmetry_tolerance, 1)
            ),
        );
    }
    if num_illegal_asymmetry != 0 {
        debug_assert!(square);
        highs_log_user(
            &options.log_options,
            HighsLogType::kError,
            &format!(
                "Square Hessian contains {} non-symmetr{} in [{}, {}] exceeding tolerance of {}\n",
                num_illegal_asymmetry,
                if num_illegal_asymmetry == 1 {
                    "y"
                } else {
                    "ies"
                },
                c_fmt_g(min_illegal_asymmetry, 2),
                c_fmt_g(max_illegal_asymmetry, 2),
                c_fmt_g(k_square_hessian_asymmetry_tolerance, 1)
            ),
        );
        error_found = true;
    }
    if num_upper_triangle != 0 {
        debug_assert!(triangular);
        highs_log_user(
            &options.log_options,
            HighsLogType::kWarning,
            &format!(
                "Triangular Hessian contains {} entr{} in upper triangle: added to lower triangle, requiring {} non-trivial summation{}\n",
                num_upper_triangle,
                if num_upper_triangle == 1 { "y" } else { "ies" },
                num_summation,
                if num_summation == 1 { "" } else { "s" }
            ),
        );
        debug_assert!(num_summation <= num_upper_triangle);
        warning_found = true;
    }
    let mut return_status = HighsStatus::kOk;
    if error_found {
        return_status = HighsStatus::kError;
    } else if warning_found {
        return_status = HighsStatus::kWarning;
    }
    return_status
}

/// HiGHS: HighsHessianUtils.cpp::completeHessian (lines 568-584)
pub fn complete_hessian(full_dim: HighsInt, hessian: &mut HighsHessian) {
    // Ensure that any non-zero Hessian of dimension less than the
    // number of columns in the model is completed with explicit zero
    // diagonal entries
    debug_assert!(hessian.dim_ <= full_dim);
    if hessian.dim_ == full_dim {
        return;
    }
    let mut nnz = hessian.num_nz();
    hessian.exact_resize();
    let mut i_col = hessian.dim_;
    while i_col < full_dim {
        hessian.index_.push(i_col);
        hessian.value_.push(0.0);
        nnz += 1;
        hessian.start_.push(nnz);
        i_col += 1;
    }
    hessian.dim_ = full_dim;
    debug_assert!(hessian.start_.len() as HighsInt == hessian.dim_ + 1);
}

/// HiGHS: HighsHessianUtils.cpp::reportHessian (lines 586-604)
pub fn report_hessian(
    log_options: &HighsLogOptions,
    dim: HighsInt,
    num_nz: HighsInt,
    start: &[HighsInt],
    index: &[HighsInt],
    value: &[f64],
) {
    if dim <= 0 {
        return;
    }
    highs_log_user(
        log_options,
        HighsLogType::kInfo,
        "Hessian Index              Value\n",
    );
    let mut col: HighsInt = 0;
    while col < dim {
        highs_log_user(
            log_options,
            HighsLogType::kInfo,
            &format!("    {:8} Start   {:10}\n", col, start[col as usize]),
        );
        let to_el = if col < dim - 1 {
            start[(col + 1) as usize]
        } else {
            num_nz
        };
        let mut el = start[col as usize];
        while el < to_el {
            highs_log_user(
                log_options,
                HighsLogType::kInfo,
                &format!(
                    "          {:8} {}\n",
                    index[el as usize],
                    c_fmt_g_width(value[el as usize], 6, 12)
                ),
            );
            el += 1;
        }
        col += 1;
    }
    highs_log_user(
        log_options,
        HighsLogType::kInfo,
        &format!("             Start   {:10}\n", num_nz),
    );
}

/// HiGHS: HighsHessianUtils.cpp::userScaleHessian (lines 606-627)
///
/// C++ default arguments: apply = true
pub fn user_scale_hessian(hessian: &mut HighsHessian, data: &mut HighsUserScaleData, apply: bool) {
    data.num_infinite_hessian_values = 0;
    if hessian.dim_ == 0 {
        return;
    }
    let user_objective_scale = data.user_objective_scale;
    let user_bound_scale = data.user_bound_scale;
    if user_objective_scale == 0 && user_bound_scale == 0 {
        return;
    }
    // If variable bounds are scaled by bound_scale_value, then linear
    // term in objective is scaled by bound_scale_value, but Hessian
    // term is scaled by bound_scale_value**2, so have to scale down the
    // Hessian values so linear and Hessian terms are both scaled by the
    // same constant value
    let objective_scale_value = 2.0_f64.powf(user_objective_scale as f64);
    let bound_scale_value = 2.0_f64.powf(-(user_bound_scale as f64));
    let mut i_el: HighsInt = 0;
    while i_el < hessian.start_[hessian.dim_ as usize] {
        let value = hessian.value_[i_el as usize] * objective_scale_value * bound_scale_value;
        if value.abs() > data.infinite_cost {
            data.num_infinite_hessian_values += 1;
        }
        if apply {
            hessian.value_[i_el as usize] = value;
        }
        i_el += 1;
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

    #[test]
    fn complete_hessian_pads_zero_diagonals() {
        let mut h = HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 2, 3],
            index_: vec![0, 1, 1],
            value_: vec![2.0, 1.0, 4.0],
        };
        complete_hessian(3, &mut h);
        assert_eq!(h.dim_, 3);
        assert_eq!(h.start_, vec![0, 2, 3, 4]);
        assert_eq!(h.index_, vec![0, 1, 1, 2]);
        assert_eq!(h.value_, vec![2.0, 1.0, 4.0, 0.0]);
    }

    #[test]
    fn normalise_keeps_triangular() {
        let options = HighsOptions::default();
        let mut h = HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 2, 3],
            index_: vec![0, 1, 1],
            value_: vec![2.0, 1.0, 4.0],
        };
        assert_eq!(normalise_hessian(&options, &mut h), HighsStatus::kOk);
        assert_eq!(h.format_, HessianFormat::kTriangular);
        assert_eq!(h.num_nz(), 3);
    }
}
