//! HiGHS source: util/HighsMatrixUtils.cpp (assessMatrix, assessMatrixDimensions)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Model validation of a packed (column- or row-wise) matrix, as called by
//! `assessLp` (via `HighsSparseMatrix::assess`) and `lpDimensionsOk`.
//! The duplicate-index test uses a per-component marker array in place of
//! `HighsHashTable<HighsInt>` / `std::unordered_map`: both are only queried
//! for membership (and, when summing, the position of the first occurrence),
//! so the result is identical.

#![allow(
    clippy::int_plus_one,
    clippy::needless_late_init,
    clippy::collapsible_match,
    clippy::assertions_on_constants
)]

use crate::prelude::*;

/// HiGHS: HighsMatrixUtils.cpp::assessMatrix (lines 50-311, 262 lines)
/// The partitioned overload; the two unpartitioned overloads (lines 20-48)
/// forward here with `partitioned = false`.
/// C++ default arguments: sum_duplicates = false
/// HiGHS 1.16 (#3169, #3382): HighsMatrixUtils.cpp::assessMatrix (lines 189-226, 260-265, 294-310 of 1.16.0; merged state of reference/h16-prs/3169.diff lines 23-106 and 3382.diff lines 48-127), behind `h16_assess_small`
/// `h16_assess_small` is the option flag; false keeps the 1.15.1 code.
#[allow(clippy::too_many_arguments)]
pub fn assess_matrix(
    log_options: &HighsLogOptions,
    matrix_name: &str,
    vec_dim: HighsInt,
    num_vec: HighsInt,
    partitioned: bool,
    matrix_start: &mut Vec<HighsInt>,
    matrix_p_end: &mut Vec<HighsInt>,
    matrix_index: &mut Vec<HighsInt>,
    matrix_value: &mut Vec<f64>,
    small_matrix_value: f64,
    large_matrix_value: f64,
    sum_duplicates: bool,
    h16_assess_small: bool,
) -> HighsStatus {
    if assess_matrix_dimensions(
        log_options,
        num_vec,
        partitioned,
        matrix_start,
        matrix_p_end,
        matrix_index,
        matrix_value,
    ) == HighsStatus::kError
    {
        return HighsStatus::kError;
    }

    let mut error_found = false;
    let mut warning_found = false;

    let num_nz: HighsInt = matrix_start[num_vec as usize];
    // Assess the starts
    //
    // Check whether the first start is zero
    if matrix_start[0] != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "{} matrix start vector begins with {} rather than 0\n",
            matrix_name,
            matrix_start[0]
        );
        return HighsStatus::kError;
    }
    // Set up previous_start for a fictitious previous empty packed vector
    let mut previous_start: HighsInt = matrix_start[0];
    // Set up this_start to be the first start in case num_vec = 0
    let mut this_start: HighsInt = matrix_start[0];
    let mut this_p_end: HighsInt = 0;
    if partitioned {
        this_p_end = matrix_p_end[0];
    }
    for ix in 0..num_vec {
        this_start = matrix_start[ix as usize];
        let this_start_too_small = this_start < previous_start;
        if this_start_too_small {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "{} matrix packed vector {} has illegal start of {} < {} = previous start\n",
                matrix_name,
                ix,
                this_start,
                previous_start
            );
            return HighsStatus::kError;
        }
        if partitioned {
            this_p_end = matrix_p_end[ix as usize];
            let this_p_end_too_small = this_p_end < this_start;
            if this_p_end_too_small {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kError,
                    "{} matrix packed vector {} has illegal partition end of {} < {} =  start\n",
                    matrix_name,
                    ix,
                    this_p_end,
                    this_start
                );
                return HighsStatus::kError;
            }
        }
        previous_start = this_start;
    }
    let this_start_too_big = this_start > num_nz;
    if this_start_too_big {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "{} matrix packed vector {} has illegal start of {} > {} = number of nonzeros\n",
            matrix_name,
            num_vec,
            this_start,
            num_nz
        );
        return HighsStatus::kError;
    }
    if partitioned {
        let this_p_end_too_big = this_p_end > num_nz;
        if this_p_end_too_big {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "{} matrix packed vector {} has illegal partition end of {} > {} = number of nonzeros\n",
                matrix_name,
                num_vec,
                this_p_end,
                num_nz
            );
            return HighsStatus::kError;
        }
    }
    // Assess the indices and values
    // Count the number of acceptable indices/values
    let mut num_new_nz: HighsInt = 0;
    let mut num_small_value: HighsInt = 0;
    let mut max_small_value: f64 = 0.0;
    let mut min_small_value: f64 = kHighsInf;
    let mut num_large_value: HighsInt = 0;
    let mut max_large_value: f64 = 0.0;
    let mut min_large_value: f64 = kHighsInf;
    let mut num_duplicate: HighsInt = 0;
    let illegal_el: HighsInt = -1;
    // Membership marker standing in for highs_hash / index_el_map: for each
    // component, the (new) entry number of its occurrence in the current
    // packed vector, or illegal_el. Entries touched are reset per vector,
    // exactly like highs_hash.clear() / index_el_map.clear().
    // HighsMatrixUtils.cpp:153 uses std::unordered_map index_el_map, so its
    // cost is O(entries seen). A dense vec_dim marker costs O(vec_dim) per call
    // (quadratic when rows are added one at a time to a wide model), so it is
    // used only when the call has enough entries to amortise it; otherwise a
    // HashMap, as HiGHS does. Same results either way.
    let use_dense = (num_nz as i64) * 4 >= (vec_dim as i64);
    let mut index_el_map: Vec<HighsInt> = if use_dense {
        vec![illegal_el; cpp_max(vec_dim, 0) as usize]
    } else {
        Vec::new()
    };
    let mut index_el_hash: std::collections::HashMap<HighsInt, HighsInt> =
        std::collections::HashMap::new();
    let mut touched: Vec<HighsInt> = Vec::new();

    for ix in 0..num_vec {
        let mut from_el: HighsInt = matrix_start[ix as usize];
        let mut to_el: HighsInt = matrix_start[(ix + 1) as usize];
        // Account for any index-value pairs removed so far
        matrix_start[ix as usize] = num_new_nz;
        for el in from_el..to_el {
            // Check the index
            let component: HighsInt = matrix_index[el as usize];
            // Check that the index is non-negative
            let mut legal_component = component >= 0;
            if !legal_component {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kError,
                    "{} matrix packed vector {}, entry {}, is illegal index {}\n",
                    matrix_name,
                    ix,
                    el,
                    component
                );
                return HighsStatus::kError;
            }
            // Check that the index does not exceed the vector dimension
            legal_component = component < vec_dim;
            if !legal_component {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kError,
                    "{} matrix packed vector {}, entry {}, is illegal index {:12} >= {} = vector dimension\n",
                    matrix_name,
                    ix,
                    el,
                    component,
                    vec_dim
                );
                return HighsStatus::kError;
            }
            // HiGHS 1.16 #3169/#3382 (h16_assess_small): HighsMatrixUtils.cpp lines 260-265 (1.16.0)
            // If not summing duplicates, can check matrix value magnitude
            // before checking for duplicates, so that an entry that is
            // ignored because its value is small is not identified as a
            // duplicate, and is not recorded
            if h16_assess_small && !sum_duplicates {
                let ok_value = h16_ok_matrix_value_magnitude(
                    matrix_value[el as usize],
                    large_matrix_value,
                    small_matrix_value,
                    &mut num_large_value,
                    &mut max_large_value,
                    &mut min_large_value,
                    &mut num_small_value,
                    &mut max_small_value,
                    &mut min_small_value,
                );
                if !ok_value {
                    continue;
                }
            }
            // Check whether the index has already occurred.
            let previous_el: HighsInt = if use_dense {
                index_el_map[component as usize]
            } else {
                index_el_hash.get(&component).copied().unwrap_or(illegal_el)
            };
            let is_duplicate = previous_el != illegal_el;
            if is_duplicate {
                if sum_duplicates {
                    num_duplicate += 1;
                    // Sum the duplicate entry, making sure that it's been
                    // assigned
                    debug_assert!(previous_el != illegal_el);
                    matrix_value[previous_el as usize] += matrix_value[el as usize];
                    continue;
                }
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kError,
                    "{} matrix packed vector {}, entry {}, is duplicate index {}\n",
                    matrix_name,
                    ix,
                    el,
                    component
                );
                return HighsStatus::kError;
            }
            // Not a duplicate
            // Shift the index and value of the OK entry to the new
            // position in the index and value vectors, and increment
            // the new number of nonzeros
            matrix_index[num_new_nz as usize] = matrix_index[el as usize];
            matrix_value[num_new_nz as usize] = matrix_value[el as usize];
            // Record where the index has occurred
            if use_dense {
                index_el_map[component as usize] = num_new_nz;
            } else {
                index_el_hash.insert(component, num_new_nz);
            }
            touched.push(component);
            num_new_nz += 1;
        }
        if h16_assess_small {
            // HiGHS 1.16 #3169 (h16_assess_small): HighsMatrixUtils.cpp lines 294-310 (1.16.0)
            // Now that any duplicates have been summed, can check matrix
            // value magnitude
            if sum_duplicates {
                from_el = matrix_start[ix as usize];
                to_el = num_new_nz;
                // Reset num_new_nz
                num_new_nz = matrix_start[ix as usize];
                for el in from_el..to_el {
                    let ok_value = h16_ok_matrix_value_magnitude(
                        matrix_value[el as usize],
                        large_matrix_value,
                        small_matrix_value,
                        &mut num_large_value,
                        &mut max_large_value,
                        &mut min_large_value,
                        &mut num_small_value,
                        &mut max_small_value,
                        &mut min_small_value,
                    );
                    if !ok_value {
                        continue;
                    }
                    // Shift the index and value of the OK entry to the new
                    // position in the index and value vectors, and increment
                    // the new number of nonzeros
                    matrix_index[num_new_nz as usize] = matrix_index[el as usize];
                    matrix_value[num_new_nz as usize] = matrix_value[el as usize];
                    num_new_nz += 1;
                }
            }
        } else {
            from_el = matrix_start[ix as usize];
            to_el = num_new_nz;
            // Reset num_new_nz
            num_new_nz = matrix_start[ix as usize];
            for el in from_el..to_el {
                // Check the value
                let abs_value = matrix_value[el as usize].abs();
                // Check that the value is not too large
                let large_value = abs_value >= large_matrix_value;
                if large_value {
                    if max_large_value < abs_value {
                        max_large_value = abs_value;
                    }
                    if min_large_value > abs_value {
                        min_large_value = abs_value;
                    }
                    num_large_value += 1;
                }
                let ok_value = abs_value > small_matrix_value;
                if !ok_value {
                    if max_small_value < abs_value {
                        max_small_value = abs_value;
                    }
                    if min_small_value > abs_value {
                        min_small_value = abs_value;
                    }
                    num_small_value += 1;
                }
                if ok_value {
                    // Shift the index and value of the OK entry to the new
                    // position in the index and value vectors, and increment
                    // the new number of nonzeros
                    matrix_index[num_new_nz as usize] = matrix_index[el as usize];
                    matrix_value[num_new_nz as usize] = matrix_value[el as usize];
                    num_new_nz += 1;
                }
            } // Loop from_el; to_el
        }
        if use_dense {
            for &component in &touched {
                index_el_map[component as usize] = illegal_el;
            }
        } else {
            index_el_hash.clear();
        }
        touched.clear();
    } // Loop 0; num_vec
    if num_duplicate != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "{} matrix packed vector contains {} duplicate entr{}: summed\n",
            matrix_name,
            num_duplicate,
            if num_duplicate == 1 { "y" } else { "ies" }
        );
    }
    if num_large_value != 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "{} matrix packed vector contains {} |value| in [{}, {}] greater than {}\n",
            matrix_name,
            num_large_value,
            min_large_value,
            max_large_value,
            large_matrix_value
        );
        error_found = true;
    }
    if num_small_value != 0 {
        if partitioned {
            // Shouldn't happen with a partitioned row-wise matrix since its
            // values should be OK and the code above doesn't handle p_end
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "{} matrix packed partitioned vector contains {} |value| in [{}, {}] less than or equal to {}: ignored\n",
                matrix_name,
                num_small_value,
                min_small_value,
                max_small_value,
                small_matrix_value
            );
            error_found = true;
            debug_assert!(num_small_value == 0);
        }
        // If explicit zeros are ignored, then no model information is
        // lost, so only report and return a warning if small nonzeros are
        // ignored
        if max_small_value > 0.0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kWarning,
                "{} matrix packed vector contains {} |value| in [{}, {}] less than or equal to {}: ignored\n",
                matrix_name,
                num_small_value,
                min_small_value,
                max_small_value,
                small_matrix_value
            );
            warning_found = true;
        }
    }
    matrix_start[num_vec as usize] = num_new_nz;
    let mut return_status = HighsStatus::kOk;
    if error_found {
        return_status = HighsStatus::kError;
    } else if warning_found {
        return_status = HighsStatus::kWarning;
    }
    return_status
}

/// HiGHS 1.16 (#3169): HighsMatrixUtils.cpp::assessMatrix lambda okMatrixValueMagnitude, as merged (reference/h16-prs/3169.diff lines 23-56); the 1.16.0 file's later version of the lambda is not ported
/// Returns whether the value is retained; updates the large/small counts.
#[allow(clippy::too_many_arguments)]
fn h16_ok_matrix_value_magnitude(
    value: f64,
    large_matrix_value: f64,
    small_matrix_value: f64,
    num_large_value: &mut HighsInt,
    max_large_value: &mut f64,
    min_large_value: &mut f64,
    num_small_value: &mut HighsInt,
    max_small_value: &mut f64,
    min_small_value: &mut f64,
) -> bool {
    // Check the value
    let abs_value = value.abs();
    // Check that the value is not too large
    let large_value = abs_value >= large_matrix_value;
    if large_value {
        if *max_large_value < abs_value {
            *max_large_value = abs_value;
        }
        if *min_large_value > abs_value {
            *min_large_value = abs_value;
        }
        *num_large_value += 1;
    }
    let ok_value = abs_value > small_matrix_value;
    if !ok_value {
        if *max_small_value < abs_value {
            *max_small_value = abs_value;
        }
        if *min_small_value > abs_value {
            *min_small_value = abs_value;
        }
        *num_small_value += 1;
    }
    ok_value
}

/// HiGHS: HighsMatrixUtils.cpp::assessMatrixDimensions (lines 313-371, 59 lines)
pub fn assess_matrix_dimensions(
    log_options: &HighsLogOptions,
    num_vec: HighsInt,
    partitioned: bool,
    matrix_start: &[HighsInt],
    matrix_p_end: &[HighsInt],
    matrix_index: &[HighsInt],
    matrix_value: &[f64],
) -> HighsStatus {
    let mut ok = true;
    // Assess main dimensions
    let legal_num_vec = num_vec >= 0;
    if !legal_num_vec {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "Matrix dimension validation fails on number of vectors = {} < 0\n",
            num_vec
        );
    }
    ok = legal_num_vec && ok;
    let legal_matrix_start_size = (matrix_start.len() as HighsInt) >= num_vec + 1;
    if !legal_matrix_start_size {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "Matrix dimension validation fails on start size = {} < {} = num vectors + 1\n",
            matrix_start.len(),
            num_vec + 1
        );
    }
    ok = legal_matrix_start_size && ok;
    if partitioned {
        let legal_matrix_p_end_size = (matrix_p_end.len() as HighsInt) >= num_vec + 1;
        if !legal_matrix_p_end_size {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Matrix dimension validation fails on p_end size = {} < {} = num vectors + 1\n",
                matrix_p_end.len(),
                num_vec + 1
            );
        }
        ok = (matrix_p_end.len() as HighsInt) >= num_vec + 1 && ok;
    }
    // Possibly check the sizes of the index and value vectors. Can only
    // do this with the number of nonzeros, and this is only known if
    // the start vector has a legal size. Setting num_nz = 0 otherwise
    // means that all tests pass, as they just check that the sizes of
    // the index and value vectors are non-negative.
    let num_nz: HighsInt = if legal_matrix_start_size {
        matrix_start[num_vec as usize]
    } else {
        0
    };
    if num_nz >= 0 {
        let legal_matrix_index_size = (matrix_index.len() as HighsInt) >= num_nz;
        if !legal_matrix_index_size {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Matrix dimension validation fails on index size = {} < {} = number of nonzeros\n",
                matrix_index.len(),
                num_nz
            );
        }
        ok = legal_matrix_index_size && ok;
        let legal_matrix_value_size = (matrix_value.len() as HighsInt) >= num_nz;
        if !legal_matrix_value_size {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Matrix dimension validation fails on value size = {} < {} = number of nonzeros\n",
                matrix_value.len(),
                num_nz
            );
        }
        ok = legal_matrix_value_size && ok;
    } else {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "Matrix dimension validation fails on number of nonzeros = {} < 0\n",
            num_nz
        );
        ok = false;
    }
    if ok {
        HighsStatus::kOk
    } else {
        HighsStatus::kError
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::lp_assess::assess_lp;

    /// One column with two entries in row 2, values as given (from check/TestLpValidation.cpp
    /// "LP-duplicate-small-value", 1.16.0 #3382).
    fn duplicate_small_lp(values: [f64; 2]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 3;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![1.0];
        lp.row_lower_ = vec![-kHighsInf; 3];
        lp.row_upper_ = vec![1.0; 3];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 2];
        lp.a_matrix_.index_ = vec![2, 2];
        lp.a_matrix_.value_ = values.to_vec();
        lp
    }

    fn options_with_flag(h16_assess_small: bool) -> HighsOptions {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.h16_assess_small = h16_assess_small;
        options
    }

    #[test]
    fn h16_assess_small_duplicate_small_value_on() {
        let options = options_with_flag(true);
        // A small entry is dropped before duplicate detection, so it is not
        // a duplicate of the retained entry, whatever the order
        let mut lp = duplicate_small_lp([3.0, 1e-12]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kWarning);
        assert_eq!(lp.a_matrix_.num_nz(), 1);
        assert_eq!(lp.a_matrix_.index_, vec![2]);
        assert_eq!(lp.a_matrix_.value_, vec![3.0]);

        let mut lp = duplicate_small_lp([1e-12, 3.0]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kWarning);
        assert_eq!(lp.a_matrix_.num_nz(), 1);
        assert_eq!(lp.a_matrix_.value_, vec![3.0]);

        let mut lp = duplicate_small_lp([1e-12, 1e-12]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kWarning);
        assert_eq!(lp.a_matrix_.num_nz(), 0);

        // Two retained entries in the same row are still an error
        let mut lp = duplicate_small_lp([3.0, 3.0]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kError);

        // check/TestLpValidation.cpp "LP-dimension-validation" (#3169):
        // duplicate index with zero values, both discarded, kOk since
        // max_small_value is zero
        let mut lp = duplicate_small_lp([0.0, 0.0]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kOk);
        assert_eq!(lp.a_matrix_.num_nz(), 0);
    }

    #[test]
    fn h16_assess_small_duplicate_small_value_off() {
        // 1.15.1: the small entry is still recorded, so it is a duplicate
        let options = options_with_flag(false);
        let mut lp = duplicate_small_lp([3.0, 1e-12]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kError);
        let mut lp = duplicate_small_lp([1e-12, 3.0]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kError);
        let mut lp = duplicate_small_lp([1e-12, 1e-12]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kError);
        // 1.15.1 TestLpValidation.cpp: duplicate zeros are an error
        let mut lp = duplicate_small_lp([0.0, 0.0]);
        assert_eq!(assess_lp(&mut lp, &options), HighsStatus::kError);
    }

    /// Column-wise matrix: column 0 = (row 0: 1e-12, row 0: 4), column 1 =
    /// (row 2: 5, row 1: 1e-12, row 2: 0). Returns (status, start, index, value).
    fn mixed_matrix(
        h16_assess_small: bool,
    ) -> (HighsStatus, Vec<HighsInt>, Vec<HighsInt>, Vec<f64>) {
        let mut start: Vec<HighsInt> = vec![0, 2, 5];
        let mut p_end: Vec<HighsInt> = Vec::new();
        let mut index: Vec<HighsInt> = vec![0, 0, 2, 1, 2];
        let mut value: Vec<f64> = vec![1e-12, 4.0, 5.0, 1e-12, 0.0];
        let status = assess_matrix(
            &HighsLogOptions::default(),
            "LP",
            3,
            2,
            false,
            &mut start,
            &mut p_end,
            &mut index,
            &mut value,
            1e-9,
            kHighsInf,
            false,
            h16_assess_small,
        );
        (status, start, index, value)
    }

    #[test]
    fn h16_assess_small_matrix_result_on() {
        let (status, start, index, value) = mixed_matrix(true);
        assert_eq!(status, HighsStatus::kWarning);
        assert_eq!(start, vec![0, 1, 2]);
        assert_eq!(&index[..2], &[0, 2]);
        assert_eq!(&value[..2], &[4.0, 5.0]);
    }

    #[test]
    fn h16_assess_small_matrix_result_off() {
        // 1.15.1: row 0 of column 0 is a duplicate of the small entry
        let (status, _, _, _) = mixed_matrix(false);
        assert_eq!(status, HighsStatus::kError);
    }

    #[test]
    fn h16_assess_small_add_rows_explicit_zero_on() {
        // addRows with an explicit zero before a nonzero of the same column:
        // the zero is dropped, leaving one nonzero with no error
        let mut highs = Highs::default();
        highs.options_.output_flag = false;
        highs.options_.h16_assess_small = true;
        highs.model_lp_.num_col_ = 1;
        // a dimensionally consistent one-column LP (addRows asserts lpDimensionsOk)
        highs.model_lp_.col_cost_ = vec![0.0];
        highs.model_lp_.col_lower_ = vec![0.0];
        highs.model_lp_.col_upper_ = vec![1.0];
        highs.model_lp_.a_matrix_.format_ = MatrixFormat::kColwise;
        highs.model_lp_.a_matrix_.num_col_ = 1;
        highs.model_lp_.a_matrix_.num_row_ = 0;
        highs.model_lp_.a_matrix_.start_ = vec![0, 0];
        let status =
            highs.add_rows_interface(1, &[-kHighsInf], &[1.0], 2, &[0], &[0, 0], &[0.0, 2.0]);
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(highs.model_lp_.num_row_, 1);
        assert_eq!(highs.model_lp_.a_matrix_.num_nz(), 1);
    }

    #[test]
    fn h16_assess_small_add_rows_explicit_zero_off() {
        // 1.15.1: the explicit zero is recorded first, so the nonzero is a duplicate
        let mut highs = Highs::default();
        highs.options_.output_flag = false;
        highs.options_.h16_assess_small = false;
        highs.model_lp_.num_col_ = 1;
        // a dimensionally consistent one-column LP (addRows asserts lpDimensionsOk)
        highs.model_lp_.col_cost_ = vec![0.0];
        highs.model_lp_.col_lower_ = vec![0.0];
        highs.model_lp_.col_upper_ = vec![1.0];
        highs.model_lp_.a_matrix_.format_ = MatrixFormat::kColwise;
        highs.model_lp_.a_matrix_.num_col_ = 1;
        highs.model_lp_.a_matrix_.num_row_ = 0;
        highs.model_lp_.a_matrix_.start_ = vec![0, 0];
        let status =
            highs.add_rows_interface(1, &[-kHighsInf], &[1.0], 2, &[0], &[0, 0], &[0.0, 2.0]);
        assert_eq!(status, HighsStatus::kError);
    }
}
