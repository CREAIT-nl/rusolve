//! HiGHS source: util/HighsSparseMatrix.cpp (lines 794-813, 826-830)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Model validation: `assess` (called by `assessLp`) and `hasLargeValue`
//! (called by `Highs::calledOptimizeModel`).

use super::*;
use crate::prelude::*;
use crate::util::highs_matrix_utils::assess_matrix;

impl HighsSparseMatrix {
    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::assess (lines 794-813, 20 lines)
    /// HiGHS 1.16 #3169 (h16_assess_small): `h16_assess_small` is option plumbing only, passed to assess_matrix
    pub fn assess(
        &mut self,
        log_options: &HighsLogOptions,
        matrix_name: &str,
        small_matrix_value: f64,
        large_matrix_value: f64,
        h16_assess_small: bool,
    ) -> HighsStatus {
        debug_assert!(self.format_ok());
        // Identify main dimensions
        let vec_dim: HighsInt;
        let num_vec: HighsInt;
        if self.is_colwise() {
            vec_dim = self.num_row_;
            num_vec = self.num_col_;
        } else {
            vec_dim = self.num_col_;
            num_vec = self.num_row_;
        }
        let partitioned = self.format_ == MatrixFormat::kRowwisePartitioned;
        assess_matrix(
            log_options,
            matrix_name,
            vec_dim,
            num_vec,
            partitioned,
            &mut self.start_,
            &mut self.p_end_,
            &mut self.index_,
            &mut self.value_,
            small_matrix_value,
            large_matrix_value,
            false,
            h16_assess_small,
        )
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::hasLargeValue (lines 826-830, 5 lines)
    pub fn has_large_value(&self, large_matrix_value: f64) -> bool {
        for i_el in 0..self.num_nz() {
            if self.value_[i_el as usize].abs() >= large_matrix_value {
                return true;
            }
        }
        false
    }
}
