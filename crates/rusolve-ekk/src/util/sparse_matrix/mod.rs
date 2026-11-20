//! HiGHS source: util/HighsSparseMatrix.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Struct, constants and `PartialEq`. Method bodies are split over `basic.rs`,
//! `modify.rs` and `product.rs` (by HighsSparseMatrix.cpp line ranges). The `assess*` family
//! and `hasLargeValue` (model validation, Highs.cpp path) are out of scope.

use crate::prelude::*;

pub mod assess;
pub mod basic;
pub mod modify;
pub mod product;

/// HiGHS: HighsSparseMatrix.h kHyperPriceDensity
pub const kHyperPriceDensity: f64 = 0.1;
/// HiGHS: HighsSparseMatrix.h kDebugReportOff
pub const kDebugReportOff: HighsInt = -2;
/// HiGHS: HighsSparseMatrix.h kDebugReportAll
pub const kDebugReportAll: HighsInt = -1;

/// HiGHS: HighsSparseMatrix.h::HighsSparseMatrix
///
/// `Default` equals the C++ constructor (`clear()`): colwise, 0x0, `start_ = [0]`.
#[derive(Clone, Debug)]
pub struct HighsSparseMatrix {
    pub format_: MatrixFormat,
    pub num_col_: HighsInt,
    pub num_row_: HighsInt,
    pub start_: Vec<HighsInt>,
    pub p_end_: Vec<HighsInt>,
    pub index_: Vec<HighsInt>,
    pub value_: Vec<f64>,
}

impl Default for HighsSparseMatrix {
    fn default() -> Self {
        Self {
            format_: MatrixFormat::kColwise,
            num_col_: 0,
            num_row_: 0,
            start_: vec![0],
            p_end_: Vec::new(),
            index_: Vec::new(),
            value_: Vec::new(),
        }
    }
}

/// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::operator== (lines 28-37, 10 lines).
/// Compares format_, num_col_, num_row_, start_, index_, value_ (not p_end_).
impl PartialEq for HighsSparseMatrix {
    fn eq(&self, matrix: &HighsSparseMatrix) -> bool {
        let mut equal: bool = true;
        equal = self.format_ == matrix.format_ && equal;
        equal = self.num_col_ == matrix.num_col_ && equal;
        equal = self.num_row_ == matrix.num_row_ && equal;
        equal = self.start_ == matrix.start_ && equal;
        equal = self.index_ == matrix.index_ && equal;
        equal = self.value_ == matrix.value_ && equal;
        equal
    }
}
