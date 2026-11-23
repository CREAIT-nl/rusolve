//! HiGHS source: util/HighsSparseMatrix.cpp (lines 39-251, 542-592, 1020-1147)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Construction, format queries/conversion, column/row extraction, slicing.

use super::*;
use crate::prelude::*;

impl HighsSparseMatrix {
    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::equivalent (lines 39-51, 13 lines)
    pub fn equivalent(&self, matrix: &HighsSparseMatrix) -> bool {
        // Use this to test when the matrices may be in opposite format
        if *self == *matrix {
            return true;
        }
        // Creat a local copy in of matrix in the same format as this and
        // test for it being equal to this
        let mut lc_matrix = matrix.clone();
        if self.is_colwise() {
            lc_matrix.ensure_colwise();
        } else {
            lc_matrix.ensure_rowwise();
        }
        *self == lc_matrix
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::clear (lines 53-62, 10 lines)
    pub fn clear(&mut self) {
        self.num_col_ = 0;
        self.num_row_ = 0;
        self.start_.clear();
        self.p_end_.clear();
        self.index_.clear();
        self.value_.clear();
        self.format_ = MatrixFormat::kColwise;
        self.start_.clear();
        self.start_.push(0);
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::exactResize (lines 64-80, 17 lines)
    pub fn exact_resize(&mut self) {
        if self.is_colwise() {
            self.start_.resize((self.num_col_ + 1) as usize, 0);
        } else {
            self.start_.resize((self.num_row_ + 1) as usize, 0);
        }
        let num_nz = if self.is_colwise() {
            self.start_[self.num_col_ as usize]
        } else {
            self.start_[self.num_row_ as usize]
        };
        if self.format_ == MatrixFormat::kRowwisePartitioned {
            self.p_end_.resize(self.num_row_ as usize, 0);
        } else {
            debug_assert!((self.p_end_.len() as HighsInt) == 0);
            self.p_end_.clear();
        }
        self.index_.resize(num_nz as usize, 0);
        self.value_.resize(num_nz as usize, 0.0);
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::isRowwise (lines 82-85, 4 lines)
    pub fn is_rowwise(&self) -> bool {
        self.format_ == MatrixFormat::kRowwise || self.format_ == MatrixFormat::kRowwisePartitioned
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::isColwise (lines 87-89, 3 lines)
    pub fn is_colwise(&self) -> bool {
        self.format_ == MatrixFormat::kColwise
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::numNz (lines 91-100, 10 lines)
    pub fn num_nz(&self) -> HighsInt {
        debug_assert!(self.format_ok());
        if self.is_colwise() {
            debug_assert!((self.start_.len() as HighsInt) > self.num_col_);
            self.start_[self.num_col_ as usize]
        } else {
            debug_assert!((self.start_.len() as HighsInt) > self.num_row_);
            self.start_[self.num_row_ as usize]
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::range (lines 102-109, 8 lines)
    pub fn range(&self, min_value: &mut f64, max_value: &mut f64) {
        debug_assert!(self.format_ok());
        for i_el in 0..self.num_nz() {
            let value = self.value_[i_el as usize].abs();
            *min_value = cpp_min(*min_value, value);
            *max_value = cpp_max(*max_value, value);
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::setFormat (lines 111-119, 9 lines)
    pub fn set_format(&mut self, desired_format: MatrixFormat) {
        debug_assert!(self.format_ok());
        if desired_format == MatrixFormat::kColwise {
            self.ensure_colwise();
        } else {
            self.ensure_rowwise();
        }
        debug_assert!(self.format_ == desired_format);
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::ensureColwise (lines 121-174, 54 lines)
    pub fn ensure_colwise(&mut self) {
        debug_assert!(self.format_ok());
        if self.is_colwise() {
            return;
        }
        let num_col = self.num_col_;
        let num_row = self.num_row_;
        let num_nz = self.num_nz();
        debug_assert!(num_nz >= 0);
        debug_assert!((self.index_.len() as HighsInt) >= num_nz);
        debug_assert!((self.value_.len() as HighsInt) >= num_nz);
        if num_nz == 0 {
            // Empty matrix, so just ensure that there are enough zero starts
            // for the new orientation
            self.start_.clear();
            self.start_.resize((num_col + 1) as usize, 0);
            self.index_.clear();
            self.value_.clear();
        } else {
            // Matrix is non-empty, so transpose it
            //
            // Take a copy of the current matrix - that is rowwise - so that
            // the current matrix is filled colwise
            let ar_start = self.start_.clone();
            let ar_index = self.index_.clone();
            let ar_value = self.value_.clone();
            self.start_.resize((num_col + 1) as usize, 0);
            self.index_.resize(num_nz as usize, 0);
            self.value_.resize(num_nz as usize, 0.0);
            let mut a_length: Vec<HighsInt> = vec![0; num_col as usize];
            for i_el in ar_start[0]..num_nz {
                a_length[ar_index[i_el as usize] as usize] += 1;
            }
            self.start_[0] = 0;
            for i_col in 0..num_col {
                self.start_[(i_col + 1) as usize] =
                    self.start_[i_col as usize] + a_length[i_col as usize];
            }
            for i_row in 0..num_row {
                for i_el in ar_start[i_row as usize]..ar_start[(i_row + 1) as usize] {
                    let i_col = ar_index[i_el as usize];
                    let i_col_el = self.start_[i_col as usize];
                    self.index_[i_col_el as usize] = i_row;
                    self.value_[i_col_el as usize] = ar_value[i_el as usize];
                    self.start_[i_col as usize] += 1;
                }
            }
            self.start_[0] = 0;
            for i_col in 0..num_col {
                self.start_[(i_col + 1) as usize] =
                    self.start_[i_col as usize] + a_length[i_col as usize];
            }
            debug_assert!(self.start_[num_col as usize] == num_nz);
        }
        self.format_ = MatrixFormat::kColwise;
        debug_assert!((self.start_.len() as HighsInt) > num_col);
        debug_assert!(self.num_nz() >= 0);
        debug_assert!((self.index_.len() as HighsInt) >= self.num_nz());
        debug_assert!((self.value_.len() as HighsInt) >= self.num_nz());
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::ensureRowwise (lines 176-228, 53 lines)
    pub fn ensure_rowwise(&mut self) {
        debug_assert!(self.format_ok());
        if self.is_rowwise() {
            return;
        }
        let num_col = self.num_col_;
        let num_row = self.num_row_;
        let num_nz = self.num_nz();
        debug_assert!(num_nz >= 0);
        debug_assert!((self.index_.len() as HighsInt) >= num_nz);
        debug_assert!((self.value_.len() as HighsInt) >= num_nz);
        if num_nz == 0 {
            // Empty matrix, so just ensure that there are enough zero starts
            // for the new orientation
            self.start_.clear();
            self.start_.resize((num_row + 1) as usize, 0);
            self.index_.clear();
            self.value_.clear();
        } else {
            // Matrix is non-empty, so transpose it
            //
            // Take a copy of the current matrix - that is colwise - so that
            // the current matrix is filled rowwise
            let a_start = self.start_.clone();
            let a_index = self.index_.clone();
            let a_value = self.value_.clone();
            self.start_.resize((num_row + 1) as usize, 0);
            self.index_.resize(num_nz as usize, 0);
            self.value_.resize(num_nz as usize, 0.0);
            let mut ar_length: Vec<HighsInt> = vec![0; num_row as usize];
            for i_el in a_start[0]..num_nz {
                ar_length[a_index[i_el as usize] as usize] += 1;
            }
            self.start_[0] = 0;
            for i_row in 0..num_row {
                self.start_[(i_row + 1) as usize] =
                    self.start_[i_row as usize] + ar_length[i_row as usize];
            }
            for i_col in 0..num_col {
                for i_el in a_start[i_col as usize]..a_start[(i_col + 1) as usize] {
                    let i_row = a_index[i_el as usize];
                    let i_row_el = self.start_[i_row as usize];
                    self.index_[i_row_el as usize] = i_col;
                    self.value_[i_row_el as usize] = a_value[i_el as usize];
                    self.start_[i_row as usize] += 1;
                }
            }
            self.start_[0] = 0;
            for i_row in 0..num_row {
                self.start_[(i_row + 1) as usize] =
                    self.start_[i_row as usize] + ar_length[i_row as usize];
            }
            debug_assert!(self.start_[num_row as usize] == num_nz);
        }
        self.format_ = MatrixFormat::kRowwise;
        debug_assert!((self.start_.len() as HighsInt) > num_row);
        debug_assert!(self.num_nz() >= 0);
        debug_assert!((self.index_.len() as HighsInt) >= self.num_nz());
        debug_assert!((self.value_.len() as HighsInt) >= self.num_nz());
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::addVec (lines 230-251, 22 lines)
    /// C++ default arguments: multiple = 1
    pub fn add_vec(&mut self, num_nz: HighsInt, index: &[HighsInt], value: &[f64], multiple: f64) {
        let num_vec = if self.is_colwise() {
            self.num_col_
        } else {
            self.num_row_
        };
        debug_assert!((self.start_.len() as HighsInt) == num_vec + 1);
        debug_assert!((self.index_.len() as HighsInt) == self.num_nz());
        debug_assert!((self.value_.len() as HighsInt) == self.num_nz());
        for i_el in 0..num_nz {
            self.index_.push(index[i_el as usize]);
            self.value_.push(multiple * value[i_el as usize]);
        }
        let new_start = self.start_[num_vec as usize] + num_nz;
        self.start_.push(new_start);
        if self.is_colwise() {
            self.num_col_ += 1;
        } else {
            self.num_row_ += 1;
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::getCol (lines 542-566, 25 lines)
    /// C++ callers may pass `nullptr` for `index`/`value` (then only `num_nz` is set).
    pub fn get_col(
        &self,
        i_col: HighsInt,
        num_nz: &mut HighsInt,
        index: Option<&mut [HighsInt]>,
        value: Option<&mut [f64]>,
    ) {
        let mut index = index;
        let mut value = value;
        debug_assert!(i_col >= 0 && i_col < self.num_col_);
        *num_nz = 0;
        if self.is_colwise() {
            for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                if let Some(ix) = index.as_deref_mut() {
                    ix[*num_nz as usize] = self.index_[i_el as usize];
                }
                if let Some(vx) = value.as_deref_mut() {
                    vx[*num_nz as usize] = self.value_[i_el as usize];
                }
                *num_nz += 1;
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    if self.index_[i_el as usize] == i_col {
                        if let Some(ix) = index.as_deref_mut() {
                            ix[*num_nz as usize] = i_row;
                        }
                        if let Some(vx) = value.as_deref_mut() {
                            vx[*num_nz as usize] = self.value_[i_el as usize];
                        }
                        *num_nz += 1;
                        break;
                    }
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::getRow (lines 568-592, 25 lines)
    /// C++ callers may pass `nullptr` for `index`/`value` (then only `num_nz` is set).
    pub fn get_row(
        &self,
        i_row: HighsInt,
        num_nz: &mut HighsInt,
        index: Option<&mut [HighsInt]>,
        value: Option<&mut [f64]>,
    ) {
        let mut index = index;
        let mut value = value;
        debug_assert!(i_row >= 0 && i_row < self.num_row_);
        *num_nz = 0;
        if self.is_rowwise() {
            for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                if let Some(ix) = index.as_deref_mut() {
                    ix[*num_nz as usize] = self.index_[i_el as usize];
                }
                if let Some(vx) = value.as_deref_mut() {
                    vx[*num_nz as usize] = self.value_[i_el as usize];
                }
                *num_nz += 1;
            }
        } else {
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    if self.index_[i_el as usize] == i_row {
                        if let Some(ix) = index.as_deref_mut() {
                            ix[*num_nz as usize] = i_col;
                        }
                        if let Some(vx) = value.as_deref_mut() {
                            vx[*num_nz as usize] = self.value_[i_el as usize];
                        }
                        *num_nz += 1;
                        break;
                    }
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::createSlice (lines 1020-1049, 30 lines)
    pub fn create_slice(
        &mut self,
        matrix: &HighsSparseMatrix,
        from_col: HighsInt,
        to_col: HighsInt,
    ) {
        debug_assert!(matrix.format_ok());
        debug_assert!(matrix.is_colwise());
        debug_assert!(self.format_ok());
        let num_row = matrix.num_row_;
        let a_start = &matrix.start_;
        let a_index = &matrix.index_;
        let a_value = &matrix.value_;
        let slice_num_col = to_col + 1 - from_col;
        let slice_num_nz = a_start[(to_col + 1) as usize] - a_start[from_col as usize];
        self.start_.resize((slice_num_col + 1) as usize, 0);
        self.index_.resize(slice_num_nz as usize, 0);
        self.value_.resize(slice_num_nz as usize, 0.0);
        let from_col_start = a_start[from_col as usize];
        for i_col in from_col..(to_col + 1) {
            self.start_[(i_col - from_col) as usize] = a_start[i_col as usize] - from_col_start;
        }
        self.start_[slice_num_col as usize] = slice_num_nz;
        for i_el in a_start[from_col as usize]..a_start[(to_col + 1) as usize] {
            self.index_[(i_el - from_col_start) as usize] = a_index[i_el as usize];
            self.value_[(i_el - from_col_start) as usize] = a_value[i_el as usize];
        }
        self.num_col_ = slice_num_col;
        self.num_row_ = num_row;
        self.format_ = MatrixFormat::kColwise;
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::createRowwise (lines 1051-1098, 48 lines)
    pub fn create_rowwise(&mut self, matrix: &HighsSparseMatrix) {
        debug_assert!(matrix.format_ok());
        debug_assert!(matrix.is_colwise());
        debug_assert!(self.format_ok());

        let num_col = matrix.num_col_;
        let num_row = matrix.num_row_;
        let num_nz = matrix.num_nz();
        let a_start = &matrix.start_;
        let a_index = &matrix.index_;
        let a_value = &matrix.value_;

        // Use ar_end to compute lengths, which are then transformed into
        // the ends of the inserted entries
        let mut ar_end: Vec<HighsInt> = Vec::new();
        self.start_.resize((num_row + 1) as usize, 0);
        ar_end.resize(num_row as usize, 0);
        // Count the nonzeros in each row
        for i_col in 0..num_col {
            for i_el in a_start[i_col as usize]..a_start[(i_col + 1) as usize] {
                let i_row = a_index[i_el as usize];
                ar_end[i_row as usize] += 1;
            }
        }
        // Compute the starts and turn the lengths into ends
        self.start_[0] = 0;
        for i_row in 0..num_row {
            self.start_[(i_row + 1) as usize] =
                self.start_[i_row as usize] + ar_end[i_row as usize];
            ar_end[i_row as usize] = self.start_[i_row as usize];
        }
        self.index_.resize(num_nz as usize, 0);
        self.value_.resize(num_nz as usize, 0.0);
        // Insert the entries
        for i_col in 0..num_col {
            for i_el in a_start[i_col as usize]..a_start[(i_col + 1) as usize] {
                let i_row = a_index[i_el as usize];
                let i_to_el = ar_end[i_row as usize];
                ar_end[i_row as usize] += 1;
                self.index_[i_to_el as usize] = i_col;
                self.value_[i_to_el as usize] = a_value[i_el as usize];
            }
        }
        self.format_ = MatrixFormat::kRowwise;
        self.num_col_ = num_col;
        self.num_row_ = num_row;
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::createColwise (lines 1100-1147, 48 lines)
    pub fn create_colwise(&mut self, matrix: &HighsSparseMatrix) {
        debug_assert!(matrix.format_ok());
        debug_assert!(matrix.is_rowwise());
        debug_assert!(self.format_ok());

        let num_col = matrix.num_col_;
        let num_row = matrix.num_row_;
        let num_nz = matrix.num_nz();
        let ar_start = &matrix.start_;
        let ar_index = &matrix.index_;
        let ar_value = &matrix.value_;

        // Use a_end to compute lengths, which are then transformed into
        // the ends of the inserted entries
        let mut a_end: Vec<HighsInt> = Vec::new();
        self.start_.resize((num_col + 1) as usize, 0);
        a_end.resize(num_col as usize, 0);
        // Count the nonzeros in each col
        for i_row in 0..num_row {
            for i_el in ar_start[i_row as usize]..ar_start[(i_row + 1) as usize] {
                let i_col = ar_index[i_el as usize];
                a_end[i_col as usize] += 1;
            }
        }
        // Compute the starts and turn the lengths into ends
        self.start_[0] = 0;
        for i_col in 0..num_col {
            self.start_[(i_col + 1) as usize] = self.start_[i_col as usize] + a_end[i_col as usize];
            a_end[i_col as usize] = self.start_[i_col as usize];
        }
        self.index_.resize(num_nz as usize, 0);
        self.value_.resize(num_nz as usize, 0.0);
        // Insert the entries
        for i_row in 0..num_row {
            for i_el in ar_start[i_row as usize]..ar_start[(i_row + 1) as usize] {
                let i_col = ar_index[i_el as usize];
                let i_to_el = a_end[i_col as usize];
                a_end[i_col as usize] += 1;
                self.index_[i_to_el as usize] = i_row;
                self.value_[i_to_el as usize] = ar_value[i_el as usize];
            }
        }
        self.format_ = MatrixFormat::kColwise;
        self.num_col_ = num_col;
        self.num_row_ = num_row;
    }

    /// HiGHS: HighsSparseMatrix.h::HighsSparseMatrix::formatOk (inline)
    pub fn format_ok(&self) -> bool {
        self.is_colwise() || self.is_rowwise()
    }
}
