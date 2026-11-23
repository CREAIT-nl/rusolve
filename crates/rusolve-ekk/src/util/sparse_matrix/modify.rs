//! HiGHS source: util/HighsSparseMatrix.cpp (lines 253-540, 594-734, 832-1018)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Adding/deleting columns and rows, scaling.

use super::*;
use crate::prelude::*;

impl HighsSparseMatrix {
    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::addCols (lines 253-375, 123 lines)
    /// C++ default arguments: in_partition = NULL
    /// C++ takes `new_cols` by value (a copy); a shared reference is equivalent.
    pub fn add_cols(&mut self, new_cols: &HighsSparseMatrix, in_partition: Option<&[i8]>) {
        debug_assert!(new_cols.is_colwise());
        let num_new_col: HighsInt = new_cols.num_col_;
        let num_new_nz: HighsInt = new_cols.num_nz();
        let new_matrix_start = &new_cols.start_;
        let new_matrix_index = &new_cols.index_;
        let new_matrix_value = &new_cols.value_;

        debug_assert!(self.format_ok());
        // Adding columns to a row-wise partitioned matrix needs the
        // partition information
        let partitioned = self.format_ == MatrixFormat::kRowwisePartitioned;
        // Cannot handle the row-wise partitioned case
        debug_assert!(!partitioned);
        if partitioned {
            debug_assert!(in_partition.is_some());
        }
        debug_assert!(num_new_col >= 0);
        debug_assert!(num_new_nz >= 0);
        if num_new_col == 0 {
            // No columns are being added, so check that no nonzeros are being
            // added
            debug_assert!(num_new_nz == 0);
            return;
        }
        // Adding a positive number of columns to a matrix
        if num_new_nz != 0 {
            // Nonzeros are being added, so ensure that non-null data are
            // being passed
            debug_assert!(!new_matrix_start.is_empty());
            debug_assert!(!new_matrix_index.is_empty());
            debug_assert!(!new_matrix_value.is_empty());
        }
        let num_col: HighsInt = self.num_col_;
        let num_row: HighsInt = self.num_row_;
        let num_nz: HighsInt = self.num_nz();
        // Check that nonzeros aren't being appended to a matrix with no rows
        debug_assert!(num_new_nz <= 0 || num_row > 0);

        // If matrix is currently a standard row-wise matrix and there are
        // more new nonzeros than current nonzeros so flip column-wise
        if self.format_ == MatrixFormat::kRowwise && num_new_nz > num_nz {
            self.ensure_colwise();
        }

        // Determine the new number of columns and nonzeros in the matrix
        let new_num_col: HighsInt = num_col + num_new_col;
        let new_num_nz: HighsInt = num_nz + num_new_nz;

        if self.is_colwise() {
            // Matrix is column-wise
            self.start_.resize((new_num_col + 1) as usize, 0);
            // Append the starts of the new columns
            if num_new_nz != 0 {
                // Nontrivial number of nonzeros being added, so use new_matrix_start
                for i_new_col in 0..num_new_col {
                    self.start_[(num_col + i_new_col) as usize] =
                        num_nz + new_matrix_start[i_new_col as usize];
                }
            } else {
                // No nonzeros being added, so new_matrix_start may be null, but entries
                // of zero are implied.
                for i_new_col in 0..num_new_col {
                    self.start_[(num_col + i_new_col) as usize] = num_nz;
                }
            }
            self.start_[(num_col + num_new_col) as usize] = new_num_nz;
            // Update the number of columns
            self.num_col_ += num_new_col;
            // If no nonzeros are being added then there's nothing else to do
            if num_new_nz <= 0 {
                return;
            }
            // Adding a non-trivial matrix: resize the column-wise matrix arrays
            // accordingly
            self.index_.resize(new_num_nz as usize, 0);
            self.value_.resize(new_num_nz as usize, 0.0);
            // Copy in the new indices and values
            for i_el in 0..num_new_nz {
                self.index_[(num_nz + i_el) as usize] = new_matrix_index[i_el as usize];
                self.value_[(num_nz + i_el) as usize] = new_matrix_value[i_el as usize];
            }
        } else {
            // Matrix is row-wise
            if num_new_nz != 0 {
                // Adding a positive number of nonzeros
                self.index_.resize(new_num_nz as usize, 0);
                self.value_.resize(new_num_nz as usize, 0.0);
                // Determine the row lengths of the new columns being added
                let mut new_row_length: Vec<HighsInt> = vec![0; num_row as usize];
                for i_el in 0..num_new_nz {
                    new_row_length[new_matrix_index[i_el as usize] as usize] += 1;
                }
                // Now shift the indices and values to make space
                let mut entry_offset: HighsInt = num_new_nz;
                let mut to_original_el: HighsInt = self.start_[num_row as usize];
                self.start_[num_row as usize] = new_num_nz;
                for i_row in (0..num_row).rev() {
                    entry_offset -= new_row_length[i_row as usize];
                    let from_original_el: HighsInt = self.start_[i_row as usize];
                    // Can now use this new_row_length to store the start for the
                    // new entries
                    new_row_length[i_row as usize] = to_original_el + entry_offset;
                    for i_el in (from_original_el..to_original_el).rev() {
                        self.index_[(i_el + entry_offset) as usize] = self.index_[i_el as usize];
                        self.value_[(i_el + entry_offset) as usize] = self.value_[i_el as usize];
                    }
                    to_original_el = from_original_el;
                    self.start_[i_row as usize] = entry_offset + from_original_el;
                }
                // Now insert the indices and values for the new columns
                for i_col in 0..num_new_col {
                    for i_el in
                        new_matrix_start[i_col as usize]..new_matrix_start[(i_col + 1) as usize]
                    {
                        let i_row: HighsInt = new_matrix_index[i_el as usize];
                        self.index_[new_row_length[i_row as usize] as usize] = num_col + i_col;
                        self.value_[new_row_length[i_row as usize] as usize] =
                            new_matrix_value[i_el as usize];
                        new_row_length[i_row as usize] += 1;
                    }
                }
            }
            // Have to increase the number of columns, even if no nonzeros are being
            // added
            self.num_col_ += num_new_col;
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::addRows (lines 377-540, 164 lines)
    /// C++ default arguments: in_partition = NULL
    /// C++ takes `new_rows` by value (a copy); a shared reference is equivalent.
    pub fn add_rows(&mut self, new_rows: &HighsSparseMatrix, in_partition: Option<&[i8]>) {
        debug_assert!(new_rows.is_rowwise());
        let num_new_row: HighsInt = new_rows.num_row_;
        let num_new_nz: HighsInt = new_rows.num_nz();
        let new_matrix_start = &new_rows.start_;
        let new_matrix_index = &new_rows.index_;
        let new_matrix_value = &new_rows.value_;

        debug_assert!(self.format_ok());
        // Adding rows to a row-wise partitioned matrix needs the
        // partition information
        let partitioned = self.format_ == MatrixFormat::kRowwisePartitioned;
        if partitioned {
            debug_assert!(!partitioned);
            debug_assert!(in_partition.is_some());
        }
        debug_assert!(num_new_row >= 0);
        debug_assert!(num_new_nz >= 0);
        if num_new_row == 0 {
            // No rows are being added, so check that no nonzeros are being
            // added
            debug_assert!(num_new_nz == 0);
            return;
        }
        // Adding a positive number of rows to a matrix
        if num_new_nz != 0 {
            // Nonzeros are being added, so ensure that non-null data are
            // being passed
            debug_assert!(!new_matrix_start.is_empty());
            debug_assert!(!new_matrix_index.is_empty());
            debug_assert!(!new_matrix_value.is_empty());
        }
        let num_col: HighsInt = self.num_col_;
        let num_row: HighsInt = self.num_row_;
        let num_nz: HighsInt = self.num_nz();
        // Check that nonzeros aren't being appended to a matrix with no columns
        debug_assert!(num_new_nz <= 0 || num_col > 0);

        if self.is_colwise() {
            // Matrix is currently a standard col-wise matrix, so flip
            // row-wise if there are more new nonzeros than current nonzeros
            if num_new_nz > num_nz {
                self.ensure_rowwise();
            }
        }
        // Determine the new number of rows and nonzeros in the matrix
        let new_num_nz: HighsInt = num_nz + num_new_nz;
        let new_num_row: HighsInt = num_row + num_new_row;

        if self.is_rowwise() {
            // Matrix is row-wise
            self.start_.resize((new_num_row + 1) as usize, 0);
            // Append the starts of the new rows
            if num_new_nz != 0 {
                // Nontrivial number of nonzeros being added, so use new_matrix_start
                for i_new_row in 0..num_new_row {
                    self.start_[(num_row + i_new_row) as usize] =
                        num_nz + new_matrix_start[i_new_row as usize];
                }
            } else {
                // No nonzeros being added, so new_matrix_start may be NULL, but entries
                // of zero are implied.
                for i_new_row in 0..num_new_row {
                    self.start_[(num_row + i_new_row) as usize] = num_nz;
                }
            }
            self.start_[new_num_row as usize] = new_num_nz;
            if num_new_nz > 0 {
                // Adding a non-trivial matrix: resize the matrix arrays accordingly
                self.index_.resize(new_num_nz as usize, 0);
                self.value_.resize(new_num_nz as usize, 0.0);
                // Copy in the new indices and values
                if partitioned {
                    // Insert the entries in the partition
                    debug_assert!(!partitioned);
                    let in_partition: &[i8] =
                        in_partition.expect("partitioned addRows needs in_partition");
                    for i_new_row in 0..num_new_row {
                        let i_row: HighsInt = num_row + i_new_row;
                        for i_new_el in new_matrix_start[i_new_row as usize]
                            ..new_matrix_start[(i_new_row + 1) as usize]
                        {
                            let i_col: HighsInt = new_matrix_index[i_new_el as usize];
                            if in_partition[i_col as usize] != 0 {
                                let i_el: HighsInt = self.start_[i_row as usize];
                                self.index_[i_el as usize] = new_matrix_index[i_new_el as usize];
                                self.value_[i_el as usize] = new_matrix_value[i_new_el as usize];
                                self.start_[i_row as usize] += 1;
                            }
                        }
                    }
                    // Use the incremented starts to initialise p_end, save these
                    // values and reset the starts
                    let mut save_p_end: Vec<HighsInt> = vec![0; num_new_row as usize];
                    for i_new_row in 0..num_new_row {
                        let i_row: HighsInt = num_row + i_new_row;
                        self.start_[i_row as usize] = num_nz + new_matrix_start[i_new_row as usize];
                        self.p_end_[i_row as usize] = self.start_[i_row as usize];
                        save_p_end[i_new_row as usize] = self.p_end_[i_row as usize];
                    }
                    // Insert the entries not in the partition
                    for i_new_row in 0..num_new_row {
                        let i_row: HighsInt = num_row + i_new_row;
                        for i_new_el in new_matrix_start[i_new_row as usize]
                            ..new_matrix_start[(i_new_row + 1) as usize]
                        {
                            let i_col: HighsInt = new_matrix_index[i_new_el as usize];
                            if in_partition[i_col as usize] == 0 {
                                let i_el: HighsInt = self.p_end_[i_row as usize];
                                self.index_[i_el as usize] = new_matrix_index[i_new_el as usize];
                                self.value_[i_el as usize] = new_matrix_value[i_new_el as usize];
                                self.p_end_[i_row as usize] += 1;
                            }
                        }
                    }
                    // Reset p_end using the saved values
                    for i_new_row in 0..num_new_row {
                        self.p_end_[(num_row + i_new_row) as usize] =
                            save_p_end[i_new_row as usize];
                    }
                } else {
                    for i_new_el in 0..num_new_nz {
                        self.index_[(num_nz + i_new_el) as usize] =
                            new_matrix_index[i_new_el as usize];
                        self.value_[(num_nz + i_new_el) as usize] =
                            new_matrix_value[i_new_el as usize];
                    }
                }
            }
        } else {
            // Storing the matrix column-wise, so have to insert the new rows
            debug_assert!(self.is_colwise());
            if num_new_nz != 0 {
                let mut length: Vec<HighsInt> = vec![0; num_col as usize];
                for i_el in 0..num_new_nz {
                    length[new_matrix_index[i_el as usize] as usize] += 1;
                }
                // Determine the new number of nonzeros and resize the column-wise matrix
                // arrays
                self.index_.resize(new_num_nz as usize, 0);
                self.value_.resize(new_num_nz as usize, 0.0);
                // Append the new rows
                // Shift the existing columns to make space for the new entries
                let mut new_i_el: HighsInt = new_num_nz;
                for i_col in (0..num_col).rev() {
                    let start_col_plus_1: HighsInt = new_i_el;
                    new_i_el -= length[i_col as usize];
                    let col_start: HighsInt = self.start_[i_col as usize];
                    let col_end: HighsInt = self.start_[(i_col + 1) as usize];
                    for i_el in (col_start..col_end).rev() {
                        new_i_el -= 1;
                        self.index_[new_i_el as usize] = self.index_[i_el as usize];
                        self.value_[new_i_el as usize] = self.value_[i_el as usize];
                    }
                    self.start_[(i_col + 1) as usize] = start_col_plus_1;
                }
                debug_assert!(new_i_el == 0);
                // Insert the new entries
                for i_new_row in 0..num_new_row {
                    let first_el: HighsInt = new_matrix_start[i_new_row as usize];
                    let last_el: HighsInt = if i_new_row < num_new_row - 1 {
                        new_matrix_start[(i_new_row + 1) as usize]
                    } else {
                        num_new_nz
                    };
                    for i_el in first_el..last_el {
                        let i_col: HighsInt = new_matrix_index[i_el as usize];
                        new_i_el = self.start_[(i_col + 1) as usize] - length[i_col as usize];
                        length[i_col as usize] -= 1;
                        self.index_[new_i_el as usize] = num_row + i_new_row;
                        self.value_[new_i_el as usize] = new_matrix_value[i_el as usize];
                    }
                }
            }
        }
        // Update the number of rows
        self.num_row_ += num_new_row;
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::deleteCols (lines 594-656, 63 lines)
    pub fn delete_cols(&mut self, index_collection: &HighsIndexCollection) {
        debug_assert!(self.format_ok());
        // Can't handle rowwise matrices yet
        debug_assert!(!self.is_rowwise());
        debug_assert!(ok(index_collection));
        let mut from_k: HighsInt = 0;
        let mut to_k: HighsInt = 0;
        limits(index_collection, &mut from_k, &mut to_k);
        if from_k > to_k {
            return;
        }

        let mut delete_from_col: HighsInt = 0;
        let mut delete_to_col: HighsInt = 0;
        let mut keep_from_col: HighsInt = 0;
        let mut keep_to_col: HighsInt = -1;
        let mut current_set_entry: HighsInt = 0;

        let col_dim: HighsInt = self.num_col_;
        let mut new_num_col: HighsInt = 0;
        let mut new_num_nz: HighsInt = 0;
        for k in from_k..=to_k {
            update_out_in_index(
                index_collection,
                &mut delete_from_col,
                &mut delete_to_col,
                &mut keep_from_col,
                &mut keep_to_col,
                &mut current_set_entry,
            );
            if k == from_k {
                // Account for the initial columns being kept
                new_num_col = delete_from_col;
                new_num_nz = self.start_[delete_from_col as usize];
            }
            // Ensure that the starts of the deleted columns are zeroed to
            // avoid redundant start information for columns whose indices
            // aren't used after the deletion takes place. In particular, if
            // all columns are deleted then something must be done to ensure
            // that the matrix isn't magically recreated by increasing the
            // number of columns from zero when there are no rows in the
            // matrix.
            for col in delete_from_col..=delete_to_col {
                self.start_[col as usize] = 0;
            }
            // Shift the starts - both in place and value - to account for the
            // columns and nonzeros removed
            let keep_from_el: HighsInt = self.start_[keep_from_col as usize];
            for col in keep_from_col..=keep_to_col {
                self.start_[new_num_col as usize] =
                    new_num_nz + self.start_[col as usize] - keep_from_el;
                new_num_col += 1;
            }
            for el in keep_from_el..self.start_[(keep_to_col + 1) as usize] {
                self.index_[new_num_nz as usize] = self.index_[el as usize];
                self.value_[new_num_nz as usize] = self.value_[el as usize];
                new_num_nz += 1;
            }
            if keep_to_col >= col_dim - 1 {
                break;
            }
        }
        // Ensure that the start of the spurious last column is zeroed so
        // that it doesn't give a positive number of matrix entries if the
        // number of columns in the matrix is increased when there are no
        // rows in the matrix.
        self.start_[self.num_col_ as usize] = 0;
        self.start_[new_num_col as usize] = new_num_nz;
        self.start_.resize((new_num_col + 1) as usize, 0);
        self.index_.resize(new_num_nz as usize, 0);
        self.value_.resize(new_num_nz as usize, 0.0);
        // Update the number of columns
        self.num_col_ = new_num_col;
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::deleteRows (lines 658-734, 77 lines)
    pub fn delete_rows(&mut self, index_collection: &HighsIndexCollection) {
        debug_assert!(self.format_ok());
        debug_assert!(ok(index_collection));
        let mut from_k: HighsInt = 0;
        let mut to_k: HighsInt = 0;
        limits(index_collection, &mut from_k, &mut to_k);
        if from_k > to_k {
            return;
        }

        let mut delete_from_row: HighsInt = 0;
        let mut delete_to_row: HighsInt = 0;
        let mut keep_from_row: HighsInt = 0;
        let row_dim: HighsInt = self.num_row_;
        // Both are set in the `!mask` branch below before first use (C++ sets them there).
        let mut keep_to_row: HighsInt;
        let mut current_set_entry: HighsInt;

        // Set up a row mask to indicate the new row index of kept rows and
        // -1 for deleted rows so that the kept entries in the column-wise
        // matrix can be identified and have their correct row index.
        let mut new_index: Vec<HighsInt> = vec![0; self.num_row_ as usize];
        let mut new_num_row: HighsInt = 0;
        let mask: bool = index_collection.is_mask_;
        let row_mask: &[HighsInt] = &index_collection.mask_;
        if !mask {
            keep_to_row = -1;
            current_set_entry = 0;
            for k in from_k..=to_k {
                update_out_in_index(
                    index_collection,
                    &mut delete_from_row,
                    &mut delete_to_row,
                    &mut keep_from_row,
                    &mut keep_to_row,
                    &mut current_set_entry,
                );
                if k == from_k {
                    // Account for any initial rows being kept
                    for row in 0..delete_from_row {
                        new_index[row as usize] = new_num_row;
                        new_num_row += 1;
                    }
                }
                for row in delete_from_row..=delete_to_row {
                    new_index[row as usize] = -1;
                }
                for row in keep_from_row..=keep_to_row {
                    new_index[row as usize] = new_num_row;
                    new_num_row += 1;
                }
                if keep_to_row >= row_dim - 1 {
                    break;
                }
            }
        } else {
            for row in 0..self.num_row_ {
                if row_mask[row as usize] != 0 {
                    new_index[row as usize] = -1;
                } else {
                    new_index[row as usize] = new_num_row;
                    new_num_row += 1;
                }
            }
        }
        let mut new_num_nz: HighsInt = 0;
        for col in 0..self.num_col_ {
            let from_el: HighsInt = self.start_[col as usize];
            self.start_[col as usize] = new_num_nz;
            for el in from_el..self.start_[(col + 1) as usize] {
                let row: HighsInt = self.index_[el as usize];
                let new_row: HighsInt = new_index[row as usize];
                if new_row >= 0 {
                    self.index_[new_num_nz as usize] = new_row;
                    self.value_[new_num_nz as usize] = self.value_[el as usize];
                    new_num_nz += 1;
                }
            }
        }
        self.start_[self.num_col_ as usize] = new_num_nz;
        self.start_.resize((self.num_col_ + 1) as usize, 0);
        self.index_.resize(new_num_nz as usize, 0);
        self.value_.resize(new_num_nz as usize, 0.0);
        // Update the number of rows
        self.num_row_ = new_num_row;
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::considerColScaling (lines 832-867, 36 lines)
    pub fn consider_col_scaling(
        &mut self,
        max_scale_factor_exponent: HighsInt,
        col_scale: &mut [f64],
    ) {
        let log2: f64 = (2.0f64).ln();
        let max_allow_scale: f64 = (2.0f64).powf(max_scale_factor_exponent as f64);
        let min_allow_scale: f64 = 1.0 / max_allow_scale;

        let min_allow_col_scale: f64 = min_allow_scale;
        let max_allow_col_scale: f64 = max_allow_scale;

        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                let mut col_max_value: f64 = 0.0;
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    col_max_value = cpp_max(self.value_[i_el as usize].abs(), col_max_value);
                }
                if col_max_value != 0.0 {
                    let mut col_scale_value: f64 = 1.0 / col_max_value;
                    // Convert the col scale factor to the nearest power of two, and
                    // ensure that it is not excessively large or small
                    col_scale_value = (2.0f64).powf((col_scale_value.ln() / log2 + 0.5).floor());
                    col_scale_value = cpp_min(
                        cpp_max(min_allow_col_scale, col_scale_value),
                        max_allow_col_scale,
                    );
                    col_scale[i_col as usize] = col_scale_value;
                    // Scale the column
                    for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                        self.value_[i_el as usize] *= col_scale[i_col as usize];
                    }
                } else {
                    // Empty column
                    col_scale[i_col as usize] = 1.0;
                }
            }
        } else {
            debug_assert!(self.is_colwise());
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::considerRowScaling (lines 869-904, 36 lines)
    pub fn consider_row_scaling(
        &mut self,
        max_scale_factor_exponent: HighsInt,
        row_scale: &mut [f64],
    ) {
        let log2: f64 = (2.0f64).ln();
        let max_allow_scale: f64 = (2.0f64).powf(max_scale_factor_exponent as f64);
        let min_allow_scale: f64 = 1.0 / max_allow_scale;

        let min_allow_row_scale: f64 = min_allow_scale;
        let max_allow_row_scale: f64 = max_allow_scale;

        if self.is_rowwise() {
            for i_row in 0..self.num_row_ {
                let mut row_max_value: f64 = 0.0;
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    row_max_value = cpp_max(self.value_[i_el as usize].abs(), row_max_value);
                }
                if row_max_value != 0.0 {
                    let mut row_scale_value: f64 = 1.0 / row_max_value;
                    // Convert the row scale factor to the nearest power of two, and
                    // ensure that it is not excessively large or small
                    row_scale_value = (2.0f64).powf((row_scale_value.ln() / log2 + 0.5).floor());
                    row_scale_value = cpp_min(
                        cpp_max(min_allow_row_scale, row_scale_value),
                        max_allow_row_scale,
                    );
                    row_scale[i_row as usize] = row_scale_value;
                    // Scale the row
                    for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                        self.value_[i_el as usize] *= row_scale[i_row as usize];
                    }
                } else {
                    // Empty row
                    row_scale[i_row as usize] = 1.0;
                }
            }
        } else {
            debug_assert!(self.is_rowwise());
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::scaleCol (lines 906-923, 18 lines)
    pub fn scale_col(&mut self, col: HighsInt, col_scale: f64) {
        debug_assert!(self.format_ok());
        debug_assert!(col >= 0);
        debug_assert!(col < self.num_col_);
        debug_assert!(col_scale != 0.0);

        if self.is_colwise() {
            for i_el in self.start_[col as usize]..self.start_[(col + 1) as usize] {
                self.value_[i_el as usize] *= col_scale;
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    if self.index_[i_el as usize] == col {
                        self.value_[i_el as usize] *= col_scale;
                    }
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::scaleRow (lines 925-942, 18 lines)
    pub fn scale_row(&mut self, row: HighsInt, row_scale: f64) {
        debug_assert!(self.format_ok());
        debug_assert!(row >= 0);
        debug_assert!(row < self.num_row_);
        debug_assert!(row_scale != 0.0);

        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    if self.index_[i_el as usize] == row {
                        self.value_[i_el as usize] *= row_scale;
                    }
                }
            }
        } else {
            for i_el in self.start_[row as usize]..self.start_[(row + 1) as usize] {
                self.value_[i_el as usize] *= row_scale;
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::applyScale (lines 944-963, 20 lines)
    pub fn apply_scale(&mut self, scale: &HighsScale) {
        debug_assert!(self.format_ok());
        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    let i_row: HighsInt = self.index_[i_el as usize];
                    self.value_[i_el as usize] *=
                        scale.col[i_col as usize] * scale.row[i_row as usize];
                }
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    let i_col: HighsInt = self.index_[i_el as usize];
                    self.value_[i_el as usize] *=
                        scale.col[i_col as usize] * scale.row[i_row as usize];
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::applyColScale (lines 965-980, 16 lines)
    pub fn apply_col_scale(&mut self, scale: &HighsScale) {
        debug_assert!(self.format_ok());
        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    self.value_[i_el as usize] *= scale.col[i_col as usize];
                }
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    self.value_[i_el as usize] *= scale.col[self.index_[i_el as usize] as usize];
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::applyRowScale (lines 982-997, 16 lines)
    pub fn apply_row_scale(&mut self, scale: &HighsScale) {
        debug_assert!(self.format_ok());
        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    self.value_[i_el as usize] *= scale.row[self.index_[i_el as usize] as usize];
                }
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    self.value_[i_el as usize] *= scale.row[i_row as usize];
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::unapplyScale (lines 999-1018, 20 lines)
    pub fn unapply_scale(&mut self, scale: &HighsScale) {
        debug_assert!(self.format_ok());
        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    let i_row: HighsInt = self.index_[i_el as usize];
                    self.value_[i_el as usize] /=
                        scale.col[i_col as usize] * scale.row[i_row as usize];
                }
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    let i_col: HighsInt = self.index_[i_el as usize];
                    self.value_[i_el as usize] /=
                        scale.col[i_col as usize] * scale.row[i_row as usize];
                }
            }
        }
    }
}
