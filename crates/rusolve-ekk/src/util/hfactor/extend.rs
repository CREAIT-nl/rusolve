//! HiGHS source: util/HFactorExtend.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Extending an existing factorization by new columns/rows.

use super::*;
use crate::prelude::*;

impl HFactor {
    /// HiGHS: HFactorExtend.cpp::HFactor::addCols (lines 18-21, 4 lines)
    pub fn add_cols(&mut self, num_new_col: HighsInt) {
        self.invalid_a_matrix_action();
        self.num_col += num_new_col;
    }

    /// HiGHS: HFactorExtend.cpp::HFactor::deleteNonbasicCols (lines 23-26, 4 lines)
    pub fn delete_nonbasic_cols(&mut self, num_deleted_col: HighsInt) {
        self.invalid_a_matrix_action();
        self.num_col -= num_deleted_col;
    }

    /// HiGHS: HFactorExtend.cpp::HFactor::addRows (lines 28-229, 202 lines)
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    pub fn add_rows(&mut self, basic_index: &mut [HighsInt], ar_matrix: &HighsSparseMatrix) {
        self.invalid_a_matrix_action();
        let extend_invert_when_adding_rows = kExtendInvertWhenAddingRows;
        debug_assert!(extend_invert_when_adding_rows);
        let num_new_row: HighsInt = ar_matrix.num_row_;
        let new_num_row: HighsInt = self.num_row + num_new_row;
        println!(
            "Adding {} new rows to HFactor instance: increasing dimension from {} to {} ",
            num_new_row, self.num_row, new_num_row
        );

        // Need to know where (if) a column is basic
        let mut in_basis: Vec<HighsInt> = vec![-1; self.num_col as usize];
        for i_row in 0..self.num_row {
            let i_var = basic_index[i_row as usize];
            if i_var < self.num_col {
                in_basis[i_var as usize] = i_row;
            }
        }
        for i_row in self.num_row..new_num_row {
            let i_var = basic_index[i_row as usize];
            debug_assert!(i_var >= self.num_col);
        }

        // Create a row-wise sparse matrix containing the new rows of the L
        // matrix
        let mut new_lr_rows = HighsSparseMatrix::default();
        new_lr_rows.format_ = MatrixFormat::kRowwise;
        new_lr_rows.num_col_ = self.num_row;
        let mut expected_density: f64 = 0.0;
        let mut rhs = HVector::default();
        rhs.setup(self.num_row);
        self.lr_start.reserve((new_num_row + 1) as usize);
        for inew_row in 0..num_new_row {
            // Prepare RHS for system U^T.v = r
            rhs.clear();
            rhs.pack_flag = true;
            let mut i_el = ar_matrix.start_[inew_row as usize];
            while i_el < ar_matrix.start_[(inew_row + 1) as usize] {
                let i_col = ar_matrix.index_[i_el as usize];
                let basis_index = in_basis[i_col as usize];
                if basis_index >= 0 {
                    rhs.array[basis_index as usize] = ar_matrix.value_[i_el as usize];
                    rhs.index[rhs.count as usize] = basis_index;
                    rhs.count += 1;
                }
                i_el += 1;
            }
            // Solve U^T.v = r
            self.btran_u(&mut rhs, expected_density, None);
            let local_density: f64 = (1.0 * rhs.count as f64) / self.num_row as f64;
            expected_density = kRunningAverageMultiplier * local_density
                + (1.0 - kRunningAverageMultiplier) * expected_density;
            rhs.tight();
            //
            // Append v to the matrix containing the new rows of L
            let rhs_num_nz: HighsInt = rhs.count;
            for i_x in 0..rhs_num_nz {
                let i_col = rhs.index[i_x as usize];
                new_lr_rows.index_.push(i_col);
                new_lr_rows.value_.push(rhs.array[i_col as usize]);
            }
            new_lr_rows
                .start_
                .push(new_lr_rows.index_.len() as HighsInt);
            new_lr_rows.num_row_ += 1;

            //
            // Append v to the L matrix
            for i_x in 0..rhs_num_nz {
                let i_col = rhs.index[i_x as usize];
                self.lr_index.push(i_col);
                self.lr_value.push(rhs.array[i_col as usize]);
            }
            self.lr_start.push(self.lr_index.len() as HighsInt);
        }
        // Now create a column-wise copy of the new rows
        let mut new_lr_cols = new_lr_rows.clone();
        new_lr_cols.ensure_colwise();
        //
        // Insert the column-wise copy into the L matrix
        //
        // Add pivot indices for the new columns
        self.l_pivot_index.resize(new_num_row as usize, 0);
        for i_col in self.num_row..new_num_row {
            self.l_pivot_index[i_col as usize] = i_col;
        }
        //
        // Add starts for the identity columns
        let l_matrix_new_num_nz: HighsInt = self.lr_index.len() as HighsInt;
        debug_assert!(
            l_matrix_new_num_nz == (self.l_index.len() + new_lr_cols.index_.len()) as HighsInt
        );
        self.l_start.resize((new_num_row + 1) as usize, 0);
        let mut to_el: HighsInt = l_matrix_new_num_nz;
        for i_col in (self.num_row + 1)..(new_num_row + 1) {
            self.l_start[i_col as usize] = l_matrix_new_num_nz;
        }
        //
        // Insert the new entries, remembering to offset the index values by
        // num_row, since new_lr_cols only has the new rows
        self.l_index.resize(l_matrix_new_num_nz as usize, 0);
        self.l_value.resize(l_matrix_new_num_nz as usize, 0.0);
        let mut i_col: HighsInt = self.num_row - 1;
        while i_col >= 0 {
            let from_el: HighsInt = self.l_start[(i_col + 1) as usize];
            self.l_start[(i_col + 1) as usize] = to_el;
            let mut i_el: HighsInt = new_lr_cols.start_[(i_col + 1) as usize] - 1;
            while i_el >= new_lr_cols.start_[i_col as usize] {
                to_el -= 1;
                self.l_index[to_el as usize] = self.num_row + new_lr_cols.index_[i_el as usize];
                self.l_value[to_el as usize] = new_lr_cols.value_[i_el as usize];
                i_el -= 1;
            }
            let mut i_el: HighsInt = from_el - 1;
            while i_el >= self.l_start[i_col as usize] {
                to_el -= 1;
                self.l_index[to_el as usize] = self.l_index[i_el as usize];
                self.l_value[to_el as usize] = self.l_value[i_el as usize];
                i_el -= 1;
            }
            i_col -= 1;
        }
        debug_assert!(to_el == 0);
        self.l_pivot_lookup.resize(new_num_row as usize, 0);
        for i_row in self.num_row..new_num_row {
            let idx = self.l_pivot_index[i_row as usize];
            self.l_pivot_lookup[idx as usize] = i_row;
        }
        // Now update the U matrix with identity rows and columns
        //
        // Now add pivots corresponding to identity columns in U.
        let u_count_x: HighsInt = self.u_index.len() as HighsInt;
        let u_pivot_lookup_offset: HighsInt = self.u_pivot_index.len() as HighsInt - self.num_row;
        for i_row in self.num_row..new_num_row {
            self.u_pivot_lookup.push(u_pivot_lookup_offset + i_row);
            self.u_pivot_index.push(i_row);
            self.u_pivot_value.push(1.0);
            self.u_start.push(u_count_x);
            self.u_last_p.push(u_count_x);
        }

        // Now, to extend UR
        //
        // Borrowing names from buildFinish()
        let ur_stuff_size: HighsInt = if self.update_method == kUpdateMethodFt {
            5
        } else {
            0
        };
        let ur_size: HighsInt = self.ur_index.len() as HighsInt;
        let ur_count_size: HighsInt = ur_size + ur_stuff_size * num_new_row;
        self.ur_index.resize(ur_count_size as usize, 0);
        self.ur_value.resize(ur_count_size as usize, 0.0);

        // Need to refer to just the new UR vectors
        let ur_cur_num_vec: HighsInt = self.ur_start.len() as HighsInt;
        let ur_new_num_vec: HighsInt = ur_cur_num_vec + num_new_row;
        println!(
            "\nUpdating UR vectors {} - {}",
            ur_cur_num_vec,
            ur_new_num_vec - 1
        );
        // UR pointer
        //
        // Allow space to the start of new rows, including the start for the
        // fictitious ur_new_num_vec'th row
        self.ur_start.resize((ur_new_num_vec + 1) as usize, 0);
        for i_row in (ur_cur_num_vec + 1)..(ur_new_num_vec + 1) {
            self.ur_start[i_row as usize] = ur_size;
        }
        // NB ur_temp plays the role of ur_lastp when it could be used as
        // temporary storage in buildFinish()
        let mut ur_temp: Vec<HighsInt> = vec![0; ur_new_num_vec as usize];
        //
        // ur_space has its new entries assigned (to ur_stuff_size) as in
        // buildFinish()
        self.ur_space.resize(ur_new_num_vec as usize, 0);
        for i_row in ur_cur_num_vec..ur_new_num_vec {
            self.ur_space[i_row as usize] = ur_stuff_size;
        }
        // Compute ur_temp exactly as in buildFinish()
        for k in 0..u_count_x {
            ur_temp[self.u_pivot_lookup[self.u_index[k as usize] as usize] as usize] += 1;
        }
        let mut i_start: HighsInt = ur_size;
        self.ur_start[ur_cur_num_vec as usize] = i_start;
        for i_row in (ur_cur_num_vec + 1)..(ur_new_num_vec + 1) {
            let gap: HighsInt = ur_temp[(i_row - 1) as usize] + ur_stuff_size;
            self.ur_start[i_row as usize] = i_start + gap;
            i_start += gap;
            println!(
                "ur_start[{}] = {}; gap = {}; iStart = {}",
                i_row, self.ur_start[i_row as usize], gap, i_start
            );
        }
        println!("ur_count_size = {}; iStart{}", ur_count_size, i_start);
        // Lose the start for the fictitious ur_new_num_vec'th row
        self.ur_start.resize(ur_new_num_vec as usize, 0);
        // Resize ur_lastp and initialise its new entries to be the ur_start
        // values since the rows are empty
        self.ur_lastp.resize(ur_new_num_vec as usize, 0);
        for i_row in ur_cur_num_vec..ur_new_num_vec {
            self.ur_lastp[i_row as usize] = self.ur_start[i_row as usize];
        }
        //
        // Increase the number of rows in HFactor
        self.num_row += num_new_row;
    }
}
