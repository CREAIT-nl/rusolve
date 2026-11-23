//! HiGHS source: util/HighsSparseMatrix.cpp (lines 1149-1745)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Matrix-vector products and PRICE (column-wise, row-wise, row-wise with switch to
//! column-wise; the switch density and hyper-sparse thresholds are numerically relevant).

use super::*;
use crate::prelude::*;

impl HighsSparseMatrix {
    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::alphaProductPlusY (lines 1149-1178, 30 lines)
    /// C++ default arguments: transpose = false
    pub fn alpha_product_plus_y(&self, alpha: f64, x: &[f64], y: &mut Vec<f64>, transpose: bool) {
        if self.is_colwise() {
            if transpose {
                for i_col in 0..self.num_col_ {
                    for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                        y[i_col as usize] += alpha
                            * self.value_[i_el as usize]
                            * x[self.index_[i_el as usize] as usize];
                    }
                }
            } else {
                for i_col in 0..self.num_col_ {
                    for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                        y[self.index_[i_el as usize] as usize] +=
                            alpha * self.value_[i_el as usize] * x[i_col as usize];
                    }
                }
            }
        } else if transpose {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    y[self.index_[i_el as usize] as usize] +=
                        alpha * self.value_[i_el as usize] * x[i_row as usize];
                }
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    y[i_row as usize] +=
                        alpha * self.value_[i_el as usize] * x[self.index_[i_el as usize] as usize];
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::product (lines 1180-1198, 19 lines)
    pub fn product(&self, result: &mut Vec<f64>, row: &[f64]) {
        result.clear();
        result.resize(self.num_row_ as usize, 0.0);
        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    result[self.index_[i_el as usize] as usize] +=
                        row[i_col as usize] * self.value_[i_el as usize];
                }
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    result[i_row as usize] +=
                        row[self.index_[i_el as usize] as usize] * self.value_[i_el as usize];
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::productTranspose (lines 1200-1218, 19 lines)
    pub fn product_transpose(&self, result: &mut Vec<f64>, col: &[f64]) {
        result.clear();
        result.resize(self.num_col_ as usize, 0.0);
        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    result[i_col as usize] +=
                        col[self.index_[i_el as usize] as usize] * self.value_[i_el as usize];
                }
            }
        } else {
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    result[self.index_[i_el as usize] as usize] +=
                        col[i_row as usize] * self.value_[i_el as usize];
                }
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::productQuad (lines 1220-1244, 25 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    pub fn product_quad(&self, result: &mut Vec<f64>, row: &[f64], debug_report: HighsInt) {
        result.clear();
        result.resize(self.num_row_ as usize, 0.0);
        if self.is_colwise() {
            let mut value: Vec<HighsCDouble> =
                vec![HighsCDouble::from(0.0); self.num_row_ as usize];
            for i_col in 0..self.num_col_ {
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    value[self.index_[i_el as usize] as usize] +=
                        row[i_col as usize] * self.value_[i_el as usize];
                }
            }
            for i_row in 0..self.num_row_ {
                result[i_row as usize] = f64::from(value[i_row as usize]);
            }
        } else {
            for i_row in 0..self.num_row_ {
                let mut value = HighsCDouble::from(0.0);
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    value += row[self.index_[i_el as usize] as usize] * self.value_[i_el as usize];
                }
                result[i_row as usize] = f64::from(value);
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::productTransposeQuad (lines 1246-1270, 25 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    /// Overload `(vector<double>& result, const vector<double>& row, debug_report)`.
    pub fn product_transpose_quad(
        &self,
        result: &mut Vec<f64>,
        row: &[f64],
        debug_report: HighsInt,
    ) {
        result.clear();
        result.resize(self.num_col_ as usize, 0.0);
        if self.is_colwise() {
            for i_col in 0..self.num_col_ {
                let mut value = HighsCDouble::from(0.0);
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    value += row[self.index_[i_el as usize] as usize] * self.value_[i_el as usize];
                }
                result[i_col as usize] = f64::from(value);
            }
        } else {
            let mut value: Vec<HighsCDouble> =
                vec![HighsCDouble::from(0.0); self.num_col_ as usize];
            for i_row in 0..self.num_row_ {
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    value[self.index_[i_el as usize] as usize] +=
                        row[i_row as usize] * self.value_[i_el as usize];
                }
            }
            for i_col in 0..self.num_col_ {
                result[i_col as usize] = f64::from(value[i_col as usize]);
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::productTransposeQuad (lines 1272-1316, 45 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    /// Overload `(result_value, result_index, const HVector& column, debug_report)`;
    /// uses `HighsSparseVectorSum` (crate::util::highs_sparse_vector_sum).
    pub fn product_transpose_quad_hvector(
        &self,
        result_value: &mut Vec<f64>,
        result_index: &mut Vec<HighsInt>,
        column: &HVector,
        debug_report: HighsInt,
    ) {
        if debug_report >= kDebugReportAll {
            println!("\nHighsSparseMatrix::productTranspose:");
        }
        if self.is_colwise() {
            result_value.reserve(self.num_col_ as usize);
            result_index.reserve(self.num_col_ as usize);
            for i_col in 0..self.num_col_ {
                let mut value = HighsCDouble::from(0.0);
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    value += column.array[self.index_[i_el as usize] as usize]
                        * self.value_[i_el as usize];
                }

                if value.abs() - kHighsTiny > 0.0 {
                    result_value.push(f64::from(value));
                    result_index.push(i_col);
                }
            }
        } else {
            let mut sum = HighsSparseVectorSum::new(self.num_col_);
            for i_row in 0..self.num_row_ {
                let multiplier = column.array[i_row as usize];
                for i_el in self.start_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                    sum.add(
                        self.index_[i_el as usize],
                        multiplier * self.value_[i_el as usize],
                    );
                }
            }
            if debug_report >= kDebugReportAll {
                let mut report_sum = HighsSparseVectorSum::new(self.num_col_);
                for i_row in 0..self.num_row_ {
                    let multiplier = column.array[i_row as usize];
                    if debug_report == kDebugReportAll || debug_report == i_row {
                        self.debug_report_row_price_sum(
                            i_row,
                            multiplier,
                            self.start_[(i_row + 1) as usize],
                            &mut report_sum,
                        );
                    }
                }
            }

            sum.cleanup(|_, x| x.abs() <= kHighsTiny);
            *result_index = std::mem::take(&mut sum.nonzeroinds);
            let result_num_nz = result_index.len() as HighsInt;
            result_value.reserve(result_num_nz as usize);
            for i in 0..result_num_nz {
                result_value.push(sum.get_value(result_index[i as usize]));
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::createRowwisePartitioned (lines 1318-1387, 70 lines)
    /// C++ default arguments: in_partition = NULL
    pub fn create_rowwise_partitioned(
        &mut self,
        matrix: &HighsSparseMatrix,
        in_partition: Option<&[i8]>,
    ) {
        let all_in_partition = in_partition.is_none();

        let num_col = matrix.num_col_;
        let num_row = matrix.num_row_;
        let num_nz = matrix.num_nz();
        let a_start = &matrix.start_;
        let a_index = &matrix.index_;
        let a_value = &matrix.value_;
        let ar_start = &mut self.start_;
        let ar_p_end = &mut self.p_end_;
        let ar_index = &mut self.index_;
        let ar_value = &mut self.value_;

        // Use ar_p_end and ar_end to compute lengths, which are then transformed into
        // the p_ends and ends of the inserted entries
        let mut ar_end: Vec<HighsInt> = Vec::new();
        ar_start.resize((num_row + 1) as usize, 0);
        ar_p_end.clear();
        ar_p_end.resize(num_row as usize, 0);
        ar_end.resize(num_row as usize, 0);
        // Count the nonzeros of nonbasic and basic columns in each row
        for i_col in 0..num_col {
            if all_in_partition || in_partition.unwrap()[i_col as usize] != 0 {
                for i_el in a_start[i_col as usize]..a_start[(i_col + 1) as usize] {
                    let i_row = a_index[i_el as usize];
                    ar_p_end[i_row as usize] += 1;
                }
            } else {
                for i_el in a_start[i_col as usize]..a_start[(i_col + 1) as usize] {
                    let i_row = a_index[i_el as usize];
                    ar_end[i_row as usize] += 1;
                }
            }
        }
        // Compute the starts and turn the lengths into ends
        ar_start[0] = 0;
        for i_row in 0..num_row as usize {
            ar_start[i_row + 1] = ar_start[i_row] + ar_p_end[i_row] + ar_end[i_row];
        }
        for i_row in 0..num_row as usize {
            ar_end[i_row] = ar_start[i_row] + ar_p_end[i_row];
            ar_p_end[i_row] = ar_start[i_row];
        }
        // Insert the entries
        ar_index.resize(num_nz as usize, 0);
        ar_value.resize(num_nz as usize, 0.0);
        for i_col in 0..num_col {
            if all_in_partition || in_partition.unwrap()[i_col as usize] != 0 {
                for i_el in a_start[i_col as usize]..a_start[(i_col + 1) as usize] {
                    let i_row = a_index[i_el as usize];
                    let i_to_el = ar_p_end[i_row as usize];
                    ar_p_end[i_row as usize] += 1;
                    ar_index[i_to_el as usize] = i_col;
                    ar_value[i_to_el as usize] = a_value[i_el as usize];
                }
            } else {
                for i_el in a_start[i_col as usize]..a_start[(i_col + 1) as usize] {
                    let i_row = a_index[i_el as usize];
                    let i_to_el = ar_end[i_row as usize];
                    ar_end[i_row as usize] += 1;
                    ar_index[i_to_el as usize] = i_col;
                    ar_value[i_to_el as usize] = a_value[i_el as usize];
                }
            }
        }
        self.format_ = MatrixFormat::kRowwisePartitioned;
        self.num_col_ = num_col;
        self.num_row_ = num_row;
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::debugPartitionOk (lines 1389-1410, 22 lines)
    pub fn debug_partition_ok(&self, in_partition: &[i8]) -> bool {
        let mut ok = true;
        for i_row in 0..self.num_row_ {
            for i_el in self.start_[i_row as usize]..self.p_end_[i_row as usize] {
                if in_partition[self.index_[i_el as usize] as usize] == 0 {
                    ok = false;
                    break;
                }
            }
            if !ok {
                break;
            }
            for i_el in self.p_end_[i_row as usize]..self.start_[(i_row + 1) as usize] {
                if in_partition[self.index_[i_el as usize] as usize] != 0 {
                    ok = false;
                    break;
                }
            }
            if !ok {
                break;
            }
        }
        ok
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::priceByColumn (lines 1412-1437, 26 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    pub fn price_by_column(
        &self,
        quad_precision: bool,
        result: &mut HVector,
        column: &HVector,
        debug_report: HighsInt,
    ) {
        if debug_report >= kDebugReportAll {
            println!("\nHighsSparseMatrix::priceByColumn:");
        }
        result.count = 0;
        for i_col in 0..self.num_col_ {
            let mut value = 0.0;
            if quad_precision {
                let mut quad_value = HighsCDouble::from(0.0);
                for i_el in self.start_[i_col as usize]..self.start_[(i_col + 1) as usize] {
                    quad_value += column.array[self.index_[i_el as usize] as usize]
                        * self.value_[i_el as usize];
                }
                value = f64::from(quad_value);
            } else {
                // Zipped subslices: one index range per column instead of three bounds
                // checks per nonzero. Same operand and summation order as C++ priceByColumn.
                let start = self.start_[i_col as usize] as usize;
                let end = self.start_[(i_col + 1) as usize] as usize;
                for (&i_row, &a) in self.index_[start..end].iter().zip(&self.value_[start..end]) {
                    value += column.array[i_row as usize] * a;
                }
            }
            if value.abs() > kHighsTiny {
                result.array[i_col as usize] = value;
                result.index[result.count as usize] = i_col;
                result.count += 1;
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::priceByRow (lines 1439-1455, 17 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    pub fn price_by_row(
        &self,
        quad_precision: bool,
        result: &mut HVector,
        column: &HVector,
        debug_report: HighsInt,
    ) {
        if debug_report >= kDebugReportAll {
            println!("\nHighsSparseMatrix::priceByRow:");
        }
        // Vanilla hyper-sparse row-wise PRICE. Set up parameters so that
        // priceByRowWithSwitch runs as vanilla hyper-sparse PRICE
        // Expected density always forces hyper-sparse PRICE
        let expected_density = -kHighsInf;
        // Always start from first index of column
        let from_index = 0;
        // Never switch to standard row-wise PRICE
        let switch_density = kHighsInf;
        self.price_by_row_with_switch(
            quad_precision,
            result,
            column,
            expected_density,
            from_index,
            switch_density,
            kDebugReportOff,
        );
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::priceByRowWithSwitch (lines 1457-1570, 114 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    pub fn price_by_row_with_switch(
        &self,
        quad_precision: bool,
        result: &mut HVector,
        column: &HVector,
        expected_density: f64,
        from_index: HighsInt,
        switch_density: f64,
        debug_report: HighsInt,
    ) {
        let mut sum = HighsSparseVectorSum::default();
        if quad_precision {
            sum.set_dimension(self.num_col_);
        }
        if debug_report >= kDebugReportAll {
            println!("\nHighsSparseMatrix::priceByRowWithSwitch");
        }
        // (Continue) hyper-sparse row-wise PRICE with possible switches to
        // standard row-wise PRICE either immediately based on historical
        // density or during hyper-sparse PRICE if there is too much fill-in
        let mut next_index = from_index;
        // Possibly don't perform hyper-sparse PRICE based on historical density
        if expected_density <= kHyperPriceDensity {
            let inv_num_col = 1.0 / f64::from(self.num_col_);
            let mut ix = next_index;
            while ix < column.count {
                let i_row = column.index[ix as usize];
                // Determine whether p_end_ or the next start_ ends the loop
                let to_i_el = if self.format_ == MatrixFormat::kRowwisePartitioned {
                    self.p_end_[i_row as usize]
                } else {
                    self.start_[(i_row + 1) as usize]
                };
                // Possibly switch to standard row-wise price
                let row_num_nz = to_i_el - self.start_[i_row as usize];
                let local_density = (1.0 * f64::from(result.count)) * inv_num_col;
                let switch_to_dense =
                    result.count + row_num_nz >= self.num_col_ || local_density > switch_density;
                if switch_to_dense {
                    break;
                }
                let multiplier = column.array[i_row as usize];
                if debug_report == kDebugReportAll || debug_report == i_row {
                    self.debug_report_row_price(i_row, multiplier, to_i_el, &result.array);
                }
                if multiplier != 0.0 {
                    if quad_precision {
                        for i_el in self.start_[i_row as usize]..to_i_el {
                            sum.add(
                                self.index_[i_el as usize],
                                multiplier * self.value_[i_el as usize],
                            );
                        }
                    } else {
                        // PORT NOTE (perf): the per-element loop over zipped subslices;
                        // same operations in the same order.
                        let (s_el, e_el) = (self.start_[i_row as usize] as usize, to_i_el as usize);
                        if s_el < e_el {
                            let array = &mut result.array[..];
                            for (&i_col, &mv) in
                                self.index_[s_el..e_el].iter().zip(&self.value_[s_el..e_el])
                            {
                                debug_assert!((i_col as usize) < array.len());
                                // SAFETY: i_col is a column index of this matrix
                                // (0 <= i_col < num_col_) and result.array has length
                                // num_col_ (HVector::setup(num_col)).
                                let slot = unsafe { array.get_unchecked_mut(i_col as usize) };
                                let value0 = *slot;
                                let value1 = value0 + multiplier * mv;
                                if value0 == 0.0 {
                                    result.index[result.count as usize] = i_col;
                                    result.count += 1;
                                }
                                *slot = if value1.abs() < kHighsTiny {
                                    kHighsZero
                                } else {
                                    value1
                                };
                            }
                        }
                    }
                }
                next_index = ix + 1;
                ix += 1;
            }
        }
        if quad_precision {
            sum.cleanup(|_, x| x.abs() <= kHighsTiny);
        }
        if next_index < column.count {
            // PRICE is not complete: finish without maintaining nonzeros of result
            if quad_precision {
                let mut result_array: Vec<HighsCDouble> = sum.values.clone();
                self.price_by_row_dense_result_cdouble(
                    &mut result_array,
                    column,
                    next_index,
                    kDebugReportOff,
                );
                // Determine indices of nonzeros in result
                result.count = 0;
                for i_col in 0..self.num_col_ {
                    let value1 = f64::from(result_array[i_col as usize]);
                    if value1.abs() < kHighsTiny {
                        result.array[i_col as usize] = 0.0;
                    } else {
                        result.array[i_col as usize] = value1;
                        result.index[result.count as usize] = i_col;
                        result.count += 1;
                    }
                }
            } else {
                self.price_by_row_dense_result(
                    &mut result.array,
                    column,
                    next_index,
                    kDebugReportOff,
                );
                // Determine indices of nonzeros in result
                result.count = 0;
                for i_col in 0..self.num_col_ {
                    let value1 = result.array[i_col as usize];
                    if value1.abs() < kHighsTiny {
                        result.array[i_col as usize] = 0.0;
                    } else {
                        result.index[result.count as usize] = i_col;
                        result.count += 1;
                    }
                }
            }
        } else if quad_precision {
            // HVector result should have result.index of size this->num_col_
            // by virtue of result.setup. However, it will generally lose
            // this property by virtue of the following move
            result.index = std::mem::take(&mut sum.nonzeroinds);
            let result_num_nz = result.index.len() as HighsInt;
            // Restore the size of result.index
            result.index.resize(self.num_col_ as usize, 0);
            result.count = result_num_nz;
            for i in 0..result_num_nz {
                let i_row = result.index[i as usize];
                result.array[i_row as usize] = sum.get_value(i_row);
            }
        } else {
            // PRICE is complete maintaining nonzeros of result
            // Remove small values
            result.tight();
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::update (lines 1572-1602, 31 lines)
    pub fn update(&mut self, var_in: HighsInt, var_out: HighsInt, matrix: &HighsSparseMatrix) {
        if var_in < self.num_col_ {
            for i_el in matrix.start_[var_in as usize]..matrix.start_[(var_in + 1) as usize] {
                let i_row = matrix.index_[i_el as usize];
                let mut i_find = self.start_[i_row as usize];
                self.p_end_[i_row as usize] -= 1;
                let i_swap = self.p_end_[i_row as usize];
                while self.index_[i_find as usize] != var_in {
                    i_find += 1;
                }
                self.index_.swap(i_find as usize, i_swap as usize);
                self.value_.swap(i_find as usize, i_swap as usize);
            }
        }
        if var_out < self.num_col_ {
            for i_el in matrix.start_[var_out as usize]..matrix.start_[(var_out + 1) as usize] {
                let i_row = matrix.index_[i_el as usize];
                let mut i_find = self.p_end_[i_row as usize];
                let i_swap = self.p_end_[i_row as usize];
                self.p_end_[i_row as usize] += 1;
                while self.index_[i_find as usize] != var_out {
                    i_find += 1;
                }
                self.index_.swap(i_find as usize, i_swap as usize);
                self.value_.swap(i_find as usize, i_swap as usize);
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::computeDot (lines 1604-1616, 13 lines)
    pub fn compute_dot(&self, array: &[f64], use_col: HighsInt) -> f64 {
        let mut result = 0.0;
        if use_col < self.num_col_ {
            for i_el in self.start_[use_col as usize]..self.start_[(use_col + 1) as usize] {
                result += array[self.index_[i_el as usize] as usize] * self.value_[i_el as usize];
            }
        } else {
            result = array[(use_col - self.num_col_) as usize];
        }
        result
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::collectAj (lines 1618-1637, 20 lines)
    pub fn collect_aj(&self, column: &mut HVector, use_col: HighsInt, multiplier: f64) {
        if use_col < self.num_col_ {
            for i_el in self.start_[use_col as usize]..self.start_[(use_col + 1) as usize] {
                let i_row = self.index_[i_el as usize];
                let value0 = column.array[i_row as usize];
                let value1 = value0 + multiplier * self.value_[i_el as usize];
                if value0 == 0.0 {
                    column.index[column.count as usize] = i_row;
                    column.count += 1;
                }
                column.array[i_row as usize] = if value1.abs() < kHighsTiny {
                    kHighsZero
                } else {
                    value1
                };
            }
        } else {
            let i_row = use_col - self.num_col_;
            let value0 = column.array[i_row as usize];
            let value1 = value0 + multiplier;
            if value0 == 0.0 {
                column.index[column.count as usize] = i_row;
                column.count += 1;
            }
            column.array[i_row as usize] = if value1.abs() < kHighsTiny {
                kHighsZero
            } else {
                value1
            };
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::priceByRowDenseResult (lines 1639-1665, 27 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    /// Overload with `std::vector<double>& result`.
    pub fn price_by_row_dense_result(
        &self,
        result: &mut Vec<f64>,
        column: &HVector,
        from_index: HighsInt,
        debug_report: HighsInt,
    ) {
        // Assumes that result is zeroed beforehand - in case continuing
        // priceByRow after switch from sparse
        for ix in from_index..column.count {
            let i_row = column.index[ix as usize];
            let multiplier = column.array[i_row as usize];
            // Determine whether p_end_ or the next start_ should be used to end the
            // loop
            let to_i_el = if self.format_ == MatrixFormat::kRowwisePartitioned {
                self.p_end_[i_row as usize]
            } else {
                self.start_[(i_row + 1) as usize]
            };
            if debug_report == kDebugReportAll || debug_report == i_row {
                self.debug_report_row_price(i_row, multiplier, to_i_el, result);
            }
            for i_el in self.start_[i_row as usize]..to_i_el {
                let i_col = self.index_[i_el as usize];
                let value0 = result[i_col as usize];
                let value1 = value0 + multiplier * self.value_[i_el as usize];
                result[i_col as usize] = if value1.abs() < kHighsTiny {
                    kHighsZero
                } else {
                    value1
                };
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::priceByRowDenseResult (lines 1667-1691, 25 lines)
    /// C++ default arguments: debug_report = kDebugReportOff
    /// Overload with `std::vector<HighsCDouble>& result` (quad precision).
    pub fn price_by_row_dense_result_cdouble(
        &self,
        result: &mut Vec<HighsCDouble>,
        column: &HVector,
        from_index: HighsInt,
        debug_report: HighsInt,
    ) {
        // Assumes that result is zeroed beforehand - in case continuing
        // priceByRow after switch from sparse
        for ix in from_index..column.count {
            let i_row = column.index[ix as usize];
            let multiplier = column.array[i_row as usize];
            // Determine whether p_end_ or the next start_ should be used to end the
            // loop
            let to_i_el = if self.format_ == MatrixFormat::kRowwisePartitioned {
                self.p_end_[i_row as usize]
            } else {
                self.start_[(i_row + 1) as usize]
            };
            for i_el in self.start_[i_row as usize]..to_i_el {
                let i_col = self.index_[i_el as usize];
                let value0 = result[i_col as usize];
                let value1 = value0 + multiplier * self.value_[i_el as usize];
                result[i_col as usize] = if f64::from(value1).abs() < kHighsTiny {
                    HighsCDouble::from(kHighsZero)
                } else {
                    value1
                };
            }
        }
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::debugReportRowPrice (lines 1693-1709, 17 lines)
    /// Debug report: out of scope, no-op.
    pub fn debug_report_row_price(
        &self,
        i_row: HighsInt,
        multiplier: f64,
        to_i_el: HighsInt,
        result: &[f64],
    ) {
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::debugReportRowPrice (lines 1711-1727, 17 lines)
    /// Debug report: out of scope, no-op.
    pub fn debug_report_row_price_cdouble(
        &self,
        i_row: HighsInt,
        multiplier: f64,
        to_i_el: HighsInt,
        result: &[HighsCDouble],
    ) {
    }

    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::debugReportRowPrice (lines 1729-1745, 17 lines)
    /// Debug report: out of scope, no-op.
    pub fn debug_report_row_price_sum(
        &self,
        i_row: HighsInt,
        multiplier: f64,
        to_i_el: HighsInt,
        sum: &mut HighsSparseVectorSum,
    ) {
    }

    /// HiGHS: HighsSparseMatrix.h::HighsSparseMatrix::computeDot(const HVector&, HighsInt) (inline)
    pub fn compute_dot_hvector(&self, column: &HVector, use_col: HighsInt) -> f64 {
        self.compute_dot(&column.array, use_col)
    }
}
