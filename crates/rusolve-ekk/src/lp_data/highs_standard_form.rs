//! HiGHS source: lp_data/Highs.cpp (getStandardFormLp),
//! lp_data/HighsInterface.cpp (formStandardFormLp)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The standard form `min c'x s.t. Ax = b, x >= 0` of the
//! incumbent LP, cached in `Highs::standard_form_*` (highs_run.rs) and
//! invalidated by `clear_standard_form_lp` (highs_api_support.rs, already
//! ported, called from `clear_derived_model_properties`).
//! `logHeader` is logging and is dropped.

use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::getStandardFormLp (lines 2063-2094).
    /// Each nullable output pointer (`double* cost = nullptr`, ...) ->
    /// `Option<&mut [T]>`, written only when `Some`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_standard_form_lp(
        &mut self,
        num_col: &mut HighsInt,
        num_row: &mut HighsInt,
        num_nz: &mut HighsInt,
        offset: &mut f64,
        mut cost: Option<&mut [f64]>,
        mut rhs: Option<&mut [f64]>,
        mut start: Option<&mut [HighsInt]>,
        mut index: Option<&mut [HighsInt]>,
        mut value: Option<&mut [f64]>,
    ) -> HighsStatus {
        // logHeader dropped (line 2067)
        if !self.standard_form_valid_ {
            let status = self.form_standard_form_lp();
            debug_assert!(status == HighsStatus::kOk);
            if status != HighsStatus::kOk {
                return status;
            }
        }
        *num_col = self.standard_form_cost_.len() as HighsInt;
        *num_row = self.standard_form_rhs_.len() as HighsInt;
        *num_nz = self.standard_form_matrix_.start_[*num_col as usize];
        *offset = self.standard_form_offset_;
        for i_col in 0..*num_col {
            if let Some(cost_out) = cost.as_deref_mut() {
                cost_out[i_col as usize] = self.standard_form_cost_[i_col as usize];
            }
            if let Some(start_out) = start.as_deref_mut() {
                start_out[i_col as usize] = self.standard_form_matrix_.start_[i_col as usize];
            }
            if index.is_some() || value.is_some() {
                for i_el in self.standard_form_matrix_.start_[i_col as usize]
                    ..self.standard_form_matrix_.start_[(i_col + 1) as usize]
                {
                    if let Some(index_out) = index.as_deref_mut() {
                        index_out[i_el as usize] = self.standard_form_matrix_.index_[i_el as usize];
                    }
                    if let Some(value_out) = value.as_deref_mut() {
                        value_out[i_el as usize] = self.standard_form_matrix_.value_[i_el as usize];
                    }
                }
            }
        }
        if let Some(start_out) = start {
            start_out[*num_col as usize] = self.standard_form_matrix_.start_[*num_col as usize];
        }
        if let Some(rhs_out) = rhs {
            for i_row in 0..*num_row {
                rhs_out[i_row as usize] = self.standard_form_rhs_[i_row as usize];
            }
        }
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::formStandardFormLp (lines 128-350).
    /// The incumbent matrix is made row-wise for `getRow` and back column-wise
    /// (`matrix.ensureColwise()`) before the column pass, as in C++. The
    /// closing "Standard form LP obtained ..." `highsLogUser` is kept.
    pub fn form_standard_form_lp(&mut self) -> HighsStatus {
        self.clear_standard_form_lp();
        // HighsLp& lp = this->model_.lp_;
        // HighsSparseMatrix& matrix = lp.a_matrix_;
        // Ensure that the incumbent LP and standard form LP matrices are rowwise
        self.model_lp_.a_matrix_.ensure_rowwise();
        // Original rows are processed before columns, so that any original
        // boxed rows can be transformed to pairs of one-sided rows,
        // requiring the standard form matrix to be row-wise. The original
        // columns are assumed to come before any new columns, so their
        // costs (as a minimization) must be defined befor costs of new
        // columns.
        // Determine the objective scaling, and apply it to any offset
        let sense: HighsInt = self.model_lp_.sense_ as HighsInt;
        self.standard_form_offset_ = (sense as f64) * self.model_lp_.offset_;
        for i_col in 0..self.model_lp_.num_col_ {
            self.standard_form_cost_
                .push((sense as f64) * self.model_lp_.col_cost_[i_col as usize]);
        }
        self.standard_form_matrix_.format_ = MatrixFormat::kRowwise;
        self.standard_form_matrix_.num_col_ = self.model_lp_.num_col_;
        // Create a HighsSparseMatrix instance to store rows extracted from
        // the original constraint matrix
        let local_row_min_nnz: HighsInt = cpp_max(self.model_lp_.num_col_, 2);
        let mut local_row = HighsSparseMatrix::default();
        local_row.ensure_rowwise();
        local_row.num_row_ = 1;
        local_row.num_col_ = self.model_lp_.num_col_;
        local_row.index_.resize(local_row_min_nnz as usize, 0);
        local_row.value_.resize(local_row_min_nnz as usize, 0.0);
        local_row.start_.resize(2, 0);
        // HighsInt& num_nz = local_row.start_[1]; (alias: getRow below
        // writes `&mut local_row.start_[1]` directly)
        local_row.start_[0] = 0;
        let mut num_fixed_row: HighsInt = 0;
        let mut num_boxed_row: HighsInt = 0;
        let mut num_lower_row: HighsInt = 0;
        let mut num_upper_row: HighsInt = 0;
        let mut num_free_row: HighsInt = 0;
        let mut num_fixed_col: HighsInt = 0;
        let mut num_boxed_col: HighsInt = 0;
        let mut num_lower_col: HighsInt = 0;
        let mut num_upper_col: HighsInt = 0;
        let mut num_free_col: HighsInt = 0;
        let mut slack_ix: Vec<HighsInt> = Vec::new();
        for i_row in 0..self.model_lp_.num_row_ {
            let lower = self.model_lp_.row_lower_[i_row as usize];
            let upper = self.model_lp_.row_upper_[i_row as usize];
            if lower <= -kHighsInf && upper >= kHighsInf {
                debug_assert!(false); // assert(0 == 1)
                // Free row
                num_free_row += 1;
                continue;
            }
            if lower == upper {
                // Equality row
                num_fixed_row += 1;
                self.model_lp_.a_matrix_.get_row(
                    i_row,
                    &mut local_row.start_[1],
                    Some(local_row.index_.as_mut_slice()),
                    Some(local_row.value_.as_mut_slice()),
                );
                self.standard_form_matrix_.add_rows(&local_row, None);
                self.standard_form_rhs_.push(upper);
                continue;
            } else if lower <= -kHighsInf {
                // Upper bounded row, so record the slack
                num_upper_row += 1;
                debug_assert!(upper < kHighsInf);
                let standard_form_row: HighsInt = self.standard_form_rhs_.len() as HighsInt;
                slack_ix.push(standard_form_row + 1);
                self.model_lp_.a_matrix_.get_row(
                    i_row,
                    &mut local_row.start_[1],
                    Some(local_row.index_.as_mut_slice()),
                    Some(local_row.value_.as_mut_slice()),
                );
                self.standard_form_matrix_.add_rows(&local_row, None);
                self.standard_form_rhs_.push(upper);
            } else if upper >= kHighsInf {
                // Lower bounded row, so record the slack
                num_lower_row += 1;
                debug_assert!(lower > -kHighsInf);
                let standard_form_row: HighsInt = self.standard_form_rhs_.len() as HighsInt;
                slack_ix.push(-(standard_form_row + 1));
                self.model_lp_.a_matrix_.get_row(
                    i_row,
                    &mut local_row.start_[1],
                    Some(local_row.index_.as_mut_slice()),
                    Some(local_row.value_.as_mut_slice()),
                );
                self.standard_form_matrix_.add_rows(&local_row, None);
                self.standard_form_rhs_.push(lower);
            } else {
                // Boxed row, so record the lower slack
                debug_assert!(lower > -kHighsInf);
                debug_assert!(upper < kHighsInf);
                num_boxed_row += 1;
                let standard_form_row: HighsInt = self.standard_form_rhs_.len() as HighsInt;
                slack_ix.push(-(standard_form_row + 1));
                self.model_lp_.a_matrix_.get_row(
                    i_row,
                    &mut local_row.start_[1],
                    Some(local_row.index_.as_mut_slice()),
                    Some(local_row.value_.as_mut_slice()),
                );
                self.standard_form_matrix_.add_rows(&local_row, None);
                self.standard_form_rhs_.push(lower);
                // .. and upper slack, adding a copy of the row
                let standard_form_row: HighsInt = self.standard_form_rhs_.len() as HighsInt;
                slack_ix.push(standard_form_row + 1);
                self.standard_form_matrix_.add_rows(&local_row, None);
                self.standard_form_rhs_.push(upper);
            }
        }
        // Add rows corresponding to boxed columns
        for i_col in 0..self.model_lp_.num_col_ {
            let lower = self.model_lp_.col_lower_[i_col as usize];
            let upper = self.model_lp_.col_upper_[i_col as usize];
            if lower > -kHighsInf && upper < kHighsInf {
                // Boxed column
                //
                // x will be replaced by x = l + X (below) with X >= 0
                //
                // Introduce variable s >= 0 so that (with x >= l still)
                //
                // x = u - s => x + s = u
                self.standard_form_cost_.push(0.0);
                self.standard_form_matrix_.num_col_ += 1;
                local_row.num_col_ += 1;
                local_row.index_[0] = i_col;
                local_row.index_[1] = self.standard_form_matrix_.num_col_ - 1;
                local_row.value_[0] = 1.0;
                local_row.value_[1] = 1.0;
                local_row.start_[1] = 2;
                self.standard_form_matrix_.add_rows(&local_row, None);
                self.standard_form_rhs_.push(upper);
            }
        }
        // Finished with both matrices, row-wise, so ensure that the
        // incumbent matrix leaves col-wise, and that the standard form
        // matrix is col-wise so RHS shifts can be applied and more columns
        // can be added
        self.model_lp_.a_matrix_.ensure_colwise();
        self.standard_form_matrix_.ensure_colwise();
        // Work through the columns, ensuring that all have non-negativity bounds
        for i_col in 0..self.model_lp_.num_col_ {
            let cost = (sense as f64) * self.model_lp_.col_cost_[i_col as usize];
            let lower = self.model_lp_.col_lower_[i_col as usize];
            let upper = self.model_lp_.col_upper_[i_col as usize];
            if lower > -kHighsInf {
                // Finite lower bound
                if upper < kHighsInf {
                    if lower == upper {
                        // Fixed column
                        num_fixed_col += 1;
                    } else {
                        // Boxed column
                        num_boxed_col += 1;
                    }
                } else {
                    // Lower column
                    num_lower_col += 1;
                }
                if lower != 0.0 {
                    // x >= l, so shift x-l = X >= 0, giving x = X + l
                    //
                    // Cost contribution c(X+l) = cX + cl
                    self.standard_form_offset_ += cost * lower;
                    // Constraint contribution a(X+l) = aX + al
                    for i_el in self.standard_form_matrix_.start_[i_col as usize]
                        ..self.standard_form_matrix_.start_[(i_col + 1) as usize]
                    {
                        self.standard_form_rhs_
                            [self.standard_form_matrix_.index_[i_el as usize] as usize] -=
                            self.standard_form_matrix_.value_[i_el as usize] * lower;
                    }
                }
            } else if upper < kHighsInf {
                // Upper column
                num_upper_col += 1;
                // Have to operate even if u=0, since cost and column values are negated
                //
                // x <= u, so shift u-x = X >= 0, giving x = u - X
                //
                // Cost contribution c(u-X) = cu - cX
                self.standard_form_offset_ += cost * upper;
                self.standard_form_cost_[i_col as usize] = -cost;
                // Constraint contribution a(u-X) = -aX + au
                for i_el in self.standard_form_matrix_.start_[i_col as usize]
                    ..self.standard_form_matrix_.start_[(i_col + 1) as usize]
                {
                    self.standard_form_rhs_
                        [self.standard_form_matrix_.index_[i_el as usize] as usize] -=
                        self.standard_form_matrix_.value_[i_el as usize] * upper;
                    self.standard_form_matrix_.value_[i_el as usize] =
                        -self.standard_form_matrix_.value_[i_el as usize];
                }
            } else {
                // Free column
                num_free_col += 1;
                // Represent as x = x+ - x-
                //
                // where original column is now x+ >= 0
                //
                // and x- >= 0 has negation of its cost and matrix column
                self.standard_form_cost_.push(-cost);
                for i_el in self.standard_form_matrix_.start_[i_col as usize]
                    ..self.standard_form_matrix_.start_[(i_col + 1) as usize]
                {
                    self.standard_form_matrix_
                        .index_
                        .push(self.standard_form_matrix_.index_[i_el as usize]);
                    self.standard_form_matrix_
                        .value_
                        .push(-self.standard_form_matrix_.value_[i_el as usize]);
                }
                self.standard_form_matrix_
                    .start_
                    .push(self.standard_form_matrix_.index_.len() as HighsInt);
            }
        }
        // Now add the slack variables
        for &i_row in &slack_ix {
            self.standard_form_cost_.push(0.0);
            if i_row > 0 {
                self.standard_form_matrix_.index_.push(i_row - 1);
                self.standard_form_matrix_.value_.push(1.0);
            } else {
                self.standard_form_matrix_.index_.push(-i_row - 1);
                self.standard_form_matrix_.value_.push(-1.0);
            }
            self.standard_form_matrix_
                .start_
                .push(self.standard_form_matrix_.index_.len() as HighsInt);
        }
        // Now set correct values for the dimensions of this->standard_form_matrix_
        self.standard_form_matrix_.num_col_ = self.standard_form_cost_.len() as HighsInt;
        self.standard_form_matrix_.num_row_ = self.standard_form_rhs_.len() as HighsInt;
        self.standard_form_valid_ = true;
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "Standard form LP obtained for LP with (free / lower / upper / boxed / fixed) variables ({} / {} / {} / {} / {}) and constraints ({} / {} / {} / {} / {}) \n",
            num_free_col,
            num_lower_col,
            num_upper_col,
            num_boxed_col,
            num_fixed_col,
            num_free_row,
            num_lower_row,
            num_upper_row,
            num_boxed_row,
            num_fixed_row
        );
        HighsStatus::kOk
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// max x0 + x1, x0 + x1 <= 4, x0 in [1, 3], x1 free.
    fn small_highs() -> Highs {
        let mut highs = Highs::default();
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.sense_ = ObjSense::kMaximize;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![1.0, -kHighsInf];
        lp.col_upper_ = vec![3.0, kHighsInf];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![4.0];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        highs.model_lp_ = lp;
        highs
    }

    #[test]
    fn standard_form_dimensions_and_validity() {
        let mut highs = small_highs();
        let (mut nc, mut nr, mut nz, mut off) = (0, 0, 0, 0.0);
        assert_eq!(
            highs.get_standard_form_lp(
                &mut nc, &mut nr, &mut nz, &mut off, None, None, None, None, None
            ),
            HighsStatus::kOk
        );
        assert!(highs.standard_form_valid_);
        // x0, x1, the slack of x0's upper-bound row, -x1 (free split) and
        // the slack of the <= row; rows: the original and x0 + s = 3
        assert_eq!((nc, nr, nz), (5, 2, 6));
        let mut cost = vec![0.0; 5];
        let mut rhs = vec![0.0; 2];
        highs.get_standard_form_lp(
            &mut nc,
            &mut nr,
            &mut nz,
            &mut off,
            Some(&mut cost),
            Some(&mut rhs),
            None,
            None,
            None,
        );
        // Maximization is negated; x0 is shifted by its lower bound 1
        assert_eq!(cost, vec![-1.0, -1.0, 0.0, 1.0, 0.0]);
        assert_eq!(rhs, vec![3.0, 2.0]);
        assert_eq!(off, -1.0);
    }

    #[test]
    fn clear_derived_model_properties_invalidates_the_cache() {
        let mut highs = small_highs();
        assert_eq!(highs.form_standard_form_lp(), HighsStatus::kOk);
        assert!(highs.standard_form_valid_);
        highs.clear_derived_model_properties();
        assert!(!highs.standard_form_valid_);
    }
}
