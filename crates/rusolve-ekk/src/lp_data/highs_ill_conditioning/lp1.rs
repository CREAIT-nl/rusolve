//! HiGHS source: lp_data/HighsInterface.cpp (formIllConditioningLp1)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The method 1 ill-conditioning LP (Klotz14 notation, see the
//! C++ comment at lines 3517-3745): B^T y (or B y for the column view) with
//! elastic variables s, t, the split |y| = u + w, e^T y = 1 and the
//! infeasibility bound e^T s + e^T t <= `ill_conditioning_bound`.

use crate::prelude::*;

impl Highs {
    /// HiGHS: HighsInterface.cpp::Highs::formIllConditioningLp1 (lines 3510-3745).
    /// `ill_conditioning_lp` is the auxiliary instance's model LP (see mod.rs).
    pub fn form_ill_conditioning_lp1(
        &mut self,
        ill_conditioning_lp: &mut HighsLp,
        basic_var: &mut Vec<HighsInt>,
        constraint: bool,
        ill_conditioning_bound: f64,
    ) {
        // HiGHS: HighsInterface.cpp::Highs::formIllConditioningLp1 (lines 3510-3745)
        // HighsLp& incumbent_lp = this->model_.lp_;
        // const HighsInt incumbent_num_row = incumbent_lp.num_row_;
        let incumbent_num_row = self.model_lp_.num_row_;

        // Using notation from Klotz14
        //
        // For constraint view, conditioning LP minimizes the
        // infeasibilities of c7
        //
        // c4: B^Ty         -   s +   t   = 0
        // c1:    y - u + w               = 0
        // c7:        u + w               = 0
        // c6: e^Ty                       = 1
        // c5:               e^Ts + e^Tt <= eps
        // y free; u, w, s, t >= 0
        //
        // Column view uses B rather than B^T
        //
        // Set up offsets
        //
        let c4_offset: HighsInt = 0;
        let c1_offset: HighsInt = incumbent_num_row;
        let c7_offset: HighsInt = 2 * incumbent_num_row;
        let c6_offset: HighsInt = 3 * incumbent_num_row;
        let c5_offset: HighsInt = 3 * incumbent_num_row + 1;

        for i_row in 0..c6_offset {
            ill_conditioning_lp.row_lower_.push(0.0);
            ill_conditioning_lp.row_upper_.push(0.0);
        }

        // HighsSparseMatrix& incumbent_matrix = incumbent_lp.a_matrix_;
        // incumbent_matrix.ensureColwise();
        self.model_lp_.a_matrix_.ensure_colwise();

        // HighsSparseMatrix& ill_conditioning_matrix = ill_conditioning_lp.a_matrix_;
        // Form the basis matrix and
        //
        // * For constraint view, add the identity matrix and vector of
        // * ones, and transpose the resulting matrix
        //
        // * For column view, add an identity matrix column and unit entry
        // * below each column
        //
        ill_conditioning_lp.num_col_ = 0;

        // Borrow incumbent_matrix fields to avoid double-borrow with self.basis_
        let incumbent_start = &self.model_lp_.a_matrix_.start_;
        let incumbent_index = &self.model_lp_.a_matrix_.index_;
        let incumbent_value = &self.model_lp_.a_matrix_.value_;
        let incumbent_num_col = self.model_lp_.num_col_;
        let col_status = &self.basis_.col_status;

        for i_col in 0..incumbent_num_col {
            if col_status[i_col as usize] != HighsBasisStatus::kBasic {
                continue;
            }
            // Basic column goes into ill-conditioning LP, possibly with
            // identity matrix column for constraint y - u + w = 0 and unit
            // entry for e^Ty = 1
            basic_var.push(i_col);
            ill_conditioning_lp
                .col_names_
                .push(format!("y_{}", ill_conditioning_lp.num_col_));
            ill_conditioning_lp.col_cost_.push(0.0);
            ill_conditioning_lp.col_lower_.push(-kHighsInf);
            ill_conditioning_lp.col_upper_.push(kHighsInf);

            let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
            for i_el in incumbent_start[i_col as usize]..incumbent_start[(i_col + 1) as usize] {
                ill_matrix.index_.push(incumbent_index[i_el as usize]);
                ill_matrix.value_.push(incumbent_value[i_el as usize]);
            }

            if !constraint {
                // Add identity matrix column for constraint y - u + w = 0
                ill_matrix
                    .index_
                    .push(c1_offset + ill_conditioning_lp.num_col_);
                ill_matrix.value_.push(1.0);
                // Add unit entry for e^Ty = 1
                ill_matrix.index_.push(c6_offset);
                ill_matrix.value_.push(1.0);
            }
            ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            ill_conditioning_lp.num_col_ += 1;
        }

        let row_status = &self.basis_.row_status;
        let incumbent_num_col = self.model_lp_.num_col_;

        for i_row in 0..incumbent_num_row {
            if row_status[i_row as usize] != HighsBasisStatus::kBasic {
                continue;
            }
            // Basic slack goes into conditioning LP
            basic_var.push(incumbent_num_col + i_row);
            ill_conditioning_lp
                .col_names_
                .push(format!("y_{}", ill_conditioning_lp.num_col_));
            ill_conditioning_lp.col_cost_.push(0.0);
            ill_conditioning_lp.col_lower_.push(-kHighsInf);
            ill_conditioning_lp.col_upper_.push(kHighsInf);

            let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
            ill_matrix.index_.push(i_row);
            ill_matrix.value_.push(-1.0);

            if !constraint {
                // Add identity matrix column for constraint y - u + w = 0
                ill_matrix
                    .index_
                    .push(c1_offset + ill_conditioning_lp.num_col_);
                ill_matrix.value_.push(1.0);
                // Add unit entry for e^Ty = 1
                ill_matrix.index_.push(c6_offset);
                ill_matrix.value_.push(1.0);
            }
            ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            ill_conditioning_lp.num_col_ += 1;
        }

        debug_assert!(ill_conditioning_lp.num_col_ == incumbent_num_row);

        if constraint {
            // Add the identity matrix for constraint y - u + w = 0
            for i_row in 0..incumbent_num_row {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(i_row);
                ill_matrix.value_.push(1.0);
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }
            // Add the square zero matrix of c7
            for _i_row in 0..incumbent_num_row {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }
            // Add the vector of ones for e^Ty = 1
            for i_row in 0..incumbent_num_row {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(i_row);
                ill_matrix.value_.push(1.0);
            }
            {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }

            // Transpose the resulting matrix
            {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.num_col_ = c6_offset + 1;
                ill_matrix.num_row_ = incumbent_num_row;
                ill_matrix.ensure_rowwise();
                ill_matrix.format_ = MatrixFormat::kColwise;
                ill_matrix.num_col_ = incumbent_num_row;
                ill_matrix.num_row_ = c6_offset + 1;
            }
        }

        debug_assert!(ill_conditioning_lp.num_col_ == incumbent_num_row);
        ill_conditioning_lp.num_row_ = 3 * incumbent_num_row + 2;

        // Now add the variables u and w
        for i_row in 0..incumbent_num_row {
            // Adding u with cost 0
            ill_conditioning_lp.col_names_.push(format!("u_{}", i_row));
            ill_conditioning_lp.col_cost_.push(0.0);
            ill_conditioning_lp.col_lower_.push(0.0);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            // Contribution to c1: y - u + w = 0
            {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(c1_offset + i_row);
                ill_matrix.value_.push(-1.0);
                // Contribution to c7: u + w = 0
                ill_matrix.index_.push(c7_offset + i_row);
                ill_matrix.value_.push(1.0);
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }
            ill_conditioning_lp.num_col_ += 1;

            // Adding w with cost 0
            ill_conditioning_lp.col_names_.push(format!("w_{}", i_row));
            ill_conditioning_lp.col_cost_.push(0.0);
            ill_conditioning_lp.col_lower_.push(0.0);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            // Contribution to c1: y - u + w = 0
            {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(c1_offset + i_row);
                ill_matrix.value_.push(1.0);
                // Contribution to c7: u + w = 0
                ill_matrix.index_.push(c7_offset + i_row);
                ill_matrix.value_.push(1.0);
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }
            ill_conditioning_lp.num_col_ += 1;
        }

        // Now add the variables s and t
        for i_row in 0..incumbent_num_row {
            // Adding s with cost 0
            ill_conditioning_lp.col_names_.push(format!("s_{}", i_row));
            ill_conditioning_lp.col_cost_.push(0.0);
            ill_conditioning_lp.col_lower_.push(0.0);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            // Contribution to c4: B^Ty - s + t = 0
            {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(c4_offset + i_row);
                ill_matrix.value_.push(-1.0);
                // Contribution to c5: e^Ts + e^Tt <= eps
                ill_matrix.index_.push(c5_offset);
                ill_matrix.value_.push(1.0);
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }
            ill_conditioning_lp.num_col_ += 1;

            // Adding t with cost 0
            ill_conditioning_lp.col_names_.push(format!("t_{}", i_row));
            ill_conditioning_lp.col_cost_.push(0.0);
            ill_conditioning_lp.col_lower_.push(0.0);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            // Contribution to c4: B^Ty - s + t = 0
            {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(c4_offset + i_row);
                ill_matrix.value_.push(1.0);
                // Contribution to c5: e^Ts + e^Tt <= eps
                ill_matrix.index_.push(c5_offset);
                ill_matrix.value_.push(1.0);
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }
            ill_conditioning_lp.num_col_ += 1;
        }

        // Add the bounds for c6: e^Ty = 1
        ill_conditioning_lp.row_lower_.push(1.0);
        ill_conditioning_lp.row_upper_.push(1.0);
        // Add the bounds for c5: e^Ts + e^Tt <= eps
        debug_assert!(ill_conditioning_bound > 0.0);
        ill_conditioning_lp.row_lower_.push(-kHighsInf);
        ill_conditioning_lp.row_upper_.push(ill_conditioning_bound);

        debug_assert!(
            ill_conditioning_lp.row_lower_.len() as HighsInt == ill_conditioning_lp.num_row_
        );
        debug_assert!(
            ill_conditioning_lp.row_upper_.len() as HighsInt == ill_conditioning_lp.num_row_
        );

        // Now add the variables to measure the infeasibilities in
        //
        // c7: u + w = r^+ - r^-
        for i_row in 0..incumbent_num_row {
            // Adding r^+ with cost 1
            ill_conditioning_lp
                .col_names_
                .push(format!("IfsPlus_{}", i_row));
            ill_conditioning_lp.col_cost_.push(1.0);
            ill_conditioning_lp.col_lower_.push(0.0);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(c7_offset + i_row);
                ill_matrix.value_.push(-1.0);
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }
            ill_conditioning_lp.num_col_ += 1;

            // Adding r^- with cost 1
            ill_conditioning_lp
                .col_names_
                .push(format!("IfsMinus_{}", i_row));
            ill_conditioning_lp.col_cost_.push(1.0);
            ill_conditioning_lp.col_lower_.push(0.0);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(c7_offset + i_row);
                ill_matrix.value_.push(1.0);
                ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            }
            ill_conditioning_lp.num_col_ += 1;
        }

        debug_assert!(ill_conditioning_lp.num_col_ == 7 * incumbent_num_row);
        debug_assert!(ill_conditioning_lp.num_row_ == 3 * incumbent_num_row + 2);

        {
            let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
            ill_matrix.num_col_ = ill_conditioning_lp.num_col_;
            ill_matrix.num_row_ = ill_conditioning_lp.num_row_;
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// One column x0 (basic, a = 2) and one row (nonbasic), column view:
    /// 7 columns (y, u, w, s, t, r+, r-) and 3 * 1 + 2 rows.
    #[test]
    fn column_view_lp1_layout() {
        let mut highs = Highs::default();
        let lp = &mut highs.model_lp_;
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![1.0];
        lp.row_lower_ = vec![0.0];
        lp.row_upper_ = vec![1.0];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![2.0];
        highs.basis_.col_status = vec![HighsBasisStatus::kBasic];
        highs.basis_.row_status = vec![HighsBasisStatus::kLower];
        let mut ill_lp = HighsLp::default();
        let mut basic_var = Vec::new();
        highs.form_ill_conditioning_lp1(&mut ill_lp, &mut basic_var, false, 1e-4);
        assert_eq!(basic_var, vec![0]);
        assert_eq!((ill_lp.num_col_, ill_lp.num_row_), (7, 5));
        assert_eq!(ill_lp.a_matrix_.start_, vec![0, 3, 5, 7, 9, 11, 12, 13]);
        assert_eq!(
            ill_lp.a_matrix_.index_,
            vec![0, 1, 3, 1, 2, 1, 2, 0, 4, 0, 4, 2, 2]
        );
        assert_eq!(ill_lp.col_cost_, vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]);
        assert_eq!(ill_lp.row_lower_, vec![0.0, 0.0, 0.0, 1.0, -kHighsInf]);
        assert_eq!(ill_lp.row_upper_, vec![0.0, 0.0, 0.0, 1.0, 1e-4]);
        assert_eq!(ill_lp.col_names_[5], "IfsPlus_0");
    }
}
