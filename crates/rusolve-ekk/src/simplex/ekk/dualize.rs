//! HiGHS source: simplex/HEkk.cpp (lines 452-1000)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Dualize/undualize and permute/unpermute of the LP.

use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkk.cpp::HEkk::dualize (lines 452-755, 304 lines)
    pub fn dualize(&mut self) -> HighsStatus {
        debug_assert!(self.lp_.a_matrix_.is_colwise());
        self.original_num_col_ = self.lp_.num_col_;
        self.original_num_row_ = self.lp_.num_row_;
        self.original_num_nz_ = self.lp_.a_matrix_.num_nz();
        self.original_offset_ = self.lp_.offset_;
        self.original_col_cost_ = self.lp_.col_cost_.clone();
        self.original_col_lower_ = self.lp_.col_lower_.clone();
        self.original_col_upper_ = self.lp_.col_upper_.clone();
        self.original_row_lower_ = self.lp_.row_lower_.clone();
        self.original_row_upper_ = self.lp_.row_upper_.clone();
        // Reserve space for simple dual
        self.lp_.col_cost_.reserve(self.original_num_row_ as usize);
        self.lp_.col_lower_.reserve(self.original_num_row_ as usize);
        self.lp_.col_upper_.reserve(self.original_num_row_ as usize);
        self.lp_.row_lower_.reserve(self.original_num_col_ as usize);
        self.lp_.row_upper_.reserve(self.original_num_col_ as usize);
        // Invalidate the original data
        self.lp_.col_cost_.resize(0, 0.0);
        self.lp_.col_lower_.resize(0, 0.0);
        self.lp_.col_upper_.resize(0, 0.0);
        self.lp_.row_lower_.resize(0, 0.0);
        self.lp_.row_upper_.resize(0, 0.0);
        // The bulk of the constraint matrix of the dual LP is the transpose
        // of the primal constraint matrix. This is obtained row-wise by
        // copying the matrix and flipping the dimensions
        let mut dual_matrix = self.lp_.a_matrix_.clone();
        dual_matrix.num_row_ = self.original_num_col_;
        dual_matrix.num_col_ = self.original_num_row_;
        dual_matrix.format_ = MatrixFormat::kRowwise;
        // The primal_bound_value vector accumulates the values of the
        // finite bounds on variables - or zero for a free variable - used
        // later to compute the offset and shift for the costs. Many of
        // these components will be zero - all for the case x>=0 - so
        // maintain a list of the nonzeros and corresponding indices. Don't
        // reserve space since they may not be needed
        let mut primal_bound_value: Vec<f64> = Vec::new();
        let mut primal_bound_index: Vec<HighsInt> = Vec::new();
        let inf: f64 = kHighsInf;
        for i_col in 0..self.original_num_col_ {
            let cost = self.original_col_cost_[i_col as usize];
            let lower = self.original_col_lower_[i_col as usize];
            let upper = self.original_col_upper_[i_col as usize];
            let primal_bound: f64;
            let row_lower: f64;
            let row_upper: f64;
            if lower == upper {
                // Fixed
                primal_bound = lower;
                // Dual activity a^Ty is free, implying dual for primal column
                // (slack for dual row) is free
                row_lower = -inf;
                row_upper = inf;
            } else if !highs_is_infinity(-lower) {
                // Finite lower bound so boxed or lower
                if !highs_is_infinity(upper) {
                    // Finite upper bound so boxed
                    //
                    // Treat as lower
                    primal_bound = lower;
                    // Dual activity a^Ty is bounded above by cost, implying dual
                    // for primal column (slack for dual row) is non-negative
                    row_lower = -inf;
                    row_upper = cost;
                    // Treat upper bound as additional constraint
                    self.upper_bound_col_.push(i_col);
                } else {
                    // Lower (since upper bound is infinite)
                    primal_bound = lower;
                    // Dual activity a^Ty is bounded above by cost, implying dual
                    // for primal column (slack for dual row) is non-negative
                    row_lower = -inf;
                    row_upper = cost;
                }
            } else if !highs_is_infinity(upper) {
                // Upper
                primal_bound = upper;
                // Dual activity a^Ty is bounded below by cost, implying dual
                // for primal column (slack for dual row) is non-positive
                row_lower = cost;
                row_upper = inf;
            } else {
                // FREE
                //
                // Dual activity a^Ty is fixed by cost, implying dual for primal
                // column (slack for dual row) is fixed at zero
                primal_bound = 0.0;
                row_lower = cost;
                row_upper = cost;
            }
            debug_assert!(row_lower < inf);
            debug_assert!(row_upper > -inf);
            debug_assert!(primal_bound < inf);
            self.lp_.row_lower_.push(row_lower);
            self.lp_.row_upper_.push(row_upper);
            if primal_bound != 0.0 {
                primal_bound_value.push(primal_bound);
                primal_bound_index.push(i_col);
            }
        }
        for i_row in 0..self.original_num_row_ {
            let lower = self.original_row_lower_[i_row as usize];
            let upper = self.original_row_upper_[i_row as usize];
            let col_cost: f64;
            let col_lower: f64;
            let col_upper: f64;
            if lower == upper {
                // Equality constraint
                //
                // Dual variable has primal RHS as cost and is free
                col_cost = lower;
                col_lower = -inf;
                col_upper = inf;
            } else if !highs_is_infinity(-lower) {
                // Finite lower bound so boxed or lower
                if !highs_is_infinity(upper) {
                    // Finite upper bound so boxed
                    //
                    // Treat as lower
                    col_cost = lower;
                    col_lower = 0.0;
                    col_upper = inf;
                    // Treat upper bound as additional constraint
                    self.upper_bound_row_.push(i_row);
                } else {
                    // Lower (since upper bound is infinite)
                    col_cost = lower;
                    col_lower = 0.0;
                    col_upper = inf;
                }
            } else if !highs_is_infinity(upper) {
                // Upper
                col_cost = upper;
                col_lower = -inf;
                col_upper = 0.0;
            } else {
                // FREE
                // Shouldn't get free rows, but handle them anyway
                col_cost = 0.0;
                col_lower = 0.0;
                col_upper = 0.0;
            }
            debug_assert!(col_lower < inf);
            debug_assert!(col_upper > -inf);
            debug_assert!(col_cost < inf);
            self.lp_.col_cost_.push(col_cost);
            self.lp_.col_lower_.push(col_lower);
            self.lp_.col_upper_.push(col_upper);
        }
        // Boxed variables and constraints yield extra columns in the dual LP
        let mut extra_columns = HighsSparseMatrix::default();
        extra_columns.ensure_colwise();
        extra_columns.num_row_ = self.original_num_col_;
        let num_upper_bound_col: HighsInt = self.upper_bound_col_.len() as HighsInt;
        let num_upper_bound_row: HighsInt = self.upper_bound_row_.len() as HighsInt;
        let one: f64 = 1.0;
        for i_x in 0..num_upper_bound_col {
            let i_col = self.upper_bound_col_[i_x as usize];
            let upper = self.original_col_upper_[i_col as usize];
            extra_columns.add_vec(
                1,
                std::slice::from_ref(&i_col),
                std::slice::from_ref(&one),
                1.0,
            );
            self.lp_.col_cost_.push(upper);
            self.lp_.col_lower_.push(-inf);
            self.lp_.col_upper_.push(0.0);
        }

        if num_upper_bound_row != 0 {
            // Need to identify the submatrix of constraint matrix rows
            // corresponding to those with a row index in
            // upper_bound_row_. When identifying numbers of entries in each
            // row of submatrix, use indirection to get corresponding row
            // index, with a dummy row for rows not in the submatrix.
            let dummy_row: HighsInt = num_upper_bound_row;
            let mut indirection: Vec<HighsInt> = vec![dummy_row; self.original_num_row_ as usize];
            let mut count: Vec<HighsInt> = vec![0; (num_upper_bound_row + 1) as usize];
            let mut extra_i_row: HighsInt = 0;
            #[allow(clippy::explicit_counter_loop)]
            for i_x in 0..num_upper_bound_row {
                let i_row = self.upper_bound_row_[i_x as usize];
                indirection[i_row as usize] = extra_i_row;
                extra_i_row += 1;
                let upper = self.original_row_upper_[i_row as usize];
                self.lp_.col_cost_.push(upper);
                self.lp_.col_lower_.push(-inf);
                self.lp_.col_upper_.push(0.0);
            }
            for i_el in 0..self.original_num_nz_ {
                let idx = self.lp_.a_matrix_.index_[i_el as usize];
                count[indirection[idx as usize] as usize] += 1;
            }
            extra_columns
                .start_
                .resize((num_upper_bound_col + num_upper_bound_row + 1) as usize, 0);
            for i_row in 0..num_upper_bound_row {
                extra_columns.start_[(num_upper_bound_col + i_row + 1) as usize] = extra_columns
                    .start_[(num_upper_bound_col + i_row) as usize]
                    + count[i_row as usize];
                count[i_row as usize] =
                    extra_columns.start_[(num_upper_bound_col + i_row) as usize];
            }
            let extra_columns_num_nz: HighsInt =
                extra_columns.start_[(num_upper_bound_col + num_upper_bound_row) as usize];
            extra_columns
                .index_
                .resize(extra_columns_num_nz as usize, 0);
            extra_columns
                .value_
                .resize(extra_columns_num_nz as usize, 0.0);
            for i_col in 0..self.original_num_col_ {
                let i_start = self.lp_.a_matrix_.start_[i_col as usize];
                let i_end = self.lp_.a_matrix_.start_[(i_col + 1) as usize];
                for i_el in i_start..i_end {
                    let i_row = indirection[self.lp_.a_matrix_.index_[i_el as usize] as usize];
                    if i_row < num_upper_bound_row {
                        let extra_columns_i_el = count[i_row as usize];
                        debug_assert!(extra_columns_i_el < extra_columns_num_nz);
                        extra_columns.index_[extra_columns_i_el as usize] = i_col;
                        extra_columns.value_[extra_columns_i_el as usize] =
                            self.lp_.a_matrix_.value_[i_el as usize];
                        count[i_row as usize] += 1;
                    }
                }
            }
            extra_columns.num_col_ += num_upper_bound_row;
        }
        // Incorporate the cost shift by subtracting A*primal_bound from the
        // cost vector; compute the objective offset
        let mut delta_offset: f64 = 0.0;
        for i_x in 0..primal_bound_index.len() {
            let i_col = primal_bound_index[i_x];
            let multiplier = primal_bound_value[i_x];
            delta_offset += multiplier * self.original_col_cost_[i_col as usize];
            let i_start = self.lp_.a_matrix_.start_[i_col as usize];
            let i_end = self.lp_.a_matrix_.start_[(i_col + 1) as usize];
            for i_el in i_start..i_end {
                let row = self.lp_.a_matrix_.index_[i_el as usize];
                self.lp_.col_cost_[row as usize] -=
                    multiplier * self.lp_.a_matrix_.value_[i_el as usize];
            }
        }
        if extra_columns.num_col_ != 0 {
            // Incorporate the cost shift by subtracting
            // extra_columns*primal_bound from the cost vector for the extra
            // dual variables
            //
            // Have to scatter the packed primal bound values into a
            // full-length vector
            //
            // ToDo Make this more efficient?
            let mut primal_bound: Vec<f64> = vec![0.0; self.original_num_col_ as usize];
            for i_x in 0..primal_bound_index.len() {
                primal_bound[primal_bound_index[i_x] as usize] = primal_bound_value[i_x];
            }

            for i_col in 0..extra_columns.num_col_ {
                let mut cost = self.lp_.col_cost_[(self.original_num_row_ + i_col) as usize];
                for i_el in
                    extra_columns.start_[i_col as usize]..extra_columns.start_[(i_col + 1) as usize]
                {
                    cost -= primal_bound[extra_columns.index_[i_el as usize] as usize]
                        * extra_columns.value_[i_el as usize];
                }
                self.lp_.col_cost_[(self.original_num_row_ + i_col) as usize] = cost;
            }
        }
        self.lp_.offset_ += delta_offset;
        // Copy the row-wise dual LP constraint matrix and transpose it.
        // ToDo Make this more efficient
        self.lp_.a_matrix_ = dual_matrix;
        self.lp_.a_matrix_.ensure_colwise();
        // Add the extra columns to the dual LP constraint matrix
        self.lp_.a_matrix_.add_cols(&extra_columns, None);

        let dual_num_col: HighsInt =
            self.original_num_row_ + num_upper_bound_col + num_upper_bound_row;
        let dual_num_row: HighsInt = self.original_num_col_;
        debug_assert!(dual_num_col == self.lp_.col_cost_.len() as HighsInt);
        debug_assert!(self.lp_.a_matrix_.num_col_ == dual_num_col);
        let ignore_scaling = true;
        if !ignore_scaling {
            // Flip any scale factors
            if self.lp_.scale_.has_scaling {
                let temp_scale = self.lp_.scale_.row.clone();
                self.lp_.scale_.row = self.lp_.scale_.col.clone();
                self.lp_.scale_.col = temp_scale;
                self.lp_.scale_.num_col = dual_num_col;
                self.lp_.scale_.num_row = dual_num_row;
            }
        }
        // Change optimization sense
        if self.lp_.sense_ == ObjSense::kMinimize {
            self.lp_.sense_ = ObjSense::kMaximize;
        } else {
            self.lp_.sense_ = ObjSense::kMinimize;
        }
        // Flip LP dimensions
        self.lp_.num_col_ = dual_num_col;
        self.lp_.num_row_ = dual_num_row;
        self.status_.is_dualized = true;
        self.status_.has_basis = false;
        self.status_.has_ar_matrix = false;
        self.status_.has_nla = false;
        // PORT FIX: edge weights and INVERT left over from a
        // previous solve of the primal LP have the primal dimensions; HiGHS keeps
        // them flagged valid (irish-electricity: index out of bounds in
        // putBacktrackingBasis on a warm EKK).
        self.status_.has_dual_steepest_edge_weights = false;
        self.status_.has_invert = false;
        self.status_.has_fresh_invert = false;
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "Solving dual LP with {} columns",
            dual_num_col
        );
        if num_upper_bound_col + num_upper_bound_row != 0 {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                " [{} extra from",
                dual_num_col - self.original_num_row_
            );
            if num_upper_bound_col != 0 {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    " {} boxed variable(s)",
                    num_upper_bound_col
                );
            }
            if num_upper_bound_col != 0 && num_upper_bound_row != 0 {
                crate::highs_log_user!(&self.options_.log_options, HighsLogType::kInfo, " and");
            }
            if num_upper_bound_row != 0 {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    " {} boxed constraint(s)",
                    num_upper_bound_row
                );
            }
            crate::highs_log_user!(&self.options_.log_options, HighsLogType::kInfo, "]");
        }
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            " and {} rows\n",
            dual_num_row
        );
        //  reportLp(options_->log_options, lp_, HighsLogType::kVerbose);
        HighsStatus::kOk
    }

    /// HiGHS: HEkk.cpp::HEkk::undualize (lines 757-989, 233 lines)
    #[allow(clippy::assertions_on_constants)]
    pub fn undualize(&mut self) -> HighsStatus {
        if !self.status_.is_dualized {
            return HighsStatus::kOk;
        }
        let dual_num_col: HighsInt = self.lp_.num_col_;
        let primal_num_tot: HighsInt = self.original_num_col_ + self.original_num_row_;
        // These two aren't used (yet)
        // vector<double>& dual_work_dual = info_.workDual_;
        // vector<double>& primal_work_value = info_.workValue_;
        // Take copies of the nonbasic information for the dual LP, since
        // its values will be over-written in constructing the corresponding
        // data for the primal problem
        let dual_nonbasic_flag: Vec<i8> = self.basis_.nonbasic_flag_.clone();
        let dual_nonbasic_move: Vec<i8> = self.basis_.nonbasic_move_.clone();
        // primal_basic_index, primal_nonbasic_flag and primal_nonbasic_move are
        // aliases of basis_.basicIndex_, basis_.nonbasicFlag_ and
        // basis_.nonbasicMove_ in C++; they are accessed through the basis below.
        self.basis_.nonbasic_flag_ = vec![kIllegalFlagValue; primal_num_tot as usize];
        self.basis_.nonbasic_move_ = vec![kIllegalMoveValue; primal_num_tot as usize];
        self.basis_.basic_index_.clear();
        // The number of dual rows is the number of primal columns, so all
        // dual basic variables are nonbasic in the primal problem.
        //
        // If there are extra dual columns due to upper bounds on boxed
        // primal variables/constraints, there will be an excess of dual
        // nonbasic variables for the required primal basic variables.
        //
        // For each pair of dual variables associated with a boxed primal
        // variable/constraint:
        //
        // * If one is basic then it yields a nonbasic primal variable, at
        //   a bound given by the basic dual
        //
        // * If both are nonbasic, then they yield a basic primal variable
        //
        // Keep track of the dual column added to handle upper bounds on
        // boxed variables/constraints
        let mut upper_bound_col: HighsInt = self.original_num_row_;
        for i_col in 0..self.original_num_col_ {
            let lower = self.original_col_lower_[i_col as usize];
            let upper = self.original_col_upper_[i_col as usize];
            let mut r#move: i8 = kIllegalMoveValue;
            let mut dual_variable: HighsInt = dual_num_col + i_col;
            let mut dual_basic: bool =
                dual_nonbasic_flag[dual_variable as usize] == kNonbasicFlagFalse;
            if lower == upper {
                // Fixed
                if dual_basic {
                    r#move = kNonbasicMoveZe;
                }
            } else if !highs_is_infinity(-lower) {
                // Finite lower bound so boxed or lower
                if !highs_is_infinity(upper) {
                    // Finite upper bound so boxed
                    if dual_basic {
                        // Primal variable is nonbasic at its lower bound
                        r#move = kNonbasicMoveUp;
                    } else {
                        // Look at the corresponding dual variable for the upper bound
                        dual_variable = upper_bound_col;
                        dual_basic =
                            dual_nonbasic_flag[dual_variable as usize] == kNonbasicFlagFalse;
                        if dual_basic {
                            // Primal variable is nonbasic at its upper bound
                            r#move = kNonbasicMoveDn;
                        }
                    }
                    upper_bound_col += 1;
                } else {
                    // Lower (since upper bound is infinite)
                    if dual_basic {
                        r#move = kNonbasicMoveUp;
                    }
                }
            } else if !highs_is_infinity(upper) {
                // Upper
                if dual_basic {
                    r#move = kNonbasicMoveDn;
                }
            } else {
                // FREE
                //
                // Dual activity a^Ty is fixed by cost, implying dual for primal
                // column (slack for dual row) is fixed at zero
                if dual_basic {
                    debug_assert!(4 == 0);
                    r#move = kNonbasicMoveZe;
                }
            }
            if dual_basic {
                // Primal nonbasic column from basic dual row
                debug_assert!(r#move != kIllegalMoveValue);
                self.basis_.nonbasic_flag_[i_col as usize] = kNonbasicFlagTrue;
                self.basis_.nonbasic_move_[i_col as usize] = r#move;
            } else {
                // Primal basic column from nonbasic dual row
                self.basis_.basic_index_.push(i_col);
                self.basis_.nonbasic_flag_[i_col as usize] = kNonbasicFlagFalse;
                self.basis_.nonbasic_move_[i_col as usize] = 0;
            }
        }
        for i_row in 0..self.original_num_row_ {
            let lower = self.original_row_lower_[i_row as usize];
            let upper = self.original_row_upper_[i_row as usize];
            let mut r#move: i8 = kIllegalMoveValue;
            let mut dual_variable: HighsInt = i_row;
            let mut dual_basic: bool =
                dual_nonbasic_flag[dual_variable as usize] == kNonbasicFlagFalse;
            if lower == upper {
                // Equality constraint
                //
                // Dual variable has primal RHS as cost and is free
                if dual_basic {
                    r#move = kNonbasicMoveZe;
                }
            } else if !highs_is_infinity(-lower) {
                // Finite lower bound so boxed or lower
                if !highs_is_infinity(upper) {
                    // Finite upper bound so boxed
                    if dual_basic {
                        // Primal variable is nonbasic at its lower bound
                        r#move = kNonbasicMoveDn;
                    } else {
                        // Look at the corresponding dual variable for the upper bound
                        dual_variable = upper_bound_col;
                        dual_basic =
                            dual_nonbasic_flag[dual_variable as usize] == kNonbasicFlagFalse;
                        if dual_basic {
                            // Primal variable is nonbasic at its upper bound
                            r#move = kNonbasicMoveUp;
                        }
                    }
                    upper_bound_col += 1;
                } else {
                    // Lower (since upper bound is infinite)
                    if dual_basic {
                        r#move = kNonbasicMoveDn;
                    }
                }
            } else if !highs_is_infinity(upper) {
                // Upper
                if dual_basic {
                    r#move = kNonbasicMoveUp;
                }
            } else {
                // FREE
                if dual_basic {
                    debug_assert!(14 == 0);
                    r#move = kNonbasicMoveZe;
                }
                // Shouldn't get free rows, but handle them anyway
            }
            if dual_basic {
                // Primal nonbasic column from basic dual row
                debug_assert!(r#move != kIllegalMoveValue);
                self.basis_.nonbasic_flag_[(self.original_num_col_ + i_row) as usize] =
                    kNonbasicFlagTrue;
                self.basis_.nonbasic_move_[(self.original_num_col_ + i_row) as usize] = r#move;
            } else {
                // Primal basic column from nonbasic dual row
                self.basis_
                    .basic_index_
                    .push(self.original_num_col_ + i_row);
                self.basis_.nonbasic_flag_[(self.original_num_col_ + i_row) as usize] =
                    kNonbasicFlagFalse;
                self.basis_.nonbasic_move_[(self.original_num_col_ + i_row) as usize] = 0;
            }
        }
        let ignore_scaling = true;
        if !ignore_scaling {
            // Flip any scale factors
            if self.lp_.scale_.has_scaling {
                let temp_scale = self.lp_.scale_.row.clone();
                self.lp_.scale_.row = self.lp_.scale_.col.clone();
                self.lp_.scale_.col = temp_scale;
                self.lp_
                    .scale_
                    .col
                    .resize(self.original_num_col_ as usize, 0.0);
                self.lp_
                    .scale_
                    .row
                    .resize(self.original_num_row_ as usize, 0.0);
                self.lp_.scale_.num_col = self.original_num_col_;
                self.lp_.scale_.num_row = self.original_num_row_;
            }
        }
        // Change optimization sense
        if self.lp_.sense_ == ObjSense::kMinimize {
            self.lp_.sense_ = ObjSense::kMaximize;
        } else {
            self.lp_.sense_ = ObjSense::kMinimize;
        }
        // Flip LP dimensions
        self.lp_.num_col_ = self.original_num_col_;
        self.lp_.num_row_ = self.original_num_row_;
        // Restore the original offset
        self.lp_.offset_ = self.original_offset_;
        // Copy back the costs and bounds
        self.lp_.col_cost_ = self.original_col_cost_.clone();
        self.lp_.col_lower_ = self.original_col_lower_.clone();
        self.lp_.col_upper_ = self.original_col_upper_.clone();
        self.lp_.row_lower_ = self.original_row_lower_.clone();
        self.lp_.row_upper_ = self.original_row_upper_.clone();
        // The primal constraint matrix is available row-wise as the first
        // original_num_row_ vectors of the dual constraint matrix
        let mut primal_matrix = HighsSparseMatrix::default();
        primal_matrix
            .start_
            .resize((self.original_num_row_ + 1) as usize, 0);
        primal_matrix
            .index_
            .resize(self.original_num_nz_ as usize, 0);
        primal_matrix
            .value_
            .resize(self.original_num_nz_ as usize, 0.0);

        for i_col in 0..self.original_num_row_ + 1 {
            primal_matrix.start_[i_col as usize] = self.lp_.a_matrix_.start_[i_col as usize];
        }
        for i_el in 0..self.original_num_nz_ {
            primal_matrix.index_[i_el as usize] = self.lp_.a_matrix_.index_[i_el as usize];
            primal_matrix.value_[i_el as usize] = self.lp_.a_matrix_.value_[i_el as usize];
        }
        primal_matrix.num_col_ = self.original_num_col_;
        primal_matrix.num_row_ = self.original_num_row_;
        primal_matrix.format_ = MatrixFormat::kRowwise;
        // Copy the row-wise primal LP constraint matrix and transpose it.
        // ToDo Make this more efficient
        self.lp_.a_matrix_ = primal_matrix;
        self.lp_.a_matrix_.ensure_colwise();
        // Some sanity checks
        debug_assert!(self.lp_.num_col_ == self.original_num_col_);
        debug_assert!(self.lp_.num_row_ == self.original_num_row_);
        debug_assert!(self.lp_.a_matrix_.num_nz() == self.original_num_nz_);
        let num_basic_variables: HighsInt = self.basis_.basic_index_.len() as HighsInt;
        let num_basic_variables_ok = num_basic_variables == self.original_num_row_;
        if !num_basic_variables_ok {
            println!(
                "HEkk::undualize: Have {} basic variables, not {}",
                num_basic_variables, self.original_num_row_
            );
        }
        debug_assert!(num_basic_variables_ok);

        // Clear the data retained when solving dual LP
        self.clear_ekk_dualize();
        self.status_.is_dualized = false;
        // Now solve with this basis. Should just be a case of reinverting
        // and re-solving for optimal primal and dual values, but
        // numerically marginal LPs will need clean-up
        self.status_.has_basis = true;
        self.status_.has_ar_matrix = false;
        self.status_.has_nla = false;
        self.status_.has_invert = false;
        // PORT FIX: HiGHS leaves has_fresh_invert set by the
        // dual solve, so computeFactor() returns early and the new NLA is never
        // inverted (highspy 1.15 logs the same stale-size warning on ns1688926).
        self.status_.has_fresh_invert = false;
        // The DSE weights left by the dual solve have the dual dimensions
        self.status_.has_dual_steepest_edge_weights = false;
        let mut primal_solve_iteration_count: HighsInt = -self.iteration_count_;
        let return_status = self.solve(false);
        primal_solve_iteration_count += self.iteration_count_;
        //  if (primal_solve_iteration_count)
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "Solving the primal LP ({}) using the optimal basis of its dual required {} simplex iterations\n",
            self.lp_.model_name_,
            primal_solve_iteration_count
        );
        return_status
    }

    /// HiGHS: HEkk.cpp::HEkk::permute (lines 991-994, 4 lines)
    #[allow(clippy::assertions_on_constants)]
    pub fn permute(&mut self) -> HighsStatus {
        debug_assert!(1 == 0);
        HighsStatus::kError
    }

    /// HiGHS: HEkk.cpp::HEkk::unpermute (lines 996-1000, 5 lines)
    #[allow(clippy::assertions_on_constants)]
    pub fn unpermute(&mut self) -> HighsStatus {
        if !self.status_.is_permuted {
            return HighsStatus::kOk;
        }
        debug_assert!(1 == 0);
        HighsStatus::kError
    }
}
