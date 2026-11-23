//! HiGHS source: lp_data/Highs.cpp (Highs::addCols, Highs::deleteCols
//! interval), lp_data/HighsInterface.cpp (addColsInterface,
//! appendNonbasicColsToBasisInterface, deleteBasisCols,
//! deleteColsInterface), lp_data/HighsLp.cpp (HighsLp::addColNames,
//! HighsLp::deleteColsFromVectors, HighsLp::deleteCols)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The column add/delete API that `HighsIis::compute` and
//! `Highs::elasticityFilter` drive. Mirrors the row versions in
//! `highs_api.rs` / `highs_interface_rows.rs` / `highs_lp_modify.rs`.
//! Dropped as there: `logHeader`, `isColDataNull` / `isMatrixDataNull`
//! (slices are never NULL), `col_hash_` (no name hash in this port).
//! `model_.hessian_.deleteCols` is `self.model_hessian_.delete_cols`.

use crate::lp_data::highs_interface_rows::delete_basis_entries;
use crate::lp_data::lp_assess::{assess_bounds, assess_costs, lp_dimensions_ok};
use crate::prelude::*;
use crate::util::highs_utils::create_interval;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::addCols (lines 2855-2870)
    ///
    /// `starts` has `num_new_col` entries (CSC without the end sentinel).
    pub fn add_cols(
        &mut self,
        num_new_col: HighsInt,
        costs: &[f64],
        lower_bounds: &[f64],
        upper_bounds: &[f64],
        num_new_nz: HighsInt,
        starts: &[HighsInt],
        indices: &[HighsInt],
        values: &[f64],
    ) -> HighsStatus {
        // logHeader dropped (stub)
        let mut return_status = HighsStatus::kOk;
        self.clear_derived_model_properties();
        let call_status = self.add_cols_interface(
            num_new_col,
            costs,
            lower_bounds,
            upper_bounds,
            num_new_nz,
            starts,
            indices,
            values,
        );
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "addCols",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: HighsInterface.cpp::Highs::addColsInterface (lines 395-524)
    /// HiGHS 1.16 (#3386): HighsInterface.cpp::addColsInterface (lines 394-526)
    pub(crate) fn add_cols_interface(
        &mut self,
        ext_num_new_col: HighsInt,
        ext_col_cost: &[f64],
        ext_col_lower: &[f64],
        ext_col_upper: &[f64],
        ext_num_new_nz: HighsInt,
        ext_a_start: &[HighsInt],
        ext_a_index: &[HighsInt],
        ext_a_value: &[f64],
    ) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;
        // HighsOptions& options = options_;
        if ext_num_new_col < 0 {
            return HighsStatus::kError;
        }
        if ext_num_new_nz < 0 {
            return HighsStatus::kError;
        }
        if ext_num_new_col == 0 {
            return HighsStatus::kOk;
        }
        // isColDataNull / isMatrixDataNull dropped: slices are never NULL.
        // HighsLp& lp = model_.lp_;
        // HighsBasis& basis = basis_;
        // HighsScale& scale = lp.scale_;
        // bool& useful_basis = basis.useful;
        // bool& lp_has_scaling = lp.scale_.has_scaling;

        // Check that if nonzeros are to be added then the model has a positive number
        // of rows
        if self.model_lp_.num_row_ <= 0 && ext_num_new_nz > 0 {
            return HighsStatus::kError;
        }

        // Record the new number of columns
        let new_num_col = self.model_lp_.num_col_ + ext_num_new_col;

        let mut index_collection = HighsIndexCollection::default();
        index_collection.dimension_ = ext_num_new_col;
        index_collection.is_interval_ = true;
        index_collection.from_ = 0;
        index_collection.to_ = ext_num_new_col - 1;

        // Take a copy of the cost and bounds that can be normalised
        // {ext_col_cost, ext_col_cost + ext_num_new_col}: exact lengths as C++.
        let mut local_col_cost: Vec<f64> = ext_col_cost[..ext_num_new_col as usize].to_vec();
        let mut local_col_lower: Vec<f64> = ext_col_lower[..ext_num_new_col as usize].to_vec();
        let mut local_col_upper: Vec<f64> = ext_col_upper[..ext_num_new_col as usize].to_vec();

        let mut local_has_infinite_cost = false;
        let call_status = assess_costs(
            &self.options_,
            self.model_lp_.num_col_,
            &index_collection,
            &mut local_col_cost,
            &mut local_has_infinite_cost,
            self.options_.infinite_cost,
        );
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "assessCosts",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        // Assess the column bounds
        let call_status = assess_bounds(
            &self.options_,
            "Col",
            self.model_lp_.num_col_,
            &index_collection,
            &mut local_col_lower,
            &mut local_col_upper,
            self.options_.infinite_bound,
            None,
        );
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "assessBounds",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        // HiGHS 1.16 #3386 (h16_add_no_modify_on_error): HighsInterface.cpp lines 452-487 (1.16.0)
        let h16_add_no_modify_on_error = self.options_.h16_add_no_modify_on_error;
        // Append the columns to the LP vectors and matrix
        if !h16_add_no_modify_on_error {
            append_cols_to_lp_vectors(
                &mut self.model_lp_,
                ext_num_new_col,
                &local_col_cost,
                &local_col_lower,
                &local_col_upper,
            );
        }
        // Form a column-wise HighsSparseMatrix of the new matrix columns so
        // that is easy to handle and, if there are nonzeros, it can be
        // normalised
        let mut local_a_matrix = HighsSparseMatrix::default();
        local_a_matrix.num_col_ = ext_num_new_col;
        local_a_matrix.num_row_ = self.model_lp_.num_row_;
        local_a_matrix.format_ = MatrixFormat::kColwise;
        if ext_num_new_nz != 0 {
            // {ext_a_start, ext_a_start + ext_num_new_col}: exact lengths as C++.
            local_a_matrix.start_ = ext_a_start[..ext_num_new_col as usize].to_vec();
            local_a_matrix
                .start_
                .resize(ext_num_new_col as usize + 1, 0);
            local_a_matrix.start_[ext_num_new_col as usize] = ext_num_new_nz;
            local_a_matrix.index_ = ext_a_index[..ext_num_new_nz as usize].to_vec();
            local_a_matrix.value_ = ext_a_value[..ext_num_new_nz as usize].to_vec();
            // Assess the matrix rows
            let call_status = local_a_matrix.assess(
                &self.options_.log_options,
                "LP",
                self.options_.small_matrix_value,
                self.options_.large_matrix_value,
                self.options_.h16_assess_small,
            );
            return_status = interpret_call_status(
                &self.options_.log_options,
                call_status,
                return_status,
                "assessMatrix",
            );
            if return_status == HighsStatus::kError {
                return return_status;
            }
        } else {
            // No nonzeros so, whether the constraint matrix is column-wise or
            // row-wise, adding the empty matrix is trivial. Complete the
            // setup of an empty column-wise HighsSparseMatrix of the new
            // matrix columns
            local_a_matrix.start_ = vec![0; ext_num_new_col as usize + 1];
        }
        // HiGHS 1.16 #3386 (h16_add_no_modify_on_error): HighsInterface.cpp lines 484-489 (1.16.0)
        if h16_add_no_modify_on_error {
            debug_assert!(return_status != HighsStatus::kError);
            // Now that data for new columns has been assessed and no errors can
            // occur, append the columns to the LP vectors and matrix
            append_cols_to_lp_vectors(
                &mut self.model_lp_,
                ext_num_new_col,
                &local_col_cost,
                &local_col_lower,
                &local_col_upper,
            );
        }
        // Append the columns to LP matrix
        self.model_lp_.a_matrix_.add_cols(&local_a_matrix, None);
        if self.model_lp_.scale_.has_scaling {
            // Extend the column scaling factors
            self.model_lp_.scale_.col.resize(new_num_col as usize, 0.0);
            for i_col in 0..ext_num_new_col {
                self.model_lp_.scale_.col[(self.model_lp_.num_col_ + i_col) as usize] = 1.0;
            }
            self.model_lp_.scale_.num_col = new_num_col;
            // Apply the existing row scaling to the new columns
            local_a_matrix.apply_row_scale(&self.model_lp_.scale_);
            // Consider applying column scaling to the new columns.
            local_a_matrix.consider_col_scaling(
                self.options_.allowed_matrix_scale_factor,
                &mut self.model_lp_.scale_.col[self.model_lp_.num_col_ as usize..],
            );
        }
        // Update the basis corresponding to new nonbasic columns
        if self.basis_.useful {
            self.append_nonbasic_cols_to_basis_interface(ext_num_new_col);
        }

        // Possibly add blank column names
        self.model_lp_.add_col_names("", ext_num_new_col);

        // Increase the number of columns in the LP
        self.model_lp_.num_col_ += ext_num_new_col;
        debug_assert!(lp_dimensions_ok(
            "addCols",
            &self.model_lp_,
            &self.options_.log_options
        ));

        // Interpret possible introduction of infinite costs
        self.model_lp_.has_infinite_cost_ =
            self.model_lp_.has_infinite_cost_ || local_has_infinite_cost;
        debug_assert!(
            self.model_lp_.has_infinite_cost_
                == self
                    .model_lp_
                    .has_infinite_cost(self.options_.infinite_cost)
        );

        // Deduce the consequences of adding new columns
        self.invalidate_model_status_solution_and_info();

        // Determine any implications for simplex data
        self.ekk_instance_
            .add_cols(&self.model_lp_, &local_a_matrix);

        // Extend any Hessian with zeros on the diagonal
        if self.model_hessian_.dim_ != 0 {
            complete_hessian(self.model_lp_.num_col_, &mut self.model_hessian_);
        }
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::appendNonbasicColsToBasisInterface
    /// (lines 1329-1404)
    pub(crate) fn append_nonbasic_cols_to_basis_interface(&mut self, ext_num_new_col: HighsInt) {
        if ext_num_new_col == 0 {
            return;
        }
        // HighsBasis& highs_basis = basis_;
        if !self.basis_.useful {
            return;
        }
        let has_simplex_basis = self.ekk_instance_.status_.has_basis;
        // SimplexBasis& simplex_basis = ekk_instance_.basis_;
        // HighsLp& lp = model_.lp_;

        debug_assert!(self.basis_.col_status.len() == self.model_lp_.num_col_ as usize);
        debug_assert!(self.basis_.row_status.len() == self.model_lp_.num_row_ as usize);

        // Add nonbasic structurals
        let new_num_col = self.model_lp_.num_col_ + ext_num_new_col;
        let new_num_tot = new_num_col + self.model_lp_.num_row_;
        self.basis_
            .col_status
            .resize(new_num_col as usize, HighsBasisStatus::default());
        if has_simplex_basis {
            self.ekk_instance_
                .basis_
                .nonbasic_flag_
                .resize(new_num_tot as usize, 0);
            self.ekk_instance_
                .basis_
                .nonbasic_move_
                .resize(new_num_tot as usize, 0);
            // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
            // necessary
            let mut i_row = self.model_lp_.num_row_ - 1;
            while i_row >= 0 {
                let i_col = self.ekk_instance_.basis_.basic_index_[i_row as usize];
                if i_col >= self.model_lp_.num_col_ {
                    // This basic variable is a row, so shift its index
                    self.ekk_instance_.basis_.basic_index_[i_row as usize] += ext_num_new_col;
                }
                self.ekk_instance_.basis_.nonbasic_flag_[(new_num_col + i_row) as usize] =
                    self.ekk_instance_.basis_.nonbasic_flag_
                        [(self.model_lp_.num_col_ + i_row) as usize];
                self.ekk_instance_.basis_.nonbasic_move_[(new_num_col + i_row) as usize] =
                    self.ekk_instance_.basis_.nonbasic_move_
                        [(self.model_lp_.num_col_ + i_row) as usize];
                i_row -= 1;
            }
        }
        // Make any new columns nonbasic
        let mut i_col = self.model_lp_.num_col_;
        while i_col < new_num_col {
            let lower = self.model_lp_.col_lower_[i_col as usize];
            let upper = self.model_lp_.col_upper_[i_col as usize];
            // Every branch below assigns both `status` and `move_val`.
            let status: HighsBasisStatus;
            let move_val: i8;
            if lower == upper {
                // Fixed
                status = HighsBasisStatus::kLower;
                move_val = kNonbasicMoveZe;
            } else if !highs_is_infinity(-lower) {
                // Finite lower bound so boxed or lower
                if !highs_is_infinity(upper) {
                    // Finite upper bound so boxed
                    if lower.abs() < upper.abs() {
                        status = HighsBasisStatus::kLower;
                        move_val = kNonbasicMoveUp;
                    } else {
                        status = HighsBasisStatus::kUpper;
                        move_val = kNonbasicMoveDn;
                    }
                } else {
                    // Lower (since upper bound is infinite)
                    status = HighsBasisStatus::kLower;
                    move_val = kNonbasicMoveUp;
                }
            } else if !highs_is_infinity(upper) {
                // Upper
                status = HighsBasisStatus::kUpper;
                move_val = kNonbasicMoveDn;
            } else {
                // FREE
                status = HighsBasisStatus::kZero;
                move_val = kNonbasicMoveZe;
            }
            debug_assert!(status != HighsBasisStatus::kNonbasic);
            self.basis_.col_status[i_col as usize] = status;
            if has_simplex_basis {
                debug_assert!(move_val != kIllegalMoveValue);
                self.ekk_instance_.basis_.nonbasic_flag_[i_col as usize] = kNonbasicFlagTrue;
                self.ekk_instance_.basis_.nonbasic_move_[i_col as usize] = move_val;
            }
            i_col += 1;
        }
    }

    /// HiGHS: Highs.cpp::Highs::deleteCols (interval overload, lines 3493-3507)
    pub fn delete_cols_interval(&mut self, from_col: HighsInt, to_col: HighsInt) -> HighsStatus {
        self.clear_derived_model_properties();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = create_interval(
            &mut index_collection,
            from_col,
            to_col,
            self.model_lp_.num_col_,
        );
        if create_error != 0 {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Interval [{}, {}] supplied to Highs::deleteCols is out of range [0, {})\n",
                from_col,
                to_col,
                self.model_lp_.num_col_
            );
            return HighsStatus::kError;
        }
        self.delete_cols_interface(&mut index_collection);
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: HighsInterface.cpp::Highs::deleteColsInterface (lines 711-767)
    ///
    /// For a mask collection, rewrites `index_collection.mask_` in place (as
    /// `delete_rows_interface`); hence `&mut`.
    pub(crate) fn delete_cols_interface(&mut self, index_collection: &mut HighsIndexCollection) {
        // HighsLp& lp = model_.lp_;
        // HighsBasis& basis = basis_;
        self.model_lp_.ensure_colwise();

        // Keep a copy of the original number of columns to check whether
        // any columns have been removed, and if there is mask to be updated
        let original_num_col = self.model_lp_.num_col_;

        self.model_lp_.delete_cols(index_collection);
        self.model_hessian_.delete_cols(index_collection);
        // Bail out if no columns were actually deleted
        if self.model_lp_.num_col_ == original_num_col {
            return;
        }

        debug_assert!(self.model_lp_.num_col_ < original_num_col);

        // Nontrivial deletion so reset the model_status and update any
        // Highs basis
        self.model_status_ = HighsModelStatus::kNotset;
        if self.basis_.useful {
            debug_assert!(self.basis_.col_status.len() == original_num_col as usize);
            // Have a full set of column basis status values, so maintain
            // them, and only invalidate the basis if a basic column has been
            // deleted
            delete_basis_cols(&mut self.basis_, index_collection, original_num_col);
        } else {
            debug_assert!(!self.basis_.valid);
        }

        if self.model_lp_.scale_.has_scaling {
            delete_scale(&mut self.model_lp_.scale_.col, index_collection);
            self.model_lp_
                .scale_
                .col
                .resize(self.model_lp_.num_col_ as usize, 0.0);
            self.model_lp_.scale_.num_col = self.model_lp_.num_col_;
        }
        // Deduce the consequences of deleting columns
        self.invalidate_model_status_solution_and_info();

        // Determine any implications for simplex data
        self.ekk_instance_.delete_cols(index_collection);

        if index_collection.is_mask_ {
            // Set the mask values to indicate the new index value of the
            // remaining columns
            let mut new_col: HighsInt = 0;
            for col in 0..original_num_col {
                if index_collection.mask_[col as usize] == 0 {
                    index_collection.mask_[col as usize] = new_col;
                    new_col += 1;
                } else {
                    index_collection.mask_[col as usize] = -1;
                }
            }
            debug_assert!(new_col == self.model_lp_.num_col_);
        }
        debug_assert!(lp_dimensions_ok(
            "deleteCols",
            &self.model_lp_,
            &self.options_.log_options
        ));
        // col_hash_.name2index.clear() dropped (no name hash in this port)
    }
}

/// HiGHS: HighsInterface.cpp::deleteBasisCols (static, lines 691-699)
fn delete_basis_cols(
    basis: &mut HighsBasis,
    index_collection: &HighsIndexCollection,
    original_num_col: HighsInt,
) {
    let mut deleted_basic = false;
    let mut deleted_nonbasic = false;
    delete_basis_entries(
        &mut basis.col_status,
        &mut deleted_basic,
        &mut deleted_nonbasic,
        index_collection,
        original_num_col,
    );
    if deleted_basic {
        basis.valid = false;
    }
}

impl HighsLp {
    /// HiGHS: HighsLp.cpp::HighsLp::addColNames (lines 305-315)
    pub fn add_col_names(&mut self, name: &str, num_new_col: HighsInt) {
        // Don't add names if there are no columns being added
        if self.num_col_ == 0 {
            return;
        }
        let col_names_size = self.col_names_.len() as HighsInt;
        if col_names_size <= 0 {
            return;
        }
        debug_assert!(col_names_size == self.num_col_ + num_new_col);
        // Handle the addition of user-defined names later
        debug_assert!(name.is_empty());
        // Blank names for the new columns were added in
        // appendColsToLpVectors
    }

    /// HiGHS: HighsLp.cpp::HighsLp::deleteColsFromVectors (lines 328-371)
    ///
    /// `col_names_` is compacted iff it is non-empty, as
    /// `delete_rows_from_vectors`.
    pub fn delete_cols_from_vectors(
        &mut self,
        new_num_col: &mut HighsInt,
        index_collection: &HighsIndexCollection,
    ) {
        debug_assert!(ok(index_collection));
        let mut from_k: HighsInt = 0;
        let mut to_k: HighsInt = 0;
        limits(index_collection, &mut from_k, &mut to_k);
        // Initialise new_num_col in case none is removed due to from_k > to_k
        *new_num_col = self.num_col_;
        if from_k > to_k {
            return;
        }

        let mut delete_from_col: HighsInt = 0;
        let mut delete_to_col: HighsInt = 0;
        let mut keep_from_col: HighsInt = 0;
        let mut keep_to_col: HighsInt = -1;
        let mut current_set_entry: HighsInt = 0;

        let col_dim = self.num_col_;
        *new_num_col = 0;
        let have_names = self.col_names_.len() != 0;
        let have_integrality = self.integrality_.len() != 0;
        let mut k = from_k;
        while k <= to_k {
            update_out_in_index(
                index_collection,
                &mut delete_from_col,
                &mut delete_to_col,
                &mut keep_from_col,
                &mut keep_to_col,
                &mut current_set_entry,
            );
            // Account for the initial columns being kept
            if k == from_k {
                *new_num_col = delete_from_col;
            }
            if delete_to_col >= col_dim - 1 {
                break;
            }
            debug_assert!(delete_to_col < col_dim);
            let mut col = keep_from_col;
            while col <= keep_to_col {
                self.col_cost_[*new_num_col as usize] = self.col_cost_[col as usize];
                self.col_lower_[*new_num_col as usize] = self.col_lower_[col as usize];
                self.col_upper_[*new_num_col as usize] = self.col_upper_[col as usize];
                if have_names {
                    self.col_names_[*new_num_col as usize] = self.col_names_[col as usize].clone();
                }
                if have_integrality {
                    self.integrality_[*new_num_col as usize] = self.integrality_[col as usize];
                }
                *new_num_col += 1;
                col += 1;
            }
            if keep_to_col >= col_dim - 1 {
                break;
            }
            k += 1;
        }
        self.col_cost_.resize(*new_num_col as usize, 0.0);
        self.col_lower_.resize(*new_num_col as usize, 0.0);
        self.col_upper_.resize(*new_num_col as usize, 0.0);
        if have_integrality {
            self.integrality_
                .resize(*new_num_col as usize, HighsVarType::kContinuous);
        }
        if have_names {
            self.col_names_.resize(*new_num_col as usize, String::new());
        }
    }

    /// HiGHS: HighsLp.cpp::HighsLp::deleteCols (lines 414-419)
    pub fn delete_cols(&mut self, index_collection: &HighsIndexCollection) {
        let mut new_num_col: HighsInt = 0;
        self.delete_cols_from_vectors(&mut new_num_col, index_collection);
        self.a_matrix_.delete_cols(index_collection);
        self.num_col_ = new_num_col;
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// Column-wise LP builder: `cols[j] = (cost, lower, upper, [(row, value)])`.
    #[allow(clippy::type_complexity)]
    fn lp(cols: &[(f64, f64, f64, &[(HighsInt, f64)])], rows: &[(f64, f64)]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = cols.len() as HighsInt;
        lp.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = vec![0];
        for &(cost, lower, upper, entries) in cols {
            lp.col_cost_.push(cost);
            lp.col_lower_.push(lower);
            lp.col_upper_.push(upper);
            for &(r, v) in entries {
                lp.a_matrix_.index_.push(r);
                lp.a_matrix_.value_.push(v);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        for &(l, u) in rows {
            lp.row_lower_.push(l);
            lp.row_upper_.push(u);
        }
        lp
    }

    fn quiet(lp: HighsLp) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs
    }

    fn iis_lp() -> HighsLp {
        lp(
            &[
                (0.0, 0.0, kHighsInf, &[(0, 1.0), (1, 1.0)]),
                (0.0, 0.0, kHighsInf, &[(0, 1.0), (1, -1.0)]),
            ],
            &[(-kHighsInf, 1.0), (2.0, kHighsInf)],
        )
    }

    #[test]
    fn add_then_delete_a_column() {
        let mut highs = quiet(iis_lp());
        assert_eq!(
            highs.add_cols(1, &[3.0], &[0.0], &[1.0], 1, &[0], &[1], &[5.0]),
            HighsStatus::kOk
        );
        let lp = &highs.model_lp_;
        assert_eq!(lp.num_col_, 3);
        assert_eq!(lp.a_matrix_.start_, vec![0, 2, 4, 5]);
        assert_eq!((lp.a_matrix_.index_[4], lp.a_matrix_.value_[4]), (1, 5.0));
        assert_eq!(lp.col_cost_[2], 3.0);
        assert_eq!(highs.delete_cols_interval(0, 0), HighsStatus::kOk);
        let lp = &highs.model_lp_;
        assert_eq!(lp.num_col_, 2);
        assert_eq!(lp.col_cost_, vec![0.0, 3.0]);
        assert_eq!(lp.a_matrix_.start_, vec![0, 2, 3]);
    }
}
