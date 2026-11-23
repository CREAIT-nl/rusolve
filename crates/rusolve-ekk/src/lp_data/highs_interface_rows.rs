//! HiGHS source: lp_data/HighsInterface.cpp (addRowsInterface,
//! deleteBasisEntries, deleteBasisRows, deleteRowsInterface,
//! appendBasicRowsToBasisInterface)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! API slice. Row addition / deletion on the
//! incumbent LP (`self.model_lp_`, C++ `model_.lp_`) with the HiGHS basis
//! kept in step. What survives (`kExtendInvertWhenAddingRows = false`,
//! HConst.h:45):
//! * `addRows`: `basis_` gets basic logicals appended (still valid);
//!   `ekk_instance_.add_rows` -> `update_status(kNewRows)` -> `HEkk::clear()`,
//!   so invert and DSE weights are gone; the next `run` refactorises from
//!   `basis_` via `HEkk::set_basis_from`.
//! * `deleteRows`: `basis_` is compacted; it stays valid iff only BASIC rows
//!   were deleted. `ekk_instance_.delete_rows` -> `kDelRows` -> `clear()`.
//!
//! Dropped: `isRowDataNull` / `isMatrixDataNull` (slices are never NULL),
//! the `kExtendInvertWhenAddingRows` debug branch (constant false),
//! `lp.addRowNames` (only asserts: the blank names for the new rows were
//! appended in `appendRowsToLpVectors`; kept as its `debug_assert!`),
//! `lp.row_hash_.name2index.clear()` (no row hash in this port).

use crate::lp_data::highs_run::Highs;
use crate::prelude::*;

/// HiGHS: HighsInterface.cpp::deleteBasisEntries (lines 648-689, 42 lines)
///
/// `static` in C++; `pub(crate)` here so that `deleteBasisCols`
/// (highs_api_cols.rs) can share it. The C++ caller passes `entry_dim` =
/// the ORIGINAL number of rows; `status.resize(new_num_entry)` is
/// `status.truncate(new_num_entry as usize)` (never grows).
pub(crate) fn delete_basis_entries(
    status: &mut Vec<HighsBasisStatus>,
    deleted_basic: &mut bool,
    deleted_nonbasic: &mut bool,
    index_collection: &HighsIndexCollection,
    entry_dim: HighsInt,
) {
    debug_assert!(ok(index_collection));
    debug_assert!(entry_dim as usize == status.len());
    let mut from_k: HighsInt = 0;
    let mut to_k: HighsInt = 0;
    limits(index_collection, &mut from_k, &mut to_k);
    if from_k > to_k {
        return;
    }

    let mut delete_from_entry: HighsInt = 0;
    let mut delete_to_entry: HighsInt = 0;
    let mut keep_from_entry: HighsInt = 0;
    let mut keep_to_entry: HighsInt = -1;
    let mut current_set_entry: HighsInt = 0;
    let mut new_num_entry: HighsInt = 0;
    *deleted_basic = false;
    *deleted_nonbasic = false;
    for k in from_k..=to_k {
        update_out_in_index(
            index_collection,
            &mut delete_from_entry,
            &mut delete_to_entry,
            &mut keep_from_entry,
            &mut keep_to_entry,
            &mut current_set_entry,
        );
        // Account for the initial entries being kept
        if k == from_k {
            new_num_entry = delete_from_entry;
        }
        // Identify whether a basic or a nonbasic entry has been deleted
        for entry in delete_from_entry..=delete_to_entry {
            if status[entry as usize] == HighsBasisStatus::kBasic {
                *deleted_basic = true;
            } else {
                *deleted_nonbasic = true;
            }
        }
        if delete_to_entry >= entry_dim - 1 {
            break;
        }
        for entry in keep_from_entry..=keep_to_entry {
            status[new_num_entry as usize] = status[entry as usize];
            new_num_entry += 1;
        }
        if keep_to_entry >= entry_dim - 1 {
            break;
        }
    }
    status.truncate(new_num_entry as usize);
}

/// HiGHS: HighsInterface.cpp::deleteBasisRows (lines 701-709, 9 lines)
///
/// Invalidates the basis iff a NONBASIC row was deleted (`deleted_nonbasic`);
/// `useful` is untouched.
fn delete_basis_rows(
    basis: &mut HighsBasis,
    index_collection: &HighsIndexCollection,
    original_num_row: HighsInt,
) {
    // C++ leaves these uninitialized on the (here unreachable) early return
    // of deleteBasisEntries; false is the only sound Rust initialisation.
    let mut deleted_basic: bool = false;
    let mut deleted_nonbasic: bool = false;
    delete_basis_entries(
        &mut basis.row_status,
        &mut deleted_basic,
        &mut deleted_nonbasic,
        index_collection,
        original_num_row,
    );
    if deleted_nonbasic {
        basis.valid = false;
    }
}

impl Highs {
    /// HiGHS: HighsInterface.cpp::Highs::addRowsInterface (lines 526-646, 121 lines)
    /// HiGHS 1.16 (#3386): HighsInterface.cpp::addRowsInterface (lines 528-650)
    ///
    /// * `assessBounds(options, "Row", lp.num_row_, index_collection, ...)`
    ///   is `crate::lp_data::lp_assess::assess_bounds(.., None)` on local
    ///   copies of the bounds; `appendRowsToLpVectors` is
    ///   `crate::lp_data::lp_utils::append_rows_to_lp_vectors`.
    /// * `local_ar_matrix` is a ROW-wise `HighsSparseMatrix`
    ///   (`format_ = MatrixFormat::kRowwise`, `num_col_ = lp.num_col_`,
    ///   `num_row_ = ext_num_new_row`); `local_ar_matrix.assess(...)` is
    ///   util/sparse_matrix/assess.rs:14; `lp.a_matrix_.addRows(local_ar_matrix)`
    ///   is `add_rows(&local_ar_matrix, None)` (modify.rs:146).
    /// * If `lp.scale_.has_scaling`: extend `scale.row` with 1.0, then
    ///   `local_ar_matrix.apply_col_scale(&lp.scale_)` and
    ///   `consider_row_scaling(options.allowed_matrix_scale_factor,
    ///   &mut scale.row[old_num_row..])`. The scaled `local_ar_matrix` is what
    ///   `ekk_instance_.add_rows(&lp, &local_ar_matrix)` receives.
    /// * Order is load-bearing: `appendBasicRowsToBasisInterface` runs BEFORE
    ///   `lp.num_row_ += ext_num_new_row` (it reads the old `num_row_`).
    /// * Ends with `invalidate_model_status_solution_and_info()`
    ///   (highs_api_support.rs) and `self.ekk_instance_.add_rows(..)`
    ///   (simplex/ekk/solve.rs:337).
    ///
    /// `ext_ar_start` has `ext_num_new_row` entries (no end sentinel), as in C++.
    pub(crate) fn add_rows_interface(
        &mut self,
        ext_num_new_row: HighsInt,
        ext_row_lower: &[f64],
        ext_row_upper: &[f64],
        ext_num_new_nz: HighsInt,
        ext_ar_start: &[HighsInt],
        ext_ar_index: &[HighsInt],
        ext_ar_value: &[f64],
    ) -> HighsStatus {
        // addRows is fundamentally different from addCols, since the new
        // matrix data are held row-wise, so we have to insert data into the
        // column-wise matrix of the LP.
        // The kExtendInvertWhenAddingRows debug branch is dropped: the
        // constant is false (HConst.h:45).
        let mut return_status = HighsStatus::kOk;
        if ext_num_new_row < 0 {
            return HighsStatus::kError;
        }
        if ext_num_new_nz < 0 {
            return HighsStatus::kError;
        }
        if ext_num_new_row == 0 {
            return HighsStatus::kOk;
        }
        // isRowDataNull / isMatrixDataNull dropped: the slices are never NULL.

        // bool& useful_basis = basis_.useful; bool& lp_has_scaling =
        // lp.scale_.has_scaling; nothing between the binding and the use
        // changes either, so they are read as values here.
        let useful_basis = self.basis_.useful;
        let lp_has_scaling = self.model_lp_.scale_.has_scaling;

        // Check that if nonzeros are to be added then the model has a positive
        // number of columns
        if self.model_lp_.num_col_ <= 0 && ext_num_new_nz > 0 {
            return HighsStatus::kError;
        }

        // Record the new number of rows
        let new_num_row = self.model_lp_.num_row_ + ext_num_new_row;

        let index_collection = HighsIndexCollection {
            dimension_: ext_num_new_row,
            is_interval_: true,
            from_: 0,
            to_: ext_num_new_row - 1,
            ..HighsIndexCollection::default()
        };
        // Take a copy of the bounds that can be normalised
        let mut local_row_lower: Vec<f64> = ext_row_lower[..ext_num_new_row as usize].to_vec();
        let mut local_row_upper: Vec<f64> = ext_row_upper[..ext_num_new_row as usize].to_vec();

        return_status = interpret_call_status(
            &self.options_.log_options,
            crate::lp_data::lp_assess::assess_bounds(
                &self.options_,
                "Row",
                self.model_lp_.num_row_,
                &index_collection,
                &mut local_row_lower,
                &mut local_row_upper,
                self.options_.infinite_bound,
                None,
            ),
            return_status,
            "assessBounds",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        // HiGHS 1.16 #3386 (h16_add_no_modify_on_error): HighsInterface.cpp lines 584-622 (1.16.0)
        let h16_add_no_modify_on_error = self.options_.h16_add_no_modify_on_error;
        // Append the rows to the LP vectors
        if !h16_add_no_modify_on_error {
            append_rows_to_lp_vectors(
                &mut self.model_lp_,
                ext_num_new_row,
                &local_row_lower,
                &local_row_upper,
            );
        }

        // Form a row-wise HighsSparseMatrix of the new matrix rows so that
        // is easy to handle and, if there are nonzeros, it can be
        // normalised
        let mut local_ar_matrix = HighsSparseMatrix::default();
        local_ar_matrix.num_col_ = self.model_lp_.num_col_;
        local_ar_matrix.num_row_ = ext_num_new_row;
        local_ar_matrix.format_ = MatrixFormat::kRowwise;
        if ext_num_new_nz != 0 {
            // {ext_ar_start, ext_ar_start + ext_num_new_row}: exact lengths as C++.
            local_ar_matrix.start_ = ext_ar_start[..ext_num_new_row as usize].to_vec();
            local_ar_matrix
                .start_
                .resize((ext_num_new_row + 1) as usize, 0);
            local_ar_matrix.start_[ext_num_new_row as usize] = ext_num_new_nz;
            local_ar_matrix.index_ = ext_ar_index[..ext_num_new_nz as usize].to_vec();
            local_ar_matrix.value_ = ext_ar_value[..ext_num_new_nz as usize].to_vec();
            // Assess the matrix columns
            return_status = interpret_call_status(
                &self.options_.log_options,
                local_ar_matrix.assess(
                    &self.options_.log_options,
                    "LP",
                    self.options_.small_matrix_value,
                    self.options_.large_matrix_value,
                    self.options_.h16_assess_small,
                ),
                return_status,
                "assessMatrix",
            );
            if return_status == HighsStatus::kError {
                return return_status;
            }
        } else {
            // No nonzeros so, whether the constraint matrix is row-wise or
            // column-wise, adding the empty matrix is trivial. Complete the
            // setup of an empty row-wise HighsSparseMatrix of the new matrix
            // rows
            local_ar_matrix.start_ = vec![0; (ext_num_new_row + 1) as usize];
        }
        // HiGHS 1.16 #3386 (h16_add_no_modify_on_error): HighsInterface.cpp lines 617-622 (1.16.0)
        if h16_add_no_modify_on_error {
            debug_assert!(return_status != HighsStatus::kError);
            // Now that data for new rows has been assessed and no errors can
            // occur, append the rows to the LP vectors and matrix
            append_rows_to_lp_vectors(
                &mut self.model_lp_,
                ext_num_new_row,
                &local_row_lower,
                &local_row_upper,
            );
        }
        // Append the rows to LP matrix
        self.model_lp_.a_matrix_.add_rows(&local_ar_matrix, None);
        if lp_has_scaling {
            // lp.num_row_ is still the OLD number of rows here (C++ reads it
            // through the `lp` reference before lp.num_row_ is increased).
            let old_num_row = self.model_lp_.num_row_;
            // Extend the row scaling factors
            self.model_lp_.scale_.row.resize(new_num_row as usize, 0.0);
            for i_row in 0..ext_num_new_row {
                self.model_lp_.scale_.row[old_num_row as usize + i_row as usize] = 1.0;
            }
            self.model_lp_.scale_.num_row = new_num_row;
            // Apply the existing column scaling to the new rows
            local_ar_matrix.apply_col_scale(&self.model_lp_.scale_);
            // Consider applying row scaling to the new rows.
            local_ar_matrix.consider_row_scaling(
                self.options_.allowed_matrix_scale_factor,
                &mut self.model_lp_.scale_.row[old_num_row as usize..],
            );
        }
        // Update the basis corresponding to new basic rows
        if useful_basis {
            self.append_basic_rows_to_basis_interface(ext_num_new_row);
        }

        // lp.addRowNames: blank names for the new rows were added in
        // appendRowsToLpVectors; only its assert remains. The MIP LP relaxation
        // keeps the model's row names (HighsLpRelaxation::loadModel copies the LP).
        debug_assert!(
            self.model_lp_.num_row_ == 0
                || self.model_lp_.row_names_.is_empty()
                || self.model_lp_.row_names_.len()
                    == (self.model_lp_.num_row_ + ext_num_new_row) as usize
        );

        // Increase the number of rows in the LP
        self.model_lp_.num_row_ += ext_num_new_row;
        debug_assert!(crate::lp_data::lp_assess::lp_dimensions_ok(
            "addRows",
            &self.model_lp_,
            &self.options_.log_options
        ));

        // Deduce the consequences of adding new rows
        self.invalidate_model_status_solution_and_info();
        // Determine any implications for simplex data
        self.ekk_instance_
            .add_rows(&self.model_lp_, &local_ar_matrix);

        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::deleteRowsInterface (lines 769-820, 52 lines)
    ///
    /// * `lp.deleteRows` is `HighsLp::delete_rows` (highs_lp_modify.rs);
    ///   `deleteScale` is `highs_lp_modify::delete_scale` followed by
    ///   `scale_.row.resize(num_row)`, `scale_.num_row = num_row`.
    /// * Bails out (no invalidation at all) when no row was deleted.
    /// * For a mask collection, rewrites `index_collection.mask_` in place to
    ///   the new row index or -1 (C++ `Highs::deleteRows(HighsInt* mask)`
    ///   copies it back to the caller); hence `&mut`.
    pub(crate) fn delete_rows_interface(&mut self, index_collection: &mut HighsIndexCollection) {
        self.model_lp_.ensure_colwise();
        // Keep a copy of the original number of rows to check whether
        // any rows have been removed, and if there is mask to be updated
        let original_num_row = self.model_lp_.num_row_;

        self.model_lp_.delete_rows(index_collection);
        // Bail out if no rows were actually deleted
        if self.model_lp_.num_row_ == original_num_row {
            return;
        }

        debug_assert!(self.model_lp_.num_row_ < original_num_row);

        // Nontrivial deletion so reset the model_status and update any
        // Highs basis
        self.model_status_ = HighsModelStatus::kNotset;
        if self.basis_.useful {
            debug_assert!(self.basis_.row_status.len() == original_num_row as usize);
            // Have a full set of row basis status values, so maintain them,
            // and only invalidate the basis if a nonbasic row has been
            // deleted
            delete_basis_rows(&mut self.basis_, index_collection, original_num_row);
        } else {
            debug_assert!(!self.basis_.valid);
        }

        if self.model_lp_.scale_.has_scaling {
            let num_row = self.model_lp_.num_row_;
            delete_scale(&mut self.model_lp_.scale_.row, index_collection);
            self.model_lp_.scale_.row.resize(num_row as usize, 0.0);
            self.model_lp_.scale_.num_row = num_row;
        }
        // Deduce the consequences of deleting rows
        self.invalidate_model_status_solution_and_info();

        // Determine any implications for simplex data
        self.ekk_instance_.delete_rows(index_collection);
        if index_collection.is_mask_ {
            let mut new_row: HighsInt = 0;
            for row in 0..original_num_row {
                if index_collection.mask_[row as usize] == 0 {
                    index_collection.mask_[row as usize] = new_row;
                    new_row += 1;
                } else {
                    index_collection.mask_[row as usize] = -1;
                }
            }
            debug_assert!(new_row == self.model_lp_.num_row_);
        }
        debug_assert!(crate::lp_data::lp_assess::lp_dimensions_ok(
            "deleteRows",
            &self.model_lp_,
            &self.options_.log_options
        ));
        // lp.row_hash_.name2index.clear() dropped: no row hash in this port.
    }

    /// HiGHS: HighsInterface.cpp::Highs::appendBasicRowsToBasisInterface (lines 1406-1436, 31 lines)
    ///
    /// Extends `basis_.row_status` with `kBasic` and, iff
    /// `ekk_instance_.status_.has_basis`, the simplex basis
    /// (`nonbasic_flag_`/`nonbasic_move_`/`basic_index_`) with the new
    /// logicals. Reads the PRE-addition `model_lp_.num_row_`.
    pub(crate) fn append_basic_rows_to_basis_interface(&mut self, ext_num_new_row: HighsInt) {
        if ext_num_new_row == 0 {
            return;
        }
        let highs_basis = &mut self.basis_;
        if !highs_basis.useful {
            return;
        }
        let has_simplex_basis = self.ekk_instance_.status_.has_basis;
        let simplex_basis = &mut self.ekk_instance_.basis_;
        let lp = &self.model_lp_;

        debug_assert!(highs_basis.col_status.len() == lp.num_col_ as usize);
        debug_assert!(highs_basis.row_status.len() == lp.num_row_ as usize);

        // Add basic logicals
        // Add the new rows to the Highs basis
        let new_num_row = lp.num_row_ + ext_num_new_row;
        highs_basis
            .row_status
            .resize(new_num_row as usize, HighsBasisStatus::default());
        for i_row in lp.num_row_..new_num_row {
            highs_basis.row_status[i_row as usize] = HighsBasisStatus::kBasic;
        }
        if has_simplex_basis {
            // Add the new rows to the simplex basis
            let new_num_tot = lp.num_col_ + new_num_row;
            simplex_basis.nonbasic_flag_.resize(new_num_tot as usize, 0);
            simplex_basis.nonbasic_move_.resize(new_num_tot as usize, 0);
            simplex_basis.basic_index_.resize(new_num_row as usize, 0);
            for i_row in lp.num_row_..new_num_row {
                simplex_basis.nonbasic_flag_[lp.num_col_ as usize + i_row as usize] =
                    kNonbasicFlagFalse;
                simplex_basis.nonbasic_move_[lp.num_col_ as usize + i_row as usize] = 0;
                simplex_basis.basic_index_[i_row as usize] = lp.num_col_ + i_row;
            }
        }
    }
}
