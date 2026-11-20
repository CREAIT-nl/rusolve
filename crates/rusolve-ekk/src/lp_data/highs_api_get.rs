//! HiGHS source: lp_data/Highs.cpp (changeColBounds, changeRowBounds,
//! changeRowsBounds interval/set, changeColsCost interval, getCols interval,
//! getRows interval), lp_data/HighsInterface.cpp (changeRowBoundsInterface,
//! getColsInterface, getRowsInterface), lp_data/HighsLpUtils.cpp
//! (changeLpRowBounds, getSubVectors, getSubVectorsTranspose)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The bound/cost change and column/row query API that
//! `HighsIis` and `Highs::elasticityFilter` use. Mirrors the column
//! versions in `highs_api.rs` / `highs_interface_bounds.rs`. Dropped as
//! there: `doubleUserDataNotNull` (slices are never NULL). Nullable output
//! arrays of `getCols` / `getRows` / `getSubVectors*` are `Option<&mut [T]>`.

use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::changeColBounds (lines 3106-3109)
    pub fn change_col_bounds(&mut self, col: HighsInt, lower: f64, upper: f64) -> HighsStatus {
        self.change_cols_bounds_set(1, &[col], &[lower], &[upper])
    }

    /// HiGHS: Highs.cpp::Highs::changeRowBounds (lines 3187-3190)
    pub fn change_row_bounds(&mut self, row: HighsInt, lower: f64, upper: f64) -> HighsStatus {
        self.change_rows_bounds_set(1, &[row], &[lower], &[upper])
    }

    /// HiGHS: Highs.cpp::Highs::changeRowsBounds (interval overload, lines 3192-3213)
    pub fn change_rows_bounds_interval(
        &mut self,
        from_row: HighsInt,
        to_row: HighsInt,
        lower: &[f64],
        upper: &[f64],
    ) -> HighsStatus {
        self.clear_derived_model_properties();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = create_interval(
            &mut index_collection,
            from_row,
            to_row,
            self.model_lp_.num_row_,
        );
        if create_error != 0 {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Interval [{}, {}] supplied to Highs::changeRowsBounds is out of range [0, {})\n",
                from_row,
                to_row,
                self.model_lp_.num_row_
            );
            return HighsStatus::kError;
        }
        let call_status = self.change_row_bounds_interface(&mut index_collection, lower, upper);
        let mut return_status = HighsStatus::kOk;
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "changeRowBounds",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::changeRowsBounds (set overload, lines 3215-3250)
    pub fn change_rows_bounds_set(
        &mut self,
        num_set_entries: HighsInt,
        set: &[HighsInt],
        lower: &[f64],
        upper: &[f64],
    ) -> HighsStatus {
        if num_set_entries == 0 {
            return HighsStatus::kOk;
        }
        // Check for NULL data in "set" version of changeRowsBounds since
        // values are sorted with set
        // PORT: doubleUserDataNotNull dropped (slices are never NULL)
        self.clear_derived_model_properties();
        // Ensure that the set and data are in ascending order
        let mut local_lower: Vec<f64> = lower[..num_set_entries as usize].to_vec();
        let mut local_upper: Vec<f64> = upper[..num_set_entries as usize].to_vec();
        let mut local_set: Vec<HighsInt> = set[..num_set_entries as usize].to_vec();
        sort_set_data(
            num_set_entries,
            &mut local_set,
            Some(lower),
            Some(upper),
            None,
            Some(&mut local_lower),
            Some(&mut local_upper),
            None,
        );
        let mut index_collection = HighsIndexCollection::default();
        let create_error = create_set(
            &mut index_collection,
            num_set_entries,
            &local_set,
            self.model_lp_.num_row_,
        );
        if create_error != 0 {
            // analyseSetCreateError is not ported: one log line instead.
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Set supplied to Highs::changeRowsBounds is illegal (create error {})\n",
                create_error
            );
            return HighsStatus::kError;
        }
        let call_status =
            self.change_row_bounds_interface(&mut index_collection, &local_lower, &local_upper);
        let mut return_status = HighsStatus::kOk;
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "changeRowBounds",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: HighsInterface.cpp::Highs::changeRowBoundsInterface (lines 1009-1058)
    /// HiGHS 1.16 (#3386): HighsInterface.cpp::changeRowBoundsInterface (lines 1012-1061)
    pub(crate) fn change_row_bounds_interface(
        &mut self,
        index_collection: &mut HighsIndexCollection,
        lower: &[f64],
        upper: &[f64],
    ) -> HighsStatus {
        let num_row_bounds = data_size(index_collection);
        // If a non-positive number of costs (may) need changing nothing needs to be
        // done
        if num_row_bounds <= 0 {
            return HighsStatus::kOk;
        }
        // PORT: doubleUserDataNotNull dropped (slices are never NULL)
        // Take a copy of the cost that can be normalised
        let mut local_row_lower: Vec<f64> = lower[..num_row_bounds as usize].to_vec();
        let mut local_row_upper: Vec<f64> = upper[..num_row_bounds as usize].to_vec();
        // If changing the bounds for a set of rows, ensure that the
        // set and data are in ascending order
        if index_collection.is_set_ {
            sort_set_data(
                index_collection.set_num_entries_,
                &mut index_collection.set_,
                Some(lower),
                Some(upper),
                None,
                Some(&mut local_row_lower),
                Some(&mut local_row_upper),
                None,
            );
        }
        let mut return_status = HighsStatus::kOk;
        return_status = interpret_call_status(
            &self.options_.log_options,
            crate::lp_data::lp_assess::assess_bounds(
                &self.options_,
                "row",
                0,
                index_collection,
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

        change_lp_row_bounds(
            &mut self.model_lp_,
            index_collection,
            &local_row_lower,
            &local_row_upper,
        );
        // Update HiGHS basis status and (any) simplex move status of
        // nonbasic variables whose bounds have changed
        self.set_nonbasic_status_interface(index_collection, false);
        // Deduce the consequences of new row bounds
        if !self.basis_.useful && self.feasible_wrt_bounds(false) {
            // Retain the solution if there's no basis, and the solution is
            // feasible
            self.invalidate_model_status_and_info();
        } else {
            // Invalidate the solution
            self.invalidate_model_status_solution_and_info();
        }
        // Determine any implications for simplex data
        self.ekk_instance_.update_status(LpAction::kNewBounds);
        // HiGHS 1.16 #3386 (h16_add_no_modify_on_error): HighsInterface.cpp line 1060 (1.16.0)
        if self.options_.h16_add_no_modify_on_error {
            return_status
        } else {
            HighsStatus::kOk
        }
    }

    /// HiGHS: Highs.cpp::Highs::changeColsCost (interval overload, lines 3042-3061)
    pub fn change_cols_cost_interval(
        &mut self,
        from_col: HighsInt,
        to_col: HighsInt,
        cost: &[f64],
    ) -> HighsStatus {
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
                "Interval [{}, {}] supplied to Highs::changeColsCost is out of range [0, {})\n",
                from_col,
                to_col,
                self.model_lp_.num_col_
            );
            return HighsStatus::kError;
        }
        let call_status = self.change_costs_interface(&index_collection, cost);
        let mut return_status = HighsStatus::kOk;
        return_status = interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "changeCosts",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::getCols (interval overload, lines 3307-3330)
    pub fn get_cols_interval(
        &self,
        from_col: HighsInt,
        to_col: HighsInt,
        num_col: &mut HighsInt,
        costs: Option<&mut [f64]>,
        lower: Option<&mut [f64]>,
        upper: Option<&mut [f64]>,
        num_nz: &mut HighsInt,
        start: Option<&mut [HighsInt]>,
        index: Option<&mut [HighsInt]>,
        value: Option<&mut [f64]>,
    ) -> HighsStatus {
        if from_col > to_col {
            // Empty interval
            *num_col = 0;
            *num_nz = 0;
            return HighsStatus::kOk;
        }
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
                "Interval [{}, {}] supplied to Highs::getCols is out of range [0, {})\n",
                from_col,
                to_col,
                self.model_lp_.num_col_
            );
            return HighsStatus::kError;
        }
        self.get_cols_interface(
            &index_collection,
            num_col,
            costs,
            lower,
            upper,
            num_nz,
            start,
            index,
            value,
        );
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::getColsInterface (lines 822-837)
    pub(crate) fn get_cols_interface(
        &self,
        index_collection: &HighsIndexCollection,
        num_col: &mut HighsInt,
        cost: Option<&mut [f64]>,
        lower: Option<&mut [f64]>,
        upper: Option<&mut [f64]>,
        num_nz: &mut HighsInt,
        start: Option<&mut [HighsInt]>,
        index: Option<&mut [HighsInt]>,
        value: Option<&mut [f64]>,
    ) {
        let lp = &self.model_lp_;
        if lp.a_matrix_.is_colwise() {
            get_sub_vectors(
                index_collection,
                lp.num_col_,
                Some(&lp.col_cost_),
                Some(&lp.col_lower_),
                Some(&lp.col_upper_),
                &lp.a_matrix_,
                num_col,
                cost,
                lower,
                upper,
                num_nz,
                start,
                index,
                value,
            );
        } else {
            get_sub_vectors_transpose(
                index_collection,
                lp.num_col_,
                Some(&lp.col_cost_),
                Some(&lp.col_lower_),
                Some(&lp.col_upper_),
                &lp.a_matrix_,
                num_col,
                cost,
                lower,
                upper,
                num_nz,
                start,
                index,
                value,
            );
        }
    }

    /// HiGHS: Highs.cpp::Highs::getRows (interval overload, lines 3396-3419)
    pub fn get_rows_interval(
        &self,
        from_row: HighsInt,
        to_row: HighsInt,
        num_row: &mut HighsInt,
        lower: Option<&mut [f64]>,
        upper: Option<&mut [f64]>,
        num_nz: &mut HighsInt,
        start: Option<&mut [HighsInt]>,
        index: Option<&mut [HighsInt]>,
        value: Option<&mut [f64]>,
    ) -> HighsStatus {
        if from_row > to_row {
            // Empty interval
            *num_row = 0;
            *num_nz = 0;
            return HighsStatus::kOk;
        }
        let mut index_collection = HighsIndexCollection::default();
        let create_error = create_interval(
            &mut index_collection,
            from_row,
            to_row,
            self.model_lp_.num_row_,
        );
        if create_error != 0 {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Interval [{}, {}] supplied to Highs::getRows is out of range [0, {})\n",
                from_row,
                to_row,
                self.model_lp_.num_row_
            );
            return HighsStatus::kError;
        }
        self.get_rows_interface(
            &index_collection,
            num_row,
            lower,
            upper,
            num_nz,
            start,
            index,
            value,
        );
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::getRowsInterface (lines 839-854)
    pub(crate) fn get_rows_interface(
        &self,
        index_collection: &HighsIndexCollection,
        num_row: &mut HighsInt,
        lower: Option<&mut [f64]>,
        upper: Option<&mut [f64]>,
        num_nz: &mut HighsInt,
        start: Option<&mut [HighsInt]>,
        index: Option<&mut [HighsInt]>,
        value: Option<&mut [f64]>,
    ) {
        let lp = &self.model_lp_;
        if lp.a_matrix_.is_colwise() {
            get_sub_vectors_transpose(
                index_collection,
                lp.num_row_,
                None,
                Some(&lp.row_lower_),
                Some(&lp.row_upper_),
                &lp.a_matrix_,
                num_row,
                None,
                lower,
                upper,
                num_nz,
                start,
                index,
                value,
            );
        } else {
            get_sub_vectors(
                index_collection,
                lp.num_row_,
                None,
                Some(&lp.row_lower_),
                Some(&lp.row_upper_),
                &lp.a_matrix_,
                num_row,
                None,
                lower,
                upper,
                num_nz,
                start,
                index,
                value,
            );
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::changeLpRowBounds (lines 1809-1815)
pub fn change_lp_row_bounds(
    lp: &mut HighsLp,
    index_collection: &HighsIndexCollection,
    new_row_lower: &[f64],
    new_row_upper: &[f64],
) {
    change_bounds(
        &mut lp.row_lower_,
        &mut lp.row_upper_,
        index_collection,
        new_row_lower,
        new_row_upper,
    );
}

/// HiGHS: HighsLpUtils.cpp::getSubVectors (lines 3471-3525)
///
/// PORT NOTE: nullable `data0..2` inputs are `Option<&[f64]>`.
pub fn get_sub_vectors(
    index_collection: &HighsIndexCollection,
    data_dim: HighsInt,
    data0: Option<&[f64]>,
    data1: Option<&[f64]>,
    data2: Option<&[f64]>,
    matrix: &HighsSparseMatrix,
    num_sub_vector: &mut HighsInt,
    mut sub_vector_data0: Option<&mut [f64]>,
    mut sub_vector_data1: Option<&mut [f64]>,
    mut sub_vector_data2: Option<&mut [f64]>,
    sub_matrix_num_nz: &mut HighsInt,
    mut sub_matrix_start: Option<&mut [HighsInt]>,
    mut sub_matrix_index: Option<&mut [HighsInt]>,
    mut sub_matrix_value: Option<&mut [f64]>,
) {
    // Ensure that if there's no data0 then it's not required in the
    // sub-vector
    if data0.is_none() {
        debug_assert!(sub_vector_data0.is_none());
    }
    debug_assert!(ok(index_collection));
    let mut from_k: HighsInt = 0;
    let mut to_k: HighsInt = 0;
    limits(index_collection, &mut from_k, &mut to_k);
    // Surely this is checked elsewhere
    debug_assert!(0 <= from_k && to_k < data_dim);
    debug_assert!(from_k <= to_k);
    let mut out_from_vector: HighsInt = 0;
    let mut out_to_vector: HighsInt = 0;
    let mut in_from_vector: HighsInt = 0;
    let mut in_to_vector: HighsInt = -1;
    let mut current_set_entry: HighsInt = 0;
    *num_sub_vector = 0;
    *sub_matrix_num_nz = 0;
    let mut k: HighsInt = from_k;
    while k <= to_k {
        update_out_in_index(
            index_collection,
            &mut out_from_vector,
            &mut out_to_vector,
            &mut in_from_vector,
            &mut in_to_vector,
            &mut current_set_entry,
        );
        debug_assert!(out_to_vector < data_dim);
        debug_assert!(in_to_vector < data_dim);
        let mut i_vector: HighsInt = out_from_vector;
        while i_vector <= out_to_vector {
            if let Some(ref mut svd0) = sub_vector_data0 {
                svd0[*num_sub_vector as usize] = data0.unwrap()[i_vector as usize];
            }
            if let Some(ref mut svd1) = sub_vector_data1 {
                svd1[*num_sub_vector as usize] = data1.unwrap()[i_vector as usize];
            }
            if let Some(ref mut svd2) = sub_vector_data2 {
                svd2[*num_sub_vector as usize] = data2.unwrap()[i_vector as usize];
            }
            if let Some(ref mut sms) = sub_matrix_start {
                sms[*num_sub_vector as usize] = *sub_matrix_num_nz
                    + matrix.start_[i_vector as usize]
                    - matrix.start_[out_from_vector as usize];
            }
            *num_sub_vector += 1;
            i_vector += 1;
        }
        let mut i_el: HighsInt = matrix.start_[out_from_vector as usize];
        while i_el < matrix.start_[(out_to_vector + 1) as usize] {
            if let Some(ref mut smi) = sub_matrix_index {
                smi[*sub_matrix_num_nz as usize] = matrix.index_[i_el as usize];
            }
            if let Some(ref mut smv) = sub_matrix_value {
                smv[*sub_matrix_num_nz as usize] = matrix.value_[i_el as usize];
            }
            *sub_matrix_num_nz += 1;
            i_el += 1;
        }
        if out_to_vector == data_dim - 1 || in_to_vector == data_dim - 1 {
            break;
        }
        k += 1;
    }
}

/// HiGHS: HighsLpUtils.cpp::getSubVectorsTranspose (lines 3527-3662)
pub fn get_sub_vectors_transpose(
    index_collection: &HighsIndexCollection,
    data_dim: HighsInt,
    data0: Option<&[f64]>,
    data1: Option<&[f64]>,
    data2: Option<&[f64]>,
    matrix: &HighsSparseMatrix,
    num_sub_vector: &mut HighsInt,
    mut sub_vector_data0: Option<&mut [f64]>,
    mut sub_vector_data1: Option<&mut [f64]>,
    mut sub_vector_data2: Option<&mut [f64]>,
    sub_matrix_num_nz: &mut HighsInt,
    mut sub_matrix_start: Option<&mut [HighsInt]>,
    mut sub_matrix_index: Option<&mut [HighsInt]>,
    mut sub_matrix_value: Option<&mut [f64]>,
) {
    // Ensure that if there's no data0 then it's not required in the
    // sub-vector
    if data0.is_none() {
        debug_assert!(sub_vector_data0.is_none());
    }
    debug_assert!(ok(index_collection));
    let mut from_k: HighsInt = 0;
    let mut to_k: HighsInt = 0;
    limits(index_collection, &mut from_k, &mut to_k);
    // Surely this is checked elsewhere
    debug_assert!(0 <= from_k && to_k < data_dim);
    debug_assert!(from_k <= to_k);
    // "Out" means not in the set to be extracted
    // "In" means in the set to be extracted
    let mut out_from_vector: HighsInt = 0;
    // PORT: C++ sets `out_to_vector = -1; current_set_entry = 0;` at the
    // top of the non-mask branch (the only branch that uses them); the
    // values are folded into the initialisers here.
    let mut out_to_vector: HighsInt = -1;
    let mut in_from_vector: HighsInt = 0;
    let mut in_to_vector: HighsInt = -1;
    let mut current_set_entry: HighsInt = 0;
    // Set up a mask so that entries to be got from the matrix can be
    // identified and have their correct index.
    let mut new_index: Vec<HighsInt> = vec![0; data_dim as usize];

    *num_sub_vector = 0;
    *sub_matrix_num_nz = 0;
    if !index_collection.is_mask_ {
        let mut k: HighsInt = from_k;
        while k <= to_k {
            update_out_in_index(
                index_collection,
                &mut in_from_vector,
                &mut in_to_vector,
                &mut out_from_vector,
                &mut out_to_vector,
                &mut current_set_entry,
            );
            if k == from_k {
                // Account for any initial vectors not being extracted
                let mut i_vector: HighsInt = 0;
                while i_vector < in_from_vector {
                    new_index[i_vector as usize] = -1;
                    i_vector += 1;
                }
            }
            let mut i_vector: HighsInt = in_from_vector;
            while i_vector <= in_to_vector {
                new_index[i_vector as usize] = *num_sub_vector;
                *num_sub_vector += 1;
                i_vector += 1;
            }
            let mut i_vector: HighsInt = out_from_vector;
            while i_vector <= out_to_vector {
                new_index[i_vector as usize] = -1;
                i_vector += 1;
            }
            if out_to_vector >= data_dim - 1 {
                break;
            }
            k += 1;
        }
    } else {
        let mut i_vector: HighsInt = 0;
        while i_vector < data_dim {
            if index_collection.mask_[i_vector as usize] != 0 {
                new_index[i_vector as usize] = *num_sub_vector;
                *num_sub_vector += 1;
            } else {
                new_index[i_vector as usize] = -1;
            }
            i_vector += 1;
        }
    }

    // Bail out if no vectors are to be extracted
    if *num_sub_vector == 0 {
        return;
    }

    let mut i_vector: HighsInt = 0;
    while i_vector < data_dim {
        let new_i_vector = new_index[i_vector as usize];
        if new_i_vector >= 0 {
            debug_assert!(new_i_vector < *num_sub_vector);
            if let Some(ref mut svd0) = sub_vector_data0 {
                svd0[new_i_vector as usize] = data0.unwrap()[i_vector as usize];
            }
            if let Some(ref mut svd1) = sub_vector_data1 {
                svd1[new_i_vector as usize] = data1.unwrap()[i_vector as usize];
            }
            if let Some(ref mut svd2) = sub_vector_data2 {
                svd2[new_i_vector as usize] = data2.unwrap()[i_vector as usize];
            }
        }
        i_vector += 1;
    }
    let extract_start = sub_matrix_start.is_some();
    let extract_index = sub_matrix_index.is_some();
    let extract_value = sub_matrix_value.is_some();
    let extract_matrix = extract_index || extract_value;
    // Allocate an array of lengths for the sub-matrix to be
    // extracted: necessary even if just the number of nonzeros is
    // required
    let mut sub_matrix_length: Vec<HighsInt> = vec![0; *num_sub_vector as usize];
    // Identify the lengths of the vectors in the sub-matrix to be extracted
    let num_vector: HighsInt = matrix.start_.len() as HighsInt - 1;
    let mut vector: HighsInt = 0;
    while vector < num_vector {
        let mut i_el: HighsInt = matrix.start_[vector as usize];
        while i_el < matrix.start_[(vector + 1) as usize] {
            let i_vector = matrix.index_[i_el as usize];
            let new_i_vector = new_index[i_vector as usize];
            if new_i_vector >= 0 {
                sub_matrix_length[new_i_vector as usize] += 1;
            }
            i_el += 1;
        }
        vector += 1;
    }
    if !extract_start {
        // bail out if no matrix starts are to be extracted, but only after
        // computing the number of nonzeros
        let mut i_vector: HighsInt = 0;
        while i_vector < *num_sub_vector {
            *sub_matrix_num_nz += sub_matrix_length[i_vector as usize];
            i_vector += 1;
        }
        return;
    }
    // Allocate an array of lengths for the sub-matrix to be extracted
    sub_matrix_start.as_mut().unwrap()[0] = 0;
    let mut i_vector: HighsInt = 0;
    while i_vector < *num_sub_vector - 1 {
        sub_matrix_start.as_mut().unwrap()[(i_vector + 1) as usize] =
            sub_matrix_start.as_ref().unwrap()[i_vector as usize]
                + sub_matrix_length[i_vector as usize];
        sub_matrix_length[i_vector as usize] =
            sub_matrix_start.as_ref().unwrap()[i_vector as usize];
        i_vector += 1;
    }
    let i_vector = *num_sub_vector - 1;
    *sub_matrix_num_nz = sub_matrix_start.as_ref().unwrap()[i_vector as usize]
        + sub_matrix_length[i_vector as usize];
    // Bail out if matrix indices and values are not required
    if !extract_matrix {
        return;
    }
    sub_matrix_length[i_vector as usize] = sub_matrix_start.as_ref().unwrap()[i_vector as usize];
    // Fill the row-wise matrix with indices and values
    let mut vector: HighsInt = 0;
    while vector < num_vector {
        let mut i_el: HighsInt = matrix.start_[vector as usize];
        while i_el < matrix.start_[(vector + 1) as usize] {
            let i_vector = matrix.index_[i_el as usize];
            let new_i_vector = new_index[i_vector as usize];
            if new_i_vector >= 0 {
                let row_i_el = sub_matrix_length[new_i_vector as usize];
                if extract_index {
                    sub_matrix_index.as_mut().unwrap()[row_i_el as usize] = vector;
                }
                if extract_value {
                    sub_matrix_value.as_mut().unwrap()[row_i_el as usize] =
                        matrix.value_[i_el as usize];
                }
                sub_matrix_length[new_i_vector as usize] += 1;
            }
            i_el += 1;
        }
        vector += 1;
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
    fn change_and_get_bounds() {
        let mut highs = quiet(iis_lp());
        assert_eq!(highs.change_col_bounds(1, -1.0, 2.0), HighsStatus::kOk);
        assert_eq!(highs.change_row_bounds(0, 0.0, 3.0), HighsStatus::kOk);
        assert_eq!(
            highs.change_cols_cost_interval(0, 1, &[5.0, 6.0]),
            HighsStatus::kOk
        );
        let (mut num_col, mut num_nz) = (0, 0);
        let (mut cost, mut lo, mut up) = (vec![0.0; 1], vec![0.0; 1], vec![0.0; 1]);
        assert_eq!(
            highs.get_cols_interval(
                1,
                1,
                &mut num_col,
                Some(&mut cost),
                Some(&mut lo),
                Some(&mut up),
                &mut num_nz,
                None,
                None,
                None
            ),
            HighsStatus::kOk
        );
        assert_eq!(
            (num_col, num_nz, cost[0], lo[0], up[0]),
            (1, 2, 6.0, -1.0, 2.0)
        );
        assert_eq!(
            (highs.model_lp_.row_lower_[0], highs.model_lp_.row_upper_[0]),
            (0.0, 3.0)
        );
    }

    /// Model of check/TestLpModification.cpp "no-modification-on-error": one
    /// integer column (cost -1, [0, 1]) in two rows, both rows [0, inf].
    fn no_modify_lp() -> HighsLp {
        let mut model = lp(
            &[(-1.0, 0.0, 1.0, &[(0, 1.0), (1, 1.0)])],
            &[(0.0, kHighsInf), (0.0, kHighsInf)],
        );
        model.integrality_ = vec![HighsVarType::kInteger];
        model.col_names_ = vec!["C".to_string()];
        model.row_names_ = vec!["R0".to_string(), "R1".to_string()];
        model
    }

    /// The fields an LP modification could touch, for "the LP is unchanged".
    #[derive(Debug, PartialEq)]
    struct LpState {
        num_col: HighsInt,
        num_row: HighsInt,
        col_cost: Vec<f64>,
        col_lower: Vec<f64>,
        col_upper: Vec<f64>,
        row_lower: Vec<f64>,
        row_upper: Vec<f64>,
        a_start: Vec<HighsInt>,
        a_index: Vec<HighsInt>,
        a_value: Vec<f64>,
        col_names: Vec<String>,
        row_names: Vec<String>,
        integrality: Vec<HighsVarType>,
    }

    fn lp_state(model: &HighsLp) -> LpState {
        LpState {
            num_col: model.num_col_,
            num_row: model.num_row_,
            col_cost: model.col_cost_.clone(),
            col_lower: model.col_lower_.clone(),
            col_upper: model.col_upper_.clone(),
            row_lower: model.row_lower_.clone(),
            row_upper: model.row_upper_.clone(),
            a_start: model.a_matrix_.start_.clone(),
            a_index: model.a_matrix_.index_.clone(),
            a_value: model.a_matrix_.value_.clone(),
            col_names: model.col_names_.clone(),
            row_names: model.row_names_.clone(),
            integrality: model.integrality_.clone(),
        }
    }

    fn no_modify_highs(flag: bool) -> Highs {
        let mut highs = quiet(no_modify_lp());
        assert_eq!(
            highs.set_option_value_bool("h16_add_no_modify_on_error", flag),
            HighsStatus::kOk
        );
        highs
    }

    #[test]
    fn h16_addnomod_on_illegal_changes_keep_lp() {
        let mut highs = no_modify_highs(true);
        let before = lp_state(highs.get_lp());
        // A cost cannot be changed to an illegal value
        assert_eq!(highs.change_col_cost(-1, kHighsInf), HighsStatus::kError);
        assert_eq!(lp_state(highs.get_lp()), before);
        assert_eq!(highs.change_col_cost(0, kHighsInf), HighsStatus::kOk);
        assert_eq!(highs.change_col_cost(0, -1.0), HighsStatus::kOk);
        assert_eq!(lp_state(highs.get_lp()), before);

        assert_eq!(
            highs.change_col_bounds(-1, kHighsInf, -kHighsInf),
            HighsStatus::kError
        );
        assert_eq!(lp_state(highs.get_lp()), before);
        assert_eq!(
            highs.change_col_bounds(0, -kHighsInf, -kHighsInf),
            HighsStatus::kError
        );
        assert_eq!(lp_state(highs.get_lp()), before);
        // Inconsistent bounds are a warning, and the bounds are restored after
        assert_eq!(highs.change_col_bounds(0, 1.0, -1.0), HighsStatus::kWarning);
        assert_eq!(highs.change_col_bounds(0, 0.0, 1.0), HighsStatus::kOk);
        assert_eq!(lp_state(highs.get_lp()), before);

        assert_eq!(
            highs.change_row_bounds(-1, kHighsInf, -kHighsInf),
            HighsStatus::kError
        );
        assert_eq!(lp_state(highs.get_lp()), before);
        assert_eq!(
            highs.change_row_bounds(0, -kHighsInf, -kHighsInf),
            HighsStatus::kError
        );
        assert_eq!(lp_state(highs.get_lp()), before);
        assert_eq!(highs.change_row_bounds(0, 1.0, -1.0), HighsStatus::kWarning);
        assert_eq!(highs.change_row_bounds(0, 0.0, kHighsInf), HighsStatus::kOk);
        assert_eq!(lp_state(highs.get_lp()), before);

        // An illegal integrality change is refused
        assert_eq!(
            highs.change_col_integrality(-1, HighsVarType::kContinuous),
            HighsStatus::kError
        );
        assert_eq!(lp_state(highs.get_lp()), before);
        assert_eq!(
            highs.change_col_integrality(0, HighsVarType::kContinuous),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.change_col_integrality(0, HighsVarType::kInteger),
            HighsStatus::kOk
        );
        assert_eq!(lp_state(highs.get_lp()), before);

        // A column with an illegal row index is refused and the LP is unchanged
        assert_eq!(
            highs.add_col(1.0, 0.0, kHighsInf, 1, &[2], &[1.0]),
            HighsStatus::kError
        );
        assert_eq!(lp_state(highs.get_lp()), before);
        // A row with an illegal column index is refused and the LP is unchanged
        assert_eq!(
            highs.add_row(0.0, kHighsInf, 1, &[1], &[1.0]),
            HighsStatus::kError
        );
        assert_eq!(lp_state(highs.get_lp()), before);

        // Legal additions still append to the LP
        assert_eq!(
            highs.add_col(1.0, 0.0, kHighsInf, 1, &[0], &[1.0]),
            HighsStatus::kOk
        );
        assert_eq!(highs.get_lp().num_col_, 2);
        assert_eq!(highs.get_lp().col_cost_.len(), 2);
        assert_eq!(
            highs.add_row(0.0, kHighsInf, 1, &[1], &[1.0]),
            HighsStatus::kOk
        );
        assert_eq!(highs.get_lp().num_row_, 3);
        assert_eq!(highs.get_lp().row_lower_.len(), 3);
    }

    #[test]
    fn h16_addnomod_off_illegal_add_col_grows_lp() {
        let mut highs = no_modify_highs(false);
        // 1.15.1 appends the column to the LP vectors before the matrix is
        // assessed, so the illegal add leaves the vectors grown while num_col_
        // stays 1 (known 1.15.1 corruption, kept on the off path)
        assert_eq!(
            highs.add_col(1.0, 0.0, kHighsInf, 1, &[2], &[1.0]),
            HighsStatus::kError
        );
        assert_eq!(highs.get_lp().num_col_, 1);
        assert_eq!(highs.get_lp().col_cost_.len(), 2);
        assert_eq!(highs.get_lp().col_lower_.len(), 2);
        assert_eq!(highs.get_lp().col_upper_.len(), 2);
        // 1.15.1 returns kOk for inconsistent bounds
        assert_eq!(highs.change_col_bounds(0, 1.0, -1.0), HighsStatus::kOk);
    }
}
