//! HiGHS source: lp_data/Highs.cpp (addCol, addVars, addRow,
//! changeObjectiveSense, changeObjectiveOffset, changeColIntegrality,
//! changeColsIntegrality (set, mask), changeColCost, changeColsCost (set),
//! changeColsBounds (mask), changeRowsBounds (mask), deleteCols (set, mask),
//! deleteRows (set), changeCoeff), lp_data/HighsInterface.cpp
//! (changeCoefficientInterface), lp_data/HighsLpUtils.cpp
//! (changeLpMatrixCoefficient), Highs.h (addVar)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The remaining model-modification overloads of the public
//! API. The other overloads and the `*Interface` methods they share are
//! already ported (highs_api.rs, highs_api_cols.rs, highs_api_get.rs,
//! highs_api_mip.rs, highs_api_misc.rs, highs_interface_bounds.rs,
//! highs_interface_rows.rs). `logHeader` is dropped (as in add_cols).
//! `sortSetData` / `create` follow the existing set overloads (e.g.
//! `change_cols_bounds_set`, highs_api.rs). The `static`
//! `analyseSetCreateError` (Highs.cpp 2957-2993) and the `HighsVarType`
//! overload of `sortSetData` (HighsSort.cpp 346-364) are ported here as
//! private helpers.

use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::addCol (lines 2846-2853).
    /// `HighsInt starts = 0` -> `&[0]`.
    pub fn add_col(
        &mut self,
        cost: f64,
        lower_bound: f64,
        upper_bound: f64,
        num_new_nz: HighsInt,
        indices: &[HighsInt],
        values: &[f64],
    ) -> HighsStatus {
        // logHeader dropped (lines 2849)
        let starts = [0 as HighsInt];
        self.add_cols(
            1,
            &[cost],
            &[lower_bound],
            &[upper_bound],
            num_new_nz,
            &starts,
            indices,
            values,
        )
    }

    /// HiGHS: Highs.h::Highs::addVar (lines 1118-1120).
    /// C++ default arguments: `lower = 0`, `upper = kHighsInf`.
    pub fn add_var(&mut self, lower: f64, upper: f64) -> HighsStatus {
        self.add_vars(1, &[lower], &[upper])
    }

    /// HiGHS: Highs.cpp::Highs::addVars (lines 2872-2882)
    pub fn add_vars(&mut self, num_new_var: HighsInt, lower: &[f64], upper: &[f64]) -> HighsStatus {
        // logHeader dropped (line 2874)
        let mut return_status = HighsStatus::kOk;
        // Avoid touching entry [0] of a vector of size 0
        if num_new_var <= 0 {
            return self.return_from_highs(return_status);
        }
        let mut cost = vec![0.0f64; num_new_var as usize];
        self.add_cols(num_new_var, &cost, lower, upper, 0, &[], &[], &[])
    }

    /// HiGHS: Highs.cpp::Highs::addRow (lines 2884-2891)
    pub fn add_row(
        &mut self,
        lower_bound: f64,
        upper_bound: f64,
        num_new_nz: HighsInt,
        indices: &[HighsInt],
        values: &[f64],
    ) -> HighsStatus {
        // logHeader dropped (line 2887)
        let starts = [0 as HighsInt];
        self.add_rows(
            1,
            &[lower_bound],
            &[upper_bound],
            num_new_nz,
            &starts,
            indices,
            values,
        )
    }

    /// HiGHS: Highs.cpp::Highs::changeObjectiveSense (lines 2910-2919)
    pub fn change_objective_sense(&mut self, sense: ObjSense) -> HighsStatus {
        if (sense == ObjSense::kMinimize) != (self.model_lp_.sense_ == ObjSense::kMinimize) {
            self.model_lp_.sense_ = sense;
            // Nontrivial change
            self.clear_derived_model_properties();
            self.invalidate_model_status_solution_and_info();
        }
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: Highs.cpp::Highs::changeObjectiveOffset (lines 2921-2927).
    /// HiGHS quirk: `presolved_model_.lp_.offset_ += offset` adds the new
    /// offset rather than the change.
    pub fn change_objective_offset(&mut self, offset: f64) -> HighsStatus {
        // Update the objective value
        self.info_.objective_function_value += offset - self.model_lp_.offset_;
        self.model_lp_.offset_ = offset;
        self.presolved_model_lp_.offset_ += offset;
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: Highs.cpp::Highs::changeColIntegrality (lines 2929-2932)
    pub fn change_col_integrality(
        &mut self,
        col: HighsInt,
        integrality: HighsVarType,
    ) -> HighsStatus {
        self.change_cols_integrality_set(1, &[col], &[integrality])
    }

    /// HiGHS: Highs.cpp::Highs::changeColsIntegrality (lines 2995-3020), the
    /// set overload.
    pub fn change_cols_integrality_set(
        &mut self,
        num_set_entries: HighsInt,
        set: &[HighsInt],
        integrality: &[HighsVarType],
    ) -> HighsStatus {
        if num_set_entries == 0 {
            return HighsStatus::kOk;
        }
        self.clear_presolve();
        // Ensure that the set and data are in ascending order
        let mut local_integrality: Vec<HighsVarType> =
            integrality[..num_set_entries as usize].to_vec();
        let mut local_set: Vec<HighsInt> = set[..num_set_entries as usize].to_vec();
        sort_set_data_var_type(
            num_set_entries,
            &mut local_set,
            Some(integrality),
            Some(&mut local_integrality),
        );
        let mut index_collection = HighsIndexCollection::default();
        let create_error = crate::util::highs_utils::create_set(
            &mut index_collection,
            num_set_entries,
            &local_set,
            self.model_lp_.num_col_,
        );
        if create_error != 0 {
            return analyse_set_create_error(
                &self.options_.log_options,
                "changeColsIntegrality",
                create_error,
                true,
                num_set_entries,
                &local_set,
                self.model_lp_.num_col_,
            );
        }
        let call_status =
            self.change_integrality_interface(&mut index_collection, &local_integrality);
        let mut return_status = HighsStatus::kOk;
        return_status = crate::lp_data::highs_status::interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "changeIntegrality",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::changeColsIntegrality (lines 3022-3036), the
    /// mask overload.
    pub fn change_cols_integrality_mask(
        &mut self,
        mask: &[HighsInt],
        integrality: &[HighsVarType],
    ) -> HighsStatus {
        self.clear_presolve();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = crate::util::highs_utils::create_mask(
            &mut index_collection,
            mask,
            self.model_lp_.num_col_,
        );
        debug_assert!(create_error == 0);
        let _ = create_error;
        let call_status = self.change_integrality_interface(&mut index_collection, integrality);
        let mut return_status = HighsStatus::kOk;
        return_status = crate::lp_data::highs_status::interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "changeIntegrality",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::changeColCost (lines 3038-3040)
    pub fn change_col_cost(&mut self, col: HighsInt, cost: f64) -> HighsStatus {
        self.change_cols_cost_set(1, &[col], &[cost])
    }

    /// HiGHS: Highs.cpp::Highs::changeColsCost (lines 3063-3090), the set
    /// overload.
    pub fn change_cols_cost_set(
        &mut self,
        num_set_entries: HighsInt,
        set: &[HighsInt],
        cost: &[f64],
    ) -> HighsStatus {
        if num_set_entries == 0 {
            return HighsStatus::kOk;
        }
        // Check for NULL data in "set" version of changeColsCost since
        // values are sorted with set
        // doubleUserDataNotNull is not ported; assume non-null as in Rust slices
        self.clear_derived_model_properties();
        // Ensure that the set and data are in ascending order
        let mut local_cost: Vec<f64> = cost[..num_set_entries as usize].to_vec();
        let mut local_set: Vec<HighsInt> = set[..num_set_entries as usize].to_vec();
        crate::util::highs_sort_set::sort_set_data(
            num_set_entries,
            &mut local_set,
            Some(cost),
            None,
            None,
            Some(&mut local_cost),
            None,
            None,
        );
        let mut index_collection = HighsIndexCollection::default();
        let create_error = crate::util::highs_utils::create_set(
            &mut index_collection,
            num_set_entries,
            &local_set,
            self.model_lp_.num_col_,
        );
        if create_error != 0 {
            return analyse_set_create_error(
                &self.options_.log_options,
                "changeColsCost",
                create_error,
                true,
                num_set_entries,
                &local_set,
                self.model_lp_.num_col_,
            );
        }
        let call_status = self.change_costs_interface(&index_collection, &local_cost);
        let mut return_status = HighsStatus::kOk;
        return_status = crate::lp_data::highs_status::interpret_call_status(
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

    /// HiGHS: Highs.cpp::Highs::changeColsBounds (lines 3171-3185), the mask
    /// overload.
    pub fn change_cols_bounds_mask(
        &mut self,
        mask: &[HighsInt],
        lower: &[f64],
        upper: &[f64],
    ) -> HighsStatus {
        self.clear_derived_model_properties();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = crate::util::highs_utils::create_mask(
            &mut index_collection,
            mask,
            self.model_lp_.num_col_,
        );
        debug_assert!(create_error == 0);
        let _ = create_error;
        let call_status = self.change_col_bounds_interface(&mut index_collection, lower, upper);
        let mut return_status = HighsStatus::kOk;
        return_status = crate::lp_data::highs_status::interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "changeColBounds",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::changeRowsBounds (lines 3252-3266), the mask
    /// overload.
    pub fn change_rows_bounds_mask(
        &mut self,
        mask: &[HighsInt],
        lower: &[f64],
        upper: &[f64],
    ) -> HighsStatus {
        self.clear_derived_model_properties();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = crate::util::highs_utils::create_mask(
            &mut index_collection,
            mask,
            self.model_lp_.num_row_,
        );
        debug_assert!(create_error == 0);
        let _ = create_error;
        let call_status = self.change_row_bounds_interface(&mut index_collection, lower, upper);
        let mut return_status = HighsStatus::kOk;
        return_status = crate::lp_data::highs_status::interpret_call_status(
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

    /// HiGHS: Highs.cpp::Highs::deleteCols (lines 3509-3522), the set
    /// overload.
    pub fn delete_cols_set(&mut self, num_set_entries: HighsInt, set: &[HighsInt]) -> HighsStatus {
        if num_set_entries == 0 {
            return HighsStatus::kOk;
        }
        self.clear_derived_model_properties();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = crate::util::highs_utils::create_set(
            &mut index_collection,
            num_set_entries,
            set,
            self.model_lp_.num_col_,
        );
        if create_error != 0 {
            return analyse_set_create_error(
                &self.options_.log_options,
                "deleteCols",
                create_error,
                false,
                num_set_entries,
                set,
                self.model_lp_.num_col_,
            );
        }
        self.delete_cols_interface(&mut index_collection);
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: Highs.cpp::Highs::deleteCols (lines 3524-3535), the mask
    /// overload: `mask` is overwritten with the new indices (-1 if deleted).
    pub fn delete_cols_mask(&mut self, mask: &mut [HighsInt]) -> HighsStatus {
        self.clear_derived_model_properties();
        let original_num_col = self.model_lp_.num_col_;
        let mut index_collection = HighsIndexCollection::default();
        let create_error =
            crate::util::highs_utils::create_mask(&mut index_collection, mask, original_num_col);
        debug_assert!(create_error == 0);
        let _ = create_error;
        self.delete_cols_interface(&mut index_collection);
        for i_col in 0..original_num_col {
            mask[i_col as usize] = index_collection.mask_[i_col as usize];
        }
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: Highs.cpp::Highs::deleteRows (lines 3553-3566), the set
    /// overload.
    pub fn delete_rows_set(&mut self, num_set_entries: HighsInt, set: &[HighsInt]) -> HighsStatus {
        if num_set_entries == 0 {
            return HighsStatus::kOk;
        }
        self.clear_derived_model_properties();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = crate::util::highs_utils::create_set(
            &mut index_collection,
            num_set_entries,
            set,
            self.model_lp_.num_row_,
        );
        if create_error != 0 {
            return analyse_set_create_error(
                &self.options_.log_options,
                "deleteRows",
                create_error,
                false,
                num_set_entries,
                set,
                self.model_lp_.num_row_,
            );
        }
        self.delete_rows_interface(&mut index_collection);
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: Highs.cpp::Highs::changeCoeff (lines 3268-3295).
    /// `"%" HIGHSINT_FORMAT` -> `{}`.
    pub fn change_coeff(&mut self, row: HighsInt, col: HighsInt, value: f64) -> HighsStatus {
        if row < 0 || row >= self.model_lp_.num_row_ {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Row {} supplied to Highs::changeCoeff is not in the range [0, {}]\n",
                row,
                self.model_lp_.num_row_
            );
            return HighsStatus::kError;
        }
        if col < 0 || col >= self.model_lp_.num_col_ {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Col {} supplied to Highs::changeCoeff is not in the range [0, {}]\n",
                col,
                self.model_lp_.num_col_
            );
            return HighsStatus::kError;
        }
        let abs_value = value.abs();
        if 0.0 < abs_value && abs_value <= self.options_.small_matrix_value {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "|Value| of {} supplied to Highs::changeCoeff is in (0, {}]: zeroes any existing coefficient, otherwise ignored\n",
                crate::util::c_printf::c_fmt_g(abs_value, 6),
                crate::util::c_printf::c_fmt_g(self.options_.small_matrix_value, 6)
            );
        }
        self.change_coefficient_interface(row, col, value);
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: HighsInterface.cpp::Highs::changeCoefficientInterface (lines 1061-1093)
    pub(crate) fn change_coefficient_interface(
        &mut self,
        ext_row: HighsInt,
        ext_col: HighsInt,
        ext_new_value: f64,
    ) {
        // HighsLp& lp = model_.lp_;
        // lp.ensureColwise();
        self.model_lp_.ensure_colwise();
        debug_assert!(0 <= ext_row && ext_row < self.model_lp_.num_row_);
        debug_assert!(0 <= ext_col && ext_col < self.model_lp_.num_col_);
        let zero_new_value = ext_new_value.abs() <= self.options_.small_matrix_value;
        change_lp_matrix_coefficient(
            &mut self.model_lp_,
            ext_row,
            ext_col,
            ext_new_value,
            zero_new_value,
        );
        // Deduce the consequences of a changed element
        //
        // ToDo: Can do something more intelligent if element is in nonbasic
        // column
        //
        let basic_column = self.basis_.col_status[ext_col as usize] == HighsBasisStatus::kBasic;
        //
        // For now, treat it as if it's a new row
        self.invalidate_model_status_solution_and_info();

        if basic_column {
            // Basis is retained, but is has to be viewed as alien, since the
            // basis matrix has changed
            self.basis_.was_alien = true;
            self.basis_.alien = true;
        }

        // Determine any implications for simplex data
        self.ekk_instance_.update_status(LpAction::kNewRows);
    }
}

/// HiGHS: Highs.cpp::analyseSetCreateError (lines 2957-2993), a `static`
/// function. `HighsLogOptions log_options` (by value) -> `&HighsLogOptions`;
/// `%s` / `%d` -> `{}`.
fn analyse_set_create_error(
    log_options: &HighsLogOptions,
    method: &str,
    create_error: HighsInt,
    ordered: bool,
    num_set_entries: HighsInt,
    set: &[HighsInt],
    dimension: HighsInt,
) -> HighsStatus {
    // HighsUtils.h: kIndexCollectionCreateIllegalSetSize = 1,
    // kIndexCollectionCreateIllegalSetDimension = 2,
    // kIndexCollectionCreateIllegalSetOrder = 3
    const kIndexCollectionCreateIllegalSetSize: HighsInt = 1;
    const kIndexCollectionCreateIllegalSetDimension: HighsInt = 2;
    const kIndexCollectionCreateIllegalSetOrder: HighsInt = 3;
    if create_error == kIndexCollectionCreateIllegalSetSize {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "Set supplied to Highs::{} has illegal size of {}\n",
            method,
            num_set_entries
        );
    } else if create_error == kIndexCollectionCreateIllegalSetOrder {
        if ordered {
            // Creating an index_collection data structure for the set
            // includes a test that the indices increase strictly. If this
            // is not the case then, since an increasing set was created
            // locally, it must contain duplicate entries. return with an
            // error
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Set supplied to Highs::{} contains duplicate entries\n",
                method
            );
        } else {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Set supplied to Highs::{} not ordered\n",
                method
            );
        }
    } else if create_error < 0 {
        let illegal_set_index: HighsInt = -1 - create_error;
        let illegal_set_entry: HighsInt = set[illegal_set_index as usize];
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "Set supplied to Highs::{} has entry {} of {} out of range [0, {})\n",
            method,
            illegal_set_index,
            illegal_set_entry,
            dimension
        );
    }
    debug_assert!(create_error != kIndexCollectionCreateIllegalSetDimension);
    HighsStatus::kError
}

/// HiGHS: HighsSort.cpp::sortSetData (lines 346-364), the `HighsVarType`
/// overload (not ported in util/highs_sort_set.rs, which has only the
/// `double` one). `const HighsVarType* data0` -> `Option<&[HighsVarType]>`,
/// `HighsVarType* sorted_data0` -> `Option<&mut [HighsVarType]>`.
fn sort_set_data_var_type(
    num_entries: HighsInt,
    set: &mut [HighsInt],
    data0: Option<&[HighsVarType]>,
    sorted_data0: Option<&mut [HighsVarType]>,
) {
    if num_entries <= 0 {
        return;
    }
    let mut sort_set_vec: Vec<HighsInt> = vec![0; (1 + num_entries) as usize];
    let mut perm_vec: Vec<HighsInt> = vec![0; (1 + num_entries) as usize];

    let sort_set: &mut [HighsInt] = &mut sort_set_vec;
    let perm: &mut [HighsInt] = &mut perm_vec;

    for ix in 0..num_entries {
        sort_set[(1 + ix) as usize] = set[ix as usize];
        perm[(1 + ix) as usize] = ix;
    }
    crate::util::highs_sort_set::maxheapsort_int(sort_set, perm, num_entries);
    let mut sorted0 = sorted_data0;
    for ix in 0..num_entries {
        set[ix as usize] = sort_set[(1 + ix) as usize];
        if let Some(d) = data0 {
            sorted0
                .as_deref_mut()
                .expect("sortSetData: sorted_data0 missing")[ix as usize] =
                d[perm[(1 + ix) as usize] as usize];
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::changeLpMatrixCoefficient (lines 1668-1713)
pub fn change_lp_matrix_coefficient(
    lp: &mut HighsLp,
    row: HighsInt,
    col: HighsInt,
    new_value: f64,
    zero_new_value: bool,
) {
    debug_assert!(0 <= row && row < lp.num_row_);
    debug_assert!(0 <= col && col < lp.num_col_);

    // Determine whether the coefficient corresponds to an existing
    // nonzero
    let mut change_el: HighsInt = -1;
    let mut el = lp.a_matrix_.start_[col as usize];
    while el < lp.a_matrix_.start_[(col + 1) as usize] {
        if lp.a_matrix_.index_[el as usize] == row {
            change_el = el;
            break;
        }
        el += 1;
    }
    if change_el < 0 {
        // Coefficient doesn't correspond to an existing nonzero
        //
        // If coefficient is small, then just ignore it
        if zero_new_value {
            return;
        }
        // New nonzero goes at the end of column "col", so have to shift
        // all index and value entries forward by 1 to accommodate it
        change_el = lp.a_matrix_.start_[(col + 1) as usize];
        let new_num_nz = lp.a_matrix_.start_[lp.num_col_ as usize] + 1;
        lp.a_matrix_.index_.resize(new_num_nz as usize, 0);
        lp.a_matrix_.value_.resize(new_num_nz as usize, 0.0);
        let mut i = col + 1;
        while i <= lp.num_col_ {
            lp.a_matrix_.start_[i as usize] += 1;
            i += 1;
        }
        let mut el = new_num_nz - 1;
        while el > change_el {
            lp.a_matrix_.index_[el as usize] = lp.a_matrix_.index_[(el - 1) as usize];
            lp.a_matrix_.value_[el as usize] = lp.a_matrix_.value_[(el - 1) as usize];
            el -= 1;
        }
    } else if zero_new_value {
        // Coefficient zeroes an existing nonzero, so shift all index and
        // value entries backward by 1 to eliminate it
        let new_num_nz = lp.a_matrix_.start_[lp.num_col_ as usize] - 1;
        let mut i = col + 1;
        while i <= lp.num_col_ {
            lp.a_matrix_.start_[i as usize] -= 1;
            i += 1;
        }
        let mut el = change_el;
        while el < new_num_nz {
            lp.a_matrix_.index_[el as usize] = lp.a_matrix_.index_[(el + 1) as usize];
            lp.a_matrix_.value_[el as usize] = lp.a_matrix_.value_[(el + 1) as usize];
            el += 1;
        }
        return;
    }
    lp.a_matrix_.index_[change_el as usize] = row;
    lp.a_matrix_.value_[change_el as usize] = new_value;
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// Two columns, one row x0 + 2 x1 in [1, 8].
    fn lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![4.0, 3.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![8.0];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 2.0];
        lp
    }

    #[test]
    fn matrix_coefficient_change_insert_and_zero() {
        let mut lp = lp();
        change_lp_matrix_coefficient(&mut lp, 0, 1, 5.0, false);
        assert_eq!(lp.a_matrix_.value_, vec![1.0, 5.0]);
        change_lp_matrix_coefficient(&mut lp, 0, 0, 0.0, true);
        assert_eq!(lp.a_matrix_.start_, vec![0, 0, 1]);
        // C++ (HighsLpUtils.cpp 1703-1708) shifts entries back but does not
        // shrink the vectors; only the first `start_[num_col]` entries are live.
        assert_eq!(lp.a_matrix_.value_[..1], [5.0]);
        change_lp_matrix_coefficient(&mut lp, 0, 0, 3.0, false);
        assert_eq!(lp.a_matrix_.start_, vec![0, 1, 2]);
        assert_eq!(lp.a_matrix_.value_, vec![3.0, 5.0]);
    }

    #[test]
    fn change_coeff_checks_indices() {
        let mut highs = Highs::default();
        highs.model_lp_ = lp();
        // HiGHS keeps basis_.col_status sized to the model; the default
        // basis is empty, so size it as the model has it.
        highs.basis_.col_status = vec![HighsBasisStatus::kNonbasic; 2];
        assert_eq!(highs.change_coeff(1, 0, 1.0), HighsStatus::kError);
        assert_eq!(highs.change_coeff(0, 2, 1.0), HighsStatus::kError);
        assert_eq!(highs.change_coeff(0, 1, -1.0), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.a_matrix_.value_, vec![1.0, -1.0]);
    }

    #[test]
    fn add_col_row_and_vars_grow_the_model() {
        let mut highs = Highs::default();
        highs.model_lp_ = lp();
        assert_eq!(
            highs.add_col(3.0, 0.0, 1.0, 1, &[0], &[4.0]),
            HighsStatus::kOk
        );
        assert_eq!(highs.model_lp_.num_col_, 3);
        assert_eq!(
            highs.add_row(-kHighsInf, 2.0, 2, &[0, 2], &[1.0, 1.0]),
            HighsStatus::kOk
        );
        assert_eq!(highs.model_lp_.num_row_, 2);
        assert_eq!(
            highs.add_vars(2, &[0.0, 0.0], &[1.0, 2.0]),
            HighsStatus::kOk
        );
        assert_eq!(highs.model_lp_.num_col_, 5);
        assert_eq!(highs.model_lp_.col_cost_[3..], [0.0, 0.0]);
        assert_eq!(highs.add_var(0.0, kHighsInf), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.num_col_, 6);
    }

    #[test]
    fn sense_offset_cost_and_integrality() {
        let mut highs = Highs::default();
        highs.model_lp_ = lp();
        assert_eq!(
            highs.change_objective_sense(ObjSense::kMaximize),
            HighsStatus::kOk
        );
        assert_eq!(highs.model_lp_.sense_, ObjSense::kMaximize);
        assert_eq!(highs.change_objective_offset(2.5), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.offset_, 2.5);
        assert_eq!(highs.change_col_cost(1, -3.0), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.col_cost_, vec![1.0, -3.0]);
        assert_eq!(
            highs.change_col_integrality(0, HighsVarType::kInteger),
            HighsStatus::kOk
        );
        assert_eq!(highs.model_lp_.integrality_[0], HighsVarType::kInteger);
        assert_eq!(
            highs.change_cols_cost_set(2, &[1, 0], &[7.0, 6.0]),
            HighsStatus::kOk
        );
        assert_eq!(highs.model_lp_.col_cost_, vec![6.0, 7.0]);
    }

    #[test]
    fn integrality_set_is_sorted_with_its_data() {
        let mut highs = Highs::default();
        highs.model_lp_ = lp();
        // sortSetData (HighsSort.cpp 346-364) reorders the data with the set
        assert_eq!(
            highs.change_cols_integrality_set(
                2,
                &[1, 0],
                &[HighsVarType::kInteger, HighsVarType::kContinuous]
            ),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.model_lp_.integrality_,
            vec![HighsVarType::kContinuous, HighsVarType::kInteger]
        );
        // Duplicates: create fails with kIndexCollectionCreateIllegalSetOrder
        assert_eq!(
            highs.change_cols_integrality_set(
                2,
                &[1, 1],
                &[HighsVarType::kInteger, HighsVarType::kInteger]
            ),
            HighsStatus::kError
        );
    }

    #[test]
    fn set_create_errors_return_error() {
        let mut highs = Highs::default();
        highs.model_lp_ = lp();
        // Out of range entry (create_error < 0)
        assert_eq!(
            highs.change_cols_cost_set(1, &[2], &[1.0]),
            HighsStatus::kError
        );
        // Unordered set for deleteCols (not sorted locally)
        assert_eq!(highs.delete_cols_set(2, &[1, 0]), HighsStatus::kError);
        assert_eq!(highs.model_lp_.num_col_, 2);
        assert_eq!(highs.delete_rows_set(1, &[1]), HighsStatus::kError);
        assert_eq!(highs.model_lp_.num_row_, 1);
    }

    #[test]
    fn mask_and_set_deletions() {
        let mut highs = Highs::default();
        highs.model_lp_ = lp();
        let mut mask = vec![1, 0];
        assert_eq!(highs.delete_cols_mask(&mut mask), HighsStatus::kOk);
        assert_eq!(mask, vec![-1, 0]);
        assert_eq!(highs.model_lp_.num_col_, 1);
        assert_eq!(highs.delete_rows_set(1, &[0]), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.num_row_, 0);
        assert_eq!(highs.delete_cols_set(1, &[0]), HighsStatus::kOk);
        assert_eq!(highs.model_lp_.num_col_, 0);
    }
}
