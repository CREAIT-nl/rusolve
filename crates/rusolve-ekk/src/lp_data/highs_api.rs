//! HiGHS source: lp_data/Highs.cpp (optimizeLp, getDualRaySparse, hasInvert,
//! getBasicVariablesArray, getBasisInverseRowSparse, setBasis, putIterate,
//! getIterate, addRows, changeColsBounds(interval, set), deleteRows(interval,
//! mask)),
//! lp_data/Highs.h (getLp, getNumCol, getNumRow, getRunTime)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! API slice: the public `Highs` methods that
//! `HighsLpRelaxation` calls on its `lpsolver` (mip/HighsLpRelaxation.cpp,
//! .h), which is what `rusolve_mip::engine::ekk::EkkEngine` drives.
//! Overloads get suffixed names (`change_cols_bounds_set`,
//! `delete_rows_mask`). `setOptionValue` is not ported: callers write the
//! `options_` field directly (the only C++ side effect, `optionChangeAction`,
//! touches IIS data, which this port does not have).
//!
//! NLA rule (simplex/ekk/clear.rs `btran` caveat): after `run` the LP has
//! been moved back out of `ekk_instance_`, so `ekk_instance_.btran` (which
//! reads `ekk_instance_.lp_`) is WRONG here. C++
//! `ekk_instance_.setNlaPointersForLpAndScale(model_.lp_); ekk_instance_.btran(v, d)`
//! translates to
//! ```text
//! let refs = NlaRefs::new(&self.model_lp_, &self.ekk_instance_.basis_.basic_index_,
//!                         &self.ekk_instance_.options_);
//! self.ekk_instance_.simplex_nla_.btran(&refs, v, d, None);
//! ```

use crate::lp_data::highs_run::Highs;
use crate::prelude::*;
use crate::simplex::nla::NlaRefs;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::optimizeLp (lines 1145-1151, 7 lines)
    /// -> `calledOptimizeModel()` ([`Highs::called_optimize_model`]), not `run()`: level 0
    /// (`reportModelStats`, `assessExcessiveObjectiveBoundScaling`) is skipped, as in C++.
    pub fn optimize_lp(&mut self) -> HighsStatus {
        self.called_optimize_model()
    }

    /// HiGHS: Highs.h::Highs::getLp (line 543)
    pub fn get_lp(&self) -> &HighsLp {
        &self.model_lp_
    }

    /// HiGHS: Highs.h::Highs::getNumCol (line 757)
    pub fn get_num_col(&self) -> HighsInt {
        self.model_lp_.num_col_
    }

    /// HiGHS: Highs.h::Highs::getNumRow (line 762)
    pub fn get_num_row(&self) -> HighsInt {
        self.model_lp_.num_row_
    }

    /// HiGHS: Highs.h::Highs::getRunTime (line 483)
    pub fn get_run_time(&self) -> f64 {
        self.timer_.read()
    }

    /// HiGHS: Highs.cpp::Highs::hasInvert (line 2251)
    pub fn has_invert(&self) -> bool {
        self.ekk_instance_.status_.has_invert
    }

    /// HiGHS: Highs.cpp::Highs::getBasicVariablesArray (lines 2253-2256, 4 lines)
    pub fn get_basic_variables_array(&self) -> &[HighsInt] {
        debug_assert!(self.ekk_instance_.status_.has_invert);
        &self.ekk_instance_.basis_.basic_index_
    }

    /// HiGHS: Highs.cpp::Highs::getBasisInverseRowSparse (lines 2272-2284, 13 lines)
    ///
    /// `row_ep_buffer` must already be `setup(num_row)` (the caller,
    /// `computeBasicDegenerateDuals`, does that). Unit vector `e_row`
    /// (`clear`, `count = 1`, `index[0] = row`, `array[row] = 1`,
    /// `pack_flag = true`), then BTRAN with density
    /// `ekk_instance_.info_.row_ep_density`, using `NlaRefs` built from
    /// `self.model_lp_` (module docs).
    pub fn get_basis_inverse_row_sparse(
        &mut self,
        row: HighsInt,
        row_ep_buffer: &mut HVector,
    ) -> HighsStatus {
        let refs = NlaRefs::new(
            &self.model_lp_,
            &self.ekk_instance_.basis_.basic_index_,
            &self.ekk_instance_.options_,
        );
        row_ep_buffer.clear();
        row_ep_buffer.count = 1;
        row_ep_buffer.index[0] = row;
        row_ep_buffer.array[row as usize] = 1.0;
        row_ep_buffer.pack_flag = true;

        let density = self.ekk_instance_.info_.row_ep_density;
        self.ekk_instance_
            .simplex_nla_
            .btran(&refs, row_ep_buffer, density, None);

        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::getDualRaySparse (lines 2146-2162, 17 lines)
    ///
    /// `has_dual_ray = ekk_instance_.dual_ray_record_.index != kNoRayIndex`;
    /// if so, `row_ep_buffer` = `sign` times the unit vector at that index,
    /// BTRAN'd as in [`Self::get_basis_inverse_row_sparse`] (`array[iRow] =
    /// sign` is an `i32 -> f64` conversion: `sign as f64`).
    pub fn get_dual_ray_sparse(
        &mut self,
        has_dual_ray: &mut bool,
        row_ep_buffer: &mut HVector,
    ) -> HighsStatus {
        *has_dual_ray = self.ekk_instance_.dual_ray_record_.index != kNoRayIndex;
        if *has_dual_ray {
            let refs = NlaRefs::new(
                &self.model_lp_,
                &self.ekk_instance_.basis_.basic_index_,
                &self.ekk_instance_.options_,
            );
            row_ep_buffer.clear();
            row_ep_buffer.count = 1;
            row_ep_buffer.pack_flag = true;
            let i_row = self.ekk_instance_.dual_ray_record_.index;
            row_ep_buffer.index[0] = i_row;
            row_ep_buffer.array[i_row as usize] = self.ekk_instance_.dual_ray_record_.sign as f64;

            let density = self.ekk_instance_.info_.row_ep_density;
            self.ekk_instance_
                .simplex_nla_
                .btran(&refs, row_ep_buffer, density, None);
        }

        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::setBasis (lines 2732-2805, 74 lines)
    /// C++ default arguments: origin = ""
    ///
    /// * Alien branch with rows: `is_basis_right_size`
    ///   (lp_data/highs_solution.rs:65), a `modifiable_basis` clone with
    ///   `was_alien = true`, a `HighsLpSolverObject::new(&mut self.model_lp_,
    ///   &mut modifiable_basis, &mut self.solution_, &mut self.info_, &mut
    ///   self.callback_, &mut self.options_, &mut self.timer_)` and
    ///   `form_simplex_lp_basis_and_factor(&mut self.ekk_instance_, &mut
    ///   solver_object, false)` (highs_basis_form.rs). Profiling dropped.
    /// * Non-alien branch: `is_basis_consistent` (presolve_support.rs:293),
    ///   else log "setBasis: invalid basis" and return kError WITHOUT
    ///   touching `basis_` or the simplex state.
    /// * `if origin != "" { basis_.debug_origin_name = origin }`; the
    ///   `highsLogDev` report is dropped; then `new_highs_basis()`
    ///   (highs_api_support.rs). Returns kOk without `returnFromHighs`.
    pub fn set_basis(&mut self, basis: &HighsBasis, origin: &str) -> HighsStatus {
        if basis.alien {
            // An alien basis needs to be checked properly, since it may be
            // singular, or even incomplete.
            if self.model_lp_.num_row_ == 0 {
                // Special case where there are no rows, so no singularity
                // issues. All columns with basic status must be set nonbasic
                for i_col in 0..self.model_lp_.num_col_ as usize {
                    self.basis_.col_status[i_col] =
                        if basis.col_status[i_col] == HighsBasisStatus::kBasic {
                            HighsBasisStatus::kNonbasic
                        } else {
                            basis.col_status[i_col]
                        };
                }
                self.basis_.alien = false;
            } else {
                // Check whether a new basis can be defined
                if !crate::lp_data::highs_solution::is_basis_right_size(&self.model_lp_, basis) {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "setBasis: User basis is rejected due to mismatch between size of column and row status vectors ({}, {}) and number of columns and rows in the model ({}, {})\n",
                        self.basis_.col_status.len(),
                        self.basis_.row_status.len(),
                        self.model_lp_.num_col_,
                        self.model_lp_.num_row_
                    );
                    return HighsStatus::kError;
                }
                let mut modifiable_basis = basis.clone();
                modifiable_basis.was_alien = true;
                let mut solver_object = HighsLpSolverObject::new(
                    &mut self.model_lp_,
                    &mut modifiable_basis,
                    &mut self.solution_,
                    &mut self.info_,
                    &mut self.callback_,
                    &mut self.options_,
                    &mut self.timer_,
                );
                let return_status = form_simplex_lp_basis_and_factor(
                    &mut self.ekk_instance_,
                    &mut solver_object,
                    false,
                );
                if return_status != HighsStatus::kOk {
                    return HighsStatus::kError;
                }
                // Update the HiGHS basis
                self.basis_ = modifiable_basis;
            }
        } else {
            // Check the user-supplied basis
            if !crate::lp_data::presolve_support::is_basis_consistent(&self.model_lp_, basis) {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "setBasis: invalid basis\n"
                );
                return HighsStatus::kError;
            }
            // Update the HiGHS basis
            self.basis_ = basis.clone();
        }
        self.basis_.valid = true;
        self.basis_.useful = true;
        if !origin.is_empty() {
            self.basis_.debug_origin_name = origin.to_string();
        }
        debug_assert!(!self.basis_.debug_origin_name.is_empty());
        debug_assert!(!self.basis_.alien);

        // Follow implications of a new HiGHS basis
        self.new_highs_basis();
        // Can't use returnFromHighs since...
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::putIterate (lines 2819-2828, 10 lines)
    ///
    /// Error unless `ekk_instance_.status_.has_invert`; then
    /// `ekk_instance_.put_iterate()` (simplex/ekk/solve_return.rs:436) and
    /// `return_from_highs(kOk)`.
    pub fn put_iterate(&mut self) -> HighsStatus {
        // Check that there is a simplex iterate to put
        if !self.ekk_instance_.status_.has_invert {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "putIterate: no simplex iterate to put\n"
            );
            return HighsStatus::kError;
        }
        self.ekk_instance_.put_iterate();
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: Highs.cpp::Highs::getIterate (lines 2830-2843, 14 lines)
    ///
    /// Error unless `ekk_instance_.status_.initialised_for_new_lp`;
    /// `ekk_instance_.get_iterate()` (solve_return.rs:451); then
    /// `basis_ = ekk_instance_.get_highs_basis(&mut self.model_lp_)`
    /// (solve.rs:449), `invalidate_model_status_solution_and_info()`,
    /// `return_from_highs(kOk)`.
    pub fn get_iterate(&mut self) -> HighsStatus {
        // Check that there is a simplex iterate to get
        if !self.ekk_instance_.status_.initialised_for_new_lp {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "getIterate: no simplex iterate to get\n"
            );
            return HighsStatus::kError;
        }
        let call_status = self.ekk_instance_.get_iterate();
        if call_status != HighsStatus::kOk {
            return call_status;
        }
        // Get the corresponding HiGHS basis
        self.basis_ = self.ekk_instance_.get_highs_basis(&mut self.model_lp_);
        // Clear everything else
        self.invalidate_model_status_solution_and_info();
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: Highs.cpp::Highs::addRows (lines 2893-2907, 15 lines)
    ///
    /// `logHeader` dropped. `clear_derived_model_properties()`, then
    /// `add_rows_interface` (highs_interface_rows.rs) through
    /// `interpret_call_status(.., "addRows")`, then `return_from_highs`.
    /// `starts` has `num_new_row` entries (CSR without the end sentinel).
    pub fn add_rows(
        &mut self,
        num_new_row: HighsInt,
        lower_bounds: &[f64],
        upper_bounds: &[f64],
        num_new_nz: HighsInt,
        starts: &[HighsInt],
        indices: &[HighsInt],
        values: &[f64],
    ) -> HighsStatus {
        self.clear_derived_model_properties();
        let call_status = self.add_rows_interface(
            num_new_row,
            lower_bounds,
            upper_bounds,
            num_new_nz,
            starts,
            indices,
            values,
        );
        let return_status = crate::lp_data::highs_status::interpret_call_status(
            &self.options_.log_options,
            call_status,
            HighsStatus::kOk,
            "addRows",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::changeColsBounds (lines 3134-3169, 36 lines)
    /// Overload `changeColsBounds(num_set_entries, const HighsInt* set,
    /// const double* lower, const double* upper)`.
    ///
    /// `num_set_entries == 0` returns kOk first. The NULL checks are dropped.
    /// Local copies, `sort_set_data` (util/highs_sort_set.rs), `create_set`
    /// (util/highs_utils.rs:83); a nonzero create error logs one kError line
    /// ("Set supplied to Highs::changeColsBounds not ordered" style;
    /// `analyseSetCreateError` is not ported) and returns kError. Then
    /// `change_col_bounds_interface` (highs_interface_bounds.rs) through
    /// `interpret_call_status(.., "changeColBounds")` and `return_from_highs`.
    pub fn change_cols_bounds_set(
        &mut self,
        num_set_entries: HighsInt,
        set: &[HighsInt],
        lower: &[f64],
        upper: &[f64],
    ) -> HighsStatus {
        if num_set_entries == 0 {
            return HighsStatus::kOk;
        }
        self.clear_derived_model_properties();
        // Ensure that the set and data are in ascending order
        let n = num_set_entries as usize;
        let mut local_lower: Vec<f64> = lower[..n].to_vec();
        let mut local_upper: Vec<f64> = upper[..n].to_vec();
        let mut local_set: Vec<HighsInt> = set[..n].to_vec();
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
            self.model_lp_.num_col_,
        );
        if create_error != 0 {
            // analyseSetCreateError is not ported: one log line instead.
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Set supplied to Highs::changeColsBounds is illegal (create error {})\n",
                create_error
            );
            return HighsStatus::kError;
        }
        let call_status =
            self.change_col_bounds_interface(&mut index_collection, &local_lower, &local_upper);
        let return_status = crate::lp_data::highs_status::interpret_call_status(
            &self.options_.log_options,
            call_status,
            HighsStatus::kOk,
            "changeColBounds",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::deleteRows (lines 3568-3579, 12 lines)
    /// Overload `deleteRows(HighsInt* mask)`.
    ///
    /// `mask[i] != 0` deletes row `i`. On return `mask[i]` is the new index
    /// of row `i`, or -1 if it was deleted (copied back from
    /// `index_collection.mask_`). `create_mask` (util/highs_utils.rs:156)
    /// cannot fail here (`debug_assert!`).
    pub fn delete_rows_mask(&mut self, mask: &mut [HighsInt]) -> HighsStatus {
        self.clear_derived_model_properties();
        let original_num_row = self.model_lp_.num_row_;
        let mut index_collection = HighsIndexCollection::default();
        let create_error = create_mask(&mut index_collection, &mask[..], original_num_row);
        debug_assert!(create_error == 0);
        self.delete_rows_interface(&mut index_collection);
        let n = original_num_row as usize;
        mask[..n].copy_from_slice(&index_collection.mask_[..n]);
        self.return_from_highs(HighsStatus::kOk)
    }

    /// HiGHS: Highs.cpp::Highs::changeColsBounds (lines 3111-3132, 22 lines)
    /// Overload `changeColsBounds(from_col, to_col, lower, upper)` (interval).
    /// `create_interval` (util/highs_utils.rs:59); a nonzero create error logs
    /// the C++ kError line and returns kError.
    pub fn change_cols_bounds_interval(
        &mut self,
        from_col: HighsInt,
        to_col: HighsInt,
        lower: &[f64],
        upper: &[f64],
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
                "Interval [{}, {}] supplied to Highs::changeColsBounds is out of range [0, {})\n",
                from_col,
                to_col,
                self.model_lp_.num_col_
            );
            return HighsStatus::kError;
        }
        let call_status = self.change_col_bounds_interface(&mut index_collection, lower, upper);
        let return_status = crate::lp_data::highs_status::interpret_call_status(
            &self.options_.log_options,
            call_status,
            HighsStatus::kOk,
            "changeColBounds",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::deleteRows (lines 3537-3551, 15 lines)
    /// Overload `deleteRows(from_row, to_row)` (interval).
    pub fn delete_rows_interval(&mut self, from_row: HighsInt, to_row: HighsInt) -> HighsStatus {
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
                "Interval [{}, {}] supplied to Highs::deleteRows is out of range [0, {})\n",
                from_row,
                to_row,
                self.model_lp_.num_row_
            );
            return HighsStatus::kError;
        }
        self.delete_rows_interface(&mut index_collection);
        self.return_from_highs(HighsStatus::kOk)
    }
}
