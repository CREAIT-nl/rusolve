//! HiGHS source: lp_data/Highs.cpp (clearDerivedModelProperties,
//! invalidateModelStatusSolutionAndInfo, invalidateModelStatusAndInfo,
//! invalidateModelStatus, newHighsBasis, forceHighsSolutionBasisSize,
//! returnFromHighs), lp_data/Highs.h (clearRayRecords)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! API slice: the bookkeeping every public
//! model-modification method of `Highs` runs before/after its interface
//! call. Ranging, IIS and the standard-form LP are not members of
//! this port, so their invalidations / clears are dropped (each listed at
//! its function).

use crate::lp_data::highs_run::Highs;
use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::clearDerivedModelProperties (lines 3897-3901, 5 lines)
    ///
    /// `clearPresolve()` is `self.clear_presolve()` (highs_presolve.rs:81);
    /// `clearStandardFormLp()` is [`Self::clear_standard_form_lp`];
    /// `clearRayRecords()` is [`Self::clear_ray_records`].
    pub fn clear_derived_model_properties(&mut self) {
        self.clear_presolve();
        self.clear_standard_form_lp();
        self.clear_ray_records();
    }

    /// HiGHS: Highs.cpp::Highs::clearStandardFormLp (lines 3909-3915, 7 lines)
    pub fn clear_standard_form_lp(&mut self) {
        self.standard_form_valid_ = false;
        self.standard_form_offset_ = 0.0;
        self.standard_form_cost_.clear();
        self.standard_form_rhs_.clear();
        self.standard_form_matrix_.clear();
    }

    /// HiGHS: Highs.h::Highs::clearRayRecords (line 1689)
    /// -> `this->ekk_instance_.clearRayRecords()` (simplex/ekk/clear.rs:54).
    pub fn clear_ray_records(&mut self) {
        self.ekk_instance_.clear_ray_records();
    }

    /// HiGHS: Highs.cpp::Highs::invalidateModelStatusSolutionAndInfo (lines 3932-3935, 4 lines)
    pub(crate) fn invalidate_model_status_solution_and_info(&mut self) {
        self.invalidate_model_status_and_info();
        self.invalidate_solution();
    }

    /// HiGHS: Highs.cpp::Highs::invalidateModelStatusAndInfo (lines 3937-3943, 7 lines)
    pub(crate) fn invalidate_model_status_and_info(&mut self) {
        self.invalidate_model_status();
        self.invalidate_ranging();
        self.info_.invalidate();
        self.run_data_.invalidate();
        self.clear_iis();
    }

    /// HiGHS: Highs.cpp::Highs::invalidateRanging (line 3970)
    pub(crate) fn invalidate_ranging(&mut self) {
        self.ranging_.invalidate();
    }

    /// HiGHS: Highs.cpp::Highs::clearIis (line 3974)
    pub(crate) fn clear_iis(&mut self) {
        self.iis_.clear();
    }

    /// HiGHS: Highs.cpp::Highs::invalidateModelStatus (lines 3945-3947, 3 lines)
    pub(crate) fn invalidate_model_status(&mut self) {
        self.model_status_ = HighsModelStatus::kNotset;
    }

    /// HiGHS: Highs.cpp::Highs::newHighsBasis (lines 4695-4698, 4 lines)
    /// -> `ekk_instance_.update_status(LpAction::kNewBasis)`.
    pub(crate) fn new_highs_basis(&mut self) {
        self.ekk_instance_.update_status(LpAction::kNewBasis);
    }

    /// HiGHS: Highs.cpp::Highs::forceHighsSolutionBasisSize (lines 4703-4737, 35 lines)
    ///
    /// `resize(n, 0)` / `resize(n, kNonbasic)` / `resize(n, kBasic)` are
    /// `Vec::resize` with the same fill values.
    pub(crate) fn force_highs_solution_basis_size(&mut self) {
        let num_col = self.model_lp_.num_col_;
        let num_row = self.model_lp_.num_row_;
        // Values
        if self.solution_.col_value.len() < num_col as usize
            || self.solution_.row_value.len() < num_row as usize
        {
            self.solution_.value_valid = false;
            self.info_.primal_solution_status = kSolutionStatusNone;
        }
        self.solution_.col_value.resize(num_col as usize, 0.0);
        self.solution_.row_value.resize(num_row as usize, 0.0);
        // Duals
        if self.solution_.col_dual.len() < num_col as usize
            || self.solution_.row_dual.len() < num_row as usize
        {
            self.solution_.dual_valid = false;
            self.info_.dual_solution_status = kSolutionStatusNone;
        }
        self.solution_.col_dual.resize(num_col as usize, 0.0);
        self.solution_.row_dual.resize(num_row as usize, 0.0);
        // Basis
        if self.basis_.col_status.len() != num_col as usize
            || self.basis_.row_status.len() != num_row as usize
        {
            self.basis_.valid = false;
            self.basis_.useful = false;
            self.info_.basis_validity = kBasisValidityInvalid;
        }
        self.basis_
            .col_status
            .resize(num_col as usize, HighsBasisStatus::kNonbasic);
        self.basis_
            .row_status
            .resize(num_row as usize, HighsBasisStatus::kBasic);
    }

    /// HiGHS: Highs.cpp::Highs::returnFromHighs (lines 4963-5017, 55 lines)
    ///
    /// * `debugHighsBasisConsistent` is
    ///   `crate::lp_data::highs_solution_debug::debug_highs_basis_consistent`;
    ///   `ekk_instance_.debugRetainedDataOk` is `debug_retained_data_ok`
    ///   (simplex/ekk/debug.rs:61). `assert` -> `debug_assert!`.
    /// * The `called_return_from_optimize_model` check is dropped (no such
    ///   member in this port; it only logs).
    /// * `if (timer_.running()) timer_.stop();`
    /// * `lpDimensionsOk` is `crate::lp_data::lp_assess::lp_dimensions_ok`.
    /// * `lpFactorRowCompatible(n)` is `lp_factor_row_compatible_num(n)`
    ///   (simplex/ekk/factor.rs:682); on a mismatch `ekk_instance_.clear()`.
    pub(crate) fn return_from_highs(&mut self, highs_return_status: HighsStatus) -> HighsStatus {
        // Applies checks before returning from HiGHS
        let mut return_status = highs_return_status;

        self.force_highs_solution_basis_size();

        let consistent =
            debug_highs_basis_consistent(&self.options_, &self.model_lp_, &self.basis_)
                != HighsDebugStatus::kLogicalError;
        if !consistent {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "returnFromHighs: Supposed to be a HiGHS basis, but not consistent\n"
            );
            debug_assert!(consistent);
            return_status = HighsStatus::kError;
        }
        // Check that any retained Ekk data - basis and NLA - are OK
        let retained_ekk_data_ok = self.ekk_instance_.debug_retained_data_ok(&self.model_lp_)
            != HighsDebugStatus::kLogicalError;
        if !retained_ekk_data_ok {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "returnFromHighs: Retained Ekk data not OK\n"
            );
            debug_assert!(retained_ekk_data_ok);
            return_status = HighsStatus::kError;
        }
        // Stop the HiGHS run clock if it is running
        if self.timer_.running() {
            self.timer_.stop();
        }
        let dimensions_ok = crate::lp_data::lp_assess::lp_dimensions_ok(
            "returnFromHighs",
            &self.model_lp_,
            &self.options_.log_options,
        );
        if !dimensions_ok {
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kError,
                "LP Dimension error in returnFromHighs()\n"
            );
            return_status = HighsStatus::kError;
        }
        debug_assert!(dimensions_ok);
        if self.ekk_instance_.status_.has_nla {
            let lp_factor_row_compatible = self
                .ekk_instance_
                .lp_factor_row_compatible_num(self.model_lp_.num_row_);
            if !lp_factor_row_compatible {
                crate::highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kWarning,
                    "Highs::returnFromHighs(): LP and HFactor have inconsistent numbers of rows\n"
                );
                debug_assert!(lp_factor_row_compatible);
                // Clear Ekk entirely
                self.ekk_instance_.clear();
            }
        }
        return_status
    }
}
