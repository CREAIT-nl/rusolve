//! HiGHS source: simplex/HSimplexNlaDebug.cpp (stubs)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Debug checks: out of scope.

use super::*;
use crate::prelude::*;

impl HSimplexNla {
    /// HiGHS: HSimplexNlaDebug.cpp::HSimplexNla::debugCheckInvert (lines 23-250, 228 lines)
    /// C++ default arguments: alt_debug_level = -1
    /// Debug check: out of scope; stub returns HighsDebugStatus::kNotChecked (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_check_invert(
        &self,
        refs: &NlaRefs,
        message: &str,
        alt_debug_level: HighsInt,
    ) -> HighsDebugStatus {
        HighsDebugStatus::kNotChecked
    }

    /// HiGHS: HSimplexNlaDebug.cpp::HSimplexNla::debugInvertResidualError (lines 252-301, 50 lines)
    /// Debug check: out of scope; stub returns 0.0 (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_invert_residual_error(
        &self,
        refs: &NlaRefs,
        transposed: bool,
        solution: &HVector,
        residual: &mut HVector,
    ) -> f64 {
        0.0
    }

    /// HiGHS: HSimplexNlaDebug.cpp::HSimplexNla::debugReportInvertSolutionError (lines 303-318, 16 lines)
    /// Debug check: out of scope; stub returns HighsDebugStatus::kNotChecked (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_report_invert_solution_error(
        &self,
        transposed: bool,
        true_solution: &HVector,
        solution: &HVector,
        residual: &mut HVector,
        force: bool,
    ) -> HighsDebugStatus {
        HighsDebugStatus::kNotChecked
    }

    /// HiGHS: HSimplexNlaDebug.cpp::HSimplexNla::debugReportInvertSolutionError (lines 320-373, 54 lines)
    /// Debug check: out of scope; stub returns HighsDebugStatus::kNotChecked (the value HiGHS returns at the default
    /// `highs_debug_level` = 0 for the callers on the simplex path).
    pub fn debug_report_invert_solution_error_norms(
        &self,
        source: &str,
        transposed: bool,
        solve_error_norm: f64,
        residual_error_norm: f64,
        force: bool,
    ) -> HighsDebugStatus {
        HighsDebugStatus::kNotChecked
    }
}
