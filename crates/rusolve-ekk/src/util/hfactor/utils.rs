//! HiGHS source: util/HFactorUtils.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Reporting helpers are out of scope (no-ops); `invalidAMatrixAction` is kept.

use super::*;
use crate::prelude::*;

impl HFactor {
    /// HiGHS: HFactorUtils.cpp::HFactor::invalidAMatrixAction (lines 13-16, 4 lines)
    pub fn invalid_a_matrix_action(&mut self) {
        self.a_matrix_valid = false;
        self.refactor_info_.clear();
    }

    /// HiGHS: HFactorUtils.cpp::HFactor::reportLu (lines 18-74, 57 lines)
    /// C++ default arguments: l_u_or_both = kReportLuBoth; full = true
    /// Reporting: out of scope, no-op.
    pub fn report_lu(&self, l_u_or_both: HighsInt, full: bool) {}

    /// HiGHS: HFactorUtils.cpp::HFactor::reportIntVector (lines 76-87, 12 lines)
    /// Reporting: out of scope, no-op.
    pub fn report_int_vector(&self, name: &str, entry: &[HighsInt]) {}

    /// HiGHS: HFactorUtils.cpp::HFactor::reportDoubleVector (lines 88-99, 12 lines)
    /// Reporting: out of scope, no-op.
    pub fn report_double_vector(&self, name: &str, entry: &[f64]) {}

    /// HiGHS: HFactorUtils.cpp::HFactor::reportAsm (lines 101-122, 22 lines)
    /// Reporting: out of scope, no-op.
    pub fn report_asm(&self) {}
}
