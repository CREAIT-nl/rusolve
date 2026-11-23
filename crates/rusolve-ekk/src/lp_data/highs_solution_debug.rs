//! HiGHS source: lp_data/HighsSolutionDebug.cpp (stubs; the subset `HEkk.cpp` calls)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Debug checks: out of scope.

use crate::prelude::*;

/// HiGHS: HighsSolutionDebug.cpp::debugHighsBasisConsistent
/// Debug check: out of scope; stub returns HighsDebugStatus::kNotChecked (the value HiGHS
/// returns at the default `highs_debug_level` = 0, i.e. below kHighsDebugLevelCheap).
pub fn debug_highs_basis_consistent(
    options: &HighsOptions,
    lp: &HighsLp,
    basis: &HighsBasis,
) -> HighsDebugStatus {
    HighsDebugStatus::kNotChecked
}
