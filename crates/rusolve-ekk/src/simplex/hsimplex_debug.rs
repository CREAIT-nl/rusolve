//! HiGHS source: simplex/HSimplexDebug.cpp (stubs; the subset `HEkkDualRow.cpp` calls)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! CHUZC failure diagnostics: out of scope. HEkkDualRow ignores their return
//! value, so the stubs only keep the call surface.

use crate::prelude::*;

/// HiGHS: HSimplexDebug.cpp::debugDualChuzcFailQuad0
/// C++ default arguments: force = false
/// Debug/report: out of scope; stub (HighsDebugStatus::kNotChecked).
pub fn debug_dual_chuzc_fail_quad0(
    options: &HighsOptions,
    work_count: HighsInt,
    work_data: &[(HighsInt, f64)],
    num_var: HighsInt,
    work_dual: &[f64],
    select_theta: f64,
    remain_theta: f64,
    force: bool,
) -> HighsDebugStatus {
    HighsDebugStatus::kNotChecked
}

/// HiGHS: HSimplexDebug.cpp::debugDualChuzcFailQuad1
/// C++ default arguments: force = false
/// Debug/report: out of scope; stub (HighsDebugStatus::kNotChecked).
pub fn debug_dual_chuzc_fail_quad1(
    options: &HighsOptions,
    work_count: HighsInt,
    work_data: &[(HighsInt, f64)],
    num_var: HighsInt,
    work_dual: &[f64],
    select_theta: f64,
    force: bool,
) -> HighsDebugStatus {
    HighsDebugStatus::kNotChecked
}

/// HiGHS: HSimplexDebug.cpp::debugDualChuzcFailHeap
/// C++ default arguments: force = false
/// Debug/report: out of scope; stub (HighsDebugStatus::kNotChecked).
pub fn debug_dual_chuzc_fail_heap(
    options: &HighsOptions,
    work_count: HighsInt,
    work_data: &[(HighsInt, f64)],
    num_var: HighsInt,
    work_dual: &[f64],
    select_theta: f64,
    force: bool,
) -> HighsDebugStatus {
    HighsDebugStatus::kNotChecked
}
