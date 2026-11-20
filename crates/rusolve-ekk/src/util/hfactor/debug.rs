//! HiGHS source: util/HFactorDebug.cpp (stubs)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Rank-deficiency reporting called from `HFactor::build*`: out of scope. Every
//! function is a no-op (HiGHS itself returns early at the default
//! `highs_debug_level` = 0, apart from log output).

use crate::prelude::*;

/// HiGHS: HFactorDebug.cpp::debugReportRankDeficiency
/// Debug/report: out of scope; stub (no-op).
pub fn debug_report_rank_deficiency(
    call_id: HighsInt,
    highs_debug_level: HighsInt,
    log_options: &HighsLogOptions,
    num_row: HighsInt,
    permute: &[HighsInt],
    iwork: &[HighsInt],
    basic_index: &[HighsInt],
    rank_deficiency: HighsInt,
    row_with_no_pivot: &[HighsInt],
    col_with_no_pivot: &[HighsInt],
) {
}

/// HiGHS: HFactorDebug.cpp::debugReportRankDeficientASM
/// Debug/report: out of scope; stub (no-op).
pub fn debug_report_rank_deficient_asm(
    highs_debug_level: HighsInt,
    log_options: &HighsLogOptions,
    num_row: HighsInt,
    mc_start: &[HighsInt],
    mc_count_a: &[HighsInt],
    mc_index: &[HighsInt],
    mc_value: &[f64],
    iwork: &[HighsInt],
    rank_deficiency: HighsInt,
    col_with_no_pivot: &[HighsInt],
    row_with_no_pivot: &[HighsInt],
) {
}

/// HiGHS: HFactorDebug.cpp::debugReportMarkSingC
/// Debug/report: out of scope; stub (no-op).
pub fn debug_report_mark_sing_c(
    call_id: HighsInt,
    highs_debug_level: HighsInt,
    log_options: &HighsLogOptions,
    num_row: HighsInt,
    iwork: &[HighsInt],
    basic_index: &[HighsInt],
) {
}

/// HiGHS: HFactorDebug.cpp::debugLogRankDeficiency
/// Debug/report: out of scope; stub (no-op).
pub fn debug_log_rank_deficiency(
    highs_debug_level: HighsInt,
    log_options: &HighsLogOptions,
    rank_deficiency: HighsInt,
    basis_matrix_num_el: HighsInt,
    invert_num_el: HighsInt,
    kernel_dim: HighsInt,
    kernel_num_el: HighsInt,
    nwork: HighsInt,
) {
}

/// HiGHS: HFactorDebug.cpp::debugPivotValueAnalysis
/// Debug/report: out of scope; stub (no-op).
pub fn debug_pivot_value_analysis(
    highs_debug_level: HighsInt,
    log_options: &HighsLogOptions,
    num_row: HighsInt,
    u_pivot_value: &[f64],
) {
}
