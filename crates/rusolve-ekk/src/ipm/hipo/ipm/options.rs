//! HiGHS source: ipm/hipo/ipm/Options.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiPO options (complete). `log_options` is a clone of the
//! `HighsLogOptions` the C++ points at (`None` = nullptr).

use super::parameters::{kIpmTolDefault, kMaxIterDefault};
use crate::io::highs_io::HighsLogOptions;
use crate::ipm::hipo::Int;
use crate::lp_data::hconst::{kHighsChooseString, kHighsOffString};
use crate::lp_data::highs_options::kHipoBothString;

/// HiGHS: Options.h `struct Options`
#[derive(Clone)]
pub struct Options {
    /// rusolve: wall-clock deadline of the HiGHS time limit, polled
    /// by the factorisation between supernodes. `None` without a finite limit.
    pub deadline_: Option<std::time::Instant>,
    pub nla: String,
    pub crossover: String,
    pub parallel: String,
    pub parallel_type: String,
    pub ordering: String,
    pub max_iter: Int,
    pub feasibility_tol: f64,
    pub optimality_tol: f64,
    pub crossover_tol: f64,
    pub refine_with_ipx: bool,
    /// rusolve: give up at analyse when the factorisation cannot
    /// finish ~20 iterations inside the time limit
    pub flops_bailout: bool,
    /// rusolve: stop when dinf stays flat while mu falls by 1e4
    pub stall_detect: bool,
    /// rusolve: AMD first; Metis only when AMD flops >= 1e9
    pub amd_first: bool,
    pub time_limit: f64,
    pub block_size: Int,
    pub display: bool,
    pub timeless_log: bool,
    pub log_options: Option<HighsLogOptions>,
    // rusolve (HiGHS 1.16 backport, default off): inert flags, default false; set in
    // Solver::setOptions (ipm/hipo/ipm/solver/mod.rs). h16_ipx_timeout is not ported.
    pub h16_ipx_stall: bool,
    pub h16_ipx_scaled: bool,
    pub h16_ipm_gap_offset: bool,
    pub h16_ipx_free_rows: bool,
    pub h16_hipo_empty_rows: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            deadline_: None,
            nla: kHighsChooseString.to_string(),
            crossover: kHighsOffString.to_string(),
            parallel: kHighsChooseString.to_string(),
            parallel_type: kHipoBothString.to_string(),
            ordering: kHighsChooseString.to_string(),
            max_iter: kMaxIterDefault,
            feasibility_tol: kIpmTolDefault,
            optimality_tol: kIpmTolDefault,
            crossover_tol: kIpmTolDefault,
            refine_with_ipx: true,
            flops_bailout: false,
            stall_detect: false,
            amd_first: false,
            time_limit: -1.0,
            block_size: 0,
            display: true,
            timeless_log: false,
            log_options: None,
            h16_ipx_stall: false,
            h16_ipx_scaled: false,
            h16_ipm_gap_offset: false,
            h16_ipx_free_rows: false,
            h16_hipo_empty_rows: false,
        }
    }
}
