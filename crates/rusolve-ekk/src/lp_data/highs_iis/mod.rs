//! HiGHS source: lp_data/HighsIis.h, lp_data/HighsIis.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HighsIis.cpp (1290 lines) is split into contiguous-range part files:
//! `part1.rs` (14-366), `part2.rs` (368-591), `part3.rs` (593-979),
//! `part4.rs` (981-1290).

pub mod part1;
pub mod part2;
pub mod part3;
pub mod part4;

use crate::prelude::*;

/// HiGHS: HighsIis.h::kIisDevReport (line 16)
pub const kIisDevReport: bool = false;

/// HiGHS: HighsIis.h::IisBoundStatus (lines 18-25)
pub type IisBoundStatus = HighsInt;
pub const kIisBoundStatusDropped: IisBoundStatus = -1;
pub const kIisBoundStatusNull: IisBoundStatus = 0;
pub const kIisBoundStatusFree: IisBoundStatus = 1;
pub const kIisBoundStatusLower: IisBoundStatus = 2;
pub const kIisBoundStatusUpper: IisBoundStatus = 3;
pub const kIisBoundStatusBoxed: IisBoundStatus = 4;

/// HiGHS: HighsIis.h::IisModelStatus (lines 27-33)
pub type IisModelStatus = HighsInt;
pub const kIisModelStatusFeasible: IisModelStatus = -1;
pub const kIisModelStatusUnknown: IisModelStatus = 0;
pub const kIisModelStatusTimeLimit: IisModelStatus = 1;
pub const kIisModelStatusReducible: IisModelStatus = 2;
pub const kIisModelStatusIrreducible: IisModelStatus = 3;

/// HiGHS: HighsIis.h::HighsIisInfo (lines 35-75)
#[derive(Clone, Debug, PartialEq)]
pub struct HighsIisInfo {
    pub num_lp_solved: HighsInt,
    pub sum_simplex_iteration_counts: HighsInt,
    pub min_simplex_iteration_count: HighsInt,
    pub max_simplex_iteration_count: HighsInt,
    pub sum_simplex_times: f64,
    pub min_simplex_time: f64,
    pub max_simplex_time: f64,
    // IIS logging state (persist across calls)
    pub iis_last_disptime: f64,
    pub iis_num_disp_lines: HighsInt,
}

impl Default for HighsIisInfo {
    fn default() -> Self {
        Self {
            num_lp_solved: 0,
            sum_simplex_iteration_counts: 0,
            min_simplex_iteration_count: 0,
            max_simplex_iteration_count: 0,
            sum_simplex_times: 0.0,
            min_simplex_time: 0.0,
            max_simplex_time: 0.0,
            iis_last_disptime: -kHighsInf,
            iis_num_disp_lines: 0,
        }
    }
}

impl HighsIisInfo {
    /// HiGHS: HighsIis.h::HighsIisInfo::clear (lines 47-58)
    pub fn clear(&mut self) {
        self.num_lp_solved = 0;
        self.sum_simplex_iteration_counts = 0;
        self.min_simplex_iteration_count = 0;
        self.max_simplex_iteration_count = 0;
        self.sum_simplex_times = 0.0;
        self.min_simplex_time = 0.0;
        self.max_simplex_time = 0.0;
        // Reset IIS logging state
        self.iis_last_disptime = -kHighsInf;
        self.iis_num_disp_lines = 0;
    }

    /// HiGHS: HighsIis.h::HighsIisInfo::update (lines 60-74)
    ///
    /// PORT NOTE: `std::min`/`std::max` -> `cpp_min`/`cpp_max`.
    pub fn update(&mut self, simplex_time: f64, simplex_iterations: HighsInt) {
        if self.num_lp_solved == 0 {
            self.min_simplex_iteration_count = kHighsIInf;
            self.min_simplex_time = kHighsInf;
        }
        self.num_lp_solved += 1;
        self.sum_simplex_times += simplex_time;
        self.min_simplex_time = cpp_min(simplex_time, self.min_simplex_time);
        self.max_simplex_time = cpp_max(simplex_time, self.max_simplex_time);
        self.sum_simplex_iteration_counts += simplex_iterations;
        self.min_simplex_iteration_count =
            cpp_min(simplex_iterations, self.min_simplex_iteration_count);
        self.max_simplex_iteration_count =
            cpp_max(simplex_iterations, self.max_simplex_iteration_count);
    }
}

/// HiGHS: HighsIis.h::HighsIis (lines 77-137)
///
/// PORT NOTE: C++ leaves `status_` uninitialised on construction; it is
/// `kIisModelStatusUnknown` (0) here, which is what `clearData` sets.
#[derive(Clone, Debug)]
pub struct HighsIis {
    pub valid_: bool,
    pub status_: HighsInt,
    pub strategy_: HighsInt,
    pub col_index_: Vec<HighsInt>,
    pub row_index_: Vec<HighsInt>,
    pub col_bound_: Vec<HighsInt>,
    pub row_bound_: Vec<HighsInt>,
    pub col_status_: Vec<HighsInt>,
    pub row_status_: Vec<HighsInt>,
    pub info_: HighsIisInfo,
    pub model_: HighsModel,
}

impl Default for HighsIis {
    fn default() -> Self {
        Self {
            valid_: false,
            status_: kIisModelStatusUnknown,
            strategy_: kIisStrategyMin,
            col_index_: Vec::new(),
            row_index_: Vec::new(),
            col_bound_: Vec::new(),
            row_bound_: Vec::new(),
            col_status_: Vec::new(),
            row_status_: Vec::new(),
            info_: HighsIisInfo::default(),
            model_: HighsModel::default(),
        }
    }
}

impl HighsIis {
    /// HiGHS: HighsIis.h::HighsIis::indexStatusOkReturn (lines 115-117)
    pub fn index_status_ok_return(&self, return_value: bool) -> bool {
        return_value
    }

    /// HiGHS: HighsIis.h::HighsIis::lpDataOkReturn (line 118)
    pub fn lp_data_ok_return(&self, return_value: bool) -> bool {
        return_value
    }

    /// HiGHS: HighsIis.h::HighsIis::lpOkReturn (line 119)
    pub fn lp_ok_return(&self, return_value: bool) -> bool {
        return_value
    }
}
