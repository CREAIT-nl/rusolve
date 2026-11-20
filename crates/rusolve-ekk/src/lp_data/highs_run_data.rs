//! HiGHS source: lp_data/HighsRunData.h, lp_data/HighsRunData.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HighsRunDataStruct fields only; the records machinery (`getRunDataValue`,
//! `writeRunDataToFile`) is dropped, as for `HighsInfo`. The struct is the same in 1.15.1
//! and 1.16.0.
//!
//! PORT NOTE: the three time fields (`presolve_time`, `solve_time`, `postsolve_time`) are
//! not kept: the port's presolve and postsolve sub-clocks are not ported (highs_presolve.rs
//! reads 0.0 for every `this_*_time`), so they would only ever hold 0.

use crate::prelude::*;

/// HiGHS: HighsRunData.h::HighsRunDataStruct (without the time fields, module doc)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HighsRunData {
    pub valid: bool,
    pub presolved_model_num_col: HighsInt,
    pub presolved_model_num_row: HighsInt,
    pub presolved_model_num_nz: HighsInt,
    pub num_simplex_iterations_after_postsolve: HighsInt,
}

impl Default for HighsRunData {
    /// C++ `HighsRunData()` runs `initRecords`, whose records set every count to its default
    /// 0. `valid` has no record and no initialiser in C++; it starts false here.
    fn default() -> Self {
        Self {
            valid: false,
            presolved_model_num_col: 0,
            presolved_model_num_row: 0,
            presolved_model_num_nz: 0,
            num_simplex_iterations_after_postsolve: 0,
        }
    }
}

impl HighsRunData {
    /// HiGHS: HighsRunData.cpp::HighsRunData::invalidate (lines 17-26)
    pub fn invalidate(&mut self) {
        self.valid = false;
        self.presolved_model_num_col = kHighsIllegalIntMeasure;
        self.presolved_model_num_row = kHighsIllegalIntMeasure;
        self.presolved_model_num_nz = kHighsIllegalIntMeasure;
        self.num_simplex_iterations_after_postsolve = kHighsIllegalIntMeasure;
    }
}
