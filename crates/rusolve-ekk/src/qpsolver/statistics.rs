//! HiGHS source: qpsolver/statistics.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete: data only. `std::chrono::high_resolution_clock::time_point`
//! becomes `Option<std::time::Instant>` (`None` = the C++ default epoch value).

use crate::prelude::*;

/// HiGHS: statistics.hpp::Statistics (lines 14-28)
#[derive(Clone, Debug, Default)]
pub struct Statistics {
    pub phase1_iterations: HighsInt,
    pub num_iterations: HighsInt,
    pub time_start: Option<std::time::Instant>,
    pub time_end: Option<std::time::Instant>,

    pub iteration: Vec<HighsInt>,
    pub nullspacedimension: Vec<HighsInt>,
    pub objval: Vec<f64>,
    pub time: Vec<f64>,
    pub sum_primal_infeasibilities: Vec<f64>,
    pub num_primal_infeasibilities: Vec<HighsInt>,
    pub density_nullspace: Vec<f64>,
    pub density_factor: Vec<f64>,
}
