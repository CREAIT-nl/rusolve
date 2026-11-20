//! HiGHS source: mip/HighsPrimalHeuristics.h (class definitions and header-inline methods)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

mod part1;
pub use part1::*;
mod part2;
pub use part2::*;

/// HiGHS: HighsPrimalHeuristics.h::HighsPrimalHeuristics (lines 21-90)
///
/// Removed C++ members:
/// * C++: `const HighsMipSolver & mipsolver`: removed; methods that use it take `ms: &mut HighsMipSolver`.
#[derive(Clone, Default, Debug)]
pub struct HighsPrimalHeuristics {
    /// C++: `std::vector<HighsInt> intcols`
    pub intcols: Vec<HighsInt>,
    /// C++: `double successObservations`
    pub success_observations: f64,
    /// C++: `HighsInt numSuccessObservations`
    pub num_success_observations: HighsInt,
    /// C++: `double infeasObservations`
    pub infeas_observations: f64,
    /// C++: `HighsInt numInfeasObservations`
    pub num_infeas_observations: HighsInt,
    /// C++: `HighsRandom randgen`
    pub randgen: HighsRandom,
    /// rusolve `mip_guided_dive`: LP iterations spent by the guided dive so far
    pub guided_lp_iterations: i64,
}

/// HiGHS: HighsPrimalHeuristics.cpp::HeuristicNeighbourhood (lines 276-312)
///
/// Removed C++ members:
/// * C++: `HighsDomain & localdom`: removed; `get_fixing_rate` takes `ms` and `localdom: &HighsDomain`.
#[derive(Clone, Default, Debug)]
pub struct HeuristicNeighbourhood {
    /// C++: `HighsInt numFixed`
    pub num_fixed: HighsInt,
    /// C++: `HighsHashTable<HighsInt> fixedCols`
    pub fixed_cols: HighsHashTable<HighsInt, ()>,
    /// C++: `size_t startCheckedChanges`
    pub start_checked_changes: usize,
    /// C++: `size_t nCheckedChanges`
    pub n_checked_changes: usize,
    /// C++: `HighsInt numTotal`
    pub num_total: HighsInt,
}

impl HighsPrimalHeuristics {
    /// HiGHS: HighsPrimalHeuristics.h::HighsPrimalHeuristics::getHeuristicRandom (lines 87-89)
    /// C++: `HighsInt getHeuristicRandom(const HighsInt sup)`
    pub fn get_heuristic_random(&mut self, sup: HighsInt) -> HighsInt {
        self.randgen.integer_sup(sup)
    }
}
