//! HiGHS source: mip/HighsSeparator.h, mip/HighsSeparator.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsSeparator.h (line 22)
/// C++: `const std::string kTableauSepaString`
pub const kTableauSepaString: &str = "Separation: Tableau";
/// HiGHS: HighsSeparator.h (line 23)
/// C++: `const std::string kPathAggrSepaString`
pub const kPathAggrSepaString: &str = "Separation: Path aggregation";
/// HiGHS: HighsSeparator.h (line 24)
/// C++: `const std::string kModKSepaString`
pub const kModKSepaString: &str = "Separation: Mod-k";
/// HiGHS latest (#2990, dc83470a9): HighsSeparator.h (lines 25-26)
/// C++: `const std::string kMachineSchedSepaString`
pub const kMachineSchedSepaString: &str = "Separation: Single Machine Scheduling";

/// HiGHS: HighsSeparator.h::HighsSeparator (lines 34-58)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsSeparator {
    /// C++: `HighsInt numCutsFound`
    pub num_cuts_found: HighsInt,
    /// C++: `HighsInt numCalls`
    pub num_calls: HighsInt,
    /// C++: `int clockIndex`
    pub clock_index: i32,
}

impl HighsSeparator {
    /// HiGHS: HighsSeparator.cpp::HighsSeparator::HighsSeparator (lines 17-30)
    /// C++: `HighsSeparator(const HighsMipSolver & mipsolver, const std::string & name)`
    // PORT NOTE: the C++ body that reads `mipsolver.analysis_` is commented out in HiGHS
    // 1.15.1; `clockIndex` is the constant 999.
    pub fn new(_mipsolver: &HighsMipSolver, _name: &str) -> Self {
        HighsSeparator {
            num_cuts_found: 0,
            num_calls: 0,
            clock_index: 999,
        }
    }

    /// HiGHS: HighsSeparator.h::HighsSeparator::getNumCutsFound (lines 51-51)
    /// C++: `HighsInt getNumCutsFound() const`
    pub fn get_num_cuts_found(&self) -> HighsInt {
        self.num_cuts_found
    }

    /// HiGHS: HighsSeparator.h::HighsSeparator::getNumCalls (lines 53-53)
    /// C++: `HighsInt getNumCalls() const`
    pub fn get_num_calls(&self) -> HighsInt {
        self.num_calls
    }

    /// HiGHS: HighsSeparator.h::HighsSeparator::getClockIndex (lines 55-55)
    /// C++: `HighsInt getClockIndex() const`
    pub fn get_clock_index(&self) -> HighsInt {
        self.clock_index as HighsInt
    }

    /// HiGHS: HighsSeparator.h::HighsSeparator::~HighsSeparator (lines 57-57)
    /// C++: `~HighsSeparator()`
    pub fn release(&mut self) {}
}

/// The virtual interface of `HighsSeparator` (C++ `separateLpSolution` is pure virtual and
/// `run` is the non-virtual driver).
///
/// Ownership, shared by every separator: `lp_relaxation` is the LP detached by
/// `HighsSeparation` (`SepaCtx::lp`); `cutpool` is an index into `ms.mipdata_.cutpools`
/// (attached; a callee that adds a cut detaches it for that call); `trans_lp` /
/// `lp_aggregator` are the separation-round locals and hold no references.
impl HighsSeparatorBox {
    /// index into `mip_stats::SEPARATOR_NAMES`.
    pub fn kind_index(&self) -> usize {
        match self {
            HighsSeparatorBox::Tableau(_) => 0,
            HighsSeparatorBox::Path(_) => 1,
            HighsSeparatorBox::Modk(_) => 2,
            HighsSeparatorBox::MachineSched(_) => 3,
            HighsSeparatorBox::KnapsackCover(_) => 4,
            HighsSeparatorBox::TiPrecedence(_) => 5,
            HighsSeparatorBox::Netdesign(_) => 6,
            HighsSeparatorBox::Dff(_) => 7,
        }
    }

    /// The C++ base-class subobject.
    pub fn base(&self) -> &HighsSeparator {
        match self {
            HighsSeparatorBox::Tableau(s) => &s.base,
            HighsSeparatorBox::Path(s) => &s.base,
            HighsSeparatorBox::Modk(s) => &s.base,
            HighsSeparatorBox::MachineSched(s) => &s.base,
            HighsSeparatorBox::KnapsackCover(s) => &s.base,
            HighsSeparatorBox::TiPrecedence(s) => &s.base,
            HighsSeparatorBox::Netdesign(s) => &s.base,
            HighsSeparatorBox::Dff(s) => &s.base,
        }
    }

    /// The C++ base-class subobject.
    pub fn base_mut(&mut self) -> &mut HighsSeparator {
        match self {
            HighsSeparatorBox::Tableau(s) => &mut s.base,
            HighsSeparatorBox::Path(s) => &mut s.base,
            HighsSeparatorBox::Modk(s) => &mut s.base,
            HighsSeparatorBox::MachineSched(s) => &mut s.base,
            HighsSeparatorBox::KnapsackCover(s) => &mut s.base,
            HighsSeparatorBox::TiPrecedence(s) => &mut s.base,
            HighsSeparatorBox::Netdesign(s) => &mut s.base,
            HighsSeparatorBox::Dff(s) => &mut s.base,
        }
    }

    /// C++ `virtual void separateLpSolution(HighsLpRelaxation&, HighsLpAggregator&,
    /// HighsTransformedLp&, HighsCutPool&) = 0` (dispatch).
    pub fn separate_lp_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &mut HighsLpRelaxation,
        lp_aggregator: &mut HighsLpAggregator,
        trans_lp: &mut HighsTransformedLp,
        cutpool: usize,
    ) {
        match self {
            HighsSeparatorBox::Tableau(s) => {
                s.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool)
            }
            HighsSeparatorBox::Path(s) => {
                s.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool)
            }
            HighsSeparatorBox::Modk(s) => {
                s.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool)
            }
            HighsSeparatorBox::MachineSched(s) => {
                s.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool)
            }
            HighsSeparatorBox::KnapsackCover(s) => {
                s.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool)
            }
            HighsSeparatorBox::TiPrecedence(s) => {
                s.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool)
            }
            HighsSeparatorBox::Netdesign(s) => {
                s.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool)
            }
            HighsSeparatorBox::Dff(s) => {
                s.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool)
            }
        }
    }

    /// HiGHS: HighsSeparator.cpp::HighsSeparator::run (lines 32-45)
    /// C++: `void run(HighsLpRelaxation & lpRelaxation, HighsLpAggregator & lpAggregator, HighsTransformedLp & transLp, HighsCutPool & cutpool)`
    // PORT NOTE: on the box (the C++ base method calls the virtual `separateLpSolution`);
    // `lpRelaxation.getMipSolver()` is `ms`; `cutpool` is an index (see the impl doc).
    pub fn run(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &mut HighsLpRelaxation,
        lp_aggregator: &mut HighsLpAggregator,
        trans_lp: &mut HighsTransformedLp,
        cutpool: usize,
    ) {
        self.base_mut().num_calls += 1;
        let curr_num_cuts = ms.mipdata_.cutpools[cutpool].get_num_cuts();

        let clock_index = self.base().get_clock_index();
        if !ms.parallel_lock_active() {
            ms.profiling_.start(clock_index, false);
        }
        let c5_start = std::time::Instant::now(); // reporting-only
        self.separate_lp_solution(ms, lp_relaxation, lp_aggregator, trans_lp, cutpool);
        let c5_kind = self.kind_index();
        ms.stats_.sepa_ns[c5_kind].add(c5_start.elapsed().as_nanos() as i64);
        if !ms.parallel_lock_active() {
            ms.profiling_.stop(clock_index);
        }

        self.base_mut().num_cuts_found +=
            ms.mipdata_.cutpools[cutpool].get_num_cuts() - curr_num_cuts;
        ms.stats_.sepa_calls[c5_kind].add(1);
        ms.stats_.sepa_cuts[c5_kind]
            .add((ms.mipdata_.cutpools[cutpool].get_num_cuts() - curr_num_cuts) as i64);
    }
}
