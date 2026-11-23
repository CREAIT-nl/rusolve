//! HiGHS source: mip/MipTimer.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: MipTimer.h::MipTimer (lines 213-635)
#[derive(Clone, Default, Debug)]
pub struct MipTimer {}

/// HiGHS: MipTimer.h::iClockMip (plain C enum: module consts).
pub type iClockMip = i32;
pub const kMipClockTotal: i32 = 0;
pub const kMipClockPresolve: i32 = 15;
pub const kMipClockSolve: i32 = 16;
pub const kMipClockPostsolve: i32 = 17;
pub const kFromMipClock: i32 = 18;
pub const kMipClockInit: i32 = 18;
pub const kMipClockRunPresolve: i32 = 19;
pub const kMipClockRunSetup: i32 = 20;
pub const kMipClockFeasibilityJump: i32 = 21;
pub const kMipClockTrivialHeuristics: i32 = 22;
pub const kMipClockEvaluateRootNode: i32 = 23;
pub const kMipClockPerformAging0: i32 = 24;
pub const kMipClockSearch: i32 = 25;
pub const kMipClockProbingPresolve: i32 = 26;
pub const kMipClockEnumerationPresolve: i32 = 27;
pub const kMipClockPerformAging1: i32 = 28;
pub const kMipClockDive: i32 = 29;
pub const kMipClockOpenNodesToQueue0: i32 = 30;
pub const kMipClockDomainPropgate: i32 = 31;
pub const kMipClockPruneInfeasibleNodes: i32 = 32;
pub const kMipClockUpdateLocalDomain: i32 = 33;
pub const kMipClockNodeSearch: i32 = 34;
pub const kMipClockStartSymmetryDetection: i32 = 35;
pub const kMipClockStartAnalyticCentreComputation: i32 = 36;
pub const kMipClockEvaluateRootLp: i32 = 37;
pub const kMipClockSeparateLpCuts: i32 = 38;
pub const kMipClockRandomizedRounding: i32 = 39;
pub const kMipClockPerformRestart: i32 = 40;
pub const kMipClockRootSeparation: i32 = 41;
pub const kMipClockFinishAnalyticCentreComputation: i32 = 42;
pub const kMipClockRootCentralRounding: i32 = 43;
pub const kMipClockRootSeparationRound0: i32 = 44;
pub const kMipClockRootHeuristicsReducedCost: i32 = 45;
pub const kMipClockRootSeparationRound1: i32 = 46;
pub const kMipClockRootHeuristicsRens: i32 = 47;
pub const kMipClockRootSeparationRound2: i32 = 48;
pub const kMipClockRootFeasibilityPump: i32 = 49;
pub const kMipClockRootSeparationRound3: i32 = 50;
pub const kMipClockEvaluateRootNode0: i32 = 51;
pub const kMipClockEvaluateRootNode1: i32 = 52;
pub const kMipClockEvaluateRootNode2: i32 = 53;
pub const kMipClockDiveEvaluateNode: i32 = 54;
pub const kMipClockDivePrimalHeuristics: i32 = 55;
pub const kMipClockTheDive: i32 = 56;
pub const kMipClockBacktrackPlunge: i32 = 57;
pub const kMipClockPerformAging2: i32 = 58;
pub const kMipClockDiveRandomizedRounding: i32 = 59;
pub const kMipClockDiveRens: i32 = 60;
pub const kMipClockDiveRins: i32 = 61;
pub const kMipClockCurrentNodeToQueue: i32 = 62;
pub const kMipClockSearchBacktrack: i32 = 63;
pub const kMipClockNodePrunedLoop: i32 = 64;
pub const kMipClockOpenNodesToQueue1: i32 = 65;
pub const kMipClockEvaluateNode1: i32 = 66;
pub const kMipClockNodeSearchSeparation: i32 = 67;
pub const kMipClockStoreBasis: i32 = 68;
pub const kMipClockRootSeparationRound: i32 = 69;
pub const kMipClockRootSeparationFinishAnalyticCentreComputation: i32 = 70;
pub const kMipClockRootSeparationCentralRounding: i32 = 71;
pub const kMipClockRootSeparationEvaluateRootLp: i32 = 72;
pub const kMipClockImplboundSepa: i32 = 73;
pub const kMipClockCliqueSepa: i32 = 74;
pub const kMipClockTableauSepa: i32 = 75;
pub const kMipClockPathAggrSepa: i32 = 76;
pub const kMipClockModKSepa: i32 = 77;
pub const kMipClockDuSimplexBasisSolveLp: i32 = 78;
pub const kMipClockDuSimplexNoBasisSolveLp: i32 = 79;
pub const kMipClockHipoSolveAnalyticCentreLp: i32 = 80;
pub const kMipClockIpxSolveAnalyticCentreLp: i32 = 81;
pub const kMipClockHipoSolveLp: i32 = 82;
pub const kMipClockIpxSolveLp: i32 = 83;
pub const kMipClockPrSimplexBasisSolveLp: i32 = 84;
pub const kMipClockPrSimplexNoBasisSolveLp: i32 = 85;
pub const kMipClockSubMipSolve: i32 = 86;
pub const kMipClockProbingImplications: i32 = 87;
pub const kLastMipClock: i32 = 87;
pub const kToMipClock: i32 = 88;

/// HiGHS: MipTimer.h::initialiseMipProfilingNames (lines 123-211)
/// C++: `static void initialiseMipProfilingNames(std::vector<std::string> & name)`
pub fn initialise_mip_profiling_names(name: &mut Vec<String>) {
    // stub (timers / debug solution are not ported)
}

impl MipTimer {
    /// HiGHS: MipTimer.h::MipTimer::initialiseMipClocks (lines 215-388)
    /// C++: `void initialiseMipClocks(HighsTimerClock & mip_timer_clock, const HighsInt thread_mip_clock_offset)`
    pub fn initialise_mip_clocks(
        &mut self,
        mip_timer_clock: &mut HighsTimerClock,
        thread_mip_clock_offset: HighsInt,
    ) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipClockList (lines 390-410)
    /// C++: `bool reportMipClockList(const char * grepStamp, const std::vector<HighsInt> mip_clock_list, const HighsTimerClock & mip_timer_clock, const HighsInt kMipClockIdeal = kMipClockTotal, const double tolerance_percent_report_ = - 1)`
    /// C++ default arguments: `k_mip_clock_ideal = kMipClockTotal`, `tolerance_percent_report_ = - 1`
    pub fn report_mip_clock_list(
        &mut self,
        grep_stamp: &[u8],
        mip_clock_list: Vec<HighsInt>,
        mip_timer_clock: &HighsTimerClock,
        k_mip_clock_ideal: HighsInt,
        tolerance_percent_report_: f64,
    ) -> bool {
        // stub (timers / debug solution are not ported)
        Default::default()
    }

    /// HiGHS: MipTimer.h::MipTimer::csvMipClockList (lines 412-446)
    /// C++: `void csvMipClockList(const std::string grep_query, const std::string model_name, const std::vector<HighsInt> mip_clock_list, const HighsTimerClock & mip_timer_clock, const HighsInt kMipClockIdeal, const bool header, const bool end_line)`
    pub fn csv_mip_clock_list(
        &mut self,
        grep_query: &str,
        model_name: &str,
        mip_clock_list: Vec<HighsInt>,
        mip_timer_clock: &HighsTimerClock,
        k_mip_clock_ideal: HighsInt,
        header: bool,
        end_line: bool,
    ) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipCoreClock (lines 448-454)
    /// C++: `void reportMipCoreClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_core_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipLevel1Clock (lines 456-468)
    /// C++: `void reportMipLevel1Clock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_level1_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipSolveLpClock (lines 470-482)
    /// C++: `void reportMipSolveLpClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_solve_lp_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipSubMipSolveClock (lines 484-488)
    /// C++: `void reportMipSubMipSolveClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_sub_mip_solve_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipPresolveClock (lines 490-495)
    /// C++: `void reportMipPresolveClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_presolve_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportAltEvaluateRootNodeClock (lines 497-504)
    /// C++: `void reportAltEvaluateRootNodeClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_alt_evaluate_root_node_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipEvaluateRootNodeClock (lines 506-530)
    /// C++: `void reportMipEvaluateRootNodeClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_evaluate_root_node_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipRootSeparationClock (lines 532-540)
    /// C++: `void reportMipRootSeparationClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_root_separation_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipSearchClock (lines 542-552)
    /// C++: `void reportMipSearchClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_search_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipDiveClock (lines 554-560)
    /// C++: `void reportMipDiveClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_dive_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipDivePrimalHeuristicsClock (lines 562-569)
    /// C++: `void reportMipDivePrimalHeuristicsClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_dive_primal_heuristics_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipNodeSearchClock (lines 571-579)
    /// C++: `void reportMipNodeSearchClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_node_search_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportMipSeparationClock (lines 581-587)
    /// C++: `void reportMipSeparationClock(const HighsTimerClock & mip_timer_clock)`
    pub fn report_mip_separation_clock(&mut self, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::csvMipClock (lines 589-597)
    /// C++: `void csvMipClock(const std::string model_name, const HighsTimerClock & mip_timer_clock, const bool header, const bool end_line)`
    pub fn csv_mip_clock(
        &mut self,
        model_name: &str,
        mip_timer_clock: &HighsTimerClock,
        header: bool,
        end_line: bool,
    ) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::csvEvaluateRootNodeClock (lines 599-621)
    /// C++: `void csvEvaluateRootNodeClock(const std::string model_name, const HighsTimerClock & mip_timer_clock, const bool header, const bool end_line)`
    pub fn csv_evaluate_root_node_clock(
        &mut self,
        model_name: &str,
        mip_timer_clock: &HighsTimerClock,
        header: bool,
        end_line: bool,
    ) {
        // stub (timers / debug solution are not ported)
    }

    /// HiGHS: MipTimer.h::MipTimer::reportFjClock (lines 623-634)
    /// C++: `void reportFjClock(std::string & model, const HighsTimerClock & mip_timer_clock)`
    pub fn report_fj_clock(&mut self, model: &mut String, mip_timer_clock: &HighsTimerClock) {
        // stub (timers / debug solution are not ported)
    }
}

/// HiGHS: HStruct.h::HighsProfiling (lines 175-205), reduced to the three
/// presolve/solve/postsolve clocks that `solvingReport` prints (clocks past
/// `kToPresolveSolvePostsolve` are ignored exactly as C++ ignores clocks at or
/// past `num_profiling_clock_`). One MIP record: the port gives every
/// (sub-)MIP its own profiler, so the sub-MIP record stays empty and sub-MIP time is
/// inside the calling MIP's Solve time. Reporting only: nothing reads these clocks.
#[derive(Clone, Debug)]
pub struct HighsProfiling {
    /// C++ `timer->read()` origin: the clocks run on a monotonic instant.
    origin: std::time::Instant,
    num_call: [HighsInt; 3],
    run_time: [f64; 3],
    /// Start time of a running clock, `None` when stopped (C++ stores `-time_start`).
    start_time: [Option<f64>; 3],
}

impl Default for HighsProfiling {
    fn default() -> Self {
        Self {
            origin: std::time::Instant::now(),
            num_call: [0; 3],
            run_time: [0.0; 3],
            start_time: [None; 3],
        }
    }
}

impl HighsProfiling {
    fn now(&self) -> f64 {
        self.origin.elapsed().as_secs_f64()
    }
    /// HiGHS: HStruct.cpp::HighsProfiling::start
    pub fn start(&mut self, profiling_clock: HighsInt, restart: bool) {
        let ci = profiling_clock as usize;
        if ci >= 3 {
            return;
        }
        self.start_time[ci] = Some(self.now());
        if restart {
            self.num_call[ci] -= 1;
        }
    }
    /// HiGHS: HStruct.cpp::HighsProfiling::stop
    pub fn stop(&mut self, profiling_clock: HighsInt) {
        let ci = profiling_clock as usize;
        if ci >= 3 {
            return;
        }
        if let Some(t0) = self.start_time[ci].take() {
            self.num_call[ci] += 1;
            self.run_time[ci] += self.now() - t0;
        }
    }
    /// HiGHS: HStruct.cpp::HighsProfiling::read (the MIP record)
    pub fn read(&self, profiling_clock: HighsInt) -> f64 {
        let ci = profiling_clock as usize;
        if ci >= 3 {
            return -kHighsInf;
        }
        let running = self.start_time[ci].map_or(0.0, |t0| self.now() - t0);
        self.run_time[ci] + running
    }
    /// HiGHS: HStruct.cpp::HighsProfiling::numCall (the MIP record)
    pub fn num_call(&self, profiling_clock: HighsInt) -> HighsInt {
        let ci = profiling_clock as usize;
        if ci >= 3 {
            return 0;
        }
        self.num_call[ci] + HighsInt::from(self.start_time[ci].is_some())
    }
    /// HiGHS: HStruct.cpp::HighsProfiling::setSubMip (stub)
    pub fn set_sub_mip(&mut self, submip: bool) {}
    /// HiGHS: HStruct.cpp::HighsProfiling::isSubMip (stub)
    pub fn is_sub_mip(&self) -> bool {
        false
    }
    /// HiGHS: HStruct.h::HighsProfiling::mip_ (line 180; a `bool` member, `false` unless
    /// MIP profiling is on, which this port never turns on). A method because the stub
    /// has no fields (HighsLpRelaxation::run reads it).
    pub fn mip_(&self) -> bool {
        false
    }
    /// HiGHS: HStruct.cpp::HighsProfiling::myThread (stub)
    pub fn my_thread(&self) -> HighsInt {
        0
    }
    /// HiGHS: HStruct.cpp::HighsProfiling::running (the MIP record)
    pub fn running(&self, profiling_clock: HighsInt) -> bool {
        (profiling_clock as usize) < 3 && self.start_time[profiling_clock as usize].is_some()
    }
    /// HiGHS: HStruct.cpp::HighsProfiling::solveCall (stub)
    pub fn solve_call(&mut self, model: &str, submip: bool) {}
}
