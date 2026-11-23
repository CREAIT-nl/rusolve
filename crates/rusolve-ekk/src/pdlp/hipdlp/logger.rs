//! HiGHS source: pdlp/hipdlp/logger.hpp, pdlp/hipdlp/logger.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiPDLP logger. Every method other than `initialise` and
//! `getLogLevel` only prints (through `highsLogUser`), and nothing it
//! changes (`iteration_stats_count_`, `iteration_stats_time_`) is read
//! outside printing, so those methods are complete no-ops here.
//! `highs_timer_p_` is read only by
//! `printIterationStats` and is not stored. `format` / `integer` / `sci` /
//! `fix` are string builders for printing and are not ported.

use super::defs::PrimalDualParams;
use super::solver_results::SolverResults;
use crate::prelude::*;

/// HiGHS: logger.hpp `enum class LogLevel`
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub enum LogLevel {
    /// No output
    #[default]
    kNone,
    /// Standard output: summary, termination, major events
    kInfo,
    /// Developer output
    kDetailed,
    /// Detailed output: iteration-level info
    kVerbose,
    /// Verbose + debug info for developers
    kDebug,
}

/// HiGHS: logger.hpp `class Logger`
#[derive(Clone, Debug, Default)]
pub struct Logger {
    pub(crate) console_level_: LogLevel,
    pub(crate) log_options_: HighsLogOptions,
    pub(crate) iteration_stats_count_: HighsInt,
    pub(crate) iteration_stats_time_: f64,
}

impl Logger {
    /// HiGHS: logger.hpp `kHipdlpIterationStatsHeaderFrequency`
    pub const kHipdlpIterationStatsHeaderFrequency: HighsInt = 50;
    /// HiGHS: logger.hpp `kHipdlpIterationStatsFrequency`
    pub const kHipdlpIterationStatsFrequency: f64 = 5.0;

    /// HiGHS: logger.cc::Logger::initialise (lines 25-41)
    ///
    /// The `HighsTimer*` parameter is dropped (print-only use).
    pub fn initialise(&mut self, log_dev_level: HighsInt, log_options: &HighsLogOptions) {
        if log_dev_level == kHighsLogDevLevelVerbose {
            self.console_level_ = LogLevel::kDebug;
        } else if log_dev_level == kHighsLogDevLevelDetailed {
            self.console_level_ = LogLevel::kVerbose;
        } else if log_dev_level == kHighsLogDevLevelInfo {
            self.console_level_ = LogLevel::kDetailed;
        } else {
            self.console_level_ = LogLevel::kInfo;
        }
        self.log_options_ = log_options.clone();
        // highs_timer_p_ = highs_timer_p; dropped (print-only use).
        self.iteration_stats_count_ = kHighsIInf;
        self.iteration_stats_time_ = -kHighsInf;
    }

    /// HiGHS: logger.hpp::Logger::getLogLevel
    pub fn getLogLevel(&self) -> LogLevel {
        self.console_level_
    }

    /// HiGHS: logger.hpp::Logger::setLogOptions
    pub fn setLogOptions(&mut self, log_options: &HighsLogOptions) {
        self.log_options_ = log_options.clone();
    }

    /// HiGHS: logger.cc::Logger::log (lines 43-47). Print only.
    fn log(&self, level: LogLevel, message: &str) {
        let _ = (level, message);
    }

    /// HiGHS: logger.cc::Logger::info (lines 49-51). Print only.
    pub fn info(&self, message: &str) {
        self.log(LogLevel::kInfo, message);
    }

    /// HiGHS: logger.cc::Logger::detailed (lines 52-54). Print only.
    pub fn detailed(&self, message: &str) {
        self.log(LogLevel::kDetailed, message);
    }

    /// HiGHS: logger.cc::Logger::verbose (lines 55-57). Print only.
    pub fn verbose(&self, message: &str) {
        self.log(LogLevel::kVerbose, message);
    }

    /// HiGHS: logger.cc::Logger::debug (lines 58-60). Print only.
    pub fn debug(&self, message: &str) {
        self.log(LogLevel::kDebug, message);
    }

    /// HiGHS: logger.cc::Logger::printHeader (lines 62-64). Print only.
    pub fn printHeader(&self) {}

    /// HiGHS: logger.cc::Logger::printParams (lines 66-129). Print only.
    pub fn printParams(&self, params: &PrimalDualParams) {
        let _ = params;
    }

    /// HiGHS: logger.cc::Logger::printIterationHeader (lines 131-136). Print only.
    pub fn printIterationHeader(&mut self) {}

    /// HiGHS: logger.cc::Logger::printIterationStats (lines 138-166). Print only.
    ///
    /// C++ default arguments: `forced = false`.
    pub fn printIterationStats(
        &mut self,
        iter: HighsInt,
        current_results: &SolverResults,
        current_eta: f64,
        forced: bool,
    ) {
        let _ = (iter, current_results, current_eta, forced);
    }

    /// HiGHS: logger.cc::Logger::printSummary (lines 168-208). Print only.
    pub fn printSummary(&self, results: &SolverResults, total_iter: HighsInt, total_time: f64) {
        let _ = (results, total_iter, total_time);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn dev_level_maps_to_console_level() {
        let mut l = Logger::default();
        l.initialise(kHighsLogDevLevelNone, &HighsLogOptions::quiet());
        assert_eq!(l.getLogLevel(), LogLevel::kInfo);
        assert_eq!(l.iteration_stats_count_, kHighsIInf);
        assert_eq!(l.iteration_stats_time_, -kHighsInf);
        l.initialise(kHighsLogDevLevelVerbose, &HighsLogOptions::quiet());
        assert_eq!(l.getLogLevel(), LogLevel::kDebug);
    }
}
