//! HiGHS source: ipm/hipo/auxiliary/Logger.h, ipm/hipo/auxiliary/Logger.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiPO logging. The C++ printf-style templates (`print`, `printw`, ...) take
//! `(const char* format, Args...)` and forward to `highsLogUser` /
//! `highsLogDev`; the Rust takes the message already formatted (`&format!(..)`
//! with the C format translated to Rust format syntax), exactly as the
//! rest of the crate calls `highs_log_user`.

use super::super::Int;
use crate::io::highs_io::{
    HighsLogOptions, highs_log_dev, highs_log_user, kHighsLogDevLevelDetailed,
    kHighsLogDevLevelInfo, kHighsLogDevLevelVerbose,
};
use crate::lp_data::hconst::HighsLogType;
use crate::util::c_printf::{c_fmt_e, c_fmt_f, c_pad_right};

/// HiGHS: Logger.h `class Logger`. `log_options_` is a clone of the
/// `HighsLogOptions` the C++ points at (`None` = nullptr).
#[derive(Clone, Default)]
pub struct Logger {
    pub log_options_: Option<HighsLogOptions>,
}

impl Logger {
    /// HiGHS: Logger.cpp::Logger::setOptions (lines 5-7)
    pub fn setOptions(&mut self, log_options: Option<&HighsLogOptions>) {
        self.log_options_ = log_options.cloned();
    }

    /// HiGHS: Logger.cpp::Logger::debug (lines 9-23)
    pub fn debug(&self, level: Int) -> bool {
        // Return true if level agrees with log_dev_level
        // PORT NOTE: C++ `log_options_->log_dev_level` is a `HighsInt*` tested for
        // non-null; the Rust field is the value itself, so only `log_options_` is tested.
        if let Some(o) = &self.log_options_ {
            if o.log_dev_level == kHighsLogDevLevelInfo {
                return level == 1;
            }

            if o.log_dev_level == kHighsLogDevLevelDetailed {
                return level == 1 || level == 2;
            }

            if o.log_dev_level == kHighsLogDevLevelVerbose {
                return level == 1 || level == 2 || level == 3;
            }
        }
        false
    }

    /// HiGHS: Logger.h::Logger::print (highsLogUser kInfo)
    pub fn print(&self, msg: &str) {
        if let Some(o) = &self.log_options_ {
            highs_log_user(o, HighsLogType::kInfo, msg);
        }
    }

    /// HiGHS: Logger.h::Logger::printw (highsLogUser kWarning)
    pub fn printw(&self, msg: &str) {
        if let Some(o) = &self.log_options_ {
            highs_log_user(o, HighsLogType::kWarning, msg);
        }
    }

    /// HiGHS: Logger.h::Logger::printe (highsLogUser kError)
    pub fn printe(&self, msg: &str) {
        if let Some(o) = &self.log_options_ {
            highs_log_user(o, HighsLogType::kError, msg);
        }
    }

    /// HiGHS: Logger.h::Logger::printInfo (highsLogDev kInfo)
    pub fn printInfo(&self, msg: &str) {
        if let Some(o) = &self.log_options_ {
            highs_log_dev(o, HighsLogType::kInfo, msg);
        }
    }

    /// HiGHS: Logger.h::Logger::printDetailed (highsLogDev kDetailed)
    pub fn printDetailed(&self, msg: &str) {
        if let Some(o) = &self.log_options_ {
            highs_log_dev(o, HighsLogType::kDetailed, msg);
        }
    }

    /// HiGHS: Logger.h::Logger::printVerbose (highsLogDev kVerbose)
    pub fn printVerbose(&self, msg: &str) {
        if let Some(o) = &self.log_options_ {
            highs_log_dev(o, HighsLogType::kVerbose, msg);
        }
    }
}

/// `std::ios_base::fmtflags` values `format` distinguishes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FloatField {
    scientific,
    fixed,
}

/// HiGHS: Logger.cpp::format (lines 25-33)
/// `std::ostringstream` with `precision(prec)`, `width(width)`, `setf(floatfield)`:
/// the libc `%.{prec}e` / `%.{prec}f` conversion, right-aligned in `width` (as IPX's
/// `format_double` in `ipx/control.rs`, which this function was taken from).
pub fn format(d: f64, width: Int, prec: Int, floatfield: FloatField) -> String {
    let s = match floatfield {
        FloatField::scientific => c_fmt_e(d, prec as usize),
        FloatField::fixed => c_fmt_f(d, prec as usize),
    };
    c_pad_right(&s, width.max(0) as usize)
}

/// HiGHS: Logger.cpp::integer (lines 35-40)
/// C++ default arguments: `width = 0`.
pub fn integer(i: Int, width: Int) -> String {
    c_pad_right(&i.to_string(), width.max(0) as usize)
}

/// HiGHS: Logger.h::sci (inline)
pub fn sci(d: f64, width: Int, prec: Int) -> String {
    format(d, width, prec, FloatField::scientific)
}

/// HiGHS: Logger.h::fix (inline)
pub fn fix(d: f64, width: Int, prec: Int) -> String {
    format(d, width, prec, FloatField::fixed)
}

/// HiGHS: Logger.h::textline (template): left-aligned, width 32.
pub fn textline(text: &str) -> String {
    format!("{text:<32}")
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn format_and_integer() {
        assert_eq!(sci(1e-8, 9, 1), "  1.0e-08");
        assert_eq!(fix(1.23456, 6, 2), "  1.23");
        assert_eq!(integer(42, 5), "   42");
        assert_eq!(sci(f64::NAN, 5, 1), "  nan");
        assert_eq!(textline("abc").len(), 32);
    }
}
