//! HiGHS source: io/HighsIO.h, io/HighsIO.cpp (minimal subset)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Logging is reduced to what the simplex path calls: `highsLogUser`,
//! `highsLogDev`, `highsFormatToString` (-> `format!`), `highsBoolToString`.
//! Translate a C++ call
//! `highsLogDev(options_->log_options, HighsLogType::kDetailed, "x = %g\n", x)`
//! as `highs_log_dev!(&self.options_.log_options, HighsLogType::kDetailed, "x = {}\n", x)`.
//! printf conversions: `%d`/`%" HIGHSINT_FORMAT "` -> `{}`, `%s` -> `{}`, and the
//! floating-point ones through `crate::util::c_printf` (`%g` -> `c_fmt_g(x, 6)`,
//! `%11.4g` -> `c_fmt_g_width(x, 4, 11)`, ...): the MIP oracle trace compares log text,
//! and Rust's `{}` float form is not printf's. (Older LP-path call sites still use
//! `{}`; their text is not compared.)
//! `FILE*` streams are only ever stdout. The `kCallbackLogging` route
//! (`user_callback` / `user_callback_active`) is ported; the C-style
//! `user_log_callback` is not (no caller in the MIP or LP path).

use crate::HighsInt;
use crate::lp_data::hconst::{HighsLogType, kCallbackLogging};
use crate::lp_data::highs_callback::{
    HighsCallbackFunctionType, HighsCallbackInput, HighsCallbackOutput,
};
use crate::util::c_printf::c_fmt_g;

/// HiGHS: HighsIO.h::kIoBufferSize
pub const kIoBufferSize: HighsInt = 1024;

/// HiGHS: HighsIO.h LogDevLevel
pub const kHighsLogDevLevelMin: HighsInt = 0;
pub const kHighsLogDevLevelNone: HighsInt = 0;
pub const kHighsLogDevLevelInfo: HighsInt = 1;
pub const kHighsLogDevLevelDetailed: HighsInt = 2;
pub const kHighsLogDevLevelVerbose: HighsInt = 3;
pub const kHighsLogDevLevelMax: HighsInt = 3;

/// HiGHS: HighsIO.h HighsLogTypeTag
pub const HighsLogTypeTag: [&str; 6] = ["", "", "", "", "WARNING: ", "ERROR:   "];

/// Replacement for C++ `FILE*`: the only stream the simplex path writes to is
/// stdout. `None` stands for a null `FILE*`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct HighsLogStream;

/// HiGHS: HighsIO.h::HighsLogOptions
///
/// In C++ the flags are pointers into `HighsOptions`; here they are owned
/// copies, kept in sync by `HighsOptions` (which owns its `log_options`).
///
/// `user_callback` is the `HighsCallbackFunctionType` that `Highs::setCallback`
/// copies in (C++ `user_callback_data` is captured by the closure), and
/// `user_callback_active` is set by `startCallback(kCallbackLogging)`.
/// `user_log_callback` / `user_log_callback_data` are not ported (always null).
#[derive(Clone)]
pub struct HighsLogOptions {
    pub log_stream: Option<HighsLogStream>,
    pub output_flag: bool,
    pub log_to_console: bool,
    pub log_dev_level: HighsInt,
    pub user_callback: Option<HighsCallbackFunctionType>,
    pub user_callback_active: bool,
}

impl std::fmt::Debug for HighsLogOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HighsLogOptions")
            .field("log_stream", &self.log_stream)
            .field("output_flag", &self.output_flag)
            .field("log_to_console", &self.log_to_console)
            .field("log_dev_level", &self.log_dev_level)
            .field("user_callback", &self.user_callback.is_some())
            .field("user_callback_active", &self.user_callback_active)
            .finish()
    }
}

impl Default for HighsLogOptions {
    /// Values as seen through the pointers of a default `HighsOptions`.
    fn default() -> Self {
        Self {
            log_stream: None,
            output_flag: true,
            log_to_console: true,
            log_dev_level: 0,
            user_callback: None,
            user_callback_active: false,
        }
    }
}

impl HighsLogOptions {
    /// Quiet options: `output_flag = false` (what tests normally want).
    pub fn quiet() -> Self {
        Self {
            output_flag: false,
            ..Self::default()
        }
    }
}

/// C++ `log_to_callback` (HighsIO.cpp highsLogUser / highsLogDev):
/// `user_log_callback || (user_callback && user_callback_active)`.
fn log_to_callback(log_options: &HighsLogOptions) -> bool {
    log_options.user_callback.is_some() && log_options.user_callback_active
}

/// C++ `std::array<char, kIoBufferSize> msgbuffer` filled by `snprintf`: at most
/// `kIoBufferSize - 1` bytes survive (the terminating NUL takes the last one).
fn truncate_to_buffer(msg: &mut String) {
    let cap = (kIoBufferSize - 1) as usize;
    if msg.len() > cap {
        let mut end = cap;
        while !msg.is_char_boundary(end) {
            end -= 1;
        }
        msg.truncate(end);
    }
}

/// The `kCallbackLogging` call of highsLogUser / highsLogDev: `data_out.log_type = type`,
/// `nullptr` callback input (a default one here: the closure signature takes `&mut`).
fn fire_log_callback(log_options: &HighsLogOptions, log_type: HighsLogType, msg: &str) {
    if !log_options.user_callback_active {
        return;
    }
    let cb = log_options
        .user_callback
        .as_ref()
        .expect("user_callback_active without user_callback");
    let data_out = HighsCallbackOutput {
        log_type: log_type as HighsInt,
        ..Default::default()
    };
    let mut data_in = HighsCallbackInput::default();
    (cb.borrow_mut())(kCallbackLogging, msg, &data_out, &mut data_in);
}

/// HiGHS: HighsIO.cpp::highsLogUser (format string already expanded).
pub fn highs_log_user(log_options: &HighsLogOptions, log_type: HighsLogType, message: &str) {
    if !log_options.output_flag {
        return;
    }
    // Log to a file if log_stream is not NULL or stdout (never, here), and to
    // console if log_to_console is true
    let log_to_console = log_options.log_to_console;
    let log_to_callback = log_to_callback(log_options);
    if !log_to_console && !log_to_callback {
        return;
    }
    debug_assert!(log_type != HighsLogType::kDetailed);
    debug_assert!(log_type != HighsLogType::kVerbose);
    let prefix = log_type == HighsLogType::kWarning || log_type == HighsLogType::kError;
    // Possibly write to console
    if log_to_console {
        if prefix {
            print!("{:<9}", HighsLogTypeTag[log_type as usize]);
        }
        print!("{message}");
    }
    // Possibly write to callback
    if log_to_callback {
        let mut msgbuffer = String::new();
        if prefix {
            msgbuffer.push_str(&format!("{:<9}", HighsLogTypeTag[log_type as usize]));
        }
        msgbuffer.push_str(message);
        truncate_to_buffer(&mut msgbuffer);
        fire_log_callback(log_options, log_type, &msgbuffer);
    }
}

/// HiGHS: HighsIO.cpp::highsLogDev (format string already expanded).
pub fn highs_log_dev(log_options: &HighsLogOptions, log_type: HighsLogType, message: &str) {
    if !log_options.output_flag || log_options.log_dev_level == 0 {
        return;
    }
    let log_to_console = log_options.log_to_console;
    let log_to_callback = log_to_callback(log_options);
    if !log_to_console && !log_to_callback {
        return;
    }
    if log_type == HighsLogType::kDetailed && log_options.log_dev_level < kHighsLogDevLevelDetailed
    {
        return;
    }
    if log_type == HighsLogType::kVerbose && log_options.log_dev_level < kHighsLogDevLevelVerbose {
        return;
    }
    // Possibly write to console
    if log_to_console {
        print!("{message}");
    }
    // Possibly write to callback
    if log_to_callback {
        let mut msgbuffer = message.to_string();
        truncate_to_buffer(&mut msgbuffer);
        fire_log_callback(log_options, log_type, &msgbuffer);
    }
}

/// Should a `highsLogDev` call with this type print? Used by the macro to
/// skip formatting entirely when it would not.
pub fn highs_log_dev_enabled(log_options: &HighsLogOptions, log_type: HighsLogType) -> bool {
    if !log_options.output_flag
        || log_options.log_dev_level == 0
        || (!log_options.log_to_console && !log_to_callback(log_options))
    {
        return false;
    }
    !(log_type == HighsLogType::kDetailed && log_options.log_dev_level < kHighsLogDevLevelDetailed
        || log_type == HighsLogType::kVerbose
            && log_options.log_dev_level < kHighsLogDevLevelVerbose)
}

/// Should a `highsLogUser` call print?
pub fn highs_log_user_enabled(log_options: &HighsLogOptions) -> bool {
    log_options.output_flag && (log_options.log_to_console || log_to_callback(log_options))
}

/// `highsLogUser(log_options, type, format, ...)` -> `highs_log_user!(log_options, type, "fmt", args...)`.
/// Formats lazily (only when the message would be printed).
#[macro_export]
macro_rules! highs_log_user {
    ($opts:expr, $ty:expr, $($arg:tt)*) => {{
        let __o: &$crate::io::highs_io::HighsLogOptions = $opts;
        if $crate::io::highs_io::highs_log_user_enabled(__o) {
            $crate::io::highs_io::highs_log_user(__o, $ty, &format!($($arg)*));
        }
    }};
}

/// `highsLogDev(log_options, type, format, ...)` -> `highs_log_dev!(log_options, type, "fmt", args...)`.
/// Formats lazily (only when the message would be printed).
#[macro_export]
macro_rules! highs_log_dev {
    ($opts:expr, $ty:expr, $($arg:tt)*) => {{
        let __o: &$crate::io::highs_io::HighsLogOptions = $opts;
        if $crate::io::highs_io::highs_log_dev_enabled(__o, $ty) {
            $crate::io::highs_io::highs_log_dev(__o, $ty, &format!($($arg)*));
        }
    }};
}

/// HiGHS: HighsIO.cpp::highsBoolToString
/// C++ default arguments: field_width = 2
pub fn highs_bool_to_string(b: bool, field_width: HighsInt) -> String {
    let abs_field_width = field_width.abs();
    if abs_field_width <= 1 {
        return if b { "T" } else { "F" }.to_string();
    }
    if abs_field_width <= 2 {
        return if b { "true" } else { "false" }.to_string();
    }
    if field_width < 0 {
        return if b { "true " } else { "false" }.to_string();
    }
    if b { " true" } else { "false" }.to_string()
}

/// HiGHS: HighsIO.cpp::highsDoubleToString (lines 39-100)
/// C++: `std::array<char, 32> highsDoubleToString(const double val, const double tolerance)`
pub fn highs_double_to_string(val: f64, tolerance: f64) -> String {
    let l: f64 = if val.abs() == crate::lp_data::hconst::kHighsInf {
        1.0
    } else {
        1.0 - tolerance
            + (crate::util::highs_utils::cpp_max(tolerance, val.abs()) / tolerance).log10()
    };
    // C++ `switch (int(l))`: truncation toward zero; `%.{n}g` for n in 1..=16.
    match l as i32 {
        0 => "0".to_string(),
        n @ 1..=15 => c_fmt_g(val, n as usize),
        _ => c_fmt_g(val, 16),
    }
}
