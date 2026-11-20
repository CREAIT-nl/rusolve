//! HiGHS source: ipm/ipx/control.h, control.cc, multistream.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `Control` is shared as `Rc<Control>` (see "Ownership" in `ipx/mod.rs`). Its C++
//! `mutable` members and the members changed through `LpSolver` while
//! clones exist (`timer_` via `setTimerOffset`/`ResetTimer`) use `Cell`/
//! `RefCell`. `parameters_` is a `RefCell` too, because C++
//! `LpSolver::SetParameters` may run after `Solve`, while `basis_` (and the
//! `Basis`'s `const Control&`) still exist; the C++ objects then see the new
//! parameters, and so do the Rust clones.
//!
//! Streams: `output_` (a `Multistream` forwarding to `std::cout` and/or the
//! logfile) is the pair of flags in [`Multistream`]; writing goes through
//! [`Control::write_output`]. `dummy_` (failbit stream) is "write nothing".
//! `Debug(level) << a << b` translates to
//! `control.debug(level, &format!(..))`: C++ evaluates the operands even
//! when the level is off, so the Rust call formats eagerly too; a caller may
//! guard an expensive operand with `if control.debug_enabled(level)` only
//! when computing it has no side effects.
//!
//! Number formatting (`Format` -> `format_int`/`format_str`/`format_double`, `Scientific` -> `scientific`, `sci2`, `fix8`, `Textline` -> `textline`, ...) is
//! implemented here on top of `util::c_printf` (libc++ ostream `scientific`/
//! `fixed` = C `%.*e`/`%.*f`; default floatfield = `%.*g` with precision 6).

use super::timer::Timer;
use super::*;
use crate::lp_data::hconst::{HighsLogType, kCallbackIpmInterrupt};
use crate::lp_data::highs_callback::HighsCallback;
use crate::util::c_printf::*;
use std::cell::{Cell, RefCell};
use std::fs::File;
use std::io::Write;

/// HiGHS: ipx/multistream.h::Multistream (as the set of sinks it forwards to)
#[derive(Debug, Clone, Copy, Default)]
pub struct Multistream {
    /// `output_.add(std::cout)`
    pub to_console: bool,
    /// `output_.add(logfile_)`
    pub to_logfile: bool,
}

/// HiGHS: ipx/control.h::Control
pub struct Control {
    parameters_: RefCell<Parameters>,
    /// C++ `HighsCallback* callback_` (pointer to the caller's callback).
    /// PORT NOTE: the callback is moved in by `LpSolver::SetCallback` and
    /// moved back out by [`Control::take_callback`] after the solve.
    callback_: RefCell<Option<HighsCallback>>,
    logfile_: RefCell<Option<File>>,
    timer_: Cell<Timer>,
    interval_: Cell<Timer>,
    output_: Cell<Multistream>,
}

impl Control {
    /// HiGHS: ipx/control.cc::Control::Control (lines 8-11)
    pub fn new() -> Self {
        Control {
            parameters_: RefCell::new(Parameters::default()),
            callback_: RefCell::new(None),
            logfile_: RefCell::new(None),
            timer_: Cell::new(Timer::default()),
            interval_: Cell::new(Timer::default()),
            output_: Cell::new(Multistream::default()),
        }
    }

    /// HiGHS: ipx/control.cc::Control::InterruptCheck (lines 13-31)
    ///
    /// C++ default arguments: `ipm_iteration_count = -1`.
    /// PORT NOTE: `HighsTaskExecutor::getThisWorkerDeque()->checkInterrupt()`
    /// is a no-op (no task executor in this port).
    pub fn interrupt_check(&self, ipm_iteration_count: Int) -> Int {
        if self.parameters_.borrow().time_limit >= 0.0
            && self.parameters_.borrow().time_limit < self.elapsed()
        {
            return IPX_ERROR_time_interrupt;
        }
        // The pointer callback_ should not be null, since that indicates
        // that it's not been set
        let mut callback = self.callback_.borrow_mut();
        debug_assert!(callback.is_some());
        if let Some(cb) = callback.as_mut() {
            if cb.user_callback.is_some() && cb.active[kCallbackIpmInterrupt as usize] {
                cb.clear_highs_callback_output();
                cb.data_out.ipm_iteration_count = ipm_iteration_count;
                if cb.callback_action(kCallbackIpmInterrupt, "IPM interrupt") {
                    return IPX_ERROR_user_interrupt;
                }
            }
        }
        0
    }

    /// HiGHS: ipx/control.cc::Control::hLog(std::stringstream&) (lines 45-56)
    ///
    /// Clears `logging` (C++ `logging.str(std::string())`).
    pub fn h_log(&self, logging: &mut String) {
        if self.parameters_.borrow().highs_logging {
            self.highs_log(logging);
        } else {
            self.write_output(logging);
        }
        logging.clear();
        // Reset interval-based logging since something has been logged
        self.reset_print_interval();
    }

    /// HiGHS: ipx/control.cc::Control::hLog(std::string) (lines 33-43)
    pub fn h_log_str(&self, str: &str) {
        if self.parameters_.borrow().highs_logging {
            self.highs_log(str);
        } else {
            self.write_output(str);
        }
        // Reset interval-based logging since something has been logged
        self.reset_print_interval();
    }

    /// HiGHS: ipx/control.cc::Control::hIntervalLog (lines 58-72)
    pub fn h_interval_log(&self, logging: &mut String) {
        if self.parameters_.borrow().print_interval >= 0.0
            && self.interval_.get().elapsed() >= self.parameters_.borrow().print_interval
        {
            self.reset_print_interval();
            if self.parameters_.borrow().highs_logging {
                self.highs_log(logging);
            } else {
                self.write_output(logging);
            }
        }
        logging.clear();
    }

    /// HiGHS: ipx/control.cc::Control::Debug (lines 74-79)
    ///
    /// `Debug(level) << text`: writes `text` to `output_` iff
    /// `parameters_.debug >= level`. C++ default arguments: `level = 1`.
    pub fn debug(&self, level: Int, text: &str) {
        if self.debug_enabled(level) {
            self.write_output(text);
        }
    }

    /// `parameters_.debug >= level` (the branch inside `Debug`).
    pub fn debug_enabled(&self, level: Int) -> bool {
        self.parameters_.borrow().debug >= level
    }

    /// `output_ << text` (Multistream: console then logfile).
    pub fn write_output(&self, text: &str) {
        let ms = self.output_.get();
        if ms.to_console {
            let _ = std::io::stdout().write_all(text.as_bytes());
        }
        if ms.to_logfile {
            if let Some(file) = self.logfile_.borrow_mut().as_mut() {
                let _ = file.write_all(text.as_bytes());
            }
        }
    }

    /// HiGHS: ipx/control.cc::Control::ResetPrintInterval (line 81)
    pub fn reset_print_interval(&self) {
        let mut t = self.interval_.get();
        t.reset(false);
        self.interval_.set(t);
    }

    /// HiGHS: ipx/control.cc::Control::Elapsed (line 85)
    pub fn elapsed(&self) -> f64 {
        let t = self.timer_.get();
        t.offset_ + t.elapsed()
    }

    pub fn dualize(&self) -> Int {
        self.parameters_.borrow().dualize
    }
    pub fn scale(&self) -> Int {
        self.parameters_.borrow().scale
    }
    pub fn ipm_maxiter(&self) -> Int {
        self.parameters_.borrow().ipm_maxiter
    }
    pub fn ipm_feasibility_tol(&self) -> f64 {
        self.parameters_.borrow().ipm_feasibility_tol
    }
    pub fn ipm_optimality_tol(&self) -> f64 {
        self.parameters_.borrow().ipm_optimality_tol
    }
    pub fn ipm_drop_primal(&self) -> f64 {
        self.parameters_.borrow().ipm_drop_primal
    }
    pub fn ipm_drop_dual_soft(&self) -> bool {
        self.parameters_.borrow().ipm_drop_dual_soft
    }
    pub fn ipm_drop_dual(&self) -> f64 {
        self.parameters_.borrow().ipm_drop_dual
    }
    pub fn kkt_tol(&self) -> f64 {
        self.parameters_.borrow().kkt_tol
    }
    /// rusolve BASICLU time limit (not in the HiGHS reference): ipx/control.h::Control::timeLimit
    ///
    /// PORT NOTE: `Parameters::time_limit` defaults to -1 (no limit) here, as
    /// in `interrupt_check`; HiGHS always sets it before calling IPX, where
    /// the C++ asserts `>= 0`. A negative value is returned as infinity so
    /// that unset parameters never make the LU time out.
    pub fn time_limit(&self) -> f64 {
        let t = self.parameters_.borrow().time_limit;
        if t < 0.0 { f64::INFINITY } else { t }
    }
    pub fn crash_basis(&self) -> Int {
        self.parameters_.borrow().crash_basis
    }
    pub fn dependency_tol(&self) -> f64 {
        self.parameters_.borrow().dependency_tol
    }
    pub fn volume_tol(&self) -> f64 {
        self.parameters_.borrow().volume_tol
    }
    pub fn rows_per_slice(&self) -> Int {
        self.parameters_.borrow().rows_per_slice
    }
    pub fn maxskip_updates(&self) -> Int {
        self.parameters_.borrow().maxskip_updates
    }
    pub fn lu_kernel(&self) -> Int {
        self.parameters_.borrow().lu_kernel
    }
    pub fn lu_pivottol(&self) -> f64 {
        self.parameters_.borrow().lu_pivottol
    }
    pub fn run_crossover(&self) -> Int {
        self.parameters_.borrow().run_crossover
    }
    pub fn start_crossover_tol(&self) -> f64 {
        self.parameters_.borrow().start_crossover_tol
    }
    pub fn pfeasibility_tol(&self) -> f64 {
        self.parameters_.borrow().pfeasibility_tol
    }
    pub fn dfeasibility_tol(&self) -> f64 {
        self.parameters_.borrow().dfeasibility_tol
    }
    pub fn switchiter(&self) -> Int {
        self.parameters_.borrow().switchiter
    }
    pub fn stop_at_switch(&self) -> Int {
        self.parameters_.borrow().stop_at_switch
    }
    pub fn update_heuristic(&self) -> Int {
        self.parameters_.borrow().update_heuristic
    }
    pub fn maxpasses(&self) -> Int {
        self.parameters_.borrow().maxpasses
    }
    /// HiGHS: ipx/control.h::Control::reportBasisData
    pub fn report_basis_data(&self) -> bool {
        self.parameters_.borrow().analyse_basis_data
    }
    /// HiGHS: ipx/control.h::Control::runCentring
    pub fn run_centring(&self) -> Int {
        self.parameters_.borrow().run_centring
    }
    /// HiGHS: ipx/control.h::Control::maxCentringSteps
    pub fn max_centring_steps(&self) -> Int {
        self.parameters_.borrow().max_centring_steps
    }
    /// HiGHS: ipx/control.h::Control::centringRatioTolerance
    pub fn centring_ratio_tolerance(&self) -> f64 {
        self.parameters_.borrow().centring_ratio_tolerance
    }
    /// HiGHS: ipx/control.h::Control::centringRatioReduction
    pub fn centring_ratio_reduction(&self) -> f64 {
        self.parameters_.borrow().centring_ratio_reduction
    }
    /// HiGHS: ipx/control.h::Control::centringAlphaScaling
    pub fn centring_alpha_scaling(&self) -> f64 {
        self.parameters_.borrow().centring_alpha_scaling
    }
    /// HiGHS: ipx/control.h::Control::badProductsTolerance
    pub fn bad_products_tolerance(&self) -> Int {
        self.parameters_.borrow().bad_products_tolerance
    }
    /// HiGHS: ipx/control.h::Control::timelessLog
    pub fn timeless_log(&self) -> bool {
        self.parameters_.borrow().timeless_log
    }
    /// rusolve (HiGHS 1.16 backport, default off): Parameters::h16_ipx_stall
    pub fn h16_ipx_stall(&self) -> bool {
        self.parameters_.borrow().h16_ipx_stall
    }
    /// rusolve `lp_crossover_batch_push`: Parameters::crossover_batch_push
    pub fn crossover_batch_push(&self) -> bool {
        self.parameters_.borrow().crossover_batch_push
    }

    /// rusolve `lp_crossover_ipm_crash`: Parameters::crossover_ipm_crash
    pub fn crossover_ipm_crash(&self) -> bool {
        self.parameters_.borrow().crossover_ipm_crash
    }
    /// rusolve: Parameters::stall_to_crossover
    pub fn stall_to_crossover(&self) -> bool {
        self.parameters_.borrow().stall_to_crossover
    }
    /// rusolve (HiGHS 1.16 backport, default off): Parameters::h16_ipx_scaled
    pub fn h16_ipx_scaled(&self) -> bool {
        self.parameters_.borrow().h16_ipx_scaled
    }
    /// rusolve (HiGHS 1.16 backport, default off): Parameters::h16_ipm_gap_offset
    pub fn h16_ipm_gap_offset(&self) -> bool {
        self.parameters_.borrow().h16_ipm_gap_offset
    }

    /// HiGHS: ipx/control.h::Control::setTimerOffset
    pub fn set_timer_offset(&self, offset: f64) {
        let mut t = self.timer_.get();
        t.offset_ = offset;
        self.timer_.set(t);
    }

    /// HiGHS: ipx/control.cc::Control::parameters() const (line 89)
    pub fn parameters(&self) -> std::cell::Ref<'_, Parameters> {
        self.parameters_.borrow()
    }

    /// HiGHS: ipx/control.cc::Control::parameters(const Parameters&) (lines 93-96)
    pub fn set_parameters(&self, new_parameters: &Parameters) {
        *self.parameters_.borrow_mut() = new_parameters.clone();
        self.make_stream();
    }

    /// HiGHS: ipx/control.cc::Control::callback (line 98)
    pub fn callback(&self, callback: Option<HighsCallback>) {
        *self.callback_.borrow_mut() = callback;
    }

    /// PORT NOTE: hands the callback moved in by [`Control::callback`] back
    /// to the caller (C++ never moved it).
    pub fn take_callback(&self) -> Option<HighsCallback> {
        self.callback_.borrow_mut().take()
    }

    /// HiGHS: ipx/control.cc::Control::OpenLogfile (lines 102-108)
    pub fn open_logfile(&self) {
        // logfile_.close()
        *self.logfile_.borrow_mut() = None;
        if let Some(filename) = self.parameters_.borrow().logfile.as_deref() {
            if !filename.is_empty() {
                // "Ignores if an error occurs; in this case no file log is written."
                *self.logfile_.borrow_mut() = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(filename)
                    .ok();
            }
        }
        self.make_stream();
    }

    /// HiGHS: ipx/control.cc::Control::CloseLogfile (lines 110-113)
    pub fn close_logfile(&self) {
        *self.logfile_.borrow_mut() = None;
        self.make_stream();
    }

    /// HiGHS: ipx/control.cc::Control::ResetTimer (line 115)
    pub fn reset_timer(&self) {
        let mut t = self.timer_.get();
        t.reset(false);
        self.timer_.set(t);
    }

    /// HiGHS: ipx/control.cc::Control::MakeStream (lines 119-125)
    fn make_stream(&self) {
        // output_.clear(); add cout if display; add logfile_ if open
        let to_logfile = self.logfile_.borrow().is_some();
        self.output_.set(Multistream {
            to_console: self.parameters_.borrow().display != 0,
            to_logfile,
        });
    }

    /// HiGHS: ipx/control.cc::Control::hLog, `highsLogUser(log_options_, kInfo, "%s", ..)` branch.
    ///
    /// PORT NOTE: C++ `assert(parameters_.log_options)`; a release build with
    /// no log options skips the message instead of dereferencing null.
    fn highs_log(&self, text: &str) {
        debug_assert!(self.parameters_.borrow().log_options.is_some());
        if let Some(log_options) = self.parameters_.borrow().log_options.as_ref() {
            crate::highs_log_user!(log_options, HighsLogType::kInfo, "{}", text);
        }
    }
}

impl Default for Control {
    fn default() -> Self {
        Self::new()
    }
}

/// HiGHS: ipx/control.cc::Format(Int, int) (lines 127-132): right-aligned in `width`.
pub fn format_int(i: Int, width: i32) -> String {
    c_pad_right(&i.to_string(), width.max(0) as usize)
}

/// HiGHS: ipx/control.cc::Format(const char*, int) (lines 134-139): right-aligned in `width`.
pub fn format_str(c: &str, width: i32) -> String {
    c_pad_right(c, width.max(0) as usize)
}

/// `std::ios_base::fmtflags floatfield` values used by IPX.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FloatField {
    scientific,
    fixed,
}

/// HiGHS: ipx/control.cc::Format(double, int, int, fmtflags) (lines 141-149)
pub fn format_double(d: f64, width: i32, prec: i32, floatfield: FloatField) -> String {
    let s = match floatfield {
        FloatField::scientific => c_fmt_e(d, prec as usize),
        FloatField::fixed => c_fmt_f(d, prec as usize),
    };
    c_pad_right(&s, width.max(0) as usize)
}

/// HiGHS: ipx/control.h::Scientific
pub fn scientific(d: f64, width: i32, prec: i32) -> String {
    format_double(d, width, prec, FloatField::scientific)
}
/// HiGHS: ipx/control.h::Fixed
pub fn fixed(d: f64, width: i32, prec: i32) -> String {
    format_double(d, width, prec, FloatField::fixed)
}
/// HiGHS: ipx/control.h::sci2
pub fn sci2(d: f64) -> String {
    scientific(d, 0, 2)
}
/// HiGHS: ipx/control.h::sci8
pub fn sci8(d: f64) -> String {
    scientific(d, 0, 8)
}
/// HiGHS: ipx/control.h::fix2
pub fn fix2(d: f64) -> String {
    fixed(d, 0, 2)
}
/// HiGHS: ipx/control.h::fix8
pub fn fix8(d: f64) -> String {
    fixed(d, 0, 8)
}
/// HiGHS: ipx/control.h::time (named `time_fmt`: `time` would shadow locals)
pub fn time_fmt(d: f64) -> String {
    fixed(d, 8, 1)
}
/// HiGHS: ipx/control.h::Textline: 4 spaces, then `text` left-justified in 52.
pub fn textline(text: &str) -> String {
    format!("    {}", c_pad_left(text, 52))
}
/// `os << d` with default ostream state (`%g`, precision 6), e.g.
/// `Debug(1) << x` or `h_logging << value`.
pub fn ostream_double(d: f64) -> String {
    c_fmt_g(d, 6)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ipx_format_helpers() {
        assert_eq!(format_int(42, 5), "   42");
        assert_eq!(sci2(12345.0), "1.23e+04");
        assert_eq!(fix2(1.005), "1.00");
        assert_eq!(time_fmt(3.0), "     3.0");
        assert_eq!(textline("abc").len(), 56);
        assert_eq!(ostream_double(0.1), "0.1");
    }
}
