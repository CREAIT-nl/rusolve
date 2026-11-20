//! HiGHS source: ipm/hipo/ipm/Control.h, ipm/hipo/ipm/Control.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Interface to the HiGHS timer and callback. The C++ holds three pointers;
//! the Rust owns the callback (moved in by `setCallback`, handed on to IPX
//! by `take_callback` in `Solver::prepareIpx` and back to `solveHipo` at the
//! end), a clone of the timer (only `read()` is used) and a clone of the
//! options (only `time_limit` is read, and it is fixed once
//! `Solver::resetOptions` has run).

use super::options::Options;
use crate::ipm::hipo::Int;
use crate::lp_data::hconst::kCallbackIpmInterrupt;
use crate::lp_data::highs_callback::HighsCallback;
use crate::util::highs_timer::HighsTimer;

/// HiGHS: Control.h `class Control`
#[derive(Default)]
pub struct Control {
    pub callback_: Option<HighsCallback>,
    pub timer_: Option<HighsTimer>,
    pub options_: Option<Options>,
}

impl Control {
    /// HiGHS: Control.cpp::Control::setCallback (lines 10-10)
    pub fn setCallback(&mut self, callback: HighsCallback) {
        self.callback_ = Some(callback);
    }

    /// Rust-only: move the callback back out (`None` if it was never set).
    pub fn take_callback(&mut self) -> Option<HighsCallback> {
        self.callback_.take()
    }

    /// HiGHS: Control.cpp::Control::setTimer (lines 11-11)
    pub fn setTimer(&mut self, timer: &HighsTimer) {
        self.timer_ = Some(timer.clone());
    }

    /// HiGHS: Control.cpp::Control::setOptions (lines 12-12)
    pub fn setOptions(&mut self, options: &Options) {
        self.options_ = Some(options.clone());
    }

    /// HiGHS: Control.cpp::Control::callback (lines 14-14)
    pub fn callback(&mut self) -> Option<&mut HighsCallback> {
        self.callback_.as_mut()
    }

    /// HiGHS: Control.cpp::Control::elapsed (lines 16-16)
    pub fn elapsed(&self) -> f64 {
        match &self.timer_ {
            Some(t) => t.read(),
            None => -1.0,
        }
    }

    /// HiGHS: Control.cpp::Control::interruptCheck (lines 18-36)
    /// C++ default argument `ipm_iteration_count = -1`. `&mut self` because
    /// `callbackAction` writes the callback's output data.
    pub fn interruptCheck(&mut self, ipm_iteration_count: Int) -> Int {
        // C++ line 19: HighsTaskExecutor::getThisWorkerDeque()->checkInterrupt();
        let deque =
            crate::parallel::highs_task_executor::HighsTaskExecutor::get_this_worker_deque();
        // SAFETY: this thread's own worker deque (thread-local); checkInterrupt is
        // owner-only and this is the owner thread (C++ protocol).
        unsafe { (*deque).check_interrupt() };

        if let Some(ref options) = self.options_ {
            if options.time_limit >= 0.0 && self.elapsed() > options.time_limit {
                return crate::ipm::hipo::ipm::status::kStatusTimeLimit;
            }
        }

        if let Some(ref mut callback) = self.callback_ {
            if callback.user_callback.is_some() && callback.active[kCallbackIpmInterrupt as usize] {
                callback.clear_highs_callback_output();
                callback.data_out.ipm_iteration_count = ipm_iteration_count;
                if callback.callback_action(kCallbackIpmInterrupt, "IPM interrupt") {
                    return crate::ipm::hipo::ipm::status::kStatusUserInterrupt;
                }
            }
        }

        0
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn interrupt_check_without_timer_or_callback_is_zero() {
        let mut c = Control::default();
        assert_eq!(c.interruptCheck(-1), 0);
        assert_eq!(c.elapsed(), -1.0);
    }

    // Unit tests.

    #[test]
    fn w5_hipo_ctrl_threads_from_scheduler() {
        // The HiPO thread count (Solver.cpp:1201 log line, `setParallel` choice) comes from `highs::parallel::num_threads()` instead of the constant 1.
        std::thread::spawn(|| {
            crate::parallel::highs_parallel::initialize_scheduler(2);
            assert_eq!(crate::parallel::highs_parallel::num_threads(), 2);
            crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
        })
        .join()
        .unwrap();
    }
}
