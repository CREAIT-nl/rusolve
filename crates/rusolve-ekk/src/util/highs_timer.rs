//! HiGHS source: util/HighsTimer.h (stub)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Out of scope: timing. A wall-clock stub with the only behaviour the simplex path
//! observes: `read()` = seconds since construction/`start()`, used for time-limit tests
//! (HEkk::returnFromSolve, HFactor::buildKernel). Clock ids are accepted and ignored.

use crate::prelude::*;
use std::time::Instant;

/// HiGHS: HighsTimer.h::HighsTimerClock (stub; only passed around as `Option<&HighsTimerClock>`)
#[derive(Clone, Copy, Debug, Default)]
pub struct HighsTimerClock;

/// HiGHS: HighsTimer.h::HighsTimer (stub)
#[derive(Clone, Debug)]
pub struct HighsTimer {
    pub start_time: Instant,
    pub running: bool,
    /// Seconds banked by earlier start/stop intervals (HiGHS `clock_time[0]`).
    pub accumulated: f64,
    /// True between a `stop()` and the next `start()`: the clock is frozen at `accumulated`.
    pub stopped: bool,
}

impl Default for HighsTimer {
    fn default() -> Self {
        Self::new()
    }
}

impl HighsTimer {
    /// HiGHS: HighsTimer.h::HighsTimer::HighsTimer
    pub fn new() -> Self {
        Self {
            start_time: Instant::now(),
            running: false,
            accumulated: 0.0,
            // HiGHS defines the run clock stopped at 0 until `start()`.
            stopped: true,
        }
    }

    /// HiGHS: HighsTimer.h::HighsTimer::start (stub: `i_clock` ignored). The run clock
    /// accumulates across start/stop pairs like HiGHS's `clock_time[0]`: restarting banks
    /// the time already counted instead of resetting it, so a long-lived LP `Highs`
    /// (the MIP's relaxation solver) keeps one cumulative run time and its LP time limit
    /// (`getRunTime() + limit - mip_timer`) stays consistent with the clock HEkk reads.
    pub fn start(&mut self) {
        if self.running {
            return;
        }
        self.start_time = Instant::now();
        self.running = true;
        self.stopped = false;
    }

    /// HiGHS: HighsTimer.h::HighsTimer::stop (stub): freezes the clock at its current reading.
    pub fn stop(&mut self) {
        if self.running {
            self.accumulated += self.start_time.elapsed().as_secs_f64();
            self.stopped = true;
        }
        self.running = false;
    }

    /// HiGHS: HighsTimer.h::HighsTimer::read (stub: banked seconds plus the open interval)
    pub fn read(&self) -> f64 {
        if self.stopped {
            self.accumulated
        } else {
            self.accumulated + self.start_time.elapsed().as_secs_f64()
        }
    }

    /// HiGHS: HighsTimer.h::HighsTimer::running (stub: `i_clock` ignored; the
    /// run clock state)
    pub fn running(&self) -> bool {
        self.running
    }

    /// HiGHS: HighsTimer.h::HighsTimer::getWallTime (stub)
    pub fn get_wall_time(&self) -> f64 {
        self.read()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread::sleep;
    use std::time::Duration;

    /// A long-lived LP `Highs` (the MIP relaxation solver) computes its LP time limit as
    /// `getRunTime() + limit - mip_timer`; the clock HEkk reads must be the same cumulative
    /// clock. Before the fix `start()` reset the clock to 0 on every run, so an LP solve began
    /// with a fresh full limit and the root of a MIP overran its time limit by the setup time.
    #[test]
    fn run_clock_is_zero_until_started_and_cumulative_across_runs() {
        let mut t = HighsTimer::new();
        sleep(Duration::from_millis(30));
        assert_eq!(t.read(), 0.0, "an unstarted run clock reads 0");
        t.start();
        sleep(Duration::from_millis(30));
        t.stop();
        let first = t.read();
        assert!(first >= 0.03, "{first}");
        sleep(Duration::from_millis(30));
        assert_eq!(t.read(), first, "a stopped clock is frozen");
        t.start();
        let at_restart = t.read();
        assert!(
            (at_restart - first).abs() < 0.02,
            "restart must not reset: {first} -> {at_restart}"
        );
        sleep(Duration::from_millis(30));
        assert!(t.read() >= first + 0.03);
        let before = t.read();
        t.start(); // starting a running clock is a no-op
        assert!(t.read() >= before);
    }
}
