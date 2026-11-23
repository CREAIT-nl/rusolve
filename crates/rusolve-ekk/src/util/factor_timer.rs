//! HiGHS source: util/FactorTimer.h (clock ids + no-op FactorTimer)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Out of scope: timing. Clock ids plus a `FactorTimer` whose start/stop do nothing.

use crate::HighsInt;
use crate::util::highs_timer::HighsTimerClock;

// HiGHS: FactorTimer.h enum iClockFactor
pub const FactorInvert: HighsInt = 0;
pub const FactorInvertSimple: HighsInt = 1;
pub const FactorInvertKernel: HighsInt = 2;
pub const FactorInvertDeficient: HighsInt = 3;
pub const FactorInvertFinish: HighsInt = 4;
pub const FactorFtran: HighsInt = 5;
pub const FactorFtranLower: HighsInt = 6;
pub const FactorFtranLowerAPF: HighsInt = 7;
pub const FactorFtranLowerDse: HighsInt = 8;
pub const FactorFtranLowerSps: HighsInt = 9;
pub const FactorFtranLowerHyper: HighsInt = 10;
pub const FactorFtranUpper: HighsInt = 11;
pub const FactorFtranUpperFT: HighsInt = 12;
pub const FactorFtranUpperMPF: HighsInt = 13;
pub const FactorFtranUpperDse: HighsInt = 14;
pub const FactorFtranUpperSps0: HighsInt = 15;
pub const FactorFtranUpperSps1: HighsInt = 16;
pub const FactorFtranUpperSps2: HighsInt = 17;
pub const FactorFtranUpperHyper0: HighsInt = 18;
pub const FactorFtranUpperHyper1: HighsInt = 19;
pub const FactorFtranUpperHyper2: HighsInt = 20;
pub const FactorFtranUpperHyper3: HighsInt = 21;
pub const FactorFtranUpperHyper4: HighsInt = 22;
pub const FactorFtranUpperHyper5: HighsInt = 23;
pub const FactorFtranUpperPF: HighsInt = 24;
pub const FactorBtran: HighsInt = 25;
pub const FactorBtranLower: HighsInt = 26;
pub const FactorBtranLowerDse: HighsInt = 27;
pub const FactorBtranLowerSps: HighsInt = 28;
pub const FactorBtranLowerHyper: HighsInt = 29;
pub const FactorBtranLowerAPF: HighsInt = 30;
pub const FactorBtranUpper: HighsInt = 31;
pub const FactorBtranUpperPF: HighsInt = 32;
pub const FactorBtranUpperDse: HighsInt = 33;
pub const FactorBtranUpperSps: HighsInt = 34;
pub const FactorBtranUpperHyper: HighsInt = 35;
pub const FactorBtranUpperFT: HighsInt = 36;
pub const FactorBtranUpperMPF: HighsInt = 37;
pub const FactorReinvert: HighsInt = 38;
pub const FactorNumClock: HighsInt = 39;

/// HiGHS: FactorTimer.h::FactorTimer (stub)
#[derive(Clone, Copy, Debug, Default)]
pub struct FactorTimer;

impl FactorTimer {
    /// HiGHS: FactorTimer.h::FactorTimer::start (stub: no-op)
    pub fn start(
        &self,
        factor_clock: HighsInt,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
    }

    /// HiGHS: FactorTimer.h::FactorTimer::stop (stub: no-op)
    pub fn stop(
        &self,
        factor_clock: HighsInt,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
    }
}
