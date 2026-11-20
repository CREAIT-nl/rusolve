//! HiGHS source: parallel/HighsSchedulerConstants.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete. The spin / sleep thresholds of the scheduler. They only decide when a waiting
//! thread stops spinning and blocks; no result depends on them.

/// HiGHS: HighsSchedulerConstants.h::HighsSchedulerConstants::kNumTryFac (line 13)
pub const kNumTryFac: i32 = 16;
/// HiGHS: HighsSchedulerConstants.h::HighsSchedulerConstants::kMicroSecsBeforeSleep (line 14)
pub const kMicroSecsBeforeSleep: i64 = 5000;
/// HiGHS: HighsSchedulerConstants.h::HighsSchedulerConstants::kMicroSecsBeforeGlobalSync (line 15)
pub const kMicroSecsBeforeGlobalSync: i64 = 1000;
