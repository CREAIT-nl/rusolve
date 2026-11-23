//! HiGHS source: lp_data/HighsStatus.h, lp_data/HighsStatus.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use crate::prelude::*;

/// HiGHS: HighsStatus.h enum class HighsStatus
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsStatus {
    kError = -1,
    #[default]
    kOk = 0,
    kWarning = 1,
}

/// HiGHS: HighsStatus.cpp::highsStatusToString (lines 12-24, 13 lines)
pub fn highs_status_to_string(status: HighsStatus) -> String {
    match status {
        HighsStatus::kOk => "OK".to_string(),
        HighsStatus::kWarning => "Warning".to_string(),
        HighsStatus::kError => "Error".to_string(),
    }
}

/// HiGHS: HighsStatus.cpp::interpretCallStatus (lines 26-37, 12 lines)
pub fn interpret_call_status(
    log_options: &HighsLogOptions,
    call_status: HighsStatus,
    from_return_status: HighsStatus,
    message: &str,
) -> HighsStatus {
    let to_return_status: HighsStatus = worse_status(call_status, from_return_status);
    if call_status != HighsStatus::kOk {
        crate::highs_log_dev!(
            log_options,
            HighsLogType::kWarning,
            "{} return of HighsStatus::{}\n",
            message,
            highs_status_to_string(call_status)
        );
    }
    to_return_status
}

/// HiGHS: HighsStatus.cpp::worseStatus (lines 39-48, 10 lines)
pub fn worse_status(status0: HighsStatus, status1: HighsStatus) -> HighsStatus {
    if status0 == HighsStatus::kError || status1 == HighsStatus::kError {
        HighsStatus::kError
    } else if status0 == HighsStatus::kWarning || status1 == HighsStatus::kWarning {
        HighsStatus::kWarning
    } else {
        HighsStatus::kOk
    }
}
