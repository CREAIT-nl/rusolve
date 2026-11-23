//! HiGHS source: ipm/hipo/ipm/Status.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiPO status codes (complete). `enum Status` is an `Int` alias with one
//! constant per enumerator, in C++ declaration order, because the IPM
//! returns and compares them as integers (`if (Int status = ...)`).

use crate::ipm::hipo::Int;
use crate::ipm::ipx::{
    IPX_STATUS_dual_infeas, IPX_STATUS_failed, IPX_STATUS_imprecise, IPX_STATUS_iter_limit,
    IPX_STATUS_no_progress, IPX_STATUS_not_run, IPX_STATUS_optimal, IPX_STATUS_primal_infeas,
    IPX_STATUS_time_limit, IPX_STATUS_user_interrupt,
};

/// HiGHS: Status.h `enum Status`
pub type Status = Int;

pub const kStatusOk: Status = 0;
pub const kStatusNotRun: Status = 1;
pub const kStatusNoProgress: Status = 2;
pub const kStatusMaxIter: Status = 3;
pub const kStatusTimeLimit: Status = 4;
pub const kStatusUserInterrupt: Status = 5;
pub const kStatusFailed: Status = 6;
pub const kStatusError: Status = 7;
pub const kStatusOverflow: Status = 8;
pub const kStatusErrorAnalyse: Status = 9;
pub const kStatusErrorFactorise: Status = 10;
pub const kStatusErrorSolve: Status = 11;
pub const kStatusBadModel: Status = 12;
pub const kStatusUnknown: Status = 13;
pub const kStatusSolved: Status = 14;
pub const kStatusPrimalInfeasible: Status = 15;
pub const kStatusDualInfeasible: Status = 16;
pub const kStatusImprecise: Status = 17;
pub const kStatusPDFeas: Status = 18;
pub const kStatusBasic: Status = 19;

/// HiGHS: Status.h::IpxToHipoStatus
pub fn IpxToHipoStatus(ipx_status: Int) -> Status {
    match ipx_status {
        IPX_STATUS_not_run => kStatusNotRun,
        IPX_STATUS_iter_limit => kStatusMaxIter,
        IPX_STATUS_no_progress => kStatusNoProgress,
        IPX_STATUS_failed => kStatusError,
        IPX_STATUS_time_limit => kStatusTimeLimit,
        IPX_STATUS_user_interrupt => kStatusUserInterrupt,
        IPX_STATUS_primal_infeas => kStatusPrimalInfeasible,
        IPX_STATUS_dual_infeas => kStatusDualInfeasible,
        IPX_STATUS_optimal => kStatusPDFeas,
        IPX_STATUS_imprecise => kStatusImprecise,
        _ => kStatusUnknown,
    }
}

/// HiGHS: Status.h::statusString
pub fn statusString(status: Status) -> &'static str {
    match status {
        kStatusNotRun => "not run",
        kStatusMaxIter => "max iterations",
        kStatusNoProgress => "no progress",
        kStatusImprecise => "imprecise",
        kStatusError => "internal error",
        kStatusOverflow => "integer overflow",
        kStatusErrorAnalyse => "error in analyse phase",
        kStatusErrorFactorise => "error in factorise phase",
        kStatusErrorSolve => "error in solve phase",
        kStatusBadModel => "invalid model",
        kStatusTimeLimit => "time limit",
        kStatusUserInterrupt => "user interrupt",
        kStatusPrimalInfeasible => "primal infeasible",
        kStatusDualInfeasible => "dual infeasible",
        kStatusPDFeas => "primal-dual feasible",
        kStatusBasic => "crossover optimal",
        _ => "unknown",
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    #[ignore = "not yet filled"]
    fn status_maps() {
        assert_eq!(IpxToHipoStatus(IPX_STATUS_optimal), kStatusPDFeas);
        assert_eq!(IpxToHipoStatus(-123), kStatusUnknown);
        assert_eq!(statusString(kStatusBasic), "crossover optimal");
        assert_eq!(statusString(kStatusOk), "unknown");
    }
}
