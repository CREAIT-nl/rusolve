//! HiGHS source: qpsolver/qpconst.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete: enums only.

/// HiGHS: qpconst.hpp::QpSolverStatus (line 12)
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum QpSolverStatus {
    #[default]
    OK = 0,
    NOTPOSITIVDEFINITE,
    DEGENERATE,
    kError,
}

/// HiGHS: qpconst.hpp::QpModelStatus (lines 14-26)
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum QpModelStatus {
    #[default]
    kNotset = 0,
    kUndetermined,
    kOptimal,
    kUnbounded,
    kInfeasible,
    kIterationLimit,
    kTimeLimit,
    kLargeNullspace,
    kNonConvex,
    kInterrupt,
    kError,
}

/// HiGHS: qpconst.hpp::BasisStatus (lines 28-33)
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum BasisStatus {
    #[default]
    kInactive = 0,
    kActiveAtLower = 1,
    kActiveAtUpper,
    kInactiveInBasis,
}
