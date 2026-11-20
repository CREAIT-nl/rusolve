//! HiGHS source: ipm/hipo/ipm/Parameters.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! IPM parameters and the static regularisation (complete).

use crate::ipm::hipo::Int;

// parameters for termination
pub const kMaxIterDefault: Int = 200;
pub const kIpmTolDefault: f64 = 1e-8;
pub const kMaxBadIter: Int = 5;
pub const kDivergeTol: f64 = 1e3;

// parameters for correctors
pub const kGammaCorrector: f64 = 0.1;
pub const kSigmaAffine: f64 = 0.01;
pub const kMaxCorrectors: Int = 10;
pub const kMccIncreaseAlpha: f64 = 0.1;
pub const kMccIncreaseMin: f64 = 0.1;
pub const kSmallProduct: f64 = 1e-3;
pub const kLargeProduct: f64 = 1e3;

// parameters for choice of AS or NE
pub const kSpopsWeight: f64 = 30.0;
pub const kRatioOpsThresh: f64 = 10.0;
pub const kRatioSnThresh: f64 = 1.5;
pub const kSymbNzMult: f64 = 5.0;

// parameters for choice of parallelism
pub const kLargeFlopsThresh: f64 = 1e7;
pub const kLargeSpeedupThresh: f64 = 1.0;
pub const kLargeSnThresh: f64 = 20.0;
pub const kSmallSnThresh: f64 = 5.0;
pub const kMinNumberSn: Int = 10;
pub const kLargeStorageGB: f64 = 20.0;
pub const kLargeFillin: f64 = 50.0;
pub const kMaxTreeDepth: f64 = 1000.0;

// parameters for choice of ordering
pub const kFlopsOrderingThresh: f64 = 1.2;

// parameters for skipping AS or NE
pub const kNzBoundsRatio: f64 = 50.0;

// parameters for iterative refinement
pub const kMaxIterRefine: Int = 3;
pub const kTolRefine: f64 = 1e-12;

// parameters for scaling
pub const kSmallScalingCoeff: f64 = 1e-4;
pub const kLargeScalingCoeff: f64 = 1e4;
pub const kSmallBoundDiff: f64 = 1e-3;

// parameters for free variables
pub const kFreeVarsInitialBound: f64 = 1e4;
pub const kFreeVarsCloseRatio: f64 = 0.5;

/// HiGHS: Parameters.h `struct Regularisation`
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Regularisation {
    pub primal: f64,
    pub dual: f64,
}

impl Default for Regularisation {
    fn default() -> Self {
        Regularisation {
            primal: 1e-12,
            dual: 1e-10,
        }
    }
}
