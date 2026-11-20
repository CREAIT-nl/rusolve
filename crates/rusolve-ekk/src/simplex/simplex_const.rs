//! HiGHS source: simplex/SimplexConst.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete translation (constants are data).

use crate::HighsInt;
use crate::lp_data::hconst::*;

/// HiGHS: SimplexConst.h enum class SimplexAlgorithm
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum SimplexAlgorithm {
    #[default]
    kNone = 0,
    kPrimal = 1,
    kDual = 2,
}

// HiGHS: SimplexConst.h enum SimplexStrategy
pub const kSimplexStrategyMin: HighsInt = 0;
pub const kSimplexStrategyChoose: HighsInt = 0;
pub const kSimplexStrategyDual: HighsInt = 1;
pub const kSimplexStrategyDualPlain: HighsInt = 1;
pub const kSimplexStrategyDualTasks: HighsInt = 2;
pub const kSimplexStrategyDualMulti: HighsInt = 3;
pub const kSimplexStrategyPrimal: HighsInt = 4;
pub const kSimplexStrategyMax: HighsInt = 4;
pub const kSimplexStrategyNum: HighsInt = 5;

// HiGHS: SimplexConst.h enum SimplexUnscaledSolutionStrategy
pub const kSimplexUnscaledSolutionStrategyMin: HighsInt = 0;
pub const kSimplexUnscaledSolutionStrategyNone: HighsInt = 0;
pub const kSimplexUnscaledSolutionStrategyRefine: HighsInt = 1;
pub const kSimplexUnscaledSolutionStrategyDirect: HighsInt = 2;
pub const kSimplexUnscaledSolutionStrategyMax: HighsInt = 2;
pub const kSimplexUnscaledSolutionStrategyNum: HighsInt = 3;

// HiGHS: SimplexConst.h enum SimplexSolvePhase
pub const kSolvePhaseMin: HighsInt = -3;
pub const kSolvePhaseError: HighsInt = -3;
pub const kSolvePhaseExit: HighsInt = -2;
pub const kSolvePhaseUnknown: HighsInt = -1;
pub const kSolvePhaseOptimal: HighsInt = 0;
pub const kSolvePhase1: HighsInt = 1;
pub const kSolvePhase2: HighsInt = 2;
pub const kSolvePhasePrimalInfeasibleCleanup: HighsInt = 3;
pub const kSolvePhaseOptimalCleanup: HighsInt = 4;
pub const kSolvePhaseTabooBasis: HighsInt = 5;
pub const kSolvePhaseMax: HighsInt = 5;

// HiGHS: SimplexConst.h enum SimplexCrashStrategy
pub const kSimplexCrashStrategyMin: HighsInt = 0;
pub const kSimplexCrashStrategyOff: HighsInt = 0;
pub const kSimplexCrashStrategyLtssfK: HighsInt = 1;
pub const kSimplexCrashStrategyLtssf: HighsInt = 1;
pub const kSimplexCrashStrategyBixby: HighsInt = 2;
pub const kSimplexCrashStrategyLtssfPri: HighsInt = 3;
pub const kSimplexCrashStrategyLtsfK: HighsInt = 4;
pub const kSimplexCrashStrategyLtsfPri: HighsInt = 5;
pub const kSimplexCrashStrategyLtsf: HighsInt = 6;
pub const kSimplexCrashStrategyBixbyNoNonzeroColCosts: HighsInt = 7;
pub const kSimplexCrashStrategyBasic: HighsInt = 8;
pub const kSimplexCrashStrategyTestSing: HighsInt = 9;
pub const kSimplexCrashStrategyMax: HighsInt = 9;

// HiGHS: SimplexConst.h enum SimplexEdgeWeightStrategy
pub const kSimplexEdgeWeightStrategyMin: HighsInt = -1;
pub const kSimplexEdgeWeightStrategyChoose: HighsInt = -1;
pub const kSimplexEdgeWeightStrategyDantzig: HighsInt = 0;
pub const kSimplexEdgeWeightStrategyDevex: HighsInt = 1;
pub const kSimplexEdgeWeightStrategySteepestEdge: HighsInt = 2;
pub const kSimplexEdgeWeightStrategyMax: HighsInt = 2;

// HiGHS: SimplexConst.h enum SimplexPriceStrategy
pub const kSimplexPriceStrategyMin: HighsInt = 0;
pub const kSimplexPriceStrategyCol: HighsInt = 0;
pub const kSimplexPriceStrategyRow: HighsInt = 1;
pub const kSimplexPriceStrategyRowSwitch: HighsInt = 2;
pub const kSimplexPriceStrategyRowSwitchColSwitch: HighsInt = 3;
pub const kSimplexPriceStrategyMax: HighsInt = 3;

// HiGHS: SimplexConst.h enum SimplexPivotalRowRefinementStrategy
pub const kSimplexInfeasibilityProofRefinementMin: HighsInt = 0;
pub const kSimplexInfeasibilityProofRefinementNo: HighsInt = 0;
pub const kSimplexInfeasibilityProofRefinementUnscaledLp: HighsInt = 1;
pub const kSimplexInfeasibilityProofRefinementAlsoScaledLp: HighsInt = 2;
pub const kSimplexInfeasibilityProofRefinementMax: HighsInt = 2;

// HiGHS: SimplexConst.h enum SimplexPrimalCorrectionStrategy
pub const kSimplexPrimalCorrectionStrategyNone: HighsInt = 0;
pub const kSimplexPrimalCorrectionStrategyInRebuild: HighsInt = 1;
pub const kSimplexPrimalCorrectionStrategyAlways: HighsInt = 2;

// HiGHS: SimplexConst.h enum RebuildReason
pub const kRebuildReasonCleanup: HighsInt = -1;
pub const kRebuildReasonNo: HighsInt = 0;
pub const kRebuildReasonUpdateLimitReached: HighsInt = 1;
pub const kRebuildReasonSyntheticClockSaysInvert: HighsInt = 2;
pub const kRebuildReasonPossiblyOptimal: HighsInt = 3;
pub const kRebuildReasonPossiblyPhase1Feasible: HighsInt = 4;
pub const kRebuildReasonPossiblyPrimalUnbounded: HighsInt = 5;
pub const kRebuildReasonPossiblyDualUnbounded: HighsInt = 6;
pub const kRebuildReasonPossiblySingularBasis: HighsInt = 7;
pub const kRebuildReasonPrimalInfeasibleInPrimalSimplex: HighsInt = 8;
pub const kRebuildReasonChooseColumnFail: HighsInt = 9;
pub const kRebuildReasonForceRefactor: HighsInt = 10;
pub const kRebuildReasonExcessivePrimalValue: HighsInt = 11;
pub const kRebuildReasonCount: HighsInt = 12;

// HiGHS: SimplexConst.h enum SimplexNlaOperation
pub const kSimplexNlaNull: HighsInt = -1;
pub const kSimplexNlaBtranFull: HighsInt = 0;
pub const kSimplexNlaPriceFull: HighsInt = 1;
pub const kSimplexNlaBtranBasicFeasibilityChange: HighsInt = 2;
pub const kSimplexNlaPriceBasicFeasibilityChange: HighsInt = 3;
pub const kSimplexNlaBtranEp: HighsInt = 4;
pub const kSimplexNlaPriceAp: HighsInt = 5;
pub const kSimplexNlaFtran: HighsInt = 6;
pub const kSimplexNlaFtranBfrt: HighsInt = 7;
pub const kSimplexNlaFtranDse: HighsInt = 8;
pub const kSimplexNlaBtranPse: HighsInt = 9;
pub const kNumSimplexNlaOperation: HighsInt = 10;

/// HiGHS: SimplexConst.h enum class EdgeWeightMode
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum EdgeWeightMode {
    #[default]
    kDantzig = 0,
    kDevex = 1,
    kSteepestEdge = 2,
    kCount = 3,
}

pub const kDualTasksMinConcurrency: HighsInt = 3;
pub const kDualMultiMinConcurrency: HighsInt = 1;
pub const kNonbasicFlagTrue: i8 = 1;
pub const kNonbasicFlagFalse: i8 = 0;
pub const kIllegalFlagValue: i8 = -99;
pub const kNonbasicMoveUp: i8 = 1;
pub const kNonbasicMoveDn: i8 = -1;
pub const kNonbasicMoveZe: i8 = 0;
pub const kIllegalMoveValue: i8 = -99;
pub const kAcceptDseWeightThreshold: f64 = 0.25;
pub const kMinDualSteepestEdgeWeight: f64 = 1e-4;
pub const kNoRowSought: HighsInt = -2;
pub const kNoRowChosen: HighsInt = -1;
pub const kNoRayIndex: HighsInt = -1;
pub const kNoRaySign: HighsInt = 0;
pub const kDebugMipNodeDualFeasible: bool = false;
/// HiGHS: SimplexConst.h enum class LpAction
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum LpAction {
    #[default]
    kScale = 0,
    kNewCosts = 1,
    kNewBounds = 2,
    kNewBasis = 3,
    kNewCols = 4,
    kNewRows = 5,
    kDelCols = 6,
    kDelNonbasicCols = 7,
    kDelRows = 8,
    kDelRowsBasisOk = 9,
    kScaledCol = 10,
    kScaledRow = 11,
    kBacktracking = 12,
}

/// HiGHS: SimplexConst.h enum class BadBasisChangeReason
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum BadBasisChangeReason {
    #[default]
    kAll = 0,
    kSingular = 1,
    kCycling = 2,
    kFailedInfeasibilityProof = 3,
}

pub const kAllowedNumBadDevexWeight: HighsInt = 3;
pub const kBadDevexWeightFactor: f64 = 3.0;

/// HiGHS: HEkkDual.h kHighsSlicedLimit
pub const kHighsSlicedLimit: HighsInt = kSimplexConcurrencyLimit;
