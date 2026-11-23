//! HiGHS source: lp_data/HConst.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete translation of HConst.h (constants are data, not bodies).

use crate::HighsInt;

pub const kHighsCopyrightStatement: &str = "Copyright (c) 2026 under MIT licence terms";
pub const kHighsSize_tInf: usize = usize::MAX;
pub const kHighsIInf: HighsInt = HighsInt::MAX;
pub const kHighsIInf32: HighsInt = i32::MAX;
pub const kHighsInf: f64 = f64::INFINITY;
pub const kHighsUndefined: f64 = kHighsInf;
pub const kHighsTiny: f64 = 1e-14;
pub const kHighsMacheps: f64 = f64::EPSILON;
/// HiGHS 1.16 (#3267): lp_data/HConst.h lines 30-31 (1.16.0): 2^52, doubles of larger magnitude cannot represent fractional parts
pub const kHighsInvMacheps: f64 = 1.0 / kHighsMacheps;
pub const kHighsZero: f64 = 1e-50;
pub const kHighsOffString: &str = "off";
pub const kHighsChooseString: &str = "choose";
pub const kHighsOnString: &str = "on";
pub const kHighsMaxStringLength: HighsInt = 512;
pub const kSimplexConcurrencyLimit: HighsInt = 8;
pub const kRunningAverageMultiplier: f64 = 0.05;
pub const kExcessivelySmallObjectiveCoefficient: f64 = 1e-4;
pub const kExcessivelyLargeObjectiveCoefficient: f64 = 1e6;
pub const kExcessivelySmallBoundValue: f64 = 1e-4;
pub const kExcessivelyLargeBoundValue: f64 = 1e6;
/// HiGHS 1.16 (h16_scaling_report): lp_data/HConst.h line 46 (1.16.0)
pub const kExcessivelyLargeIntegerBoundValue: f64 = 1e5;
pub const kNoThreadInstance: HighsInt = -1;
pub const kAllowDeveloperAssert: bool = false;
pub const kExtendInvertWhenAddingRows: bool = false;
/// HiGHS: HConst.h enum class HighsLogType
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsLogType {
    #[default]
    kInfo = 1,
    kDetailed = 2,
    kVerbose = 3,
    kWarning = 4,
    kError = 5,
}

// HiGHS: HConst.h enum SimplexScaleStrategy
pub const kSimplexScaleStrategyMin: HighsInt = 0;
pub const kSimplexScaleStrategyOff: HighsInt = 0;
pub const kSimplexScaleStrategyChoose: HighsInt = 1;
pub const kSimplexScaleStrategyEquilibration: HighsInt = 2;
pub const kSimplexScaleStrategyForcedEquilibration: HighsInt = 3;
pub const kSimplexScaleStrategyMaxValue: HighsInt = 4;
pub const kSimplexScaleStrategyMaxValue015: HighsInt = 4;
pub const kSimplexScaleStrategyMaxValue0157: HighsInt = 4;
pub const kSimplexScaleStrategyMax: HighsInt = 4;

// HiGHS: HConst.h enum HighsDebugLevel
pub const kHighsDebugLevelNone: HighsInt = 0;
pub const kHighsDebugLevelCheap: HighsInt = 1;
pub const kHighsDebugLevelCostly: HighsInt = 2;
pub const kHighsDebugLevelExpensive: HighsInt = 3;
pub const kHighsDebugLevelMin: HighsInt = 0;
pub const kHighsDebugLevelMax: HighsInt = 3;

/// HiGHS: HConst.h enum class HighsDebugStatus
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsDebugStatus {
    kNotChecked = -1,
    #[default]
    kOk = 0,
    kSmallError = 1,
    kWarning = 2,
    kLargeError = 3,
    kError = 4,
    kExcessiveError = 5,
    kLogicalError = 6,
}

// HiGHS: HConst.h enum HighsAnalysisLevel
pub const kHighsAnalysisLevelNone: HighsInt = 0;
pub const kHighsAnalysisLevelModelData: HighsInt = 1;
pub const kHighsAnalysisLevelSolverSummaryData: HighsInt = 2;
pub const kHighsAnalysisLevelSolverRuntimeData: HighsInt = 4;
pub const kHighsAnalysisLevelSolverTime: HighsInt = 8;
pub const kHighsAnalysisLevelNlaData: HighsInt = 16;
pub const kHighsAnalysisLevelNlaTime: HighsInt = 32;
pub const kHighsAnalysisLevelMipData: HighsInt = 64;
pub const kHighsAnalysisLevelMipTime: HighsInt = 128;
pub const kHighsAnalysisLevelPresolveTime: HighsInt = 256;
pub const kHighsAnalysisLevelMin: HighsInt = 0;
pub const kHighsAnalysisLevelMax: HighsInt = 511;

/// HiGHS: HConst.h enum class HighsVarType
#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsVarType {
    #[default]
    kContinuous = 0,
    kInteger = 1,
    kSemiContinuous = 2,
    kSemiInteger = 3,
    kImplicitInteger = 4,
}

/// HiGHS: HConst.h enum class HighsOptionType
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsOptionType {
    #[default]
    kBool = 0,
    kInt = 1,
    kDouble = 2,
    kString = 3,
}

/// HiGHS: HConst.h enum class HighsInfoType
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsInfoType {
    #[default]
    kInt64 = -1,
    kInt = 1,
    kDouble = 2,
}

/// HiGHS: HConst.h enum class HighsRunDataType
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsRunDataType {
    #[default]
    kInt64 = -1,
    kInt = 1,
    kDouble = 2,
}

// HiGHS: HConst.h enum OptionOffChooseOn
pub const kHighsOptionOff: HighsInt = -1;
pub const kHighsOptionChoose: HighsInt = 0;
pub const kHighsOptionOn: HighsInt = 1;

// HiGHS: HConst.h enum IpxDualizeStrategy
pub const kIpxDualizeStrategyOff: HighsInt = -1;
pub const kIpxDualizeStrategyChoose: HighsInt = 0;
pub const kIpxDualizeStrategyOn: HighsInt = 1;
pub const kIpxDualizeStrategyLukas: HighsInt = 2;
pub const kIpxDualizeStrategyFilippo: HighsInt = 3;
pub const kIpxDualizeStrategyMin: HighsInt = -1;
pub const kIpxDualizeStrategyMax: HighsInt = 3;

/// HiGHS: HConst.h enum class ObjSense
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum ObjSense {
    #[default]
    kMinimize = 1,
    kMaximize = -1,
}

/// HiGHS: HConst.h enum class MatrixFormat
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum MatrixFormat {
    #[default]
    kColwise = 1,
    kRowwise = 2,
    kRowwisePartitioned = 3,
}

/// HiGHS: HConst.h enum class HessianFormat
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HessianFormat {
    #[default]
    kTriangular = 1,
    kSquare = 2,
}

// HiGHS: HConst.h enum SolutionStatus
pub const kSolutionStatusNone: HighsInt = 0;
pub const kSolutionStatusInfeasible: HighsInt = 1;
pub const kSolutionStatusFeasible: HighsInt = 2;
pub const kSolutionStatusMin: HighsInt = 0;
pub const kSolutionStatusMax: HighsInt = 2;

// HiGHS: HConst.h enum BasisValidity
pub const kBasisValidityInvalid: HighsInt = 0;
pub const kBasisValidityValid: HighsInt = 1;
pub const kBasisValidityMin: HighsInt = 0;
pub const kBasisValidityMax: HighsInt = 1;

// HiGHS 1.16 #3456 (h16_ac_centring): HConst.h enum CentringStatus (1.16.0, appended; constants only)
pub const kCentringStatusNone: HighsInt = 0;
pub const kCentringStatusFailed: HighsInt = 1;
pub const kCentringStatusSuccess: HighsInt = 2;
pub const kCentringStatusMin: HighsInt = 0;
pub const kCentringStatusMax: HighsInt = 2;

pub const kHighsBasisFileV1: &str = "v1";
pub const kHighsBasisFileV2: &str = "v2";
// HiGHS: HConst.h enum SolutionStyle
pub const kSolutionStyleOldRaw: HighsInt = -1;
pub const kSolutionStyleRaw: HighsInt = 0;
pub const kSolutionStylePretty: HighsInt = 1;
pub const kSolutionStyleGlpsolRaw: HighsInt = 2;
pub const kSolutionStyleGlpsolPretty: HighsInt = 3;
pub const kSolutionStyleSparse: HighsInt = 4;
pub const kSolutionStyleMin: HighsInt = -1;
pub const kSolutionStyleMax: HighsInt = 4;

// HiGHS: HConst.h enum GlpsolCostRowLocation
pub const kGlpsolCostRowLocationLast: HighsInt = -2;
pub const kGlpsolCostRowLocationNone: HighsInt = -1;
pub const kGlpsolCostRowLocationNoneIfEmpty: HighsInt = 0;
pub const kGlpsolCostRowLocationMin: HighsInt = -2;

pub const kHighsFilenameDefault: &str = "";
pub const kHighsMinimalColNamePrefix: &str = "c";
pub const kHighsMinimalrowNamePrefix: &str = "r";
pub const kHighsUniqueColNamePrefix: &str = "c_ekk";
pub const kHighsUniquerowNamePrefix: &str = "r_ekk";
/// HiGHS: HConst.h enum class HighsPresolveStatus
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsPresolveStatus {
    kNotPresolved = -1,
    #[default]
    kNotReduced = 0,
    kInfeasible = 1,
    kUnboundedOrInfeasible = 2,
    kReduced = 3,
    kReducedToEmpty = 4,
    kTimeout = 5,
    kNullError = 6,
    kOptionsError = 7,
    kNotSet = 8,
    kOutOfMemory = 9,
}

/// HiGHS: HConst.h enum class HighsPostsolveStatus
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsPostsolveStatus {
    kNotPresolved = -1,
    #[default]
    kNoPrimalSolutionError = 0,
    kSolutionRecovered = 1,
    kBasisError = 2,
}

/// HiGHS: HConst.h enum class HighsModelStatus
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsModelStatus {
    #[default]
    kNotset = 0,
    kLoadError = 1,
    kModelError = 2,
    kPresolveError = 3,
    kSolveError = 4,
    kPostsolveError = 5,
    kModelEmpty = 6,
    kOptimal = 7,
    kInfeasible = 8,
    kUnboundedOrInfeasible = 9,
    kUnbounded = 10,
    kObjectiveBound = 11,
    kObjectiveTarget = 12,
    kTimeLimit = 13,
    kIterationLimit = 14,
    kUnknown = 15,
    kSolutionLimit = 16,
    kInterrupt = 17,
    kMemoryLimit = 18,
    kHighsInterrupt = 19,
}

/// HiGHS: enumerators of `HighsModelStatus` that repeat an earlier value.
impl HighsModelStatus {
    pub const kMin: HighsModelStatus = HighsModelStatus::kNotset;
    pub const kMax: HighsModelStatus = HighsModelStatus::kHighsInterrupt;
}

// HiGHS: HConst.h enum HighsCallbackType
pub const kCallbackMin: HighsInt = 0;
pub const kCallbackLogging: HighsInt = 0;
pub const kCallbackSimplexInterrupt: HighsInt = 1;
pub const kCallbackIpmInterrupt: HighsInt = 2;
pub const kCallbackMipSolution: HighsInt = 3;
pub const kCallbackMipImprovingSolution: HighsInt = 4;
pub const kCallbackMipLogging: HighsInt = 5;
pub const kCallbackMipInterrupt: HighsInt = 6;
pub const kCallbackMipGetCutPool: HighsInt = 7;
pub const kCallbackMipDefineLazyConstraints: HighsInt = 8;
pub const kCallbackMipUserSolution: HighsInt = 9;
pub const kCallbackMax: HighsInt = 9;
pub const kNumCallbackType: HighsInt = 10;

/// HiGHS: HConst.h enum class HighsBasisStatus
#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsBasisStatus {
    #[default]
    kLower = 0,
    kBasic = 1,
    kUpper = 2,
    kZero = 3,
    kNonbasic = 4,
}

// HiGHS: HConst.h enum PresolveRuleType
pub const kPresolveRuleIllegal: HighsInt = -1;
pub const kPresolveRuleMin: HighsInt = 0;
pub const kPresolveRuleEmptyRow: HighsInt = 0;
pub const kPresolveRuleSingletonRow: HighsInt = 1;
pub const kPresolveRuleRedundantRow: HighsInt = 2;
pub const kPresolveRuleEmptyCol: HighsInt = 3;
pub const kPresolveRuleFixedCol: HighsInt = 4;
pub const kPresolveRuleDominatedCol: HighsInt = 5;
pub const kPresolveRuleFirstAllowOff: HighsInt = 6;
pub const kPresolveRuleForcingRow: HighsInt = 6;
pub const kPresolveRuleForcingCol: HighsInt = 7;
pub const kPresolveRuleFreeColSubstitution: HighsInt = 8;
pub const kPresolveRuleDoubletonEquation: HighsInt = 9;
pub const kPresolveRuleDependentEquations: HighsInt = 10;
pub const kPresolveRuleDependentFreeCols: HighsInt = 11;
pub const kPresolveRuleAggregator: HighsInt = 12;
pub const kPresolveRuleParallelRowsAndCols: HighsInt = 13;
pub const kPresolveRuleSparsify: HighsInt = 14;
pub const kPresolveRuleProbing: HighsInt = 15;
pub const kPresolveRuleEnumeration: HighsInt = 16;
pub const kPresolveRuleDualFixing: HighsInt = 17;
pub const kPresolveRuleColStuffing: HighsInt = 18;
pub const kPresolveRuleInitialSweep: HighsInt = 19;
pub const kPresolveRuleMax: HighsInt = 19;
pub const kPresolveRuleLastAllowOff: HighsInt = 19;
pub const kPresolveRuleCount: HighsInt = 20;

// HiGHS 1.16 (#3536): HConst.h enum PresolveMode (1.16.0 lines 278-285), read only under
// h16_presolve_rules.
pub const kPresolveModeMin: HighsInt = 0;
pub const kPresolveModeDefault: HighsInt = 0;
pub const kPresolveModeInitialSweep: HighsInt = 1;
pub const kPresolveModeLight: HighsInt = 2;
pub const kPresolveModeFull: HighsInt = 3;
pub const kPresolveModeMax: HighsInt = 3;

// HiGHS 1.16 (#3167 #3218 #3329): HConst.h enum PresolveRuleType (1.16.0 lines 288-318)
// adds ZeroCostSingleton, DualFixProbing, FourierMotzkin and WeaklyDominatedCol and
// renumbers ColStuffing / InitialSweep. The port keeps the 1.15.1 numbering above (rule
// F) and, under h16_presolve_rules, appends the new rules as extra allow_rule_ slots;
// h16_presolve_rule_from_116 maps a 1.16 rule index (a presolve_rule_off bit) to its slot.
pub const kH16PresolveRuleZeroCostSingleton: HighsInt = 20;
pub const kH16PresolveRuleDualFixProbing: HighsInt = 21;
pub const kH16PresolveRuleFourierMotzkin: HighsInt = 22;
pub const kH16PresolveRuleWeaklyDominatedCol: HighsInt = 23;
pub const kH16PresolveRuleCount: HighsInt = 24;

/// 1.16 rule index (HConst.h 1.16.0 lines 290-314) -> port allow_rule_ slot.
pub fn h16_presolve_rule_from_116(rule: HighsInt) -> HighsInt {
    match rule {
        18 => kH16PresolveRuleZeroCostSingleton,
        19 => kPresolveRuleColStuffing,
        20 => kH16PresolveRuleDualFixProbing,
        21 => kPresolveRuleInitialSweep,
        22 => kH16PresolveRuleFourierMotzkin,
        23 => kH16PresolveRuleWeaklyDominatedCol,
        r => r,
    }
}

// HiGHS: HConst.h enum IisStrategy
pub const kIisStrategyMin: HighsInt = 0;
pub const kIisStrategyLight: HighsInt = 0;
pub const kIisStrategyFromRay: HighsInt = 1;
pub const kIisStrategyFromLp: HighsInt = 2;
pub const kIisStrategyIrreducible: HighsInt = 4;
pub const kIisStrategyColPriority: HighsInt = 8;
pub const kIisStrategyRelaxation: HighsInt = 16;
pub const kIisStrategyDefault: HighsInt = 0;
pub const kIisStrategyMax: HighsInt = 31;

// HiGHS: HConst.h enum IisStatus
pub const kIisStatusMin: HighsInt = -1;
pub const kIisStatusNotInConflict: HighsInt = -1;
pub const kIisStatusMaybeInConflict: HighsInt = 0;
pub const kIisStatusInConflict: HighsInt = 1;
pub const kIisStatusMax: HighsInt = 1;

// HiGHS: HConst.h enum MipChooseSubMipRecord
pub const kMipRecord: HighsInt = -1;
pub const kChooseRecord: HighsInt = 0;
pub const kSubMipRecord: HighsInt = 1;

// HiGHS: HConst.h enum PresolveSolvePostsolveIndex
pub const kPresolveTime: HighsInt = 0;
pub const kSolveTime: HighsInt = 1;
pub const kPostsolveTime: HighsInt = 2;
pub const kToPresolveSolvePostsolve: HighsInt = 3;

// HiGHS: HConst.h enum SubSolverIndex
pub const kFromSubSolver: HighsInt = 3;
pub const kSubSolverMip: HighsInt = 3;
pub const kSubSolverDuSimplexBasis: HighsInt = 4;
pub const kSubSolverDuSimplexNoBasis: HighsInt = 5;
pub const kSubSolverPrSimplexBasis: HighsInt = 6;
pub const kSubSolverPrSimplexNoBasis: HighsInt = 7;
pub const kSubSolverHipo: HighsInt = 8;
pub const kSubSolverIpx: HighsInt = 9;
pub const kSubSolverHipoAc: HighsInt = 10;
pub const kSubSolverIpxAc: HighsInt = 11;
pub const kSubSolverPdlp: HighsInt = 12;
pub const kSubSolverQpAsm: HighsInt = 13;
pub const kSubSolverSubMip: HighsInt = 14;
pub const kLastSubSolver: HighsInt = 14;
pub const kToSubSolver: HighsInt = 15;

pub const kMinimumKktTolerance: f64 = 1e-10;
pub const kDefaultKktTolerance: f64 = 1e-7;
pub const kMinimumMipTolerance: f64 = 1e-10;
pub const kDefaultMipTolerance: f64 = 1e-6;
pub const kMinimumIpmTolerance: f64 = 1e-12;
pub const kDefaultIpmTolerance: f64 = 1e-1 * kDefaultKktTolerance;
pub const kHessianRegularizationValue: f64 = 1e-7;
pub const kDefaultAllowedMatrixPow2Scale: HighsInt = 20;
pub const kMaxAllowedMatrixPow2Scale: HighsInt = 30;
pub const kHighsIllegalInfeasibilityMeasure: f64 = kHighsInf;
pub const kHighsIllegalInfeasibilityCount: HighsInt = -1;
pub const kHighsIllegalResidualMeasure: f64 = kHighsInf;
pub const kHighsIllegalResidualCount: HighsInt = -1;
pub const kHighsIllegalErrorValue: f64 = kHighsInf;
pub const kHighsIllegalErrorIndex: HighsInt = -1;
pub const kHighsIllegalComplementarityViolation: f64 = kHighsInf;
pub const kHighsIllegalComplementarityCount: HighsInt = -1;
pub const kHighsIllegalDoubleMeasure: f64 = -kHighsInf;
pub const kHighsIllegalIntMeasure: HighsInt = -1;
pub const kMaxSemiVariableUpper: f64 = 1e5;
pub const kExcessivePrimalValue: f64 = 1e25;
pub const kHashIsDuplicate: HighsInt = -1;
pub const kModelValueToStringTolerance: f64 = 1e-15;
pub const kRangingValueToStringTolerance: f64 = 1e-13;
pub const kHighsSolutionValueToStringTolerance: f64 = 1e-13;
pub const kGlpsolSolutionValueToStringTolerance: f64 = 1e-12;
pub const kNoLink: HighsInt = -1;
pub const kPivotIllegal: i8 = -1;
pub const kPivotLogical: i8 = 0;
pub const kPivotUnit: i8 = 1;
pub const kPivotRowSingleton: i8 = 2;
pub const kPivotColSingleton: i8 = 3;
pub const kPivotMarkowitz: i8 = 4;
// HiGHS: HConst.h enum ConstraintType
pub const EQ: HighsInt = 0;
pub const LEQ: HighsInt = 1;
pub const GEQ: HighsInt = 2;
pub const BOUND: HighsInt = 3;
pub const FREE: HighsInt = 4;

// HiGHS: HConst.h enum PdlpFeaturesOff
pub const kPdlpAllFeaturesOn: HighsInt = 0;
pub const kPdlpScalingOff: HighsInt = 1;
pub const kPdlpRestartOff: HighsInt = 2;
pub const kPdlpAdaptiveStepSizeOff: HighsInt = 4;
pub const kPdlpAllFeaturesOff: HighsInt = 7;

// HiGHS: HConst.h enum PdlpScalingBit
pub const kPdlpScalingMin: HighsInt = 0;
pub const kPdlpScalingRuiz: HighsInt = 1;
pub const kPdlpScalingL2: HighsInt = 2;
pub const kPdlpScalingPC: HighsInt = 4;
pub const kPdlpScalingMax: HighsInt = 7;

// HiGHS: HConst.h enum PdlpStepSizeStrategy
pub const kPdlpStepSizeStrategyFixed: HighsInt = 0;
pub const kPdlpStepSizeStrategyMin: HighsInt = 0;
pub const kPdlpStepSizeStrategyAdaptive: HighsInt = 1;
pub const kPdlpStepSizeStrategyMalitskyPock: HighsInt = 2;
pub const kPdlpStepSizeStrategyPid: HighsInt = 3;
pub const kPdlpStepSizeStrategyMax: HighsInt = 3;

// HiGHS: HConst.h enum PdlpRestartStrategy
pub const kPdlpRestartStrategyOff: HighsInt = 0;
pub const kPdlpRestartStrategyMin: HighsInt = 0;
pub const kPdlpRestartStrategyFixed: HighsInt = 1;
pub const kPdlpRestartStrategyAdaptive: HighsInt = 2;
pub const kPdlpRestartStrategyHalpern: HighsInt = 3;
pub const kPdlpRestartStrategyMax: HighsInt = 3;
