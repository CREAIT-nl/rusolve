//! HiGHS source: simplex/SimplexTimer.h (clock ids only)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Out of scope: timing. Only the clock-id constants, so timer calls translate verbatim
//! (`analysis->simplexTimerStart(InvertClock)` -> `ekk.analysis_.simplex_timer_start(InvertClock, 0)`).

use crate::HighsInt;

// HiGHS: SimplexTimer.h enum iClockSimplex
pub const SimplexTotalClock: HighsInt = 0;
pub const SimplexIzDseWtClock: HighsInt = 1;
pub const SimplexDualPhase1Clock: HighsInt = 2;
pub const SimplexDualPhase2Clock: HighsInt = 3;
pub const SimplexPrimalPhase1Clock: HighsInt = 4;
pub const SimplexPrimalPhase2Clock: HighsInt = 5;
pub const Group1Clock: HighsInt = 6;
pub const IterateClock: HighsInt = 7;
pub const IterateDualRebuildClock: HighsInt = 8;
pub const IteratePrimalRebuildClock: HighsInt = 9;
pub const IterateChuzrClock: HighsInt = 10;
pub const IterateChuzcClock: HighsInt = 11;
pub const IterateFtranClock: HighsInt = 12;
pub const IterateVerifyClock: HighsInt = 13;
pub const IterateDualClock: HighsInt = 14;
pub const IteratePrimalClock: HighsInt = 15;
pub const IterateDevexIzClock: HighsInt = 16;
pub const IteratePivotsClock: HighsInt = 17;
pub const initialiseSimplexLpBasisAndFactorClock: HighsInt = 18;
pub const ScaleClock: HighsInt = 19;
pub const CrashClock: HighsInt = 20;
pub const BasisConditionClock: HighsInt = 21;
pub const matrixSetupClock: HighsInt = 22;
pub const setNonbasicMoveClock: HighsInt = 23;
pub const allocateSimplexArraysClock: HighsInt = 24;
pub const initialiseSimplexCostBoundsClock: HighsInt = 25;
pub const DseIzClock: HighsInt = 26;
pub const InvertClock: HighsInt = 27;
pub const PermWtClock: HighsInt = 28;
pub const ComputeDualClock: HighsInt = 29;
pub const CorrectDualClock: HighsInt = 30;
pub const CollectPrIfsClock: HighsInt = 31;
pub const ComputePrIfsClock: HighsInt = 32;
pub const ComputeDuIfsClock: HighsInt = 33;
pub const ComputePrimalClock: HighsInt = 34;
pub const ComputeDuObjClock: HighsInt = 35;
pub const ComputePrObjClock: HighsInt = 36;
pub const ReportRebuildClock: HighsInt = 37;
pub const ChuzrDualClock: HighsInt = 38;
pub const Chuzr1Clock: HighsInt = 39;
pub const Chuzr2Clock: HighsInt = 40;
pub const ChuzcPrimalClock: HighsInt = 41;
pub const ChuzcHyperInitialiselClock: HighsInt = 42;
pub const ChuzcHyperBasicFeasibilityChangeClock: HighsInt = 43;
pub const ChuzcHyperDualClock: HighsInt = 44;
pub const ChuzcHyperClock: HighsInt = 45;
pub const Chuzc0Clock: HighsInt = 46;
pub const PriceChuzc1Clock: HighsInt = 47;
pub const Chuzc1Clock: HighsInt = 48;
pub const Chuzc2Clock: HighsInt = 49;
pub const Chuzc3Clock: HighsInt = 50;
pub const Chuzc4Clock: HighsInt = 51;
pub const Chuzc4a0Clock: HighsInt = 52;
pub const Chuzc4a1Clock: HighsInt = 53;
pub const Chuzc4bClock: HighsInt = 54;
pub const Chuzc4cClock: HighsInt = 55;
pub const Chuzc4dClock: HighsInt = 56;
pub const Chuzc4eClock: HighsInt = 57;
pub const Chuzc5Clock: HighsInt = 58;
pub const DevexWtClock: HighsInt = 59;
pub const BtranClock: HighsInt = 60;
pub const BtranBasicFeasibilityChangeClock: HighsInt = 61;
pub const BtranFullClock: HighsInt = 62;
pub const PriceClock: HighsInt = 63;
pub const PriceBasicFeasibilityChangeClock: HighsInt = 64;
pub const PriceFullClock: HighsInt = 65;
pub const FtranClock: HighsInt = 66;
pub const FtranDseClock: HighsInt = 67;
pub const BtranPseClock: HighsInt = 68;
pub const FtranMixParClock: HighsInt = 69;
pub const FtranMixFinalClock: HighsInt = 70;
pub const FtranBfrtClock: HighsInt = 71;
pub const UpdateRowClock: HighsInt = 72;
pub const UpdateDualClock: HighsInt = 73;
pub const UpdateDualBasicFeasibilityChangeClock: HighsInt = 74;
pub const UpdatePrimalClock: HighsInt = 75;
pub const DevexIzClock: HighsInt = 76;
pub const DevexUpdateWeightClock: HighsInt = 77;
pub const DseUpdateWeightClock: HighsInt = 78;
pub const UpdatePivotsClock: HighsInt = 79;
pub const UpdateFactorClock: HighsInt = 80;
pub const UpdateMatrixClock: HighsInt = 81;
pub const UpdateRowEpClock: HighsInt = 82;
pub const SimplexNumClock: HighsInt = 83;
