//! HiGHS source: util/HFactorConst.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete translation (constants are data).

use crate::HighsInt;

// HiGHS: HFactorConst.h enum UPDATE_METHOD
pub const kUpdateMethodFt: HighsInt = 1;
pub const kUpdateMethodPf: HighsInt = 2;
pub const kUpdateMethodMpf: HighsInt = 3;
pub const kUpdateMethodApf: HighsInt = 4;

pub const kMinPivotThreshold: f64 = 8e-4;
pub const kDefaultPivotThreshold: f64 = 0.1;
pub const kPivotThresholdChangeFactor: f64 = 5.0;
pub const kMaxPivotThreshold: f64 = 0.5;
pub const kMinPivotTolerance: f64 = 0.0;
pub const kDefaultPivotTolerance: f64 = 1e-10;
pub const kMaxPivotTolerance: f64 = 1.0;
pub const kHyperFtranL: f64 = 0.15;
pub const kHyperFtranU: f64 = 0.10;
pub const kHyperBtranL: f64 = 0.10;
pub const kHyperBtranU: f64 = 0.15;
pub const kHyperCancel: f64 = 0.05;
pub const kHyperResult: f64 = 0.10;
pub const kMultiBuildSyntheticTickMu: f64 = 1.0;
pub const kNumericalTroubleTolerance: f64 = 1e-7;
pub const kMultiNumericalTroubleTolerance: f64 = 1e-7;
pub const kSyntheticTickReinversionMinUpdateCount: HighsInt = 50;
pub const kMultiSyntheticTickReinversionMinUpdateCount: HighsInt =
    kSyntheticTickReinversionMinUpdateCount;
pub const kMCExtraEntriesMultiplier: HighsInt = 2;
pub const kMRExtraEntriesMultiplier: HighsInt = 2;
pub const kLFactorExtraEntriesMultiplier: HighsInt = 3;
pub const kUFactorExtraVectors: HighsInt = 1000;
pub const kUFactorExtraEntriesMultiplier: HighsInt = 3;
pub const kPFFPivotEntries: HighsInt = 1000;
pub const kPFVectors: HighsInt = 2000;
pub const kPFEntriesMultiplier: HighsInt = 4;
pub const kNewLRRowsExtraNz: HighsInt = 100;
// HiGHS: HFactorConst.h enum ReportLuOption
pub const kReportLuJustL: HighsInt = 1;
pub const kReportLuJustU: HighsInt = 2;
pub const kReportLuBoth: HighsInt = 3;
