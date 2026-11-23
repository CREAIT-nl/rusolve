//! HiGHS source: ipm/hipo/factorhighs/FactorHiGHSSettings.h, ipm/hipo/factorhighs/ReturnValues.h, ipm/hipo/factorhighs/Timing.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! FactorHiGHS, the supernodal multifrontal LDL^T solver of HiPO: settings,
//! return values and timing items. `HIPO_TIMING_LEVEL` is 0 and
//! `HIPO_COLLECT_EXPENSIVE_DATA` is undefined in the HiGHS build, so the
//! `HIPO_CLOCK_*` macros expand to nothing and are not ported.

pub mod analyse;
pub mod call_and_time_blas;
pub mod clique_stack;
pub mod data_collector;
pub mod dense_fact_hybrid;
pub mod dense_fact_kernel;
pub mod dgemm_parallel;
pub mod factor_highs;
pub mod factorise;
pub mod format_handler;
pub mod hybrid_hybrid_format_handler;
pub mod hybrid_solve_handler;
pub mod numeric;
pub mod solve_handler;
pub mod swaps;
pub mod symbolic;

use super::Int;

// HiGHS: FactorHiGHSSettings.h
pub const kStartThreshRelax: Int = 256;
pub const kUpperRatioRelax: f64 = 0.02;
pub const kLowerRatioRelax: f64 = 0.01;
pub const kMaxIterRelax: Int = 20;
pub const kSnSizeRelax: Int = 16;
pub const kSpopsWeightSn: f64 = 50.0;
pub const kAlphaBK: f64 = 0.01;
pub const kBlockGrainSize: Int = 1;
pub const kBlockParallelThreshold: Int = 5;
pub const kMinConsecutiveSums: Int = 1;
pub const kDynamicDiagCoeff: f64 = 1e-24;
pub const kMetisSeed: Int = 42;

/// HiGHS: FactorHiGHSSettings.h `struct Regul`
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Regul {
    pub primal: f64,
    pub dual: f64,
}

// HiGHS: ReturnValues.h `enum RetValue`
pub const kRetOk: Int = 0;
pub const kRetInvalidInput: Int = 1;
pub const kRetOutOfMemory: Int = 2;
pub const kRetInvalidPivot: Int = 3;
pub const kRetOrderingError: Int = 4;
pub const kRetIntOverflow: Int = 5;
pub const kRetInvalidPointer: Int = 6;
pub const kRetGeneric: Int = 7;

// HiGHS: Timing.h `enum TimeItems`
pub const kTimeAnalyse: Int = 0;
pub const kTimeAnalyseTree: Int = 1;
pub const kTimeAnalyseCount: Int = 2;
pub const kTimeAnalysePattern: Int = 3;
pub const kTimeAnalyseSn: Int = 4;
pub const kTimeAnalyseReorder: Int = 5;
pub const kTimeAnalyseRelInd: Int = 6;
pub const kTimeFactorise: Int = 7;
pub const kTimeFactorisePrepare: Int = 8;
pub const kTimeFactoriseAssembleOriginal: Int = 9;
pub const kTimeFactoriseAssembleChildrenFrontal: Int = 10;
pub const kTimeFactoriseAssembleChildrenClique: Int = 11;
pub const kTimeFactoriseDenseFact: Int = 12;
pub const kTimeDenseFact_main: Int = 13;
pub const kTimeDenseFact_schur: Int = 14;
pub const kTimeDenseFact_kernel: Int = 15;
pub const kTimeDenseFact_convert: Int = 16;
pub const kTimeDenseFact_pivoting: Int = 17;
pub const kTimeFactoriseTerminate: Int = 18;
pub const kTimeSolve: Int = 19;
pub const kTimeSolvePrepare: Int = 20;
pub const kTimeSolveSolve: Int = 21;
pub const kTimeSolveSolve_dense: Int = 22;
pub const kTimeSolveSolve_sparse: Int = 23;
pub const kTimeSolveSolve_swap: Int = 24;
pub const kTimeBlasStart: Int = 25;
pub const kTimeBlas_copy: Int = 25;
pub const kTimeBlas_axpy: Int = 26;
pub const kTimeBlas_scal: Int = 27;
pub const kTimeBlas_swap: Int = 28;
pub const kTimeBlas_gemv: Int = 29;
pub const kTimeBlas_trsv: Int = 30;
pub const kTimeBlas_tpsv: Int = 31;
pub const kTimeBlas_ger: Int = 32;
pub const kTimeBlas_trsm: Int = 33;
pub const kTimeBlas_syrk: Int = 34;
pub const kTimeBlas_gemm: Int = 35;
pub const kTimeBlasEnd: Int = 35;
pub const kTimeSize: Int = 36;

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn time_items_layout() {
        assert_eq!(kTimeBlasEnd - kTimeBlasStart + 1, 11);
        assert_eq!(kTimeSize, kTimeBlasEnd + 1);
        assert_eq!(
            Regul::default(),
            Regul {
                primal: 0.0,
                dual: 0.0
            }
        );
    }
}
