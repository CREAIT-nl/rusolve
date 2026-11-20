//! HiGHS source: pdlp/cupdlp/glbopts.h, pdlp/cupdlp/cupdlp_defs.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Types, constants, macros, enums and structs of cuPDLP-C as HiGHS compiles
//! them without CUDA (`CUPDLP_CPU` defined, `DLONG` / `SFLOAT` undefined,
//! `PDLP_DEBUG_LOG == 0`).
//!
//! # Macro expansions (defined once here)
//!
//! | C macro (glbopts.h / cupdlp_defs.h)                 | Rust                                           |
//! |-----------------------------------------------------|------------------------------------------------|
//! | `cupdlp_int` / `cupdlp_float` / `cupdlp_bool`       | [`cupdlp_int`] = `i32`, [`cupdlp_float`] = `f64`, [`cupdlp_bool`] = `bool` |
//! | `cupdlp_retcode`                                    | [`cupdlp_retcode`] = `cupdlp_int`              |
//! | `RETCODE_OK` / `RETCODE_FAILED`                     | consts `0` / `1`                               |
//! | `MAX(a,b)` / `MIN(a,b)` = `((a)>(b))?(a):(b)`       | [`MAX`] / [`MIN`] (C ternary: NaN / signed-zero semantics kept) |
//! | `ABS(x)` = `((x)<0)?-(x):(x)`                       | [`ABS`] (`ABS(-0.0)` is `-0.0`, unlike `f64::abs`) |
//! | `POWF` = `pow`, `SQRTF` = `sqrt`                    | `x.powf(y)` (never `powi`), `x.sqrt()`         |
//! | `SAFEDIV_POS(X,Y)`                                  | [`SAFEDIV_POS`]                                |
//! | `EPS_TOL`, `EPS`, `CONVERGED_INTERVAL`, `INDETERMINATE_TOL`, `OUR_DBL_MAX` | consts             |
//! | `cupdlp_copy(dst, src, type, size)` / `CUPDLP_COPY_VEC` / `cupdlp_copy_vec` (memcpy) | [`cupdlp_copy`]: `dst[..size].copy_from_slice(&src[..size])` |
//! | `cupdlp_zero(var, type, size)` / `CUPDLP_ZERO_VEC` (memset 0) | [`cupdlp_zero`]                     |
//! | `cupdlp_init_*`, `CUPDLP_INIT_*`, `CUPDLP_INIT_ZERO_*`, `cupdlp_init_zero_vec_*` (malloc / calloc + NULL check) | `vec![0.0; n]` / `vec![0; n]` / `T::default()`; allocation cannot fail, so the `goto exit_cleanup` on NULL disappears |
//! | `CUPDLP_FREE`, `cupdlp_free`, `*_clear`              | `Drop` (nothing to write)                      |
//! | `CUPDLP_CALL(f)` (`if ((retcode = (f)) != RETCODE_OK) { retcode = RETCODE_FAILED; goto exit_cleanup; }`) | `retcode = f; if retcode != RETCODE_OK { retcode = RETCODE_FAILED; /* jump to the exit_cleanup code */ }` written out at each site, keeping the code after the label |
//! | `CUPDLP_CHECK_TIMEOUT(pdhg)` (solver.h)             | `PDHG_Compute_SolvingTime(pdhg); if pdhg.timers.dSolvingTime > pdhg.settings.dTimeLim { retcode = RETCODE_FAILED; /* exit_cleanup */ }` written out |
//! | `CUPDLP_ASSERT` = `assert`                           | `debug_assert!`                                |
//! | `cupdlp_printf`, `cupdlp_snprintf`, `CUPDLP_*_PRINT*`, `DEBUG_FUNC`, MATLAB / PYTHON variants | dropped: print only (logging rule) |
//! | `CHECK_CUBLAS`, `CHECK_CUSPARSE`, `CHECK_CUDA`, `cupdlp_*_cuda`, `#ifndef CUPDLP_CPU` blocks | dropped (CUDA only) |
//!
//! # C enums
//!
//! cuPDLP's `typedef enum {...} name;` are plain C enums whose values are
//! assigned from and to `cupdlp_int` (`stepsize->eLineSearchMethod =
//! intParam[E_LINE_SEARCH_METHOD]`, `*model_status = (cupdlp_int)termCode`).
//! Per the crate rule for plain C enums each becomes a type alias of
//! [`cupdlp_int`] plus one `pub const` per enumerator.
//!
//! # Structs
//!
//! Every `T *` member that cuPDLP mallocs once becomes an owned value
//! (`Vec<cupdlp_float>` for arrays, the struct itself for struct pointers).
//! `CUPDLPdata`'s three matrix pointers, which are NULL unless the matrix
//! format asks for them, are `Option<Box<..>>`. GPU-only members are
//! dropped. Fields are zero by `Default` (C `malloc` leaves them
//! uninitialised; every one is written before it is read).

use crate::HighsInt;

/// HiGHS: glbopts.h `#define cupdlp_int int` (no `DLONG`).
pub type cupdlp_int = i32;
/// HiGHS: glbopts.h `#define cupdlp_float double` (no `SFLOAT`).
pub type cupdlp_float = f64;
/// HiGHS: glbopts.h `#define cupdlp_bool bool`.
pub type cupdlp_bool = bool;
/// HiGHS: cupdlp_defs.h `typedef cupdlp_int cupdlp_retcode;`
pub type cupdlp_retcode = cupdlp_int;

/// HiGHS: glbopts.h `RETCODE_OK`.
pub const RETCODE_OK: cupdlp_retcode = 0;
/// HiGHS: glbopts.h `RETCODE_FAILED`.
pub const RETCODE_FAILED: cupdlp_retcode = 1;
/// HiGHS: glbopts.h `BOOL`.
pub const BOOL: cupdlp_int = 1;

/// HiGHS: glbopts.h `EPS_TOL`.
pub const EPS_TOL: cupdlp_float = 1E-18;
/// HiGHS: glbopts.h `EPS`.
pub const EPS: cupdlp_float = 1E-8;
/// HiGHS: glbopts.h `CONVERGED_INTERVAL`.
pub const CONVERGED_INTERVAL: cupdlp_int = 1;
/// HiGHS: glbopts.h `INDETERMINATE_TOL`.
pub const INDETERMINATE_TOL: cupdlp_float = 1e-9;
/// HiGHS: glbopts.h `OUR_DBL_MAX`.
pub const OUR_DBL_MAX: cupdlp_float = 1E+20;

/// HiGHS: cupdlp_defs.h `PDHG_USE_TIMERS`.
pub const PDHG_USE_TIMERS: cupdlp_int = 1;
/// HiGHS: cupdlp_defs.h `USE_MY_BLAS`.
pub const USE_MY_BLAS: cupdlp_int = 1;
/// HiGHS: cupdlp_defs.h `USE_KERNELS`.
pub const USE_KERNELS: cupdlp_int = 1;
/// HiGHS: cupdlp_defs.h `PDHG_STEPSIZE_REDUCTION_EXP`.
pub const PDHG_STEPSIZE_REDUCTION_EXP: cupdlp_float = 0.3;
/// HiGHS: cupdlp_defs.h `PDHG_STEPSIZE_GROWTH_EXP`.
pub const PDHG_STEPSIZE_GROWTH_EXP: cupdlp_float = 0.6;
/// HiGHS: cupdlp_defs.h `PDHG_DISPLAY_TERMINATION_CHECK`.
pub const PDHG_DISPLAY_TERMINATION_CHECK: cupdlp_int = 1;
/// HiGHS: cupdlp_defs.h `PDHG_PROJECT_INITIAL_X`.
pub const PDHG_PROJECT_INITIAL_X: cupdlp_int = 0;
/// HiGHS: cupdlp_defs.h `CUPDLP_DEBUG_INTERVAL`.
pub const CUPDLP_DEBUG_INTERVAL: cupdlp_int = 40;
/// HiGHS: cupdlp_defs.h `CUPDLP_RELEASE_INTERVAL`.
pub const CUPDLP_RELEASE_INTERVAL: cupdlp_int = 40;
/// HiGHS: cupdlp_defs.h `CUPDLP_DUMP_ITERATES_STATS`.
pub const CUPDLP_DUMP_ITERATES_STATS: cupdlp_int = 1;
/// HiGHS: cupdlp_defs.h `CUPDLP_DUMP_LINESEARCH_STATS`.
pub const CUPDLP_DUMP_LINESEARCH_STATS: cupdlp_int = 1;
/// HiGHS: cupdlp_defs.h `CUPDLP_INEXACT_EPS`.
pub const CUPDLP_INEXACT_EPS: cupdlp_float = 1e-4;
/// HiGHS: cupdlp_defs.h `PDLP_DEBUG_LOG` (0: every `#if PDLP_DEBUG_LOG` block is dead).
pub const PDLP_DEBUG_LOG: cupdlp_int = 0;

/// HiGHS: glbopts.h `MAX(a, b)` = `((a) > (b)) ? (a) : (b)`.
#[inline]
pub fn MAX<T: PartialOrd>(a: T, b: T) -> T {
    if a > b { a } else { b }
}

/// HiGHS: glbopts.h `MIN(a, b)` = `((a) < (b)) ? (a) : (b)`.
#[inline]
pub fn MIN<T: PartialOrd>(a: T, b: T) -> T {
    if a < b { a } else { b }
}

/// HiGHS: glbopts.h `ABS(x)` = `((x) < 0) ? -(x) : (x)`.
#[inline]
pub fn ABS(x: cupdlp_float) -> cupdlp_float {
    if x < 0.0 { -x } else { x }
}

/// HiGHS: glbopts.h `SAFEDIV_POS(X, Y)` = `(Y) < EPS_TOL ? ((X) / EPS_TOL) : (X) / (Y)`.
#[inline]
pub fn SAFEDIV_POS(x: cupdlp_float, y: cupdlp_float) -> cupdlp_float {
    if y < EPS_TOL { x / EPS_TOL } else { x / y }
}

/// HiGHS: glbopts.h `cupdlp_copy(dst, src, type, size)` / `CUPDLP_COPY_VEC` /
/// CupdlpWrapper.h `cupdlp_copy_vec` (`memcpy` of `size` elements).
#[inline]
pub fn cupdlp_copy<T: Copy>(dst: &mut [T], src: &[T], size: cupdlp_int) {
    let n = size as usize;
    dst[..n].copy_from_slice(&src[..n]);
}

/// HiGHS: glbopts.h `cupdlp_zero(var, type, size)` / `CUPDLP_ZERO_VEC` (`memset` 0).
#[inline]
pub fn cupdlp_zero<T: Copy + Default>(var: &mut [T], size: cupdlp_int) {
    for v in var[..size as usize].iter_mut() {
        *v = T::default();
    }
}

/// HiGHS: cupdlp_defs.h `termination_code`.
pub type termination_code = cupdlp_int;
pub const OPTIMAL: termination_code = 0;
pub const INFEASIBLE: termination_code = 1;
pub const UNBOUNDED: termination_code = 2;
pub const INFEASIBLE_OR_UNBOUNDED: termination_code = 3;
pub const TIMELIMIT_OR_ITERLIMIT: termination_code = 4;
pub const FEASIBLE: termination_code = 5;

/// HiGHS: cupdlp_defs.h `termination_iterate`.
pub type termination_iterate = cupdlp_int;
pub const LAST_ITERATE: termination_iterate = 0;
pub const AVERAGE_ITERATE: termination_iterate = 1;

/// HiGHS: cupdlp_defs.h `pdhg_linesearch`.
pub type pdhg_linesearch = cupdlp_int;
pub const PDHG_FIXED_LINESEARCH: pdhg_linesearch = 0;
pub const PDHG_MALITSKY_POCK_LINESEARCH: pdhg_linesearch = 1;
pub const PDHG_ADAPTIVE_LINESEARCH: pdhg_linesearch = 2;

/// HiGHS: cupdlp_defs.h `pdhg_restart`.
pub type pdhg_restart = cupdlp_int;
pub const PDHG_WITHOUT_RESTART: pdhg_restart = 0;
pub const PDHG_GPU_RESTART: pdhg_restart = 1;
pub const PDHG_CPU_RESTART: pdhg_restart = 2;

/// HiGHS: cupdlp_defs.h `CUPDLP_DEVICE`.
pub type CUPDLP_DEVICE = cupdlp_int;
pub const CPU: CUPDLP_DEVICE = 0;
pub const SINGLE_GPU: CUPDLP_DEVICE = 1;
pub const MULTI_GPU: CUPDLP_DEVICE = 2;

/// HiGHS: cupdlp_defs.h `CUPDLP_MATRIX_FORMAT`.
pub type CUPDLP_MATRIX_FORMAT = cupdlp_int;
pub const DENSE: CUPDLP_MATRIX_FORMAT = 0;
pub const CSR: CUPDLP_MATRIX_FORMAT = 1;
pub const CSC: CUPDLP_MATRIX_FORMAT = 2;
pub const CSR_CSC: CUPDLP_MATRIX_FORMAT = 3;

/// HiGHS: cupdlp_defs.h `CUPDLP_INT_USER_PARAM_INDEX` (indices into `intParam`).
pub type CUPDLP_INT_USER_PARAM_INDEX = cupdlp_int;
pub const N_ITER_LIM: usize = 0;
pub const IF_SCALING: usize = 1;
pub const I_SCALING_METHOD: usize = 2;
pub const E_LINE_SEARCH_METHOD: usize = 3;
pub const E_RESTART_METHOD: usize = 4;
pub const IF_RUIZ_SCALING: usize = 5;
pub const IF_L2_SCALING: usize = 6;
pub const IF_PC_SCALING: usize = 7;
pub const N_LOG_LEVEL: usize = 8;
pub const N_LOG_INTERVAL: usize = 9;
pub const IF_PRESOLVE: usize = 10;
pub const I_INF_NORM_ABS_LOCAL_TERMINATION: usize = 11;
pub const N_INT_USER_PARAM: usize = 12;

/// HiGHS: cupdlp_defs.h `CUPDLP_FLOAT_USER_PARAM_INDEX` (indices into `floatParam`).
pub type CUPDLP_FLOAT_USER_PARAM_INDEX = cupdlp_int;
pub const D_SCALING_LIMIT: usize = 0;
pub const D_PRIMAL_TOL: usize = 1;
pub const D_DUAL_TOL: usize = 2;
pub const D_GAP_TOL: usize = 3;
pub const D_FEAS_TOL: usize = 4;
pub const D_TIME_LIM: usize = 5;
pub const N_FLOAT_USER_PARAM: usize = 6;

/// HiGHS: cupdlp_restart.h `PDHG_restart_choice`.
pub type PDHG_restart_choice = cupdlp_int;
pub const PDHG_NO_RESTART: PDHG_restart_choice = 0;
pub const PDHG_RESTART_TO_CURRENT: PDHG_restart_choice = 1;
pub const PDHG_RESTART_TO_AVERAGE: PDHG_restart_choice = 2;

/// HiGHS: cupdlp_defs.h `struct CUPDLP_CUDA_DENSE_VEC` (`CUPDLPvec`; `cuda_vec` dropped).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPvec {
    pub len: cupdlp_int,
    pub data: Vec<cupdlp_float>,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_DENSE_MATRIX` (`CUPDLPdense`).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPdense {
    pub nRows: cupdlp_int,
    pub nCols: cupdlp_int,
    pub data: Vec<cupdlp_float>,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_CSR_MATRIX` (`CUPDLPcsr`; `cuda_csr` dropped).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPcsr {
    pub nRows: cupdlp_int,
    pub nCols: cupdlp_int,
    pub nMatElem: cupdlp_int,
    pub rowMatBeg: Vec<cupdlp_int>,
    pub rowMatIdx: Vec<cupdlp_int>,
    pub rowMatElem: Vec<cupdlp_float>,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_CSC_MATRIX` (`CUPDLPcsc`; `cuda_csc` dropped).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPcsc {
    pub nRows: cupdlp_int,
    pub nCols: cupdlp_int,
    pub nMatElem: cupdlp_int,
    pub colMatBeg: Vec<cupdlp_int>,
    pub colMatIdx: Vec<cupdlp_int>,
    pub colMatElem: Vec<cupdlp_float>,
    pub MatElemNormInf: cupdlp_float,
}

/// The `void *src` + `CUPDLP_MATRIX_FORMAT src_matrix_format` pair of
/// `*_alloc_matrix` / `data_alloc` / `problem_alloc`: the format tag picks the
/// cast, so the Rust parameter is the tagged reference.
#[derive(Clone, Copy, Debug)]
pub enum CupdlpSrcMatrix<'a> {
    Dense(&'a CUPDLPdense),
    Csr(&'a CUPDLPcsr),
    Csc(&'a CUPDLPcsc),
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_DATA` (`CUPDLPdata`).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPdata {
    pub nRows: cupdlp_int,
    pub nCols: cupdlp_int,
    pub matrix_format: CUPDLP_MATRIX_FORMAT,
    pub dense_matrix: Option<Box<CUPDLPdense>>,
    pub csr_matrix: Option<Box<CUPDLPcsr>>,
    pub csc_matrix: Option<Box<CUPDLPcsc>>,
    pub device: CUPDLP_DEVICE,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_SETTINGS` (`CUPDLPsettings`).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPsettings {
    pub ifScaling: cupdlp_int,
    pub iScalingMethod: cupdlp_int,
    pub dScalingLimit: cupdlp_float,
    pub dPrimalTol: cupdlp_float,
    pub dDualTol: cupdlp_float,
    pub dGapTol: cupdlp_float,
    pub nIterLim: cupdlp_int,
    pub dTimeLim: cupdlp_float,
    pub nLogLevel: cupdlp_int,
    pub nLogInterval: cupdlp_int,
    pub eRestartMethod: pdhg_restart,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_PROBLEM` (`CUPDLPproblem`).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPproblem {
    pub data: CUPDLPdata,
    pub lower: Vec<cupdlp_float>,
    pub upper: Vec<cupdlp_float>,
    pub cost: Vec<cupdlp_float>,
    pub rhs: Vec<cupdlp_float>,
    pub dMaxCost: cupdlp_float,
    pub dMaxRhs: cupdlp_float,
    pub dMaxRowBound: cupdlp_float,
    pub nRows: cupdlp_int,
    pub nCols: cupdlp_int,
    pub nEqs: cupdlp_int,
    pub hasLower: Vec<cupdlp_float>,
    pub hasUpper: Vec<cupdlp_float>,
    pub offset: cupdlp_float,
    pub sense_origin: cupdlp_float,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_RES_OBJ` (`CUPDLPresobj`).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPresobj {
    pub dFeasTol: cupdlp_float,
    pub dPrimalObj: cupdlp_float,
    pub dDualObj: cupdlp_float,
    pub dDualityGap: cupdlp_float,
    pub dComplementarity: cupdlp_float,
    pub dPrimalFeas: cupdlp_float,
    pub dDualFeas: cupdlp_float,
    pub dRelObjGap: cupdlp_float,
    pub primalResidual: Vec<cupdlp_float>,
    pub dualResidual: Vec<cupdlp_float>,
    pub dSlackPos: Vec<cupdlp_float>,
    pub dSlackNeg: Vec<cupdlp_float>,
    pub dSlackPosAverage: Vec<cupdlp_float>,
    pub dSlackNegAverage: Vec<cupdlp_float>,
    pub dLowerFiltered: Vec<cupdlp_float>,
    pub dUpperFiltered: Vec<cupdlp_float>,

    /* for infeasibility detection */
    pub primalCode: termination_code,
    pub dualCode: termination_code,
    pub termInfeasIterate: termination_iterate,
    pub dPrimalInfeasObj: cupdlp_float,
    pub dDualInfeasObj: cupdlp_float,
    pub dPrimalInfeasRes: cupdlp_float,
    pub dDualInfeasRes: cupdlp_float,
    pub dPrimalInfeasObjAverage: cupdlp_float,
    pub dDualInfeasObjAverage: cupdlp_float,
    pub dPrimalInfeasResAverage: cupdlp_float,
    pub dDualInfeasResAverage: cupdlp_float,
    pub primalInfeasRay: Vec<cupdlp_float>,
    pub primalInfeasConstr: Vec<cupdlp_float>,
    pub primalInfeasBound: Vec<cupdlp_float>,
    pub dualInfeasRay: Vec<cupdlp_float>,
    pub dualInfeasLbRay: Vec<cupdlp_float>,
    pub dualInfeasUbRay: Vec<cupdlp_float>,
    pub dualInfeasConstr: Vec<cupdlp_float>,

    pub dPrimalObjAverage: cupdlp_float,
    pub dDualObjAverage: cupdlp_float,
    pub dDualityGapAverage: cupdlp_float,
    pub dComplementarityAverage: cupdlp_float,
    pub dPrimalFeasAverage: cupdlp_float,
    pub dDualFeasAverage: cupdlp_float,
    pub dRelObjGapAverage: cupdlp_float,
    pub primalResidualAverage: Vec<cupdlp_float>,
    pub dualResidualAverage: Vec<cupdlp_float>,

    pub dPrimalFeasLastRestart: cupdlp_float,
    pub dDualFeasLastRestart: cupdlp_float,
    pub dDualityGapLastRestart: cupdlp_float,

    pub dPrimalFeasLastCandidate: cupdlp_float,
    pub dDualFeasLastCandidate: cupdlp_float,
    pub dDualityGapLastCandidate: cupdlp_float,

    pub termCode: termination_code,
    pub termIterate: termination_iterate,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_ITERATES` (`CUPDLPiterates`).
///
/// `x[k % 2]` holds `x^k`, `x[(k + 1) % 2]` receives `x^{k+1}` (C comment);
/// the four `CUPDLPvec *name[2]` arrays become `[CUPDLPvec; 2]`.
#[derive(Clone, Debug, Default)]
pub struct CUPDLPiterates {
    pub nRows: cupdlp_int,
    pub nCols: cupdlp_int,

    pub iLastRestartIter: cupdlp_int,
    pub dLastRestartDualityGap: cupdlp_float,
    pub dLastRestartBeta: cupdlp_float,
    pub xSum: Vec<cupdlp_float>,
    pub ySum: Vec<cupdlp_float>,

    pub xLastRestart: Vec<cupdlp_float>,
    pub yLastRestart: Vec<cupdlp_float>,

    pub x: [CUPDLPvec; 2],
    pub y: [CUPDLPvec; 2],
    pub ax: [CUPDLPvec; 2],
    pub aty: [CUPDLPvec; 2],
    pub xAverage: CUPDLPvec,
    pub yAverage: CUPDLPvec,
    pub axAverage: CUPDLPvec,
    pub atyAverage: CUPDLPvec,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_STEPSIZE` (`CUPDLPstepsize`).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPstepsize {
    pub eLineSearchMethod: pdhg_linesearch,
    pub dPrimalStep: cupdlp_float,
    pub dDualStep: cupdlp_float,
    pub dSumPrimalStep: cupdlp_float,
    pub dSumDualStep: cupdlp_float,
    pub dBeta: cupdlp_float,
    pub dTheta: cupdlp_float,
    pub nStepSizeIter: cupdlp_int,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_SCALING` (`CUPDLPscaling`).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPscaling {
    pub ifScaled: cupdlp_int,
    pub rowScale: Vec<cupdlp_float>,
    pub colScale: Vec<cupdlp_float>,
    pub ifRuizScaling: cupdlp_int,
    pub ifL2Scaling: cupdlp_int,
    pub ifPcScaling: cupdlp_int,
    pub RuizTimes: cupdlp_int,
    pub RuizNorm: cupdlp_float,
    pub PcAlpha: cupdlp_float,
    pub dNormCost: cupdlp_float,
    pub dNormRhs: cupdlp_float,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_TIMERS` (`CUPDLPtimers`; `PDHG_USE_TIMERS == 1`
/// members kept, `#ifndef CUPDLP_CPU` members dropped).
#[derive(Clone, Debug, Default)]
pub struct CUPDLPtimers {
    pub nIter: cupdlp_int,
    pub dSolvingTime: cupdlp_float,
    pub dSolvingBeg: cupdlp_float,
    pub dScalingTime: cupdlp_float,
    pub dPresolveTime: cupdlp_float,
    pub dAtyTime: cupdlp_float,
    pub dAxTime: cupdlp_float,
    pub dComputeResidualsTime: cupdlp_float,
    pub dUpdateIterateTime: cupdlp_float,
    pub nAtyCalls: cupdlp_int,
    pub nAxCalls: cupdlp_int,
    pub nComputeResidualsCalls: cupdlp_int,
    pub nUpdateIterateCalls: cupdlp_int,
}

/// HiGHS: cupdlp_defs.h `struct CUPDLP_WORK` (`CUPDLPwork`).
///
/// Owns what its pointers point at: `CupdlpWrapper.cpp::solveLpCupdlp`
/// moves the problem and scaling in (`w->problem = prob; w->scaling =
/// scaling;`), `PDHG_Alloc` creates the rest. The CUDA handles and the
/// `PDLP_DEBUG_LOG` members are dropped.
#[derive(Clone, Debug, Default)]
pub struct CUPDLPwork {
    pub problem: CUPDLPproblem,
    pub settings: CUPDLPsettings,
    pub resobj: CUPDLPresobj,
    pub iterates: CUPDLPiterates,
    pub stepsize: CUPDLPstepsize,
    pub scaling: CUPDLPscaling,
    pub timers: CUPDLPtimers,
    pub buffer: CUPDLPvec,
    pub buffer2: Vec<cupdlp_float>,
    pub buffer3: Vec<cupdlp_float>,

    pub rowScale: Vec<cupdlp_float>,
    pub colScale: Vec<cupdlp_float>,
}

#[cfg(test)]
mod tests {
    //! Unit tests: the macro equivalents.
    use super::*;

    #[test]
    fn macros_keep_c_ternary_semantics() {
        assert_eq!(MAX(1.0, 2.0), 2.0);
        assert_eq!(MIN(1, 2), 1);
        // ((x) < 0) ? -(x) : (x): -0.0 is not < 0, so it is returned as is.
        assert_eq!(ABS(-0.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(ABS(-3.0), 3.0);
        // MAX(NaN, 1) = (NaN > 1) ? NaN : 1 = 1; MAX(1, NaN) = NaN.
        assert_eq!(MAX(f64::NAN, 1.0), 1.0);
        assert!(MAX(1.0, f64::NAN).is_nan());
        assert_eq!(SAFEDIV_POS(1.0, 0.0), 1.0 / EPS_TOL);
        assert_eq!(SAFEDIV_POS(1.0, 4.0), 0.25);
    }

    #[test]
    fn copy_and_zero_touch_size_elements() {
        let mut d = vec![9.0; 4];
        cupdlp_copy(&mut d, &[1.0, 2.0, 3.0], 2);
        assert_eq!(d, vec![1.0, 2.0, 9.0, 9.0]);
        cupdlp_zero(&mut d, 3);
        assert_eq!(d, vec![0.0, 0.0, 0.0, 9.0]);
    }

    #[test]
    fn enum_values_match_c() {
        assert_eq!(N_INT_USER_PARAM, 12);
        assert_eq!(N_FLOAT_USER_PARAM, 6);
        assert_eq!(TIMELIMIT_OR_ITERLIMIT, 4);
        assert_eq!(PDHG_ADAPTIVE_LINESEARCH, 2);
        assert_eq!(CSR_CSC, 3);
        let _: HighsInt = OPTIMAL;
    }
}
