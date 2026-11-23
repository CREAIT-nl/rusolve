//! HiGHS source: pdlp/CupdlpWrapper.cpp (lines 280-717, 839-848)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiGHS LP -> cuPDLP standard form, the cuPDLP problem constructors and
//! the option translation.
//!
//! `formulateLP_highs` returns its eleven `double**` / `int**` arrays
//! through out-parameters that it `malloc`s; the port takes `&mut Vec<_>`
//! for each and assigns a fresh vector (the caller frees them; here they
//! drop). Scalar out-parameters are `&mut`.

use super::super::cupdlp::utils::*;
use super::super::cupdlp::*;
use crate::prelude::*;

/// HiGHS: CupdlpWrapper.cpp::formulateLP_highs (lines 280-448)
///
/// Row classification (`hconst.rs` `ConstraintType`, threshold `1e20`):
/// `EQ` (finite equal bounds), `GEQ` (lower only), `LEQ` (upper only,
/// negated), `BOUND` (both finite, or neither: a slack column `x_j` with
/// the row bounds is appended and the row becomes `a'x - x_j = 0`, `nnz`
/// and `nEqs` grow by one). The new row order is `EQ` / `BOUND` rows first,
/// then `LEQ` / `GEQ`, recorded in `constraint_new_idx`. Bounds beyond
/// `1e20` become `+-INFINITY`. `local_log_level` is unused upstream; the
/// free-row `printf` warning (line 338) is dropped. `cupdlp_init_*` allocates without initialising; every entry is
/// written before it is read, so `vec![0; n]` is equivalent.
#[allow(clippy::too_many_arguments)]
pub fn formulate_lp_highs(
    local_log_level: cupdlp_int,
    lp: &HighsLp,
    cost: &mut Vec<f64>,
    nCols: &mut i32,
    nRows: &mut i32,
    nnz: &mut i32,
    nEqs: &mut i32,
    csc_beg: &mut Vec<i32>,
    csc_idx: &mut Vec<i32>,
    csc_val: &mut Vec<f64>,
    rhs: &mut Vec<f64>,
    lower: &mut Vec<f64>,
    upper: &mut Vec<f64>,
    offset: &mut f64,
    sense_origin: &mut f64,
    nCols_origin: &mut i32,
    constraint_new_idx: &mut Vec<i32>,
    constraint_type: &mut Vec<i32>,
) -> i32 {
    let retcode = 0;

    // problem size for malloc
    let nCols_clp = lp.num_col_;
    let nRows_clp = lp.num_row_;
    let nnz_clp = lp.a_matrix_.start_[lp.num_col_ as usize];
    *nCols_origin = nCols_clp;
    *nRows = nRows_clp; // need not recalculate
    *nCols = nCols_clp; // need recalculate
    *nEqs = 0; // need recalculate
    *nnz = nnz_clp; // need recalculate
    *offset = lp.offset_; // need not recalculate
    if lp.sense_ == ObjSense::kMinimize {
        *sense_origin = 1.0;
    } else if lp.sense_ == ObjSense::kMaximize {
        *sense_origin = -1.0;
    }

    let lhs_clp = &lp.row_lower_;
    let rhs_clp = &lp.row_upper_;
    let A_csc_beg = &lp.a_matrix_.start_;
    let A_csc_idx = &lp.a_matrix_.index_;
    let A_csc_val = &lp.a_matrix_.value_;
    let mut has_lower: bool;
    let mut has_upper: bool;

    // cupdlp_init_int(*constraint_type, *nRows);
    *constraint_type = vec![0; *nRows as usize];
    // cupdlp_init_int(*constraint_new_idx, *nRows);
    *constraint_new_idx = vec![0; *nRows as usize];

    // recalculate nRows and nnz for Ax - z = 0
    for i in 0..nRows_clp {
        has_lower = lhs_clp[i as usize] > -1e20;
        has_upper = rhs_clp[i as usize] < 1e20;

        // count number of equations and rows
        if has_lower && has_upper && lhs_clp[i as usize] == rhs_clp[i as usize] {
            constraint_type[i as usize] = EQ;
            *nEqs += 1;
        } else if has_lower && !has_upper {
            constraint_type[i as usize] = GEQ;
        } else if !has_lower && has_upper {
            constraint_type[i as usize] = LEQ;
        } else if has_lower && has_upper {
            constraint_type[i as usize] = BOUND;
            *nCols += 1;
            *nnz += 1;
            *nEqs += 1;
        } else {
            // printf("Warning: constraint %d has no lower and upper bound\n", i); // line 338
            constraint_type[i as usize] = BOUND;
            *nCols += 1;
            *nnz += 1;
            *nEqs += 1;
        }
    }

    // allocate memory
    // cupdlp_init_double(*cost, *nCols);
    *cost = vec![0.0; *nCols as usize];
    // cupdlp_init_double(*lower, *nCols);
    *lower = vec![0.0; *nCols as usize];
    // cupdlp_init_double(*upper, *nCols);
    *upper = vec![0.0; *nCols as usize];
    // cupdlp_init_int(*csc_beg, *nCols + 1);
    *csc_beg = vec![0; (*nCols + 1) as usize];
    // cupdlp_init_int(*csc_idx, *nnz);
    *csc_idx = vec![0; *nnz as usize];
    // cupdlp_init_double(*csc_val, *nnz);
    *csc_val = vec![0.0; *nnz as usize];
    // cupdlp_init_double(*rhs, *nRows);
    *rhs = vec![0.0; *nRows as usize];

    // cost, lower, upper
    for i in 0..nCols_clp {
        cost[i as usize] = lp.col_cost_[i as usize] * (*sense_origin);
        lower[i as usize] = lp.col_lower_[i as usize];

        upper[i as usize] = lp.col_upper_[i as usize];
    }
    // slack costs
    for i in nCols_clp..*nCols {
        cost[i as usize] = 0.0;
    }
    // slack bounds
    {
        let mut j = nCols_clp;
        for i in 0..*nRows {
            if constraint_type[i as usize] == BOUND {
                lower[j as usize] = lhs_clp[i as usize];
                upper[j as usize] = rhs_clp[i as usize];
                j += 1;
            }
        }
    }

    for i in 0..*nCols {
        if lower[i as usize] < -1e20 {
            lower[i as usize] = -f64::INFINITY;
        }
        if upper[i as usize] > 1e20 {
            upper[i as usize] = f64::INFINITY;
        }
    }

    // permute LP rhs
    // EQ or BOUND first
    {
        let mut j = 0;
        for i in 0..*nRows {
            if constraint_type[i as usize] == EQ {
                rhs[j as usize] = lhs_clp[i as usize];
                constraint_new_idx[i as usize] = j;
                j += 1;
            } else if constraint_type[i as usize] == BOUND {
                rhs[j as usize] = 0.0;
                constraint_new_idx[i as usize] = j;
                j += 1;
            }
        }
    }
    // then LEQ or GEQ
    {
        let mut j = *nEqs;
        for i in 0..*nRows {
            if constraint_type[i as usize] == LEQ {
                rhs[j as usize] = -rhs_clp[i as usize]; // multiply -1
                constraint_new_idx[i as usize] = j;
                j += 1;
            } else if constraint_type[i as usize] == GEQ {
                rhs[j as usize] = lhs_clp[i as usize];
                constraint_new_idx[i as usize] = j;
                j += 1;
            }
        }
    }

    // formulate and permute LP matrix
    // beg remains the same
    for i in 0..nCols_clp + 1 {
        csc_beg[i as usize] = A_csc_beg[i as usize];
    }
    for i in nCols_clp + 1..*nCols + 1 {
        csc_beg[i as usize] = csc_beg[(i - 1) as usize] + 1;
    }

    // row idx changes
    {
        let mut k = 0;
        for i in 0..nCols_clp {
            // same order as in rhs
            // EQ or BOUND first
            for j in csc_beg[i as usize]..csc_beg[(i + 1) as usize] {
                if constraint_type[A_csc_idx[j as usize] as usize] == EQ
                    || constraint_type[A_csc_idx[j as usize] as usize] == BOUND
                {
                    csc_idx[k as usize] = constraint_new_idx[A_csc_idx[j as usize] as usize];
                    csc_val[k as usize] = A_csc_val[j as usize];
                    k += 1;
                }
            }
            // then LEQ or GEQ
            for j in csc_beg[i as usize]..csc_beg[(i + 1) as usize] {
                if constraint_type[A_csc_idx[j as usize] as usize] == LEQ {
                    csc_idx[k as usize] = constraint_new_idx[A_csc_idx[j as usize] as usize];
                    csc_val[k as usize] = -A_csc_val[j as usize]; // multiply -1
                    k += 1;
                } else if constraint_type[A_csc_idx[j as usize] as usize] == GEQ {
                    csc_idx[k as usize] = constraint_new_idx[A_csc_idx[j as usize] as usize];
                    csc_val[k as usize] = A_csc_val[j as usize];
                    k += 1;
                }
            }
        }
    }

    // slacks for BOUND
    {
        let mut j = nCols_clp;
        for i in 0..*nRows {
            if constraint_type[i as usize] == BOUND {
                csc_idx[csc_beg[j as usize] as usize] = constraint_new_idx[i as usize];
                csc_val[csc_beg[j as usize] as usize] = -1.0;
                j += 1;
            }
        }
    }

    retcode
}

/// HiGHS: CupdlpWrapper.cpp::problem_create (lines 450-456)
pub fn problem_create(prob: &mut Option<Box<CUPDLPproblem>>) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    // cupdlp_init_problem(*prob, 1);
    *prob = Some(Box::new(CUPDLPproblem::default()));

    retcode
}

/// HiGHS: CupdlpWrapper.cpp::data_alloc (lines 466-515)
///
/// CPU build: `data->device = CPU`; `dst_matrix_format` CSR_CSC builds both
/// the CSC copy and the CSR via `csc2csr`.
pub fn data_alloc(
    data: &mut CUPDLPdata,
    nRows: cupdlp_int,
    nCols: cupdlp_int,
    matrix: CupdlpSrcMatrix<'_>,
    dst_matrix_format: CUPDLP_MATRIX_FORMAT,
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    data.nRows = nRows;
    data.nCols = nCols;
    data.matrix_format = dst_matrix_format;
    data.dense_matrix = None;
    data.csr_matrix = None;
    data.csc_matrix = None;
    data.device = CPU;

    match dst_matrix_format {
        DENSE => {
            dense_create(&mut data.dense_matrix);
            dense_alloc_matrix(
                data.dense_matrix.as_deref_mut().unwrap(),
                nRows,
                nCols,
                matrix,
            );
        }
        CSR => {
            csr_create(&mut data.csr_matrix);
            csr_alloc_matrix(
                data.csr_matrix.as_deref_mut().unwrap(),
                nRows,
                nCols,
                matrix,
            );
        }
        CSC => {
            csc_create(&mut data.csc_matrix);
            csc_alloc_matrix(
                data.csc_matrix.as_deref_mut().unwrap(),
                nRows,
                nCols,
                matrix,
            );
        }
        CSR_CSC => {
            csc_create(&mut data.csc_matrix);
            csc_alloc_matrix(
                data.csc_matrix.as_deref_mut().unwrap(),
                nRows,
                nCols,
                matrix,
            );
            csr_create(&mut data.csr_matrix);
            csr_alloc_matrix(
                data.csr_matrix.as_deref_mut().unwrap(),
                nRows,
                nCols,
                matrix,
            );
        }
        _ => {}
    }
    // currently, only supprot that input matrix is CSC, and store both CSC and
    // CSR data->csc_matrix = matrix;

    retcode
}

/// HiGHS: CupdlpWrapper.cpp::problem_alloc (lines 517-585)
///
/// Copies cost / rhs / lower / upper, sets `hasLower` / `hasUpper` with
/// `cupdlp_hasLower(.., -INFINITY, ..)` / `cupdlp_hasUpper(.., +INFINITY,
/// ..)` and `data->csc_matrix->MatElemNormInf` with `infNorm` of the source
/// CSC values. `dMaxCost` / `dMaxRhs` / `dMaxRowBound` are not set (upstream
/// TODO, line 581). The two timing out-parameters are written with
/// `getTimeStamp` differences.
#[allow(clippy::too_many_arguments)]
pub fn problem_alloc(
    prob: &mut CUPDLPproblem,
    nRows: cupdlp_int,
    nCols: cupdlp_int,
    nEqs: cupdlp_int,
    cost: &[cupdlp_float],
    offset: cupdlp_float,
    sign_origin: cupdlp_float,
    matrix: CupdlpSrcMatrix<'_>,
    dst_matrix_format: CUPDLP_MATRIX_FORMAT,
    rhs: &[cupdlp_float],
    lower: &[cupdlp_float],
    upper: &[cupdlp_float],
    alloc_matrix_time: &mut cupdlp_float,
    copy_vec_time: &mut cupdlp_float,
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    prob.nRows = nRows;
    prob.nCols = nCols;
    prob.nEqs = nEqs;
    // prob.data = cupdlp_NULL;  // Default already
    // prob.cost = cupdlp_NULL;
    prob.offset = offset;
    prob.sense_origin = sign_origin;
    // prob.rhs = cupdlp_NULL;
    // prob.lower = cupdlp_NULL;
    // prob.upper = cupdlp_NULL;

    let begin = getTimeStamp();

    // cupdlp_init_data(prob->data, 1);
    prob.data = CUPDLPdata::default();

    // cupdlp_init_vec_double(prob->cost, nCols);
    prob.cost = vec![0.0; nCols as usize];
    // cupdlp_init_vec_double(prob->rhs, nRows);
    prob.rhs = vec![0.0; nRows as usize];
    // cupdlp_init_vec_double(prob->lower, nCols);
    prob.lower = vec![0.0; nCols as usize];
    // cupdlp_init_vec_double(prob->upper, nCols);
    prob.upper = vec![0.0; nCols as usize];
    // cupdlp_init_zero_vec_double(prob->hasLower, nCols);
    prob.hasLower = vec![0.0; nCols as usize];
    // cupdlp_init_zero_vec_double(prob->hasUpper, nCols);
    prob.hasUpper = vec![0.0; nCols as usize];

    data_alloc(&mut prob.data, nRows, nCols, matrix, dst_matrix_format);
    *alloc_matrix_time = getTimeStamp() - begin;

    // C++ casts matrix to CUPDLPcsc* for MatElemNormInf
    let csc_src = match matrix {
        CupdlpSrcMatrix::Csc(csc) => csc,
        _ => unreachable!(),
    };
    prob.data.csc_matrix.as_mut().unwrap().MatElemNormInf =
        inf_norm(&csc_src.colMatElem, csc_src.nMatElem);

    let begin = getTimeStamp();
    // CPU branch:
    cupdlp_copy(&mut prob.cost, cost, nCols);
    cupdlp_copy(&mut prob.rhs, rhs, nRows);
    cupdlp_copy(&mut prob.lower, lower, nCols);
    cupdlp_copy(&mut prob.upper, upper, nCols);

    *copy_vec_time = getTimeStamp() - begin;

    // Keep
    cupdlp_hasLower(&mut prob.hasLower, &prob.lower, -f64::INFINITY, nCols);
    cupdlp_hasUpper(&mut prob.hasUpper, &prob.upper, f64::INFINITY, nCols);

    // TODO: cal dMaxCost, dMaxRhs, dMaxRowBound

    // exit_cleanup: (line 583) — no goto reaches it in C; label dropped.
    retcode
}

/// HiGHS: CupdlpWrapper.cpp::infNorm (lines 590-594)
///
/// `norm = std::max(fabs(x[i]), norm)` from `norm = 0` (named `inf_norm` here:
/// `linalg::infNorm` is a different function).
pub fn inf_norm(x: &[f64], n: cupdlp_int) -> f64 {
    let mut norm = 0.0;
    for iX in 0..n {
        norm = cpp_max(x[iX as usize].abs(), norm);
    }
    norm
}

/// HiGHS: CupdlpWrapper.cpp::cupdlp_hasLower (lines 610-615)
pub fn cupdlp_hasLower(
    haslb: &mut [cupdlp_float],
    lb: &[cupdlp_float],
    bound: cupdlp_float,
    len: cupdlp_int,
) {
    for i in 0..len {
        haslb[i as usize] = if lb[i as usize] > bound { 1.0 } else { 0.0 };
    }
}

/// HiGHS: CupdlpWrapper.cpp::cupdlp_hasUpper (lines 617-622)
pub fn cupdlp_hasUpper(
    hasub: &mut [cupdlp_float],
    ub: &[cupdlp_float],
    bound: cupdlp_float,
    len: cupdlp_int,
) {
    for i in 0..len {
        hasub[i as usize] = if ub[i as usize] < bound { 1.0 } else { 0.0 };
    }
}

/// HiGHS: CupdlpWrapper.cpp::cupdlp_haslb (lines 624-631)
///
/// Declared and defined but never called in HiGHS.
pub fn cupdlp_haslb(
    haslb: &mut [cupdlp_float],
    lb: &[cupdlp_float],
    bound: cupdlp_float,
    len: cupdlp_int,
) {
    // CPU branch (CUPDLP_CPU defined):
    cupdlp_hasLower(haslb, lb, bound, len);
}

/// HiGHS: CupdlpWrapper.cpp::cupdlp_hasub (lines 633-640)
///
/// Declared and defined but never called in HiGHS.
pub fn cupdlp_hasub(
    hasub: &mut [cupdlp_float],
    ub: &[cupdlp_float],
    bound: cupdlp_float,
    len: cupdlp_int,
) {
    // CPU branch (CUPDLP_CPU defined):
    cupdlp_hasUpper(hasub, ub, bound, len);
}

/// HiGHS: CupdlpWrapper.cpp::getUserParamsFromOptions (lines 642-717, `static`)
/// HiGHS 1.16 (#3450): CupdlpWrapper.cpp::getUserParamsFromOptions (lines 643-718)
///
/// `N_ITER_LIM = min(pdlp_iteration_limit, kHighsIInf32)`; `N_LOG_LEVEL =
/// getCupdlpLogLevel`; `IF_SCALING` / `E_LINE_SEARCH_METHOD` /
/// `E_RESTART_METHOD` from `pdlp_features_off` (restart also off when
/// `pdlp_cupdlpc_restart_method == 0`); `D_PRIMAL_TOL` / `D_DUAL_TOL` /
/// `D_GAP_TOL` from the feasibility and `pdlp_optimality_tolerance`
/// options, all three overwritten by `kkt_tolerance` when it is not
/// `kDefaultKktTolerance`. `D_TIME_LIM = options.time_limit` on the 1.15.1
/// path: upstream computes the remaining time with `timer.read()` into a
/// local and then does not use it; the `timer.read()` call is kept (it has
/// no effect on state). With `h16_cupdlp_time` set (1.16.0), `D_TIME_LIM` is
/// that remaining time, `max(0, time_limit - timer.read())` when the limit
/// is finite. The `highsLogUser` lines are emitted as upstream.
pub fn get_user_params_from_options(
    options: &HighsOptions,
    timer: &mut HighsTimer,
    ifChangeIntParam: &mut [cupdlp_bool],
    intParam: &mut [cupdlp_int],
    ifChangeFloatParam: &mut [cupdlp_bool],
    floatParam: &mut [cupdlp_float],
) {
    for i in 0..N_INT_USER_PARAM {
        ifChangeIntParam[i] = false;
    }
    for i in 0..N_FLOAT_USER_PARAM {
        ifChangeFloatParam[i] = false;
    }
    // Assume all PDLP-related options in HiGHS cause changes
    ifChangeIntParam[N_ITER_LIM] = true;
    // If HiGHS is using 64-bit integers, then the default value of
    // options.pdlp_iteration_limit is kHighsIInf, so copying this to
    // intParam[N_ITER_LIM] will overflow.
    // HighsInt is i32 in this port, so the 64-bit clamp (C++ lines 655-657)
    // can never fire: the value is copied as is.
    intParam[N_ITER_LIM] = options.pdlp_iteration_limit;
    //
    ifChangeIntParam[N_LOG_LEVEL] = true;
    intParam[N_LOG_LEVEL] = get_cupdlp_log_level(options);
    //
    ifChangeIntParam[IF_SCALING] = true;
    let scaling_on: cupdlp_int = if (options.pdlp_features_off & kPdlpScalingOff) == 0 {
        1
    } else {
        0
    };
    intParam[IF_SCALING] = scaling_on;
    if scaling_on == 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "PDLP: Scaling off\n"
        );
    }
    //
    ifChangeIntParam[E_LINE_SEARCH_METHOD] = true;
    let adaptive_lineasearch: cupdlp_int =
        if (options.pdlp_features_off & kPdlpAdaptiveStepSizeOff) == 0 {
            1
        } else {
            0
        };
    intParam[E_LINE_SEARCH_METHOD] = adaptive_lineasearch;
    if adaptive_lineasearch == 1 {
        intParam[E_LINE_SEARCH_METHOD] = PDHG_ADAPTIVE_LINESEARCH;
    } else {
        intParam[E_LINE_SEARCH_METHOD] = PDHG_FIXED_LINESEARCH;
    }
    if adaptive_lineasearch == 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "PDLP: Adaptive line search off\n"
        );
    }
    //
    ifChangeFloatParam[D_PRIMAL_TOL] = true;
    floatParam[D_PRIMAL_TOL] = options.primal_feasibility_tolerance;
    //
    ifChangeFloatParam[D_DUAL_TOL] = true;
    floatParam[D_DUAL_TOL] = options.dual_feasibility_tolerance;
    //
    ifChangeFloatParam[D_GAP_TOL] = true;
    floatParam[D_GAP_TOL] = options.pdlp_optimality_tolerance;
    // Possibly over-write with uniform KKT tolerance
    if options.kkt_tolerance != kDefaultKktTolerance {
        floatParam[D_PRIMAL_TOL] = options.kkt_tolerance;
        floatParam[D_DUAL_TOL] = options.kkt_tolerance;
        floatParam[D_GAP_TOL] = options.kkt_tolerance;
    }
    //
    // cuPDLP-C has its own timer, so set its time limit according to
    // the time remaining with respect to the HiGHS time limit (if
    // finite)
    let mut time_limit = options.time_limit;
    if time_limit < kHighsInf {
        time_limit -= timer.read();
        time_limit = cpp_max(0.0, time_limit);
    }
    // 1.15.1 C++ computes the remaining time above and never reads it (lines
    // 701-705); 1.16.0 (#3450) reads it at line 708 (gated below).
    let _ = time_limit;
    ifChangeFloatParam[D_TIME_LIM] = true;
    // HiGHS 1.16 #3450 (h16_cupdlp_time): CupdlpWrapper.cpp line 708 (1.16.0) uses the
    // remaining time; the 1.15.1 path (below) passes options.time_limit.
    floatParam[D_TIME_LIM] = if options.h16_cupdlp_time {
        time_limit
    } else {
        options.time_limit
    };
    //
    ifChangeIntParam[E_RESTART_METHOD] = true;
    let mut restart_on: cupdlp_int = if (options.pdlp_features_off & kPdlpRestartOff) == 0 {
        1
    } else {
        0
    };
    if options.pdlp_cupdlpc_restart_method == 0 {
        restart_on = 0;
    }
    intParam[E_RESTART_METHOD] = restart_on;
    if restart_on == 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "PDLP: Restart off\n"
        );
    }
}

/// HiGHS: CupdlpWrapper.cpp::getCupdlpLogLevel (lines 839-848)
pub fn get_cupdlp_log_level(options: &HighsOptions) -> cupdlp_int {
    if options.output_flag {
        if options.log_dev_level != 0 {
            return 2;
        } else {
            return 1;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn log_level_follows_output_flag() {
        let mut o = HighsOptions::default();
        o.output_flag = false;
        assert_eq!(get_cupdlp_log_level(&o), 0);
        o.output_flag = true;
        o.log_dev_level = 0;
        assert_eq!(get_cupdlp_log_level(&o), 1);
        o.log_dev_level = 1;
        assert_eq!(get_cupdlp_log_level(&o), 2);
    }

    #[test]
    fn user_params_from_default_options() {
        let mut o = HighsOptions::default();
        o.output_flag = false;
        let mut timer = HighsTimer::default();
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        get_user_params_from_options(&o, &mut timer, &mut ici, &mut ip, &mut icf, &mut fp);
        assert!(ici[N_ITER_LIM]);
        assert_eq!(ip[N_ITER_LIM], kHighsIInf32);
        assert_eq!(ip[IF_SCALING], 1);
        assert_eq!(ip[E_LINE_SEARCH_METHOD], PDHG_ADAPTIVE_LINESEARCH);
        assert_eq!(ip[E_RESTART_METHOD], 1);
        assert_eq!(fp[D_GAP_TOL], o.pdlp_optimality_tolerance);
    }

    /// One ranged row 1 <= x0 + x1 <= 3 becomes a `>=` row with a slack
    /// column: nCols = 3, constraint_type = [3].
    #[test]
    fn ranged_row_gets_slack_column() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![3.0];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        let (mut cost, mut beg, mut idx, mut val) = (vec![], vec![], vec![], vec![]);
        let (mut rhs, mut lo, mut up, mut cni, mut ct) = (vec![], vec![], vec![], vec![], vec![]);
        let (mut nc, mut nr, mut nnz, mut ne, mut nco) = (0, 0, 0, 0, 0);
        let (mut off, mut sense) = (0.0, 0.0);
        let rc = formulate_lp_highs(
            0, &lp, &mut cost, &mut nc, &mut nr, &mut nnz, &mut ne, &mut beg, &mut idx, &mut val,
            &mut rhs, &mut lo, &mut up, &mut off, &mut sense, &mut nco, &mut cni, &mut ct,
        );
        assert_eq!(rc, 0);
        assert_eq!((nc, nr, nco, ne), (3, 1, 2, 1));
        assert_eq!(ct, vec![3]);
        assert_eq!(sense, 1.0);
    }

    #[test]
    fn inf_norm_basic() {
        let x = vec![-3.0, 1.0, -5.0, 2.0];
        assert_eq!(inf_norm(&x, 4), 5.0);
        assert_eq!(inf_norm(&x, 0), 0.0);
    }

    #[test]
    fn has_lower_upper() {
        let lb = vec![-f64::INFINITY, 0.0, 1.0];
        let mut haslb = vec![0.0; 3];
        cupdlp_hasLower(&mut haslb, &lb, -f64::INFINITY, 3);
        assert_eq!(haslb, vec![0.0, 1.0, 1.0]);

        let ub = vec![f64::INFINITY, 5.0, 10.0];
        let mut hasub = vec![0.0; 3];
        cupdlp_hasUpper(&mut hasub, &ub, f64::INFINITY, 3);
        assert_eq!(hasub, vec![0.0, 1.0, 1.0]);
    }

    #[test]
    fn haslb_hasub_cpu_branch() {
        // Dead CPU code: cupdlp_haslb / cupdlp_hasub call hasLower / hasUpper.
        let lb = vec![-f64::INFINITY, 0.0, 1.0];
        let mut haslb = vec![0.0; 3];
        cupdlp_haslb(&mut haslb, &lb, -f64::INFINITY, 3);
        assert_eq!(haslb, vec![0.0, 1.0, 1.0]);

        let ub = vec![f64::INFINITY, 5.0, 10.0];
        let mut hasub = vec![0.0; 3];
        cupdlp_hasub(&mut hasub, &ub, f64::INFINITY, 3);
        assert_eq!(hasub, vec![0.0, 1.0, 1.0]);
    }

    #[test]
    fn problem_create_sets_default() {
        let mut prob: Option<Box<CUPDLPproblem>> = None;
        let rc = problem_create(&mut prob);
        assert_eq!(rc, RETCODE_OK);
        assert!(prob.is_some());
    }

    #[test]
    fn data_alloc_csc() {
        let csc = CUPDLPcsc {
            nRows: 2,
            nCols: 3,
            nMatElem: 4,
            colMatBeg: vec![0, 2, 3, 4],
            colMatIdx: vec![0, 1, 0, 1],
            colMatElem: vec![1.0, 2.0, 3.0, 4.0],
            MatElemNormInf: 0.0,
        };
        let mut data = CUPDLPdata::default();
        let rc = data_alloc(&mut data, 2, 3, CupdlpSrcMatrix::Csc(&csc), CSC);
        assert_eq!(rc, RETCODE_OK);
        assert_eq!(data.nRows, 2);
        assert_eq!(data.nCols, 3);
        assert_eq!(data.device, CPU);
        assert!(data.csc_matrix.is_some());
        assert!(data.csr_matrix.is_none());
    }

    #[test]
    fn problem_alloc_basic() {
        let csc = CUPDLPcsc {
            nRows: 2,
            nCols: 3,
            nMatElem: 4,
            colMatBeg: vec![0, 2, 3, 4],
            colMatIdx: vec![0, 1, 0, 1],
            colMatElem: vec![1.0, 2.0, 3.0, 4.0],
            MatElemNormInf: 0.0,
        };
        let cost = vec![1.0, 2.0, 3.0];
        let rhs = vec![10.0, 20.0];
        let lower = vec![0.0, 0.0, 0.0];
        let upper = vec![f64::INFINITY, f64::INFINITY, f64::INFINITY];
        let mut prob = CUPDLPproblem::default();
        let mut alloc_time = 0.0;
        let mut copy_time = 0.0;
        let rc = problem_alloc(
            &mut prob,
            2,
            3,
            0,
            &cost,
            0.0,
            1.0,
            CupdlpSrcMatrix::Csc(&csc),
            CSC,
            &rhs,
            &lower,
            &upper,
            &mut alloc_time,
            &mut copy_time,
        );
        assert_eq!(rc, RETCODE_OK);
        assert_eq!(prob.nRows, 2);
        assert_eq!(prob.nCols, 3);
        assert_eq!(prob.cost, vec![1.0, 2.0, 3.0]);
        assert_eq!(prob.rhs, vec![10.0, 20.0]);
        assert_eq!(prob.lower, vec![0.0, 0.0, 0.0]);
        assert_eq!(
            prob.upper,
            vec![f64::INFINITY, f64::INFINITY, f64::INFINITY]
        );
        assert_eq!(prob.hasLower, vec![1.0, 1.0, 1.0]);
        assert_eq!(prob.hasUpper, vec![0.0, 0.0, 0.0]);
        assert_eq!(prob.data.csc_matrix.as_ref().unwrap().MatElemNormInf, 4.0);
    }

    #[test]
    fn formulate_eq_row() {
        // EQ row: lhs == rhs, both finite
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![kHighsInf];
        lp.row_lower_ = vec![5.0];
        lp.row_upper_ = vec![5.0];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let (mut cost, mut beg, mut idx, mut val) = (vec![], vec![], vec![], vec![]);
        let (mut rhs, mut lo, mut up, mut cni, mut ct) = (vec![], vec![], vec![], vec![], vec![]);
        let (mut nc, mut nr, mut nnz, mut ne, mut nco) = (0, 0, 0, 0, 0);
        let (mut off, mut sense) = (0.0, 0.0);
        let rc = formulate_lp_highs(
            0, &lp, &mut cost, &mut nc, &mut nr, &mut nnz, &mut ne, &mut beg, &mut idx, &mut val,
            &mut rhs, &mut lo, &mut up, &mut off, &mut sense, &mut nco, &mut cni, &mut ct,
        );
        assert_eq!(rc, 0);
        assert_eq!(ct, vec![EQ]);
        assert_eq!(nc, 1); // no extra column for EQ
        assert_eq!(ne, 1);
        assert_eq!(rhs, vec![5.0]);
    }

    #[test]
    fn formulate_leq_row() {
        // LEQ row: no lower, finite upper
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![kHighsInf];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![10.0];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let (mut cost, mut beg, mut idx, mut val) = (vec![], vec![], vec![], vec![]);
        let (mut rhs, mut lo, mut up, mut cni, mut ct) = (vec![], vec![], vec![], vec![], vec![]);
        let (mut nc, mut nr, mut nnz, mut ne, mut nco) = (0, 0, 0, 0, 0);
        let (mut off, mut sense) = (0.0, 0.0);
        let rc = formulate_lp_highs(
            0, &lp, &mut cost, &mut nc, &mut nr, &mut nnz, &mut ne, &mut beg, &mut idx, &mut val,
            &mut rhs, &mut lo, &mut up, &mut off, &mut sense, &mut nco, &mut cni, &mut ct,
        );
        assert_eq!(rc, 0);
        assert_eq!(ct, vec![LEQ]);
        assert_eq!(nc, 1);
        assert_eq!(ne, 0);
        assert_eq!(rhs, vec![-10.0]); // negated
    }

    #[test]
    fn formulate_geq_row() {
        // GEQ row: finite lower, no upper
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![kHighsInf];
        lp.row_lower_ = vec![3.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let (mut cost, mut beg, mut idx, mut val) = (vec![], vec![], vec![], vec![]);
        let (mut rhs, mut lo, mut up, mut cni, mut ct) = (vec![], vec![], vec![], vec![], vec![]);
        let (mut nc, mut nr, mut nnz, mut ne, mut nco) = (0, 0, 0, 0, 0);
        let (mut off, mut sense) = (0.0, 0.0);
        let rc = formulate_lp_highs(
            0, &lp, &mut cost, &mut nc, &mut nr, &mut nnz, &mut ne, &mut beg, &mut idx, &mut val,
            &mut rhs, &mut lo, &mut up, &mut off, &mut sense, &mut nco, &mut cni, &mut ct,
        );
        assert_eq!(rc, 0);
        assert_eq!(ct, vec![GEQ]);
        assert_eq!(nc, 1);
        assert_eq!(ne, 0);
        assert_eq!(rhs, vec![3.0]);
    }

    #[test]
    fn formulate_maximize_sense() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.sense_ = ObjSense::kMaximize;
        lp.col_cost_ = vec![2.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![kHighsInf];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let (mut cost, mut beg, mut idx, mut val) = (vec![], vec![], vec![], vec![]);
        let (mut rhs, mut lo, mut up, mut cni, mut ct) = (vec![], vec![], vec![], vec![], vec![]);
        let (mut nc, mut nr, mut nnz, mut ne, mut nco) = (0, 0, 0, 0, 0);
        let (mut off, mut sense) = (0.0, 0.0);
        formulate_lp_highs(
            0, &lp, &mut cost, &mut nc, &mut nr, &mut nnz, &mut ne, &mut beg, &mut idx, &mut val,
            &mut rhs, &mut lo, &mut up, &mut off, &mut sense, &mut nco, &mut cni, &mut ct,
        );
        assert_eq!(sense, -1.0);
        assert_eq!(cost[0], -2.0); // 2.0 * -1.0
    }

    #[test]
    fn user_params_kkt_override() {
        let mut o = HighsOptions::default();
        o.output_flag = false;
        o.kkt_tolerance = 1e-9; // not kDefaultKktTolerance
        let mut timer = HighsTimer::default();
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        get_user_params_from_options(&o, &mut timer, &mut ici, &mut ip, &mut icf, &mut fp);
        assert_eq!(fp[D_PRIMAL_TOL], 1e-9);
        assert_eq!(fp[D_DUAL_TOL], 1e-9);
        assert_eq!(fp[D_GAP_TOL], 1e-9);
    }

    #[test]
    fn user_params_restart_off() {
        let mut o = HighsOptions::default();
        o.output_flag = false;
        o.pdlp_cupdlpc_restart_method = 0;
        let mut timer = HighsTimer::default();
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        get_user_params_from_options(&o, &mut timer, &mut ici, &mut ip, &mut icf, &mut fp);
        assert_eq!(ip[E_RESTART_METHOD], 0);
    }

    #[test]
    fn user_params_scaling_off() {
        let mut o = HighsOptions::default();
        o.output_flag = false;
        o.pdlp_features_off = kPdlpScalingOff;
        let mut timer = HighsTimer::default();
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        get_user_params_from_options(&o, &mut timer, &mut ici, &mut ip, &mut icf, &mut fp);
        assert_eq!(ip[IF_SCALING], 0);
    }

    #[test]
    fn user_params_adaptive_off() {
        let mut o = HighsOptions::default();
        o.output_flag = false;
        o.pdlp_features_off = kPdlpAdaptiveStepSizeOff;
        let mut timer = HighsTimer::default();
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        get_user_params_from_options(&o, &mut timer, &mut ici, &mut ip, &mut icf, &mut fp);
        assert_eq!(ip[E_LINE_SEARCH_METHOD], PDHG_FIXED_LINESEARCH);
    }

    /// A timer whose clock started 4 seconds ago (`HighsTimer`'s fields are
    /// public, so the start instant can be set directly).
    fn timer_started_seconds_ago(secs: u64) -> HighsTimer {
        HighsTimer {
            start_time: std::time::Instant::now() - std::time::Duration::from_secs(secs),
            running: true,
            accumulated: 0.0,
            stopped: false,
        }
    }

    /// HiGHS 1.16 #3450 (h16_cupdlp_time): with the flag on, D_TIME_LIM is the
    /// remaining time, time_limit - timer.read() = 10 - 4 = 6 (1.16.0
    /// CupdlpWrapper.cpp line 708).
    #[test]
    fn h16_cupdlp_time_remaining_limit_when_on() {
        let mut o = HighsOptions::default();
        o.output_flag = false;
        o.time_limit = 10.0;
        o.h16_cupdlp_time = true;
        let mut timer = timer_started_seconds_ago(4);
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        get_user_params_from_options(&o, &mut timer, &mut ici, &mut ip, &mut icf, &mut fp);
        assert!(icf[D_TIME_LIM]);
        // timer.read() is slightly more than 4 s, so allow a small tolerance.
        assert!(
            (fp[D_TIME_LIM] - 6.0).abs() < 1e-2,
            "D_TIME_LIM = {}",
            fp[D_TIME_LIM]
        );
    }

    /// HiGHS 1.16 #3450 (h16_cupdlp_time): with the flag off (1.15.1 path), D_TIME_LIM
    /// stays options.time_limit = 10 even though the timer has run for 4 s.
    #[test]
    fn h16_cupdlp_time_full_limit_when_off() {
        let mut o = HighsOptions::default();
        o.output_flag = false;
        o.time_limit = 10.0;
        o.h16_cupdlp_time = false;
        let mut timer = timer_started_seconds_ago(4);
        let mut ici = [false; N_INT_USER_PARAM];
        let mut ip = [0; N_INT_USER_PARAM];
        let mut icf = [false; N_FLOAT_USER_PARAM];
        let mut fp = [0.0; N_FLOAT_USER_PARAM];
        get_user_params_from_options(&o, &mut timer, &mut ici, &mut ip, &mut icf, &mut fp);
        assert_eq!(fp[D_TIME_LIM], 10.0);
    }
}
