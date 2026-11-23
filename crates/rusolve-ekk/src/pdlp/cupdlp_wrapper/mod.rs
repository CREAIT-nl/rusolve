//! HiGHS source: pdlp/CupdlpWrapper.h, pdlp/CupdlpWrapper.cpp (lines 1-278)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiGHS -> cuPDLP-C bridge (`solver=pdlp`). `CupdlpWrapper.cpp` is 848
//! lines, split by contiguous line range:
//!
//! | Rust module    | lines           | contents |
//! |----------------|-----------------|----------|
//! | this file      | 1-278           | `solveLpCupdlp` (both overloads) |
//! | [`formulate`]  | 280-717, 839-848 | `formulateLP_highs`, `problem_create`, `data_alloc`, `problem_alloc`, the bound helpers, `getUserParamsFromOptions`, `getCupdlpLogLevel` |
//!
//! Not ported: `analysePdlpSolution` (719-837) is called only under
//! `#if CUPDLP_DEBUG` (undefined); the `cupdlp_init_*` malloc macros of
//! `CupdlpWrapper.h` are `Vec` / `Box` allocation.
//!
//! These are HiGHS C++ functions, so they take snake_case names (repo
//! convention); the cuPDLP C functions they call keep their C names.
//!
//! # Ownership in `solveLpCupdlp`
//!
//! `cupdlp_init_work(w, 1)` + `w->problem = prob; w->scaling = scaling;
//! PDHG_Alloc(w);` becomes a `Box<CUPDLPwork>` that takes the problem and
//! the scaling by value. `highs_solution.col_value.data()` etc. reach
//! `LP_SolvePDHG` as `vec_data_mut` (`None` for an empty vector, as libc++
//! `data()` is `nullptr`), see [`crate::ipm::ipx_wrapper::vec_data_mut`].
//! The `free` / `*_clear` calls at the end are drops.

pub mod formulate;

pub use formulate::*;

use super::cupdlp::scaling::*;
use super::cupdlp::solver::*;
use super::cupdlp::utils::*;
use super::cupdlp::*;
use crate::ipm::ipx_wrapper::vec_data_mut;
use crate::lp_data::highs_solution::reset_model_status_and_highs_info;
use crate::prelude::*;

/// HiGHS: CupdlpWrapper.cpp::solveLpCupdlp(HighsLpSolverObject&) (lines 23-28)
pub fn solve_lp_cupdlp_so(solver_object: &mut HighsLpSolverObject) -> HighsStatus {
    solve_lp_cupdlp(
        &*solver_object.options_,
        &mut *solver_object.timer_,
        &*solver_object.lp_,
        &mut *solver_object.basis_,
        &mut *solver_object.solution_,
        &mut solver_object.model_status_,
        &mut *solver_object.highs_info_,
        &mut *solver_object.callback_,
    )
}

/// HiGHS: CupdlpWrapper.cpp::solveLpCupdlp (lines 30-278)
///
/// Order of state-changing calls: `resetModelStatusAndHighsInfo`,
/// `getUserParamsFromOptions`, `getCupdlpLogLevel`, `formulateLP_highs`,
/// `Init_Scaling`, `problem_create`, `csc_create` + copy of the CSC arrays,
/// `PDHG_Scale_Data`, `problem_alloc`, `PDHG_Alloc`, the
/// `dScalingTime` / `dPresolveTime` timers, copy of `scaling->rowScale` /
/// `colScale` into `w->rowScale` / `colScale`, resize of the four solution
/// vectors to `lp.num_col_` / `lp.num_row_`, `LP_SolvePDHG`,
/// `highs_info.pdlp_iteration_count = pdlp_num_iter`, `PDHG_Destroy`, then
/// the status mapping: `kUnknown`; `retcode == RETCODE_OK` maps OPTIMAL ->
/// `kOptimal`, INFEASIBLE -> `kInfeasible`, UNBOUNDED -> `kUnbounded`,
/// INFEASIBLE_OR_UNBOUNDED -> `kUnboundedOrInfeasible`,
/// TIMELIMIT_OR_ITERLIMIT -> `pdlp_num_iter >= intParam[N_ITER_LIM] - 1 ?
/// kIterationLimit : kTimeLimit`, FEASIBLE -> `assert(111 == 666)` then
/// `kUnknown` (a `debug_assert!`), anything else `assert(111 == 777)`;
/// otherwise `kSolveError`. Returns `kOk` iff `retcode == RETCODE_OK`.
/// `highs_basis.valid = false` always.
#[allow(clippy::too_many_arguments)]
pub fn solve_lp_cupdlp(
    options: &HighsOptions,
    timer: &mut HighsTimer,
    lp: &HighsLp,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
    model_status: &mut HighsModelStatus,
    highs_info: &mut HighsInfo,
    callback: &mut HighsCallback,
) -> HighsStatus {
    // Indicate that no imprecise solution has (yet) been found
    reset_model_status_and_highs_info(model_status, highs_info);

    let fp: Option<&str> = None; // char* fp = nullptr
    let fp_sol: Option<&str> = None; // char* fp_sol = nullptr

    let mut nCols: cupdlp_int = 0;
    let mut nRows: cupdlp_int = 0;
    let mut nEqs: cupdlp_int = 0;
    let mut nCols_origin: cupdlp_int = 0;
    let ifSaveSol: cupdlp_bool = false;
    // cupdlp_bool ifPresolve = false; (declared in C++, never read)

    let mut nnz: cupdlp_int = 0;

    let mut rhs: Vec<cupdlp_float> = Vec::new(); // double* rhs = NULL
    let mut cost: Vec<cupdlp_float> = Vec::new(); // double* cost = NULL

    let mut lower: Vec<cupdlp_float> = Vec::new(); // cupdlp_float* lower = NULL
    let mut upper: Vec<cupdlp_float> = Vec::new(); // cupdlp_float* upper = NULL

    // -------------------------
    let mut csc_beg: Vec<cupdlp_int> = Vec::new(); // int* csc_beg = NULL
    let mut csc_idx: Vec<cupdlp_int> = Vec::new(); // int* csc_idx = NULL
    let mut csc_val: Vec<cupdlp_float> = Vec::new(); // double* csc_val = NULL
    let mut offset: cupdlp_float = 0.0; // true objVal = sig * c'x - offset, sig = 1 (min) or -1 (max)
    let mut sense_origin: cupdlp_float = 1.0; // 1 (min) or -1 (max)
    let mut constraint_new_idx: Vec<cupdlp_int> = Vec::new(); // int* constraint_new_idx = NULL
    let mut constraint_type: Vec<cupdlp_int> = Vec::new(); // int* constraint_type = NULL

    // CUPDLPscaling* scaling = (CUPDLPscaling*)cupdlp_malloc(sizeof(CUPDLPscaling));
    // The WIP printf / ifRuizScaling lines (C++ 68-78) are dropped.
    let mut scaling = CUPDLPscaling::default();

    // claim solvers variables
    // prepare pointers
    // CUPDLP_MATRIX_FORMAT src_matrix_format = CSC; the tag travels with the
    // CupdlpSrcMatrix::Csc argument of problem_alloc.
    let dst_matrix_format: CUPDLP_MATRIX_FORMAT = CSR_CSC;
    let mut csc_cpu: Option<Box<CUPDLPcsc>> = None; // CUPDLPcsc* csc_cpu = cupdlp_NULL;
    let mut prob: Option<Box<CUPDLPproblem>> = None; // CUPDLPproblem* prob = cupdlp_NULL;

    // load parameters

    // set solver parameters
    let mut ifChangeIntParam = [false; N_INT_USER_PARAM];
    let mut intParam = [0; N_INT_USER_PARAM];
    let mut ifChangeFloatParam = [false; N_FLOAT_USER_PARAM];
    let mut floatParam = [0.0; N_FLOAT_USER_PARAM];

    // Transfer from options
    get_user_params_from_options(
        options,
        timer,
        &mut ifChangeIntParam,
        &mut intParam,
        &mut ifChangeFloatParam,
        &mut floatParam,
    );

    // std::vector<int> constraint_type(lp.num_row_);

    let local_log_level: cupdlp_int = get_cupdlp_log_level(options);
    // if (local_log_level) cupdlp_printf("Solving with cuPDLP-C\n"); dropped: print

    formulate_lp_highs(
        local_log_level,
        lp,
        &mut cost,
        &mut nCols,
        &mut nRows,
        &mut nnz,
        &mut nEqs,
        &mut csc_beg,
        &mut csc_idx,
        &mut csc_val,
        &mut rhs,
        &mut lower,
        &mut upper,
        &mut offset,
        &mut sense_origin,
        &mut nCols_origin,
        &mut constraint_new_idx,
        &mut constraint_type,
    );

    // H_Init_Scaling(local_log_level, scaling, nCols, nRows, cost, rhs); (commented out in C++)
    Init_Scaling(local_log_level, &mut scaling, nCols, nRows, &cost, &rhs);

    let ifScaling: cupdlp_int = intParam[IF_SCALING];

    // cupdlp_init_work(w, 1): the work struct is boxed; PDHG_Destroy takes the Option.
    let mut w: Option<Box<CUPDLPwork>> = Some(Box::new(CUPDLPwork::default())); // CUPDLPwork* w = cupdlp_NULL;

    // #ifdef CUPDLP_GPU block (cuSPARSE / cuBLAS handles and timers) dropped: CPU build.

    problem_create(&mut prob);

    // currently, only supprot that input matrix is CSC, and store both CSC and
    // CSR
    csc_create(&mut csc_cpu);
    let csc: &mut CUPDLPcsc = csc_cpu.as_deref_mut().unwrap();
    csc.nRows = nRows;
    csc.nCols = nCols;
    csc.nMatElem = nnz;
    // csc_cpu->colMatBeg = malloc((1 + nCols) * sizeof(int)); memcpy(..., csc_beg, (nCols + 1) * sizeof(int));
    csc.colMatBeg = csc_beg[..(1 + nCols) as usize].to_vec();
    // csc_cpu->colMatIdx = malloc(nnz * sizeof(int)); memcpy(..., csc_idx, nnz * sizeof(int));
    csc.colMatIdx = csc_idx[..nnz as usize].to_vec();
    // csc_cpu->colMatElem = malloc(nnz * sizeof(double)); memcpy(..., csc_val, nnz * sizeof(double));
    csc.colMatElem = csc_val[..nnz as usize].to_vec();

    let mut scaling_time: cupdlp_float = getTimeStamp();

    PDHG_Scale_Data(
        local_log_level,
        csc,
        ifScaling,
        &mut scaling,
        &mut cost,
        &mut lower,
        &mut upper,
        &mut rhs,
    );

    scaling_time = getTimeStamp() - scaling_time;

    let mut alloc_matrix_time: cupdlp_float = 0.0;
    let mut copy_vec_time: cupdlp_float = 0.0;

    problem_alloc(
        prob.as_deref_mut().unwrap(),
        nRows,
        nCols,
        nEqs,
        &cost,
        offset,
        sense_origin,
        CupdlpSrcMatrix::Csc(csc),
        dst_matrix_format,
        &rhs,
        &lower,
        &upper,
        &mut alloc_matrix_time,
        &mut copy_vec_time,
    );

    let wk: &mut CUPDLPwork = w.as_deref_mut().unwrap();
    wk.problem = *prob.take().unwrap(); // w->problem = prob;
    wk.scaling = scaling; // w->scaling = scaling;
    PDHG_Alloc(wk);
    wk.timers.dScalingTime = scaling_time;
    wk.timers.dPresolveTime = 0.0; // presolve_time (0 in C++: no presolve)
    // CUPDLP_CPU branch: cupdlp_copy_vec(w->rowScale, scaling->rowScale, cupdlp_float, nRows);
    cupdlp_copy(&mut wk.rowScale, &wk.scaling.rowScale, nRows);
    cupdlp_copy(&mut wk.colScale, &wk.scaling.colScale, nCols);

    // #ifdef CUPDLP_GPU timer / transfer block dropped: CPU build.

    // CUPDLP_CALL(LP_SolvePDHG(...)) commented out in C++ (lines 184-187).

    // Resize the highs_solution so cuPDLP-c can use it internally
    highs_solution.col_value.resize(lp.num_col_ as usize, 0.0);
    highs_solution.row_value.resize(lp.num_row_ as usize, 0.0);
    highs_solution.col_dual.resize(lp.num_col_ as usize, 0.0);
    highs_solution.row_dual.resize(lp.num_row_ as usize, 0.0);
    let mut value_valid: cupdlp_int = highs_solution.value_valid as cupdlp_int;
    let mut dual_valid: cupdlp_int = highs_solution.dual_valid as cupdlp_int;
    let mut pdlp_model_status: cupdlp_int = 0;
    let mut pdlp_num_iter: cupdlp_int = 0;

    let retcode: cupdlp_retcode = LP_SolvePDHG(
        wk,
        &ifChangeIntParam,
        &intParam,
        &ifChangeFloatParam,
        &floatParam,
        fp,
        nCols_origin,
        vec_data_mut(highs_solution.col_value.as_mut_slice()),
        vec_data_mut(highs_solution.col_dual.as_mut_slice()),
        vec_data_mut(highs_solution.row_value.as_mut_slice()),
        vec_data_mut(highs_solution.row_dual.as_mut_slice()),
        &mut value_valid,
        &mut dual_valid,
        ifSaveSol,
        fp_sol,
        Some(constraint_new_idx.as_slice()),
        Some(constraint_type.as_slice()),
        &mut pdlp_model_status,
        &mut pdlp_num_iter,
    );
    highs_info.pdlp_iteration_count = pdlp_num_iter;

    // #if PDLP_DEBUG_LOG debugPdlpFinalSolutionLog / fclose dropped (PDLP_DEBUG_LOG == 0).
    PDHG_Destroy(&mut w);

    *model_status = HighsModelStatus::kUnknown;
    highs_solution.value_valid = value_valid != 0;
    highs_solution.dual_valid = dual_valid != 0;
    highs_basis.valid = false;

    if retcode == RETCODE_OK {
        if pdlp_model_status == OPTIMAL {
            *model_status = HighsModelStatus::kOptimal;
        } else if pdlp_model_status == INFEASIBLE {
            *model_status = HighsModelStatus::kInfeasible;
        } else if pdlp_model_status == UNBOUNDED {
            *model_status = HighsModelStatus::kUnbounded;
        } else if pdlp_model_status == INFEASIBLE_OR_UNBOUNDED {
            *model_status = HighsModelStatus::kUnboundedOrInfeasible;
        } else if pdlp_model_status == TIMELIMIT_OR_ITERLIMIT {
            *model_status = if pdlp_num_iter >= intParam[N_ITER_LIM] - 1 {
                HighsModelStatus::kIterationLimit
            } else {
                HighsModelStatus::kTimeLimit
            };
        } else if pdlp_model_status == FEASIBLE {
            debug_assert_eq!(111, 666); // C++ assert(111 == 666)
            *model_status = HighsModelStatus::kUnknown;
        } else {
            debug_assert_eq!(111, 777); // C++ assert(111 == 777)
        }
        // #if CUPDLP_DEBUG analysePdlpSolution(...) dropped (CUPDLP_DEBUG undefined).
    } else {
        // Failure return from LP_SolvePDHG
        *model_status = HighsModelStatus::kSolveError;
    }
    // free problem: scaling_clear, the cupdlp_free calls, csc_clear_host and
    // problem_clear are Drop in Rust. #if !defined(CUPDLP_CPU) cudaDeviceReset dropped.

    if retcode == RETCODE_OK {
        HighsStatus::kOk
    } else {
        HighsStatus::kError
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests. End-to-end parity is the
    //! oracle test `tests/t_w3_cupdlp.rs`.
    use super::*;

    /// min x s.t. x >= 1, 0 <= x <= 10: cuPDLP reaches kOptimal near x = 1.
    #[test]
    fn one_variable_lp_is_optimal() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let mut options = HighsOptions::default();
        options.output_flag = false;
        let mut timer = HighsTimer::default();
        let mut basis = HighsBasis::default();
        let mut solution = HighsSolution::default();
        let mut status = HighsModelStatus::kNotset;
        let mut info = HighsInfo::default();
        let mut callback = HighsCallback::default();
        let rc = solve_lp_cupdlp(
            &options,
            &mut timer,
            &lp,
            &mut basis,
            &mut solution,
            &mut status,
            &mut info,
            &mut callback,
        );
        assert_eq!(rc, HighsStatus::kOk);
        assert_eq!(status, HighsModelStatus::kOptimal);
        assert!((solution.col_value[0] - 1.0).abs() < 1e-4);
        assert!(!basis.valid);
        assert!(info.pdlp_iteration_count > 0);
    }
}
