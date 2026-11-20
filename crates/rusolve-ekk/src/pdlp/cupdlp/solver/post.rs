//! HiGHS source: pdlp/cupdlp/cupdlp_solver.c (lines 1217-1498)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Warm-start mapping, solution recovery and the cuPDLP entry point.
//!
//! The nullable C pointers become `Option`s: `col_value` / `col_dual` /
//! `row_value` / `row_dual` are `Option<&mut [cupdlp_float]>` (`if
//! (col_value)` tests), `constraint_new_idx` / `constraint_type` are
//! `Option<&[cupdlp_int]>`, `value_valid` / `dual_valid` are `&mut
//! cupdlp_int` (HiGHS always passes them; `PDHG_PreSolve` dereferences them
//! unconditionally). Where `PDHG_PreSolve` dereferences an optional pointer
//! unconditionally, `None` is read as an empty slice, so it panics only on
//! an index C would read through NULL. The `CUPDLP_INIT_DOUBLE` buffers are local `Vec`s;
//! their allocation failure path (`goto exit_cleanup`) has no Rust
//! equivalent.

use super::super::defs::*;
use super::super::linalg::*;
use super::super::utils::*;
use super::*;

/// HiGHS: cupdlp_solver.c::PDHG_PreSolve (lines 1217-1279)
///
/// Returns at once unless both `*value_valid` and `*dual_valid`. Fills
/// `x[0]` / `y[0]`: `row_buffer[constraint_new_idx[iRow]] = sense * mu *
/// row_dual[iRow]` with `mu = constraint_type[iRow] == 1 ? -1 : 1`, the slack
/// of a `constraint_type == 3` row is `col_buffer[iCol++] = row_value[iRow]`;
/// then `edot` by `colScale` / `rowScale` when `scaling->ifScaled`.
#[allow(clippy::too_many_arguments)]
pub fn PDHG_PreSolve(
    pdhg: &mut CUPDLPwork,
    nCols_origin: cupdlp_int,
    constraint_new_idx: Option<&[cupdlp_int]>,
    constraint_type: Option<&[cupdlp_int]>,
    col_value: Option<&mut [cupdlp_float]>,
    col_dual: Option<&mut [cupdlp_float]>,
    row_value: Option<&mut [cupdlp_float]>,
    row_dual: Option<&mut [cupdlp_float]>,
    value_valid: &mut cupdlp_int,
    dual_valid: &mut cupdlp_int,
) -> cupdlp_retcode {
    let mut retcode: cupdlp_retcode = RETCODE_OK;

    if *value_valid == 0 || *dual_valid == 0 {
        // Nothing to presolve
        return retcode;
    }

    let problem = &pdhg.problem;
    let iterates = &mut pdhg.iterates;
    let scaling = &pdhg.scaling;
    // CUPDLPresobj *resobj = pdhg->resobj; // unused in this function
    let sense = problem.sense_origin;

    let nCols = problem.nCols;
    let nRows = problem.nRows;

    // Borrow x[0] and y[0] together
    let [x0, _x1] = &mut iterates.x;
    let [y0, _y1] = &mut iterates.y;

    debug_assert!(constraint_new_idx.is_some());
    debug_assert!(constraint_type.is_some());

    // Have to build cuPDLP-C solution in CPU row and column buffers, so
    // that only copy operations are required to put the cuPDLP-C
    // solution into GPU buffers
    // CUPDLP_INIT_DOUBLE(col_buffer, problem->nCols);
    // CUPDLP_INIT_DOUBLE(row_buffer, problem->nRows);
    let mut col_buffer = vec![0.0f64; nCols as usize];
    let mut row_buffer = vec![0.0f64; nRows as usize];

    // CUPDLP_ZERO_VEC(x->data, cupdlp_float, problem->nCols);
    cupdlp_zero(&mut x0.data, nCols);
    // CUPDLP_ZERO_VEC(y->data, cupdlp_float, problem->nRows);
    cupdlp_zero(&mut y0.data, nRows);

    // C dereferences these pointers unconditionally; a NULL pointer (an
    // empty HiGHS vector) becomes an empty slice, so it fails only where C
    // would actually read through it.
    let cv: &[cupdlp_float] = col_value.map_or(&[], |v| &*v);
    let cni = constraint_new_idx.unwrap_or_default();
    let ct = constraint_type.unwrap_or_default();
    let rd: &[cupdlp_float] = row_dual.map_or(&[], |v| &*v);
    let rv: &[cupdlp_float] = row_value.map_or(&[], |v| &*v);

    let mut iCol: cupdlp_int = 0;
    // for (; iCol < nCols_origin; iCol++)
    //   col_buffer[iCol] = col_value[iCol];
    while iCol < nCols_origin {
        col_buffer[iCol as usize] = cv[iCol as usize];
        iCol += 1;
    }

    // for (int iRow = 0; iRow < problem->nRows; iRow++) {
    //   const double mu = constraint_type[iRow] == 1 ? -1 : 1;
    //   row_buffer[constraint_new_idx[iRow]] = sense * mu * row_dual[iRow];
    //   if (constraint_type[iRow] == 3)
    //     col_buffer[iCol++] = row_value[iRow];
    // }
    for iRow in 0..nRows {
        let mu: cupdlp_float = if ct[iRow as usize] == 1 { -1.0 } else { 1.0 };
        row_buffer[cni[iRow as usize] as usize] = sense * mu * rd[iRow as usize];
        if ct[iRow as usize] == 3 {
            col_buffer[iCol as usize] = rv[iRow as usize];
            iCol += 1;
        }
    }

    // Copy into GPU buffers
    // CUPDLP_COPY_VEC(x->data, col_buffer, cupdlp_float, problem->nCols);
    cupdlp_copy(&mut x0.data, &col_buffer, nCols);
    // CUPDLP_COPY_VEC(y->data, row_buffer, cupdlp_float, problem->nRows);
    cupdlp_copy(&mut y0.data, &row_buffer, nRows);

    // Scale
    if scaling.ifScaled != 0 {
        cupdlp_edot(&mut x0.data, &pdhg.colScale, nCols);
        cupdlp_edot(&mut y0.data, &pdhg.rowScale, nRows);
    }

    // exit_cleanup:
    // CUPDLP_FREE(col_buffer);
    // CUPDLP_FREE(row_buffer);
    // Buffers are Vecs, dropped automatically.

    retcode
}

/// HiGHS: cupdlp_solver.c::PDHG_PostSolve (lines 1281-1435)
///
/// Unscales the iterate `k = nIter % 2` in place when `scaling->ifScaled`
/// (`ediv` x / y, `edot` dSlackPos / dSlackNeg / ax / aty), then maps back
/// through `constraint_new_idx` / `constraint_type`; `ScaleVector(sense,
/// ...)` on the duals; sets `*value_valid` / `*dual_valid` from the flags.
#[allow(clippy::too_many_arguments)]
pub fn PDHG_PostSolve(
    pdhg: &mut CUPDLPwork,
    nCols_origin: cupdlp_int,
    constraint_new_idx: Option<&[cupdlp_int]>,
    constraint_type: Option<&[cupdlp_int]>,
    col_value: Option<&mut [cupdlp_float]>,
    col_dual: Option<&mut [cupdlp_float]>,
    row_value: Option<&mut [cupdlp_float]>,
    row_dual: Option<&mut [cupdlp_float]>,
    value_valid: &mut cupdlp_int,
    dual_valid: &mut cupdlp_int,
) -> cupdlp_retcode {
    let mut retcode: cupdlp_retcode = RETCODE_OK;

    let problem = &pdhg.problem;
    let iterates = &mut pdhg.iterates;
    let scaling = &pdhg.scaling;
    let resobj = &mut pdhg.resobj;
    let sense = problem.sense_origin;

    let nCols = problem.nCols;
    let nRows = problem.nRows;

    // flag
    let mut col_value_flag: cupdlp_int = 0;
    let mut col_dual_flag: cupdlp_int = 0;
    let mut row_value_flag: cupdlp_int = 0;
    let mut row_dual_flag: cupdlp_int = 0;

    // allocate buffer
    // CUPDLP_INIT_DOUBLE(col_buffer, problem->nCols);
    // CUPDLP_INIT_DOUBLE(row_buffer, problem->nRows);
    // CUPDLP_INIT_DOUBLE(col_buffer2, problem->nCols);
    let mut col_buffer = vec![0.0f64; nCols as usize];
    let mut row_buffer = vec![0.0f64; nRows as usize];
    let mut col_buffer2 = vec![0.0f64; nCols as usize];

    let iter = pdhg.timers.nIter;
    let k = (iter % 2) as usize;

    // x = iterates->x[iter % 2], and likewise y / ax / aty (split borrows).
    let [x0, x1] = &mut iterates.x;
    let [y0, y1] = &mut iterates.y;
    let [ax0, ax1] = &mut iterates.ax;
    let [aty0, aty1] = &mut iterates.aty;

    let x = if k == 0 { x0 } else { x1 };
    let y = if k == 0 { y0 } else { y1 };
    let ax = if k == 0 { ax0 } else { ax1 };
    let aty = if k == 0 { aty0 } else { aty1 };

    // unscale
    if scaling.ifScaled != 0 {
        cupdlp_ediv(&mut x.data, &pdhg.colScale, nCols);
        cupdlp_ediv(&mut y.data, &pdhg.rowScale, nRows);
        cupdlp_edot(&mut resobj.dSlackPos, &pdhg.colScale, nCols);
        cupdlp_edot(&mut resobj.dSlackNeg, &pdhg.colScale, nCols);
        cupdlp_edot(&mut ax.data, &pdhg.rowScale, nRows);
        cupdlp_edot(&mut aty.data, &pdhg.colScale, nCols);
    }

    // col value: extract x from (x, z)
    if let Some(cv) = col_value {
        // CUPDLP_COPY_VEC(col_value, x->data, cupdlp_float, nCols_origin);
        cupdlp_copy(cv, &x.data, nCols_origin);
        col_value_flag = 1;
    }

    // row value
    if let Some(rv) = row_value {
        if let Some(cni) = constraint_new_idx {
            // CUPDLP_COPY_VEC(row_buffer, ax->data, cupdlp_float, problem->nRows);
            cupdlp_copy(&mut row_buffer, &ax.data, nRows);

            // un-permute row value
            // for (int i = 0; i < problem->nRows; i++) {
            //   row_value[i] = row_buffer[constraint_new_idx[i]];
            // }
            for i in 0..nRows {
                rv[i as usize] = row_buffer[cni[i as usize] as usize];
            }
        } else {
            // CUPDLP_COPY_VEC(row_value, ax->data, cupdlp_float, problem->nRows);
            cupdlp_copy(rv, &ax.data, nRows);
        }

        if let Some(ct) = constraint_type {
            // CUPDLP_COPY_VEC(col_buffer, x->data, cupdlp_float, problem->nCols);
            cupdlp_copy(&mut col_buffer, &x.data, nCols);

            // EQ = 0, LEQ = 1, GEQ = 2, BOUND = 3
            // for (int i = 0, j = 0; i < problem->nRows; i++) {
            //   if (constraint_type[i] == 1) {  // LEQ: multiply -1
            //     row_value[i] = -row_value[i];
            //   } else if (constraint_type[i] == 3) {  // BOUND: get Ax from Ax - z
            //     row_value[i] = row_value[i] + col_buffer[nCols_origin + j];
            //     j++;
            //   }
            // }
            let mut j: cupdlp_int = 0;
            for i in 0..nRows {
                if ct[i as usize] == 1 {
                    // LEQ: multiply -1
                    rv[i as usize] = -rv[i as usize];
                } else if ct[i as usize] == 3 {
                    // BOUND: get Ax from Ax - z
                    rv[i as usize] = rv[i as usize] + col_buffer[(nCols_origin + j) as usize];
                    j += 1;
                }
            }
        }

        row_value_flag = 1;
    }

    // col duals of l <= x <= u
    if let Some(cd) = col_dual {
        // CUPDLP_COPY_VEC(col_buffer, resobj->dSlackPos, cupdlp_float, nCols_origin);
        cupdlp_copy(&mut col_buffer, &resobj.dSlackPos, nCols_origin);
        // CUPDLP_COPY_VEC(col_buffer2, resobj->dSlackNeg, cupdlp_float, nCols_origin);
        cupdlp_copy(&mut col_buffer2, &resobj.dSlackNeg, nCols_origin);

        // for (int i = 0; i < nCols_origin; i++) {
        //   col_dual[i] = col_buffer[i] - col_buffer2[i];
        // }
        for i in 0..nCols_origin {
            cd[i as usize] = col_buffer[i as usize] - col_buffer2[i as usize];
        }

        ScaleVector(sense, cd, nCols_origin);

        col_dual_flag = 1;
    }

    // row dual: recover y
    if let Some(rd) = row_dual {
        if let Some(cni) = constraint_new_idx {
            // CUPDLP_COPY_VEC(row_buffer, y->data, cupdlp_float, problem->nRows);
            cupdlp_copy(&mut row_buffer, &y.data, nRows);
            // un-permute row dual
            // for (int i = 0; i < problem->nRows; i++) {
            //   row_dual[i] = row_buffer[constraint_new_idx[i]];
            // }
            for i in 0..nRows {
                rd[i as usize] = row_buffer[cni[i as usize] as usize];
            }
        } else {
            // CUPDLP_COPY_VEC(row_dual, y->data, cupdlp_float, problem->nRows);
            cupdlp_copy(rd, &y.data, nRows);
        }

        ScaleVector(sense, rd, nRows);

        if let Some(ct) = constraint_type {
            // EQ = 0, LEQ = 1, GEQ = 2, BOUND = 3
            // for (int i = 0; i < problem->nRows; i++) {
            //   if (constraint_type[i] == 1) {  // LEQ: multiply -1
            //     row_dual[i] = -row_dual[i];
            //   }
            // }
            for i in 0..nRows {
                if ct[i as usize] == 1 {
                    // LEQ: multiply -1
                    rd[i as usize] = -rd[i as usize];
                }
            }
        }

        row_dual_flag = 1;
    }

    // valid
    // if (value_valid) { *value_valid = col_value_flag && row_value_flag; }
    // value_valid is always present in Rust (`&mut`), so the guard is dropped.
    *value_valid = if col_value_flag != 0 && row_value_flag != 0 {
        1
    } else {
        0
    };

    // if (dual_valid) { *dual_valid = col_dual_flag && row_dual_flag; }
    *dual_valid = if col_dual_flag != 0 && row_dual_flag != 0 {
        1
    } else {
        0
    };

    // exit_cleanup:
    // CUPDLP_FREE(col_buffer);
    // CUPDLP_FREE(row_buffer);
    // CUPDLP_FREE(col_buffer2);
    // Buffers are Vecs, dropped automatically.

    retcode
}

/// HiGHS: cupdlp_solver.c::LP_SolvePDHG (lines 1437-1498)
///
/// `CUPDLP_CALL(PDHG_SetUserParam)`, `PDHG_PrintHugeCUPDHG` (print-only,
/// dropped), `CUPDLP_CALL(PDHG_PreSolve)`, `has_variables = (*value_valid +
/// *dual_valid) != 0`, `CUPDLP_CALL(PDHG_Solve)`, `*model_status =
/// termCode; *num_iter = nIter`, `CUPDLP_CALL(PDHG_PostSolve)`. `fp` /
/// `fp_sol` select `writeJson` / `writeSol` (file writers, deferred:
/// HiGHS passes `fp = NULL`, `ifSaveSol = false`), so they are accepted and
/// ignored. The solution slices are reborrowed for the two calls that take
/// them.
#[allow(clippy::too_many_arguments)]
pub fn LP_SolvePDHG(
    pdhg: &mut CUPDLPwork,
    ifChangeIntParam: &[cupdlp_bool],
    intParam: &[cupdlp_int],
    ifChangeFloatParam: &[cupdlp_bool],
    floatParam: &[cupdlp_float],
    fp: Option<&str>,
    nCols_origin: cupdlp_int,
    mut col_value: Option<&mut [cupdlp_float]>,
    mut col_dual: Option<&mut [cupdlp_float]>,
    mut row_value: Option<&mut [cupdlp_float]>,
    mut row_dual: Option<&mut [cupdlp_float]>,
    value_valid: &mut cupdlp_int,
    dual_valid: &mut cupdlp_int,
    ifSaveSol: cupdlp_bool,
    fp_sol: Option<&str>,
    constraint_new_idx: Option<&[cupdlp_int]>,
    constraint_type: Option<&[cupdlp_int]>,
    model_status: &mut cupdlp_int,
    num_iter: &mut cupdlp_int,
) -> cupdlp_retcode {
    let mut retcode: cupdlp_retcode = RETCODE_OK;

    // Set the parameters first - which is silent
    // CUPDLP_CALL(PDHG_SetUserParam(pdhg, ifChangeIntParam, intParam,
    //                               ifChangeFloatParam, floatParam));
    let rc = PDHG_SetUserParam(
        pdhg,
        ifChangeIntParam,
        intParam,
        ifChangeFloatParam,
        floatParam,
    );
    if rc != RETCODE_OK {
        retcode = rc;
        // goto exit_cleanup;
        return retcode;
    }

    // Call PDHG_PrintHugeCUPDHG() if logging level (set in
    // PDHG_SetUserParam) is verbose
    if pdhg.settings.nLogLevel > 1 {
        // PDHG_PrintHugeCUPDHG(); // Print-only, dropped
    }

    // #if !defined(CUPDLP_CPU) ... #endif // GPU code, dropped

    // CUPDLP_CALL(PDHG_PreSolve(pdhg, nCols_origin, constraint_new_idx,
    //                           constraint_type, col_value, col_dual, row_value,
    //                           row_dual, value_valid, dual_valid));
    let rc = PDHG_PreSolve(
        pdhg,
        nCols_origin,
        constraint_new_idx,
        constraint_type,
        col_value.as_deref_mut(),
        col_dual.as_deref_mut(),
        row_value.as_deref_mut(),
        row_dual.as_deref_mut(),
        value_valid,
        dual_valid,
    );
    if rc != RETCODE_OK {
        retcode = rc;
        // goto exit_cleanup;
        return retcode;
    }

    let has_variables: cupdlp_int = if (*value_valid + *dual_valid) != 0 {
        1
    } else {
        0
    };

    // if (pdhg->settings->nLogLevel>0 && has_variables)
    //   cupdlp_printf("Hot starting with given column primal values and row dual values\n");
    // Print-only, dropped

    // CUPDLP_CALL(PDHG_Solve(&has_variables, pdhg));
    let rc = PDHG_Solve(&has_variables, pdhg);
    if rc != RETCODE_OK {
        retcode = rc;
        // goto exit_cleanup;
        return retcode;
    }

    *model_status = pdhg.resobj.termCode;
    *num_iter = pdhg.timers.nIter;

    // CUPDLP_CALL(PDHG_PostSolve(pdhg, nCols_origin, constraint_new_idx,
    //                            constraint_type, col_value, col_dual, row_value,
    //                            row_dual, value_valid, dual_valid));
    let rc = PDHG_PostSolve(
        pdhg,
        nCols_origin,
        constraint_new_idx,
        constraint_type,
        col_value,
        col_dual,
        row_value,
        row_dual,
        value_valid,
        dual_valid,
    );
    if rc != RETCODE_OK {
        retcode = rc;
        // goto exit_cleanup;
        return retcode;
    }

    // if (fp)
    //   writeJson(fp, pdhg);
    // File writer dropped (cupdlp_solver.c:1480-1481)

    // if (ifSaveSol && fp_sol) {
    //   if (strcmp(fp, fp_sol) != 0) {
    //     writeSol(fp_sol, nCols_origin, pdhg->problem->nRows, col_value, col_dual,
    //              row_value, row_dual);
    //   } else {
    //     if (pdhg->settings->nLogLevel>0)
    //       cupdlp_printf("Warning: fp and fp_sol are the same, stop saving solution.\n");
    //   }
    // }
    // File writer dropped (cupdlp_solver.c:1483-1491)

    // exit_cleanup:
    // Move this to CupdlpWrapper so pdhg->debug_pdlp_log_file_ can
    // still be used
    //
    //  PDHG_Destroy(&pdhg);

    retcode
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn presolve_is_a_no_op_without_valid_start() {
        let mut w = CUPDLPwork::default();
        let (mut vv, mut dv) = (0, 1);
        assert_eq!(
            PDHG_PreSolve(
                &mut w, 0, None, None, None, None, None, None, &mut vv, &mut dv
            ),
            RETCODE_OK
        );
    }

    /// Unscaled one-column one-row LEQ problem: the row value and dual are
    /// negated back (`constraint_type == 1`).
    #[test]
    fn postsolve_negates_leq_rows() {
        let mut w = CUPDLPwork::default();
        w.problem.nRows = 1;
        w.problem.nCols = 1;
        w.problem.sense_origin = 1.0;
        let v1 = |v: f64| CUPDLPvec {
            len: 1,
            data: vec![v],
        };
        w.iterates.x = [v1(3.0), v1(0.0)];
        w.iterates.y = [v1(2.0), v1(0.0)];
        w.iterates.ax = [v1(-3.0), v1(0.0)];
        w.iterates.aty = [v1(0.0), v1(0.0)];
        w.resobj.dSlackPos = vec![0.5];
        w.resobj.dSlackNeg = vec![0.25];
        let (mut cv, mut cd, mut rv, mut rd) = ([0.0], [0.0], [0.0], [0.0]);
        let (mut vv, mut dv) = (0, 0);
        let idx = [0];
        let ty = [1];
        let rc = PDHG_PostSolve(
            &mut w,
            1,
            Some(&idx),
            Some(&ty),
            Some(&mut cv),
            Some(&mut cd),
            Some(&mut rv),
            Some(&mut rd),
            &mut vv,
            &mut dv,
        );
        assert_eq!(rc, RETCODE_OK);
        assert_eq!((cv[0], cd[0], rv[0], rd[0]), (3.0, 0.25, 3.0, -2.0));
        assert_eq!((vv, dv), (1, 1));
    }
}
