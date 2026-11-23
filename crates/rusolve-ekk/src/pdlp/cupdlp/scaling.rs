//! HiGHS source: pdlp/cupdlp/cupdlp_scaling.h, pdlp/cupdlp/cupdlp_scaling.c
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Ruiz, L2 and Pock-Chambolle scaling, applied in that order by
//! `PDHG_Scale_Data`. The `#if CUPDLP_DEBUG` blocks (statistics printed
//! before and after rescaling) are not compiled (`CUPDLP_DEBUG`
//! undefined); the `if (log_level) cupdlp_printf(...)` lines print only and
//! are dropped, the computations next to them stay.

use super::defs::*;
use super::linalg::*;

/// HiGHS: cupdlp_scaling.c::scale_problem (lines 17-45) (`static`)
pub fn scale_problem(
    csc: &mut CUPDLPcsc,
    cost: &mut [cupdlp_float],
    lower: &mut [cupdlp_float],
    upper: &mut [cupdlp_float],
    rhs: &mut [cupdlp_float],
    col_scaling: &[cupdlp_float],
    row_scaling: &[cupdlp_float],
) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;
    let nRows = csc.nRows;
    let nCols = csc.nCols;

    cupdlp_cdiv(cost, col_scaling, nCols);
    cupdlp_cdot(lower, col_scaling, nCols);
    cupdlp_cdot(upper, col_scaling, nCols);
    cupdlp_cdiv(rhs, row_scaling, nRows);

    // row scaling
    for i in 0..csc.colMatBeg[nCols as usize] {
        csc.colMatElem[i as usize] /= row_scaling[csc.colMatIdx[i as usize] as usize];
    }
    // col scaling
    for i in 0..nCols {
        for j in csc.colMatBeg[i as usize]..csc.colMatBeg[(i + 1) as usize] {
            csc.colMatElem[j as usize] /= col_scaling[i as usize];
        }
    }

    // w->scaling->ifScaled = 1;

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_scaling.c::cupdlp_ruiz_scaling (lines 47-120)
pub fn cupdlp_ruiz_scaling(
    csc: &mut CUPDLPcsc,
    cost: &mut [cupdlp_float],
    lower: &mut [cupdlp_float],
    upper: &mut [cupdlp_float],
    rhs: &mut [cupdlp_float],
    scaling: &mut CUPDLPscaling,
) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;

    let nRows = csc.nRows;
    let nCols = csc.nCols;

    // CUPDLP_INIT_ZERO_DOUBLE(current_col_scaling, nCols);
    let mut current_col_scaling = vec![0.0; nCols as usize];
    // CUPDLP_INIT_ZERO_DOUBLE(current_row_scaling, nRows);
    let mut current_row_scaling = vec![0.0; nRows as usize];

    for i in 0..scaling.RuizTimes {
        cupdlp_zero(&mut current_col_scaling, nCols);
        cupdlp_zero(&mut current_row_scaling, nRows);

        // if (csc != NULL) is always true in Rust
        {
            for j in 0..nCols {
                if csc.colMatBeg[j as usize] == csc.colMatBeg[(j + 1) as usize] {
                    current_col_scaling[j as usize] = 0.0;
                } else {
                    let start = csc.colMatBeg[j as usize] as usize;
                    let len = csc.colMatBeg[(j + 1) as usize] - csc.colMatBeg[j as usize];
                    current_col_scaling[j as usize] =
                        GenNorm(&csc.colMatElem[start..], len, scaling.RuizNorm).sqrt();
                }
            }
        }
        for j in 0..nCols {
            if current_col_scaling[j as usize] == 0.0 {
                current_col_scaling[j as usize] = 1.0;
            }
        }

        if scaling.RuizNorm == f64::INFINITY {
            if nRows > 0 {
                // inf norm of rows of csc
                for j in 0..csc.colMatBeg[nCols as usize] {
                    if current_row_scaling[csc.colMatIdx[j as usize] as usize]
                        < ABS(csc.colMatElem[j as usize])
                    {
                        current_row_scaling[csc.colMatIdx[j as usize] as usize] =
                            ABS(csc.colMatElem[j as usize]);
                    }
                }
                for j in 0..nRows {
                    if current_row_scaling[j as usize] == 0.0 {
                        current_row_scaling[j as usize] = 1.0;
                    } else {
                        current_row_scaling[j as usize] = current_row_scaling[j as usize].sqrt();
                    }
                }
            }
        } else {
            // C lines 103-104: cupdlp_printf dropped, exit(1) kept.
            std::process::exit(1);
        }

        // apply scaling
        // scale_problem(work, current_col_scaling, current_row_scaling);
        scale_problem(
            csc,
            cost,
            lower,
            upper,
            rhs,
            &current_col_scaling,
            &current_row_scaling,
        );

        // update scaling
        cupdlp_cdot(&mut scaling.colScale, &current_col_scaling, nCols);
        cupdlp_cdot(&mut scaling.rowScale, &current_row_scaling, nRows);
    }

    // exit_cleanup:
    // cupdlp_free(current_col_scaling);
    // cupdlp_free(current_row_scaling);
    retcode
}

/// HiGHS: cupdlp_scaling.c::cupdlp_l2norm_scaling (lines 122-172)
pub fn cupdlp_l2norm_scaling(
    csc: &mut CUPDLPcsc,
    cost: &mut [cupdlp_float],
    lower: &mut [cupdlp_float],
    upper: &mut [cupdlp_float],
    rhs: &mut [cupdlp_float],
    scaling: &mut CUPDLPscaling,
) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;
    let nRows = csc.nRows;
    let nCols = csc.nCols;

    // CUPDLP_INIT_ZERO_DOUBLE(current_col_scaling, nCols);
    let mut current_col_scaling = vec![0.0; nCols as usize];
    // CUPDLP_INIT_ZERO_DOUBLE(current_row_scaling, nRows);
    let mut current_row_scaling = vec![0.0; nRows as usize];

    if nRows > 0 {
        for j in 0..nCols {
            if csc.colMatBeg[j as usize] == csc.colMatBeg[(j + 1) as usize] {
                current_col_scaling[j as usize] = 1.0;
            } else {
                let start = csc.colMatBeg[j as usize] as usize;
                let len = csc.colMatBeg[(j + 1) as usize] - csc.colMatBeg[j as usize];
                current_col_scaling[j as usize] =
                    GenNorm(&csc.colMatElem[start..], len, 2.0).sqrt();
            }
        }

        for i in 0..csc.colMatBeg[nCols as usize] {
            current_row_scaling[csc.colMatIdx[i as usize] as usize] +=
                csc.colMatElem[i as usize].powf(2.0);
        }
        for i in 0..nRows {
            current_row_scaling[i as usize] = current_row_scaling[i as usize].sqrt().sqrt();
            if current_row_scaling[i as usize] == 0.0 {
                current_row_scaling[i as usize] = 1.0;
            }
        }
    }
    // apply scaling
    // scale_problem(work, current_col_scaling, current_row_scaling);
    scale_problem(
        csc,
        cost,
        lower,
        upper,
        rhs,
        &current_col_scaling,
        &current_row_scaling,
    );

    // update scaling
    cupdlp_cdot(&mut scaling.colScale, &current_col_scaling, nCols);
    cupdlp_cdot(&mut scaling.rowScale, &current_row_scaling, nRows);

    // exit_cleanup:
    // cupdlp_free(current_col_scaling);
    // cupdlp_free(current_row_scaling);
    retcode
}

/// HiGHS: cupdlp_scaling.c::cupdlp_pc_scaling (lines 174-231)
pub fn cupdlp_pc_scaling(
    csc: &mut CUPDLPcsc,
    cost: &mut [cupdlp_float],
    lower: &mut [cupdlp_float],
    upper: &mut [cupdlp_float],
    rhs: &mut [cupdlp_float],
    scaling: &mut CUPDLPscaling,
) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;
    let nRows = csc.nRows;
    let nCols = csc.nCols;
    let alpha = scaling.PcAlpha;

    // CUPDLP_INIT_ZERO_DOUBLE(current_col_scaling, nCols);
    let mut current_col_scaling = vec![0.0; nCols as usize];
    // CUPDLP_INIT_ZERO_DOUBLE(current_row_scaling, nRows);
    let mut current_row_scaling = vec![0.0; nRows as usize];

    if alpha > 2.0 || alpha < 0.0 {
        // C lines 190-191: cupdlp_printf dropped, exit(1) kept.
        std::process::exit(1);
    }

    if nRows > 0 {
        for i in 0..nCols {
            for j in csc.colMatBeg[i as usize]..csc.colMatBeg[(i + 1) as usize] {
                current_col_scaling[i as usize] += ABS(csc.colMatElem[j as usize]).powf(alpha);
            }
            current_col_scaling[i as usize] =
                current_col_scaling[i as usize].powf(1.0 / alpha).sqrt();
            if current_col_scaling[i as usize] == 0.0 {
                current_col_scaling[i as usize] = 1.0;
            }
        }

        for i in 0..csc.colMatBeg[nCols as usize] {
            current_row_scaling[csc.colMatIdx[i as usize] as usize] +=
                ABS(csc.colMatElem[i as usize]).powf(2.0 - alpha);
        }
        for i in 0..nRows {
            current_row_scaling[i as usize] = current_row_scaling[i as usize]
                .powf(1.0 / (2.0 - alpha))
                .sqrt();
            if current_row_scaling[i as usize] == 0.0 {
                current_row_scaling[i as usize] = 1.0;
            }
        }
    }

    // apply scaling
    // scale_problem(work, current_col_scaling, current_row_scaling);
    scale_problem(
        csc,
        cost,
        lower,
        upper,
        rhs,
        &current_col_scaling,
        &current_row_scaling,
    );

    // update scaling
    cupdlp_cdot(&mut scaling.colScale, &current_col_scaling, nCols);
    cupdlp_cdot(&mut scaling.rowScale, &current_row_scaling, nRows);

    // exit_cleanup:
    // cupdlp_free(current_col_scaling);
    // cupdlp_free(current_row_scaling);
    retcode
}

/// HiGHS: cupdlp_scaling.c::PDHG_Scale_Data (lines 233-393)
pub fn PDHG_Scale_Data(
    log_level: cupdlp_int,
    csc: &mut CUPDLPcsc,
    ifScaling: cupdlp_int,
    scaling: &mut CUPDLPscaling,
    cost: &mut [cupdlp_float],
    lower: &mut [cupdlp_float],
    upper: &mut [cupdlp_float],
    rhs: &mut [cupdlp_float],
) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;
    // scaling->dObjScale = 1.0;

    // #if CUPDLP_DEBUG block (lines 240-304) is not compiled.

    'exit_cleanup: {
        if ifScaling != 0 {
            // C lines 307-308: log_level prints dropped.

            if scaling.ifRuizScaling != 0 {
                // C line 311: log_level print dropped.
                // CUPDLP_CALL(cupdlp_ruiz_scaling(...)) (line 312)
                retcode = cupdlp_ruiz_scaling(csc, cost, lower, upper, rhs, scaling);
                if retcode != RETCODE_OK {
                    retcode = RETCODE_FAILED;
                    break 'exit_cleanup;
                }
                scaling.ifScaled = 1;
            }
            if scaling.ifL2Scaling != 0 {
                // C line 316: log_level print dropped.
                // CUPDLP_CALL(cupdlp_l2norm_scaling(...)) (line 317)
                retcode = cupdlp_l2norm_scaling(csc, cost, lower, upper, rhs, scaling);
                if retcode != RETCODE_OK {
                    retcode = RETCODE_FAILED;
                    break 'exit_cleanup;
                }
                scaling.ifScaled = 1;
            }
            if scaling.ifPcScaling != 0 {
                // C line 321: log_level print dropped.
                // CUPDLP_CALL(cupdlp_pc_scaling(...)) (line 322)
                retcode = cupdlp_pc_scaling(csc, cost, lower, upper, rhs, scaling);
                if retcode != RETCODE_OK {
                    retcode = RETCODE_FAILED;
                    break 'exit_cleanup;
                }
                scaling.ifScaled = 1;
            }

            // C line 326: log_level print dropped.
        }

        // C line 330: commented-out csc2csr dropped.
        // #if CUPDLP_DEBUG block (lines 331-388) is not compiled.
    }

    // exit_cleanup:

    retcode
}

/// HiGHS: cupdlp_scaling.c::Init_Scaling (lines 395-425)
pub fn Init_Scaling(
    log_level: cupdlp_int,
    scaling: &mut CUPDLPscaling,
    ncols: cupdlp_int,
    nrows: cupdlp_int,
    cost: &[cupdlp_float],
    rhs: &[cupdlp_float],
) -> cupdlp_retcode {
    let mut retcode = RETCODE_OK;

    // C lines 400-401, 411-415: NULL, CUPDLP_INIT_DOUBLE, then every entry set to 1.0.
    scaling.rowScale = vec![1.0; nrows as usize];
    scaling.colScale = vec![1.0; ncols as usize];

    scaling.ifRuizScaling = 1;
    scaling.ifL2Scaling = 0;
    scaling.ifPcScaling = 1;

    // todo, read these paras
    scaling.RuizTimes = 10;
    scaling.RuizNorm = f64::INFINITY;
    scaling.PcAlpha = 1.0;

    scaling.dNormCost = twoNorm(cost, ncols);
    scaling.dNormRhs = twoNorm(rhs, nrows);
    // C lines 419-421: log_level print dropped.

    // exit_cleanup:
    retcode
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn init_scaling_defaults() {
        let mut s = CUPDLPscaling::default();
        assert_eq!(
            Init_Scaling(0, &mut s, 2, 1, &[3.0, 4.0], &[2.0]),
            RETCODE_OK
        );
        assert_eq!((s.ifRuizScaling, s.ifL2Scaling, s.ifPcScaling), (1, 0, 1));
        assert_eq!(s.RuizTimes, 10);
        assert_eq!(s.RuizNorm, f64::INFINITY);
        assert_eq!(s.PcAlpha, 1.0);
        assert_eq!(s.colScale, vec![1.0, 1.0]);
        assert_eq!(s.rowScale, vec![1.0]);
        assert_eq!(s.dNormCost, 5.0);
        assert_eq!(s.dNormRhs, 2.0);
    }

    #[test]
    fn scale_problem_divides_by_row_then_col() {
        let mut csc = CUPDLPcsc {
            nRows: 1,
            nCols: 1,
            nMatElem: 1,
            colMatBeg: vec![0, 1],
            colMatIdx: vec![0],
            colMatElem: vec![8.0],
            MatElemNormInf: 0.0,
        };
        let (mut c, mut l, mut u, mut r) = (vec![6.0], vec![1.0], vec![3.0], vec![10.0]);
        assert_eq!(
            scale_problem(&mut csc, &mut c, &mut l, &mut u, &mut r, &[2.0], &[4.0]),
            RETCODE_OK
        );
        assert_eq!(csc.colMatElem, vec![1.0]);
        assert_eq!((c[0], l[0], u[0], r[0]), (3.0, 2.0, 6.0, 2.5));
    }

    fn csc_from(
        nRows: cupdlp_int,
        nCols: cupdlp_int,
        beg: Vec<cupdlp_int>,
        idx: Vec<cupdlp_int>,
        elem: Vec<cupdlp_float>,
    ) -> CUPDLPcsc {
        CUPDLPcsc {
            nRows,
            nCols,
            nMatElem: elem.len() as cupdlp_int,
            colMatBeg: beg,
            colMatIdx: idx,
            colMatElem: elem,
            MatElemNormInf: 0.0,
        }
    }

    fn ones_scaling(ncols: usize, nrows: usize) -> CUPDLPscaling {
        CUPDLPscaling {
            colScale: vec![1.0; ncols],
            rowScale: vec![1.0; nrows],
            RuizTimes: 10,
            RuizNorm: f64::INFINITY,
            PcAlpha: 1.0,
            ..Default::default()
        }
    }

    // C 64-115: iteration 1 gives col sqrt(4) = 2 and row sqrt(4) = 2 (row inf norm read
    // before the column division), so A becomes 1; later iterations scale by 1. The empty
    // column (C 70-71, 79-83) and the empty row (C 95-96) get 1.
    #[test]
    fn ruiz_scaling_with_empty_col_and_row() {
        let mut csc = csc_from(2, 2, vec![0, 1, 1], vec![0], vec![4.0]);
        let mut s = ones_scaling(2, 2);
        let (mut c, mut l, mut u, mut r) = (
            vec![6.0, 5.0],
            vec![1.0, 1.0],
            vec![3.0, 3.0],
            vec![10.0, 7.0],
        );
        assert_eq!(
            cupdlp_ruiz_scaling(&mut csc, &mut c, &mut l, &mut u, &mut r, &mut s),
            RETCODE_OK
        );
        assert_eq!(csc.colMatElem, vec![1.0]);
        assert_eq!(s.colScale, vec![2.0, 1.0]);
        assert_eq!(s.rowScale, vec![2.0, 1.0]);
        assert_eq!(
            (c, l, u, r),
            (
                vec![3.0, 5.0],
                vec![2.0, 1.0],
                vec![6.0, 3.0],
                vec![5.0, 7.0]
            )
        );
    }

    // C 138-157: col sqrt(twoNorm([4])) = 2, row sqrt(sqrt(pow(4, 2))) = 2.
    #[test]
    fn l2norm_scaling_one_entry() {
        let mut csc = csc_from(1, 1, vec![0, 1], vec![0], vec![4.0]);
        let mut s = ones_scaling(1, 1);
        let (mut c, mut l, mut u, mut r) = (vec![6.0], vec![1.0], vec![3.0], vec![10.0]);
        assert_eq!(
            cupdlp_l2norm_scaling(&mut csc, &mut c, &mut l, &mut u, &mut r, &mut s),
            RETCODE_OK
        );
        assert_eq!(csc.colMatElem, vec![1.0]);
        assert_eq!((s.colScale[0], s.rowScale[0]), (2.0, 2.0));
        assert_eq!((c[0], l[0], u[0], r[0]), (3.0, 2.0, 6.0, 5.0));
    }

    // C 138: with nRows == 0 the zero-initialised column scaling is applied as is, so the
    // cost is divided by 0 and colScale multiplied by 0. Same for PC scaling (C 194).
    #[test]
    fn l2_and_pc_scaling_without_rows_apply_zero_col_scaling() {
        for pc in [false, true] {
            let mut csc = csc_from(0, 1, vec![0, 0], vec![], vec![]);
            let mut s = ones_scaling(1, 0);
            let (mut c, mut l, mut u, mut r) = (vec![2.0], vec![1.0], vec![3.0], vec![]);
            let ret = if pc {
                cupdlp_pc_scaling(&mut csc, &mut c, &mut l, &mut u, &mut r, &mut s)
            } else {
                cupdlp_l2norm_scaling(&mut csc, &mut c, &mut l, &mut u, &mut r, &mut s)
            };
            assert_eq!(ret, RETCODE_OK);
            assert_eq!(c[0], f64::INFINITY);
            assert_eq!((l[0], u[0], s.colScale[0]), (0.0, 0.0, 0.0));
        }
    }

    // C 194-215 with alpha = 1 on one column [4 (row 0), -12 (row 1)]: col
    // sqrt(pow(|4| + |-12|, 1)) = 4; row 0 sqrt(pow(4, 1)) = 2, row 1 sqrt(pow(12, 1)),
    // and the empty row 2 gives sqrt(0) = 0, reset to 1 (C 212-213).
    #[test]
    fn pc_scaling_alpha_one() {
        let mut csc = csc_from(3, 1, vec![0, 2], vec![0, 1], vec![4.0, -12.0]);
        let mut s = ones_scaling(1, 3);
        let (mut c, mut l, mut u, mut r) = (vec![8.0], vec![1.0], vec![3.0], vec![10.0, 9.0, 7.0]);
        assert_eq!(
            cupdlp_pc_scaling(&mut csc, &mut c, &mut l, &mut u, &mut r, &mut s),
            RETCODE_OK
        );
        assert_eq!(s.colScale[0], 4.0);
        assert_eq!(s.rowScale[0], 2.0);
        assert_eq!(s.rowScale[1], 12.0_f64.powf(1.0).sqrt());
        assert_eq!(s.rowScale[2], 1.0);
        assert_eq!(csc.colMatElem[0], 4.0 / 2.0 / 4.0);
        assert_eq!(csc.colMatElem[1], -12.0 / s.rowScale[1] / 4.0);
        assert_eq!((c[0], l[0], u[0]), (2.0, 4.0, 12.0));
        assert_eq!(r, vec![5.0, 9.0 / s.rowScale[1], 7.0]);
    }

    // C 306-327: Ruiz then PC (Init_Scaling defaults); ifScaling == 0 leaves everything.
    #[test]
    fn pdhg_scale_data_ruiz_then_pc() {
        let mut csc = csc_from(1, 1, vec![0, 1], vec![0], vec![4.0]);
        let mut s = CUPDLPscaling::default();
        Init_Scaling(0, &mut s, 1, 1, &[6.0], &[10.0]);
        let (mut c, mut l, mut u, mut r) = (vec![6.0], vec![1.0], vec![3.0], vec![10.0]);
        assert_eq!(
            PDHG_Scale_Data(0, &mut csc, 0, &mut s, &mut c, &mut l, &mut u, &mut r),
            RETCODE_OK
        );
        assert_eq!((s.ifScaled, csc.colMatElem[0], c[0]), (0, 4.0, 6.0));
        assert_eq!(
            PDHG_Scale_Data(0, &mut csc, 1, &mut s, &mut c, &mut l, &mut u, &mut r),
            RETCODE_OK
        );
        assert_eq!(s.ifScaled, 1);
        assert_eq!(csc.colMatElem, vec![1.0]);
        assert_eq!((s.colScale[0], s.rowScale[0]), (2.0, 2.0));
        assert_eq!((c[0], l[0], u[0], r[0]), (3.0, 2.0, 6.0, 5.0));
    }
}
