//! HiGHS source: pdlp/cupdlp/cupdlp_utils.c (lines 1160-1525)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Matrix copies, format conversions (through [`super::super::cs`]),
//! `*_create` and `*_alloc` / `*_alloc_matrix`. The `void *src` +
//! `CUPDLP_MATRIX_FORMAT` pairs are [`CupdlpSrcMatrix`]; `CUPDLPx **` out
//! parameters are `&mut Option<Box<CUPDLPx>>`.
//!
//! Upstream quirk kept: `csr2csc` sizes its CSparse matrix from
//! `csc->nRows` / `csc->nMatElem` (the destination), not from the source
//! CSR. HiGHS never reaches `csr2csc` (`problem_alloc` is called with a
//! CSC source), but the port keeps the expression.

use super::super::cs::*;
use super::super::defs::*;

/// HiGHS: cupdlp_utils.c::dense_copy (lines 1160-1166)
pub fn dense_copy(dst: &mut CUPDLPdense, src: &CUPDLPdense) {
    dst.nRows = src.nRows;
    dst.nCols = src.nCols;
    cupdlp_copy(&mut dst.data, &src.data, src.nRows * src.nCols);
}

/// HiGHS: cupdlp_utils.c::csr_copy (lines 1168-1177)
pub fn csr_copy(dst: &mut CUPDLPcsr, src: &CUPDLPcsr) {
    dst.nRows = src.nRows;
    dst.nCols = src.nCols;
    dst.nMatElem = src.nMatElem;
    cupdlp_copy(&mut dst.rowMatBeg, &src.rowMatBeg, src.nRows + 1);
    cupdlp_copy(&mut dst.rowMatIdx, &src.rowMatIdx, src.nMatElem);
    cupdlp_copy(&mut dst.rowMatElem, &src.rowMatElem, src.nMatElem);
}

/// HiGHS: cupdlp_utils.c::csc_copy (lines 1179-1197)
pub fn csc_copy(dst: &mut CUPDLPcsc, src: &CUPDLPcsc) -> cupdlp_int {
    dst.nRows = src.nRows;
    dst.nCols = src.nCols;
    dst.nMatElem = src.nMatElem;
    cupdlp_copy(&mut dst.colMatBeg, &src.colMatBeg, src.nCols + 1);
    cupdlp_copy(&mut dst.colMatIdx, &src.colMatIdx, src.nMatElem);
    cupdlp_copy(&mut dst.colMatElem, &src.colMatElem, src.nMatElem);
    // GPU branch (lines 1188-1194) dropped: CUPDLP_CPU defined.
    0
}

/// HiGHS: cupdlp_utils.c::csr2csc (lines 1199-1220)
pub fn csr2csc(csc: &mut CUPDLPcsc, csr: &CUPDLPcsr) {
    let mut cs_csr = cupdlp_dcs_spalloc(csr.nCols, csc.nRows, csc.nMatElem, 1, 0)
        .expect("cupdlp_dcs_spalloc returned NULL");
    cupdlp_copy(&mut cs_csr.p, &csr.rowMatBeg, csr.nRows + 1);
    cupdlp_copy(&mut cs_csr.i, &csr.rowMatIdx, csr.nMatElem);
    cupdlp_copy(
        cs_csr.x.as_mut().expect("cupdlp_dcs values"),
        &csr.rowMatElem,
        csr.nMatElem,
    );

    let mut cs_csc = cupdlp_dcs_transpose(&cs_csr, 1).expect("cupdlp_dcs_transpose returned NULL");
    csc.nCols = cs_csc.m;
    csc.nRows = cs_csc.n;
    csc.nMatElem = cs_csc.nzmax;
    cupdlp_copy(&mut csc.colMatBeg, &cs_csc.p, cs_csc.n + 1);
    cupdlp_copy(&mut csc.colMatIdx, &cs_csc.i, cs_csc.nzmax);
    cupdlp_copy(
        &mut csc.colMatElem,
        cs_csc.x.as_ref().expect("cupdlp_dcs values"),
        cs_csc.nzmax,
    );

    // cupdlp_dcs_spfree(cs_csc); // Drop
    // cupdlp_dcs_spfree(cs_csr); // Drop
}

/// HiGHS: cupdlp_utils.c::csc2csr (lines 1222-1254)
pub fn csc2csr(csr: &mut CUPDLPcsr, csc: &CUPDLPcsc) -> cupdlp_int {
    // The transpose may need to be done on the GPU
    // Currently, it is done on the CPU

    let mut cs_csc = cupdlp_dcs_spalloc(csc.nRows, csc.nCols, csc.nMatElem, 1, 0)
        .expect("cupdlp_dcs_spalloc returned NULL");
    cupdlp_copy(&mut cs_csc.p, &csc.colMatBeg, csc.nCols + 1);
    cupdlp_copy(&mut cs_csc.i, &csc.colMatIdx, csc.nMatElem);
    cupdlp_copy(
        cs_csc.x.as_mut().expect("cupdlp_dcs values"),
        &csc.colMatElem,
        csc.nMatElem,
    );

    let mut cs_csr = cupdlp_dcs_transpose(&cs_csc, 1).expect("cupdlp_dcs_transpose returned NULL");
    csr.nCols = cs_csr.m;
    csr.nRows = cs_csr.n;
    csr.nMatElem = cs_csr.nzmax;
    cupdlp_copy(&mut csr.rowMatBeg, &cs_csr.p, cs_csr.n + 1);
    cupdlp_copy(&mut csr.rowMatIdx, &cs_csr.i, cs_csr.nzmax);
    cupdlp_copy(
        &mut csr.rowMatElem,
        cs_csr.x.as_ref().expect("cupdlp_dcs values"),
        cs_csr.nzmax,
    );

    // cupdlp_dcs_spfree(cs_csc); // Drop
    // cupdlp_dcs_spfree(cs_csr); // Drop

    // GPU branch (lines 1245-1251) dropped: CUPDLP_CPU defined.

    0
}

/// HiGHS: cupdlp_utils.c::dense2csr (lines 1256-1277)
pub fn dense2csr(csr: &mut CUPDLPcsr, dense: &CUPDLPdense) {
    csr.nRows = dense.nRows;
    csr.nCols = dense.nCols;

    let mut nnz: cupdlp_int = 0;
    let mut iCol: cupdlp_int;
    let mut iRow: cupdlp_int;
    csr.rowMatBeg[0] = 0;
    iRow = 0;
    while iRow < csr.nRows {
        iCol = 0;
        while iCol < csr.nCols {
            if dense.data[(iCol * csr.nRows + iRow) as usize] != 0.0 {
                csr.rowMatIdx[nnz as usize] = iCol;
                csr.rowMatElem[nnz as usize] = dense.data[(iCol * csr.nRows + iRow) as usize];
                nnz += 1;
            }
            iCol += 1;
        }
        csr.rowMatBeg[(iRow + 1) as usize] = nnz;
        iRow += 1;
    }
    csr.nMatElem = nnz;
}

/// HiGHS: cupdlp_utils.c::dense2csc (lines 1279-1301)
pub fn dense2csc(csc: &mut CUPDLPcsc, dense: &CUPDLPdense) {
    csc.nRows = dense.nRows;
    csc.nCols = dense.nCols;

    let mut nnz: cupdlp_int = 0;
    let mut iCol: cupdlp_int;
    let mut iRow: cupdlp_int;
    csc.colMatBeg[0] = 0;
    iCol = 0;
    while iCol < csc.nCols {
        iRow = 0;
        while iRow < csc.nRows {
            if dense.data[(iCol * csc.nRows + iRow) as usize] != 0.0 {
                csc.colMatIdx[nnz as usize] = iRow;
                csc.colMatElem[nnz as usize] = dense.data[(iCol * csc.nRows + iRow) as usize];
                nnz += 1;
            }
            iRow += 1;
        }
        csc.colMatBeg[(iCol + 1) as usize] = nnz;
        iCol += 1;
    }

    csc.nMatElem = nnz;
}

/// HiGHS: cupdlp_utils.c::csr2dense (lines 1303-1321)
pub fn csr2dense(dense: &mut CUPDLPdense, csr: &CUPDLPcsr) {
    dense.nRows = csr.nRows;
    dense.nCols = csr.nCols;

    let mut iRow: cupdlp_int;
    let mut iCol: cupdlp_int;
    let mut iMatElem: cupdlp_int = 0;
    iRow = 0;
    while iRow < dense.nRows {
        iCol = 0;
        while iCol < dense.nCols {
            if Some(iCol) == csr.rowMatIdx.get(iMatElem as usize).copied() {
                dense.data[(iRow * dense.nCols + iCol) as usize] =
                    csr.rowMatElem[iMatElem as usize];
                iMatElem += 1;
            } else {
                dense.data[(iRow * dense.nCols + iCol) as usize] = 0.0;
            }
            iCol += 1;
        }
        iRow += 1;
    }
}

/// HiGHS: cupdlp_utils.c::csc2dense (lines 1323-1341)
pub fn csc2dense(dense: &mut CUPDLPdense, csc: &CUPDLPcsc) {
    dense.nRows = csc.nRows;
    dense.nCols = csc.nCols;

    let mut iRow: cupdlp_int;
    let mut iCol: cupdlp_int;
    let mut iMatElem: cupdlp_int = 0;
    iCol = 0;
    while iCol < dense.nCols {
        iRow = 0;
        while iRow < dense.nRows {
            if Some(iRow) == csc.colMatIdx.get(iMatElem as usize).copied() {
                dense.data[(iRow * dense.nCols + iCol) as usize] =
                    csc.colMatElem[iMatElem as usize];
                iMatElem += 1;
            } else {
                dense.data[(iRow * dense.nCols + iCol) as usize] = 0.0;
            }
            iRow += 1;
        }
        iCol += 1;
    }
}

/// HiGHS: cupdlp_utils.c::dense_create (lines 1343-1349)
pub fn dense_create(dense: &mut Option<Box<CUPDLPdense>>) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    *dense = Some(Box::new(CUPDLPdense::default()));

    retcode
}

/// HiGHS: cupdlp_utils.c::csr_create (lines 1351-1357)
pub fn csr_create(csr: &mut Option<Box<CUPDLPcsr>>) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    *csr = Some(Box::new(CUPDLPcsr::default()));

    retcode
}

/// HiGHS: cupdlp_utils.c::csc_create (lines 1359-1366)
pub fn csc_create(csc: &mut Option<Box<CUPDLPcsc>>) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    *csc = Some(Box::new(CUPDLPcsc::default()));

    retcode
}

/// HiGHS: cupdlp_utils.c::dense_alloc_matrix (lines 1368-1389)
pub fn dense_alloc_matrix(
    dense: &mut CUPDLPdense,
    nRows: cupdlp_int,
    nCols: cupdlp_int,
    src: CupdlpSrcMatrix<'_>,
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    dense.data = vec![0.0; (nRows * nCols) as usize];

    match src {
        CupdlpSrcMatrix::Dense(s) => {
            dense_copy(dense, s);
        }
        CupdlpSrcMatrix::Csr(s) => {
            csr2dense(dense, s);
        }
        CupdlpSrcMatrix::Csc(s) => {
            csc2dense(dense, s);
        }
    }

    retcode
}

/// HiGHS: cupdlp_utils.c::csr_alloc_matrix (lines 1391-1429)
pub fn csr_alloc_matrix(
    csr: &mut CUPDLPcsr,
    nRows: cupdlp_int,
    nCols: cupdlp_int,
    src: CupdlpSrcMatrix<'_>,
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    // C switch on src_matrix_format; `default: break` cannot occur (enum is exhaustive).
    let nnz: cupdlp_int = match src {
        CupdlpSrcMatrix::Dense(_) => nRows * nCols,
        CupdlpSrcMatrix::Csr(s) => s.nMatElem,
        CupdlpSrcMatrix::Csc(s) => s.nMatElem,
    };
    // todo make sure this is right
    csr.rowMatBeg = vec![0; (nRows + 1) as usize];
    csr.rowMatIdx = vec![0; nnz as usize];
    csr.rowMatElem = vec![0.0; nnz as usize];

    match src {
        CupdlpSrcMatrix::Dense(s) => {
            dense2csr(csr, s);
        }
        CupdlpSrcMatrix::Csr(s) => {
            csr_copy(csr, s);
        }
        CupdlpSrcMatrix::Csc(s) => {
            csc2csr(csr, s);
        }
    }

    retcode
}

/// HiGHS: cupdlp_utils.c::csc_alloc_matrix (lines 1431-1469)
pub fn csc_alloc_matrix(
    csc: &mut CUPDLPcsc,
    nRows: cupdlp_int,
    nCols: cupdlp_int,
    src: CupdlpSrcMatrix<'_>,
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    // C switch on src_matrix_format; `default: break` cannot occur (enum is exhaustive).
    let nnz: cupdlp_int = match src {
        CupdlpSrcMatrix::Dense(_) => nRows * nCols,
        CupdlpSrcMatrix::Csr(s) => s.nMatElem,
        CupdlpSrcMatrix::Csc(s) => s.nMatElem,
    };

    csc.colMatBeg = vec![0; (nCols + 1) as usize];
    csc.colMatIdx = vec![0; nnz as usize];
    csc.colMatElem = vec![0.0; nnz as usize];

    match src {
        CupdlpSrcMatrix::Dense(s) => {
            dense2csc(csc, s);
        }
        CupdlpSrcMatrix::Csr(s) => {
            csr2csc(csc, s);
        }
        CupdlpSrcMatrix::Csc(s) => {
            csc_copy(csc, s);
        }
    }

    retcode
}

/// HiGHS: cupdlp_utils.c::dense_alloc (lines 1471-1482)
pub fn dense_alloc(
    dense: &mut CUPDLPdense,
    nRows: cupdlp_int,
    nCols: cupdlp_int,
    val: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    dense.nRows = nRows;
    dense.nCols = nCols;
    dense.data = vec![0.0; (nRows * nCols) as usize];

    cupdlp_copy(&mut dense.data, val, nRows * nCols);

    retcode
}

/// HiGHS: cupdlp_utils.c::csr_alloc (lines 1484-1504)
pub fn csr_alloc(
    csr: &mut CUPDLPcsr,
    nRows: cupdlp_int,
    nCols: cupdlp_int,
    nnz: cupdlp_int,
    row_ptr: &[cupdlp_int],
    col_ind: &[cupdlp_int],
    val: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    csr.nRows = nRows;
    csr.nCols = nCols;
    csr.nMatElem = nnz;

    csr.rowMatBeg = vec![0; (nRows + 1) as usize];
    csr.rowMatIdx = vec![0; nnz as usize];
    csr.rowMatElem = vec![0.0; nnz as usize];

    cupdlp_copy(&mut csr.rowMatBeg, row_ptr, nRows + 1);
    cupdlp_copy(&mut csr.rowMatIdx, col_ind, nnz);
    cupdlp_copy(&mut csr.rowMatElem, val, nnz);

    retcode
}

/// HiGHS: cupdlp_utils.c::csc_alloc (lines 1506-1525)
pub fn csc_alloc(
    csc: &mut CUPDLPcsc,
    nRows: cupdlp_int,
    nCols: cupdlp_int,
    nnz: cupdlp_int,
    col_ptr: &[cupdlp_int],
    row_ind: &[cupdlp_int],
    val: &[cupdlp_float],
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    csc.nRows = nRows;
    csc.nCols = nCols;
    csc.nMatElem = nnz;
    csc.colMatBeg = vec![0; (nCols + 1) as usize];
    csc.colMatIdx = vec![0; nnz as usize];
    csc.colMatElem = vec![0.0; nnz as usize];

    cupdlp_copy(&mut csc.colMatBeg, col_ptr, nCols + 1);
    cupdlp_copy(&mut csc.colMatIdx, row_ind, nnz);
    cupdlp_copy(&mut csc.colMatElem, val, nnz);

    retcode
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// 2x3 matrix [[1 0 2] [0 3 0]] in CSC.
    fn csc_2x3() -> CUPDLPcsc {
        CUPDLPcsc {
            nRows: 2,
            nCols: 3,
            nMatElem: 3,
            colMatBeg: vec![0, 1, 2, 3],
            colMatIdx: vec![0, 1, 0],
            colMatElem: vec![1.0, 3.0, 2.0],
            MatElemNormInf: 0.0,
        }
    }

    #[test]
    fn csc_to_csr_via_alloc_matrix() {
        let csc = csc_2x3();
        let mut csr = CUPDLPcsr::default();
        assert_eq!(
            csr_alloc_matrix(&mut csr, 2, 3, CupdlpSrcMatrix::Csc(&csc)),
            RETCODE_OK
        );
        assert_eq!((csr.nRows, csr.nCols, csr.nMatElem), (2, 3, 3));
        assert_eq!(csr.rowMatBeg, vec![0, 2, 3]);
        assert_eq!(csr.rowMatIdx, vec![0, 2, 1]);
        assert_eq!(csr.rowMatElem, vec![1.0, 2.0, 3.0]);
        let mut copy = CUPDLPcsc::default();
        assert_eq!(
            csc_alloc_matrix(&mut copy, 2, 3, CupdlpSrcMatrix::Csc(&csc)),
            RETCODE_OK
        );
        assert_eq!(copy.colMatBeg, csc.colMatBeg);
        assert_eq!(copy.colMatElem, csc.colMatElem);
    }

    #[test]
    fn dense_round_trip() {
        // dense2csc reads the dense buffer column-major (cupdlp_utils.c 1279-1301).
        // Column-major [[1 0 2] [0 3 0]] is [1, 0, 0, 3, 2, 0].
        let d = CUPDLPdense {
            nRows: 2,
            nCols: 3,
            data: vec![1.0, 0.0, 0.0, 3.0, 2.0, 0.0],
        };
        let mut csc = CUPDLPcsc::default();
        assert_eq!(
            csc_alloc_matrix(&mut csc, 2, 3, CupdlpSrcMatrix::Dense(&d)),
            RETCODE_OK
        );
        let expect = csc_2x3();
        assert_eq!(csc.colMatBeg, expect.colMatBeg);
        // csc_alloc_matrix allocates nRows * nCols entries (C: nnz = nRows * nCols);
        // only the first nMatElem are written.
        let nnz = csc.nMatElem as usize;
        assert_eq!(csc.colMatIdx[..nnz], expect.colMatIdx[..]);
        assert_eq!(csc.colMatElem[..nnz], expect.colMatElem[..]);
    }

    #[test]
    fn csc2dense_matches_c_sequential_scan() {
        // csc2dense (cupdlp_utils.c 1323-1341) walks iCol x iRow and consumes
        // colMatIdx[iMatElem] sequentially without checking column boundaries,
        // so on [[1 0 2] [0 3 0]] it writes row-major [1, 2, 0, 3, 0, 0].
        // The skeleton expected the mathematically correct layout; the C
        // (ported faithfully) does not produce it. See commit message.
        let csc = csc_2x3();
        let mut d = CUPDLPdense::default();
        assert_eq!(
            dense_alloc_matrix(&mut d, 2, 3, CupdlpSrcMatrix::Csc(&csc)),
            RETCODE_OK
        );
        assert_eq!(d.data, vec![1.0, 2.0, 0.0, 3.0, 0.0, 0.0]);
    }
}
