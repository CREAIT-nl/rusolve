//! HiGHS source: ipm/hipo/ipm/KktMatrix.h, ipm/hipo/ipm/KktMatrix.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Augmented system (AS) and normal equations (NE) matrices. The four
//! references the C++ caches (`model`, `regul`, `info`, `logger`) are passed
//! to the methods that use them. `NE_nz_limit` is the C++ `std::atomic<Int64>`
//! (`AtomicI64`, relaxed): `chooseNla`'s AS task stores it while the NE task
//! reads it (FactorHiGHSSolver.cpp 266, KktMatrix.cpp 171-172).
//! The two `chooseNla` tasks share `kkt_` in C++ and touch disjoint members;
//! [`KktMatrix::split_AS_NE`] borrows those members apart ([`KktAsPart`],
//! [`KktNePart`]) so each task holds its own `&mut` (borrowed apart).

use super::info::Info;
use super::model::Model;
use super::parameters::Regularisation;
use super::status::{kStatusOk, kStatusOverflow};
use crate::ipm::hipo::auxiliary::Clock;
use crate::ipm::hipo::auxiliary::freeVector;
use crate::ipm::hipo::auxiliary::logger::Logger;
use crate::ipm::hipo::{Int, Int64};
use crate::lp_data::hconst::kHighsIInf;
use crate::util::highs_utils::cpp_min;
use std::sync::atomic::{AtomicI64, Ordering};

/// HiGHS: KktMatrix.h `struct KktMatrix`
#[derive(Debug)]
pub struct KktMatrix {
    pub ptrNE: Vec<Int>,
    pub rowsNE: Vec<Int>,
    pub valNE: Vec<f64>,
    pub ptrA_rw: Vec<Int>,
    pub idxA_rw: Vec<Int>,
    pub corr_A: Vec<Int>,
    pub NE_nz_limit: AtomicI64,
    pub ptrAS: Vec<Int>,
    pub rowsAS: Vec<Int>,
    pub valAS: Vec<f64>,
    pub iperm: Vec<Int>,
}

/// Ours (no C++ counterpart): the members of [`KktMatrix`] that the AS task of
/// `chooseNla` touches (`buildASstructure`, `analyseAS`).
pub struct KktAsPart<'a> {
    pub ptrAS: &'a mut Vec<Int>,
    pub rowsAS: &'a mut Vec<Int>,
    pub valAS: &'a mut Vec<f64>,
    pub NE_nz_limit: &'a AtomicI64,
}

/// Ours (no C++ counterpart): the members of [`KktMatrix`] that the NE task of
/// `chooseNla` touches (`buildNEstructure`).
pub struct KktNePart<'a> {
    pub ptrNE: &'a mut Vec<Int>,
    pub rowsNE: &'a mut Vec<Int>,
    pub ptrA_rw: &'a mut Vec<Int>,
    pub idxA_rw: &'a mut Vec<Int>,
    pub corr_A: &'a mut Vec<Int>,
    pub NE_nz_limit: &'a AtomicI64,
}

/// Largest number of dense columns the Schur-complement treatment takes on.
pub const kDenseMax: usize = 128;

/// rusolve `ipm_hipo_dense_cols` (Andersen 1996): columns of A with
/// more nonzeros than this are left out of A D A^T and recovered by preconditioned CG in
/// `solveNE`. `Int::MAX` when the lever is off or the model has a Hessian.
pub fn dense_col_threshold(A: &crate::util::sparse_matrix::HighsSparseMatrix, qp: bool) -> Int {
    if qp || !crate::ipm::hipo::blas::DENSE_COLS.load(Ordering::Relaxed) || A.num_col_ <= 0 {
        return Int::MAX;
    }
    let avg = A.num_nz() as f64 / A.num_col_ as f64;
    let thr = cpp_min(Int::MAX as f64, (10.0 * avg).max(40.0)) as Int;
    // Woodbury costs one solve per dense column: with more than kDenseMax of them the columns
    // are not "dense" in the sense of Andersen 1996 and the lever stays inert
    if crate::ipm::hipo::blas::DENSE_SCHUR.load(Ordering::Relaxed) {
        let k = (0..A.num_col_ as usize)
            .filter(|&j| A.start_[j + 1] - A.start_[j] > thr)
            .count();
        if k > kDenseMax {
            return Int::MAX;
        }
    }
    thr
}

impl KktMatrix {
    /// Ours (no C++ counterpart): borrow the AS and NE members apart; both parts
    /// share the atomic `NE_nz_limit`, as the C++ tasks do.
    pub fn split_AS_NE(&mut self) -> (KktAsPart<'_>, KktNePart<'_>) {
        (
            KktAsPart {
                ptrAS: &mut self.ptrAS,
                rowsAS: &mut self.rowsAS,
                valAS: &mut self.valAS,
                NE_nz_limit: &self.NE_nz_limit,
            },
            KktNePart {
                ptrNE: &mut self.ptrNE,
                rowsNE: &mut self.rowsNE,
                ptrA_rw: &mut self.ptrA_rw,
                idxA_rw: &mut self.idxA_rw,
                corr_A: &mut self.corr_A,
                NE_nz_limit: &self.NE_nz_limit,
            },
        )
    }

    /// HiGHS: KktMatrix.cpp::KktMatrix::KktMatrix (lines 7-9)
    pub fn new() -> Self {
        KktMatrix {
            ptrNE: Vec::new(),
            rowsNE: Vec::new(),
            valNE: Vec::new(),
            ptrA_rw: Vec::new(),
            idxA_rw: Vec::new(),
            corr_A: Vec::new(),
            NE_nz_limit: AtomicI64::new(kHighsIInf as Int64),
            ptrAS: Vec::new(),
            rowsAS: Vec::new(),
            valAS: Vec::new(),
            iperm: Vec::new(),
        }
    }

    /// HiGHS: KktMatrix.cpp::KktMatrix::buildASstructure (lines 11-71); the body is
    /// [`KktAsPart::buildASstructure`].
    pub fn buildASstructure(&mut self, model: &Model, info: &mut Info, logger: &Logger) -> Int {
        let (mut kkt_AS, _) = self.split_AS_NE();
        kkt_AS.buildASstructure(model, &mut info.AS_structure_time, logger)
    }

    /// HiGHS: KktMatrix.cpp::KktMatrix::buildASvalues (lines 73-90)
    pub fn buildASvalues(&mut self, model: &Model, info: &mut Info, scaling: &[f64]) -> Int {
        // build AS values that change during iterations.

        debug_assert!(!self.ptrAS.is_empty() && !self.rowsAS.is_empty());
        let clock = Clock::new();

        let n = model.A().num_col_;

        for i in 0..n {
            self.valAS[self.ptrAS[i as usize] as usize] = if scaling.is_empty() {
                -1.0
            } else {
                -scaling[i as usize]
            };
            if model.qp() {
                self.valAS[self.ptrAS[i as usize] as usize] -= model.sense() * model.Q().diag(i);
            }
        }

        info.matrix_time += clock.stop();

        kStatusOk
    }

    /// HiGHS: KktMatrix.cpp::KktMatrix::buildNEstructure (lines 92-189); the body is
    /// [`KktNePart::buildNEstructure`].
    pub fn buildNEstructure(&mut self, model: &Model, info: &mut Info, logger: &Logger) -> Int {
        let (_, mut kkt_NE) = self.split_AS_NE();
        kkt_NE.buildNEstructure(model, &mut info.NE_structure_time, logger)
    }

    /// HiGHS: KktMatrix.cpp::KktMatrix::buildNEvalues (lines 191-246)
    pub fn buildNEvalues(
        &mut self,
        model: &Model,
        regul: &Regularisation,
        info: &mut Info,
        scaling: &[f64],
    ) -> Int {
        // given the NE structure already computed, fill in the NE values

        debug_assert!(!self.ptrNE.is_empty() && !self.rowsNE.is_empty());
        let clock = Clock::new();

        let A = model.A();
        let Q = model.Q();
        let m = A.num_row_;

        self.valNE.resize(self.rowsNE.len(), 0.0);

        // rusolve lever `ipm_hipo_par_ne`: every column of A D A^T
        // is computed from scratch by one row of A with its own scatter array, so blocks of
        // rows can run on different threads and each entry keeps its summation order
        // (bit-identical to the serial loop).
        let threads = crate::parallel::highs_parallel::num_threads();
        let par = crate::ipm::hipo::blas::PAR_NE.load(Ordering::Relaxed)
            && threads > 1
            && self.rowsNE.len() >= 100_000;
        let nchunks: Int = if par {
            cpp_min(m, 8 * threads as Int)
        } else {
            1
        };
        let ptrA_rw = &self.ptrA_rw;
        let idxA_rw = &self.idxA_rw;
        let corr_A = &self.corr_A;
        let ptrNE = &self.ptrNE;
        let rowsNE = &self.rowsNE;
        let (a_val, a_start, a_index) = (&A.value_, &A.start_, &A.index_);
        let qp = model.qp();
        let sense = model.sense();
        let qdiag: Vec<f64> = if qp {
            (0..A.num_col_).map(|c| Q.diag(c)).collect()
        } else {
            Vec::new()
        };
        let reg_primal = regul.primal;
        let dense_thr = dense_col_threshold(A, qp);
        struct OutPtr(*mut f64);
        // SAFETY: each chunk writes only the entries valNE[ptrNE[row]..ptrNE[row+1]] of its own
        // rows, and the row ranges of different chunks are disjoint.
        unsafe impl Sync for OutPtr {}
        let out = OutPtr(self.valNE.as_mut_ptr());
        let out = &out;
        let run = |c0: Int, c1: Int| {
            let mut work: Vec<f64> = vec![0.0; m as usize];
            let r0 = (c0 as i64 * m as i64 / nchunks as i64) as Int;
            let r1 = (c1 as i64 * m as i64 / nchunks as i64) as Int;
            for row in r0..r1 {
                // go along the entries of the row, and then down each column.
                // this builds the lower triangular part of the row-th column of AAt.

                let mut dense_diag = 0.0f64;
                let mut has_sparse = false;
                let mut el = ptrA_rw[row as usize];
                while el < ptrA_rw[(row + 1) as usize] {
                    let col = idxA_rw[el as usize];
                    let corr = corr_A[el as usize];

                    let mut denom: f64 = if scaling.is_empty() {
                        1.0
                    } else {
                        scaling[col as usize]
                    };
                    denom += reg_primal;
                    if qp {
                        denom += sense * qdiag[col as usize];
                    }

                    let mult: f64 = 1.0 / denom;
                    let row_value: f64 = mult * a_val[corr as usize];
                    if a_start[(col + 1) as usize] - a_start[col as usize] > dense_thr {
                        // dense column: only its diagonal contribution enters the factor, so
                        // the preconditioner stays positive definite (ipm_hipo_dense_cols)
                        dense_diag += row_value * a_val[corr as usize];
                        el += 1;
                        continue;
                    }
                    has_sparse = true;

                    // for each nonzero in the row, go down corresponding column, starting
                    // from current position
                    let mut colEl = corr;
                    while colEl < a_start[(col + 1) as usize] {
                        let row2 = a_index[colEl as usize];

                        // row2 is guaranteed to be larger or equal than row
                        // (provided that the columns of A are sorted)

                        // compute and accumulate value
                        let value: f64 = row_value * a_val[colEl as usize];
                        work[row2 as usize] += value;
                        colEl += 1;
                    }
                    el += 1;
                }
                // intersection of row with rows below finished.
                // with the Schur complement the dense part is added back exactly, and the row's
                // own regularisation keeps the pivot positive; the PCG lever needs its diagonal
                if !has_sparse && !crate::ipm::hipo::blas::DENSE_SCHUR.load(Ordering::Relaxed) {
                    work[row as usize] += dense_diag;
                }

                // read from work, using indices of column "row" of AAt
                let mut el = ptrNE[row as usize];
                while el < ptrNE[(row + 1) as usize] {
                    let index = rowsNE[el as usize];
                    // SAFETY: see OutPtr; el is inside this row's range of valNE.
                    unsafe { *out.0.add(el as usize) = work[index as usize] };
                    work[index as usize] = 0.0;
                    el += 1;
                }
            }
        };
        if par {
            crate::parallel::highs_parallel::for_each(0, nchunks, &run, 1);
        } else {
            run(0, 1);
        }

        info.matrix_time += clock.stop();

        kStatusOk
    }

    /// HiGHS: KktMatrix.cpp::KktMatrix::freeASmemory (lines 248-253)
    pub fn freeASmemory(&mut self) {
        // Give up memory used for AS.
        freeVector(&mut self.ptrAS);
        freeVector(&mut self.rowsAS);
        freeVector(&mut self.valAS);
    }

    /// HiGHS: KktMatrix.cpp::KktMatrix::freeNEmemory (lines 255-263)
    pub fn freeNEmemory(&mut self) {
        // Give up memory used for NE.
        freeVector(&mut self.ptrNE);
        freeVector(&mut self.rowsNE);
        freeVector(&mut self.valNE);
        freeVector(&mut self.ptrA_rw);
        freeVector(&mut self.idxA_rw);
        freeVector(&mut self.corr_A);
    }
}

impl Default for KktMatrix {
    fn default() -> Self {
        Self::new()
    }
}

impl KktAsPart<'_> {
    /// HiGHS: KktMatrix.cpp::KktMatrix::buildASstructure (lines 11-71). C++ `info.AS_structure_time`
    /// is `AS_structure_time`.
    pub fn buildASstructure(
        &mut self,
        model: &Model,
        AS_structure_time: &mut f64,
        logger: &Logger,
    ) -> Int {
        // Build lower triangular structure of the augmented system.
        // Build values of AS that will not change during the iterations.
        let clock = Clock::new();

        let A = model.A();
        let Q = model.Q();
        let n = A.num_col_;
        let m = A.num_row_;
        let nzA = A.num_nz();
        let nzQ = Q.num_nz();

        logger.printInfo("Building AS structure\n");

        let nzBlock11: Int = if model.qp() { nzQ } else { n };

        // AS matrix must fit into HighsInt
        if (nzBlock11 as Int64) + (m as Int64) + (nzA as Int64) >= kHighsIInf as Int64 {
            return kStatusOverflow;
        }

        self.ptrAS.resize((n + m + 1) as usize, 0);
        self.rowsAS.resize((nzBlock11 + nzA + m) as usize, 0);
        self.valAS.resize((nzBlock11 + nzA + m) as usize, 0.0);

        let mut next: Int = 0;

        for i in 0..n {
            // diagonal element
            self.rowsAS[next as usize] = i;
            next += 1;

            // column of Q
            if model.qp() {
                debug_assert!(Q.index_[Q.start_[i as usize] as usize] == i);
                let mut el = Q.start_[i as usize] + 1;
                while el < Q.start_[(i + 1) as usize] {
                    self.rowsAS[next as usize] = Q.index_[el as usize];
                    self.valAS[next as usize] = -Q.value_[el as usize]; // values of AS that will not change
                    next += 1;
                    el += 1;
                }
            }

            // column of A
            let mut el = A.start_[i as usize];
            while el < A.start_[(i + 1) as usize] {
                self.rowsAS[next as usize] = n + A.index_[el as usize];
                self.valAS[next as usize] = A.value_[el as usize]; // values of AS that will not change
                next += 1;
                el += 1;
            }

            self.ptrAS[(i + 1) as usize] = next;
        }

        // 2,2 block
        for i in 0..m {
            self.rowsAS[next as usize] = n + i;
            next += 1;
            self.ptrAS[(n + i + 1) as usize] = self.ptrAS[(n + i) as usize] + 1;
        }

        *AS_structure_time = clock.stop();

        kStatusOk
    }
}

impl KktNePart<'_> {
    /// HiGHS: KktMatrix.cpp::KktMatrix::buildNEstructure (lines 92-189). C++ `info.NE_structure_time`
    /// is `NE_structure_time`.
    pub fn buildNEstructure(
        &mut self,
        model: &Model,
        NE_structure_time: &mut f64,
        logger: &Logger,
    ) -> Int {
        // Build lower triangular structure of AAt.
        // This approach uses a column-wise copy of A, a partial row-wise copy and a
        // vector of corresponding indices.
        // NB: A must have sorted columns for this to work
        let clock = Clock::new();

        let A = model.A();
        let n = A.num_col_;
        let m = A.num_row_;
        let nzA = A.num_nz();

        logger.printInfo("Building NE structure\n");

        // create partial row-wise representation without values, and array or
        // corresponding indices between cw and rw representation

        self.ptrA_rw.clear();
        self.ptrA_rw.resize((m + 1) as usize, 0);
        self.idxA_rw.clear();
        self.idxA_rw.resize(nzA as usize, 0);

        // pointers of row-start
        for el in 0..nzA {
            self.ptrA_rw[(A.index_[el as usize] + 1) as usize] += 1;
        }
        for i in 0..m {
            self.ptrA_rw[(i + 1) as usize] += self.ptrA_rw[i as usize];
        }

        let mut temp: Vec<Int> = self.ptrA_rw.clone();
        self.corr_A.clear();
        self.corr_A.resize(nzA as usize, 0);

        // rw-indices and corresponding indices created together
        for col in 0..n {
            let mut el = A.start_[col as usize];
            while el < A.start_[(col + 1) as usize] {
                let row = A.index_[el as usize];

                self.corr_A[temp[row as usize] as usize] = el;
                self.idxA_rw[temp[row as usize] as usize] = col;
                temp[row as usize] += 1;
                el += 1;
            }
        }

        self.ptrNE.clear();
        self.rowsNE.clear();

        // ptr is allocated its exact size
        self.ptrNE.resize((m + 1) as usize, 0);

        // keep track if given entry is nonzero, in column considered
        let mut is_nz: Vec<bool> = vec![false; m as usize];

        // temporary storage of indices
        let mut temp_index: Vec<Int> = vec![0; m as usize];

        let dense_thr = dense_col_threshold(A, model.qp());
        for row in 0..m {
            // go along the entries of the row, and then down each column.
            // this builds the lower triangular part of the row-th column of AAt.

            let mut nz_in_col: Int = 0;
            if dense_thr != Int::MAX {
                // a row covered only by dense columns would have no diagonal entry
                is_nz[row as usize] = true;
                temp_index[0] = row;
                nz_in_col = 1;
            }

            let mut el = self.ptrA_rw[row as usize];
            while el < self.ptrA_rw[(row + 1) as usize] {
                let col = self.idxA_rw[el as usize];
                if A.start_[(col + 1) as usize] - A.start_[col as usize] > dense_thr {
                    el += 1;
                    continue;
                }
                let corr = self.corr_A[el as usize];

                // for each nonzero in the row, go down corresponding column, starting
                // from current position
                let mut colEl = corr;
                while colEl < A.start_[(col + 1) as usize] {
                    let row2 = A.index_[colEl as usize];

                    // row2 is guaranteed to be larger or equal than row
                    // (provided that the columns of A are sorted)

                    // save information that there is nonzero in position (row2,row).
                    if !is_nz[row2 as usize] {
                        is_nz[row2 as usize] = true;
                        temp_index[nz_in_col as usize] = row2;
                        nz_in_col += 1;
                    }
                    colEl += 1;
                }
                el += 1;
            }
            // intersection of row with rows below finished.

            // if the total number of nonzeros exceeds the maximum, return error.
            if (self.ptrNE[row as usize] as Int64) + (nz_in_col as Int64)
                >= self.NE_nz_limit.load(Ordering::Relaxed)
            {
                return kStatusOverflow;
            }

            // update pointers
            self.ptrNE[(row + 1) as usize] = self.ptrNE[row as usize] + nz_in_col;

            // now assign indices
            for i in 0..nz_in_col {
                let index = temp_index[i as usize];
                // push_back is better then reserve, because the final length is not known
                self.rowsNE.push(index);
                is_nz[index as usize] = false;
            }
        }

        *NE_structure_time = clock.stop();
        kStatusOk
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn free_memory_empties_vectors() {
        let mut k = KktMatrix::new();
        k.ptrAS = vec![0, 1];
        k.rowsAS = vec![0];
        k.valAS = vec![1.0];
        k.freeASmemory();
        assert!(k.ptrAS.is_empty() && k.rowsAS.is_empty() && k.valAS.is_empty());
        assert_eq!(k.NE_nz_limit.load(Ordering::Relaxed), kHighsIInf as Int64);
    }
}
