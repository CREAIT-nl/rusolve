//! HiGHS source: ipm/ipx/sparse_matrix.h, sparse_matrix.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The IPX `SparseMatrix` (CSC with a "queue" for building column by
//! column), not `util::sparse_matrix::HighsSparseMatrix`. `SortIndices`
//! sorts `(Int, double)` pairs with HiGHS pdqsort and the default pair `<`
//! (non-branchless instantiation: `util::highs_pdqsort::pdqsort_by`).
//! Header inline members and the `.cc` bodies are implemented.

use super::utils::{infnorm, onenorm, pair_less};
use super::*;
use crate::util::highs_pdqsort::pdqsort_by;
use crate::util::highs_utils::cpp_max;

/// HiGHS: ipx/sparse_matrix.h::SparseMatrix
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SparseMatrix {
    nrow_: Int,
    colptr_: Vec<Int>,
    rowidx_: Vec<Int>,
    values_: Vec<f64>,
    rowidx_queue_: Vec<Int>,
    values_queue_: Vec<f64>,
}

impl SparseMatrix {
    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::SparseMatrix() (lines 11-13)
    pub fn new() -> Self {
        let mut a = Self::default();
        a.resize(0, 0, 0);
        a
    }
    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::SparseMatrix(Int, Int) (lines 15-17)
    pub fn with_dims(nrow: Int, ncol: Int) -> Self {
        let mut a = Self::default();
        a.resize(nrow, ncol, 0);
        a
    }
    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::SparseMatrix(Int, Int, Int) (lines 19-21)
    pub fn with_capacity(nrow: Int, ncol: Int, min_capacity: Int) -> Self {
        let mut a = Self::default();
        a.resize(nrow, ncol, min_capacity);
        a
    }

    pub fn rows(&self) -> Int {
        self.nrow_
    }
    pub fn cols(&self) -> Int {
        self.colptr_.len() as Int - 1
    }
    pub fn entries(&self) -> Int {
        *self.colptr_.last().unwrap()
    }
    /// HiGHS: ipx/sparse_matrix.h::SparseMatrix::entries(Int j)
    pub fn entries_col(&self, j: Int) -> Int {
        self.end(j) - self.begin(j)
    }
    pub fn capacity(&self) -> Int {
        self.rowidx_.len() as Int
    }

    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::reserve (lines 23-28)
    pub fn reserve(&mut self, min_capacity: Int) {
        if min_capacity > self.capacity() {
            self.rowidx_.resize(min_capacity as usize, 0);
            self.values_.resize(min_capacity as usize, 0.0);
        }
    }
    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::resize (lines 30-42)
    ///
    /// C++ default arguments: `min_capacity = 0`.
    pub fn resize(&mut self, nrow: Int, ncol: Int, min_capacity: Int) {
        debug_assert!(nrow >= 0);
        debug_assert!(ncol >= 0);
        debug_assert!(min_capacity >= 0);
        self.nrow_ = nrow;
        // colptr_.resize(ncol+1); then std::fill(.., 0): every entry becomes 0.
        self.colptr_.resize((ncol + 1) as usize, 0);
        self.colptr_.fill(0);
        self.rowidx_.resize(min_capacity as usize, 0);
        self.values_.resize(min_capacity as usize, 0.0);
    }
    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::clear (lines 44-46)
    pub fn clear(&mut self) {
        self.resize(0, 0, 0);
    }
    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::LoadFromArrays (lines 48-68)
    pub fn load_from_arrays(
        &mut self,
        nrow: Int,
        ncol: Int,
        Abegin: &[Int],
        Aend: &[Int],
        Ai: &[Int],
        Ax: &[f64],
    ) {
        let mut nz: Int = 0;
        for j in 0..ncol as usize {
            nz += Aend[j] - Abegin[j];
        }
        self.resize(nrow, ncol, nz);
        let mut put: Int = 0;
        for j in 0..ncol as usize {
            self.colptr_[j] = put;
            for p in Abegin[j] as usize..Aend[j] as usize {
                if Ax[p] != 0.0 {
                    self.rowidx_[put as usize] = Ai[p];
                    self.values_[put as usize] = Ax[p];
                    put += 1;
                }
            }
        }
        self.colptr_[ncol as usize] = put;
        self.sort_indices();
    }

    #[inline]
    pub fn begin(&self, j: Int) -> Int {
        self.colptr_[j as usize]
    }
    #[inline]
    pub fn end(&self, j: Int) -> Int {
        self.colptr_[j as usize + 1]
    }
    #[inline]
    pub fn index(&self, pos: Int) -> Int {
        self.rowidx_[pos as usize]
    }
    #[inline]
    pub fn value(&self, pos: Int) -> f64 {
        self.values_[pos as usize]
    }
    /// C++ non-const `index(pos)` (a reference); not `std::ops::IndexMut`.
    #[allow(clippy::should_implement_trait)]
    pub fn index_mut(&mut self, pos: Int) -> &mut Int {
        &mut self.rowidx_[pos as usize]
    }
    pub fn value_mut(&mut self, pos: Int) -> &mut f64 {
        &mut self.values_[pos as usize]
    }
    pub fn colptr(&self) -> &[Int] {
        &self.colptr_
    }
    pub fn rowidx(&self) -> &[Int] {
        &self.rowidx_
    }
    pub fn values(&self) -> &[f64] {
        &self.values_
    }
    pub fn colptr_mut(&mut self) -> &mut [Int] {
        &mut self.colptr_
    }
    pub fn rowidx_mut(&mut self) -> &mut [Int] {
        &mut self.rowidx_
    }
    pub fn values_mut(&mut self) -> &mut [f64] {
        &mut self.values_
    }
    /// The three arrays at once (C++ code that holds `Ap`, `Ai`, `Ax`
    /// pointers together, e.g. `RemoveDiagonal`, BASICLU `get_factors`).
    pub fn arrays_mut(&mut self) -> (&mut [Int], &mut [Int], &mut [f64]) {
        (&mut self.colptr_, &mut self.rowidx_, &mut self.values_)
    }

    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::SortIndices (lines 70-87)
    pub fn sort_indices(&mut self) {
        if self.is_sorted() {
            return;
        }
        let mut work: Vec<(Int, f64)> = vec![(0, 0.0); self.rows() as usize];
        for j in 0..self.cols() {
            let mut nz: Int = 0; // # entries in column j
            for p in self.begin(j)..self.end(j) {
                work[nz as usize] = (self.index(p), self.value(p));
                nz += 1;
            }
            pdqsort_by(&mut work[..nz as usize], |a, b| pair_less(*a, *b));
            for (k, p) in (self.begin(j)..self.end(j)).enumerate() {
                *self.index_mut(p) = work[k].0;
                *self.value_mut(p) = work[k].1;
            }
        }
    }

    pub fn push_back(&mut self, i: Int, x: f64) {
        self.rowidx_queue_.push(i);
        self.values_queue_.push(x);
    }
    pub fn queue_size(&self) -> Int {
        self.rowidx_queue_.len() as Int
    }
    pub fn qindex(&self, pos: Int) -> Int {
        self.rowidx_queue_[pos as usize]
    }
    pub fn qvalue(&self, pos: Int) -> f64 {
        self.values_queue_[pos as usize]
    }
    pub fn qindex_mut(&mut self, pos: Int) -> &mut Int {
        &mut self.rowidx_queue_[pos as usize]
    }
    pub fn qvalue_mut(&mut self, pos: Int) -> &mut f64 {
        &mut self.values_queue_[pos as usize]
    }

    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::add_column (lines 89-97)
    pub fn add_column(&mut self) {
        let nz = self.entries();
        let nznew = nz + self.queue_size();
        self.reserve(nznew);
        for (k, &i) in self.rowidx_queue_.iter().enumerate() {
            self.rowidx_[nz as usize + k] = i;
        }
        for (k, &x) in self.values_queue_.iter().enumerate() {
            self.values_[nz as usize + k] = x;
        }
        self.colptr_.push(nznew);
        self.clear_queue();
    }
    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::clear_queue (lines 99-102)
    pub fn clear_queue(&mut self) {
        self.rowidx_queue_.clear();
        self.values_queue_.clear();
    }
    /// HiGHS: ipx/sparse_matrix.cc::SparseMatrix::IsSorted (lines 104-111)
    fn is_sorted(&self) -> bool {
        for j in 0..self.cols() {
            for p in self.begin(j)..self.end(j) - 1 {
                if self.index(p) > self.index(p + 1) {
                    return false;
                }
            }
        }
        true
    }
}

/// HiGHS: ipx/sparse_matrix.cc::Transpose(const SparseMatrix&) (lines 113-117)
pub fn transpose(A: &SparseMatrix) -> SparseMatrix {
    let mut AT = SparseMatrix::new();
    transpose_into(A, &mut AT);
    AT
}

/// HiGHS: ipx/sparse_matrix.cc::Transpose(const SparseMatrix&, SparseMatrix&) (lines 119-150)
pub fn transpose_into(A: &SparseMatrix, AT: &mut SparseMatrix) {
    let m = A.rows();
    let n = A.cols();
    let nz = A.entries();
    AT.resize(n, m, nz);

    // Compute row counts of A in workspace.
    let mut work: Vec<Int> = vec![0; m as usize];
    for p in 0..nz {
        work[A.index(p) as usize] += 1;
    }

    // Set column pointers for AT.
    let (ATp, ATi, ATx) = AT.arrays_mut();
    let mut sum: Int = 0;
    for i in 0..m as usize {
        ATp[i] = sum;
        sum += work[i];
        work[i] = ATp[i];
    }
    debug_assert_eq!(sum, nz);
    ATp[m as usize] = sum;

    // Fill AT with one column of A at a time.
    // work[i] is the next free slot in column i of AT.
    for j in 0..n {
        for p in A.begin(j)..A.end(j) {
            let row = A.index(p) as usize;
            let put = work[row];
            work[row] += 1;
            ATi[put as usize] = j;
            ATx[put as usize] = A.value(p);
        }
    }
}

/// HiGHS: ipx/sparse_matrix.cc::CopyColumns (lines 152-160)
pub fn copy_columns(A: &SparseMatrix, cols: &[Int]) -> SparseMatrix {
    let mut A2 = SparseMatrix::with_dims(A.rows(), 0);
    for &j in cols {
        for p in A.begin(j)..A.end(j) {
            A2.push_back(A.index(p), A.value(p));
        }
        A2.add_column();
    }
    A2
}

/// HiGHS: ipx/sparse_matrix.cc::PermuteRows (lines 162-165)
pub fn permute_rows(A: &mut SparseMatrix, perm: &[Int]) {
    for p in 0..A.entries() {
        let row = A.index(p) as usize;
        *A.index_mut(p) = perm[row];
    }
}

/// HiGHS: ipx/sparse_matrix.h::ScaleColumn (inline)
pub fn scale_column(A: &mut SparseMatrix, j: Int, s: f64) {
    let p1 = A.begin(j);
    let p2 = A.end(j);
    for p in p1..p2 {
        *A.value_mut(p) *= s;
    }
}

/// HiGHS: ipx/sparse_matrix.cc::RemoveDiagonal (lines 167-191)
///
/// `diag` may be `nullptr`.
pub fn remove_diagonal(A: &mut SparseMatrix, diag: Option<&mut [f64]>) -> Int {
    let mut diag = diag;
    let ncol = A.cols();
    let (Ap, Ai, Ax) = A.arrays_mut();
    let mut get: Int = 0;
    let mut put: Int = 0;
    for j in 0..ncol {
        if let Some(d) = diag.as_deref_mut() {
            d[j as usize] = 0.0; // if no diagonal entry in column j
        }
        Ap[j as usize] = put;
        while get < Ap[(j + 1) as usize] {
            if Ai[get as usize] == j {
                if let Some(d) = diag.as_deref_mut() {
                    d[j as usize] = Ax[get as usize];
                }
            } else {
                Ai[put as usize] = Ai[get as usize];
                Ax[put as usize] = Ax[get as usize];
                put += 1;
            }
            get += 1;
        }
    }
    Ap[ncol as usize] = put;
    get - put
}

/// HiGHS: ipx/sparse_matrix.h::DotColumn (inline)
///
/// PORT NOTE (perf): the column's index/value
/// ranges are sliced once (one bounds check each) instead of per entry; the
/// summation order is the C++ loop's.
///
/// PORT NOTE (perf): `inline(always)`, as the C++ header function
/// is. LLVM kept the plain-hint version out of line in per-column loops such as
/// `add_normal_product` and the KKT solvers, where the call and the column-range setup cost
/// more than the short column loops.
#[inline(always)]
pub fn dot_column(A: &SparseMatrix, j: Int, rhs: &[f64]) -> f64 {
    let p1 = A.begin(j) as usize;
    let p2 = A.end(j) as usize;
    let idx = &A.rowidx_[p1..p2];
    let val = &A.values_[p1..p2];
    let mut d = 0.0;
    for (&i, &x) in idx.iter().zip(val) {
        d += rhs[i as usize] * x;
    }
    d
}

/// HiGHS: ipx/sparse_matrix.h::ScatterColumn (inline)
///
/// PORT NOTE (perf): sliced as in [`dot_column`]; same update order. `inline(always)` as
/// [`dot_column`].
#[inline(always)]
pub fn scatter_column(A: &SparseMatrix, j: Int, alpha: f64, lhs: &mut [f64]) {
    let p1 = A.begin(j) as usize;
    let p2 = A.end(j) as usize;
    let idx = &A.rowidx_[p1..p2];
    let val = &A.values_[p1..p2];
    for (&i, &x) in idx.iter().zip(val) {
        lhs[i as usize] += alpha * x;
    }
}

/// HiGHS: ipx/sparse_matrix.cc::MultiplyAdd (lines 193-208)
pub fn multiply_add(A: &SparseMatrix, rhs: &[f64], alpha: f64, lhs: &mut [f64], trans: u8) {
    let m = A.rows();
    let n = A.cols();
    if trans == b't' || trans == b'T' {
        debug_assert_eq!(rhs.len() as Int, m);
        debug_assert_eq!(lhs.len() as Int, n);
        for j in 0..n {
            lhs[j as usize] += alpha * dot_column(A, j, rhs);
        }
    } else {
        debug_assert_eq!(rhs.len() as Int, n);
        debug_assert_eq!(lhs.len() as Int, m);
        for j in 0..n {
            scatter_column(A, j, alpha * rhs[j as usize], lhs);
        }
    }
}

/// HiGHS: ipx/sparse_matrix.cc::AddNormalProduct (lines 210-221)
///
/// `D` may be `nullptr`.
pub fn add_normal_product(A: &SparseMatrix, D: Option<&[f64]>, rhs: &[f64], lhs: &mut [f64]) {
    let m = A.rows();
    let n = A.cols();
    debug_assert_eq!(rhs.len() as Int, m);
    debug_assert_eq!(lhs.len() as Int, m);
    for j in 0..n {
        let mut temp = dot_column(A, j, rhs);
        if let Some(d) = D {
            temp *= d[j as usize] * d[j as usize];
        }
        scatter_column(A, j, temp, lhs);
    }
}

/// HiGHS: ipx/sparse_matrix.cc::TriangularSolve (lines 223-300)
///
/// `uplo` is the C string `"u"`/`"l"` (only `uplo[0]` is read).
///
/// PORT NOTE (perf): each column's `Ai`/`Ax` range is
/// sliced once; entry and summation order are the C++ loops'.
pub fn triangular_solve(
    A: &SparseMatrix,
    x: &mut [f64],
    trans: u8,
    uplo: &str,
    unitdiag: i32,
) -> Int {
    let ncol = A.cols();
    let Ap = A.colptr();
    let Ai = A.rowidx();
    let Ax = A.values();
    // C: `*uplo == 'u' || *uplo == 'U'` (a NUL first byte is "lower")
    let upper = matches!(uplo.as_bytes().first(), Some(b'u') | Some(b'U'));
    let unit = unitdiag != 0;
    let mut nz: Int = 0;
    if trans == b't' || trans == b'T' {
        if upper {
            // transposed solve with upper triangular matrix
            for i in 0..ncol {
                let begin = Ap[i as usize];
                let end = Ap[(i + 1) as usize] - (if unit { 0 } else { 1 });
                let mut d = 0.0;
                let (b, e) = (begin as usize, end as usize);
                for (&i, &v) in Ai[b..e].iter().zip(&Ax[b..e]) {
                    d += x[i as usize] * v;
                }
                x[i as usize] -= d;
                if !unit {
                    debug_assert_eq!(Ai[end as usize], i);
                    x[i as usize] /= Ax[end as usize];
                }
                if x[i as usize] != 0.0 {
                    nz += 1;
                }
            }
        } else {
            // transposed solve with lower triangular matrix
            for i in (0..ncol).rev() {
                let begin = Ap[i as usize] + (if unit { 0 } else { 1 });
                let end = Ap[(i + 1) as usize];
                let mut d = 0.0;
                let (b, e) = (begin as usize, end as usize);
                for (&i, &v) in Ai[b..e].iter().zip(&Ax[b..e]) {
                    d += x[i as usize] * v;
                }
                x[i as usize] -= d;
                if !unit {
                    debug_assert_eq!(Ai[(begin - 1) as usize], i);
                    x[i as usize] /= Ax[(begin - 1) as usize];
                }
                if x[i as usize] != 0.0 {
                    nz += 1;
                }
            }
        }
    } else if upper {
        // forward solve with upper triangular matrix
        for j in (0..ncol).rev() {
            let begin = Ap[j as usize];
            let end = Ap[(j + 1) as usize] - (if unit { 0 } else { 1 });
            if !unit {
                debug_assert_eq!(Ai[end as usize], j);
                x[j as usize] /= Ax[end as usize];
            }
            let temp = x[j as usize];
            if temp != 0.0 {
                let (b, e) = (begin as usize, end as usize);
                for (&i, &v) in Ai[b..e].iter().zip(&Ax[b..e]) {
                    x[i as usize] -= v * temp;
                }
                nz += 1;
            }
        }
    } else {
        // forward solve with lower triangular matrix
        for j in 0..ncol {
            let begin = Ap[j as usize] + (if unit { 0 } else { 1 });
            let end = Ap[(j + 1) as usize];
            if !unit {
                debug_assert_eq!(Ai[(begin - 1) as usize], j);
                x[j as usize] /= Ax[(begin - 1) as usize];
            }
            let temp = x[j as usize];
            if temp != 0.0 {
                let (b, e) = (begin as usize, end as usize);
                for (&i, &v) in Ai[b..e].iter().zip(&Ax[b..e]) {
                    x[i as usize] -= v * temp;
                }
                nz += 1;
            }
        }
    }
    nz
}

/// HiGHS: ipx/sparse_matrix.cc::ForwardSolve (lines 302-305)
pub fn forward_solve(L: &SparseMatrix, U: &SparseMatrix, x: &mut [f64]) {
    triangular_solve(L, x, b'n', "lower", 1);
    triangular_solve(U, x, b'n', "upper", 0);
}

/// HiGHS: ipx/sparse_matrix.cc::BackwardSolve (lines 307-310)
pub fn backward_solve(L: &SparseMatrix, U: &SparseMatrix, x: &mut [f64]) {
    triangular_solve(U, x, b't', "upper", 0);
    triangular_solve(L, x, b't', "lower", 1);
}

/// HiGHS: ipx/sparse_matrix.cc::Onenorm(const SparseMatrix&) (lines 312-321)
pub fn onenorm_matrix(A: &SparseMatrix) -> f64 {
    let mut norm = 0.0;
    for j in 0..A.cols() {
        let mut colsum = 0.0;
        for p in A.begin(j)..A.end(j) {
            colsum += A.value(p).abs();
        }
        norm = cpp_max(norm, colsum);
    }
    norm
}

/// HiGHS: ipx/sparse_matrix.cc::Infnorm(const SparseMatrix&) (lines 323-330)
pub fn infnorm_matrix(A: &SparseMatrix) -> f64 {
    let mut rowsum = Vector::new(A.rows());
    for j in 0..A.cols() {
        for p in A.begin(j)..A.end(j) {
            rowsum[A.index(p) as usize] += A.value(p).abs();
        }
    }
    infnorm(&rowsum)
}

/// HiGHS: ipx/sparse_matrix.cc::NormestInverse (lines 332-380)
pub fn normest_inverse(A: &SparseMatrix, uplo: &str, unitdiag: i32) -> f64 {
    let m = A.rows();
    let mut x = Vector::new(m);
    debug_assert_eq!(A.rows(), A.cols());
    let upper = matches!(uplo.as_bytes().first(), Some(b'u') | Some(b'U'));
    let unit = unitdiag != 0;

    // Solve A'x=b, where the entries of b are +/-1 chosen dynamically to make x
    // large.
    if upper {
        for j in 0..m {
            let begin = A.begin(j);
            let mut end = A.end(j);
            if !unit {
                end -= 1;
            }
            let mut temp = 0.0;
            for p in begin..end {
                temp -= x[A.index(p) as usize] * A.value(p);
            }
            temp += if temp >= 0.0 { 1.0 } else { -1.0 }; // choose b[j] = 1 or b[j] = -1
            if !unit {
                debug_assert_eq!(A.index(end), j);
                temp /= A.value(end);
            }
            x[j as usize] = temp;
        }
    } else {
        for j in (0..m).rev() {
            let mut begin = A.begin(j);
            let end = A.end(j);
            if !unit {
                begin += 1;
            }
            let mut temp = 0.0;
            for p in begin..end {
                temp -= x[A.index(p) as usize] * A.value(p);
            }
            temp += if temp >= 0.0 { 1.0 } else { -1.0 }; // choose b[j] = 1 or b[j] = -1
            if !unit {
                debug_assert_eq!(A.index(begin - 1), j);
                temp /= A.value(begin - 1);
            }
            x[j as usize] = temp;
        }
    }
    let x1norm = onenorm(&x);
    let xinfnorm = infnorm(&x);

    // Solve Ay=x, solution overwrites x.
    triangular_solve(A, &mut x, b'n', uplo, unitdiag);
    let y1norm = onenorm(&x);

    cpp_max(y1norm / x1norm, xinfnorm)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2x2 matrix [[1, 0], [2, 3]] built through the queue, as IPX builds bases.
    fn small() -> SparseMatrix {
        let mut a = SparseMatrix::with_dims(2, 0);
        a.push_back(0, 1.0);
        a.push_back(1, 2.0);
        a.add_column();
        a.push_back(1, 3.0);
        a.add_column();
        a
    }

    #[test]
    fn transpose_and_queue() {
        let a = small();
        assert_eq!(a.colptr(), &[0, 2, 3]);
        let at = transpose(&a);
        assert_eq!(at.colptr(), &[0, 1, 3]);
        assert_eq!(at.rowidx(), &[0, 0, 1]);
        assert_eq!(at.values(), &[1.0, 2.0, 3.0]);
    }

    #[test]
    fn remove_diagonal_reports_and_compacts() {
        // [[5, 1], [2, 7]]: diagonal 5 and 7 are removed, two entries remain.
        let mut a = SparseMatrix::with_dims(2, 0);
        a.push_back(0, 5.0);
        a.push_back(1, 2.0);
        a.add_column();
        a.push_back(0, 1.0);
        a.push_back(1, 7.0);
        a.add_column();
        let mut diag = [0.0; 2];
        let removed = remove_diagonal(&mut a, Some(&mut diag[..]));
        assert_eq!(removed, 2);
        assert_eq!(diag, [5.0, 7.0]);
        assert_eq!(a.colptr(), &[0, 1, 2]);
        assert_eq!(a.rowidx()[..2], [1, 0]);
        assert_eq!(a.values()[..2], [2.0, 1.0]);
    }

    #[test]
    fn upper_triangular_solve_hand_computed() {
        // U = [[2, 1], [0, 4]] (column 0: row 0; column 1: rows 0 and 1).
        let mut u = SparseMatrix::with_dims(2, 0);
        u.push_back(0, 2.0);
        u.add_column();
        u.push_back(0, 1.0);
        u.push_back(1, 4.0);
        u.add_column();
        // U x = [4, 8]  =>  x1 = 8 / 4 = 2, x0 = (4 - 1 * 2) / 2 = 1.
        let mut x = [4.0, 8.0];
        let nz = triangular_solve(&u, &mut x, b'n', "upper", 0);
        assert_eq!(x, [1.0, 2.0]);
        assert_eq!(nz, 2);
    }
}
