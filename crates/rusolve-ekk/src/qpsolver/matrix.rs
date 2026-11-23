//! HiGHS source: qpsolver/matrix.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Overloads: `mat_vec(other, target)` -> `mat_vec_into`,
//! `vec_mat(other, target)` -> `vec_mat_into`, `vec_mat(idx, val, nnz)` ->
//! `vec_mat_raw`, `extractcol(col, target)` -> `extractcol_into`; the three
//! `Matrix::append` overloads -> `append` (QpVector), `append_idx_val`
//! (`idx, val, nnz`), `append_nz` (`num_nz, index, value`).
//! `MatrixBase::num_row/num_col` are uninitialised in C++ for `MatrixBase res;`
//! (always assigned before use); `Default` gives 0.
//! `Matrix(nr, nc)` leaves `mat.start` EMPTY (the first `append` pushes 0).

use crate::prelude::*;
use crate::qpsolver::*;
use crate::util::c_printf::c_fmt_f;

/// HiGHS: matrix.hpp::MatrixBase (lines 16-155)
#[derive(Clone, Debug, Default)]
pub struct MatrixBase {
    pub num_row: HighsInt,
    pub num_col: HighsInt,
    pub start: Vec<HighsInt>,
    pub index: Vec<HighsInt>,
    pub value: Vec<f64>,
}

impl MatrixBase {
    /// HiGHS: matrix.hpp::MatrixBase::mat_vec(const QpVector&, QpVector&) (lines 23-25)
    pub fn mat_vec_into<'a>(&self, other: &QpVector, target: &'a mut QpVector) -> &'a mut QpVector {
        self.mat_vec_seq(other, target)
    }

    /// HiGHS: matrix.hpp::MatrixBase::mat_vec_seq (lines 27-38)
    pub fn mat_vec_seq<'a>(&self, other: &QpVector, target: &'a mut QpVector) -> &'a mut QpVector {
        target.reset();
        for i in 0..other.num_nz {
            let col = other.index[i as usize];
            let mut idx = self.start[col as usize];
            while idx < self.start[(col + 1) as usize] {
                let row = self.index[idx as usize];
                target.value[row as usize] += self.value[idx as usize] * other.value[col as usize];
                idx += 1;
            }
        }
        target.resparsify();
        target
    }

    /// HiGHS: matrix.hpp::MatrixBase::mat_vec(const QpVector&) (lines 40-44)
    pub fn mat_vec(&self, other: &QpVector) -> QpVector {
        let mut result = QpVector::new(self.num_row);
        self.mat_vec_into(other, &mut result);
        result
    }

    /// HiGHS: matrix.hpp::MatrixBase::vec_mat(HighsInt*, double*, HighsInt) (lines 46-82)
    pub fn vec_mat_raw(&self, idx: &[HighsInt], val: &[f64], nnz: HighsInt) -> QpVector {
        let mut result = QpVector::new(self.num_col);
        for i in 0..self.num_col {
            let mut dot = 0.0;
            // C++ lines 50-60: commented-out merge loop, dropped.
            let mut j = self.start[i as usize];
            while j < self.start[(i + 1) as usize] {
                // does the vector have an entry for index index[j]?
                let mut other_value = 0.0;
                for k in 0..nnz {
                    if idx[k as usize] == self.index[j as usize] {
                        other_value = val[k as usize];
                        break;
                    }
                }
                dot += other_value * self.value[j as usize];
                j += 1;
            }
            if dot != 0.0 {
                result.value[i as usize] = dot;
                result.index[result.num_nz as usize] = i;
                result.num_nz += 1;
            }
        }
        result
    }

    /// HiGHS: matrix.hpp::MatrixBase::vec_mat(const QpVector&, QpVector&) (lines 84-86)
    pub fn vec_mat_into<'a>(&self, other: &QpVector, target: &'a mut QpVector) -> &'a mut QpVector {
        self.vec_mat_1(other, target)
    }

    /// HiGHS: matrix.hpp::MatrixBase::vec_mat_1 (lines 88-100)
    pub fn vec_mat_1<'a>(&self, other: &QpVector, target: &'a mut QpVector) -> &'a mut QpVector {
        target.reset();
        for col in 0..self.num_col {
            let mut dot = 0.0;
            let mut j = self.start[col as usize];
            while j < self.start[(col + 1) as usize] {
                dot += other.value[self.index[j as usize] as usize] * self.value[j as usize];
                j += 1;
            }
            target.value[col as usize] = dot;
        }
        target.resparsify();
        target
    }

    /// HiGHS: matrix.hpp::MatrixBase::vec_mat(const QpVector&) (lines 102-106)
    pub fn vec_mat(&self, other: &QpVector) -> QpVector {
        let mut result = QpVector::new(self.num_col);
        self.vec_mat_into(other, &mut result);
        result
    }

    /// HiGHS: matrix.hpp::MatrixBase::tran_mat_ (lines 108-129)
    pub fn tran_mat_(&self, other: &MatrixBase) -> MatrixBase {
        let mut res = MatrixBase::default();
        res.num_row = self.num_col;
        res.num_col = other.num_col;

        res.start.push(0);
        let mut buffer_col = QpVector::new(other.num_row);
        let mut buffer_col_res = QpVector::new(self.num_col);
        for r in 0..other.num_col {
            other.extractcol_into(r, &mut buffer_col);

            self.vec_mat_into(&buffer_col, &mut buffer_col_res);
            for i in 0..buffer_col_res.num_nz {
                res.index.push(buffer_col_res.index[i as usize]);
                res.value
                    .push(buffer_col_res.value[buffer_col_res.index[i as usize] as usize]);
            }
            res.start
                .push(res.start[r as usize] + buffer_col_res.num_nz);
        }

        res
    }

    /// HiGHS: matrix.hpp::MatrixBase::extractcol(HighsInt, QpVector&) (lines 131-148)
    pub fn extractcol_into<'a>(&self, col: HighsInt, target: &'a mut QpVector) -> &'a mut QpVector {
        debug_assert!(target.dim == self.num_row);
        target.reset();

        if col >= self.num_col {
            target.index[0] = col - self.num_col;
            target.value[(col - self.num_col) as usize] = 1.0;
            target.num_nz = 1;
        } else {
            for i in 0..self.start[(col + 1) as usize] - self.start[col as usize] {
                target.index[i as usize] =
                    self.index[self.start[col as usize] as usize + i as usize];
                target.value[target.index[i as usize] as usize] =
                    self.value[self.start[col as usize] as usize + i as usize];
            }
            target.num_nz = self.start[(col + 1) as usize] - self.start[col as usize];
        }

        target
    }

    /// HiGHS: matrix.hpp::MatrixBase::extractcol(HighsInt) (lines 150-154)
    pub fn extractcol(&self, col: HighsInt) -> QpVector {
        let mut res = QpVector::new(self.num_row);
        self.extractcol_into(col, &mut res);
        res
    }

    /// HiGHS 1.16 (#3153): matrix.hpp::MatrixBase::diagonal(HighsInt) (lines 97-114),
    /// non-oracle branch only (oracle_ is not ported, so isOracle() is false)
    pub fn diagonal(&self, i_col: HighsInt) -> f64 {
        let mut diagonal_value: f64 = 0.0;
        if self.num_col > 0 {
            let mut i_el = self.start[i_col as usize];
            while i_el < self.start[(i_col + 1) as usize] {
                if self.index[i_el as usize] == i_col {
                    diagonal_value = self.value[i_el as usize];
                    break;
                }
                i_el += 1;
            }
        }
        diagonal_value
    }

    /// HiGHS 1.16 (#3153): matrix.hpp::MatrixBase::getColumn(HighsInt, HighsInt&, HighsInt*,
    /// double*) (lines 116-143), non-oracle branch only
    pub fn get_column(
        &self,
        i_col: HighsInt,
        num_entries: &mut HighsInt,
        index: &mut [HighsInt],
        value: &mut [f64],
    ) {
        if self.num_col > 0 {
            *num_entries = 0;
            let mut i_el = self.start[i_col as usize];
            while i_el < self.start[(i_col + 1) as usize] {
                index[*num_entries as usize] = self.index[i_el as usize];
                value[*num_entries as usize] = self.value[i_el as usize];
                *num_entries += 1;
                i_el += 1;
            }
        }
    }

    /// HiGHS 1.16 (#3153): matrix.hpp::MatrixBase::numNz() (lines 348-361),
    /// non-oracle branch only
    pub fn num_nz(&self) -> HighsInt {
        debug_assert!(self.start.len() >= self.num_col as usize);
        self.start[self.num_col as usize]
    }

    /// HiGHS 1.16 (#3153): matrix.hpp::MatrixBase::isDiagonal() (lines 363-375),
    /// non-oracle branch only. As in C++, this compares the stored entry count
    /// (numNz) with num_col.
    pub fn is_diagonal(&self) -> bool {
        debug_assert!(self.num_col == self.num_row);
        self.num_nz() == self.num_col
    }
}

/// HiGHS: matrix.hpp::Matrix (lines 157-340)
#[derive(Clone, Debug, Default)]
pub struct Matrix {
    tran: MatrixBase,
    has_transpose: bool,
    pub mat: MatrixBase,
}

impl Matrix {
    /// HiGHS: matrix.hpp::Matrix::transpose (lines 162-196)
    fn transpose(&mut self) {
        if !self.has_transpose {
            let mut row_indices: Vec<Vec<HighsInt>> = vec![Vec::new(); self.mat.num_row as usize];
            let mut row_values: Vec<Vec<f64>> = vec![Vec::new(); self.mat.num_row as usize];

            for col in 0..self.mat.num_col {
                let mut entry = self.mat.start[col as usize];
                while entry < self.mat.start[(col + 1) as usize] {
                    let row = self.mat.index[entry as usize];
                    let val = self.mat.value[entry as usize];
                    row_indices[row as usize].push(col);
                    row_values[row as usize].push(val);
                    entry += 1;
                }
            }
            self.tran.start.clear();
            self.tran.index.clear();
            self.tran.value.clear();
            self.tran.start.reserve(self.mat.num_row as usize + 1);
            self.tran.index.reserve(self.mat.index.len());
            self.tran.value.reserve(self.mat.value.len());

            self.tran.start.push(0);
            for row in 0..self.mat.num_row {
                self.tran
                    .index
                    .extend_from_slice(&row_indices[row as usize]);
                self.tran.value.extend_from_slice(&row_values[row as usize]);

                self.tran.start.push(
                    self.tran.start[row as usize] + row_indices[row as usize].len() as HighsInt,
                );
            }

            self.tran.num_col = self.mat.num_row;
            self.tran.num_row = self.mat.num_col;
        }
    }

    /// HiGHS: matrix.hpp::Matrix::Matrix(HighsInt, HighsInt) (lines 201-204)
    pub fn new(nr: HighsInt, nc: HighsInt) -> Self {
        let mut m = Matrix::default();
        m.mat.num_row = nr;
        m.mat.num_col = nc;
        m
    }

    /// HiGHS: matrix.hpp::Matrix::Matrix(const MatrixBase&, bool) (lines 206-211)
    pub fn from_base(m: &MatrixBase, needstran: bool) -> Self {
        let mut result = Matrix::default();
        result.mat = m.clone();
        // C++ lines 208-210: commented-out transpose() call, dropped.
        let _ = needstran;
        result
    }

    /// HiGHS: matrix.hpp::Matrix::append(const QpVector&) (lines 213-224)
    pub fn append(&mut self, vec: &QpVector) {
        if self.mat.num_col == 0 && self.mat.start.len() == 0 {
            self.mat.start.push(0);
        }
        for i in 0..vec.num_nz {
            self.mat.index.push(vec.index[i as usize]);
            self.mat
                .value
                .push(vec.value[vec.index[i as usize] as usize]);
        }
        self.mat
            .start
            .push(self.mat.start[self.mat.num_col as usize] + vec.num_nz);
        self.mat.num_col += 1;
        self.has_transpose = false;
    }

    /// HiGHS: matrix.hpp::Matrix::append(HighsInt*, double*, HighsInt) (lines 226-237)
    pub fn append_idx_val(&mut self, idx: &[HighsInt], val: &[f64], nnz: HighsInt) {
        if self.mat.num_col == 0 && self.mat.start.len() == 0 {
            self.mat.start.push(0);
        }
        for i in 0..nnz {
            self.mat.index.push(idx[i as usize]);
            self.mat.value.push(val[i as usize]);
        }
        self.mat
            .start
            .push(self.mat.start[self.mat.num_col as usize] + nnz);
        self.mat.num_col += 1;
        self.has_transpose = false;
    }

    /// HiGHS: matrix.hpp::Matrix::append(HighsInt, HighsInt*, double*) (lines 239-250)
    pub fn append_nz(&mut self, num_nz: HighsInt, index: &[HighsInt], value: &[f64]) {
        if self.mat.num_col == 0 && self.mat.start.len() == 0 {
            self.mat.start.push(0);
        }
        for i in 0..num_nz {
            self.mat.index.push(index[i as usize]);
            self.mat.value.push(value[i as usize]);
        }
        self.mat
            .start
            .push(self.mat.start[self.mat.num_col as usize] + num_nz);
        self.mat.num_col += 1;
        self.has_transpose = false;
    }

    /// HiGHS: matrix.hpp::Matrix::dropcol (lines 252-267)
    pub fn dropcol(&mut self, col: HighsInt) {
        debug_assert!(col < self.mat.num_col);
        self.has_transpose = false;

        self.mat.index.drain(
            self.mat.start[col as usize] as usize..self.mat.start[(col + 1) as usize] as usize,
        );
        self.mat.value.drain(
            self.mat.start[col as usize] as usize..self.mat.start[(col + 1) as usize] as usize,
        );

        let num_elements_in_col = self.mat.start[(col + 1) as usize] - self.mat.start[col as usize];
        let mut cc = col;
        while cc < self.mat.num_col {
            self.mat.start[cc as usize] = self.mat.start[(cc + 1) as usize] - num_elements_in_col;
            cc += 1;
        }
        self.mat.start.pop();
        self.mat.num_col -= 1;
    }

    /// HiGHS: matrix.hpp::Matrix::t (lines 269-275)
    pub fn t(&mut self) -> &mut MatrixBase {
        if !self.has_transpose {
            self.transpose();
            self.has_transpose = true;
        }
        &mut self.tran
    }

    /// HiGHS: matrix.hpp::Matrix::mat_mat (lines 277-287)
    pub fn mat_mat(&mut self, other: &mut Matrix) -> Matrix {
        let mut res = Matrix::new(self.mat.num_row, 0);

        let mut buffer = QpVector::new(other.mat.num_row);
        let mut buffer2 = QpVector::new(self.mat.num_col);
        for col in 0..other.mat.num_col {
            other.mat.extractcol_into(col, &mut buffer);
            self.vec_mat_into(&buffer, &mut buffer2);
            res.append(&buffer2);
        }

        res
    }

    /// HiGHS: matrix.hpp::Matrix::tran_mat (lines 289-298)
    pub fn tran_mat(&mut self, other: &mut Matrix) -> Matrix {
        let mut res = Matrix::new(self.mat.num_col, 0);

        let mut buffer = QpVector::new(other.mat.num_row);
        let mut buffer2 = QpVector::new(self.mat.num_row);
        for col in 0..other.mat.num_col {
            other.mat.extractcol_into(col, &mut buffer);
            self.mat_vec_into(&buffer, &mut buffer2);
            res.append(&buffer2);
        }
        res
    }

    /// HiGHS: matrix.hpp::Matrix::mat_vec(const QpVector&, QpVector&) (lines 300-302)
    pub fn mat_vec_into<'a>(&self, other: &QpVector, target: &'a mut QpVector) -> &'a mut QpVector {
        self.mat.mat_vec_into(other, target)
    }

    /// HiGHS: matrix.hpp::Matrix::mat_vec(const QpVector&) (line 304)
    pub fn mat_vec(&self, other: &QpVector) -> QpVector {
        self.mat.mat_vec(other)
    }

    /// HiGHS: matrix.hpp::Matrix::vec_mat(const QpVector&) (line 306)
    pub fn vec_mat(&self, other: &QpVector) -> QpVector {
        self.mat.vec_mat(other)
    }

    /// HiGHS: matrix.hpp::Matrix::vec_mat(const QpVector&, QpVector&) (lines 308-310)
    pub fn vec_mat_into<'a>(&self, other: &QpVector, target: &'a mut QpVector) -> &'a mut QpVector {
        self.mat.vec_mat_into(other, target)
    }

    /// HiGHS: matrix.hpp::Matrix::vec_mat(HighsInt*, double*, HighsInt) (lines 312-314)
    pub fn vec_mat_raw(&self, index: &[HighsInt], value: &[f64], num_nz: HighsInt) -> QpVector {
        self.mat.vec_mat_raw(index, value, num_nz)
    }

    /// HiGHS: matrix.hpp::Matrix::report (lines 316-339)
    ///
    /// C++ default arguments: name = ""
    pub fn report(&self, name: &str) {
        if !name.is_empty() {
            print!("{}:", name);
        }
        println!("[{} x {}]", self.mat.num_row, self.mat.num_col);
        print!("start: ");
        for i in self.mat.start.iter() {
            print!("{} ", i);
        }
        println!();

        print!("index: ");
        for i in self.mat.index.iter() {
            print!("{} ", i);
        }
        println!();

        print!("value: ");
        for d in self.mat.value.iter() {
            // `%lf` is `%f`: through `c_fmt_f` so values print as libc does.
            print!("{} ", c_fmt_f(*d, 6));
        }
        println!();
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use crate::qpsolver::*;

    #[test]
    fn mat_vec_and_vec_mat() {
        // [[1, 2], [0, 3]] column-wise
        let base = MatrixBase {
            num_row: 2,
            num_col: 2,
            start: vec![0, 1, 3],
            index: vec![0, 0, 1],
            value: vec![1.0, 2.0, 3.0],
        };
        let m = Matrix::from_base(&base, true);
        let mut x = QpVector::new(2);
        x.value = vec![1.0, 1.0];
        x.resparsify();
        assert_eq!(m.mat_vec(&x).value, vec![3.0, 3.0]);
        assert_eq!(m.vec_mat(&x).value, vec![1.0, 5.0]);
        assert_eq!(base.extractcol(1).value, vec![2.0, 3.0]);
    }

    // [[1, 2, 0], [0, 3, 4]] column-wise
    fn base_2x3() -> MatrixBase {
        MatrixBase {
            num_row: 2,
            num_col: 3,
            start: vec![0, 1, 3, 4],
            index: vec![0, 0, 1, 1],
            value: vec![1.0, 2.0, 3.0, 4.0],
        }
    }

    #[test]
    fn append_transpose_dropcol() {
        // Matrix(nr, 0) leaves start empty; the first append pushes 0 (C++ lines 214-216).
        let mut m = Matrix::new(2, 0);
        assert!(m.mat.start.is_empty());
        let mut v = QpVector::new(2);
        v.value = vec![1.0, 0.0];
        v.resparsify();
        m.append(&v);
        m.append_idx_val(&[0, 1], &[2.0, 3.0], 2);
        m.append_nz(1, &[1], &[4.0]);
        assert_eq!(m.mat.num_col, 3);
        assert_eq!(m.mat.start, vec![0, 1, 3, 4]);
        assert_eq!(m.mat.index, vec![0, 0, 1, 1]);
        assert_eq!(m.mat.value, vec![1.0, 2.0, 3.0, 4.0]);

        let t = m.t();
        assert_eq!((t.num_row, t.num_col), (3, 2));
        assert_eq!(t.start, vec![0, 2, 4]);
        assert_eq!(t.index, vec![0, 1, 1, 2]);
        assert_eq!(t.value, vec![1.0, 2.0, 3.0, 4.0]);

        // dropcol clears has_transpose, so t() rebuilds (C++ lines 254, 270).
        m.dropcol(1);
        assert_eq!(m.mat.num_col, 2);
        assert_eq!(m.mat.start, vec![0, 1, 2]);
        assert_eq!(m.mat.index, vec![0, 1]);
        assert_eq!(m.mat.value, vec![1.0, 4.0]);
        let t = m.t();
        assert_eq!((t.num_row, t.num_col), (2, 2));
        assert_eq!(t.start, vec![0, 1, 2]);
        assert_eq!(t.index, vec![0, 1]);
        assert_eq!(t.value, vec![1.0, 4.0]);
    }

    #[test]
    fn extractcol_slack_and_vec_mat_raw() {
        let base = base_2x3();
        // col >= num_col: unit vector at col - num_col (C++ lines 135-138).
        let e = base.extractcol(4);
        assert_eq!(e.num_nz, 1);
        assert_eq!(e.index[0], 1);
        assert_eq!(e.value, vec![0.0, 1.0]);

        // x = (0, 2): column 0 gives dot 0 and is skipped (C++ line 75).
        let r = base.vec_mat_raw(&[1], &[2.0], 1);
        assert_eq!(r.num_nz, 2);
        assert_eq!(&r.index[..2], &[1, 2]);
        assert_eq!(r.value, vec![0.0, 6.0, 8.0]);
    }

    #[test]
    fn tran_mat_and_products() {
        let base = base_2x3();
        let ones2 = MatrixBase {
            num_row: 2,
            num_col: 1,
            start: vec![0, 2],
            index: vec![0, 1],
            value: vec![1.0, 1.0],
        };
        let r = base.tran_mat_(&ones2);
        assert_eq!((r.num_row, r.num_col), (3, 1));
        assert_eq!(r.start, vec![0, 3]);
        assert_eq!(r.index, vec![0, 1, 2]);
        assert_eq!(r.value, vec![1.0, 5.0, 4.0]);

        let mut m = Matrix::from_base(&base, false);
        let mut o2 = Matrix::from_base(&ones2, false);
        let mm = m.mat_mat(&mut o2);
        assert_eq!((mm.mat.num_row, mm.mat.num_col), (2, 1));
        assert_eq!(mm.mat.start, vec![0, 3]);
        assert_eq!(mm.mat.index, vec![0, 1, 2]);
        assert_eq!(mm.mat.value, vec![1.0, 5.0, 4.0]);

        let ones3 = MatrixBase {
            num_row: 3,
            num_col: 1,
            start: vec![0, 3],
            index: vec![0, 1, 2],
            value: vec![1.0, 1.0, 1.0],
        };
        let mut o3 = Matrix::from_base(&ones3, false);
        let tm = m.tran_mat(&mut o3);
        assert_eq!((tm.mat.num_row, tm.mat.num_col), (3, 1));
        assert_eq!(tm.mat.start, vec![0, 2]);
        assert_eq!(tm.mat.index, vec![0, 1]);
        assert_eq!(tm.mat.value, vec![3.0, 7.0]);
    }
}
