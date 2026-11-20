//! HiGHS source: qpsolver/qpvector.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! C++ operators become named methods (`operator+` -> `op_add`,
//! `operator-` -> `op_sub`, unary `operator-` -> `op_neg`, `operator*(double)`
//! -> `op_mul_scalar`, `operator*(QpVector)` -> `dot`, `operator+=` ->
//! `op_add_assign`, `operator*=` -> `op_mul_assign`). The copy constructor
//! and `operator=` are both wholesale copies: `Clone` / `clone_from`
//! (`op_assign` keeps the C++ name visible at call sites).
//! `Default` (dim 0) is Rust-only, for `std::mem::take`.

use crate::prelude::*;
use crate::util::c_printf::c_fmt_f;

/// HiGHS: qpvector.hpp::QpVector (lines 18-240)
#[derive(Clone, Debug, Default, PartialEq)]
pub struct QpVector {
    pub num_nz: HighsInt,
    pub dim: HighsInt,
    pub index: Vec<HighsInt>,
    pub value: Vec<f64>,
}

impl QpVector {
    /// HiGHS: qpvector.hpp::QpVector::QpVector(HighsInt) (lines 24-28)
    pub fn new(d: HighsInt) -> Self {
        QpVector {
            num_nz: 0,
            dim: d,
            index: vec![0; d as usize],
            value: vec![0.0; d as usize],
        }
    }

    /// HiGHS: qpvector.hpp::QpVector::reset (lines 33-39)
    pub fn reset(&mut self) {
        for i in 0..self.num_nz {
            self.value[self.index[i as usize] as usize] = 0.0;
            self.index[i as usize] = 0;
        }
        self.num_nz = 0;
    }

    /// HiGHS: qpvector.hpp::QpVector::repopulate (lines 41-49)
    pub fn repopulate(&mut self, other: &QpVector) -> &mut QpVector {
        self.reset();
        for i in 0..other.num_nz {
            self.index[i as usize] = other.index[i as usize];
            self.value[self.index[i as usize] as usize] =
                other.value[other.index[i as usize] as usize];
        }
        self.num_nz = other.num_nz;
        self
    }

    /// HiGHS: qpvector.hpp::QpVector::operator= (lines 51-57)
    pub fn op_assign(&mut self, other: &QpVector) -> &mut QpVector {
        self.num_nz = other.num_nz;
        self.dim = other.dim;
        self.index = other.index.clone();
        self.value = other.value.clone();
        self
    }

    /// HiGHS: qpvector.hpp::QpVector::unit(HighsInt, HighsInt, QpVector&) (lines 59-65)
    pub fn unit_into(dim: HighsInt, u: HighsInt, target: &mut QpVector) -> &mut QpVector {
        target.reset();
        target.index[0] = u;
        target.value[u as usize] = 1.0;
        target.num_nz = 1;
        target
    }

    /// HiGHS: qpvector.hpp::QpVector::unit(HighsInt, HighsInt) (lines 67-73)
    pub fn unit(dim: HighsInt, u: HighsInt) -> QpVector {
        let mut vec = QpVector::new(dim);
        vec.index[0] = u;
        vec.value[u as usize] = 1.0;
        vec.num_nz = 1;
        vec
    }

    /// HiGHS: qpvector.hpp::QpVector::report (lines 75-83)
    ///
    /// C++ default arguments: name = ""
    pub fn report(&self, name: &str) {
        if !name.is_empty() {
            print!("{}: ", name);
        }
        for i in 0..self.num_nz {
            // `%lf` is `%f`: through `c_fmt_f` so non-finite values print as libc does.
            print!(
                "[{}] {} ",
                self.index[i as usize],
                c_fmt_f(self.value[self.index[i as usize] as usize], 6)
            );
        }
        println!();
    }

    /// HiGHS: qpvector.hpp::QpVector::norm2 (lines 85-93)
    pub fn norm2(&self) -> f64 {
        let mut val = 0.0;
        for i in 0..self.num_nz {
            val += self.value[self.index[i as usize] as usize]
                * self.value[self.index[i as usize] as usize];
        }
        val
    }

    /// HiGHS: qpvector.hpp::QpVector::sanitize (lines 95-107)
    ///
    /// C++ default arguments: threshold = 1e-14
    pub fn sanitize(&mut self, threshold: f64) {
        let mut new_idx: HighsInt = 0;
        for i in 0..self.num_nz {
            if self.value[self.index[i as usize] as usize].abs() > threshold {
                self.index[new_idx as usize] = self.index[i as usize];
                new_idx += 1;
            } else {
                self.value[self.index[i as usize] as usize] = 0.0;
                self.index[i as usize] = 0;
            }
        }
        self.num_nz = new_idx;
    }

    /// HiGHS: qpvector.hpp::QpVector::resparsify (lines 109-116)
    pub fn resparsify(&mut self) {
        self.num_nz = 0;
        for i in 0..self.dim {
            if self.value[i as usize] != 0.0 {
                self.index[self.num_nz as usize] = i;
                self.num_nz += 1;
            }
        }
    }

    /// HiGHS: qpvector.hpp::QpVector::scale (lines 118-123)
    pub fn scale(&mut self, a: f64) -> &mut QpVector {
        for i in 0..self.num_nz {
            self.value[self.index[i as usize] as usize] *= a;
        }
        self
    }

    /// HiGHS: qpvector.hpp::QpVector::saxpy(double, double, const QpVector&) (lines 125-129)
    pub fn saxpy_ab(&mut self, a: f64, b: f64, x: &QpVector) -> &mut QpVector {
        self.scale(a);
        self.saxpy(b, x);
        self
    }

    /// HiGHS: qpvector.hpp::QpVector::saxpy(double, const QpVector&) (lines 131-142)
    pub fn saxpy(&mut self, a: f64, x: &QpVector) -> &mut QpVector {
        self.sanitize(0.0);
        for i in 0..x.num_nz {
            if self.value[x.index[i as usize] as usize] == 0.0 {
                self.index[self.num_nz as usize] = x.index[i as usize];
                self.num_nz += 1;
            }
            self.value[x.index[i as usize] as usize] += a * x.value[x.index[i as usize] as usize];
        }
        self.resparsify();
        // sanitize(0.0); // commented out in C++ (line 140)
        self
    }

    /// HiGHS: qpvector.hpp::QpVector::operator+ (lines 151-162)
    pub fn op_add(&self, other: &QpVector) -> QpVector {
        let mut result = QpVector::new(self.dim);
        for i in 0..self.dim {
            result.value[i as usize] = self.value[i as usize] + other.value[i as usize];
            if result.value[i as usize] != 0.0 {
                result.index[result.num_nz as usize] = i;
                result.num_nz += 1;
            }
        }
        result
    }

    /// HiGHS: qpvector.hpp::QpVector::operator-(const QpVector&) (lines 164-175)
    pub fn op_sub(&self, other: &QpVector) -> QpVector {
        let mut result = QpVector::new(self.dim);
        for i in 0..self.dim {
            result.value[i as usize] = self.value[i as usize] - other.value[i as usize];
            if result.value[i as usize] != 0.0 {
                result.index[result.num_nz as usize] = i;
                result.num_nz += 1;
            }
        }
        result
    }

    /// HiGHS: qpvector.hpp::QpVector::operator-() (lines 177-187)
    pub fn op_neg(&self) -> QpVector {
        let mut result = QpVector::new(self.dim);
        for i in 0..self.num_nz {
            result.index[i as usize] = self.index[i as usize];
            result.value[self.index[i as usize] as usize] =
                -self.value[self.index[i as usize] as usize];
        }
        result.num_nz = self.num_nz;
        result
    }

    /// HiGHS: qpvector.hpp::QpVector::operator*(const double) (lines 189-199)
    pub fn op_mul_scalar(&self, d: f64) -> QpVector {
        let mut result = QpVector::new(self.dim);
        for i in 0..self.num_nz {
            result.index[i as usize] = self.index[i as usize];
            result.value[self.index[i as usize] as usize] =
                d * self.value[self.index[i as usize] as usize];
        }
        result.num_nz = self.num_nz;
        result
    }

    /// HiGHS: qpvector.hpp::QpVector::dot(const QpVector&) (lines 201-208);
    /// also `operator*(const QpVector&)` (line 210)
    pub fn dot(&self, other: &QpVector) -> f64 {
        let mut dot = 0.0;
        for i in 0..self.num_nz {
            dot += self.value[self.index[i as usize] as usize]
                * other.value[self.index[i as usize] as usize];
        }
        dot
    }

    /// HiGHS: qpvector.hpp::QpVector::dot(const HighsInt*, const double*, HighsInt) (lines 212-219)
    pub fn dot_raw(&self, idx: &[HighsInt], val: &[f64], nnz: HighsInt) -> f64 {
        let mut dot = 0.0;
        for i in 0..nnz {
            dot += self.value[idx[i as usize] as usize] * val[i as usize];
        }
        dot
    }

    /// HiGHS: qpvector.hpp::QpVector::operator+= (lines 221-231)
    pub fn op_add_assign(&mut self, other: &QpVector) -> &mut QpVector {
        // sanitize(); // commented out in C++ (line 222)
        for i in 0..other.num_nz {
            // if (value[other.index[i]] == 0.0) { // commented out in C++ (lines 224-226)
            //    index[num_nz++] = other.index[i];
            // }
            self.value[other.index[i as usize] as usize] +=
                other.value[other.index[i as usize] as usize];
        }
        self.resparsify();
        self
    }

    /// HiGHS: qpvector.hpp::QpVector::operator*= (lines 233-239)
    pub fn op_mul_assign(&mut self, d: f64) -> &mut QpVector {
        for i in 0..self.num_nz {
            self.value[self.index[i as usize] as usize] *= d;
        }
        self
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
    fn vector_algebra() {
        let v = QpVector::new(3);
        assert_eq!((v.dim, v.num_nz), (3, 0));
        assert_eq!(v.value, vec![0.0; 3]);
        let mut a = QpVector::unit(3, 1);
        assert_eq!((a.num_nz, a.value[1]), (1, 1.0));
        let mut b = QpVector::new(3);
        b.value = vec![3.0, 4.0, 0.0];
        b.resparsify();
        assert_eq!(b.num_nz, 2);
        // norm2 is the squared norm (qpvector.hpp 85-93)
        assert_eq!(b.norm2(), 25.0);
        assert_eq!(a.dot(&b), 4.0);
        a.saxpy(2.0, &b);
        assert_eq!(a.value, vec![6.0, 9.0, 0.0]);
        assert_eq!(a.op_sub(&b).value, vec![3.0, 5.0, 0.0]);
    }

    /// Branches `vector_algebra` does not reach; values hand-derived from
    /// qpvector.hpp (sanitize 161-173, repopulate 107-115, unit 125-131,
    /// saxpy 191-208, unary - 243-253, * 255-265, + 217-228, += 287-297,
    /// *= 299-305, dot(idx, val, nnz) 278-285, = 117-123).
    #[test]
    fn vector_algebra_more() {
        let mut x = QpVector::new(4);
        x.value = vec![0.0, 1e-15, 2.0, 0.0];
        x.resparsify();
        assert_eq!((x.num_nz, x.index.clone()), (2, vec![1, 2, 0, 0]));
        // index[0] is zeroed, then overwritten by the kept entry; index[1] stays.
        x.sanitize(1e-14);
        assert_eq!(x.num_nz, 1);
        assert_eq!(x.index, vec![2, 2, 0, 0]);
        assert_eq!(x.value, vec![0.0, 0.0, 2.0, 0.0]);

        let mut y = QpVector::new(4);
        y.value[3] = 5.0;
        y.index[0] = 3;
        y.num_nz = 1;
        y.repopulate(&x);
        assert_eq!((y.num_nz, y.index[0]), (1, 2));
        assert_eq!(y.value, vec![0.0, 0.0, 2.0, 0.0]);

        QpVector::unit_into(4, 3, &mut y);
        assert_eq!((y.num_nz, y.index[0]), (1, 3));
        assert_eq!(y.value, vec![0.0, 0.0, 0.0, 1.0]);

        y.saxpy_ab(2.0, 3.0, &x);
        assert_eq!(y.num_nz, 2);
        assert_eq!(&y.index[..2], &[2, 3]);
        assert_eq!(y.value, vec![0.0, 0.0, 6.0, 2.0]);

        let n = y.op_neg();
        assert_eq!((n.num_nz, &n.index[..2]), (2, &[2, 3][..]));
        assert_eq!(n.value, vec![0.0, 0.0, -6.0, -2.0]);
        let m = y.op_mul_scalar(0.5);
        assert_eq!(m.value, vec![0.0, 0.0, 3.0, 1.0]);
        let s = m.op_add(&x);
        assert_eq!((s.num_nz, &s.index[..2]), (2, &[2, 3][..]));
        assert_eq!(s.value, vec![0.0, 0.0, 5.0, 1.0]);

        y.op_add_assign(&x);
        assert_eq!(y.value, vec![0.0, 0.0, 8.0, 2.0]);
        y.op_mul_assign(0.5);
        assert_eq!(y.value, vec![0.0, 0.0, 4.0, 1.0]);
        assert_eq!(y.dot_raw(&[3, 2], &[2.0, 1.0], 2), 6.0);

        // Exact cancellation drops the entry through resparsify.
        y.saxpy(-2.0, &x);
        assert_eq!((y.num_nz, y.index[0]), (1, 3));
        assert_eq!(y.value, vec![0.0, 0.0, 0.0, 1.0]);

        let mut z = QpVector::new(1);
        z.op_assign(&y);
        assert_eq!(z, y);
    }
}
