//! HiGHS source: ipm/ipx/indexed_vector.h, indexed_vector.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::*;

/// HiGHS: ipx/indexed_vector.h::IndexedVector
///
/// `nnz_ < 0` means "pattern unknown" (`InvalidatePattern`).
#[derive(Debug, Clone, Default)]
pub struct IndexedVector {
    elements_: Vector,
    pattern_: Vec<Int>,
    nnz_: Int,
}

impl IndexedVector {
    /// HiGHS: ipx/indexed_vector.cc::IndexedVector::IndexedVector (line 5)
    ///
    /// C++ default arguments: `dim = 0`.
    pub fn new(dim: Int) -> Self {
        IndexedVector {
            elements_: Vector::new(dim),
            pattern_: vec![0; dim as usize],
            nnz_: 0,
        }
    }
    /// HiGHS: ipx/indexed_vector.h::IndexedVector::dim
    pub fn dim(&self) -> Int {
        self.elements_.size()
    }
    /// HiGHS: ipx/indexed_vector.h::IndexedVector::operator[] const
    pub fn get(&self, i: Int) -> f64 {
        self.elements_[i as usize]
    }
    /// HiGHS: ipx/indexed_vector.h::IndexedVector::operator[]
    pub fn get_mut(&mut self, i: Int) -> &mut f64 {
        &mut self.elements_[i as usize]
    }
    /// HiGHS: ipx/indexed_vector.cc::IndexedVector::sparse (lines 7-9)
    pub fn sparse(&self) -> bool {
        self.nnz() >= 0 && (self.nnz() as f64) <= kHypersparseThreshold * (self.dim() as f64)
    }
    pub fn nnz(&self) -> Int {
        self.nnz_
    }
    /// HiGHS: ipx/indexed_vector.h::IndexedVector::InvalidatePattern
    pub fn invalidate_pattern(&mut self) {
        self.nnz_ = -1;
    }
    pub fn set_nnz(&mut self, new_nnz: Int) {
        self.nnz_ = new_nnz;
    }
    pub fn elements(&self) -> &[f64] {
        &self.elements_
    }
    pub fn elements_mut(&mut self) -> &mut [f64] {
        &mut self.elements_
    }
    pub fn pattern(&self) -> &[Int] {
        &self.pattern_
    }
    pub fn pattern_mut(&mut self) -> &mut [Int] {
        &mut self.pattern_
    }
    /// Both arrays at once (callers that fill pattern and values together,
    /// e.g. the `basiclu_solve_for_update` output of `BasicLu`).
    pub fn elements_and_pattern_mut(&mut self) -> (&mut [f64], &mut [Int]) {
        (&mut self.elements_, &mut self.pattern_)
    }
    /// HiGHS: ipx/indexed_vector.cc::IndexedVector::set_to_zero (lines 11-19)
    pub fn set_to_zero(&mut self) {
        if self.sparse() {
            for p in 0..self.nnz_ {
                let i = self.pattern_[p as usize];
                self.elements_[i as usize] = 0.0;
            }
        } else {
            // C++ `elements_ = 0.0` (valarray scalar assignment)
            self.elements_.fill(0.0);
        }
        self.nnz_ = 0;
    }
}

/// HiGHS: ipx/indexed_vector.h::for_each_nonzero(const IndexedVector&, C&)
pub fn for_each_nonzero<C: FnMut(Int, f64)>(v: &IndexedVector, mut c: C) {
    if v.sparse() {
        let nnz = v.nnz();
        for p in 0..nnz {
            let i = v.pattern()[p as usize];
            c(i, v.get(i));
        }
    } else {
        let dim = v.dim();
        for i in 0..dim {
            c(i, v.get(i));
        }
    }
}

/// HiGHS: ipx/indexed_vector.h::for_each_nonzero(IndexedVector&, C&)
/// (the closure may write the element: C++ passes `double&`)
pub fn for_each_nonzero_mut<C: FnMut(Int, &mut f64)>(v: &mut IndexedVector, mut c: C) {
    if v.sparse() {
        let nnz = v.nnz();
        for p in 0..nnz {
            let i = v.pattern_[p as usize];
            c(i, &mut v.elements_[i as usize]);
        }
    } else {
        let dim = v.dim();
        for i in 0..dim {
            c(i, &mut v.elements_[i as usize]);
        }
    }
}

/// HiGHS: ipx/indexed_vector.cc::Dot(const IndexedVector&, const Vector&) (lines 21-28)
pub fn dot_indexed(x: &IndexedVector, y: &[f64]) -> f64 {
    let mut d = 0.0;
    for_each_nonzero(x, |p, f| {
        d += y[p as usize] * f;
    });
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexed_vector_sparse_threshold_and_dot() {
        let mut v = IndexedVector::new(10);
        // nnz 1 <= 0.1 * 10: sparse, so only pattern entries are visited.
        v.set_nnz(1);
        v.pattern_mut()[0] = 3;
        v.elements_mut()[3] = 2.0;
        assert!(v.sparse());
        let y: Vec<f64> = (0..10).map(|i| i as f64).collect();
        assert_eq!(dot_indexed(&v, &y), 6.0);
        // nnz 2 > 0.1 * 10: dense path, entries outside the pattern are visited too.
        v.set_nnz(2);
        assert!(!v.sparse());
        v.set_to_zero();
        assert_eq!(v.nnz(), 0);
        assert_eq!(v.elements()[3], 0.0);
    }
}
