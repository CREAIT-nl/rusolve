//! HiGHS source: util/HighsSparseVectorSum.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Header-only class. `partition`/`cleanup` take closures (C++ template functors).
//! `clear` assigns `false` in C++, i.e. `HighsCDouble::from(0.0)`.

use crate::prelude::*;
use crate::util::highs_cdouble::HighsCDouble;

/// HiGHS: HighsSparseVectorSum.h::HighsSparseVectorSum
#[derive(Clone, Debug, Default)]
pub struct HighsSparseVectorSum {
    pub values: Vec<HighsCDouble>,
    pub nonzeroinds: Vec<HighsInt>,
}

impl HighsSparseVectorSum {
    /// HiGHS: HighsSparseVectorSum.h::HighsSparseVectorSum(HighsInt dimension)
    pub fn new(dimension: HighsInt) -> Self {
        let mut s = Self::default();
        s.set_dimension(dimension);
        s
    }

    /// HiGHS: HighsSparseVectorSum.h::setDimension
    pub fn set_dimension(&mut self, dimension: HighsInt) {
        self.values
            .resize(dimension as usize, HighsCDouble::default());
        let len = self.nonzeroinds.len();
        self.nonzeroinds
            .reserve((dimension as usize).saturating_sub(len));
    }

    /// HiGHS: HighsSparseVectorSum.h::add(HighsInt index, double value)
    /// A sum that cancels to exactly 0 is stored as `f64::MIN_POSITIVE`
    /// (C++ `std::numeric_limits<double>::min()`), so the index stays listed.
    // PORT NOTE (perf): `#[inline]`, as the C++ `HighsSparseVectorSum`
    // header inline it translates; `transform`/`untransform` call it per nonzero.
    #[inline]
    pub fn add(&mut self, index: HighsInt, value: f64) {
        debug_assert!(index >= 0 && (index as usize) < self.values.len());
        let i = index as usize;
        if self.values[i] != 0.0 {
            self.values[i] += value;
        } else {
            self.values[i] = HighsCDouble::from(value);
            self.nonzeroinds.push(index);
        }

        if self.values[i] == 0.0 {
            self.values[i] = HighsCDouble::from(f64::MIN_POSITIVE);
        }
    }

    /// HiGHS: HighsSparseVectorSum.h::add(HighsInt index, HighsCDouble value)
    pub fn add_cdouble(&mut self, index: HighsInt, value: HighsCDouble) {
        let i = index as usize;
        if self.values[i] != 0.0 {
            self.values[i] += value;
        } else {
            self.values[i] = value;
            self.nonzeroinds.push(index);
        }

        if self.values[i] == 0.0 {
            self.values[i] = HighsCDouble::from(f64::MIN_POSITIVE);
        }
    }

    /// HiGHS: HighsSparseVectorSum.h::getNonzeros
    pub fn get_nonzeros(&self) -> &[HighsInt] {
        &self.nonzeroinds
    }

    /// HiGHS: HighsSparseVectorSum.h::getValue
    // PORT NOTE (perf): `#[inline]`, as the C++ header inlines.
    #[inline]
    pub fn get_value(&self, index: HighsInt) -> f64 {
        f64::from(self.values[index as usize])
    }

    /// HiGHS latest (#3223, 7db1c8180): HighsSparseVectorSum.h::getValue<HighsCDouble> (lines 59-62)
    pub fn get_value_cd(&self, index: HighsInt) -> HighsCDouble {
        self.values[index as usize]
    }

    /// HiGHS: HighsSparseVectorSum.h::clear
    pub fn clear(&mut self) {
        if 10 * self.nonzeroinds.len() < 3 * self.values.len() {
            for &i in self.nonzeroinds.iter() {
                self.values[i as usize] = HighsCDouble::from(0.0);
            }
        } else {
            self.values.fill(HighsCDouble::from(0.0));
        }

        self.nonzeroinds.clear();
    }

    /// HiGHS: HighsSparseVectorSum.h::partition<Pred> (`std::partition` on
    /// `nonzeroinds`; returns the size of the true-part). `std::partition` is
    /// not stable; reproduce libstdc++'s forward/backward swap order exactly
    /// if a caller depends on the order.
    pub fn partition<F: FnMut(&HighsInt) -> bool>(&mut self, pred: F) -> HighsInt {
        // libstdc++ `std::partition` on random-access iterators dispatches to the
        // bidirectional overload: forward scan for a false element, backward scan
        // for a true element, swap, repeat. `last` is exclusive.
        let mut pred = pred;
        let v = &mut self.nonzeroinds;
        let mut first: usize = 0;
        let mut last: usize = v.len();
        loop {
            loop {
                if first == last {
                    return first as HighsInt;
                } else if pred(&v[first]) {
                    first += 1;
                } else {
                    break;
                }
            }
            last -= 1;
            loop {
                if first == last {
                    return first as HighsInt;
                } else if !pred(&v[last]) {
                    last -= 1;
                } else {
                    break;
                }
            }
            v.swap(first, last);
            first += 1;
        }
    }

    /// HiGHS: HighsSparseVectorSum.h::cleanup<IsZero>
    pub fn cleanup<F: FnMut(HighsInt, f64) -> bool>(&mut self, is_zero: F) {
        let mut is_zero = is_zero;
        let mut num_nz: HighsInt = self.nonzeroinds.len() as HighsInt;

        let mut i: HighsInt = num_nz - 1;
        while i >= 0 {
            let pos: HighsInt = self.nonzeroinds[i as usize];
            let val: f64 = f64::from(self.values[pos as usize]);

            if is_zero(pos, val) {
                self.values[pos as usize] = HighsCDouble::from(0.0);
                num_nz -= 1;
                self.nonzeroinds.swap(num_nz as usize, i as usize);
            }
            i -= 1;
        }

        self.nonzeroinds.truncate(num_nz as usize);
    }
}
