//! HiGHS source: util/HVectorBase.h, util/HVectorBase.cpp, util/HVector.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `HVectorBase<Real>` is instantiated only as `HVector = HVectorBase<double>`
//! on the simplex path, so the port is the `double` instantiation, named
//! `HVector`. `HVectorQuad` (`HVectorBase<HighsCDouble>`) is unused there and
//! omitted; so are the `copy`/`saxpy` instantiations that mix in
//! `HighsCDouble` (template parameters `FromReal`, `RealPivX`, `RealPiv` are
//! all `double` here).
//!
//! The `next` member (the PAMI linked list of `HVector*`) is
//! [`HVector::next`]: an [`HVectorNext`] raw link, reset by `setup` and
//! `clear_scalars` exactly where C++ resets it, and copied by `Clone` (C++ copy).
//! `cwork` is `vector<char>` in C++ and is `Vec<i8>` here.

use crate::HighsInt;
use crate::lp_data::hconst::{kHighsTiny, kHighsZero};

/// HiGHS: HVectorBase.h::HVectorBase<double> (alias `HVector`)
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HVector {
    /// dimension of the vector
    pub size: HighsInt,
    /// number of nonzeros (-1 = index not valid)
    pub count: HighsInt,
    /// packed indices
    pub index: Vec<HighsInt>,
    /// full (dense) values
    pub array: Vec<f64>,
    /// synthetic clock for operations with this vector
    pub synthetic_tick: f64,
    /// character work space (C++ `vector<char>`)
    pub cwork: Vec<i8>,
    /// integer work space
    pub iwork: Vec<HighsInt>,
    /// packed data flag: kept by `pack()`
    pub pack_flag: bool,
    pub pack_count: HighsInt,
    pub pack_index: Vec<HighsInt>,
    pub pack_value: Vec<f64>,
    /// HiGHS: HVectorBase.h `HVectorBase<Real>* next` (line 50): "Allows vectors to be
    /// linked for PAMI". `None` is C++ `0`. Set only by
    /// `HEkkDual::majorUpdateFactor` (HEkkDual.cpp 906-933) and walked by
    /// `HFactor::updateCFT`; see [`HVectorNext`].
    pub next: Option<HVectorNext>,
}

/// HiGHS: HVectorBase.h `next` link. A raw `HVector*`: the PAMI multi-update chains
/// the `col_aq` / `row_ep` vectors of `HEkkDual::multi_finish[]`, which live in one
/// array whose elements never move while the chain is in use. Like the C++ pointer the
/// link may go stale after the chain is consumed; it is only dereferenced inside
/// `HFactor::update_cft` (an `unsafe fn` whose caller guarantees the chain is live).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HVectorNext(pub std::ptr::NonNull<HVector>);

// SAFETY: the link is plain data; it is only dereferenced by `HFactor::update_cft` on the
// thread that built the chain (C++ `HVector*` has the same rule).
unsafe impl Send for HVectorNext {}
// SAFETY: as above; sharing the address does not dereference it.
unsafe impl Sync for HVectorNext {}

/// HiGHS: HVector.h `HVector_ptr` (= `HVector*`). In Rust a pointer
/// parameter is `&mut HVector`; this alias exists only for documentation.
pub type HVectorPtr<'a> = &'a mut HVector;

impl HVector {
    /// HiGHS: HVectorBase.cpp::HVectorBase::setup (lines 22-44)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn setup(&mut self, size_: HighsInt) {
        self.size = size_;
        self.count = 0;
        self.index.resize(size_ as usize, 0);
        self.array.clear();
        self.array.resize(size_ as usize, 0.0);
        self.cwork.clear();
        self.cwork.resize((size_ + 6400) as usize, 0);
        self.iwork.clear();
        self.iwork.resize((size_ * 4) as usize, 0);

        self.pack_count = 0;
        self.pack_index.resize(size_ as usize, 0);
        self.pack_value.resize(size_ as usize, 0.0);

        self.pack_flag = false;
        self.synthetic_tick = 0.0;
        self.next = None;
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::clear (lines 47-64)
    pub fn clear(&mut self) {
        // Standard HVector to clear
        let dense_clear = self.count < 0 || (self.count as f64) > (self.size as f64) * 0.3;
        if dense_clear {
            // Treat the array as full if there are no indices or too many
            // indices
            let size = self.size as usize;
            self.array.clear();
            self.array.resize(size, 0.0);
        } else {
            // Zero according to the indices of (possible) nonzeros
            // PORT NOTE (perf): the C++ loop `array[index[i]] = 0` over the first
            // `count` indices, written as a slice walk so only the scatter target is
            // bounds-checked. Same stores in the same order.
            let count = self.count.max(0) as usize;
            let array = &mut self.array;
            for &idx in &self.index[..count] {
                array[idx as usize] = 0.0;
            }
        }
        self.clear_scalars();
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::clearScalars (lines 67-79)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn clear_scalars(&mut self) {
        // Reset the flag to indicate when to pack
        self.pack_flag = false;
        // Zero the number of stored indices
        self.count = 0;
        // Zero the synthetic clock for operations with this vector
        self.synthetic_tick = 0.0;
        // Initialise the next value
        self.next = None;
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::tight (lines 82-104)
    pub fn tight(&mut self) {
        if self.count < 0 {
            for val in self.array.iter_mut() {
                if val.abs() < kHighsTiny {
                    *val = 0.0;
                }
            }
        } else {
            // PORT NOTE (perf): the C++ loop over the first `count` indices, with the
            // index list and the value array borrowed as disjoint fields and `index`
            // sliced to `count`. The in-place compaction `index[totalCount++]` writes
            // at or below the read position, as in C++.
            let count = self.count.max(0) as usize;
            let index = &mut self.index[..count];
            let array = &mut self.array;
            let mut total = 0usize;
            for i in 0..count {
                let my_index = index[i];
                let value = array[my_index as usize];
                if value.abs() >= kHighsTiny {
                    index[total] = my_index;
                    total += 1;
                } else {
                    array[my_index as usize] = 0.0;
                }
            }
            self.count = total as HighsInt;
        }
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::pack (lines 107-121)
    pub fn pack(&mut self) {
        if !self.pack_flag {
            return;
        }
        self.pack_flag = false;
        self.pack_count = 0;
        for i in 0..(self.count.max(0) as usize) {
            let ipack = self.index[i];
            self.pack_index[self.pack_count as usize] = ipack;
            self.pack_value[self.pack_count as usize] = self.array[ipack as usize];
            self.pack_count += 1;
        }
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::reIndex (lines 124-134)
    pub fn re_index(&mut self) {
        // Don't do it if there are relatively few nonzeros
        if self.count >= 0 && (self.count as f64) <= (self.size as f64) * 0.1 {
            return;
        }
        // PORT NOTE (perf): the C++ dense scan over `array[0..size)`, written as a slice
        // walk so the per-element read needs no bounds check. Same order, same stores.
        let size = self.size as usize;
        let index = &mut self.index;
        let mut count = 0usize;
        for (i, &v) in self.array[..size].iter().enumerate() {
            // (double)array[i]: truthy iff nonzero (NaN is truthy)
            if v != 0.0 {
                index[count] = i as HighsInt;
                count += 1;
            }
        }
        self.count = count as HighsInt;
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::copy (lines 138-155), `FromReal = double`
    /// MUST-PORT (determinism, synthetic ticks drive reinversion): copies `synthetic_tick`.
    pub fn copy(&mut self, from: &HVector) {
        self.clear();
        self.synthetic_tick = from.synthetic_tick;
        let from_count = from.count;
        self.count = from_count;
        for i in 0..(from_count.max(0) as usize) {
            let i_from = from.index[i];
            let x_from = from.array[i_from as usize];
            self.index[i] = i_from;
            self.array[i_from as usize] = x_from;
        }
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::norm2 (lines 157-171)
    pub fn norm2(&self) -> f64 {
        let work_count = self.count;
        let mut result = 0.0_f64;
        for i in 0..(work_count.max(0) as usize) {
            let value = self.array[self.index[i] as usize];
            result += value * value;
        }
        result
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::saxpy (lines 176-200), `RealPivX = RealPiv = double`
    pub fn saxpy(&mut self, pivot_x: f64, pivot: &HVector) {
        let mut work_count = self.count;
        let pivot_count = pivot.count;
        for k in 0..(pivot_count.max(0) as usize) {
            let i_row = pivot.index[k] as usize;
            let x0 = self.array[i_row];
            let x1 = x0 + pivot_x * pivot.array[i_row];
            if x0 == 0.0 {
                self.index[work_count as usize] = i_row as HighsInt;
                work_count += 1;
            }
            self.array[i_row] = if x1.abs() < kHighsTiny {
                kHighsZero
            } else {
                x1
            };
        }
        self.count = work_count;
    }

    /// HiGHS: HVectorBase.cpp::HVectorBase::isEqual (lines 203-213)
    pub fn is_equal(&self, v0: &HVector) -> bool {
        if self.size != v0.size {
            return false;
        }
        if self.count != v0.count {
            return false;
        }
        if self.index != v0.index {
            return false;
        }
        if self.array != v0.array {
            return false;
        }
        if self.synthetic_tick != v0.synthetic_tick {
            return false;
        }
        true
    }
}
