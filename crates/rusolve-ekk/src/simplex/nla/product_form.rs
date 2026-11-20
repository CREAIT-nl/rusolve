//! HiGHS source: simplex/HSimplexNlaProductForm.cpp, simplex/HSimplexNlaFreeze.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Product-form update (used after `HEkk::getIterate` restores a frozen basis) and the
//! freeze/unfreeze helpers.

use super::*;
use crate::prelude::*;

const kProductFormExtraEntries: HighsInt = 1000;
const kProductFormMaxUpdates: HighsInt = 50;
const kProductFormPivotTolerance: f64 = 1e-8;

impl ProductFormUpdate {
    /// HiGHS: HSimplexNlaProductForm.cpp::ProductFormUpdate::clear (lines 22-31, 10 lines)
    pub fn clear(&mut self) {
        self.valid_ = false;
        self.num_row_ = 0;
        self.update_count_ = 0;
        self.pivot_index_.clear();
        self.pivot_value_.clear();
        self.start_.clear();
        self.index_.clear();
        self.value_.clear();
    }

    /// HiGHS: HSimplexNlaProductForm.cpp::ProductFormUpdate::setup (lines 33-44, 12 lines)
    pub fn setup(&mut self, num_row: HighsInt, expected_density: f64) {
        self.valid_ = true;
        self.num_row_ = num_row;
        self.update_count_ = 0;
        self.start_.push(0);
        let reserve_entry_space: HighsInt = (kProductFormExtraEntries as f64
            + kProductFormMaxUpdates as f64 * num_row as f64 * expected_density)
            as HighsInt;
        self.index_.reserve(reserve_entry_space as usize);
        self.value_.reserve(reserve_entry_space as usize);
    }

    /// HiGHS: HSimplexNlaProductForm.cpp::ProductFormUpdate::update (lines 46-64, 19 lines)
    pub fn update(&mut self, aq: &mut HVector, pivot_row: &mut HighsInt) -> HighsInt {
        debug_assert!(0 <= *pivot_row && *pivot_row < self.num_row_);
        if self.update_count_ >= kProductFormMaxUpdates {
            return kRebuildReasonUpdateLimitReached;
        }
        let pivot: f64 = aq.array[*pivot_row as usize];
        if pivot.abs() < kProductFormPivotTolerance {
            return kRebuildReasonPossiblySingularBasis;
        }
        self.pivot_index_.push(*pivot_row);
        self.pivot_value_.push(pivot);
        for i_x in 0..aq.count {
            let i_row: HighsInt = aq.index[i_x as usize];
            if i_row == *pivot_row {
                continue;
            }
            self.index_.push(i_row);
            self.value_.push(aq.array[i_row as usize]);
        }
        self.start_.push(self.index_.len() as HighsInt);
        self.update_count_ += 1;
        kRebuildReasonNo
    }

    /// HiGHS: HSimplexNlaProductForm.cpp::ProductFormUpdate::btran (lines 66-80, 15 lines)
    pub fn btran(&self, rhs: &mut HVector) {
        if !self.valid_ {
            return;
        }
        debug_assert!(rhs.size == self.num_row_);
        debug_assert!(self.start_.len() as HighsInt == self.update_count_ + 1);
        let mut i_x: HighsInt = self.update_count_ - 1;
        while i_x >= 0 {
            let pivot_index: HighsInt = self.pivot_index_[i_x as usize];
            let mut pivot_value: f64 = rhs.array[pivot_index as usize];
            for i_el in self.start_[i_x as usize]..self.start_[(i_x + 1) as usize] {
                pivot_value -=
                    self.value_[i_el as usize] * rhs.array[self.index_[i_el as usize] as usize];
            }
            pivot_value /= self.pivot_value_[i_x as usize];
            if rhs.array[pivot_index as usize] == 0.0 {
                rhs.index[rhs.count as usize] = pivot_index;
                rhs.count += 1;
            }
            rhs.array[pivot_index as usize] = if pivot_value.abs() < kHighsTiny {
                1e-100
            } else {
                pivot_value
            };
            i_x -= 1;
        }
    }

    /// HiGHS: HSimplexNlaProductForm.cpp::ProductFormUpdate::ftran (lines 82-113, 32 lines)
    pub fn ftran(&self, rhs: &mut HVector) {
        if !self.valid_ {
            return;
        }
        debug_assert!(rhs.size == self.num_row_);
        debug_assert!(self.start_.len() as HighsInt == self.update_count_ + 1);
        // Use the zeroed rhs.cwork to record whether a row is in the index
        // list. If RHS fill-in occurs in a row, then we have to add it to
        // the list. We're not tracking cancellation, so we don't need to
        // know where a row appears in the list
        for i_x in 0..rhs.count {
            let r = rhs.index[i_x as usize];
            rhs.cwork[r as usize] = 1;
        }

        for i_x in 0..self.update_count_ {
            let pivot_index: HighsInt = self.pivot_index_[i_x as usize];
            let mut pivot_value: f64 = rhs.array[pivot_index as usize];
            if pivot_value.abs() > kHighsTiny {
                debug_assert!(rhs.cwork[pivot_index as usize] != 0);
                pivot_value /= self.pivot_value_[i_x as usize];
                rhs.array[pivot_index as usize] = pivot_value;
                for i_el in self.start_[i_x as usize]..self.start_[(i_x + 1) as usize] {
                    let i_row: HighsInt = self.index_[i_el as usize];
                    rhs.array[i_row as usize] -= pivot_value * self.value_[i_el as usize];
                    if rhs.cwork[i_row as usize] != 0 {
                        continue;
                    }
                    rhs.cwork[i_row as usize] = 1;
                    rhs.index[rhs.count as usize] = i_row;
                    rhs.count += 1;
                }
            } else {
                rhs.array[pivot_index as usize] = 0.0;
            }
        }
        // Zero the in_index entries used to point into the index list
        for i_x in 0..rhs.count {
            let r = rhs.index[i_x as usize];
            rhs.cwork[r as usize] = 0;
        }
    }
}

impl SimplexIterate {
    /// HiGHS: HSimplexNlaFreeze.cpp::SimplexIterate::clear (lines 17-22, 6 lines)
    pub fn clear(&mut self) {
        self.valid_ = false;
        self.basis_.clear();
        self.invert_.clear();
        self.dual_edge_weight_.clear();
    }
}

impl HSimplexNla {
    /// HiGHS: HSimplexNlaFreeze.cpp::HSimplexNla::putInvert (lines 24-27, 4 lines)
    pub fn put_invert(&mut self) {
        self.simplex_iterate_.valid_ = true;
        self.simplex_iterate_.invert_ = self.factor_.get_invert();
    }

    /// HiGHS: HSimplexNlaFreeze.cpp::HSimplexNla::getInvert (lines 28-28, 1 lines)
    pub fn get_invert(&mut self) {
        self.factor_.set_invert(&self.simplex_iterate_.invert_);
    }
}
