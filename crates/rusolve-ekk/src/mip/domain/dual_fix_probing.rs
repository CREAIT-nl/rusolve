//! HiGHS source: mip/HighsDomain.cpp (lines 641-883) and mip/HighsDomain.h (lines 239-363,
//! 830-836) at the #3244 merge 8a5f9c1d8
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Port of HiGHS PR #3244 (merge 8a5f9c1d8).
//!
//! Borrow story (re-entrant objects and HighsDomain propagators):
//! * The C++ back-pointers `domain` and `mipsolver` of `DualFixProbingPropagation` are gone.
//!   Methods that only read/write the struct's own fields are `impl HighsDomain_DualFixProbingPropagation`
//!   methods. Methods that read the domain or the solver are `impl HighsDomain` methods taking
//!   `ms` (and `&mut self` for the domain).
//! * The C++ `std::function<void()> storeLiftingOpportunity` is
//!   `Option<Arc<dyn Fn(&mut HighsDomain, &mut HighsMipSolver) + Send + Sync>>` (nullptr =
//!   None). The C++ closure (HighsImplications.cpp, `[&]() { storeLiftingOpportunities(col,
//!   val); }`) reads the global domain and the solver through references; with the
//!   back-pointers gone, the domain and `ms` are passed at the call instead.

use crate::prelude::*;

/// C++ `std::vector::assign(n, val)`: clear, then fill with `n` copies of `val`.
trait CppAssign<T> {
    fn assign(&mut self, n: usize, val: T);
}

impl<T: Clone> CppAssign<T> for Vec<T> {
    fn assign(&mut self, n: usize, val: T) {
        self.clear();
        self.resize(n, val);
    }
}

/// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::DualFixProbingFixDirection (lines 256-260)
/// C++: `enum DualFixProbingFixDirection { FixUndecided, FixLowerBound, FixUpperBound }`
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HighsDomain_DualFixProbingFixDirection {
    #[default]
    FixUndecided,
    FixLowerBound,
    FixUpperBound,
}

/// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::RowSide (lines 243-247)
/// C++: `struct RowSide { HighsInt row; bool isRhs; operator size_t() const { return 2 * row + isRhs; } }`
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HighsDomain_DualFixProbingRowSide {
    pub row: HighsInt,
    pub is_rhs: bool,
}

impl HighsDomain_DualFixProbingRowSide {
    /// C++: `operator size_t() const { return 2 * row + isRhs; }`
    pub fn index(&self) -> usize {
        2 * self.row as usize + self.is_rhs as usize
    }
}

/// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::FixedZeroCostColumn (lines 262-265)
/// C++: `struct FixedZeroCostColumn { HighsInt col; DualFixProbingFixDirection direction; }`
#[derive(Clone, Copy, Debug, Default)]
pub struct HighsDomain_DualFixProbingFixedZeroCostColumn {
    pub col: HighsInt,
    pub direction: HighsDomain_DualFixProbingFixDirection,
}

/// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation (lines 239-363)
/// C++: `struct DualFixProbingPropagation { ... }`
///
/// Removed C++ members:
/// * `HighsDomain* domain`: removed; methods take the owning domain as `self`.
/// * `HighsMipSolver* mipsolver`: removed; methods take `ms`.
/// * `std::function<void()> storeLiftingOpportunity`: `Option<Arc<dyn Fn(&mut HighsDomain,
///   &mut HighsMipSolver) + Send + Sync>>` (the closure gets the domain and `ms` at the call).
#[derive(Clone)]
pub struct HighsDomain_DualFixProbingPropagation {
    /// C++: `std::vector<HighsBool> redundantRowFlags_`
    pub redundant_row_flags_: Vec<bool>,
    /// C++: `std::vector<RowSide> redundantRowInds_`
    pub redundant_row_inds_: Vec<HighsDomain_DualFixProbingRowSide>,
    /// C++: `HighsInt previousRedundantRowSize = 0`
    pub previous_redundant_row_size: HighsInt,
    /// C++: `HighsInt numGlobalRedundantRows_ = 0`
    pub num_global_redundant_rows_: HighsInt,
    /// C++: `HighsInt numCachedGlobalRedundantRows_ = 0`
    pub num_cached_global_redundant_rows_: HighsInt,
    /// C++: `std::vector<DualFixProbingFixDirection> zeroCostDirections_`
    pub zero_cost_directions_: Vec<HighsDomain_DualFixProbingFixDirection>,
    /// C++: `std::vector<FixedZeroCostColumn> fixedZeroCostColumns_`
    pub fixed_zero_cost_columns_: Vec<HighsDomain_DualFixProbingFixedZeroCostColumn>,
    /// C++: `bool applyingZeroCostFixings_ = false`
    pub applying_zero_cost_fixings_: bool,
    /// C++: `HighsInt zeroCostStartPos_ = kHighsIInf`
    pub zero_cost_start_pos_: HighsInt,
    /// C++: `std::function<void()> storeLiftingOpportunity`
    pub store_lifting_opportunity:
        Option<std::sync::Arc<dyn Fn(&mut HighsDomain, &mut HighsMipSolver) + Send + Sync>>,
    /// C++: `bool enabled_ = false`
    pub enabled_: bool,
    /// C++: `std::vector<HighsInt> colLowerLocksOriginal_`
    pub col_lower_locks_original_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> colUpperLocksOriginal_`
    pub col_upper_locks_original_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> colLowerReducedNumLocks_`
    pub col_lower_reduced_num_locks_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> colUpperReducedNumLocks_`
    pub col_upper_reduced_num_locks_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> clearColNumReducedLocks_`
    pub clear_col_num_reduced_locks_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> candidateFixedCols_`
    pub candidate_fixed_cols_: Vec<HighsInt>,
    /// C++: `std::vector<HighsBool> candidateColFixedFlags_`
    pub candidate_col_fixed_flags_: Vec<bool>,
    /// C++: `std::vector<HighsInt> globalCandidateFixedCols_`
    pub global_candidate_fixed_cols_: Vec<HighsInt>,
}

impl Default for HighsDomain_DualFixProbingPropagation {
    fn default() -> Self {
        Self {
            redundant_row_flags_: Default::default(),
            redundant_row_inds_: Default::default(),
            previous_redundant_row_size: 0,
            num_global_redundant_rows_: 0,
            num_cached_global_redundant_rows_: 0,
            zero_cost_directions_: Default::default(),
            fixed_zero_cost_columns_: Default::default(),
            applying_zero_cost_fixings_: false,
            zero_cost_start_pos_: kHighsIInf,
            store_lifting_opportunity: None,
            enabled_: false,
            col_lower_locks_original_: Default::default(),
            col_upper_locks_original_: Default::default(),
            col_lower_reduced_num_locks_: Default::default(),
            col_upper_reduced_num_locks_: Default::default(),
            clear_col_num_reduced_locks_: Default::default(),
            candidate_fixed_cols_: Default::default(),
            candidate_col_fixed_flags_: Default::default(),
            global_candidate_fixed_cols_: Default::default(),
        }
    }
}

impl std::fmt::Debug for HighsDomain_DualFixProbingPropagation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HighsDomain_DualFixProbingPropagation")
            .field("enabled_", &self.enabled_)
            .field(
                "previous_redundant_row_size",
                &self.previous_redundant_row_size,
            )
            .field(
                "num_global_redundant_rows_",
                &self.num_global_redundant_rows_,
            )
            .field(
                "num_cached_global_redundant_rows_",
                &self.num_cached_global_redundant_rows_,
            )
            .field(
                "applying_zero_cost_fixings_",
                &self.applying_zero_cost_fixings_,
            )
            .field("zero_cost_start_pos_", &self.zero_cost_start_pos_)
            .finish_non_exhaustive()
    }
}

impl HighsDomain_DualFixProbingPropagation {
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::setEnabled (line 288)
    /// C++: `void setEnabled(const bool val) { enabled_ = val; }`
    pub fn set_enabled(&mut self, val: bool) {
        self.enabled_ = val;
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::isEnabled (line 290)
    /// C++: `bool isEnabled() const { return enabled_; }`
    pub fn is_enabled(&self) -> bool {
        self.enabled_
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::isActive (lines 293-297)
    /// C++: `bool isActive() const { return enabled_ && (!candidateFixedCols_.empty() || static_cast<HighsInt>(redundantRowInds_.size()) > previousRedundantRowSize); }`
    pub fn is_active(&self) -> bool {
        self.enabled_
            && (!self.candidate_fixed_cols_.is_empty()
                || self.redundant_row_inds_.len() as HighsInt > self.previous_redundant_row_size)
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::isZeroCostFixingActive (line 299)
    /// C++: `bool isZeroCostFixingActive() const { return applyingZeroCostFixings_; }`
    pub fn is_zero_cost_fixing_active(&self) -> bool {
        self.applying_zero_cost_fixings_
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::setZeroCostFixingPosition (line 303)
    /// C++: `void setZeroCostFixingPosition(HighsInt v) { zeroCostStartPos_ = v; }`
    pub fn set_zero_cost_fixing_position(&mut self, v: HighsInt) {
        self.zero_cost_start_pos_ = v;
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::getZeroCostFixingPosition (line 305)
    /// C++: `HighsInt getZeroCostFixingPosition() const { return zeroCostStartPos_; }`
    pub fn get_zero_cost_fixing_position(&self) -> HighsInt {
        self.zero_cost_start_pos_
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::endProbing (lines 334-348)
    /// C++: `void endProbing() { ... }`
    pub fn end_probing(&mut self) {
        self.set_enabled(false);
        for i in self.num_global_redundant_rows_..self.redundant_row_inds_.len() as HighsInt {
            let idx = self.redundant_row_inds_[i as usize];
            self.redundant_row_flags_[idx.index()] = false;
        }
        self.redundant_row_inds_
            .truncate(self.num_global_redundant_rows_ as usize);
        self.fixed_zero_cost_columns_.clear();
        self.applying_zero_cost_fixings_ = false;
        self.store_lifting_opportunity = None;
        for &col in &self.candidate_fixed_cols_ {
            self.candidate_col_fixed_flags_[col as usize] = false;
        }
        self.candidate_fixed_cols_.clear();
    }
}

impl HighsDomain {
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::HighsDomain::getDualFixProbingPropagation (lines 830-832)
    /// C++: `DualFixProbingPropagation& getDualFixProbingPropagation()`
    pub fn get_dual_fix_probing_propagation(
        &mut self,
    ) -> &mut HighsDomain_DualFixProbingPropagation {
        &mut self.dual_fix_probing_
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::HighsDomain::setDualFixProbingActive (line 834)
    /// C++: `void setDualFixProbingActive(const bool val) { dualFixProbingActive_ = val; }`
    pub fn set_dual_fix_probing_active(&mut self, val: bool) {
        self.dual_fix_probing_active = val;
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::HighsDomain::getDualFixProbingActive (line 836)
    /// C++: `bool getDualFixProbingActive() const { return dualFixProbingActive_; }`
    pub fn get_dual_fix_probing_active(&self) -> bool {
        self.dual_fix_probing_active
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::ableToFixToLb (lines 307-310)
    /// C++: `bool ableToFixToLb(const HighsInt col) const { return domain->col_lower_[col] != -kHighsInf && mipsolver->model_->col_cost_[col] >= 0; }`
    pub fn dfp_able_to_fix_to_lb(&self, ms: &HighsMipSolver, col: HighsInt) -> bool {
        self.col_lower_[col as usize] != -kHighsInf && mip_model!(ms).col_cost_[col as usize] >= 0.0
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::ableToFixToUb (lines 312-315)
    /// C++: `bool ableToFixToUb(const HighsInt col) const { return domain->col_upper_[col] != kHighsInf && mipsolver->model_->col_cost_[col] <= 0; }`
    pub fn dfp_able_to_fix_to_ub(&self, ms: &HighsMipSolver, col: HighsInt) -> bool {
        self.col_upper_[col as usize] != kHighsInf && mip_model!(ms).col_cost_[col as usize] <= 0.0
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.h::DualFixProbingPropagation::beginProbing (lines 317-332)
    /// C++: `void beginProbing() { ... }`
    pub fn dfp_begin_probing(&mut self, ms: &HighsMipSolver) {
        self.dual_fix_probing_.previous_redundant_row_size =
            self.dual_fix_probing_.num_global_redundant_rows_;
        debug_assert!(
            self.dual_fix_probing_.redundant_row_inds_.len() as HighsInt
                == self.dual_fix_probing_.num_global_redundant_rows_
        );
        self.dual_fix_probing_.fixed_zero_cost_columns_.clear();
        self.dual_fix_probing_
            .set_zero_cost_fixing_position(kHighsIInf);
        self.dual_fix_probing_.applying_zero_cost_fixings_ = false;
        self.dual_fix_probing_.set_enabled(true);

        for &col in &self.dual_fix_probing_.clear_col_num_reduced_locks_ {
            self.dual_fix_probing_.col_lower_reduced_num_locks_[col as usize] = 0;
            self.dual_fix_probing_.col_upper_reduced_num_locks_[col as usize] = 0;
        }
        self.dual_fix_probing_.clear_col_num_reduced_locks_.clear();
        self.dfp_cache_global_redundant_rows(ms);
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.cpp::DualFixProbingPropagation::recomputeLocks (lines 641-689)
    /// C++: `void recomputeLocks() { ... }`
    pub fn dfp_recompute_locks(&mut self, ms: &HighsMipSolver) {
        // mipsolver = domain->mipsolver; (dropped: ms is passed)
        self.dual_fix_probing_
            .redundant_row_flags_
            .assign(2 * ms.num_row() as usize, false);
        self.dual_fix_probing_.redundant_row_inds_.clear();
        self.dual_fix_probing_
            .redundant_row_inds_
            .reserve(2 * ms.num_row() as usize);
        self.dual_fix_probing_.zero_cost_directions_.assign(
            ms.num_col() as usize,
            HighsDomain_DualFixProbingFixDirection::FixUndecided,
        );
        self.dual_fix_probing_.fixed_zero_cost_columns_.clear();
        self.dual_fix_probing_
            .fixed_zero_cost_columns_
            .reserve(ms.num_col() as usize);

        self.dual_fix_probing_.applying_zero_cost_fixings_ = false;
        self.dual_fix_probing_.previous_redundant_row_size = 0;
        self.dual_fix_probing_.num_global_redundant_rows_ = 0;
        self.dual_fix_probing_.num_cached_global_redundant_rows_ = 0;

        self.dual_fix_probing_
            .col_lower_locks_original_
            .assign(ms.num_col() as usize, 0);
        self.dual_fix_probing_
            .col_upper_locks_original_
            .assign(ms.num_col() as usize, 0);
        self.dual_fix_probing_
            .col_lower_reduced_num_locks_
            .assign(ms.num_col() as usize, 0);
        self.dual_fix_probing_
            .col_upper_reduced_num_locks_
            .assign(ms.num_col() as usize, 0);

        self.dual_fix_probing_.candidate_fixed_cols_.clear();
        self.dual_fix_probing_
            .candidate_fixed_cols_
            .reserve(ms.num_col() as usize);
        self.dual_fix_probing_
            .candidate_col_fixed_flags_
            .assign(ms.num_col() as usize, false);
        self.dual_fix_probing_.clear_col_num_reduced_locks_.clear();
        self.dual_fix_probing_
            .clear_col_num_reduced_locks_
            .reserve(ms.num_col() as usize);
        self.dual_fix_probing_.global_candidate_fixed_cols_.clear();
        self.dual_fix_probing_
            .global_candidate_fixed_cols_
            .reserve(ms.num_col() as usize);

        // compute the locks for each variable
        let model = mip_model!(ms);
        for col in 0..model.a_matrix_.num_col_ {
            for k in
                model.a_matrix_.start_[col as usize]..model.a_matrix_.start_[(col + 1) as usize]
            {
                let row = model.a_matrix_.index_[k as usize];
                let val = model.a_matrix_.value_[k as usize];
                let lhs = model.row_lower_[row as usize];
                let rhs = model.row_upper_[row as usize];
                if (val > 0.0 && rhs != kHighsInf) || (val < 0.0 && lhs != -kHighsInf) {
                    self.dual_fix_probing_.col_upper_locks_original_[col as usize] += 1;
                }
                if (val > 0.0 && lhs != -kHighsInf) || (val < 0.0 && rhs != kHighsInf) {
                    self.dual_fix_probing_.col_lower_locks_original_[col as usize] += 1;
                }
            }
        }

        // Check for any pre-redundant rows
        for row in 0..ms.num_row() {
            self.dfp_update_rhs_redundant(ms, row);
            self.dfp_update_lhs_redundant(ms, row);
        }
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.cpp::DualFixProbingPropagation::cacheGlobalRedundantRows (lines 691-725)
    /// C++: `void cacheGlobalRedundantRows() { ... }`
    pub fn dfp_cache_global_redundant_rows(&mut self, ms: &HighsMipSolver) {
        // Go through rows made redundant via domain changes made while
        // dualFixProbingPropagation was not enabled
        while self.dual_fix_probing_.num_cached_global_redundant_rows_
            < self.dual_fix_probing_.num_global_redundant_rows_
        {
            let side = self.dual_fix_probing_.redundant_row_inds_
                [self.dual_fix_probing_.num_cached_global_redundant_rows_ as usize];
            let start = ms.mipdata_.a_rstart_[side.row as usize];
            let end = ms.mipdata_.a_rstart_[(side.row + 1) as usize];
            for i in start..end {
                let val = ms.mipdata_.a_rvalue_[i as usize];
                if val == 0.0 {
                    continue;
                }
                let col = ms.mipdata_.a_rindex_[i as usize];
                let is_upper_lock = (val > 0.0) == side.is_rhs;
                let locks_ref = if is_upper_lock {
                    &mut self.dual_fix_probing_.col_upper_locks_original_[col as usize]
                } else {
                    &mut self.dual_fix_probing_.col_lower_locks_original_[col as usize]
                };
                *locks_ref -= 1;
                let locks = *locks_ref;
                let cost = mip_model!(ms).col_cost_[col as usize];
                if locks == 0 && ((is_upper_lock && cost <= 0.0) || (!is_upper_lock && cost >= 0.0))
                {
                    self.dual_fix_probing_
                        .global_candidate_fixed_cols_
                        .push(col);
                }
            }
            self.dual_fix_probing_.num_cached_global_redundant_rows_ += 1;
        }

        // Schedule candidate fixes. These will naturally get fixed in the
        // global domain over time.
        let mut num_candidates: HighsInt = 0;
        for i in 0..self.dual_fix_probing_.global_candidate_fixed_cols_.len() {
            let col = self.dual_fix_probing_.global_candidate_fixed_cols_[i];
            if self.is_fixed(col) || self.dual_fix_probing_.candidate_col_fixed_flags_[col as usize]
            {
                continue;
            }
            self.dual_fix_probing_.global_candidate_fixed_cols_[num_candidates as usize] = col;
            num_candidates += 1;
            self.dual_fix_probing_.candidate_col_fixed_flags_[col as usize] = true;
            self.dual_fix_probing_.candidate_fixed_cols_.push(col);
        }
        self.dual_fix_probing_
            .global_candidate_fixed_cols_
            .truncate(num_candidates as usize);
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.cpp::DualFixProbingPropagation::updateRhsRedundant (lines 727-739)
    /// HiGHS 1.16 (#3575): HighsDomain.cpp::DualFixProbingPropagation::updateRhsRedundant (lines 732-743)
    /// C++: `void updateRhsRedundant(HighsInt row) { ... }`
    pub fn dfp_update_rhs_redundant(&mut self, ms: &HighsMipSolver, row: HighsInt) {
        let idx = HighsDomain_DualFixProbingRowSide { row, is_rhs: true };
        if self.activitymaxinf_[row as usize] != 0
            || self.dual_fix_probing_.redundant_row_flags_[idx.index()]
            || mip_model!(ms).row_upper_[row as usize] == kHighsInf
        {
            return;
        }

        // HiGHS 1.16 #3575 (h16_row_redundant): HighsDomain.cpp lines 738 (1.16.0)
        let h16_row_redundant = ms.options_mip_.h16_row_redundant;
        let redundant = if h16_row_redundant {
            self.is_row_upper_redundant(ms, row)
        } else {
            self.get_max_activity(row)
                <= mip_model!(ms).row_upper_[row as usize] + ms.mipdata_.feastol
        };
        if redundant {
            self.dual_fix_probing_.redundant_row_inds_.push(idx);
            self.dual_fix_probing_.redundant_row_flags_[idx.index()] = true;
            if !self.dual_fix_probing_.is_enabled() {
                self.dual_fix_probing_.num_global_redundant_rows_ += 1;
            }
        }
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.cpp::DualFixProbingPropagation::updateLhsRedundant (lines 741-753)
    /// HiGHS 1.16 (#3575): HighsDomain.cpp::DualFixProbingPropagation::updateLhsRedundant (lines 745-756)
    /// C++: `void updateLhsRedundant(HighsInt row) { ... }`
    pub fn dfp_update_lhs_redundant(&mut self, ms: &HighsMipSolver, row: HighsInt) {
        let idx = HighsDomain_DualFixProbingRowSide { row, is_rhs: false };
        if self.activitymininf_[row as usize] != 0
            || self.dual_fix_probing_.redundant_row_flags_[idx.index()]
            || mip_model!(ms).row_lower_[row as usize] == -kHighsInf
        {
            return;
        }

        // HiGHS 1.16 #3575 (h16_row_redundant): HighsDomain.cpp lines 751 (1.16.0)
        let h16_row_redundant = ms.options_mip_.h16_row_redundant;
        let redundant = if h16_row_redundant {
            self.is_row_lower_redundant(ms, row)
        } else {
            self.get_min_activity(row)
                >= mip_model!(ms).row_lower_[row as usize] - ms.mipdata_.feastol
        };
        if redundant {
            self.dual_fix_probing_.redundant_row_inds_.push(idx);
            self.dual_fix_probing_.redundant_row_flags_[idx.index()] = true;
            if !self.dual_fix_probing_.is_enabled() {
                self.dual_fix_probing_.num_global_redundant_rows_ += 1;
            }
        }
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.cpp::DualFixProbingPropagation::propagate (lines 755-847)
    /// C++: `void propagate() { ... }`
    pub fn dfp_propagate(&mut self, ms: &mut HighsMipSolver) {
        if !self.dual_fix_probing_.is_active() {
            return;
        }

        for i in self.dual_fix_probing_.previous_redundant_row_size
            ..self.dual_fix_probing_.redundant_row_inds_.len() as HighsInt
        {
            let idx = self.dual_fix_probing_.redundant_row_inds_[i as usize];
            let row = idx.row;
            let is_lhs = !idx.is_rhs;
            let start = ms.mipdata_.a_rstart_[row as usize];
            let end = ms.mipdata_.a_rstart_[(row + 1) as usize];
            for j in start..end {
                let col = ms.mipdata_.a_rindex_[j as usize];
                if self.is_fixed(col) {
                    continue;
                }
                let val = ms.mipdata_.a_rvalue_[j as usize];
                let cost = mip_model!(ms).col_cost_[col as usize];
                if val != 0.0 {
                    // if LHS: Positive val removes a lower lock, negative an upper
                    // if RHS: Negative val removes a lower lock, positive an upper
                    let reduces_lower_lock = (val > 0.0) == is_lhs;
                    if reduces_lower_lock && cost >= 0.0 {
                        if self.dual_fix_probing_.col_lower_reduced_num_locks_[col as usize] == 0
                            && self.dual_fix_probing_.col_upper_reduced_num_locks_[col as usize]
                                == 0
                        {
                            self.dual_fix_probing_
                                .clear_col_num_reduced_locks_
                                .push(col);
                        }
                        self.dual_fix_probing_.col_lower_reduced_num_locks_[col as usize] += 1;
                        if self.dual_fix_probing_.col_lower_reduced_num_locks_[col as usize]
                            == self.dual_fix_probing_.col_lower_locks_original_[col as usize]
                        {
                            self.dfp_add_candidate_fixing(col);
                        }
                    } else if !reduces_lower_lock && cost <= 0.0 {
                        if self.dual_fix_probing_.col_lower_reduced_num_locks_[col as usize] == 0
                            && self.dual_fix_probing_.col_upper_reduced_num_locks_[col as usize]
                                == 0
                        {
                            self.dual_fix_probing_
                                .clear_col_num_reduced_locks_
                                .push(col);
                        }
                        self.dual_fix_probing_.col_upper_reduced_num_locks_[col as usize] += 1;
                        if self.dual_fix_probing_.col_upper_reduced_num_locks_[col as usize]
                            == self.dual_fix_probing_.col_upper_locks_original_[col as usize]
                        {
                            self.dfp_add_candidate_fixing(col);
                        }
                    }
                }
            }
        }

        self.dual_fix_probing_.previous_redundant_row_size =
            self.dual_fix_probing_.redundant_row_inds_.len() as HighsInt;

        // C++ range-for over candidateFixedCols_; the body never pushes to it.
        for ci in 0..self.dual_fix_probing_.candidate_fixed_cols_.len() {
            let col = self.dual_fix_probing_.candidate_fixed_cols_[ci];
            if self.is_fixed(col) {
                continue;
            }
            let can_be_fixed_to_lower = self.dfp_able_to_fix_to_lb(ms, col)
                && self.dual_fix_probing_.col_lower_reduced_num_locks_[col as usize]
                    == self.dual_fix_probing_.col_lower_locks_original_[col as usize];
            let can_be_fixed_to_upper = self.dfp_able_to_fix_to_ub(ms, col)
                && self.dual_fix_probing_.col_upper_reduced_num_locks_[col as usize]
                    == self.dual_fix_probing_.col_upper_locks_original_[col as usize];
            if !can_be_fixed_to_lower && !can_be_fixed_to_upper {
                continue;
            }
            let cost = mip_model!(ms).col_cost_[col as usize];
            if cost == 0.0 {
                let mut direction = self.dual_fix_probing_.zero_cost_directions_[col as usize];
                if direction == HighsDomain_DualFixProbingFixDirection::FixUndecided {
                    direction = if can_be_fixed_to_lower {
                        HighsDomain_DualFixProbingFixDirection::FixLowerBound
                    } else {
                        HighsDomain_DualFixProbingFixDirection::FixUpperBound
                    };
                }
                if direction == HighsDomain_DualFixProbingFixDirection::FixLowerBound
                    && can_be_fixed_to_lower
                {
                    self.dual_fix_probing_.fixed_zero_cost_columns_.push(
                        HighsDomain_DualFixProbingFixedZeroCostColumn {
                            col,
                            direction: HighsDomain_DualFixProbingFixDirection::FixLowerBound,
                        },
                    );
                } else if direction == HighsDomain_DualFixProbingFixDirection::FixUpperBound
                    && can_be_fixed_to_upper
                {
                    self.dual_fix_probing_.fixed_zero_cost_columns_.push(
                        HighsDomain_DualFixProbingFixedZeroCostColumn {
                            col,
                            direction: HighsDomain_DualFixProbingFixDirection::FixUpperBound,
                        },
                    );
                }
            } else {
                if can_be_fixed_to_lower {
                    let boundval = self.col_lower_[col as usize];
                    self.change_bound_2(
                        ms,
                        HighsBoundType::kUpper,
                        col,
                        boundval,
                        HighsDomain_Reason::unspecified(),
                    );
                } else if can_be_fixed_to_upper {
                    let boundval = self.col_upper_[col as usize];
                    self.change_bound_2(
                        ms,
                        HighsBoundType::kLower,
                        col,
                        boundval,
                        HighsDomain_Reason::unspecified(),
                    );
                }
                if self.infeasible() {
                    break;
                }
            }
        }

        for &x in &self.dual_fix_probing_.candidate_fixed_cols_ {
            self.dual_fix_probing_.candidate_col_fixed_flags_[x as usize] = false;
        }
        self.dual_fix_probing_.candidate_fixed_cols_.clear();
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.cpp::DualFixProbingPropagation::propagateZeroCosts (lines 849-883)
    /// C++: `void propagateZeroCosts() { ... }`
    pub fn dfp_propagate_zero_costs(&mut self, ms: &mut HighsMipSolver) {
        if self.dual_fix_probing_.fixed_zero_cost_columns_.is_empty() {
            return;
        }

        if let Some(f) = self.dual_fix_probing_.store_lifting_opportunity.clone() {
            f(self, ms);
            self.dual_fix_probing_.store_lifting_opportunity = None;
        }

        self.dual_fix_probing_.applying_zero_cost_fixings_ = true;
        if self.dual_fix_probing_.zero_cost_start_pos_ == kHighsIInf {
            let pos = self.get_domain_change_stack().len() as HighsInt;
            self.dual_fix_probing_.set_zero_cost_fixing_position(pos);
        }

        // C++ range-for over fixedZeroCostColumns_; the body never pushes to it.
        for fi in 0..self.dual_fix_probing_.fixed_zero_cost_columns_.len() {
            let fixing = self.dual_fix_probing_.fixed_zero_cost_columns_[fi];
            if self.is_fixed(fixing.col) {
                continue;
            }
            debug_assert!(
                self.dual_fix_probing_.zero_cost_directions_[fixing.col as usize]
                    == HighsDomain_DualFixProbingFixDirection::FixUndecided
                    || self.dual_fix_probing_.zero_cost_directions_[fixing.col as usize]
                        == fixing.direction
            );
            if fixing.direction == HighsDomain_DualFixProbingFixDirection::FixLowerBound {
                let boundval = self.col_lower_[fixing.col as usize];
                self.change_bound_2(
                    ms,
                    HighsBoundType::kUpper,
                    fixing.col,
                    boundval,
                    HighsDomain_Reason::unspecified(),
                );
            } else {
                let boundval = self.col_upper_[fixing.col as usize];
                self.change_bound_2(
                    ms,
                    HighsBoundType::kLower,
                    fixing.col,
                    boundval,
                    HighsDomain_Reason::unspecified(),
                );
            }
            if !self.infeasible() {
                self.dual_fix_probing_.zero_cost_directions_[fixing.col as usize] =
                    fixing.direction;
            } else {
                break;
            }
        }

        self.dual_fix_probing_.fixed_zero_cost_columns_.clear();
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsDomain.cpp::DualFixProbingPropagation::propagate (lambda addCandidateFixing, lines 758-763)
    /// C++: `auto addCandidateFixing = [&](const HighsInt col) { ... };`
    fn dfp_add_candidate_fixing(&mut self, col: HighsInt) {
        if !self.dual_fix_probing_.candidate_col_fixed_flags_[col as usize] {
            self.dual_fix_probing_.candidate_fixed_cols_.push(col);
            self.dual_fix_probing_.candidate_col_fixed_flags_[col as usize] = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x1_dfp_d_nonzero_cost() {
        // Probing x=1 makes rows 0,1 redundant, removing all lower locks on y.
        // Since y has positive cost, dual fix probing fixes y=0 during x=1 probe.
        // Probing x=0 propagates y=2 from rows 0,1,2.  Merging gives y=2-2x
        // substitution, incrementing numReductions.
        //
        // Rows:  x + y >= 1    (lower lock on y)
        //       2x + y >= 2    (lower lock on y)
        //      -2x + y <= 2    (upper lock on y)
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 3;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_ = vec![0.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![1.0, 10.0];
        lp.row_lower_ = vec![1.0, 2.0, -kHighsInf];
        lp.row_upper_ = vec![kHighsInf, kHighsInf, 2.0];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kContinuous];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 3, 6];
        lp.a_matrix_.index_ = vec![0, 1, 2, 0, 1, 2];
        lp.a_matrix_.value_ = vec![1.0, 2.0, -2.0, 1.0, 1.0, 1.0];

        let mut options = HighsOptions::default();
        options.mip_dual_fix_probing = true;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.new_data();
        ms.mipdata_.feastol = 1e-6;
        ms.setup_domain_propagation();

        crate::detach!(ms.mipdata_.domains[0], |d| d.dfp_recompute_locks(&ms));
        ms.mipdata_.domains[0].set_dual_fix_probing_active(true);

        let mut num_bound_chgs: HighsInt = 0;
        crate::detach!(ms.mipdata_.implications, |imp| imp.run_probing(
            &mut ms,
            0,
            &mut num_bound_chgs
        ));

        ms.mipdata_.domains[0].set_dual_fix_probing_active(false);
        crate::detach!(ms.mipdata_.domains[0], |d| d.propagate(&mut ms));

        assert!(!ms.mipdata_.domains[0].infeasible());
        // Merging up/down implications yields substitution y = 2 - 2x
        assert!(num_bound_chgs > 0);
        assert_eq!(ms.mipdata_.implications.substitutions.len(), 1);
        assert_eq!(ms.mipdata_.implications.substitutions[0].substcol, 1);
        assert_eq!(ms.mipdata_.implications.substitutions[0].staycol, 0);
        assert_eq!(ms.mipdata_.implications.substitutions[0].offset, 2.0);
        assert_eq!(ms.mipdata_.implications.substitutions[0].scale, -2.0);
    }

    #[test]
    fn x1_dfp_d_tentative_global_fix() {
        // Probing x in either direction forces k=0 via propagation.
        // The tentative clique mechanism should detect this and globally fix k=0.
        //
        // Rows:  x + k + z <= 1   (x=1 forces k+z <= 0 => k=0)
        //       -x + k + z <= 0   (x=0 forces k+z <= 0 => k=0)
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_ = vec![0.0, 0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0, 1.0];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![1.0, 0.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
        ];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        // col 0 (x): row 0 coeff 1, row 1 coeff -1
        // col 1 (k): row 0 coeff 1, row 1 coeff 1
        // col 2 (z): row 0 coeff 1, row 1 coeff 1
        lp.a_matrix_.start_ = vec![0, 2, 4, 6];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![1.0, -1.0, 1.0, 1.0, 1.0, 1.0];

        let mut options = HighsOptions::default();
        options.mip_dual_fix_probing = true;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.new_data();
        ms.mipdata_.feastol = 1e-6;
        ms.setup_domain_propagation();

        crate::detach!(ms.mipdata_.domains[0], |d| d.dfp_recompute_locks(&ms));
        ms.mipdata_.domains[0].set_dual_fix_probing_active(true);

        let mut num_bound_chgs: HighsInt = 0;
        crate::detach!(ms.mipdata_.implications, |imp| imp.run_probing(
            &mut ms,
            0,
            &mut num_bound_chgs
        ));

        ms.mipdata_.domains[0].set_dual_fix_probing_active(false);
        crate::detach!(ms.mipdata_.domains[0], |d| d.propagate(&mut ms));

        assert!(!ms.mipdata_.domains[0].infeasible());
        // k (col 1) globally fixed to 0 by tentative clique mechanism
        assert!(ms.mipdata_.domains[0].is_fixed(1));
        assert_eq!(ms.mipdata_.domains[0].col_lower_[1], 0.0);
        assert_eq!(ms.mipdata_.domains[0].col_upper_[1], 0.0);
    }

    #[test]
    fn x1_dfp_d_tentative_clique_complement() {
        // Complement case: x=1 fixes k=0 via propagation (safe implication),
        // x=0 fixes k=1 via zero-cost lock removal (unsafe, only tentative).
        // The tentative clique mechanism discovers x + k = 1.
        //
        // Cols: x binary (0), k binary (1, cost 0), y continuous [0,10] (2)
        // Rows: x + k <= 1        (upper lock on k; propagation: x=1 => k=0)
        //       k + y >= 0.5      (lower lock on k; not redundant before zero-cost
        //       fix)
        //
        // x=1: propagation from R0 fixes k=0 (safe).
        // x=0: R0 max activity = 0+1 = 1 <= 1 => RHS redundant, upper lock removed.
        //      All upper locks gone, cost=0 => zero-cost fix k=1 (unsafe).
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_ = vec![0.0, 0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0, 10.0];
        lp.row_lower_ = vec![-kHighsInf, 0.5];
        lp.row_upper_ = vec![1.0, kHighsInf];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
        ];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        // col 0 (x): row 0 coeff 1
        // col 1 (k): row 0 coeff 1, row 1 coeff 1
        // col 2 (y): row 1 coeff 1
        lp.a_matrix_.start_ = vec![0, 1, 3, 4];
        lp.a_matrix_.index_ = vec![0, 0, 1, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 1.0];

        let mut options = HighsOptions::default();
        options.mip_dual_fix_probing = true;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.new_data();
        ms.mipdata_.feastol = 1e-6;
        ms.setup_domain_propagation();

        crate::detach!(ms.mipdata_.domains[0], |d| d.dfp_recompute_locks(&ms));
        ms.mipdata_.domains[0].set_dual_fix_probing_active(true);

        let mut num_bound_chgs: HighsInt = 0;
        crate::detach!(ms.mipdata_.implications, |imp| imp.run_probing(
            &mut ms,
            0,
            &mut num_bound_chgs
        ));

        ms.mipdata_.domains[0].set_dual_fix_probing_active(false);
        crate::detach!(ms.mipdata_.domains[0], |d| d.propagate(&mut ms));

        assert!(!ms.mipdata_.domains[0].infeasible());
        // k (col 1) substituted by complement of x: k = 1 - x
        let subst = ms.mipdata_.cliquetable.get_substitution(1);
        assert!(subst.is_some());
        let subst = subst.unwrap();
        assert_eq!(subst.substcol, 1);
        assert_eq!(subst.replace.col(), 0);
        assert_eq!(subst.replace.val(), 0);
    }

    #[test]
    fn x1_dfp_d_tentative_clique_equivalence() {
        // Equivalence case: x=1 fixes k=1 via propagation (safe implication),
        // x=0 fixes k=0 via zero-cost lock removal (unsafe, only tentative).
        // The tentative clique mechanism discovers x = k.
        //
        // Cols: x binary (0), k binary (1, cost 0), y continuous [0,10] (2)
        // Rows: -x + k >= 0       (lower lock on k; propagation: x=1 => k>=1 => k=1)
        //        k + y <= 10.5    (upper lock on k; not redundant before zero-cost
        //        fix)
        //
        // x=1: propagation from R0 fixes k=1 (safe).
        // x=0: R0 min activity = 0+0 = 0 >= 0 => LHS redundant, lower lock removed.
        //      All lower locks gone, cost=0 => zero-cost fix k=0 (unsafe).
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_ = vec![0.0, 0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0, 10.0];
        lp.row_lower_ = vec![0.0, -kHighsInf];
        lp.row_upper_ = vec![kHighsInf, 10.5];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
        ];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        // col 0 (x): row 0 coeff -1
        // col 1 (k): row 0 coeff 1, row 1 coeff 1
        // col 2 (y): row 1 coeff 1
        lp.a_matrix_.start_ = vec![0, 1, 3, 4];
        lp.a_matrix_.index_ = vec![0, 0, 1, 1];
        lp.a_matrix_.value_ = vec![-1.0, 1.0, 1.0, 1.0];

        let mut options = HighsOptions::default();
        options.mip_dual_fix_probing = true;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.new_data();
        ms.mipdata_.feastol = 1e-6;
        ms.setup_domain_propagation();

        crate::detach!(ms.mipdata_.domains[0], |d| d.dfp_recompute_locks(&ms));
        ms.mipdata_.domains[0].set_dual_fix_probing_active(true);

        let mut num_bound_chgs: HighsInt = 0;
        crate::detach!(ms.mipdata_.implications, |imp| imp.run_probing(
            &mut ms,
            0,
            &mut num_bound_chgs
        ));

        ms.mipdata_.domains[0].set_dual_fix_probing_active(false);
        crate::detach!(ms.mipdata_.domains[0], |d| d.propagate(&mut ms));

        assert!(!ms.mipdata_.domains[0].infeasible());
        // k (col 1) substituted by x: k = x
        let subst = ms.mipdata_.cliquetable.get_substitution(1);
        assert!(subst.is_some());
        let subst = subst.unwrap();
        assert_eq!(subst.substcol, 1);
        assert_eq!(subst.replace.col(), 0);
        assert_eq!(subst.replace.val(), 1);
    }

    #[test]
    fn x1_dfp_a_recompute_locks_counts() {
        // Model: 3 columns x0,x1,x2 with col_lower=[0,0,0], col_upper=[1,1,1], col_cost=[0,0,0],
        // integrality continuous; 3 rows, row-wise CSR:
        //   row0: x0 coef 1, x1 coef 1 ; lhs = -inf, rhs = 1
        //   row1: x1 coef 1, x2 coef -1 ; lhs = 0, rhs = +inf
        //   row2: x0 coef 2, x2 coef 1 ; lhs = 1, rhs = 3

        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 3;
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0, 1.0];
        lp.col_cost_ = vec![0.0, 0.0, 0.0];
        lp.integrality_ = vec![
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
        ];
        lp.row_lower_ = vec![-kHighsInf, 0.0, 1.0];
        lp.row_upper_ = vec![1.0, kHighsInf, 3.0];

        // The model matrix is column-wise (C++ `model->a_matrix_` is read by column in
        // recomputeLocks, lines 670-681). Column 0: rows 0, 2; column 1: rows 0, 1;
        // column 2: rows 1, 2.
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 2, 4, 6];
        lp.a_matrix_.index_ = vec![0, 2, 0, 1, 1, 2];
        lp.a_matrix_.value_ = vec![1.0, 2.0, 1.0, 1.0, -1.0, 1.0];

        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);

        // Row-wise AR arrays for the global-redundant-row cache (C++ ARstart_ / ARindex_ /
        // ARvalue_). Row 0: cols 0, 1; row 1: cols 1, 2; row 2: cols 0, 2.
        ms.mipdata_.a_rstart_ = vec![0, 2, 4, 6];
        ms.mipdata_.a_rindex_ = vec![0, 1, 1, 2, 0, 2];
        ms.mipdata_.a_rvalue_ = vec![1.0, 1.0, 1.0, -1.0, 2.0, 1.0];

        let mut dom = HighsDomain::default();
        dom.col_lower_ = vec![0.0; 3];
        dom.col_upper_ = vec![1.0; 3];
        // Activities set to the values of the boxed bounds:
        // row0: min = 0*1+0*1 = 0, max = 1*1+1*1 = 2
        // row1: min = 0*1+1*(-1) = -1, max = 1*1+0*(-1) = 1
        // row2: min = 0*2+0*1 = 0, max = 1*2+1*1 = 3
        dom.activitymin_ = vec![
            HighsCDouble::from(0.0),
            HighsCDouble::from(-1.0),
            HighsCDouble::from(0.0),
        ];
        dom.activitymax_ = vec![
            HighsCDouble::from(2.0),
            HighsCDouble::from(1.0),
            HighsCDouble::from(3.0),
        ];
        dom.activitymininf_ = vec![0, 0, 0];
        dom.activitymaxinf_ = vec![0, 0, 0];

        dom.dfp_recompute_locks(&ms);

        assert_eq!(
            dom.dual_fix_probing_.col_lower_locks_original_,
            vec![1, 1, 1]
        );
        assert_eq!(
            dom.dual_fix_probing_.col_upper_locks_original_,
            vec![2, 1, 2]
        );
        assert_eq!(dom.dual_fix_probing_.num_global_redundant_rows_, 1);
        assert_eq!(
            dom.dual_fix_probing_.redundant_row_inds_,
            vec![HighsDomain_DualFixProbingRowSide {
                row: 2,
                is_rhs: true
            }]
        );

        // beginProbing (HighsDomain.h 317-332) caches the global redundant row 2 (rhs):
        // cacheGlobalRedundantRows (HighsDomain.cpp 691-725) removes one upper lock from
        // col 0 (coef 2) and col 2 (coef 1). No lock count reaches 0, so no global
        // candidate is scheduled.
        dom.dfp_begin_probing(&ms);
        assert!(dom.dual_fix_probing_.is_enabled());
        assert_eq!(dom.dual_fix_probing_.previous_redundant_row_size, 1);
        assert_eq!(dom.dual_fix_probing_.num_cached_global_redundant_rows_, 1);
        assert_eq!(
            dom.dual_fix_probing_.col_upper_locks_original_,
            vec![1, 1, 1]
        );
        assert!(
            dom.dual_fix_probing_
                .global_candidate_fixed_cols_
                .is_empty()
        );
        assert!(!dom.dual_fix_probing_.is_active());

        // endProbing (HighsDomain.h 334-348) keeps the global redundant rows.
        dom.dual_fix_probing_.end_probing();
        assert!(!dom.dual_fix_probing_.is_enabled());
        assert_eq!(dom.dual_fix_probing_.redundant_row_inds_.len(), 1);
        assert!(dom.dual_fix_probing_.redundant_row_flags_[5]);
    }
}

#[cfg(test)]
mod h16_row_redundant_dfp_tests {
    use super::*;

    // One continuous column, one row with feastol 1e-6; the row sides are given by the caller.
    fn dfp_mip(row_lower: f64, row_upper: f64, h16_row_redundant: bool) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![kHighsInf];
        lp.row_lower_ = vec![row_lower];
        lp.row_upper_ = vec![row_upper];
        lp.integrality_ = vec![HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        ms.mipdata_.a_rstart_ = vec![0, 1];
        ms.mipdata_.a_rindex_ = vec![0];
        ms.mipdata_.a_rvalue_ = vec![1.0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.h16_row_redundant = h16_row_redundant;
        ms
    }

    // d1 = round(2.0 + 1e-6) is the double 1.15.1 compares against; d1 - 2.0 is above feastol in
    // CDouble, so the rhs row is redundant in 1.15.1 but not in 1.16.
    #[test]
    fn h16_row_redundant_on_dfp_rhs_not_redundant() {
        let ms = dfp_mip(-kHighsInf, 2.0, true);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymaxinf_ = vec![0];
        dom.activitymax_ = vec![HighsCDouble::from(2.0 + 1e-6)];
        dom.dual_fix_probing_.redundant_row_flags_ = vec![false; 8];
        dom.dfp_update_rhs_redundant(&ms, 0);
        assert!(dom.dual_fix_probing_.redundant_row_inds_.is_empty());
    }

    #[test]
    fn h16_row_redundant_off_dfp_rhs_redundant() {
        let ms = dfp_mip(-kHighsInf, 2.0, false);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymaxinf_ = vec![0];
        dom.activitymax_ = vec![HighsCDouble::from(2.0 + 1e-6)];
        dom.dual_fix_probing_.redundant_row_flags_ = vec![false; 8];
        dom.dfp_update_rhs_redundant(&ms, 0);
        assert_eq!(dom.dual_fix_probing_.redundant_row_inds_.len(), 1);
    }

    // Mirror on the lhs: -2.0 - feastol rounds to -d1 in 1.15.1, while -2.0 - (-d1) is above feastol.
    #[test]
    fn h16_row_redundant_on_dfp_lhs_not_redundant() {
        let ms = dfp_mip(-2.0, kHighsInf, true);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymininf_ = vec![0];
        dom.activitymin_ = vec![HighsCDouble::from(-(2.0 + 1e-6))];
        dom.dual_fix_probing_.redundant_row_flags_ = vec![false; 8];
        dom.dfp_update_lhs_redundant(&ms, 0);
        assert!(dom.dual_fix_probing_.redundant_row_inds_.is_empty());
    }

    #[test]
    fn h16_row_redundant_off_dfp_lhs_redundant() {
        let ms = dfp_mip(-2.0, kHighsInf, false);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymininf_ = vec![0];
        dom.activitymin_ = vec![HighsCDouble::from(-(2.0 + 1e-6))];
        dom.dual_fix_probing_.redundant_row_flags_ = vec![false; 8];
        dom.dfp_update_lhs_redundant(&ms, 0);
        assert_eq!(dom.dual_fix_probing_.redundant_row_inds_.len(), 1);
    }
}
