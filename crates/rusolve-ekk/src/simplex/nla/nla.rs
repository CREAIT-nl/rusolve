//! HiGHS source: simplex/HSimplexNla.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Scaling-aware FTRAN/BTRAN and INVERT wrapper around `HFactor`.

use super::*;
use crate::prelude::*;

impl HSimplexNla {
    /// HiGHS: HSimplexNla.cpp::HSimplexNla::setup (lines 23-41, 19 lines)
    /// `timer`/`analysis` dropped; `lp`/`basic_index`/`options` are used immediately (build
    /// `NlaRefs::new(lp, basic_index, options)` for the `debugCheckData` assertion, a `debug_assert!`).
    /// The `factor_.setupGeneral(...)` call is `self.factor_.setup_general_raw(...)` (no basic_index).
    pub fn setup(
        &mut self,
        lp: &HighsLp,
        basic_index: &[HighsInt],
        options: &HighsOptions,
        factor_a_matrix: &HighsSparseMatrix,
        factor_pivot_threshold: f64,
    ) {
        self.report_ = false;
        self.factor_.setup_general_raw(
            lp.num_col_,
            lp.num_row_,
            lp.num_row_,
            &factor_a_matrix.start_,
            &factor_a_matrix.index_,
            &factor_a_matrix.value_,
            factor_pivot_threshold,
            options.factor_pivot_tolerance,
            options.highs_debug_level,
            &options.log_options,
            true,
            kUpdateMethodFt,
        );
        // rusolve (HiGHS 1.16 backport, default off): the h16_hfactor_pivot flag
        self.factor_.h16_hfactor_pivot = options.h16_hfactor_pivot;
        let refs = NlaRefs::new(lp, basic_index, options);
        debug_assert!(
            self.debug_check_data(&refs, "After HSimplexNla::setup") == HighsDebugStatus::kOk
        );
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::setLpAndScalePointers (lines 43-48, 6 lines)
    /// No-op: the pointers it sets are supplied per call via `NlaRefs::new`.
    pub fn set_lp_and_scale_pointers(&mut self, for_lp: &HighsLp) {}

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::setBasicIndexPointers (lines 50-53, 4 lines)
    /// No-op: `basic_index` is supplied per call (`NlaRefs` / `basic_index` params).
    pub fn set_basic_index_pointers(&mut self, basic_index: &mut [HighsInt]) {}

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::setPointers (lines 55-66, 12 lines)
    /// C++ default arguments: factor_a_matrix = NULL
    /// Only the `factor_a_matrix` branch has an effect (`self.factor_.setup_matrix(m)`);
    /// the other C++ parameters are pointer refreshes with no Rust counterpart.
    pub fn set_pointers(&mut self, for_lp: &HighsLp, factor_a_matrix: Option<&HighsSparseMatrix>) {
        if let Some(m) = factor_a_matrix {
            self.factor_.setup_matrix(m);
        }
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::clear (lines 68-77, 10 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn clear(&mut self) {
        self.report_ = false;
        self.build_synthetic_tick_ = 0.0;
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::invert (lines 79-94, 16 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion): `build_synthetic_tick_ = factor_.build_synthetic_tick`.
    /// `analysis_->analyse_factor_time` is false: `factor_timer_clock_pointer = None`.
    /// `factor_.build()` -> `self.factor_.build(basic_index, None)`.
    pub fn invert(&mut self, basic_index: &mut [HighsInt]) -> HighsInt {
        let factor_timer_clock_pointer: Option<&HighsTimerClock> = None;
        let rank_deficiency = self.factor_.build(basic_index, factor_timer_clock_pointer);
        // Must not have timed out
        debug_assert!(rank_deficiency >= 0);
        self.build_synthetic_tick_ = self.factor_.build_synthetic_tick;
        rank_deficiency
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::btran (lines 96-101, 6 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    pub fn btran(
        &self,
        refs: &NlaRefs,
        rhs: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        self.apply_basis_matrix_col_scale(refs, rhs);
        self.btran_in_scaled_space(rhs, expected_density, factor_timer_clock_pointer);
        self.apply_basis_matrix_row_scale(refs, rhs);
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::ftran (lines 103-108, 6 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    pub fn ftran(
        &self,
        refs: &NlaRefs,
        rhs: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        self.apply_basis_matrix_row_scale(refs, rhs);
        self.ftran_in_scaled_space(rhs, expected_density, factor_timer_clock_pointer);
        self.apply_basis_matrix_col_scale(refs, rhs);
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::btranInScaledSpace (lines 110-114, 5 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    pub fn btran_in_scaled_space(
        &self,
        rhs: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        self.factor_
            .btran_call(rhs, expected_density, factor_timer_clock_pointer);
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::ftranInScaledSpace (lines 116-120, 5 lines)
    /// C++ default arguments: factor_timer_clock_pointer = NULL
    pub fn ftran_in_scaled_space(
        &self,
        rhs: &mut HVector,
        expected_density: f64,
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) {
        self.factor_
            .ftran_call(rhs, expected_density, factor_timer_clock_pointer);
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::update (lines 122-132, 11 lines)
    /// `factor_.update(aq, ep, iRow, hint)` -> `self.factor_.update(basic_index, aq, ep, i_row, hint)`.
    pub fn update(
        &mut self,
        basic_index: &mut [HighsInt],
        aq: &mut HVector,
        ep: &mut HVector,
        i_row: &mut HighsInt,
        hint: &mut HighsInt,
    ) {
        self.report_pack_value("  pack: aq Bf ", aq, false);
        self.report_pack_value("  pack: ep Bf ", ep, false);
        self.factor_.refactor_info_.clear();
        if !self.update_.valid_ {
            self.factor_.update(basic_index, aq, ep, i_row, hint);
        } else {
            *hint = self.update_.update(aq, i_row);
        }
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::update (lines 122-132), PAMI chained form:
    /// `factor_.update(aq, ep, iRow, hint)` with `aq->next` set ->
    /// [`HFactor::update_multi`]; the product-form branch (`update_.valid_`) uses
    /// `aq` / `i_row[0]` as in C++.
    ///
    /// # Safety
    /// As [`HFactor::update_cft`].
    pub unsafe fn update_multi(
        &mut self,
        basic_index: &mut [HighsInt],
        aq: *mut HVector,
        ep: *mut HVector,
        i_row: &mut [HighsInt],
        hint: &mut HighsInt,
    ) {
        // SAFETY: caller guarantees aq/ep are valid HVector pointers (PAMI chain from
        // HEkkDual::majorUpdateFactor); i_row[0] is the first pivot row.
        let aq_ref = unsafe { &mut *aq };
        let ep_ref = unsafe { &mut *ep };
        self.report_pack_value("  pack: aq Bf ", aq_ref, false);
        self.report_pack_value("  pack: ep Bf ", ep_ref, false);
        self.factor_.refactor_info_.clear();
        if !self.update_.valid_ {
            // SAFETY: caller guarantees the HVector chain is live and i_row has at least
            // num_update entries (PAMI protocol).
            unsafe {
                self.factor_.update_multi(basic_index, aq, ep, i_row, hint);
            }
        } else {
            *hint = self.update_.update(aq_ref, &mut i_row[0]);
        }
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::rowEp2NormInScaledSpace (lines 134-184, 51 lines)
    pub fn row_ep2_norm_in_scaled_space(
        &self,
        refs: &NlaRefs,
        i_row: HighsInt,
        row_ep: &HVector,
    ) -> f64 {
        let scale = match refs.scale_ {
            None => return row_ep.norm2(),
            Some(scale) => scale,
        };
        let row_scale = &scale.row;
        // Get the 2-norm of row_ep in the scaled space
        //
        // Determine the scaling that was applied to the unit RHS before
        // scaled BTRAN. This must be unapplied to all components of the
        // result.
        let col_scale_value = self.basic_col_scale_factor(refs, i_row);
        // Now compute the 2-norm of row_ep in the scaled space, unapplying
        // the scaling that was applied after BTRAN
        let mut row_ep_2norm = 0.0;
        let mut to_entry: HighsInt = 0;
        let use_row_indices =
            self.sparse_loop_style(row_ep.count, refs.lp_.num_row_, &mut to_entry);
        for i_entry in 0..to_entry {
            let i_row = if use_row_indices {
                row_ep.index[i_entry as usize]
            } else {
                i_entry
            };
            let value_in_scaled_space =
                row_ep.array[i_row as usize] / (row_scale[i_row as usize] * col_scale_value);
            row_ep_2norm += value_in_scaled_space * value_in_scaled_space;
        }
        row_ep_2norm
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::transformForUpdate (lines 186-221, 36 lines)
    pub fn transform_for_update(
        &mut self,
        refs: &NlaRefs,
        aq: &mut HVector,
        ep: &mut HVector,
        variable_in: HighsInt,
        row_out: HighsInt,
    ) {
        if refs.scale_.is_none() {
            return;
        }
        self.report_pack_value("pack aq Bf ", aq, false);
        let cq_scale_factor = self.variable_scale_factor(refs, variable_in);

        for ix in 0..aq.pack_count {
            aq.pack_value[ix as usize] *= cq_scale_factor;
        }
        self.report_pack_value("pack aq Af ", aq, false);
        //
        // Now focus on the pivot value, aq->array[row_out]
        let pivot_in_scaled_space = self.pivot_in_scaled_space(refs, aq, variable_in, row_out);
        // First scale by cq
        aq.array[row_out as usize] *= cq_scale_factor;
        //
        // Also have to unscale by cp
        let cp_scale_factor = self.basic_col_scale_factor(refs, row_out);
        aq.array[row_out as usize] /= cp_scale_factor;
        debug_assert!(pivot_in_scaled_space == aq.array[row_out as usize]);
        // Hence packValue needs to be unscaled by cp
        for ix in 0..ep.pack_count {
            ep.pack_value[ix as usize] /= cp_scale_factor;
        }
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::variableScaleFactor (lines 223-227, 5 lines)
    pub fn variable_scale_factor(&self, refs: &NlaRefs, i_var: HighsInt) -> f64 {
        let scale = match refs.scale_ {
            None => return 1.0,
            Some(scale) => scale,
        };
        if i_var < refs.lp_.num_col_ {
            scale.col[i_var as usize]
        } else {
            1.0 / scale.row[(i_var - refs.lp_.num_col_) as usize]
        }
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::basicColScaleFactor (lines 229-232, 4 lines)
    pub fn basic_col_scale_factor(&self, refs: &NlaRefs, i_col: HighsInt) -> f64 {
        if refs.scale_.is_none() {
            return 1.0;
        }
        self.variable_scale_factor(refs, refs.basic_index_[i_col as usize])
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::pivotInScaledSpace (lines 234-239, 6 lines)
    pub fn pivot_in_scaled_space(
        &self,
        refs: &NlaRefs,
        aq: &HVector,
        variable_in: HighsInt,
        row_out: HighsInt,
    ) -> f64 {
        aq.array[row_out as usize] * self.variable_scale_factor(refs, variable_in)
            / self.variable_scale_factor(refs, refs.basic_index_[row_out as usize])
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::setPivotThreshold (lines 241-243, 3 lines)
    pub fn set_pivot_threshold(&mut self, new_pivot_threshold: f64) {
        self.factor_.set_pivot_threshold(new_pivot_threshold);
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::applyBasisMatrixRowScale (lines 245-255, 11 lines)
    pub fn apply_basis_matrix_row_scale(&self, refs: &NlaRefs, rhs: &mut HVector) {
        let scale = match refs.scale_ {
            None => return,
            Some(scale) => scale,
        };
        let row_scale = &scale.row;
        let mut to_entry: HighsInt = 0;
        let use_row_indices = self.sparse_loop_style(rhs.count, refs.lp_.num_row_, &mut to_entry);
        for i_entry in 0..to_entry {
            let i_row = if use_row_indices {
                rhs.index[i_entry as usize]
            } else {
                i_entry
            };
            rhs.array[i_row as usize] *= row_scale[i_row as usize];
        }
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::applyBasisMatrixColScale (lines 257-273, 17 lines)
    pub fn apply_basis_matrix_col_scale(&self, refs: &NlaRefs, rhs: &mut HVector) {
        let scale = match refs.scale_ {
            None => return,
            Some(scale) => scale,
        };
        let col_scale = &scale.col;
        let row_scale = &scale.row;
        let mut to_entry: HighsInt = 0;
        let use_row_indices = self.sparse_loop_style(rhs.count, refs.lp_.num_row_, &mut to_entry);
        for i_entry in 0..to_entry {
            let i_col = if use_row_indices {
                rhs.index[i_entry as usize]
            } else {
                i_entry
            };
            let i_var = refs.basic_index_[i_col as usize];
            if i_var < refs.lp_.num_col_ {
                rhs.array[i_col as usize] *= col_scale[i_var as usize];
            } else {
                rhs.array[i_col as usize] /= row_scale[(i_var - refs.lp_.num_col_) as usize];
            }
        }
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::unapplyBasisMatrixRowScale (lines 275-285, 11 lines)
    pub fn unapply_basis_matrix_row_scale(&self, refs: &NlaRefs, rhs: &mut HVector) {
        let scale = match refs.scale_ {
            None => return,
            Some(scale) => scale,
        };
        let row_scale = &scale.row;
        let mut to_entry: HighsInt = 0;
        let use_row_indices = self.sparse_loop_style(rhs.count, refs.lp_.num_row_, &mut to_entry);
        for i_entry in 0..to_entry {
            let i_row = if use_row_indices {
                rhs.index[i_entry as usize]
            } else {
                i_entry
            };
            rhs.array[i_row as usize] /= row_scale[i_row as usize];
        }
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::addCols (lines 287-295, 9 lines)
    /// No-op: C++ only calls `setLpAndScalePointers(updated_lp)`.
    pub fn add_cols(&mut self, updated_lp: &HighsLp) {}

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::addRows (lines 297-307, 11 lines)
    /// `setLpAndScalePointers`/pointer stores have no counterpart; translate the
    /// `factor_.addRows(scaled_ar_matrix)` call as `self.factor_.add_rows(basic_index, scaled_ar_matrix)`.
    pub fn add_rows(
        &mut self,
        updated_lp: &HighsLp,
        basic_index: &mut [HighsInt],
        scaled_ar_matrix: &HighsSparseMatrix,
    ) {
        self.factor_.add_rows(basic_index, scaled_ar_matrix);
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::sparseLoopStyle (lines 309-320, 12 lines)
    pub fn sparse_loop_style(
        &self,
        count: HighsInt,
        dim: HighsInt,
        to_entry: &mut HighsInt,
    ) -> bool {
        // Parameter to decide whether to use just the values in a HVector, or
        // use the indices of their nonzeros
        let use_indices = count >= 0 && (count as f64) < kDensityForIndexing * (dim as f64);
        if use_indices {
            *to_entry = count;
        } else {
            *to_entry = dim;
        }
        use_indices
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::reportArray (lines 322-325, 4 lines)
    /// C++ default arguments: force = false
    /// Reporting: out of scope, no-op.
    pub fn report_array(&self, message: &str, vector: &HVector, force: bool) {}

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::reportArray (lines 327-342, 16 lines)
    /// C++ default arguments: force = false
    /// Reporting: out of scope, no-op.
    pub fn report_array_offset(
        &self,
        message: &str,
        offset: HighsInt,
        vector: &HVector,
        force: bool,
    ) {
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::reportVector (lines 344-363, 20 lines)
    /// Reporting: out of scope, no-op.
    pub fn report_vector(
        &self,
        message: &str,
        num_index: HighsInt,
        vector_value: &[f64],
        vector_index: &[HighsInt],
        force: bool,
    ) {
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::reportArraySparse (lines 365-369, 5 lines)
    /// C++ default arguments: force = false
    /// Reporting: out of scope, no-op.
    pub fn report_array_sparse(&self, message: &str, vector: &HVector, force: bool) {}

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::reportArraySparse (lines 371-402, 32 lines)
    /// C++ default arguments: force = false
    /// Reporting: out of scope, no-op.
    pub fn report_array_sparse_offset(
        &self,
        message: &str,
        offset: HighsInt,
        vector: &HVector,
        force: bool,
    ) {
    }

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::reportPackValue (lines 404-422, 19 lines)
    /// C++ default arguments: force = false
    /// Reporting: out of scope, no-op.
    pub fn report_pack_value(&self, message: &str, vector: &HVector, force: bool) {}

    /// HiGHS: HSimplexNla.cpp::HSimplexNla::debugCheckData (lines 424-517, 94 lines)
    /// HiGHS 1.16 (#3515): HSimplexNla.cpp::HSimplexNla::debugCheckData (lines 424-426)
    /// C++ default arguments: message = ""
    /// With `h16_nla_debug_gate` on, `highs_debug_level < kHighsDebugLevelCheap` returns kOk
    /// before any check (1.16.0 early return). With it off, the check runs at every level,
    /// as in 1.15.1 (HEkk::solve fails if this is not kOk). Translate fully;
    /// C++ `factor_.getAstart()` etc. are the owned copies in `self.factor_`.
    pub fn debug_check_data(&self, refs: &NlaRefs, message: &str) -> HighsDebugStatus {
        let h16_nla_debug_gate = refs.options_.h16_nla_debug_gate;
        if h16_nla_debug_gate && refs.options_.highs_debug_level < kHighsDebugLevelCheap {
            return HighsDebugStatus::kOk;
        }
        let scale_status = if refs.scale_.is_none() {
            "NULL"
        } else {
            "non-NULL"
        };
        // PORT NOTE (perf): C++ copies the whole LP (`HighsLp check_lp = *lp_`) and scales
        // it, but only `check_lp.a_matrix_` is read. Scaling leaves `start_`/`index_`
        // unchanged, so those are compared in place, and only the matrix values are
        // copied and scaled with the `HighsSparseMatrix::applyScale` expression
        // (`v *= col*row`), so every comparison sees the values C++ compares.
        let lp_matrix = &refs.lp_.a_matrix_;
        let scaled_value: Vec<f64>;
        let check_value: &[f64] = match refs.scale_ {
            Some(scale) if lp_matrix.is_colwise() => {
                let mut value = lp_matrix.value_.clone();
                for i_col in 0..lp_matrix.num_col_ as usize {
                    let (begin, end) = (
                        lp_matrix.start_[i_col] as usize,
                        lp_matrix.start_[i_col + 1] as usize,
                    );
                    for i_el in begin..end {
                        let i_row = lp_matrix.index_[i_el] as usize;
                        value[i_el] *= scale.col[i_col] * scale.row[i_row];
                    }
                }
                scaled_value = value;
                &scaled_value
            }
            Some(scale) => {
                let mut m = lp_matrix.clone();
                m.apply_scale(scale);
                scaled_value = m.value_;
                &scaled_value
            }
            None => &lp_matrix.value_,
        };
        let factor_a_start = self.factor_.get_astart();
        let factor_a_index = self.factor_.get_aindex();
        let factor_a_value = self.factor_.get_avalue();

        // The C++ pointer-identity checks (factor_ matrix pointers vs lp_ matrix pointers when
        // scale_ == NULL) have no counterpart: HFactor owns copies, so only the content checks
        // below apply.
        let mut error_col: HighsInt = -1;
        for i_col in 0..refs.lp_.num_col_ + 1 {
            if lp_matrix.start_[i_col as usize] != factor_a_start[i_col as usize] {
                error_col = i_col;
                break;
            }
        }
        let mut error_found = error_col >= 0;
        if error_found {
            crate::highs_log_user!(
                &refs.options_.log_options,
                HighsLogType::kError,
                "CheckNlaData: ({}) scale_ is {} check_lp.a_matrix_.start_ != factor_Astart for col {} ({} != {})\n",
                message,
                scale_status,
                error_col,
                lp_matrix.start_[error_col as usize],
                factor_a_start[error_col as usize]
            );
            debug_assert!(!error_found);
            return HighsDebugStatus::kLogicalError;
        }
        let nnz = lp_matrix.num_nz();
        let mut error_el: HighsInt = -1;
        for i_el in 0..nnz {
            if lp_matrix.index_[i_el as usize] != factor_a_index[i_el as usize] {
                error_el = i_el;
                break;
            }
        }
        error_found = error_el >= 0;
        if error_found {
            crate::highs_log_user!(
                &refs.options_.log_options,
                HighsLogType::kError,
                "CheckNlaData: ({}) scale_ is {} check_lp.a_matrix_.index_ != factor_Aindex for el {} ({} != {})\n",
                message,
                scale_status,
                error_el,
                lp_matrix.index_[error_el as usize],
                factor_a_index[error_el as usize]
            );
            debug_assert!(!error_found);
            return HighsDebugStatus::kLogicalError;
        }
        for i_el in 0..nnz {
            if check_value[i_el as usize] != factor_a_value[i_el as usize] {
                error_el = i_el;
                break;
            }
        }
        error_found = error_el >= 0;
        if error_found {
            crate::highs_log_user!(
                &refs.options_.log_options,
                HighsLogType::kError,
                "CheckNlaData: ({}) scale_ is {} check_lp.a_matrix_.value_ != factor_Avalue for el {} ({} != {})\n",
                message,
                scale_status,
                error_el,
                check_value[error_el as usize],
                factor_a_value[error_el as usize]
            );
            debug_assert!(!error_found);
            return HighsDebugStatus::kLogicalError;
        }
        HighsDebugStatus::kOk
    }
}

#[cfg(test)]
mod h16_nla_debug_gate_tests {
    use super::*;
    use crate::lp_data::highs_options::HighsOptions;
    use crate::lp_data::highs_run::Highs;
    use crate::prelude::*;

    /// Same LP as the `h16_*` option test fixture: optimum 6 - 8/7 + 0.5.
    fn tiny_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 3,
            col_cost_: vec![3.0, -4.0],
            col_lower_: vec![0.0, -kHighsInf],
            col_upper_: vec![10.0, kHighsInf],
            row_lower_: vec![1.0, -kHighsInf, 2.0],
            row_upper_: vec![4.0, kHighsInf, 2.0],
            sense_: ObjSense::kMaximize,
            offset_: 0.5,
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 1, 2];
        lp.a_matrix_.value_ = vec![2.0, 3.0, 5.0, 7.0];
        lp
    }

    /// The factor holds the copy taken at setup; the returned LP differs from it in one
    /// value, so the check reports kLogicalError unless the 1.16.0 early return is taken.
    fn nla_with_bad_lp() -> (HSimplexNla, HighsLp, Vec<HighsInt>) {
        let lp = tiny_lp();
        let basic_index: Vec<HighsInt> = vec![2, 3, 4];
        let mut options = HighsOptions::default();
        options.highs_debug_level = kHighsDebugLevelNone;
        let mut nla = HSimplexNla::default();
        nla.setup(&lp, &basic_index, &options, &lp.a_matrix_, 0.1);
        let mut bad = lp;
        bad.a_matrix_.value_[0] = 99.0;
        (nla, bad, basic_index)
    }

    fn check_bad(gate: bool, debug_level: HighsInt) -> HighsDebugStatus {
        let (nla, bad, basic_index) = nla_with_bad_lp();
        let mut options = HighsOptions::default();
        options.highs_debug_level = debug_level;
        options.h16_nla_debug_gate = gate;
        let refs = NlaRefs::new(&bad, &basic_index, &options);
        nla.debug_check_data(&refs, "gate test")
    }

    /// Flag on at level 0 (below kHighsDebugLevelCheap): returns kOk without checking.
    #[test]
    fn h16_nla_debug_gate_skips_check_at_level_zero() {
        assert_eq!(check_bad(true, kHighsDebugLevelNone), HighsDebugStatus::kOk);
    }

    /// Flag off (1.15.1: the check runs at every debug level), and flag on at
    /// kHighsDebugLevelCheap: the check runs and fails. A failing check logs and then
    /// `assert(!error_found)`s (HSimplexNla.cpp:513), so the kLogicalError return is only
    /// observable with debug assertions off.
    #[test]
    #[cfg(not(debug_assertions))]
    fn h16_nla_debug_gate_runs_check_otherwise() {
        assert_eq!(
            check_bad(false, kHighsDebugLevelNone),
            HighsDebugStatus::kLogicalError
        );
        assert_eq!(
            check_bad(true, kHighsDebugLevelCheap),
            HighsDebugStatus::kLogicalError
        );
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "!error_found")]
    fn h16_nla_debug_gate_off_trips_cpp_assert() {
        check_bad(false, kHighsDebugLevelNone);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "!error_found")]
    fn h16_nla_debug_gate_on_cheap_trips_cpp_assert() {
        check_bad(true, kHighsDebugLevelCheap);
    }

    /// The simplex result on a consistent LP is bit-identical with the flag off and on.
    #[test]
    fn h16_nla_debug_gate_solve_unchanged() {
        let expected = 6.0 - 8.0 / 7.0 + 0.5;
        let mut objectives = Vec::new();
        for gate in [false, true] {
            let mut highs = Highs::default();
            assert_eq!(
                highs.set_option_value_bool("log_to_console", false),
                HighsStatus::kOk
            );
            assert_eq!(
                highs.set_option_value_str("presolve", "off"),
                HighsStatus::kOk
            );
            assert_eq!(
                highs.set_option_value_str("solver", "simplex"),
                HighsStatus::kOk
            );
            assert_eq!(
                highs.set_option_value_bool("h16_nla_debug_gate", gate),
                HighsStatus::kOk
            );
            assert_eq!(highs.pass_model(tiny_lp()), HighsStatus::kOk);
            assert_eq!(highs.run(), HighsStatus::kOk, "gate {gate}");
            assert_eq!(
                highs.get_model_status(),
                HighsModelStatus::kOptimal,
                "gate {gate}"
            );
            let objective = highs.get_info().objective_function_value;
            assert!(
                (objective - expected).abs() < 1e-6,
                "gate {gate}: {objective}"
            );
            objectives.push(objective.to_bits());
        }
        assert_eq!(objectives[0], objectives[1]);
    }
}
