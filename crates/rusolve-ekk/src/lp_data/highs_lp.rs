//! HiGHS source: lp_data/HighsLp.h, lp_data/HighsLp.cpp (subset)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Subset: the members the simplex path touches. Dropped: col_hash_/row_hash_ (name hashing).
//! `mods_` (semi-variable / infinite-cost modifications) is ported too;
//! `unapply_mods` is in lp_utils_semi.rs.
//! `Default` reproduces `HighsLp() { clear(); }`.

use crate::lp_data::hstruct::{HighsLpMods, HighsScale};
use crate::prelude::*;
use crate::util::highs_cdouble::HighsCDouble;
use crate::util::sparse_matrix::HighsSparseMatrix;

/// HiGHS: HighsLp.h::HighsLp
#[derive(Clone, Debug)]
pub struct HighsLp {
    pub num_col_: HighsInt,
    pub num_row_: HighsInt,
    pub col_cost_: Vec<f64>,
    pub col_lower_: Vec<f64>,
    pub col_upper_: Vec<f64>,
    pub row_lower_: Vec<f64>,
    pub row_upper_: Vec<f64>,
    pub a_matrix_: HighsSparseMatrix,
    pub sense_: ObjSense,
    pub offset_: f64,
    pub model_name_: String,
    pub origin_name_: String,
    pub objective_name_: String,
    pub col_name_prefix_: String,
    pub row_name_prefix_: String,
    pub col_name_suffix_: HighsInt,
    pub row_name_suffix_: HighsInt,
    pub col_names_: Vec<String>,
    pub row_names_: Vec<String>,
    pub integrality_: Vec<HighsVarType>,
    pub scale_: HighsScale,
    pub is_scaled_: bool,
    pub is_moved_: bool,
    pub cost_row_location_: HighsInt,
    pub has_infinite_cost_: bool,
    /// HiGHS 1.16 (#3255): HighsLp.h `HighsInt fme_obj_col_ = -1` (1.16.0 line 61), the
    /// objective column z added by the Fourier-Motzkin objective reformulation. Only FM
    /// presolve (`mip_fm_presolve`) sets it; -1 everywhere else.
    pub fme_obj_col_: HighsInt,
    /// HiGHS: HighsLp.h `HighsLpMods mods_` (line 58)
    pub mods_: HighsLpMods,
}

impl Default for HighsLp {
    /// HiGHS: HighsLp.h::HighsLp::HighsLp() { clear(); } -- values written out from HighsLp::clear.
    fn default() -> Self {
        Self {
            num_col_: 0,
            num_row_: 0,
            col_cost_: Vec::new(),
            col_lower_: Vec::new(),
            col_upper_: Vec::new(),
            row_lower_: Vec::new(),
            row_upper_: Vec::new(),
            a_matrix_: HighsSparseMatrix::default(),
            sense_: ObjSense::kMinimize,
            offset_: 0.0,
            model_name_: String::new(),
            origin_name_: String::new(),
            objective_name_: String::new(),
            col_name_prefix_: String::new(),
            row_name_prefix_: String::new(),
            col_name_suffix_: 0,
            row_name_suffix_: 0,
            col_names_: Vec::new(),
            row_names_: Vec::new(),
            integrality_: Vec::new(),
            scale_: HighsScale::default(),
            is_scaled_: false,
            is_moved_: false,
            cost_row_location_: -1,
            has_infinite_cost_: false,
            fme_obj_col_: -1,
            mods_: HighsLpMods::default(),
        }
    }
}

impl HighsLp {
    /// HiGHS: HighsLp.cpp::HighsLp::isMip (lines 18-26, 9 lines)
    pub fn is_mip(&self) -> bool {
        let integrality_size = self.integrality_.len();
        if integrality_size != 0 {
            debug_assert!(integrality_size as HighsInt == self.num_col_);
            for i_col in 0..self.num_col_ {
                if self.integrality_[i_col as usize] != HighsVarType::kContinuous {
                    return true;
                }
            }
        }
        false
    }

    /// HiGHS: HighsLp.cpp::HighsLp::hasInfiniteCost (lines 28-34, 7 lines)
    pub fn has_infinite_cost(&self, infinite_cost: f64) -> bool {
        for i_col in 0..self.num_col_ {
            if self.col_cost_[i_col as usize] >= infinite_cost {
                return true;
            }
            if self.col_cost_[i_col as usize] <= -infinite_cost {
                return true;
            }
        }
        false
    }

    /// HiGHS: HighsLp.cpp::HighsLp::hasSemiVariables (lines 36-46, 11 lines)
    pub fn has_semi_variables(&self) -> bool {
        let integrality_size = self.integrality_.len();
        if integrality_size != 0 {
            debug_assert!(integrality_size as HighsInt == self.num_col_);
            for i_col in 0..self.num_col_ {
                if self.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
                    || self.integrality_[i_col as usize] == HighsVarType::kSemiInteger
                {
                    return true;
                }
            }
        }
        false
    }

    /// HiGHS: HighsLp.cpp::HighsLp::hasMods (lines 48-54, 7 lines)
    pub fn has_mods(&self) -> bool {
        !self.mods_.save_non_semi_variable_index.is_empty()
            || !self.mods_.save_inconsistent_semi_variable_index.is_empty()
            || !self
                .mods_
                .save_relaxed_semi_variable_lower_bound_index
                .is_empty()
            || !self
                .mods_
                .save_tightened_semi_variable_upper_bound_index
                .is_empty()
            || !self.mods_.save_inf_cost_variable_index.is_empty()
    }

    /// HiGHS: HighsLp.cpp::HighsLp::needsMods (lines 56-59, 4 lines)
    pub fn needs_mods(&self, infinite_cost: f64) -> bool {
        debug_assert!(self.has_infinite_cost_ == self.has_infinite_cost(infinite_cost));
        self.has_infinite_cost_ || self.has_semi_variables()
    }

    /// HiGHS: HighsLp.cpp::HighsLp::objectiveValue (lines 124-130, 7 lines)
    /// HiGHS 1.16 (#3513): HighsLp.cpp::HighsLp::objectiveValue (lines 124-127)
    pub fn objective_value(&self, solution: &[f64], h16_obj_cdouble: bool) -> f64 {
        debug_assert!(solution.len() as HighsInt >= self.num_col_);
        // HiGHS 1.16 #3513 (h16_obj_cdouble): HighsLp.cpp line 126 (1.16.0)
        if h16_obj_cdouble {
            return f64::from(self.objective_c_double_value(solution));
        }
        let mut objective_function_value = self.offset_;
        for i_col in 0..self.num_col_ {
            objective_function_value += self.col_cost_[i_col as usize] * solution[i_col as usize];
        }
        objective_function_value
    }

    /// HiGHS: HighsLp.cpp::HighsLp::objectiveCDoubleValue (lines 132-140, 9 lines)
    pub fn objective_c_double_value(&self, solution: &[f64]) -> HighsCDouble {
        debug_assert!(solution.len() as HighsInt >= self.num_col_);
        let mut objective_function_value = HighsCDouble::from(self.offset_);
        for i_col in 0..self.num_col_ {
            objective_function_value +=
                HighsCDouble::from(self.col_cost_[i_col as usize]) * solution[i_col as usize];
        }
        objective_function_value
    }

    /// HiGHS: HighsLp.cpp::HighsLp::setMatrixDimensions (lines 142-145, 4 lines)
    pub fn set_matrix_dimensions(&mut self) {
        self.a_matrix_.num_col_ = self.num_col_;
        self.a_matrix_.num_row_ = self.num_row_;
    }

    /// HiGHS: HighsLp.cpp::HighsLp::resetScale (lines 147-153, 7 lines)
    pub fn reset_scale(&mut self) {
        // The C++ dimension check is commented out upstream; scaling is always cleared.
        self.clear_scale();
    }

    /// HiGHS: HighsLp.cpp::HighsLp::setFormat (lines 155-157, 3 lines)
    pub fn set_format(&mut self, format: MatrixFormat) {
        self.a_matrix_.set_format(format);
    }

    /// HiGHS: HighsLp.cpp::HighsLp::exactResize (lines 159-170, 12 lines)
    pub fn exact_resize(&mut self) {
        self.col_cost_.resize(self.num_col_ as usize, 0.0);
        self.col_lower_.resize(self.num_col_ as usize, 0.0);
        self.col_upper_.resize(self.num_col_ as usize, 0.0);
        self.row_lower_.resize(self.num_row_ as usize, 0.0);
        self.row_upper_.resize(self.num_row_ as usize, 0.0);
        self.a_matrix_.exact_resize();

        if !self.col_names_.is_empty() {
            self.col_names_
                .resize(self.num_col_ as usize, String::new());
        }
        if !self.row_names_.is_empty() {
            self.row_names_
                .resize(self.num_row_ as usize, String::new());
        }
        if !self.integrality_.is_empty() {
            self.integrality_
                .resize(self.num_col_ as usize, HighsVarType::kContinuous);
        }
    }

    /// HiGHS: HighsLp.cpp::HighsLp::clear (lines 193-230, 38 lines)
    /// Also clears what C++ clears in col_hash_/row_hash_ (dropped here).
    pub fn clear(&mut self) {
        self.num_col_ = 0;
        self.num_row_ = 0;

        self.col_cost_.clear();
        self.col_lower_.clear();
        self.col_upper_.clear();
        self.row_lower_.clear();
        self.row_upper_.clear();

        self.a_matrix_.clear();

        self.sense_ = ObjSense::kMinimize;
        self.offset_ = 0.0;

        self.model_name_.clear();
        self.origin_name_.clear();
        self.objective_name_.clear();

        self.col_name_prefix_.clear();
        self.row_name_prefix_.clear();
        self.col_name_suffix_ = 0;
        self.row_name_suffix_ = 0;
        self.col_names_.clear();
        self.row_names_.clear();

        self.integrality_.clear();

        // C++ col_hash_.clear() / row_hash_.clear(): dropped (see module doc).

        self.clear_scale();
        self.is_scaled_ = false;
        self.is_moved_ = false;
        self.cost_row_location_ = -1;
        self.has_infinite_cost_ = false;
        // HiGHS 1.16 (#3255): HighsLp.cpp::clear line 226 (1.16.0)
        self.fme_obj_col_ = -1;
        self.mods_.clear();
    }

    /// HiGHS: HighsLp.cpp::HighsLp::clearScale (lines 232-240, 9 lines)
    pub fn clear_scale(&mut self) {
        self.scale_.strategy = kSimplexScaleStrategyOff;
        self.scale_.has_scaling = false;
        self.scale_.num_col = 0;
        self.scale_.num_row = 0;
        self.scale_.cost = 0.0;
        self.scale_.col.clear();
        self.scale_.row.clear();
    }

    /// HiGHS: HighsLp.cpp::HighsLp::clearScaling (lines 242-245, 4 lines)
    pub fn clear_scaling(&mut self) {
        self.unapply_scale();
        self.clear_scale();
    }

    /// HiGHS: HighsLp.cpp::HighsLp::applyScale (lines 247-272, 26 lines)
    pub fn apply_scale(&mut self) {
        // Ensure that any scaling is applied
        if self.is_scaled_ {
            // Already scaled - so check that there is scaling and return
            debug_assert!(self.scale_.has_scaling);
            return;
        }
        // No scaling currently applied
        self.is_scaled_ = false;
        if self.scale_.has_scaling {
            // Apply the scaling to the bounds, costs and matrix, and record
            // that it has been applied
            for i_col in 0..self.num_col_ {
                let ic = i_col as usize;
                self.col_lower_[ic] /= self.scale_.col[ic];
                self.col_upper_[ic] /= self.scale_.col[ic];
                self.col_cost_[ic] *= self.scale_.col[ic];
            }
            for i_row in 0..self.num_row_ {
                let ir = i_row as usize;
                self.row_lower_[ir] *= self.scale_.row[ir];
                self.row_upper_[ir] *= self.scale_.row[ir];
            }
            self.a_matrix_.apply_scale(&self.scale_);
            self.is_scaled_ = true;
        }
    }

    /// HiGHS: HighsLp.cpp::HighsLp::unapplyScale (lines 274-296, 23 lines)
    pub fn unapply_scale(&mut self) {
        // Ensure that any scaling is not applied
        if !self.is_scaled_ {
            // Not scaled so return
            return;
        }
        // Already scaled - so check that there is scaling
        debug_assert!(self.scale_.has_scaling);
        // Unapply the scaling to the bounds, costs and matrix, and record
        // that it has been unapplied
        for i_col in 0..self.num_col_ {
            let ic = i_col as usize;
            self.col_lower_[ic] *= self.scale_.col[ic];
            self.col_upper_[ic] *= self.scale_.col[ic];
            self.col_cost_[ic] /= self.scale_.col[ic];
        }
        for i_row in 0..self.num_row_ {
            let ir = i_row as usize;
            self.row_lower_[ir] /= self.scale_.row[ir];
            self.row_upper_[ir] /= self.scale_.row[ir];
        }
        self.a_matrix_.unapply_scale(&self.scale_);
        self.is_scaled_ = false;
    }

    /// HiGHS: HighsLp.cpp::HighsLp::moveBackLpAndUnapplyScaling (lines 298-303, 6 lines)
    pub fn move_back_lp_and_unapply_scaling(&mut self, lp: &mut HighsLp) {
        debug_assert!(self.is_moved_);
        // C++ `*this = std::move(lp)`; std::mem::take leaves Default in `lp`.
        *self = std::mem::take(lp);
        self.unapply_scale();
        debug_assert!(!self.is_moved_);
    }

    /// HiGHS: HighsLp.h::HighsLp::ensureColwise (inline)
    pub fn ensure_colwise(&mut self) {
        self.a_matrix_.ensure_colwise();
    }

    /// HiGHS: HighsLp.h::HighsLp::ensureRowwise (inline)
    pub fn ensure_rowwise(&mut self) {
        self.a_matrix_.ensure_rowwise();
    }
}

#[cfg(test)]
mod h16_obj_cdouble_tests {
    use super::*;

    /// offset 1e16, costs [1, 1], x = [1, 1]: in double each `+ 1` rounds
    /// back to 1e16 (tie to even), while HighsCDouble keeps the exact 1e16 + 2.
    fn cancel_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf];
        lp.offset_ = 1e16;
        lp
    }

    #[test]
    fn h16_obj_cdouble_objective_value_on() {
        // HighsLp.cpp line 126 (1.16.0): double(objectiveCDoubleValue(solution))
        let lp = cancel_lp();
        assert_eq!(lp.objective_value(&[1.0, 1.0], true), 10000000000000002.0);
    }

    #[test]
    fn h16_obj_cdouble_objective_value_off() {
        // HighsLp.cpp lines 124-130 (1.15.1): plain double accumulation
        let lp = cancel_lp();
        assert_eq!(lp.objective_value(&[1.0, 1.0], false), 1e16);
    }
}
