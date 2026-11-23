//! HiGHS source: lp_data/HStruct.h, lp_data/HighsSolution.cpp (HighsSolution/HighsBasis members)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Structs from HStruct.h used by the simplex path. `RefactorInfo::clear` lives in
//! `util/hfactor/refactor.rs` (HFactorRefactor.cpp); `HighsSimplexStats::report/initialise`
//! live in `simplex/ekk/proof.rs` (HEkk.cpp).

use crate::prelude::*;

/// HiGHS: HStruct.h::HighsSolution
#[derive(Clone, Debug)]
pub struct HighsSolution {
    pub value_valid: bool,
    pub dual_valid: bool,
    pub col_value: Vec<f64>,
    pub col_dual: Vec<f64>,
    pub row_value: Vec<f64>,
    pub row_dual: Vec<f64>,
}

impl Default for HighsSolution {
    fn default() -> Self {
        Self {
            value_valid: false,
            dual_valid: false,
            col_value: Default::default(),
            col_dual: Default::default(),
            row_value: Default::default(),
            row_dual: Default::default(),
        }
    }
}

/// HiGHS: HStruct.h::RefactorInfo
#[derive(Clone, Debug)]
pub struct RefactorInfo {
    pub r#use: bool,
    pub pivot_row: Vec<HighsInt>,
    pub pivot_var: Vec<HighsInt>,
    pub pivot_type: Vec<i8>,
    pub build_synthetic_tick: f64,
}

impl Default for RefactorInfo {
    fn default() -> Self {
        Self {
            r#use: false,
            pivot_row: Default::default(),
            pivot_var: Default::default(),
            pivot_type: Default::default(),
            build_synthetic_tick: Default::default(),
        }
    }
}

/// HiGHS: HStruct.h::HotStart
#[derive(Clone, Debug)]
pub struct HotStart {
    pub valid: bool,
    pub refactor_info: RefactorInfo,
    pub nonbasic_move: Vec<i8>,
}

impl Default for HotStart {
    fn default() -> Self {
        Self {
            valid: false,
            refactor_info: Default::default(),
            nonbasic_move: Default::default(),
        }
    }
}

/// HiGHS: HStruct.h::HighsBasis
#[derive(Clone, Debug)]
pub struct HighsBasis {
    pub valid: bool,
    pub alien: bool,
    pub useful: bool,
    pub was_alien: bool,
    pub debug_id: HighsInt,
    pub debug_update_count: HighsInt,
    pub debug_origin_name: String,
    pub col_status: Vec<HighsBasisStatus>,
    pub row_status: Vec<HighsBasisStatus>,
}

impl Default for HighsBasis {
    fn default() -> Self {
        Self {
            valid: false,
            alien: true,
            useful: false,
            was_alien: true,
            debug_id: -1,
            debug_update_count: -1,
            debug_origin_name: "None".to_string(),
            col_status: Default::default(),
            row_status: Default::default(),
        }
    }
}

/// HiGHS: HStruct.h::HighsScale
/// C++ leaves the scalars uninitialised; `HighsLp::clear` -> `clearScale` sets them.
/// Default here = the `HighsLp::clearScale` values (strategy kSimplexScaleStrategyOff, cost 0, empty vectors).
#[derive(Clone, Debug)]
pub struct HighsScale {
    pub strategy: HighsInt,
    pub has_scaling: bool,
    pub num_col: HighsInt,
    pub num_row: HighsInt,
    pub cost: f64,
    pub col: Vec<f64>,
    pub row: Vec<f64>,
}

impl Default for HighsScale {
    fn default() -> Self {
        Self {
            strategy: kSimplexScaleStrategyOff,
            has_scaling: Default::default(),
            num_col: Default::default(),
            num_row: Default::default(),
            cost: 0.0,
            col: Default::default(),
            row: Default::default(),
        }
    }
}

/// HiGHS: HStruct.h::HighsSimplexStats
#[derive(Clone, Debug)]
pub struct HighsSimplexStats {
    pub valid: bool,
    pub iteration_count: HighsInt,
    pub num_invert: HighsInt,
    pub last_invert_num_el: HighsInt,
    pub last_factored_basis_num_el: HighsInt,
    pub col_aq_density: f64,
    pub row_ep_density: f64,
    pub row_ap_density: f64,
    pub row_dse_density: f64,
}

impl Default for HighsSimplexStats {
    fn default() -> Self {
        Self {
            valid: Default::default(),
            iteration_count: Default::default(),
            num_invert: Default::default(),
            last_invert_num_el: Default::default(),
            last_factored_basis_num_el: Default::default(),
            col_aq_density: Default::default(),
            row_ep_density: Default::default(),
            row_ap_density: Default::default(),
            row_dse_density: Default::default(),
        }
    }
}

impl HighsSolution {
    /// HiGHS: HighsSolution.cpp::HighsSolution::invalidate (lines 2154-2157, 4 lines)
    pub fn invalidate(&mut self) {
        self.value_valid = false;
        self.dual_valid = false;
    }

    /// HiGHS: HighsSolution.cpp::HighsSolution::clear (lines 2159-2165, 7 lines)
    pub fn clear(&mut self) {
        self.invalidate();
        self.col_value.clear();
        self.row_value.clear();
        self.col_dual.clear();
        self.row_dual.clear();
    }
}

impl HighsBasis {
    /// HiGHS: HighsSolution.cpp::HighsBasis::invalidate (lines 2223-2231, 9 lines)
    pub fn invalidate(&mut self) {
        self.valid = false;
        self.alien = true;
        self.useful = false;
        self.was_alien = true;
        self.debug_id = -1;
        self.debug_update_count = -1;
        self.debug_origin_name = "None".to_string();
    }

    /// HiGHS: HighsSolution.cpp::HighsBasis::clear (lines 2233-2237, 5 lines)
    pub fn clear(&mut self) {
        self.invalidate();
        self.row_status.clear();
        self.col_status.clear();
    }
}

/// HiGHS: HStruct.h::HighsPresolveRuleLog (C++ fields uninitialised;
/// always set by `HighsPresolveLog::clear` before use).
#[derive(Clone, Copy, Debug, Default)]
pub struct HighsPresolveRuleLog {
    pub call: HighsInt,
    pub col_removed: HighsInt,
    pub row_removed: HighsInt,
}

/// HiGHS: HStruct.h::HighsPresolveLog
#[derive(Clone, Debug, Default)]
pub struct HighsPresolveLog {
    pub rule: Vec<HighsPresolveRuleLog>,
}

impl HighsPresolveLog {
    /// HiGHS: HPresolveAnalysis.cpp::HighsPresolveLog::clear (lines 76-83, 8 lines)
    pub fn clear(&mut self) {
        self.rule
            .resize(kPresolveRuleCount as usize, HighsPresolveRuleLog::default());
        for rule_type in 0..kPresolveRuleCount as usize {
            self.rule[rule_type].call = 0;
            self.rule[rule_type].col_removed = 0;
            self.rule[rule_type].row_removed = 0;
        }
    }
}

/// HiGHS: HStruct.h::HighsUserScaleData (lines 224-250)
///
/// The C++ leaves the non-defaulted members uninitialised until `initialise`, which only
/// `Highs::userScale` calls for non-zero user scales; the port zero-initialises them.
/// Its methods are in highs_user_scale.rs.
#[derive(Clone, Debug, Default)]
pub struct HighsUserScaleData {
    pub user_objective_scale: HighsInt,
    pub user_bound_scale: HighsInt,
    pub infinite_cost: f64,
    pub infinite_bound: f64,
    pub small_matrix_value: f64,
    pub large_matrix_value: f64,
    pub num_infinite_costs: HighsInt,
    pub num_infinite_hessian_values: HighsInt,
    pub num_infinite_col_bounds: HighsInt,
    pub num_infinite_row_bounds: HighsInt,
    pub num_small_matrix_values: HighsInt,
    pub num_large_matrix_values: HighsInt,
    pub suggested_user_objective_scale: HighsInt,
    pub suggested_user_bound_scale: HighsInt,
    pub applied: bool,
}

/// HiGHS: HStruct.h::HighsLpMods (lines 98-127): the semi-variable and
/// infinite-cost modifications `Highs::optimizeModel` makes to `model_.lp_`
/// and undoes in `returnFromOptimizeModel` (`HighsLp::unapplyMods`,
/// `Highs::restoreInfCost`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HighsLpMods {
    // Semi-variables with zero lower bound that are treated as non-semi
    pub save_non_semi_variable_index: Vec<HighsInt>,
    // Semi-variables with inconsistent bounds that are fixed at zero
    pub save_inconsistent_semi_variable_index: Vec<HighsInt>,
    pub save_inconsistent_semi_variable_lower_bound_value: Vec<f64>,
    pub save_inconsistent_semi_variable_upper_bound_value: Vec<f64>,
    pub save_inconsistent_semi_variable_type: Vec<HighsVarType>,
    // Semi-variables whose lower bound is ignored when solving the relaxation
    pub save_relaxed_semi_variable_lower_bound_index: Vec<HighsInt>,
    pub save_relaxed_semi_variable_lower_bound_value: Vec<f64>,
    // Semi-variables whose upper bound is too large to be used as a big-M
    pub save_tightened_semi_variable_upper_bound_index: Vec<HighsInt>,
    pub save_tightened_semi_variable_upper_bound_value: Vec<f64>,
    // Variables with infinite costs that are fixed during solve
    pub save_inf_cost_variable_index: Vec<HighsInt>,
    pub save_inf_cost_variable_cost: Vec<f64>,
    pub save_inf_cost_variable_lower: Vec<f64>,
    pub save_inf_cost_variable_upper: Vec<f64>,
}

impl HighsLpMods {
    /// HiGHS: HighsLp.cpp::HighsLpMods::clear (lines 495-509)
    pub fn clear(&mut self) {
        self.save_non_semi_variable_index.clear();
        self.save_inconsistent_semi_variable_index.clear();
        self.save_inconsistent_semi_variable_lower_bound_value
            .clear();
        self.save_inconsistent_semi_variable_upper_bound_value
            .clear();
        self.save_inconsistent_semi_variable_type.clear();
        self.save_relaxed_semi_variable_lower_bound_index.clear();
        self.save_relaxed_semi_variable_lower_bound_value.clear();
        self.save_tightened_semi_variable_upper_bound_index.clear();
        self.save_tightened_semi_variable_upper_bound_value.clear();
        self.save_inf_cost_variable_index.clear();
        self.save_inf_cost_variable_cost.clear();
        self.save_inf_cost_variable_lower.clear();
        self.save_inf_cost_variable_upper.clear();
    }

    /// HiGHS: HighsLp.cpp::HighsLpMods::isClear (lines 511-524). HiGHS quirk: the
    /// four `save_inf_cost_*` vectors are not tested.
    pub fn is_clear(&self) -> bool {
        if !self.save_non_semi_variable_index.is_empty() {
            return false;
        }
        if !self.save_inconsistent_semi_variable_index.is_empty() {
            return false;
        }
        if !self
            .save_inconsistent_semi_variable_lower_bound_value
            .is_empty()
        {
            return false;
        }
        if !self
            .save_inconsistent_semi_variable_upper_bound_value
            .is_empty()
        {
            return false;
        }
        if !self.save_inconsistent_semi_variable_type.is_empty() {
            return false;
        }
        if !self.save_relaxed_semi_variable_lower_bound_index.is_empty() {
            return false;
        }
        if !self.save_relaxed_semi_variable_lower_bound_value.is_empty() {
            return false;
        }
        if !self
            .save_tightened_semi_variable_upper_bound_index
            .is_empty()
        {
            return false;
        }
        if !self
            .save_tightened_semi_variable_upper_bound_value
            .is_empty()
        {
            return false;
        }
        true
    }
}

/// HiGHS: HStruct.h::HighsIllConditioningRecord (lines 149-152). The C++ members are
/// uninitialised; the port zero-initialises them.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HighsIllConditioningRecord {
    pub index: HighsInt,
    pub multiplier: f64,
}

/// HiGHS: HStruct.h::HighsIllConditioning (lines 154-157)
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HighsIllConditioning {
    pub record: Vec<HighsIllConditioningRecord>,
}

impl HighsIllConditioning {
    /// HiGHS: HighsInterface.cpp::HighsIllConditioning::clear (line 3204)
    pub fn clear(&mut self) {
        self.record.clear();
    }
}

/// HiGHS: HStruct.h::HighsLinearObjective (lines 159-167)
#[derive(Clone, Debug, PartialEq)]
pub struct HighsLinearObjective {
    pub weight: f64,
    pub offset: f64,
    pub coefficients: Vec<f64>,
    pub abs_tolerance: f64,
    pub rel_tolerance: f64,
    pub priority: HighsInt,
}

impl Default for HighsLinearObjective {
    /// HiGHS: HStruct.h default member initialisers (lines 160-165)
    fn default() -> Self {
        Self {
            weight: 0.0,
            offset: 0.0,
            coefficients: Vec::new(),
            abs_tolerance: -1.0,
            rel_tolerance: -1.0,
            priority: 0,
        }
    }
}

impl HighsLinearObjective {
    /// HiGHS: HighsInterface.cpp::HighsLinearObjective::clear (lines 4278-4285). HiGHS
    /// quirk: the tolerances are cleared to 0, not to their -1 initialisers.
    pub fn clear(&mut self) {
        self.weight = 0.0;
        self.offset = 0.0;
        self.coefficients.clear();
        self.abs_tolerance = 0.0;
        self.rel_tolerance = 0.0;
        self.priority = 0;
    }
}

#[cfg(test)]
mod g1_struct_tests {
    use super::*;

    #[test]
    fn lp_mods_clear_and_is_clear() {
        let mut mods = HighsLpMods::default();
        assert!(mods.is_clear());
        mods.save_inf_cost_variable_index.push(1);
        // HiGHS quirk: isClear ignores the infinite-cost vectors
        assert!(mods.is_clear());
        mods.save_non_semi_variable_index.push(0);
        assert!(!mods.is_clear());
        mods.clear();
        assert!(mods.is_clear() && mods.save_inf_cost_variable_index.is_empty());
    }

    #[test]
    fn linear_objective_default_and_clear() {
        let mut o = HighsLinearObjective::default();
        assert_eq!((o.abs_tolerance, o.rel_tolerance), (-1.0, -1.0));
        o.coefficients.push(1.0);
        o.clear();
        assert_eq!((o.abs_tolerance, o.rel_tolerance), (0.0, 0.0));
        assert!(o.coefficients.is_empty());
    }
}
