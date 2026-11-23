//! HiGHS source: lp_data/HighsInterface.cpp (userScale, userUnscale,
//! userScaleModel, userScaleSolution), lp_data/HighsLpUtils.cpp (userScaleLp,
//! userScaleCosts, userScaleColBounds, userScaleRowBounds, userScaleMatrix,
//! userScaleStatus, initialiseUserScaleData, HighsUserScaleData methods),
//! lp_data/Highs.cpp (getObjectiveBoundScaling)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The user objective and bound scaling that `Highs::run` applies
//! before `optimizeHighs` (`user_scale`) and removes after it
//! (`user_unscale`), driven by the options `user_objective_scale` and
//! `user_bound_scale` (both 0 by default, so `user_scale` returns at once and
//! `user_unscale` sees `data.applied == false`). The Hessian part,
//! `userScaleHessian`, is already ported (`model/highs_hessian_utils/normalise.rs`).
//! `logHeader` (in `getObjectiveBoundScaling`) is logging and is dropped.

use crate::lp_data::highs_solve::assess_excessive_objective_bound_scaling;
use crate::lp_data::hstruct::HighsUserScaleData;
use crate::lp_data::presolve_support::report_kkt_failures;
use crate::prelude::*;
use crate::util::sparse_matrix::HighsSparseMatrix;

impl HighsUserScaleData {
    /// HiGHS: HighsLpUtils.cpp::HighsUserScaleData::initialise (lines 3697-3718)
    pub fn initialise(
        &mut self,
        user_objective_scale_: HighsInt,
        user_bound_scale_: HighsInt,
        infinite_cost_: f64,
        infinite_bound_: f64,
        small_matrix_value_: f64,
        large_matrix_value_: f64,
    ) {
        self.user_objective_scale = user_objective_scale_;
        self.user_bound_scale = user_bound_scale_;
        self.infinite_cost = infinite_cost_;
        self.infinite_bound = infinite_bound_;
        self.small_matrix_value = small_matrix_value_;
        self.large_matrix_value = large_matrix_value_;
        self.num_infinite_costs = 0;
        self.num_infinite_hessian_values = 0;
        self.num_infinite_col_bounds = 0;
        self.num_infinite_row_bounds = 0;
        self.num_small_matrix_values = 0;
        self.num_large_matrix_values = 0;
        self.suggested_user_objective_scale = 0;
        self.suggested_user_bound_scale = 0;
        self.applied = false;
    }

    /// HiGHS: HighsLpUtils.cpp::HighsUserScaleData::scaleError (lines 3720-3783).
    /// The `std::stringstream` text is built with `format!`/`push_str`;
    /// `ss << int` prints the plain decimal value.
    pub fn scale_error(&self, message: &mut String) -> bool {
        if self.num_infinite_costs
            + self.num_infinite_hessian_values
            + self.num_infinite_col_bounds
            + self.num_infinite_row_bounds
            + self.num_large_matrix_values
            == 0
        {
            return false;
        }
        debug_assert!(self.user_objective_scale != 0 || self.user_bound_scale != 0);
        // std::stringstream ss; ss.str(std::string());
        let mut ss = String::new();
        ss.push_str("User scaling of");
        if self.user_objective_scale != 0 {
            ss.push_str(&format!(" 2**({}) for costs", self.user_objective_scale));
        }
        if self.user_bound_scale != 0 {
            if self.user_objective_scale != 0 {
                ss.push_str(" and");
            }
            ss.push_str(&format!(" 2**({}) for bounds", self.user_bound_scale));
        }
        ss.push_str(" yields");
        if self.num_infinite_costs != 0 {
            ss.push_str(&format!(" {} infinite cost", self.num_infinite_costs));
            if self.num_infinite_costs > 1 {
                ss.push('s');
            }
        }
        if self.num_infinite_hessian_values != 0 {
            if self.num_infinite_costs != 0 {
                if self.num_infinite_col_bounds != 0 || self.num_infinite_row_bounds != 0 {
                    ss.push(',');
                } else {
                    ss.push_str(" and");
                }
            }
            ss.push_str(&format!(
                " {} infinite Hessian values",
                self.num_infinite_hessian_values
            ));
            if self.num_infinite_hessian_values > 1 {
                ss.push('s');
            }
        }
        if self.num_infinite_col_bounds != 0 {
            if self.num_infinite_costs != 0 || self.num_infinite_hessian_values != 0 {
                if self.num_infinite_row_bounds != 0 {
                    ss.push(',');
                } else {
                    ss.push_str(" and");
                }
            }
            ss.push_str(&format!(
                " {} infinite column bound",
                self.num_infinite_col_bounds
            ));
            if self.num_infinite_col_bounds > 1 {
                ss.push('s');
            }
        }
        if self.num_infinite_row_bounds != 0 {
            if self.num_infinite_costs != 0
                || self.num_infinite_hessian_values != 0
                || self.num_infinite_col_bounds != 0
            {
                ss.push_str(" and");
            }
            ss.push_str(&format!(
                " {} infinite row bound",
                self.num_infinite_row_bounds
            ));
            if self.num_infinite_row_bounds > 1 {
                ss.push('s');
            }
        }
        if self.num_large_matrix_values != 0 {
            if self.num_infinite_costs
                + self.num_infinite_hessian_values
                + self.num_infinite_col_bounds
                + self.num_infinite_row_bounds
                > 0
            {
                ss.push_str(", and");
            }
            ss.push_str(&format!(
                " {} large matrix value",
                self.num_large_matrix_values
            ));
            if self.num_large_matrix_values > 1 {
                ss.push('s');
            }
        }
        ss.push('\n');
        *message = ss;
        true
    }

    /// HiGHS: HighsLpUtils.cpp::HighsUserScaleData::scaleWarning (lines 3785-3797)
    pub fn scale_warning(&self, message: &mut String) -> bool {
        if self.num_small_matrix_values == 0 {
            return false;
        }
        debug_assert!(self.user_bound_scale != 0);
        // std::stringstream ss; ss.str(std::string());
        let mut ss = String::new();
        ss.push_str(&format!(
            "User scaling of 2**({}) for bounds yields {} small matrix value",
            self.user_bound_scale, self.num_small_matrix_values
        ));
        if self.num_small_matrix_values > 1 {
            ss.push('s');
        }
        ss.push('\n');
        *message = ss;
        true
    }
}

/// HiGHS: HighsLpUtils.cpp::initialiseUserScaleData (lines 3689-3695)
pub fn initialise_user_scale_data(
    options: &HighsOptions,
    user_scale_data: &mut HighsUserScaleData,
) {
    user_scale_data.initialise(
        options.user_objective_scale,
        options.user_bound_scale,
        options.infinite_cost,
        options.infinite_bound,
        options.small_matrix_value,
        options.large_matrix_value,
    );
}

/// HiGHS: HighsLpUtils.cpp::userScaleLp (lines 772-779), the overload
/// `userScaleLp(HighsLp&, HighsUserScaleData&, const HighsLogOptions&)`:
/// assess, report through `userScaleStatus`, then apply unless that is `kError`.
pub fn user_scale_lp_log(
    lp: &mut HighsLp,
    data: &mut HighsUserScaleData,
    log_options: &HighsLogOptions,
) -> HighsStatus {
    // This overload has no options in reach, so the `h16_user_scale_offset` flag is false
    // here. The 1.16.0 C++ overload (HighsLpUtils.cpp lines 626-633) would
    // scale lp.offset_ on its apply step. It has no Rust caller.
    user_scale_lp(lp, data, false, false);
    let return_status = user_scale_status(log_options, data);
    if return_status == HighsStatus::kError {
        return HighsStatus::kError;
    }
    user_scale_lp(lp, data, true, false);
    return_status
}

/// HiGHS: HighsLpUtils.cpp::userScaleLp (lines 781-787), the overload
/// `userScaleLp(HighsLp&, HighsUserScaleData&, const bool apply = true)`.
/// HiGHS 1.16 (#3186): lp_data/HighsLpUtils.cpp::userScaleLp (lines 635-641)
/// passes `lp.offset_` to `user_scale_costs`; `h16_user_scale_offset` is the
/// option that enables that.
pub fn user_scale_lp(
    lp: &mut HighsLp,
    data: &mut HighsUserScaleData,
    apply: bool,
    h16_user_scale_offset: bool,
) {
    user_scale_costs(
        &lp.integrality_,
        &mut lp.offset_,
        &mut lp.col_cost_,
        data,
        apply,
        h16_user_scale_offset,
    );
    user_scale_col_bounds(
        &lp.integrality_,
        &mut lp.col_lower_,
        &mut lp.col_upper_,
        data,
        apply,
    );
    user_scale_matrix(&lp.integrality_, &mut lp.a_matrix_, data, apply);
    user_scale_row_bounds(&mut lp.row_lower_, &mut lp.row_upper_, data, apply);
}

/// HiGHS: HighsLpUtils.cpp::userScaleCosts (lines 789-811)
/// HiGHS 1.16 (#3186): lp_data/HighsLpUtils.cpp::userScaleCosts (lines 643-666)
/// takes `double& offset` and, at the end, `if (apply) offset *= (bound_scale_value
/// * objective_scale_value)`. That step is gated by `h16_user_scale_offset`.
pub fn user_scale_costs(
    integrality: &[HighsVarType],
    offset: &mut f64,
    cost: &mut [f64],
    data: &mut HighsUserScaleData,
    apply: bool,
    h16_user_scale_offset: bool,
) {
    data.num_infinite_costs = 0;
    let user_bound_scale = data.user_bound_scale;
    let user_objective_scale = data.user_objective_scale;
    if user_bound_scale == 0 && user_objective_scale == 0 {
        return;
    }
    let num_col = cost.len() as HighsInt;
    if num_col <= 0 {
        return;
    }
    let integrality_size = integrality.len() as HighsInt;
    let has_integrality = integrality_size > 0;
    let bound_scale_value = 2f64.powf(user_bound_scale as f64);
    let objective_scale_value = 2f64.powf(user_objective_scale as f64);
    debug_assert!(!has_integrality || integrality_size >= num_col);
    for i_col in 0..num_col {
        let mut value = cost[i_col as usize];
        if has_integrality && integrality[i_col as usize] != HighsVarType::kContinuous {
            value *= bound_scale_value;
        }
        value *= objective_scale_value;
        if value.abs() > data.infinite_cost {
            data.num_infinite_costs += 1;
        }
        if apply {
            cost[i_col as usize] = value;
        }
    }
    // HiGHS 1.16 #3186 (h16_user_scale_offset): HighsLpUtils.cpp line 665 (1.16.0)
    if h16_user_scale_offset && apply {
        *offset *= bound_scale_value * objective_scale_value;
    }
}

/// HiGHS: HighsLpUtils.cpp::userScaleColBounds (lines 813-841)
pub fn user_scale_col_bounds(
    integrality: &[HighsVarType],
    lower: &mut [f64],
    upper: &mut [f64],
    data: &mut HighsUserScaleData,
    apply: bool,
) {
    data.num_infinite_col_bounds = 0;
    let user_bound_scale = data.user_bound_scale;
    if user_bound_scale == 0 {
        return;
    }
    let num_col = lower.len() as HighsInt;
    if num_col <= 0 {
        return;
    }
    let integrality_size = integrality.len() as HighsInt;
    let has_integrality = integrality_size > 0;
    debug_assert!(!has_integrality || integrality_size >= num_col);
    let bound_scale_value = 2f64.powf(user_bound_scale as f64);
    for i_col in 0..num_col {
        if !has_integrality || integrality[i_col as usize] == HighsVarType::kContinuous {
            if lower[i_col as usize] > -kHighsInf {
                let value = lower[i_col as usize] * bound_scale_value;
                if value.abs() > data.infinite_bound {
                    data.num_infinite_col_bounds += 1;
                }
                if apply {
                    lower[i_col as usize] = value;
                }
            }
            if upper[i_col as usize] < kHighsInf {
                let value = upper[i_col as usize] * bound_scale_value;
                if value.abs() > data.infinite_bound {
                    data.num_infinite_col_bounds += 1;
                }
                if apply {
                    upper[i_col as usize] = value;
                }
            }
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::userScaleRowBounds (lines 843-863)
pub fn user_scale_row_bounds(
    lower: &mut [f64],
    upper: &mut [f64],
    data: &mut HighsUserScaleData,
    apply: bool,
) {
    data.num_infinite_row_bounds = 0;
    let user_bound_scale = data.user_bound_scale;
    if user_bound_scale == 0 {
        return;
    }
    let num_row = lower.len() as HighsInt;
    if num_row <= 0 {
        return;
    }
    let bound_scale_value = 2f64.powf(user_bound_scale as f64);
    for i_row in 0..num_row {
        if lower[i_row as usize] > -kHighsInf {
            let value = lower[i_row as usize] * bound_scale_value;
            if value.abs() > data.infinite_bound {
                data.num_infinite_row_bounds += 1;
            }
            if apply {
                lower[i_row as usize] = value;
            }
        }
        if upper[i_row as usize] < kHighsInf {
            let value = upper[i_row as usize] * bound_scale_value;
            if value.abs() > data.infinite_bound {
                data.num_infinite_row_bounds += 1;
            }
            if apply {
                upper[i_row as usize] = value;
            }
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::userScaleMatrix (lines 865-909)
pub fn user_scale_matrix(
    integrality: &[HighsVarType],
    matrix: &mut HighsSparseMatrix,
    data: &mut HighsUserScaleData,
    apply: bool,
) {
    data.num_small_matrix_values = 0;
    data.num_large_matrix_values = 0;
    let user_bound_scale = data.user_bound_scale;
    if user_bound_scale == 0 {
        return;
    }
    if integrality.is_empty() {
        return;
    }
    let num_col = matrix.num_col_;
    if num_col <= 0 {
        return;
    }
    let num_row = matrix.num_row_;
    if num_row <= 0 {
        return;
    }
    debug_assert!(integrality.len() as HighsInt >= num_col);
    let bound_scale_value = 2f64.powf(user_bound_scale as f64);
    if matrix.is_colwise() {
        for i_col in 0..num_col {
            if integrality[i_col as usize] == HighsVarType::kContinuous {
                continue;
            }
            let mut i_el = matrix.start_[i_col as usize];
            while i_el < matrix.start_[(i_col + 1) as usize] {
                let value = matrix.value_[i_el as usize] * bound_scale_value;
                let abs_value = value.abs();
                if abs_value <= data.small_matrix_value {
                    data.num_small_matrix_values += 1;
                } else if abs_value >= data.large_matrix_value {
                    data.num_large_matrix_values += 1;
                }
                if apply {
                    matrix.value_[i_el as usize] = value;
                }
                i_el += 1;
            }
        }
    } else {
        for i_row in 0..num_row {
            let mut i_el = matrix.start_[i_row as usize];
            while i_el < matrix.start_[(i_row + 1) as usize] {
                let i_col = matrix.index_[i_el as usize];
                if integrality[i_col as usize] == HighsVarType::kContinuous {
                    i_el += 1;
                    continue;
                }
                let value = matrix.value_[i_el as usize] * bound_scale_value;
                let abs_value = value.abs();
                if abs_value <= data.small_matrix_value {
                    data.num_small_matrix_values += 1;
                } else if abs_value >= data.large_matrix_value {
                    data.num_large_matrix_values += 1;
                }
                if apply {
                    matrix.value_[i_el as usize] = value;
                }
                i_el += 1;
            }
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::userScaleStatus (lines 911-924)
pub fn user_scale_status(log_options: &HighsLogOptions, data: &HighsUserScaleData) -> HighsStatus {
    let mut return_status = HighsStatus::kOk;
    let mut message = String::new();
    if data.scale_warning(&mut message) {
        crate::highs_log_user!(log_options, HighsLogType::kWarning, "{}\n", message);
        return_status = HighsStatus::kWarning;
    }
    if data.scale_error(&mut message) {
        crate::highs_log_user!(log_options, HighsLogType::kError, "{}\n", message);
        return_status = HighsStatus::kError;
    }
    return_status
}

impl Highs {
    /// HiGHS: HighsInterface.cpp::Highs::userScale (lines 3093-3112)
    ///
    /// With the default options (both scales 0) the opening guard (lines
    /// 3094-3095) returns at once.
    pub fn user_scale(&mut self, data: &mut HighsUserScaleData) -> HighsStatus {
        if self.options_.user_objective_scale == 0 && self.options_.user_bound_scale == 0 {
            return HighsStatus::kOk;
        }
        // User objective and bound scaling data are accumulated in the
        // HighsUserScaleData struct, in particular, there is a local copy
        // of the user objective and bound scaling options values, and
        // records of resulting extreme data values that prevent the user
        // objective and bound scaling from being applied.
        initialise_user_scale_data(&self.options_, data);
        // Determine whether user scaling yields excessively large cost,
        // Hessian values, column/row bounds or matrix values. If not,
        // then apply the user scaling to the model...
        if self.user_scale_model(data) == HighsStatus::kError {
            return HighsStatus::kError;
        }
        // ... and the solution
        self.user_scale_solution(data, false);
        // Indicate that the scaling has been applied
        data.applied = true;
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::userUnscale (lines 3114-3145)
    ///
    /// The opening guard (line 3115) is in place: `data.applied` is only set
    /// by a `user_scale` that scaled.
    pub fn user_unscale(&mut self, data: &mut HighsUserScaleData) -> HighsStatus {
        if !data.applied {
            return HighsStatus::kOk;
        }
        // Unscale the incumbent model and solution
        let mut status = HighsStatus::kOk;
        // Flip the scaling sign
        data.user_objective_scale *= -1;
        data.user_bound_scale *= -1;
        let mut unscale_status = self.user_scale_model(data);
        if unscale_status == HighsStatus::kError {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Unexpected error removing user scaling from the incumbent model\n"
            );
            debug_assert!(unscale_status != HighsStatus::kError);
        }
        let update_kkt = true;
        unscale_status = self.user_scale_solution(data, update_kkt);
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "After solving the user-scaled model, the unscaled solution has objective value {}\n",
            crate::util::c_printf::c_fmt_g(self.info_.objective_function_value, 12)
        );
        if self.model_status_ == HighsModelStatus::kOptimal && unscale_status != HighsStatus::kOk {
            // KKT errors in the unscaled optimal solution, so log a warning and
            // return
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "User scaled problem solved to optimality, but unscaled solution does not satisfy feasibility and optimality tolerances\n"
            );
            status = HighsStatus::kWarning;
        }
        status
    }

    /// HiGHS: HighsInterface.cpp::Highs::userScaleModel (lines 3147-3159).
    /// `userScaleHessian(model_.hessian_, data, apply)` ->
    /// `user_scale_hessian(&mut self.model_hessian_, data, apply)`.
    pub fn user_scale_model(&mut self, data: &mut HighsUserScaleData) -> HighsStatus {
        // Consider applying user objective and bound scaling to the model
        // by first identifying whether it causes any errors due to creating
        // extreme data values...
        let h16_user_scale_offset = self.options_.h16_user_scale_offset;
        user_scale_lp(&mut self.model_lp_, data, false, h16_user_scale_offset);
        user_scale_hessian(&mut self.model_hessian_, data, false);
        let return_status = user_scale_status(&self.options_.log_options, data);
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        // ... and, if not, actually apply the scaling
        user_scale_lp(&mut self.model_lp_, data, true, h16_user_scale_offset);
        user_scale_hessian(&mut self.model_hessian_, data, true);
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::userScaleSolution (lines 3161-3202),
    /// `update_kkt` defaulting to false in Highs.h. `getKktFailures(options_,
    /// model_, solution_, basis_, info_)` -> `get_kkt_failures_model`;
    /// `reportKktFailures` -> `report_kkt_failures` (presolve_support.rs).
    /// `std::pow(2, n)` -> `2f64.powf(n as f64)` (exact for every `HighsInt` n).
    /// HiGHS 1.16 (#3186): HighsInterface.cpp::Highs::userScaleSolution (lines
    /// 3189-3193) scales the objective with no offset round trip, and the MIP
    /// dual bound; gated by `h16_user_scale_offset`.
    pub fn user_scale_solution(
        &mut self,
        data: &mut HighsUserScaleData,
        update_kkt: bool,
    ) -> HighsStatus {
        let mut return_status = HighsStatus::kOk;
        if data.user_objective_scale == 0 && data.user_bound_scale == 0 {
            return HighsStatus::kOk;
        }
        let objective_scale_value = 2f64.powf(data.user_objective_scale as f64);
        let bound_scale_value = 2f64.powf(data.user_bound_scale as f64);
        let lp = &self.model_lp_;
        let has_integrality = !lp.integrality_.is_empty();
        if self.info_.primal_solution_status != kSolutionStatusNone {
            if data.user_bound_scale != 0 {
                for i_col in 0..lp.num_col_ {
                    if has_integrality
                        && lp.integrality_[i_col as usize] != HighsVarType::kContinuous
                    {
                        continue;
                    }
                    self.solution_.col_value[i_col as usize] *= bound_scale_value;
                }
                for i_row in 0..lp.num_row_ {
                    self.solution_.row_value[i_row as usize] *= bound_scale_value;
                }
            }
        }
        if self.info_.dual_solution_status != kSolutionStatusNone {
            if data.user_objective_scale != 0 {
                for i_col in 0..lp.num_col_ {
                    self.solution_.col_dual[i_col as usize] *= objective_scale_value;
                }
                for i_row in 0..lp.num_row_ {
                    self.solution_.row_dual[i_row as usize] *= objective_scale_value;
                }
            }
        }
        if !update_kkt {
            return return_status;
        }
        if self.options_.h16_user_scale_offset {
            // HiGHS 1.16 #3186 (h16_user_scale_offset): HighsInterface.cpp lines 3191-3193 (1.16.0).
            // The offset is scaled at apply time, so the objective is scaled without
            // the offset round trip. The dual bound is scaled for MIPs.
            self.info_.objective_function_value *= bound_scale_value * objective_scale_value;
            if has_integrality {
                self.info_.mip_dual_bound *= bound_scale_value * objective_scale_value;
            }
        } else {
            // In scaling the objective function value, have to consider the offset
            let mut objective_function_value =
                self.info_.objective_function_value - self.model_lp_.offset_;
            objective_function_value *= bound_scale_value * objective_scale_value;
            objective_function_value += self.model_lp_.offset_;
            self.info_.objective_function_value = objective_function_value;
        }
        let lp = &self.model_lp_;
        let hessian = &self.model_hessian_;
        let solution = &self.solution_;
        let basis = &self.basis_;
        let info = &mut self.info_;
        get_kkt_failures_model(&self.options_, lp, hessian, solution, basis, info);
        let lp = &self.model_lp_;
        let options = &self.options_;
        let info = &self.info_;
        if report_kkt_failures(lp, options, info, "After removing user scaling") {
            HighsStatus::kWarning
        } else {
            return_status
        }
    }

    /// HiGHS: Highs.cpp::Highs::getObjectiveBoundScaling (lines 2226-2236)
    pub fn get_objective_bound_scaling(
        &mut self,
        suggested_objective_scale: &mut HighsInt,
        suggested_bound_scale: &mut HighsInt,
    ) -> HighsStatus {
        // logHeader() dropped (logging only, line 2228)
        let mut data = HighsUserScaleData::default();
        initialise_user_scale_data(&self.options_, &mut data);
        // HiGHS 1.16 (h16_scaling_report): Highs.cpp line 2443 (1.16.0) makes the matrix
        // column-wise first
        if self.options_.h16_scaling_report {
            self.model_lp_.a_matrix_.ensure_colwise();
        }
        let lp = &self.model_lp_;
        assess_excessive_objective_bound_scaling(&self.options_, lp, &mut data);
        *suggested_objective_scale = data.suggested_user_objective_scale;
        *suggested_bound_scale = data.suggested_user_bound_scale;
        HighsStatus::kOk
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    fn data(obj: HighsInt, bnd: HighsInt) -> HighsUserScaleData {
        let mut d = HighsUserScaleData::default();
        d.initialise(obj, bnd, 1e20, 1e20, 1e-9, 1e15);
        d
    }

    #[test]
    fn initialise_sets_scales_and_zeroes_counts() {
        let mut d = HighsUserScaleData::default();
        d.num_infinite_costs = 3;
        d.applied = true;
        d.initialise(2, -1, 1e20, 1e21, 1e-9, 1e15);
        assert_eq!((d.user_objective_scale, d.user_bound_scale), (2, -1));
        assert_eq!((d.infinite_cost, d.infinite_bound), (1e20, 1e21));
        assert_eq!(d.num_infinite_costs, 0);
        assert!(!d.applied);
    }

    #[test]
    fn scale_error_message_text() {
        let mut d = data(3, 2);
        let mut m = String::new();
        assert!(!d.scale_error(&mut m));
        d.num_infinite_costs = 2;
        d.num_infinite_row_bounds = 1;
        assert!(d.scale_error(&mut m));
        assert_eq!(
            m,
            "User scaling of 2**(3) for costs and 2**(2) for bounds yields 2 infinite costs and 1 infinite row bound\n"
        );
    }

    #[test]
    fn scale_warning_message_text() {
        let mut d = data(0, -4);
        let mut m = String::new();
        assert!(!d.scale_warning(&mut m));
        d.num_small_matrix_values = 1;
        assert!(d.scale_warning(&mut m));
        assert_eq!(
            m,
            "User scaling of 2**(-4) for bounds yields 1 small matrix value\n"
        );
    }

    #[test]
    fn costs_scale_by_power_of_two() {
        let integrality: Vec<HighsVarType> = Vec::new();
        let mut cost = vec![1.5, -3.0, 0.0];
        let mut d = data(2, 0);
        let mut offset = 0.0;
        user_scale_costs(&integrality, &mut offset, &mut cost, &mut d, true, false);
        assert_eq!(cost, vec![6.0, -12.0, 0.0]);
        assert_eq!(d.num_infinite_costs, 0);
    }

    /// HiGHS 1.16 #3186: `offset` is scaled only on apply, only with the flag on.
    #[test]
    fn h16_uscale_costs_offset_scaled_on_apply() {
        let integrality: Vec<HighsVarType> = Vec::new();
        let mut cost = vec![1.5, -3.0];
        let mut d = data(1, 0);
        let mut offset = 3.0;
        user_scale_costs(&integrality, &mut offset, &mut cost, &mut d, false, true);
        assert_eq!(offset, 3.0);
        user_scale_costs(&integrality, &mut offset, &mut cost, &mut d, true, false);
        assert_eq!(offset, 3.0);
        user_scale_costs(&integrality, &mut offset, &mut cost, &mut d, true, true);
        assert_eq!(offset, 6.0);
        assert_eq!(cost, vec![6.0, -12.0]);
    }

    /// HiGHS 1.16 #3186 LP test model: min x + 2y + 3 (offset 3), s.t. x + y >= 2.
    /// Optimum x = 2, y = 0, objective 5. `integer` makes x integer and asks
    /// for x + y >= 2.5, optimum x = 3 (or x = 2, y = 0.5), objective 6.
    fn h16_uscale_lp(integer: bool) -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 1,
            col_cost_: vec![1.0, 2.0],
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![10.0, 10.0],
            row_lower_: vec![if integer { 2.5 } else { 2.0 }],
            row_upper_: vec![kHighsInf],
            sense_: ObjSense::kMinimize,
            offset_: 3.0,
            ..Default::default()
        };
        if integer {
            lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kContinuous];
        }
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        lp
    }

    /// Runs the HiGHS 1.16 #3186 LP with user objective scale 1 and the given flag.
    fn h16_uscale_run(integer: bool, h16_user_scale_offset: bool) -> Highs {
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
            highs.set_option_value_bool("h16_user_scale_offset", h16_user_scale_offset),
            HighsStatus::kOk
        );
        highs.options_.user_objective_scale = 1;
        assert_eq!(highs.pass_model(h16_uscale_lp(integer)), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        highs
    }

    /// HiGHS 1.16 #3186 on path: the reported LP objective is the unscaled objective
    /// including the offset (the offset is scaled at apply time and unscaled on
    /// removal). With the flag on, the scaled objective is 2*(x + 2y) + 2*3 = 10.
    #[test]
    fn h16_uscale_lp_objective_on() {
        let highs = h16_uscale_run(false, true);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        let objective = highs.get_info().objective_function_value;
        assert!((objective - 5.0).abs() < 1e-9, "on: {objective}");
    }

    /// HiGHS 1.16 #3186 off path: the 1.15.1 offset round trip gives the same LP
    /// objective, 5.
    #[test]
    fn h16_uscale_lp_objective_off() {
        let highs = h16_uscale_run(false, false);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        let objective = highs.get_info().objective_function_value;
        assert!((objective - 5.0).abs() < 1e-9, "off: {objective}");
    }

    /// HiGHS 1.16 #3186 MIP: with the flag on, mip_dual_bound is scaled with the
    /// objective, so the gap is closed at the optimum (objective 6).
    #[test]
    fn h16_uscale_mip_dual_bound_on() {
        let highs = h16_uscale_run(true, true);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        let info = highs.get_info();
        assert!(
            (info.objective_function_value - 6.0).abs() < 1e-6,
            "objective {}",
            info.objective_function_value
        );
        assert!(
            (info.mip_dual_bound - info.objective_function_value).abs() < 1e-6,
            "dual bound {} vs objective {}",
            info.mip_dual_bound,
            info.objective_function_value
        );
    }

    /// HiGHS 1.16 #3186 MIP off path (1.15.1): the offset is not scaled and
    /// mip_dual_bound is left in scaled units, 2 * 3 + 3 = 9, while the
    /// objective is unscaled through the offset round trip to 6.
    #[test]
    fn h16_uscale_mip_dual_bound_off() {
        let highs = h16_uscale_run(true, false);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        let info = highs.get_info();
        assert!(
            (info.objective_function_value - 6.0).abs() < 1e-6,
            "objective {}",
            info.objective_function_value
        );
        assert!(
            (info.mip_dual_bound - 9.0).abs() < 1e-6,
            "dual bound {}",
            info.mip_dual_bound
        );
    }

    #[test]
    fn user_scale_is_a_no_op_with_default_options() {
        let mut highs = Highs::default();
        let mut d = HighsUserScaleData::default();
        assert_eq!(highs.user_scale(&mut d), HighsStatus::kOk);
        assert!(!d.applied);
        assert_eq!(highs.user_unscale(&mut d), HighsStatus::kOk);
    }
}
