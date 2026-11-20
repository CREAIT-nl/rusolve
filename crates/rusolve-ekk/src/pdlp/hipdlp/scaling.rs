//! HiGHS source: pdlp/hipdlp/scaling.hpp, pdlp/hipdlp/scaling.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Ruiz / Pock-Chambolle / L2 scaling of the processed LP. Back-pointers
//! are dropped: `lp_` and `params_` become the
//! `lp` / `params` parameters; `highsLogDev` / `highsLogUser` go through
//! `params.log_options_`.

use super::defs::*;
use super::linalg;
use crate::prelude::*;

/// HiGHS: scaling.hpp `class Scaling`
#[derive(Clone, Debug, Default)]
pub struct Scaling {
    pub(crate) col_scale_: Vec<f64>,
    pub(crate) row_scale_: Vec<f64>,
    pub(crate) is_scaled_: bool,
    pub(crate) norm_cost_: f64,
    pub(crate) norm_rhs_: f64,
}

impl Scaling {
    /// HiGHS: scaling.cc::Scaling::initialize (lines 21-29)
    pub fn initialize(&mut self, lp: &HighsLp) {
        self.col_scale_ = vec![1.0; lp.num_col_ as usize];
        self.row_scale_ = vec![1.0; lp.num_row_ as usize];
        self.is_scaled_ = false;

        // use linalg to compute norms
        self.norm_cost_ = linalg::computeCostNorm(lp, 2.0);
        self.norm_rhs_ = linalg::computeRhsNorm(lp, 2.0);
    }

    /// HiGHS: scaling.cc::Scaling::scaleProblem (lines 31-54)
    pub fn scaleProblem(&mut self, lp: &mut HighsLp, params: &PrimalDualParams) {
        self.is_scaled_ = false;

        if params.use_ruiz_scaling {
            highs_log_dev(
                &params.log_options_,
                HighsLogType::kInfo,
                "Applying Ruiz scaling...\n",
            );
            self.applyRuizScaling(lp, params);
            self.is_scaled_ = true;
        }

        if params.use_pc_scaling {
            highs_log_dev(
                &params.log_options_,
                HighsLogType::kInfo,
                "Applying Pock-Chambolle scaling...\n",
            );
            self.applyPockChambolleScaling(lp, params);
            self.is_scaled_ = true;
        }

        if params.use_l2_scaling {
            highs_log_dev(
                &params.log_options_,
                HighsLogType::kInfo,
                "Applying L2 norm scaling...\n",
            );
            self.applyL2Scaling(lp);
            self.is_scaled_ = true;
        }
    }

    /// HiGHS: scaling.cc::Scaling::unscaleSolution (lines 270-283)
    pub fn unscaleSolution(&self, x: &mut [f64], y: &mut [f64]) {
        if !self.is_scaled_ {
            return;
        }

        // Unscale primal variables: x_original = x_scaled / col_scale
        for i in 0..x.len() {
            x[i] /= self.col_scale_[i];
        }

        // Unscale dual variables: y_original = y_scaled / row_scale
        for i in 0..y.len() {
            y[i] /= self.row_scale_[i];
        }
    }

    /// HiGHS: scaling.hpp::Scaling::isScaled
    pub fn isScaled(&self) -> bool {
        self.is_scaled_
    }

    /// HiGHS: scaling.hpp::Scaling::getColScaling
    pub fn getColScaling(&self) -> &[f64] {
        &self.col_scale_
    }

    /// HiGHS: scaling.hpp::Scaling::getRowScaling
    pub fn getRowScaling(&self) -> &[f64] {
        &self.row_scale_
    }

    /// HiGHS: scaling.hpp::Scaling::getNormCost
    pub fn getNormCost(&self) -> f64 {
        self.norm_cost_
    }

    /// HiGHS: scaling.hpp::Scaling::getNormRhs
    pub fn getNormRhs(&self) -> f64 {
        self.norm_rhs_
    }

    /// HiGHS: scaling.cc::Scaling::applyRuizScaling (lines 56-122)
    pub(crate) fn applyRuizScaling(&mut self, lp: &mut HighsLp, params: &PrimalDualParams) {
        let mut current_col_scaling = vec![0.0f64; lp.num_col_ as usize];
        let mut current_row_scaling = vec![0.0f64; lp.num_row_ as usize];

        let mut iter: HighsInt = 0;
        while iter < params.ruiz_iterations {
            // Reset current scaling factors
            current_col_scaling.fill(0.0);
            current_row_scaling.fill(0.0);

            // Compute column norms (norm of each column)
            let mut col: HighsInt = 0;
            while col < lp.num_col_ {
                let mut col_values: Vec<f64> = Vec::new();
                let mut el: HighsInt = lp.a_matrix_.start_[col as usize];
                while el < lp.a_matrix_.start_[(col + 1) as usize] {
                    col_values.push(lp.a_matrix_.value_[el as usize]);
                    el += 1;
                }

                if !col_values.is_empty() {
                    current_col_scaling[col as usize] = self
                        .computeNorm(&col_values, col_values.len() as HighsInt, params.ruiz_norm)
                        .sqrt();
                }

                if current_col_scaling[col as usize] == 0.0 {
                    current_col_scaling[col as usize] = 1.0;
                }
                col += 1;
            }

            // Compute row norms (infinity norm for rows)
            let ruiz_norm_ok = params.ruiz_norm == f64::INFINITY;
            if !ruiz_norm_ok {
                highs_log_user(
                    &params.log_options_,
                    HighsLogType::kError,
                    "Currently only support infinity norm for Ruiz scaling\n",
                );
                debug_assert!(ruiz_norm_ok);
            }
            if params.ruiz_norm == f64::INFINITY {
                // For infinity norm, find max absolute value in each row
                let mut col: HighsInt = 0;
                while col < lp.num_col_ {
                    let mut el: HighsInt = lp.a_matrix_.start_[col as usize];
                    while el < lp.a_matrix_.start_[(col + 1) as usize] {
                        let row = lp.a_matrix_.index_[el as usize];
                        let abs_val = lp.a_matrix_.value_[el as usize].abs();
                        current_row_scaling[row as usize] =
                            cpp_max(current_row_scaling[row as usize], abs_val);
                        el += 1;
                    }
                    col += 1;
                }

                let mut row: HighsInt = 0;
                while row < lp.num_row_ {
                    if current_row_scaling[row as usize] == 0.0 {
                        current_row_scaling[row as usize] = 1.0;
                    } else {
                        current_row_scaling[row as usize] =
                            current_row_scaling[row as usize].sqrt();
                    }
                    row += 1;
                }
            }

            // Apply the scaling
            Scaling::applyScaling(lp, &current_col_scaling, &current_row_scaling);

            // Update cumulative scaling factors
            let mut i: HighsInt = 0;
            while i < lp.num_col_ {
                self.col_scale_[i as usize] *= current_col_scaling[i as usize];
                i += 1;
            }
            let mut i: HighsInt = 0;
            while i < lp.num_row_ {
                self.row_scale_[i as usize] *= current_row_scaling[i as usize];
                i += 1;
            }
            iter += 1;
        }
    }

    /// HiGHS: scaling.cc::Scaling::applyPockChambolleScaling (lines 124-180)
    pub(crate) fn applyPockChambolleScaling(
        &mut self,
        lp: &mut HighsLp,
        params: &PrimalDualParams,
    ) {
        let pc_alpha_ok = 0.0 <= params.pc_alpha && params.pc_alpha <= 2.0;
        if !pc_alpha_ok {
            highs_log_user(
                &params.log_options_,
                HighsLogType::kError,
                "PC alpha must be in [0, 2]\n",
            );
            debug_assert!(pc_alpha_ok);
        }

        let mut current_col_scaling = vec![0.0f64; lp.num_col_ as usize];
        let mut current_row_scaling = vec![0.0f64; lp.num_row_ as usize];

        // Compute column scaling: (sum |A_ij|^alpha)^(1/alpha)
        let mut col: HighsInt = 0;
        while col < lp.num_col_ {
            let mut el: HighsInt = lp.a_matrix_.start_[col as usize];
            while el < lp.a_matrix_.start_[(col + 1) as usize] {
                current_col_scaling[col as usize] +=
                    lp.a_matrix_.value_[el as usize].abs().powf(params.pc_alpha);
                el += 1;
            }

            if current_col_scaling[col as usize] > 0.0 {
                current_col_scaling[col as usize] = current_col_scaling[col as usize]
                    .powf(1.0 / params.pc_alpha)
                    .sqrt();
            } else {
                current_col_scaling[col as usize] = 1.0;
            }
            col += 1;
        }

        // Compute row scaling: (sum |A_ij|^(2-alpha))^(1/(2-alpha))
        let mut col: HighsInt = 0;
        while col < lp.num_col_ {
            let mut el: HighsInt = lp.a_matrix_.start_[col as usize];
            while el < lp.a_matrix_.start_[(col + 1) as usize] {
                let row = lp.a_matrix_.index_[el as usize];
                current_row_scaling[row as usize] += lp.a_matrix_.value_[el as usize]
                    .abs()
                    .powf(2.0 - params.pc_alpha);
                el += 1;
            }
            col += 1;
        }

        let mut row: HighsInt = 0;
        while row < lp.num_row_ {
            if current_row_scaling[row as usize] > 0.0 {
                current_row_scaling[row as usize] = current_row_scaling[row as usize]
                    .powf(1.0 / (2.0 - params.pc_alpha))
                    .sqrt();
            } else {
                current_row_scaling[row as usize] = 1.0;
            }
            row += 1;
        }

        // Apply the scaling
        Scaling::applyScaling(lp, &current_col_scaling, &current_row_scaling);

        // Update cumulative scaling factors
        let mut i: HighsInt = 0;
        while i < lp.num_col_ {
            self.col_scale_[i as usize] *= current_col_scaling[i as usize];
            i += 1;
        }
        let mut i: HighsInt = 0;
        while i < lp.num_row_ {
            self.row_scale_[i as usize] *= current_row_scaling[i as usize];
            i += 1;
        }
    }

    /// HiGHS: scaling.cc::Scaling::applyL2Scaling (lines 182-229)
    pub(crate) fn applyL2Scaling(&mut self, lp: &mut HighsLp) {
        let mut current_col_scaling = vec![1.0f64; lp.num_col_ as usize];
        let mut current_row_scaling = vec![0.0f64; lp.num_row_ as usize];

        // Compute L2 norm of each column
        let mut col: HighsInt = 0;
        while col < lp.num_col_ {
            let mut sum_sq: f64 = 0.0;
            let mut el: HighsInt = lp.a_matrix_.start_[col as usize];
            while el < lp.a_matrix_.start_[(col + 1) as usize] {
                sum_sq += lp.a_matrix_.value_[el as usize] * lp.a_matrix_.value_[el as usize];
                el += 1;
            }

            if sum_sq > 0.0 {
                current_col_scaling[col as usize] = sum_sq.sqrt().sqrt();
            } else {
                current_col_scaling[col as usize] = 1.0;
            }
            col += 1;
        }

        // Compute L2 norm of each row
        let mut col: HighsInt = 0;
        while col < lp.num_col_ {
            let mut el: HighsInt = lp.a_matrix_.start_[col as usize];
            while el < lp.a_matrix_.start_[(col + 1) as usize] {
                let row = lp.a_matrix_.index_[el as usize];
                current_row_scaling[row as usize] +=
                    lp.a_matrix_.value_[el as usize] * lp.a_matrix_.value_[el as usize];
                el += 1;
            }
            col += 1;
        }

        let mut row: HighsInt = 0;
        while row < lp.num_row_ {
            if current_row_scaling[row as usize] > 0.0 {
                current_row_scaling[row as usize] = current_row_scaling[row as usize].sqrt().sqrt();
            } else {
                current_row_scaling[row as usize] = 1.0;
            }
            row += 1;
        }

        // Apply the scaling
        Scaling::applyScaling(lp, &current_col_scaling, &current_row_scaling);

        // Update cumulative scaling factors
        let mut i: HighsInt = 0;
        while i < lp.num_col_ {
            self.col_scale_[i as usize] *= current_col_scaling[i as usize];
            i += 1;
        }
        let mut i: HighsInt = 0;
        while i < lp.num_row_ {
            self.row_scale_[i as usize] *= current_row_scaling[i as usize];
            i += 1;
        }
    }

    /// HiGHS: scaling.cc::Scaling::applyScaling (lines 231-268)
    pub(crate) fn applyScaling(lp: &mut HighsLp, col_scaling: &[f64], row_scaling: &[f64]) {
        // Scale cost vector: c_scaled = c / col_scaling
        let mut i: HighsInt = 0;
        while i < lp.num_col_ {
            lp.col_cost_[i as usize] /= col_scaling[i as usize];
            i += 1;
        }

        // Scale column bounds: l_scaled = l * col_scaling, u_scaled = u * col_scaling
        let mut i: HighsInt = 0;
        while i < lp.num_col_ {
            if lp.col_lower_[i as usize] > -kHighsInf {
                lp.col_lower_[i as usize] *= col_scaling[i as usize];
            }
            if lp.col_upper_[i as usize] < kHighsInf {
                lp.col_upper_[i as usize] *= col_scaling[i as usize];
            }
            i += 1;
        }

        // Scale row bounds: b_scaled = b / row_scaling
        let mut i: HighsInt = 0;
        while i < lp.num_row_ {
            if lp.row_lower_[i as usize] > -kHighsInf {
                lp.row_lower_[i as usize] /= row_scaling[i as usize];
            }
            if lp.row_upper_[i as usize] < kHighsInf {
                lp.row_upper_[i as usize] /= row_scaling[i as usize];
            }
            i += 1;
        }

        // Scale matrix: A_scaled = diag(1/row_scaling) * A * diag(1/col_scaling)
        let mut col: HighsInt = 0;
        while col < lp.num_col_ {
            let mut el: HighsInt = lp.a_matrix_.start_[col as usize];
            while el < lp.a_matrix_.start_[(col + 1) as usize] {
                let row = lp.a_matrix_.index_[el as usize];
                lp.a_matrix_.value_[el as usize] /=
                    row_scaling[row as usize] * col_scaling[col as usize];
                el += 1;
            }
            col += 1;
        }

        // TO DO: Scale offset
    }

    /// HiGHS: scaling.cc::Scaling::computeNorm (lines 285-307)
    pub(crate) fn computeNorm(&self, values: &[f64], size: HighsInt, norm_type: f64) -> f64 {
        if norm_type == f64::INFINITY {
            let mut max_val: f64 = 0.0;
            let mut i: HighsInt = 0;
            while i < size {
                max_val = cpp_max(max_val, values[i as usize].abs());
                i += 1;
            }
            max_val
        } else if norm_type == 2.0 {
            let mut sum_sq: f64 = 0.0;
            let mut i: HighsInt = 0;
            while i < size {
                sum_sq += values[i as usize] * values[i as usize];
                i += 1;
            }
            sum_sq.sqrt()
        } else {
            // General p-norm
            let mut sum: f64 = 0.0;
            let mut i: HighsInt = 0;
            while i < size {
                sum += values[i as usize].abs().powf(norm_type);
                i += 1;
            }
            sum.powf(1.0 / norm_type)
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn compute_norm_and_unscale() {
        let s = Scaling {
            col_scale_: vec![2.0],
            row_scale_: vec![4.0],
            is_scaled_: true,
            ..Default::default()
        };
        assert_eq!(s.computeNorm(&[3.0, -4.0], 2, f64::INFINITY), 4.0);
        assert_eq!(s.computeNorm(&[3.0, -4.0], 2, 2.0), 5.0);
        let (mut x, mut y) = (vec![6.0], vec![8.0]);
        s.unscaleSolution(&mut x, &mut y);
        assert_eq!((x[0], y[0]), (3.0, 2.0));
    }
}
