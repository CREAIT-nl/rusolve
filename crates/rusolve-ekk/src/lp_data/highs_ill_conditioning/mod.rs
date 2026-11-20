//! HiGHS source: lp_data/Highs.cpp (getIllConditioning),
//! lp_data/HighsInterface.cpp (computeIllConditioning, formIllConditioningLp0;
//! formIllConditioningLp1 is in lp1.rs)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The ill-conditioning analysis of the
//! current basis. An auxiliary `Highs conditioning` (output off) solves an
//! LP built from the basis matrix (method 0: `formIllConditioningLp0`;
//! method 1: `formIllConditioningLp1` with `ill_conditioning_bound`); the
//! normalised multipliers above `kZeroMultiplier = 1e-6` become the
//! `HighsIllConditioning` records, reported with `highsLogUser`.
//! The `dev_conditioning` (always false) writes are dropped.

pub mod lp1;

use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::getIllConditioning (lines 2213-2224).
    /// C++ default arguments: `method = 0`, `ill_conditioning_bound = 1e-4`.
    pub fn get_ill_conditioning(
        &mut self,
        ill_conditioning: &mut HighsIllConditioning,
        constraint: bool,
        method: HighsInt,
        ill_conditioning_bound: f64,
    ) -> HighsStatus {
        if !self.basis_.valid {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot get ill-conditioning without a valid basis\n"
            );
            return HighsStatus::kError;
        }
        self.compute_ill_conditioning(ill_conditioning, constraint, method, ill_conditioning_bound)
    }

    /// HiGHS: HighsInterface.cpp::Highs::computeIllConditioning (lines 3206-3404).
    /// `Highs conditioning;` -> `let mut conditioning = Highs::default()`, then
    /// `set_option_value("output_flag", false)`; `assert(assessLp(...))` ->
    /// `debug_assert!`; `conditioning.run()` keeps its position.
    pub fn compute_ill_conditioning(
        &mut self,
        ill_conditioning: &mut HighsIllConditioning,
        constraint: bool,
        method: HighsInt,
        ill_conditioning_bound: f64,
    ) -> HighsStatus {
        const k_zero_multiplier: f64 = 1e-6;
        ill_conditioning.clear();
        // HighsLp& incumbent_lp = this->model_.lp_;
        let mut conditioning = Highs::default();
        // const bool dev_conditioning = false;
        // conditioning.setOptionValue("output_flag", false);
        conditioning.set_option_value_bool("output_flag", false);
        let mut basic_var: Vec<HighsInt> = Vec::new();
        // HighsLp& ill_conditioning_lp = conditioning.model_.lp_;
        // Form the ill-conditioning LP according to method
        if method == 0 {
            self.form_ill_conditioning_lp0(&mut conditioning.model_lp_, &mut basic_var, constraint);
        } else {
            self.form_ill_conditioning_lp1(
                &mut conditioning.model_lp_,
                &mut basic_var,
                constraint,
                ill_conditioning_bound,
            );
            // conditioning.writeModel(""); // dropped (writeModel)
        }

        // if (dev_conditioning) conditioning.writeModel(""); // dropped (dev_conditioning always false)

        debug_assert!(
            crate::lp_data::lp_assess::assess_lp(&mut conditioning.model_lp_, &self.options_)
                == HighsStatus::kOk
        );
        // Solve the ill-conditioning analysis LP
        let return_status = conditioning.run();
        let model_status = conditioning.get_model_status();
        let type_str = if constraint { "Constraint" } else { "Column" };
        let failed = return_status != HighsStatus::kOk
            || (method == 0 && model_status != HighsModelStatus::kOptimal)
            || (method == 1
                && (model_status != HighsModelStatus::kOptimal
                    && model_status != HighsModelStatus::kInfeasible));
        if failed {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "\n{} view ill-conditioning analysis has failed\n",
                type_str
            );
            return HighsStatus::kError;
        }
        if method == 1 && model_status == HighsModelStatus::kInfeasible {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "\n{} view ill-conditioning bound of {} is insufficient for analysis: try {}\n",
                type_str,
                crate::util::c_printf::c_fmt_g(ill_conditioning_bound, 6),
                crate::util::c_printf::c_fmt_g(1e1 * ill_conditioning_bound, 6)
            );
            return HighsStatus::kOk;
        }
        // if (dev_conditioning) conditioning.writeSolution("", 1); // dropped (dev_conditioning always false)
        // Extract and normalise the multipliers
        // HighsSolution& solution = conditioning.solution_;
        let mut multiplier_norm: f64 = 0.0;
        for i_row in 0..self.model_lp_.num_row_ {
            multiplier_norm += conditioning.solution_.col_value[i_row as usize].abs();
        }
        debug_assert!(multiplier_norm > 0.0);
        let ill_conditioning_measure = if method == 0 {
            conditioning.get_info().objective_function_value
        } else {
            conditioning.solution_.row_value[(conditioning.get_num_row() - 1) as usize]
        } / multiplier_norm;
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "\n{} view ill-conditioning analysis: 1-norm distance of basis matrix from singularity is estimated to be {}\n",
            type_str,
            crate::util::c_printf::c_fmt_g(ill_conditioning_measure, 6)
        );
        let mut abs_list: Vec<(f64, HighsInt)> = Vec::new();
        for i_row in 0..self.model_lp_.num_row_ {
            let abs_multiplier =
                conditioning.solution_.col_value[i_row as usize].abs() / multiplier_norm;
            if abs_multiplier <= k_zero_multiplier {
                continue;
            }
            abs_list.push((abs_multiplier, i_row));
        }
        // PORT NOTE: std::sort, keys distinct
        abs_list.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        // Report on ill-conditioning multipliers
        let mut ss = String::new();
        let has_row_names = self.model_lp_.row_names_.len() as HighsInt == self.model_lp_.num_row_;
        let has_col_names = self.model_lp_.col_names_.len() as HighsInt == self.model_lp_.num_col_;
        const coefficient_zero_tolerance: f64 = 1e-8;

        // HiGHS: HighsInterface.cpp::Highs::computeIllConditioning (lines 3283-3300)
        // auto printCoefficient = [&](const double multiplier, const bool first) { ... };
        let print_coefficient = |ss: &mut String, multiplier: f64, first: bool| {
            if multiplier.abs() < coefficient_zero_tolerance {
                ss.push_str("+ 0");
            } else if (multiplier - 1.0).abs() < coefficient_zero_tolerance {
                let str = if first { "" } else { "+ " };
                ss.push_str(str);
            } else if (multiplier + 1.0).abs() < coefficient_zero_tolerance {
                let str = if first { "-" } else { "- " };
                ss.push_str(str);
            } else if multiplier < 0.0 {
                let str = if first { "-" } else { "- " };
                ss.push_str(str);
                ss.push_str(&crate::util::c_printf::c_fmt_g(-multiplier, 6));
                ss.push(' ');
            } else {
                let str = if first { "" } else { "+ " };
                ss.push_str(str);
                ss.push_str(&crate::util::c_printf::c_fmt_g(multiplier, 6));
                ss.push(' ');
            }
        };

        let mut i_x = abs_list.len() as HighsInt - 1;
        while i_x >= 0 {
            let i_row = abs_list[i_x as usize].1;
            let mut record = HighsIllConditioningRecord::default();
            record.index = i_row;
            record.multiplier = conditioning.solution_.col_value[i_row as usize] / multiplier_norm;
            ill_conditioning.record.push(record);
            i_x -= 1;
        }
        // HighsSparseMatrix& incumbent_matrix = incumbent_lp.a_matrix_;
        if constraint {
            let mut num_nz: HighsInt;
            let mut index: Vec<HighsInt> = vec![0; self.model_lp_.num_col_ as usize];
            let mut value: Vec<f64> = vec![0.0; self.model_lp_.num_col_ as usize];
            for i_x in 0..ill_conditioning.record.len() as HighsInt {
                ss.clear();
                let mut newline = false;
                let i_row = ill_conditioning.record[i_x as usize].index;
                let multiplier = ill_conditioning.record[i_x as usize].multiplier;
                // Extract the row corresponding to this constraint
                num_nz = 0;
                self.model_lp_.a_matrix_.get_row(
                    i_row,
                    &mut num_nz,
                    Some(&mut index),
                    Some(&mut value),
                );
                let row_name = if has_row_names {
                    self.model_lp_.row_names_[i_row as usize].clone()
                } else {
                    format!("R{}", i_row)
                };
                ss.push_str(&format!(
                    "(Mu={})",
                    crate::util::c_printf::c_fmt_g(multiplier, 6)
                ));
                ss.push_str(&row_name);
                ss.push_str(": ");
                let lower = self.model_lp_.row_lower_[i_row as usize];
                let upper = self.model_lp_.row_upper_[i_row as usize];
                if lower > -kHighsInf && lower != upper {
                    ss.push_str(&crate::util::c_printf::c_fmt_g(
                        self.model_lp_.row_lower_[i_row as usize],
                        6,
                    ));
                    ss.push_str(" <= ");
                }
                for i_el in 0..num_nz {
                    if newline {
                        ss.push_str("  ");
                        newline = false;
                    }
                    let i_col = index[i_el as usize];
                    print_coefficient(&mut ss, value[i_el as usize], i_el == 0);
                    let col_name = if has_col_names {
                        self.model_lp_.col_names_[i_col as usize].clone()
                    } else {
                        format!("C{}", i_col)
                    };
                    ss.push_str(&col_name);
                    ss.push(' ');
                    let length_ss = ss.len() as HighsInt;
                    if length_ss > 72 && i_el < num_nz - 1 {
                        crate::highs_log_user!(
                            &self.options_.log_options,
                            HighsLogType::kInfo,
                            "{}\n",
                            ss
                        );
                        ss.clear();
                        newline = true;
                    }
                }
                if upper < kHighsInf {
                    if lower == upper {
                        ss.push_str("= ");
                        ss.push_str(&crate::util::c_printf::c_fmt_g(upper, 6));
                    } else {
                        ss.push_str("<= ");
                        ss.push_str(&crate::util::c_printf::c_fmt_g(upper, 6));
                    }
                }
                if !ss.is_empty() {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "{}\n",
                        ss
                    );
                }
            }
        } else {
            for rec in &ill_conditioning.record {
                ss.clear();
                let mut newline = false;
                let multiplier = rec.multiplier;
                let i_col = basic_var[rec.index as usize];
                if i_col < self.model_lp_.num_col_ {
                    let col_name = if has_col_names {
                        self.model_lp_.col_names_[i_col as usize].clone()
                    } else {
                        format!("C{}", i_col)
                    };
                    ss.push_str(&format!(
                        "(Mu={})",
                        crate::util::c_printf::c_fmt_g(multiplier, 6)
                    ));
                    ss.push_str(&col_name);
                    ss.push_str(": ");
                    for i_el in self.model_lp_.a_matrix_.start_[i_col as usize]
                        ..self.model_lp_.a_matrix_.start_[(i_col + 1) as usize]
                    {
                        if newline {
                            ss.push_str("  ");
                            newline = false;
                        } else {
                            if i_el > self.model_lp_.a_matrix_.start_[i_col as usize] {
                                ss.push_str(" | ");
                            }
                        }
                        let i_row = self.model_lp_.a_matrix_.index_[i_el as usize];
                        print_coefficient(
                            &mut ss,
                            self.model_lp_.a_matrix_.value_[i_el as usize],
                            true,
                        );
                        let row_name = if has_row_names {
                            self.model_lp_.row_names_[i_row as usize].clone()
                        } else {
                            format!("R{}", i_row)
                        };
                        ss.push_str(&row_name);
                        let length_ss = ss.len() as HighsInt;
                        if length_ss > 72
                            && i_el < self.model_lp_.a_matrix_.start_[(i_col + 1) as usize] - 1
                        {
                            ss.push_str(" | ");
                            crate::highs_log_user!(
                                &self.options_.log_options,
                                HighsLogType::kInfo,
                                "{}\n",
                                ss
                            );
                            ss.clear();
                            newline = true;
                        }
                    }
                } else {
                    let i_row = i_col - self.model_lp_.num_col_;
                    let col_name = if has_row_names {
                        format!("Slack_{}", self.model_lp_.row_names_[i_row as usize])
                    } else {
                        format!("Slack_R{}", i_row)
                    };
                    ss.push_str(&format!(
                        "(Mu={})",
                        crate::util::c_printf::c_fmt_g(multiplier, 6)
                    ));
                    ss.push_str(&col_name);
                    ss.push_str(": ");
                }
                if !ss.is_empty() {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        "{}\n",
                        ss
                    );
                }
            }
        }
        HighsStatus::kOk
    }

    /// HiGHS: HighsInterface.cpp::Highs::formIllConditioningLp0 (lines 3406-3508).
    /// `ill_conditioning_lp` is `conditioning.model_.lp_` of the auxiliary
    /// instance, so it is a separate `&mut HighsLp`.
    pub fn form_ill_conditioning_lp0(
        &mut self,
        ill_conditioning_lp: &mut HighsLp,
        basic_var: &mut Vec<HighsInt>,
        constraint: bool,
    ) {
        // HighsLp& incumbent_lp = this->model_.lp_;
        // Conditioning LP minimizes the infeasibilities of
        //
        // [B^T]y = [0]; y free - for constraint view
        // [e^T]    [1]
        //
        // [ B ]y = [0]; y free - for column view
        // [e^T]    [1]
        //
        ill_conditioning_lp.num_row_ = self.model_lp_.num_row_ + 1;
        for i_row in 0..self.model_lp_.num_row_ {
            ill_conditioning_lp.row_lower_.push(0.0);
            ill_conditioning_lp.row_upper_.push(0.0);
        }
        ill_conditioning_lp.row_lower_.push(1.0);
        ill_conditioning_lp.row_upper_.push(1.0);
        // HighsSparseMatrix& incumbent_matrix = incumbent_lp.a_matrix_;
        self.model_lp_.a_matrix_.ensure_colwise();
        // HighsSparseMatrix& ill_conditioning_matrix = ill_conditioning_lp.a_matrix_;
        ill_conditioning_lp.a_matrix_.num_row_ = ill_conditioning_lp.num_row_;
        // Form the basis matrix and
        //
        // * For constraint view, add the column e, and transpose the
        // * resulting matrix
        //
        // * For column view, add a unit entry to each column
        //
        let ill_conditioning_lp_e_row = ill_conditioning_lp.num_row_ - 1;
        for i_col in 0..self.model_lp_.num_col_ {
            if self.basis_.col_status[i_col as usize] != HighsBasisStatus::kBasic {
                continue;
            }
            // Basic column goes into conditioning LP, possibly with unit
            // coefficient for constraint e^Ty=1
            basic_var.push(i_col);
            ill_conditioning_lp.col_cost_.push(0.0);
            ill_conditioning_lp.col_lower_.push(-kHighsInf);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
            for i_el in self.model_lp_.a_matrix_.start_[i_col as usize]
                ..self.model_lp_.a_matrix_.start_[(i_col + 1) as usize]
            {
                ill_matrix
                    .index_
                    .push(self.model_lp_.a_matrix_.index_[i_el as usize]);
                ill_matrix
                    .value_
                    .push(self.model_lp_.a_matrix_.value_[i_el as usize]);
            }
            if !constraint {
                ill_matrix.index_.push(ill_conditioning_lp_e_row);
                ill_matrix.value_.push(1.0);
            }
            ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
        }
        for i_row in 0..self.model_lp_.num_row_ {
            if self.basis_.row_status[i_row as usize] != HighsBasisStatus::kBasic {
                continue;
            }
            // Basic slack goes into conditioning LP
            basic_var.push(self.model_lp_.num_col_ + i_row);
            ill_conditioning_lp.col_cost_.push(0.0);
            ill_conditioning_lp.col_lower_.push(-kHighsInf);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
            ill_matrix.index_.push(i_row);
            ill_matrix.value_.push(-1.0);
            if !constraint {
                ill_matrix.index_.push(ill_conditioning_lp_e_row);
                ill_matrix.value_.push(1.0);
            }
            ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
        }
        if constraint {
            // Add the column e, and transpose the resulting matrix
            for i_row in 0..self.model_lp_.num_row_ {
                let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
                ill_matrix.index_.push(i_row);
                ill_matrix.value_.push(1.0);
            }
            let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
            ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            ill_matrix.num_row_ = self.model_lp_.num_row_;
            ill_matrix.num_col_ = self.model_lp_.num_row_ + 1;
            ill_matrix.ensure_rowwise();
            ill_matrix.format_ = MatrixFormat::kColwise;
        }
        // Now add the variables to measure the infeasibilities
        for i_row in 0..self.model_lp_.num_row_ {
            // Adding x_+ with cost 1
            ill_conditioning_lp.col_cost_.push(1.0);
            ill_conditioning_lp.col_lower_.push(0.0);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
            ill_matrix.index_.push(i_row);
            ill_matrix.value_.push(1.0);
            ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
            // Subtracting x_- with cost 1
            ill_conditioning_lp.col_cost_.push(1.0);
            ill_conditioning_lp.col_lower_.push(0.0);
            ill_conditioning_lp.col_upper_.push(kHighsInf);
            let ill_matrix = &mut ill_conditioning_lp.a_matrix_;
            ill_matrix.index_.push(i_row);
            ill_matrix.value_.push(-1.0);
            ill_matrix.start_.push(ill_matrix.index_.len() as HighsInt);
        }
        ill_conditioning_lp.num_col_ = 3 * self.model_lp_.num_row_;
        ill_conditioning_lp.a_matrix_.num_col_ = ill_conditioning_lp.num_col_;
        ill_conditioning_lp.a_matrix_.num_row_ = ill_conditioning_lp.num_row_;
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn ill_conditioning_needs_a_valid_basis() {
        let mut highs = Highs::default();
        let mut ill = HighsIllConditioning::default();
        ill.record.push(HighsIllConditioningRecord {
            index: 0,
            multiplier: 1.0,
        });
        assert_eq!(
            highs.get_ill_conditioning(&mut ill, true, 0, 1e-4),
            HighsStatus::kError
        );
        // getIllConditioning returns before computeIllConditioning clears
        assert_eq!(ill.record.len(), 1);
    }
}
