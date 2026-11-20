//! HiGHS source: lp_data/HighsInterface.cpp (Highs::elasticityFilterReturn,
//! Highs::elasticityFilter), lp_data/Highs.cpp (Highs::feasibilityRelaxation)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `elasticityFilter` is a single 613-line function and is kept
//! whole. Its local `solveLp` lambda (lines 2543-2554) becomes a nested
//! closure or a private helper taking `&mut self` and `&mut HighsIis`; it
//! must keep the C++ order of the `iis.info_.update` and `output_flag`
//! restore. `writeModel` / `writeSolution` calls (under `write_model` /
//! `kIisDevReport`, both false) are file-writer code: not ported.
//! `passModelName` / `passColName` / `passRowName` are ported.

use crate::prelude::*;
use crate::util::c_printf::c_fmt_g_width;

impl Highs {
    /// HiGHS: HighsInterface.cpp::Highs::elasticityFilterReturn (lines 2155-2225)
    /// HiGHS 1.16 (#3370): HighsInterface.cpp::elasticityFilterReturn (lines 2186-2259)
    pub(crate) fn elasticity_filter_return(
        &mut self,
        return_status: HighsStatus,
        original_model_name: &str,
        original_model_status: HighsModelStatus,
        original_num_col: HighsInt,
        original_num_row: HighsInt,
        original_sense: ObjSense,
        original_col_cost: &[f64],
        original_col_lower: &[f64],
        original_col_upper: &[f64],
        original_integrality: &[HighsVarType],
    ) -> HighsStatus {
        // C++ `const HighsLp& lp = model_.lp_;` is read through self.model_lp_ below
        // The IIS is cleared by restoring the original LP, so save it
        let iis = self.iis_.clone();
        let objective_function_value = self.info_.objective_function_value;
        // Delete any additional rows and columns, and restore the original
        // column costs and bounds
        let mut run_status;
        run_status = self.delete_rows_interval(original_num_row, self.model_lp_.num_row_ - 1);
        debug_assert!(run_status == HighsStatus::kOk);

        run_status = self.delete_cols_interval(original_num_col, self.model_lp_.num_col_ - 1);
        debug_assert!(run_status == HighsStatus::kOk);

        // Now that deleteRows and deleteCols may yield a valid basis, the
        // lack of dual values triggers an assert in
        // getKktFailures. Ultimately (#2081) the dual values will be
        // available but, for now, make the basis invalid.
        self.basis_.valid = false;

        run_status = self.change_cols_cost_interval(0, original_num_col - 1, original_col_cost);
        debug_assert!(run_status == HighsStatus::kOk);

        // HiGHS 1.16 #3370 (h16_feasrelax_sense): HighsInterface.cpp lines 2217-2218 (1.16.0)
        let h16_feasrelax_sense = self.options_.h16_feasrelax_sense;
        if h16_feasrelax_sense {
            run_status = self.change_objective_sense(original_sense);
            debug_assert!(run_status == HighsStatus::kOk);
        }

        run_status = self.change_cols_bounds_interval(
            0,
            original_num_col - 1,
            original_col_lower,
            original_col_upper,
        );
        debug_assert!(run_status == HighsStatus::kOk);

        if kIisStrategyRelaxation & self.options_.iis_strategy != 0
            && !original_integrality.is_empty()
        {
            self.change_cols_integrality_interval(0, original_num_col - 1, original_integrality);
            debug_assert!(run_status == HighsStatus::kOk);
        }
        // Revert the model name
        self.pass_model_name(original_model_name);
        debug_assert!(self.model_lp_.num_col_ == original_num_col);
        debug_assert!(self.model_lp_.num_row_ == original_num_row);

        // Clear solver data
        self.invalidate_solver_data();

        if return_status == HighsStatus::kOk {
            // Solution is invalidated by deleting rows and columns, but
            // primal values are correct. Have to recompute row activities,
            // though
            self.model_lp_.a_matrix_.product_quad(
                &mut self.solution_.row_value,
                &self.solution_.col_value,
                kDebugReportOff,
            );
            self.solution_.value_valid = true;
            // Set the feasibility objective and any KKT failures
            self.info_.objective_function_value = objective_function_value;
            crate::model::highs_model::get_kkt_failures_model(
                &self.options_,
                &self.model_lp_,
                &self.model_hessian_,
                &self.solution_,
                &self.basis_,
                &mut self.info_,
            );
            self.info_.valid = true;
        }

        // Revert model status
        self.model_status_ = original_model_status;

        // Restore IIS
        self.iis_ = iis;

        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::elasticityFilter (lines 2227-2839)
    /// HiGHS 1.16 (#3370): HighsInterface.cpp::elasticityFilter (lines 2261-2877)
    ///
    /// C++ default arguments: get_iis = false
    ///
    /// PORT NOTE: nullable penalty arrays are `Option<&[f64]>`.
    // evar_ix++ at the end of the last block is dead in the C++ too (kept for fidelity)
    #[allow(unused_assignments)]
    pub(crate) fn elasticity_filter(
        &mut self,
        global_lower_penalty: f64,
        global_upper_penalty: f64,
        global_rhs_penalty: f64,
        local_lower_penalty: Option<&[f64]>,
        local_upper_penalty: Option<&[f64]>,
        local_rhs_penalty: Option<&[f64]>,
        get_iis: bool,
    ) -> HighsStatus {
        // Line 2234: commented-out writeModel("infeasible.mps") dropped.
        //
        // Solve the feasibility relaxation problem for the given penalties,
        // continuing to act as the elasticity filter if get_iis
        // is true, resulting in an infeasibility subset for further
        // refinement as an IIS
        //
        // Construct the e-LP:
        //
        // Constraints L <= Ax <= U; l <= x <= u
        //
        // Transformed to
        //
        // L <= Ax + e_L - e_U <= U,
        //
        // l <=  x + e_l - e_u <= u,
        //
        // where the elastic variables are not used if the corresponding
        // bound is infinite or the local/global penalty is negative.
        //
        // x is free, and the objective is the linear function of the
        // elastic variables given by the local/global penalties
        //
        // col_of_ecol lists the column indices corresponding to the entries in
        // bound_of_col_of_ecol_is_lower so that the results can be interpreted
        //
        // row_of_ecol lists the row indices corresponding to the entries in
        // bound_of_row_of_ecol_is_lower so that the results can be interpreted
        let mut col_of_ecol: Vec<HighsInt> = Vec::new();
        let mut row_of_ecol: Vec<HighsInt> = Vec::new();
        let mut bound_of_row_of_ecol_is_lower: Vec<bool> = Vec::new();
        let mut bound_of_col_of_ecol_is_lower: Vec<bool> = Vec::new();
        let mut erow_lower: Vec<f64> = Vec::new();
        let mut erow_upper: Vec<f64> = Vec::new();
        let mut erow_start: Vec<HighsInt> = Vec::new();
        let mut erow_index: Vec<HighsInt> = Vec::new();
        let mut erow_value: Vec<f64> = Vec::new();
        // Accumulate names for ecols and erows, re-using ecol_name for the
        // names of row ecols after defining the names of col ecols
        let mut ecol_name: Vec<String> = Vec::new();
        let mut erow_name: Vec<String> = Vec::new();
        let mut ecol_cost: Vec<f64> = Vec::new();
        let mut ecol_lower: Vec<f64>;
        let mut ecol_upper: Vec<f64>;
        let lp = &self.model_lp_;
        let mut evar_ix: HighsInt = lp.num_col_;
        let mut run_status: HighsStatus;
        let write_model = false;
        // Take copies of the original model name, dimensions and column
        // data vectors, as they will be modified in forming the e-LP
        let original_model_name = lp.model_name_.clone();
        let original_model_status = self.get_model_status();
        let original_sense = lp.sense_;
        let original_num_col = lp.num_col_;
        let original_num_row = lp.num_row_;
        let original_col_cost = lp.col_cost_.clone();
        let original_col_lower = lp.col_lower_.clone();
        let original_col_upper = lp.col_upper_.clone();
        let original_integrality = lp.integrality_.clone();
        let mut iis = self.iis_.clone();
        // Clear solver data
        self.invalidate_solver_data();
        // Give the model a new name to avoid confusing logging when the
        // elastic LP is solved
        self.pass_model_name(&(original_model_name.clone() + "_elastic"));
        // Zero the column costs
        let zero_costs: Vec<f64> = vec![0.0; original_num_col as usize];
        run_status = self.change_cols_cost_interval(0, self.model_lp_.num_col_ - 1, &zero_costs);
        debug_assert!(run_status == HighsStatus::kOk);
        // HiGHS 1.16 #3370 (h16_feasrelax_sense): HighsInterface.cpp lines 2338-2340 (1.16.0)
        let h16_feasrelax_sense = self.options_.h16_feasrelax_sense;
        if h16_feasrelax_sense {
            run_status = self.change_objective_sense(ObjSense::kMinimize);
            debug_assert!(run_status == HighsStatus::kOk);
        }

        let mip_relaxation_iis = get_iis
            && self.model_lp_.is_mip()
            && (kIisStrategyRelaxation & self.options_.iis_strategy != 0);
        if mip_relaxation_iis {
            // Set any integrality to continuous
            let all_continuous: Vec<HighsVarType> =
                vec![HighsVarType::kContinuous; original_num_col as usize];
            run_status = self.change_cols_integrality_interval(
                0,
                self.model_lp_.num_col_ - 1,
                &all_continuous,
            );
            debug_assert!(run_status == HighsStatus::kOk);
        }

        // For the columns
        let has_local_lower_penalty = local_lower_penalty.is_some();
        let has_global_elastic_lower = global_lower_penalty >= 0.0;
        let has_elastic_lower = has_local_lower_penalty || has_global_elastic_lower;
        let has_local_upper_penalty = local_upper_penalty.is_some();
        let has_global_elastic_upper = global_upper_penalty >= 0.0;
        let has_elastic_upper = has_local_upper_penalty || has_global_elastic_upper;
        let has_elastic_columns = has_elastic_lower || has_elastic_upper;
        // For the rows
        let has_local_rhs_penalty = local_rhs_penalty.is_some();
        let has_global_elastic_rhs = global_rhs_penalty >= 0.0;
        let has_elastic_rows = has_local_rhs_penalty || has_global_elastic_rhs;
        debug_assert!(has_elastic_columns || has_elastic_rows);

        let col_ecol_offset = self.model_lp_.num_col_;
        if has_elastic_columns {
            // Accumulate bounds to be used for columns
            let mut col_lower: Vec<f64> = Vec::new();
            let mut col_upper: Vec<f64> = Vec::new();
            // When defining names, need to know the column number
            let has_col_names = self.model_lp_.col_names_.len() > 0;
            erow_start.push(0);
            let lp = &self.model_lp_;
            let mut i_col: HighsInt = 0;
            while i_col < lp.num_col_ {
                let lower = lp.col_lower_[i_col as usize];
                let upper = lp.col_upper_[i_col as usize];
                // Original bounds used unless e-variable introduced
                col_lower.push(lower);
                col_upper.push(upper);
                // Free columns have no erow
                if lower <= -kHighsInf && upper >= kHighsInf {
                    i_col += 1;
                    continue;
                }

                // Get the penalty for violating the lower bounds on this column
                let lower_penalty = if has_local_lower_penalty {
                    local_lower_penalty.unwrap()[i_col as usize]
                } else {
                    global_lower_penalty
                };
                // Negative lower penalty and infinite upper bound implies that the
                // bounds cannot be violated
                if lower_penalty < 0.0 && upper >= kHighsInf {
                    i_col += 1;
                    continue;
                }

                // Get the penalty for violating the upper bounds on this column
                let upper_penalty = if has_local_upper_penalty {
                    local_upper_penalty.unwrap()[i_col as usize]
                } else {
                    global_upper_penalty
                };
                // Infinite upper bound and negative lower penalty implies that the
                // bounds cannot be violated
                if lower <= -kHighsInf && upper_penalty < 0.0 {
                    i_col += 1;
                    continue;
                }
                // rusolve fix (not in HiGHS 1.16 HighsInterface.cpp): with both
                // bounds finite and both penalties negative no e-variable is made, so
                // the erow would hold x[iCol] alone and the assertion below
                // (row_nz == 2 || 3) would fire in a debug build. HiGHS only avoids
                // that when one bound is infinite (the two `continue`s above). The
                // column keeps both bounds, so the skipped row is redundant.
                if !((lower > -kHighsInf && lower_penalty >= 0.0)
                    || (upper < kHighsInf && upper_penalty >= 0.0))
                {
                    i_col += 1;
                    continue;
                }
                erow_lower.push(lower);
                erow_upper.push(upper);
                if has_col_names {
                    erow_name.push(format!(
                        "row_{}_{}_erow",
                        i_col, lp.col_names_[i_col as usize]
                    ));
                }
                // Define the entry for x[iCol]
                erow_index.push(i_col);
                erow_value.push(1.0);
                if lower > -kHighsInf && lower_penalty >= 0.0 {
                    // New e_l variable
                    col_of_ecol.push(i_col);
                    if has_col_names {
                        ecol_name.push(format!(
                            "col_{}_{}_lower",
                            i_col, lp.col_names_[i_col as usize]
                        ));
                    }
                    // Save the original lower bound on this column and free its
                    // lower bound
                    bound_of_col_of_ecol_is_lower.push(true);
                    col_lower[i_col as usize] = -kHighsInf;
                    erow_index.push(evar_ix);
                    erow_value.push(1.0);
                    ecol_cost.push(lower_penalty);
                    evar_ix += 1;
                }
                if upper < kHighsInf && upper_penalty >= 0.0 {
                    // New e_u variable
                    col_of_ecol.push(i_col);
                    if has_col_names {
                        ecol_name.push(format!(
                            "col_{}_{}_upper",
                            i_col, lp.col_names_[i_col as usize]
                        ));
                    }
                    // Save the original upper bound on this column and free its
                    // upper bound
                    bound_of_col_of_ecol_is_lower.push(false);
                    col_upper[i_col as usize] = kHighsInf;
                    erow_index.push(evar_ix);
                    erow_value.push(-1.0);
                    ecol_cost.push(upper_penalty);
                    evar_ix += 1;
                }
                erow_start.push(erow_index.len() as HighsInt);
                let row_nz = erow_start[erow_start.len() - 1] - erow_start[erow_start.len() - 2];
                // Lines 2405-2406: commented-out printf dropped.
                debug_assert!(row_nz == 2 || row_nz == 3);
                i_col += 1;
            }
            let num_new_col = col_of_ecol.len() as HighsInt;
            let num_new_row = (erow_start.len() as HighsInt) - 1;
            let num_new_nz = erow_start[num_new_row as usize];
            if kIisDevReport {
                println!(
                    "Elasticity filter: For columns there are {} variables and {} constraints",
                    num_new_col, num_new_row
                );
            }
            // Apply the original column bound changes
            debug_assert!(col_lower.len() == self.model_lp_.num_col_ as usize);
            debug_assert!(col_upper.len() == self.model_lp_.num_col_ as usize);
            run_status = self.change_cols_bounds_interval(
                0,
                self.model_lp_.num_col_ - 1,
                &col_lower,
                &col_upper,
            );
            debug_assert!(run_status == HighsStatus::kOk);
            // Add the new columns
            ecol_lower = vec![0.0; num_new_col as usize];
            ecol_upper = vec![kHighsInf; num_new_col as usize];
            run_status = self.add_cols(
                num_new_col,
                &ecol_cost,
                &ecol_lower,
                &ecol_upper,
                0,
                &[],
                &[],
                &[],
            );
            debug_assert!(run_status == HighsStatus::kOk);
            // Add the new rows
            debug_assert!(erow_start.len() == (num_new_row + 1) as usize);
            debug_assert!(erow_index.len() == num_new_nz as usize);
            debug_assert!(erow_value.len() == num_new_nz as usize);
            run_status = self.add_rows(
                num_new_row,
                &erow_lower,
                &erow_upper,
                num_new_nz,
                &erow_start,
                &erow_index,
                &erow_value,
            );
            debug_assert!(run_status == HighsStatus::kOk);
            if has_col_names {
                let mut i_col: HighsInt = 0;
                while i_col < num_new_col {
                    self.pass_col_name(col_ecol_offset + i_col, &ecol_name[i_col as usize]);
                    i_col += 1;
                }
                let mut i_row: HighsInt = 0;
                while i_row < num_new_row {
                    self.pass_row_name(original_num_row + i_row, &erow_name[i_row as usize]);
                    i_row += 1;
                }
            }
            debug_assert!(ecol_cost.len() == num_new_col as usize);
            debug_assert!(ecol_lower.len() == num_new_col as usize);
            debug_assert!(ecol_upper.len() == num_new_col as usize);
            if write_model {
                // Lines 2446-2453: writeModel dropped (file writer).
            }
        }
        let row_ecol_offset = self.model_lp_.num_col_;
        if has_elastic_rows {
            // Add the columns corresponding to the e_L and e_U variables for
            // the constraints
            ecol_name.clear();
            ecol_cost.clear();
            let mut ecol_start: Vec<HighsInt> = Vec::new();
            let mut ecol_index: Vec<HighsInt> = Vec::new();
            let mut ecol_value: Vec<f64> = Vec::new();
            ecol_start.push(0);
            let has_row_names = self.model_lp_.row_names_.len() > 0;
            let lp = &self.model_lp_;
            let mut i_row: HighsInt = 0;
            while i_row < original_num_row {
                // Get the penalty for violating the bounds on this row
                let penalty = if has_local_rhs_penalty {
                    local_rhs_penalty.unwrap()[i_row as usize]
                } else {
                    global_rhs_penalty
                };
                // Negative penalty implies that the bounds cannot be violated
                if penalty < 0.0 {
                    i_row += 1;
                    continue;
                }
                let lower = lp.row_lower_[i_row as usize];
                let upper = lp.row_upper_[i_row as usize];
                if lower > -kHighsInf {
                    // Create an e-var for the row lower bound
                    row_of_ecol.push(i_row);
                    if has_row_names {
                        ecol_name.push(format!(
                            "row_{}_{}_lower",
                            i_row, lp.row_names_[i_row as usize]
                        ));
                    }
                    bound_of_row_of_ecol_is_lower.push(true);
                    // Define the sub-matrix column
                    ecol_index.push(i_row);
                    ecol_value.push(1.0);
                    ecol_start.push(ecol_index.len() as HighsInt);
                    ecol_cost.push(penalty);
                    evar_ix += 1;
                }
                if upper < kHighsInf {
                    // Create an e-var for the row upper bound
                    row_of_ecol.push(i_row);
                    if has_row_names {
                        ecol_name.push(format!(
                            "row_{}_{}_upper",
                            i_row, lp.row_names_[i_row as usize]
                        ));
                    }
                    bound_of_row_of_ecol_is_lower.push(false);
                    // Define the sub-matrix column
                    ecol_index.push(i_row);
                    ecol_value.push(-1.0);
                    ecol_start.push(ecol_index.len() as HighsInt);
                    ecol_cost.push(penalty);
                    evar_ix += 1;
                }
                i_row += 1;
            }
            let num_new_col = (ecol_start.len() as HighsInt) - 1;
            let num_new_nz = ecol_start[num_new_col as usize];
            ecol_lower = vec![0.0; num_new_col as usize];
            ecol_upper = vec![kHighsInf; num_new_col as usize];
            debug_assert!(ecol_cost.len() == num_new_col as usize);
            debug_assert!(ecol_lower.len() == num_new_col as usize);
            debug_assert!(ecol_upper.len() == num_new_col as usize);
            debug_assert!(ecol_start.len() == (num_new_col + 1) as usize);
            debug_assert!(ecol_index.len() == num_new_nz as usize);
            debug_assert!(ecol_value.len() == num_new_nz as usize);
            run_status = self.add_cols(
                num_new_col,
                &ecol_cost,
                &ecol_lower,
                &ecol_upper,
                num_new_nz,
                &ecol_start,
                &ecol_index,
                &ecol_value,
            );
            debug_assert!(run_status == HighsStatus::kOk);
            if has_row_names {
                let mut i_col: HighsInt = 0;
                while i_col < num_new_col {
                    self.pass_col_name(row_ecol_offset + i_col, &ecol_name[i_col as usize]);
                    i_col += 1;
                }
            }

            if write_model {
                // Lines 2522-2529: writeModel dropped (file writer).
            }
        }

        if write_model {
            // Line 2532: writeModel("elastic.mps") dropped (file writer).
        }
        // Initial logging
        if get_iis {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Running elasticity filter to identify an infeasible subset of rows\n"
            );
            iis.report_iteration(&self.options_, 0, 0, true);
        }
        // Working with a copy of the IIS so clear this->iis_
        self.iis_.clear();
        // Lambda for gathering data when solving an LP
        // Lines 2543-2554: solveLp lambda becomes a private helper.

        // Solve the elastic LP
        run_status = self.solve_lp_for_elasticity(&mut iis);
        if run_status != HighsStatus::kOk {
            // Failed to establish feasibility or infeasibility due to unknown error
            // or hitting a time limit
            if self.model_status_ == HighsModelStatus::kTimeLimit {
                iis.status_ = kIisModelStatusTimeLimit;
                if get_iis {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Elasticity filter failed because time limit was reached\n"
                    );
                }
            } else {
                iis.status_ = kIisModelStatusUnknown;
                if get_iis {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Elasticity filter failed because it encountered an unknown model status\n"
                    );
                }
            }
            iis.valid_ = false;
            self.iis_ = iis;
            return self.elasticity_filter_return(
                HighsStatus::kError,
                &original_model_name,
                original_model_status,
                original_num_col,
                original_num_row,
                original_sense,
                &original_col_cost,
                &original_col_lower,
                &original_col_upper,
                &original_integrality,
            );
        }
        if kIisDevReport {
            // Line 2581: writeSolution dropped (file writer).
        }
        // Model status should be optimal, unless model is unbounded
        //
        // rusolve: HiGHS 1.16 (HighsInterface.cpp) has the narrower
        //   assert(this->model_status_ == HighsModelStatus::kOptimal ||
        //          this->model_status_ == HighsModelStatus::kUnbounded);
        // and it only exists in debug builds. A release HiGHS returns kOk with
        // model status kInfeasible when the elastic LP is infeasible. A feasibility
        // relaxation (not get_iis) with negative penalties can legitimately have an
        // infeasible elastic LP, when the bounds that may not be relaxed conflict,
        // so kInfeasible is accepted there. Release behaviour is unchanged.
        debug_assert!(
            self.model_status_ == HighsModelStatus::kOptimal
                || self.model_status_ == HighsModelStatus::kUnbounded
                || (!get_iis && self.model_status_ == HighsModelStatus::kInfeasible)
        );
        iis.valid_ = true;
        iis.status_ = if self.info_.objective_function_value > 0.0 {
            kIisModelStatusReducible
        } else {
            kIisModelStatusFeasible
        };
        if !get_iis {
            self.iis_ = iis;
            return self.elasticity_filter_return(
                HighsStatus::kOk,
                &original_model_name,
                original_model_status,
                original_num_col,
                original_num_row,
                original_sense,
                &original_col_cost,
                &original_col_lower,
                &original_col_upper,
                &original_integrality,
            );
        }
        // If getting an (I)IS, assume no elastic columns, so no additional rows
        debug_assert!(!has_elastic_columns);
        debug_assert!(has_elastic_rows);
        debug_assert!(original_num_row == self.model_lp_.num_row_);

        // Get solution
        // Now fix e-variables that are positive and re-solve until e-LP is
        // infeasible
        let mut loop_k: HighsInt = 0;
        let mut row_set: std::collections::BTreeSet<HighsInt> = std::collections::BTreeSet::new();

        loop {
            loop_k += 1;
            if kIisDevReport {
                println!("\nElasticity filter pass {}\n==============", loop_k);
            }
            let mut num_fixed: HighsInt = 0;
            if has_elastic_columns {
                let mut e_col: usize = 0;
                while e_col < col_of_ecol.len() {
                    let i_col = col_of_ecol[e_col];
                    if self.solution_.col_value[(col_ecol_offset + e_col as HighsInt) as usize]
                        > self.options_.primal_feasibility_tolerance
                    {
                        if kIisDevReport {
                            println!(
                                "E-col {:2} (column {:2}) corresponds to column {:2} with {} bound {} and has solution value {}",
                                e_col,
                                col_ecol_offset + e_col as HighsInt,
                                i_col,
                                if bound_of_col_of_ecol_is_lower[e_col] {
                                    "lower"
                                } else {
                                    "upper"
                                },
                                c_fmt_g_width(
                                    if bound_of_col_of_ecol_is_lower[e_col] {
                                        self.model_lp_.col_lower_[i_col as usize]
                                    } else {
                                        self.model_lp_.col_upper_[i_col as usize]
                                    },
                                    4,
                                    11,
                                ),
                                c_fmt_g_width(
                                    self.solution_.col_value
                                        [(col_ecol_offset + e_col as HighsInt) as usize],
                                    4,
                                    11,
                                )
                            );
                        }
                        self.change_col_bounds(col_ecol_offset + e_col as HighsInt, 0.0, 0.0);
                        num_fixed += 1;
                    }
                    e_col += 1;
                }
            }
            if has_elastic_rows {
                let mut e_col: usize = 0;
                while e_col < row_of_ecol.len() {
                    let i_row = row_of_ecol[e_col];
                    if self.solution_.col_value[(row_ecol_offset + e_col as HighsInt) as usize]
                        > self.options_.primal_feasibility_tolerance
                    {
                        if kIisDevReport {
                            println!(
                                "E-row {:2} (column {:2}) corresponds to    row {:2} with {} bound {} and has solution value {}",
                                e_col,
                                row_ecol_offset + e_col as HighsInt,
                                i_row,
                                if bound_of_row_of_ecol_is_lower[e_col] {
                                    "lower"
                                } else {
                                    "upper"
                                },
                                c_fmt_g_width(
                                    if bound_of_row_of_ecol_is_lower[e_col] {
                                        self.model_lp_.row_lower_[i_row as usize]
                                    } else {
                                        self.model_lp_.row_upper_[i_row as usize]
                                    },
                                    4,
                                    11,
                                ),
                                c_fmt_g_width(
                                    self.solution_.col_value
                                        [(row_ecol_offset + e_col as HighsInt) as usize],
                                    4,
                                    11,
                                )
                            );
                        }
                        self.change_col_bounds(row_ecol_offset + e_col as HighsInt, 0.0, 0.0);
                        num_fixed += 1;
                        row_set.insert(i_row);
                    }
                    e_col += 1;
                }
            }

            if num_fixed == 0 {
                // No elastic variables were positive, so problem is feasible
                iis.status_ = kIisModelStatusFeasible;
                iis.report_iteration(&self.options_, loop_k, row_set.len() as HighsInt, true);
                break;
            }
            let run_status = self.solve_lp_for_elasticity(&mut iis);
            let model_status = self.get_model_status();
            if run_status != HighsStatus::kOk {
                // Solve failed
                if model_status == HighsModelStatus::kTimeLimit {
                    iis.status_ = kIisModelStatusTimeLimit;

                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Elasticity filter failed because time limit was reached\n"
                    );
                } else {
                    iis.status_ = kIisModelStatusUnknown;

                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Elasticity filter failed because it encountered an unknown model status\n"
                    );
                }
                iis.valid_ = false;
                self.iis_ = iis;
                return self.elasticity_filter_return(
                    HighsStatus::kError,
                    &original_model_name,
                    original_model_status,
                    original_num_col,
                    original_num_row,
                    original_sense,
                    &original_col_cost,
                    &original_col_lower,
                    &original_col_upper,
                    &original_integrality,
                );
            }
            if kIisDevReport {
                // Line 2687: writeSolution dropped (file writer).
            }
            let terminate = model_status == HighsModelStatus::kInfeasible;
            iis.report_iteration(&self.options_, loop_k, row_set.len() as HighsInt, terminate);
            if terminate {
                break;
            }
        }
        let mut num_enforced_col_ecol: HighsInt = 0;
        let mut num_enforced_row_ecol: HighsInt = 0;
        if has_elastic_columns {
            let mut e_col: usize = 0;
            while e_col < col_of_ecol.len() {
                let i_col = col_of_ecol[e_col];
                if self.model_lp_.col_upper_[(col_ecol_offset + e_col as HighsInt) as usize] == 0.0
                {
                    num_enforced_col_ecol += 1;
                    if kIisDevReport {
                        println!(
                            "Col e-col {:2} (column {:2}) corresponds to column {:2} with {} bound {} and is enforced",
                            e_col,
                            col_ecol_offset + e_col as HighsInt,
                            i_col,
                            if bound_of_col_of_ecol_is_lower[e_col] {
                                "lower"
                            } else {
                                "upper"
                            },
                            c_fmt_g_width(
                                if bound_of_col_of_ecol_is_lower[e_col] {
                                    self.model_lp_.col_lower_[i_col as usize]
                                } else {
                                    self.model_lp_.col_upper_[i_col as usize]
                                },
                                4,
                                11,
                            )
                        );
                    }
                }
                e_col += 1;
            }
        }
        if has_elastic_rows {
            let mut e_col: usize = 0;
            while e_col < row_of_ecol.len() {
                let i_row = row_of_ecol[e_col];
                if self.model_lp_.col_upper_[(row_ecol_offset + e_col as HighsInt) as usize] == 0.0
                {
                    num_enforced_row_ecol += 1;
                    if iis.row_index_.is_empty()
                        || iis.row_index_[iis.row_index_.len() - 1] != i_row
                    {
                        iis.row_index_.push(i_row);
                    }

                    // Lines 2720-2724: commented-out code dropped.
                    if kIisDevReport {
                        println!(
                            "Row e-col {:2} (column {:2}) corresponds to    row {:2} with {} bound {} and is enforced",
                            e_col,
                            row_ecol_offset + e_col as HighsInt,
                            i_row,
                            if bound_of_row_of_ecol_is_lower[e_col] {
                                "lower"
                            } else {
                                "upper"
                            },
                            c_fmt_g_width(
                                if bound_of_row_of_ecol_is_lower[e_col] {
                                    self.model_lp_.row_lower_[i_row as usize]
                                } else {
                                    self.model_lp_.row_upper_[i_row as usize]
                                },
                                4,
                                11,
                            )
                        );
                    }
                }
                e_col += 1;
            }
        }
        let num_iis_row = iis.row_index_.len() as HighsInt;
        if iis.status_ == kIisModelStatusFeasible {
            debug_assert!(num_enforced_col_ecol == 0 && num_enforced_row_ecol == 0);
            debug_assert!(num_iis_row == 0);
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Elasticity filter failed to reproduce infeasibility\n"
            );
        } else {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Elasticity filter after {} passes found an infeasible subset of {} rows\n",
                loop_k,
                row_set.len()
            );
        }

        iis.valid_ = true;
        iis.strategy_ = self.options_.iis_strategy;
        if iis.status_ == kIisModelStatusFeasible {
            self.iis_ = iis;
            return self.elasticity_filter_return(
                HighsStatus::kOk,
                &original_model_name,
                original_model_status,
                original_num_col,
                original_num_row,
                original_sense,
                &original_col_cost,
                &original_col_lower,
                &original_col_upper,
                &original_integrality,
            );
        }

        // Model is infeasible because there are (at least) a positive
        // number of rows in the infeasibility set. Hence the IIS status is
        // reducible
        debug_assert!(iis.status_ == kIisModelStatusReducible);
        debug_assert!(num_iis_row > 0);

        // Have to be able to map from original row indices to position in
        // iis.row_index_
        let mut in_row_index: Vec<HighsInt> = vec![-1; original_num_row as usize];
        let mut i_x: HighsInt = 0;
        while i_x < num_iis_row {
            in_row_index[iis.row_index_[i_x as usize] as usize] = i_x;
            i_x += 1;
        }

        // Determine the columns with nonzeros in the row subset
        let mut nonzero_in_row_index: Vec<bool> = vec![false; original_num_col as usize];
        if self.model_lp_.a_matrix_.is_colwise() {
            let lp = &self.model_lp_;
            let mut i_col: HighsInt = 0;
            while i_col < original_num_col {
                let mut i_el = lp.a_matrix_.start_[i_col as usize];
                while i_el < lp.a_matrix_.start_[(i_col + 1) as usize] {
                    if in_row_index[lp.a_matrix_.index_[i_el as usize] as usize] >= 0 {
                        nonzero_in_row_index[i_col as usize] = true;
                    }
                    i_el += 1;
                }
                i_col += 1;
            }
        } else {
            let lp = &self.model_lp_;
            let mut i_x: HighsInt = 0;
            while i_x < num_iis_row {
                let i_row = iis.row_index_[i_x as usize];
                let mut i_el = lp.a_matrix_.start_[i_row as usize];
                while i_el < lp.a_matrix_.start_[(i_row + 1) as usize] {
                    nonzero_in_row_index[lp.a_matrix_.index_[i_el as usize] as usize] = true;
                    i_el += 1;
                }
                i_x += 1;
            }
        }
        // Now form iis.col_index_ and iis.col_bound_
        debug_assert!(iis.col_index_.len() == 0);
        let lp = &self.model_lp_;
        let mut i_col: HighsInt = 0;
        while i_col < original_num_col {
            if nonzero_in_row_index[i_col as usize] {
                iis.col_index_.push(i_col);
                if lp.col_lower_[i_col as usize] > -kHighsInf {
                    if lp.col_upper_[i_col as usize] < kHighsInf {
                        // Boxed
                        iis.col_bound_.push(kIisBoundStatusBoxed);
                    } else {
                        // LB
                        iis.col_bound_.push(kIisBoundStatusLower);
                    }
                } else {
                    if lp.col_upper_[i_col as usize] < kHighsInf {
                        // UB
                        iis.col_bound_.push(kIisBoundStatusUpper);
                    } else {
                        // Free
                        iis.col_bound_.push(kIisBoundStatusFree);
                    }
                }
            }
            i_col += 1;
        }
        // Now, go through the elastic columns for the rows, identifying the bound
        // that is required
        iis.row_bound_ = vec![-1; num_iis_row as usize];
        let mut e_col: usize = 0;
        while e_col < row_of_ecol.len() {
            let i_row = row_of_ecol[e_col];
            if self.model_lp_.col_upper_[(row_ecol_offset + e_col as HighsInt) as usize] == 0.0 {
                let i_x = in_row_index[i_row as usize];
                if iis.row_bound_[i_x as usize] == -1 {
                    // Unassigned
                    iis.row_bound_[i_x as usize] = if bound_of_row_of_ecol_is_lower[e_col] {
                        kIisBoundStatusLower
                    } else {
                        kIisBoundStatusUpper
                    };
                } else {
                    // Mark row as boxed when encountered twice, indicating both elastic
                    // variables were enforced
                    iis.row_bound_[i_x as usize] = kIisBoundStatusBoxed;
                }
            }
            e_col += 1;
        }
        debug_assert!(iis.row_bound_.len() == iis.row_index_.len());
        let mut i_x: HighsInt = 0;
        while i_x < num_iis_row {
            debug_assert!(iis.row_bound_[i_x as usize] >= 0);
            i_x += 1;
        }
        self.iis_ = iis;
        self.elasticity_filter_return(
            HighsStatus::kOk,
            &original_model_name,
            original_model_status,
            original_num_col,
            original_num_row,
            original_sense,
            &original_col_cost,
            &original_col_lower,
            &original_col_upper,
            &original_integrality,
        )
    }

    /// HiGHS: HighsInterface.cpp::Highs::elasticityFilter (lines 2543-2554)
    ///
    /// The `solveLp` lambda from `elasticityFilter`.
    fn solve_lp_for_elasticity(&mut self, iis: &mut HighsIis) -> HighsStatus {
        let output_flag = self.options_.output_flag;
        self.options_.output_flag = false;
        let mut simplex_time = -self.get_run_time();
        let mut simplex_iterations = -self.info_.simplex_iteration_count;
        let run_status = self.optimize_model();
        simplex_time += self.get_run_time();
        simplex_iterations += self.info_.simplex_iteration_count;
        iis.info_.update(simplex_time, simplex_iterations);
        self.options_.output_flag = output_flag;
        run_status
    }

    /// HiGHS: Highs.cpp::Highs::feasibilityRelaxation (lines 2201-2211)
    ///
    /// C++ calls `elasticityFilter` with its default `get_iis = false`.
    pub fn feasibility_relaxation(
        &mut self,
        global_lower_penalty: f64,
        global_upper_penalty: f64,
        global_rhs_penalty: f64,
        local_lower_penalty: Option<&[f64]>,
        local_upper_penalty: Option<&[f64]>,
        local_rhs_penalty: Option<&[f64]>,
    ) -> HighsStatus {
        self.elasticity_filter(
            global_lower_penalty,
            global_upper_penalty,
            global_rhs_penalty,
            local_lower_penalty,
            local_upper_penalty,
            local_rhs_penalty,
            false,
        )
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// Column-wise LP builder: `cols[j] = (cost, lower, upper, [(row, value)])`.
    #[allow(clippy::type_complexity)]
    fn lp(cols: &[(f64, f64, f64, &[(HighsInt, f64)])], rows: &[(f64, f64)]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = cols.len() as HighsInt;
        lp.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = vec![0];
        for &(cost, lower, upper, entries) in cols {
            lp.col_cost_.push(cost);
            lp.col_lower_.push(lower);
            lp.col_upper_.push(upper);
            for &(r, v) in entries {
                lp.a_matrix_.index_.push(r);
                lp.a_matrix_.value_.push(v);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        for &(l, u) in rows {
            lp.row_lower_.push(l);
            lp.row_upper_.push(u);
        }
        lp
    }

    fn quiet(lp: HighsLp) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs
    }

    fn iis_lp() -> HighsLp {
        lp(
            &[
                (0.0, 0.0, kHighsInf, &[(0, 1.0), (1, 1.0)]),
                (0.0, 0.0, kHighsInf, &[(0, 1.0), (1, -1.0)]),
            ],
            &[(-kHighsInf, 1.0), (2.0, kHighsInf)],
        )
    }

    #[test]
    fn feasibility_relaxation_unit_penalties() {
        // Cheapest repair: x1 = -0.5 (violation 0.5), x0 = 1.5.
        let mut highs = quiet(iis_lp());
        assert_eq!(
            highs.feasibility_relaxation(1.0, 1.0, 1.0, None, None, None),
            HighsStatus::kOk
        );
        assert!((highs.get_info().objective_function_value - 0.5).abs() <= 1e-9);
        let x = &highs.get_solution().col_value;
        assert!((x[0] - 1.5).abs() <= 1e-9 && (x[1] + 0.5).abs() <= 1e-9);
        // the original LP is restored
        assert_eq!((highs.model_lp_.num_col_, highs.model_lp_.num_row_), (2, 2));
        assert_eq!(highs.model_lp_.col_cost_, vec![0.0, 0.0]);
    }

    /// HiGHS check/TestIis.cpp issue-3368 LP (PR #3370 test), row-wise matrix:
    /// row 0: -x + y = 0; row 1: y >= 2; x in [0, 1] with cost 0; y free with cost 1.
    fn issue_3368_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![0.0, 1.0];
        lp.col_lower_ = vec![0.0, -kHighsInf];
        lp.col_upper_ = vec![1.0, kHighsInf];
        lp.row_lower_ = vec![0.0, 2.0];
        lp.row_upper_ = vec![0.0, kHighsInf];
        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 3];
        lp.a_matrix_.index_ = vec![0, 1, 1];
        lp.a_matrix_.value_ = vec![-1.0, 1.0, 1.0];
        lp
    }

    /// PR #3370 (issue-3368): with h16_feasrelax_sense on, the elastic LP is
    /// solved as a minimisation for both LP senses and the sense is restored.
    #[test]
    fn h16_feasrelax_on_solves_as_minimize_for_both_senses() {
        let mut highs = quiet(issue_3368_lp());
        assert_eq!(
            highs.set_option_value_bool("h16_feasrelax_sense", true),
            HighsStatus::kOk
        );
        for sense in [ObjSense::kMinimize, ObjSense::kMaximize] {
            assert_eq!(highs.change_objective_sense(sense), HighsStatus::kOk);
            assert_eq!(
                highs.feasibility_relaxation(1.0, 1.0, 1.0, None, None, None),
                HighsStatus::kOk
            );
            assert!((highs.get_info().objective_function_value - 1.0).abs() <= 1e-9);
            let x = &highs.get_solution().col_value;
            assert!((x[0] - 2.0).abs() <= 1e-9 && (x[1] - 2.0).abs() <= 1e-9);
            // the LP's sense is restored after the elastic solve
            assert_eq!(highs.model_lp_.sense_, sense);
        }
    }
    /// PR #3370 off path: with h16_feasrelax_sense off (the default), the elastic
    /// LP keeps the 1.15.1 behaviour. For a maximisation LP the elastic costs
    /// are maximised as well, giving objective 2 (not 1) at x = y = 0.
    #[test]
    fn h16_feasrelax_off_keeps_1_15_1_sense_behaviour() {
        let mut highs = quiet(issue_3368_lp());
        assert_eq!(
            highs.change_objective_sense(ObjSense::kMinimize),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.feasibility_relaxation(1.0, 1.0, 1.0, None, None, None),
            HighsStatus::kOk
        );
        assert!((highs.get_info().objective_function_value - 1.0).abs() <= 1e-9);
        assert_eq!(highs.model_lp_.sense_, ObjSense::kMinimize);

        let mut highs = quiet(issue_3368_lp());
        assert_eq!(
            highs.change_objective_sense(ObjSense::kMaximize),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.feasibility_relaxation(1.0, 1.0, 1.0, None, None, None),
            HighsStatus::kOk
        );
        assert!((highs.get_info().objective_function_value - 2.0).abs() <= 1e-9);
        let x = &highs.get_solution().col_value;
        assert!(x[0].abs() <= 1e-9 && x[1].abs() <= 1e-9);
        assert_eq!(highs.model_lp_.sense_, ObjSense::kMaximize);
    }
}
