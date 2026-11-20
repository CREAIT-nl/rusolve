//! HiGHS source: lp_data/HighsIis.cpp (lines 593-979)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `compute` drives a nested `Highs` instance (passModel, changeColsCost,
//! run, getBasicVariablesArray, getBasisTransposeSolve, getPrimalPhase1Dual,
//! changeColBounds/changeRowBounds, getCols/getRows). Those API methods are
//! ported by jobs IA1/IA2. PORT NOTE: `setLogCallback` / user-callback
//! propagation (lines 620-636) is interface code: the port propagates
//! `callback` with `Highs::set_callback` only when it is active, matching
//! the C++ conditions, and skips `user_log_callback` (always unset here).

use super::*;
use crate::prelude::*;

impl HighsIis {
    /// HiGHS: HighsIis.cpp::HighsIis::compute (lines 593-922)
    /// C++ default arguments: basis = nullptr
    pub fn compute(
        &mut self,
        lp: &HighsLp,
        options: &HighsOptions,
        callback: &HighsCallback,
        basis: Option<&HighsBasis>,
    ) -> HighsStatus {
        let log_options = &options.log_options;
        let col_priority = kIisStrategyColPriority & options.iis_strategy != 0;
        let row_priority = !col_priority;
        // Initially all columns and rows are candidates for the IIS
        let mut num_rows: HighsInt = 0;
        for i_col in 0..lp.num_col_ {
            let col_status = self.determine_bound_status(
                lp.col_lower_[i_col as usize],
                lp.col_upper_[i_col as usize],
                false,
            );
            self.add_col(i_col, col_status);
        }
        for i_row in 0..lp.num_row_ {
            let row_status = self.determine_bound_status(
                lp.row_lower_[i_row as usize],
                lp.row_upper_[i_row as usize],
                true,
            );
            self.add_row(i_row, row_status);
            if row_status != kIisBoundStatusDropped {
                num_rows += 1;
            }
        }

        let mut highs = Highs::default();
        highs.pass_options(options);
        highs.set_option_value_bool("output_flag", kIisDevReport);
        highs.set_option_value_str("presolve", kHighsOffString);
        highs.set_option_value_double(
            "time_limit",
            cpp_max(options.iis_time_limit - self.info_.sum_simplex_times, 0.0),
        );
        // Handle the callback propagation after setting output_flag false,
        // otherwise deprecation message for setLogCallback is echoed
        // PORT NOTE: user_log_callback skipped (always unset here).
        if callback.active[kCallbackLogging as usize]
            || callback.active[kCallbackSimplexInterrupt as usize]
        {
            // Setting the logging callbacks currently serves no purpose since
            // output_flag is set to kIisDevReport which is false (unless
            // developing) so that solves with this Highs instance are silent
            // C++ lines 627-630: commented-out setLogCallback block dropped.
            // PORT NOTE: C++ setCallback(callback.user_callback, data) also runs
            // with a null user_callback; Rust set_callback takes a non-Option
            // closure, so the call is made only when one is set.
            if let Some(ref cb) = callback.user_callback {
                highs.set_callback(cb.clone());
            }
            if callback.active[kCallbackLogging as usize] {
                highs.start_callback(kCallbackLogging);
            }
            if callback.active[kCallbackSimplexInterrupt as usize] {
                highs.start_callback(kCallbackSimplexInterrupt);
            }
        }
        // C++ lines 637-640: const references into highs. In Rust these are
        // read fresh at each use site to avoid holding an immutable borrow
        // across the mutable calls that follow.
        let mut run_status = highs.pass_model(lp.clone());
        debug_assert!(run_status == HighsStatus::kOk);
        if let Some(b) = basis {
            highs.set_basis(b, "");
        }
        // Zero the objective
        let mut cost = vec![0.0f64; lp.num_col_ as usize];
        run_status = highs.change_cols_cost_interval(0, lp.num_col_ - 1, &cost);
        debug_assert!(run_status == HighsStatus::kOk);
        // Solve the LP
        if let Some(b) = basis {
            highs.set_basis(b, "");
        }
        let use_sensitivity_filter = false;
        let mut _primal_phase1_dual: Vec<f64> = Vec::new();
        let mut row_deletion = false;
        let mut i_x: HighsInt = -1;
        let _drop_lower = false;

        // Lambda for gathering data when solving an LP
        // HiGHS: HighsIis.cpp::HighsIis::compute::solveLp (lines 657-758)
        // The C++ lambda captures by reference; here the mutable state is
        // threaded through explicit parameters.
        #[allow(clippy::assertions_on_constants)]
        let solve_lp = |highs: &mut Highs,
                        info_: &mut HighsIisInfo,
                        run_status: &mut HighsStatus,
                        _row_deletion: &bool,
                        _i_x: &HighsInt,
                        _drop_lower: &bool,
                        _use_sensitivity_filter: &bool,
                        _primal_phase1_dual: &mut Vec<f64>,
                        _options: &HighsOptions,
                        _lp: &HighsLp|
         -> HighsStatus {
            let mut simplex_time = -highs.get_run_time();
            let mut simplex_iterations = -highs.get_info().simplex_iteration_count;
            *run_status = highs.optimize_model();
            simplex_time += highs.get_run_time();
            simplex_iterations += highs.get_info().simplex_iteration_count;
            info_.update(simplex_time, simplex_iterations);
            if *run_status != HighsStatus::kOk {
                return *run_status;
            }
            let model_status = highs.get_model_status();
            if *_use_sensitivity_filter && model_status == HighsModelStatus::kInfeasible {
                // C++ lines 668-670: debug printf behind the compile-time
                // use_sensitivity_filter (false), not a runtime option — dropped.
                let output_flag = highs.get_option_value_bool("output_flag").unwrap_or(false);
                highs.set_option_value_bool("output_flag", true);
                let simplex_strategy = highs.get_option_value_int("simplex_strategy").unwrap_or(0);
                highs.set_option_value_int("simplex_strategy", kSimplexStrategyPrimal);
                // Solve the LP
                *run_status = highs.optimize_model();
                if *run_status != HighsStatus::kOk {
                    return *run_status;
                }
                // C++ line 680: highs.writeSolution("", kSolutionStylePretty) —
                // file writer, not ported.
                let basic_index = highs.get_basic_variables_array();
                let mut rhs = vec![0.0f64; _lp.num_row_ as usize];
                // Get duals for nonbasic rows, and initialise duals so that basic duals
                // are zero
                debug_assert!(101 == 202);

                {
                    let solution = highs.get_solution();
                    for i_row in 0.._lp.num_row_ {
                        let i_var = basic_index[i_row as usize];
                        let lower = if i_var < _lp.num_col_ {
                            _lp.col_lower_[i_var as usize]
                        } else {
                            _lp.row_lower_[(i_var - _lp.num_col_) as usize]
                        };
                        let upper = if i_var < _lp.num_col_ {
                            _lp.col_upper_[i_var as usize]
                        } else {
                            _lp.row_upper_[(i_var - _lp.num_col_) as usize]
                        };
                        let value = if i_var < _lp.num_col_ {
                            solution.col_value[i_var as usize]
                        } else {
                            solution.row_value[(i_var - _lp.num_col_) as usize]
                        };
                        if value < lower - _options.primal_feasibility_tolerance {
                            rhs[i_row as usize] = -1.0;
                        } else if value > upper + _options.primal_feasibility_tolerance {
                            rhs[i_row as usize] = 1.0;
                        }
                    }
                }
                let mut pi = HVector::default();
                pi.setup(_lp.num_row_);
                highs.get_basis_transpose_solve(&rhs, &mut pi.array, None, None);
                pi.count = _lp.num_row_;
                let mut reduced_costs_value: Vec<f64> = Vec::new();
                let mut reduced_costs_index: Vec<HighsInt> = Vec::new();
                _lp.a_matrix_.product_transpose_quad_hvector(
                    &mut reduced_costs_value,
                    &mut reduced_costs_index,
                    &pi,
                    kDebugReportOff,
                );

                *_primal_phase1_dual = highs.get_primal_phase1_dual().to_vec();
                let mut num_zero_dual: HighsInt = 0;
                {
                    let incumbent_basis = highs.get_basis();
                    let solution = highs.get_solution();
                    for i_col in 0.._lp.num_col_ {
                        let status = incumbent_basis.col_status[i_col as usize];
                        let dual = _primal_phase1_dual[i_col as usize];
                        let lower = _lp.col_lower_[i_col as usize];
                        let upper = _lp.col_upper_[i_col as usize];
                        let value = solution.col_value[i_col as usize];
                        if status != HighsBasisStatus::kBasic
                            && dual.abs() < _options.dual_feasibility_tolerance
                        {
                            num_zero_dual += 1;
                            // Small dual for nonbasic variable
                            // C++ lines 726-730: printf — dev report, dropped.
                            // C++ line 731: commented-out assert — dropped.
                        }
                    }
                    for i_row in 0.._lp.num_row_ {
                        let status = incumbent_basis.row_status[i_row as usize];
                        let dual = _primal_phase1_dual[(_lp.num_col_ + i_row) as usize];
                        let lower = _lp.row_lower_[i_row as usize];
                        let upper = _lp.row_upper_[i_row as usize];
                        let value = solution.row_value[i_row as usize];
                        if status != HighsBasisStatus::kBasic
                            && dual.abs() < _options.dual_feasibility_tolerance
                        {
                            num_zero_dual += 1;
                            // Small dual for nonbasic variable
                            // C++ lines 744-748: printf — dev report, dropped.
                            // C++ line 749: commented-out assert — dropped.
                        }
                    }
                }
                highs.set_option_value_bool("output_flag", output_flag);
                highs.set_option_value_int("simplex_strategy", simplex_strategy);
                debug_assert!(num_zero_dual == 0);
            }

            *run_status
        };

        solve_lp(
            &mut highs,
            &mut self.info_,
            &mut run_status,
            &row_deletion,
            &i_x,
            &_drop_lower,
            &use_sensitivity_filter,
            &mut _primal_phase1_dual,
            options,
            lp,
        );
        // If we fail to establish infeasibility, return the initial subset
        if run_status != HighsStatus::kOk {
            self.valid_ = true;
            self.strategy_ = options.iis_strategy;
            if highs.get_model_status() == HighsModelStatus::kTimeLimit {
                self.status_ = kIisModelStatusTimeLimit;
            } else {
                self.status_ = kIisModelStatusReducible;
            }
            return HighsStatus::kWarning;
        }

        debug_assert!(highs.get_model_status() == HighsModelStatus::kInfeasible);
        let mut iis_status: IisModelStatus = kIisModelStatusIrreducible;
        let mut search_return_status = HighsStatus::kOk;

        // Initial logging
        self.clear_log_info();
        let mut iter: HighsInt = 0;
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "\nRunning deletion filter to identify an IIS\n"
        );
        self.report_iteration(options, iter, num_rows, true);

        // Pass twice: rows before columns, or columns before rows, according to
        // row_priority
        for k in 0..2 {
            row_deletion = (row_priority && k == 0) || (!row_priority && k == 1);
            let type_str = if row_deletion { "Row" } else { "Col" };
            // Perform deletion pass
            let num_index: HighsInt = if row_deletion {
                lp.num_row_
            } else {
                lp.num_col_
            };
            i_x = 0;
            while i_x < num_index {
                // Get logging info
                iter += 1;
                let force = row_deletion && (i_x == (num_index - 1));
                // Get status
                let ix_status: HighsInt = if row_deletion {
                    self.row_bound_[i_x as usize]
                } else {
                    self.col_bound_[i_x as usize]
                };
                // Skip if status is already free or dropped
                if ix_status == kIisBoundStatusDropped || ix_status == kIisBoundStatusFree {
                    // Possibly report
                    self.report_iteration(options, iter, num_rows, force);
                    i_x += 1;
                    continue;
                }

                let mut lower = if row_deletion {
                    lp.row_lower_[i_x as usize]
                } else {
                    lp.col_lower_[i_x as usize]
                };
                let mut upper = if row_deletion {
                    lp.row_upper_[i_x as usize]
                } else {
                    lp.col_upper_[i_x as usize]
                };
                // Record whether the upper bound has been dropped due to the lower
                // bound being kept
                if lower > -kHighsInf && iis_status != kIisModelStatusTimeLimit {
                    // Drop the lower bound temporarily
                    // C++ `bool drop_lower = true;` shadows the outer variable;
                    // the solveLp lambda still captures the outer one (false).
                    let _shadow_drop_lower = true;
                    run_status = if row_deletion {
                        highs.change_row_bounds(i_x, -kHighsInf, upper)
                    } else {
                        highs.change_col_bounds(i_x, -kHighsInf, upper)
                    };
                    debug_assert!(run_status == HighsStatus::kOk);
                    // Solve the LP
                    solve_lp(
                        &mut highs,
                        &mut self.info_,
                        &mut run_status,
                        &row_deletion,
                        &i_x,
                        &_drop_lower,
                        &use_sensitivity_filter,
                        &mut _primal_phase1_dual,
                        options,
                        lp,
                    );
                    self.process_bound_relaxation(
                        &mut highs,
                        row_deletion,
                        true,
                        i_x,
                        &mut lower,
                        &mut upper,
                        &mut iis_status,
                        &mut search_return_status,
                    );
                }
                if upper < kHighsInf && iis_status != kIisModelStatusTimeLimit {
                    // Drop the upper bound temporarily
                    run_status = if row_deletion {
                        highs.change_row_bounds(i_x, lower, kHighsInf)
                    } else {
                        highs.change_col_bounds(i_x, lower, kHighsInf)
                    };
                    debug_assert!(run_status == HighsStatus::kOk);
                    // Solve the LP
                    solve_lp(
                        &mut highs,
                        &mut self.info_,
                        &mut run_status,
                        &row_deletion,
                        &i_x,
                        &_drop_lower,
                        &use_sensitivity_filter,
                        &mut _primal_phase1_dual,
                        options,
                        lp,
                    );
                    self.process_bound_relaxation(
                        &mut highs,
                        row_deletion,
                        false,
                        i_x,
                        &mut lower,
                        &mut upper,
                        &mut iis_status,
                        &mut search_return_status,
                    );
                }
                let debug_bound_change = true;
                if debug_bound_change {
                    // Check bounds have been changed correctly
                    let mut check_lower = 0.0f64;
                    let mut check_upper = 0.0f64;
                    let mut check_cost = 0.0f64;
                    let mut check_num_ix: HighsInt = 0;
                    let mut check_num_nz: HighsInt = 0;
                    run_status = if row_deletion {
                        highs.get_rows_interval(
                            i_x,
                            i_x,
                            &mut check_num_ix,
                            Some(std::slice::from_mut(&mut check_lower)),
                            Some(std::slice::from_mut(&mut check_upper)),
                            &mut check_num_nz,
                            None,
                            None,
                            None,
                        )
                    } else {
                        highs.get_cols_interval(
                            i_x,
                            i_x,
                            &mut check_num_ix,
                            Some(std::slice::from_mut(&mut check_cost)),
                            Some(std::slice::from_mut(&mut check_lower)),
                            Some(std::slice::from_mut(&mut check_upper)),
                            &mut check_num_nz,
                            None,
                            None,
                            None,
                        )
                    };
                    debug_assert!(run_status == HighsStatus::kOk);
                    debug_assert!(check_lower == lower);
                    debug_assert!(check_upper == upper);
                }
                let iss_bound_status = self.determine_bound_status(lower, upper, row_deletion);
                if row_deletion {
                    self.row_bound_[i_x as usize] = iss_bound_status;
                    if iss_bound_status == kIisBoundStatusDropped {
                        num_rows -= 1;
                    }
                } else {
                    self.col_bound_[i_x as usize] = iss_bound_status;
                }
                // Possibly report on iteration
                self.report_iteration(options, iter, num_rows, force);
                if kIisDevReport {
                    crate::highs_log_user!(
                        log_options,
                        HighsLogType::kInfo,
                        "{} {} has status {}\n",
                        type_str,
                        i_x,
                        self.iis_bound_status_to_string(iss_bound_status)
                    );
                }
                i_x += 1;
            }
            if k == 1 {
                continue;
            }
            // End of first pass: look to simplify second pass
            if kIisDevReport {
                self.report("End of deletion", highs.get_lp());
            }
            if row_deletion {
                // Mark empty columns as dropped
                for i_col in 0..lp.num_col_ {
                    let mut empty_col = true;
                    let mut i_el = lp.a_matrix_.start_[i_col as usize];
                    while i_el < lp.a_matrix_.start_[(i_col + 1) as usize] {
                        if self.row_bound_[lp.a_matrix_.index_[i_el as usize] as usize]
                            != kIisBoundStatusDropped
                        {
                            empty_col = false;
                            break;
                        }
                        i_el += 1;
                    }
                    if empty_col {
                        if kIisDevReport {
                            crate::highs_log_user!(
                                log_options,
                                HighsLogType::kInfo,
                                "Col {} has status Dropped: Empty\n",
                                i_col
                            );
                        }
                        self.col_bound_[i_col as usize] = kIisBoundStatusDropped;
                        run_status = highs.change_col_bounds(i_col, -kHighsInf, kHighsInf);
                        debug_assert!(run_status == HighsStatus::kOk);
                    }
                }
            }
            if kIisDevReport {
                self.report("End of pass 1", highs.get_lp());
            }
        }
        if kIisDevReport {
            self.report("End of pass 2", highs.get_lp());
        }
        let mut iss_num_col: HighsInt = 0;
        for i_col in 0..lp.num_col_ {
            if self.col_bound_[i_col as usize] != kIisBoundStatusDropped {
                self.col_index_[iss_num_col as usize] = self.col_index_[i_col as usize];
                self.col_bound_[iss_num_col as usize] = self.col_bound_[i_col as usize];
                iss_num_col += 1;
            }
        }
        let mut iss_num_row: HighsInt = 0;
        for i_row in 0..lp.num_row_ {
            if self.row_bound_[i_row as usize] != kIisBoundStatusDropped {
                self.row_index_[iss_num_row as usize] = self.row_index_[i_row as usize];
                self.row_bound_[iss_num_row as usize] = self.row_bound_[i_row as usize];
                iss_num_row += 1;
            }
        }
        // Return final result
        self.valid_ = true;
        self.status_ = iis_status;
        self.strategy_ = options.iis_strategy;
        self.col_index_.resize(iss_num_col as usize, 0);
        self.col_bound_.resize(iss_num_col as usize, 0);
        self.row_index_.resize(iss_num_row as usize, 0);
        self.row_bound_.resize(iss_num_row as usize, 0);
        search_return_status
    }

    /// HiGHS: HighsIis.cpp::HighsIis::processBoundRelaxation (lines 924-979)
    pub fn process_bound_relaxation(
        &mut self,
        highs: &mut Highs,
        row_deletion: bool,
        drop_lower: bool,
        i_x: HighsInt,
        lower: &mut f64,
        upper: &mut f64,
        iis_status: &mut IisModelStatus,
        search_return_status: &mut HighsStatus,
    ) {
        // Declare run_status
        let mut run_status;
        // Get model status
        let model_status = highs.get_model_status();
        if model_status == HighsModelStatus::kOptimal {
            // Now feasible, so restore the bound
            run_status = if row_deletion {
                highs.change_row_bounds(i_x, *lower, *upper)
            } else {
                highs.change_col_bounds(i_x, *lower, *upper)
            };
            debug_assert!(run_status == HighsStatus::kOk);
            // If the lower bound must be kept, then any finite upper bound
            // must be dropped
            let apply_reciprocal_rule = drop_lower && (*upper < kHighsInf);
            if apply_reciprocal_rule {
                // Drop the upper bound permanently
                *upper = kHighsInf;
                run_status = if row_deletion {
                    highs.change_row_bounds(i_x, *lower, *upper)
                } else {
                    highs.change_col_bounds(i_x, *lower, *upper)
                };
                debug_assert!(run_status == HighsStatus::kOk);
            }
            if drop_lower {
                debug_assert!(*upper >= kHighsInf);
            }
            // Since upper = kHighsInf, allow the loop to run so that
            // bound status is set as if upper were set to kHighsInf
            // by relaxing it and finding that the LP was still
            // infeasible
        } else if model_status == HighsModelStatus::kInfeasible {
            // Bound can be dropped permanently
            if drop_lower {
                *lower = -kHighsInf;
            } else {
                *upper = kHighsInf;
            }
        } else if model_status == HighsModelStatus::kTimeLimit {
            // Time limit reached, so restore the bound and set iis status to
            // kIisModelStatusTimeLimit
            run_status = if row_deletion {
                highs.change_row_bounds(i_x, *lower, *upper)
            } else {
                highs.change_col_bounds(i_x, *lower, *upper)
            };
            debug_assert!(run_status == HighsStatus::kOk);
            *iis_status = kIisModelStatusTimeLimit;
            *search_return_status = HighsStatus::kWarning;
        } else {
            // Unknown failure, so restore the bound and set iis status to
            // kIisModelStatusReducible
            run_status = if row_deletion {
                highs.change_row_bounds(i_x, *lower, *upper)
            } else {
                highs.change_col_bounds(i_x, *lower, *upper)
            };
            debug_assert!(run_status == HighsStatus::kOk);
            *iis_status = kIisModelStatusReducible;
            *search_return_status = HighsStatus::kWarning;
        }
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

    fn _quiet(lp: HighsLp) -> Highs {
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
    fn compute_irreducible_from_lp() {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.iis_strategy = kIisStrategyFromLp | kIisStrategyIrreducible;
        let callback = HighsCallback::default();
        let mut iis = HighsIis::default();
        assert_eq!(
            iis.compute(&iis_lp(), &options, &callback, None),
            HighsStatus::kOk
        );
        let mut rows = iis.row_index_.clone();
        rows.sort();
        assert_eq!(rows, vec![0, 1]);
        // C++ semantics: dropping column 0's lower bound leaves it with both
        // bounds infinite, which determineBoundStatus maps to kIisBoundStatusFree
        // (not Dropped), so column 0 stays in the index with status Free.
        assert_eq!(iis.col_index_, vec![0, 1]);
        assert_eq!(
            iis.col_bound_,
            vec![kIisBoundStatusFree, kIisBoundStatusLower]
        );
    }
}
