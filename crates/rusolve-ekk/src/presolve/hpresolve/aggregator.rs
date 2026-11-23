//! HiGHS source: presolve/HPresolve.cpp (lines 6476-7021)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Dependent rows/columns, the aggregator and removal helpers. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::removeDependentEquations (lines 6476-6595, 120 lines)
    ///
    /// Uses `HFactor` (`crate::util::hfactor`): `factor.setup(matrix, colSet)` ->
    /// `factor.setup(&matrix, &mut col_set, kDefaultPivotThreshold, kDefaultPivotTolerance, kHighsDebugLevelMin, None)`
    /// (the C++ default arguments); `factor.build()` takes `basic_index` explicitly (see `HFactor::build`).
    pub fn remove_dependent_equations(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::DEPENDENT_EQ);
        debug_assert!(self.analysis_.allow_rule_[kPresolveRuleDependentEquations as usize]);
        let logging_on = self.analysis_.logging_on_;
        if self.equations.is_empty() {
            return HPresolveResult::kOk;
        }

        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleDependentEquations);
        }

        let mut matrix = HighsSparseMatrix::default();
        matrix.num_col_ = self.equations.len() as HighsInt;
        matrix.num_row_ = self.model.num_col_ + 1;
        matrix.start_.resize(matrix.num_col_ as usize + 1, 0);
        matrix.start_[0] = 0;
        let max_capacity = self.num_nonzeros() + matrix.num_col_;
        matrix.value_.reserve(max_capacity as usize);
        matrix.index_.reserve(max_capacity as usize);

        let mut eq_set: Vec<HighsInt> = vec![0; matrix.num_col_ as usize];
        let mut i: usize = 0;
        for &(_, eq) in self.equations.iter() {
            eq_set[i] = eq;
            i += 1;

            // add entries of equation
            for nonz in self.get_row_vector(eq).iter(self) {
                matrix.value_.push(nonz.value());
                matrix.index_.push(nonz.index());
            }

            // add entry for artificial rhs column
            if self.model.row_lower_[eq as usize] != 0.0 {
                matrix.value_.push(self.model.row_lower_[eq as usize]);
                matrix.index_.push(self.model.num_col_);
            }

            matrix.start_[i] = matrix.value_.len() as HighsInt;
        }
        let mut col_set: Vec<HighsInt> = (0..matrix.num_col_).collect();
        let mut factor = HFactor::default();
        // HiGHS 1.16 #3354 (h16_hfactor_pivot): the C++ local `HFactor factor;` runs the 1.16.0
        // kernel
        factor.h16_hfactor_pivot = self.options.h16_hfactor_pivot;
        factor.setup(
            &matrix,
            &mut col_set,
            kDefaultPivotThreshold,
            kDefaultPivotTolerance,
            kHighsDebugLevelMin,
            None,
        );
        // Set up a time limit to prevent the redundant rows factorization
        // taking forever.
        //
        // Allow no more than 1% of the time limit to be spent on removing
        // dependent equations, but ensure that there is some limit since
        // options->time_limit is infinity by default
        //
        // ToDo: This is strictly non-deterministic, but so conservative
        // that it'll only reap the cases when factor.build never finishes
        // HiGHS 1.16 #3167/#3398 (h16_dep_eq_time_cap): HPresolveDependent.cpp lines 91-93 (1.16.0),
        // kMaxDependentEquationsTime = 100; 1.15.1 HPresolve.cpp:6528 used 1000
        // (ff138f533 introduced the constant at 10, settled at 100 in 31b27b7eb)
        // rusolve presolve_dep_eq_fixed_budget: the budget no longer
        // scales with options.time_limit, so presolve does not change with the limit
        let time_limit: f64 = dependent_equations_time_limit(
            if self.options.presolve_dep_eq_fixed_budget {
                f64::INFINITY
            } else {
                self.options.time_limit
            },
            self.options.h16_dep_eq_time_cap,
        );
        factor.set_time_limit(time_limit);
        let silent = self.silent_log();
        // Determine rank deficiency of the equations
        if !silent {
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Dependent equations search running on {} equations with time limit of {:.2}s\n",
                matrix.num_col_,
                time_limit
            );
        }
        let mut time_taken: f64 = -self.timer.as_ref().unwrap().read();
        let build_return: HighsInt = factor.build(&mut col_set, None);
        time_taken += self.timer.as_ref().unwrap().read();
        if build_return == kBuildKernelReturnTimeout {
            // HFactor::build has timed out, so just return
            if !silent {
                crate::highs_log_user!(
                    &self.options.log_options,
                    HighsLogType::kInfo,
                    "Dependent equations search terminated after {:.3}s due to expected time exceeding limit\n",
                    time_taken
                );
            }
            self.analysis_.logging_on_ = logging_on;
            if logging_on {
                // Sic: the C++ stops kPresolveRuleDependentFreeCols on this path.
                self.analysis_stop_presolve_rule_log(kPresolveRuleDependentFreeCols);
            }
            return HPresolveResult::kOk;
        } else {
            let pct_off_timeout = 1e2 * (time_taken - time_limit).abs() / time_limit;
            if !silent && pct_off_timeout < 1.0 {
                crate::highs_log_user!(
                    &self.options.log_options,
                    HighsLogType::kWarning,
                    "Dependent equations search finished within {:.2}% of limit of {:.2}s: risk of non-deterministic behaviour if solve is repeated\n",
                    pct_off_timeout,
                    time_limit
                );
            }
        }
        // build_return as rank_deficiency must be valid
        debug_assert!(build_return >= 0);
        let rank_deficiency: HighsInt = build_return;
        // Analyse what's been removed
        let mut num_removed_row: HighsInt = 0;
        let mut num_removed_nz: HighsInt = 0;
        let mut num_fictitious_rows_skipped: HighsInt = 0;
        for k in 0..rank_deficiency as usize {
            if factor.var_with_no_pivot[k] >= 0 {
                let redundant_row = eq_set[factor.var_with_no_pivot[k] as usize];
                num_removed_row += 1;
                num_removed_nz += self.rowsize[redundant_row as usize];
                postsolve_stack.redundant_row(redundant_row);
                self.remove_row(redundant_row);
            } else {
                num_fictitious_rows_skipped += 1;
            }
        }
        if !silent {
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Dependent equations search removed {} rows and {} nonzeros in {:.2}s (limit = {:.2}s)\n",
                num_removed_row,
                num_removed_nz,
                time_taken,
                time_limit
            );
        }
        if num_fictitious_rows_skipped != 0 {
            crate::highs_log_dev!(
                &self.options.log_options,
                HighsLogType::kInfo,
                ", avoiding {} fictitious rows",
                num_fictitious_rows_skipped
            );
        }
        crate::highs_log_dev!(&self.options.log_options, HighsLogType::kInfo, "\n");

        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleDependentEquations);
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::removeDependentFreeCols (lines 6597-6686, 90 lines)
    pub fn remove_dependent_free_cols(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        // The C++ body is commented out and returns kOk.
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::aggregator (lines 6688-6806, 119 lines)
    ///
    /// MUST-PORT (determinism): `pdqsort(..., comp)` (C++ 6703) -> `pdqsort_by`; the comparator hashes
    /// `HighsHashHelpers::hash(std::make_pair(uint32_t(a), uint32_t(b)))` -> `hash_u32_pair(a as u32, b as u32)`;
    /// `std::make_tuple` comparison -> Rust tuple `<`.
    /// `HighsHashTable<HighsInt, HighsInt> fillinCache` -> `HighsHashTable::new()`.
    pub fn aggregator(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::AGGREGATOR);
        debug_assert!(self.analysis_.allow_rule_[kPresolveRuleAggregator as usize]);
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleAggregator);
        }
        {
            let mut opportunities = std::mem::take(&mut self.substitution_opportunities);
            // `erase(remove_if(pred))` keeps the elements where `pred` is false.
            opportunities.retain(|p| {
                let row = p.0;
                let col = p.1;
                !(self.row_deleted[row as usize] != 0
                    || self.col_deleted[col as usize] != 0
                    || !self.is_implied_free(col)
                    || !self.is_dual_implied_free(row))
            });
            self.substitution_opportunities = opportunities;
        }

        {
            let mut opportunities = std::mem::take(&mut self.substitution_opportunities);
            let this: &HPresolve = &*self;
            pdqsort_by(
                &mut opportunities,
                |nz1: &(HighsInt, HighsInt), nz2: &(HighsInt, HighsInt)| {
                    let min_len1 =
                        cpp_min(this.rowsize[nz1.0 as usize], this.colsize[nz1.1 as usize]);
                    let min_len2 =
                        cpp_min(this.rowsize[nz2.0 as usize], this.colsize[nz2.1 as usize]);
                    if min_len1 == 2 && min_len2 != 2 {
                        return true;
                    }
                    if min_len2 == 2 && min_len1 != 2 {
                        return false;
                    }

                    let size_prod1: i64 = (this.rowsize[nz1.0 as usize] as i64)
                        * (this.colsize[nz1.1 as usize] as i64);
                    let size_prod2: i64 = (this.rowsize[nz2.0 as usize] as i64)
                        * (this.colsize[nz2.1 as usize] as i64);
                    if size_prod1 < size_prod2 {
                        return true;
                    }
                    if size_prod2 < size_prod1 {
                        return false;
                    }
                    if min_len1 < min_len2 {
                        return true;
                    }
                    if min_len2 < min_len1 {
                        return false;
                    }

                    (hash_u32_pair(nz1.0 as u32, nz1.1 as u32), nz1.0, nz1.1)
                        < (hash_u32_pair(nz2.0 as u32, nz2.1 as u32), nz2.0, nz2.1)
                },
            );
            self.substitution_opportunities = opportunities;
        }

        let mut nfail: HighsInt = 0;
        let mut i: usize = 0;
        // `i < size()` is re-evaluated each pass, as in the C++ loop.
        while i < self.substitution_opportunities.len() {
            let cur = i;
            i += 1;
            let row = self.substitution_opportunities[cur].0;
            let col = self.substitution_opportunities[cur].1;

            if self.row_deleted[row as usize] != 0
                || self.col_deleted[col as usize] != 0
                || !self.is_implied_free(col)
                || !self.is_dual_implied_free(row)
            {
                self.substitution_opportunities[cur].0 = -1;
                continue;
            }

            let nz_pos = self.find_nonzero(row, col);
            if nz_pos == -1 {
                self.substitution_opportunities[cur].0 = -1;
                continue;
            }
            if self.model.integrality_[col as usize] == HighsVarType::kInteger {
                let implied_integral = self.is_implied_integral(col);
                hpresolve_checked_call!(implied_integral.result());
                if !implied_integral.flag() {
                    continue;
                }
            }

            // in the case where the row has length two or the column has length two
            // we always do the substitution since the fillin can never be problematic
            if self.rowsize[row as usize] == 2 || self.colsize[col as usize] == 2 {
                self.store_row(row);
                self.substitute_free_col(postsolve_stack, row, col, true);
                self.substitution_opportunities[cur].0 = -1;
                hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));
                hpresolve_checked_call!(self.check_limits(postsolve_stack));
                continue;
            }

            let mut max_val: f64 = if self.rowsize[row as usize] < self.colsize[col as usize] {
                self.get_max_abs_row_val(row)
            } else {
                self.get_max_abs_col_val(col)
            };
            if self.a_value[nz_pos as usize].abs() < max_val * self.options.presolve_pivot_threshold
            {
                max_val = if self.rowsize[row as usize] < self.colsize[col as usize] {
                    self.get_max_abs_col_val(col)
                } else {
                    self.get_max_abs_row_val(row)
                };
                if self.a_value[nz_pos as usize].abs()
                    < max_val * self.options.presolve_pivot_threshold
                {
                    self.substitution_opportunities[cur].0 = -1;
                    continue;
                }
            }

            self.store_row(row);
            let mut fillin: HighsInt =
                -(self.rowsize[row as usize] + self.colsize[col as usize] - 1);
            let mut it = self.get_column_vector(col);
            while let Some(nz) = it.next_nz(&*self) {
                if nz.index() == row {
                    continue;
                }
                fillin += self.count_fillin(nz.index());

                if fillin > self.options.presolve_substitution_maxfillin {
                    break;
                }
            }

            if fillin > self.options.presolve_substitution_maxfillin {
                nfail += 1;
                // if the fill in is too much for multiple tries, then we stop
                // as this indicates that the rows/columns are becoming too dense
                // for substitutions
                if nfail == 3 {
                    break;
                }
                continue;
            }

            nfail = 0;
            self.substitute_free_col(postsolve_stack, row, col, true);
            self.substitution_opportunities[cur].0 = -1;
            hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));
            hpresolve_checked_call!(self.check_limits(postsolve_stack));
        }

        self.substitution_opportunities.retain(|p| p.0 != -1);

        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleAggregator);
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::substitute (lines 6808-6848, 41 lines)
    /// HiGHS 1.16 (#3268): HPresolve.cpp::substitute (lines 4716-4725)
    ///
    /// C++ overload `substitute(substcol, staycol, offset, scale)`.
    pub fn substitute_col(
        &mut self,
        substcol: HighsInt,
        staycol: HighsInt,
        offset: f64,
        scale: f64,
    ) {
        let h16_subst_int = self.options.h16_subst_int;
        // HiGHS 1.16 #3268 (h16_subst_int): HPresolve.cpp lines 4716-4725 (1.16.0)
        // Preserve explicit integrality, i.e., upgrade implied integral
        // column if it is substituting an integral column
        if h16_subst_int
            && self.model.integrality_[substcol as usize] == HighsVarType::kInteger
            && self.model.integrality_[staycol as usize] == HighsVarType::kImplicitInteger
        {
            self.model.integrality_[staycol as usize] = HighsVarType::kInteger;
            let mut it = self.get_column_vector(staycol);
            while let Some(nonzero) = it.next_nz(self) {
                self.rowsize_integer[nonzero.index() as usize] += 1;
                self.rowsize_impl_int[nonzero.index() as usize] -= 1;
            }
        }

        let h16_row_duals = self.options.h16_row_duals;
        // substitute the column in each row where it occurs
        let mut coliter = self.colhead[substcol as usize];
        while coliter != -1 {
            let colrow = self.a_row[coliter as usize];
            let colval = self.a_value[coliter as usize];
            // walk to the next position before doing any modifications, because
            // the current position will be deleted in the loop below
            debug_assert!(self.a_col[coliter as usize] == substcol);
            let colpos = coliter;
            coliter = self.a_next[coliter as usize];
            debug_assert!(self.row_deleted[colrow as usize] == 0);
            self.unlink(colpos);

            // adjust the sides
            if h16_row_duals {
                // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp lines 4740-4741 (1.16.0)
                self.h16_add_to_row_lower(colrow, -HighsCDouble::from(colval) * offset);
                self.h16_add_to_row_upper(colrow, -HighsCDouble::from(colval) * offset);
            } else {
                if self.model.row_lower_[colrow as usize] != -kHighsInf {
                    self.model.row_lower_[colrow as usize] -= colval * offset;
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 8493 at f2ae2b320 (1.15.1 6824)
                    self.x1_row_modified(colrow);
                }

                if self.model.row_upper_[colrow as usize] != kHighsInf {
                    self.model.row_upper_[colrow as usize] -= colval * offset;
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 8494 at f2ae2b320 (1.15.1 6827)
                    self.x1_row_modified(colrow);
                }
            }

            self.add_to_matrix(colrow, staycol, scale * colval);

            // check if this is an equation row and it now has a different size
            self.reinsert_equation(colrow);
        }

        // substitute column in the objective function
        if self.model.col_cost_[substcol as usize] != 0.0 {
            self.model.offset_ += self.model.col_cost_[substcol as usize] * offset;
            debug_assert!(self.model.offset_.is_finite());

            self.model.col_cost_[staycol as usize] +=
                scale * self.model.col_cost_[substcol as usize];

            if self.model.col_cost_[staycol as usize].abs() <= self.options.small_matrix_value {
                self.model.col_cost_[staycol as usize] = 0.0;
            }
            self.model.col_cost_[substcol as usize] = 0.0;
        }
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::substitute (lines 8506-8562)
    /// C++: `Result substitute(HighsInt substcol, HighsInt staycol, double offset, double scale, HighsInt row = -1)`
    /// C++ default arguments: `row = -1`
    ///
    /// HiGHS #3312 (mip_presolve_clique_table): called only from gated sites. #3312 moves the
    /// callers' `markColDeleted` / `removeRow` into `substitute`, after the clique-table update;
    /// the rest of the body is [`HPresolve::substitute_col`] unchanged.
    pub fn substitute_x1b(
        &mut self,
        substcol: HighsInt,
        staycol: HighsInt,
        offset: f64,
        scale: f64,
        row: HighsInt,
    ) -> HPresolveResult {
        hpresolve_checked_call!(
            self.update_clique_table_substitute_col(substcol, staycol, offset, scale)
        );
        self.mark_col_deleted(substcol);
        if row != -1 {
            self.remove_row(row);
        }

        self.substitute_col(substcol, staycol, offset, scale);

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::fixColToLower (lines 6850-6865, 16 lines)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::fixColToLower (lines 8564-8579)
    pub fn fix_col_to_lower(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        let fixval = self.model.col_lower_[col as usize];
        if fixval == -kHighsInf {
            return HPresolveResult::kDualInfeasible;
        }

        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleFixedCol);
        }
        postsolve_stack.fixed_col_at_lower(
            col,
            fixval,
            self.model.col_cost_[col as usize],
            self.get_column_vector(col).iter(self),
        );
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 8575-8578 at 7b890c255
        if self.options.mip_presolve_clique_table {
            let result = self.remove_fixed_col_value_x1b(col, fixval);
            self.analysis_.logging_on_ = logging_on;
            if logging_on {
                self.analysis_stop_presolve_rule_log(kPresolveRuleFixedCol);
            }
            return result;
        }
        self.remove_fixed_col_value(col, fixval);
        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleFixedCol);
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::fixColToUpper (lines 6867-6882, 16 lines)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::fixColToUpper (lines 8581-8596)
    pub fn fix_col_to_upper(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        let fixval = self.model.col_upper_[col as usize];
        if fixval == kHighsInf {
            return HPresolveResult::kDualInfeasible;
        }

        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleFixedCol);
        }
        postsolve_stack.fixed_col_at_upper(
            col,
            fixval,
            self.model.col_cost_[col as usize],
            self.get_column_vector(col).iter(self),
        );
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 8592-8595 at 7b890c255
        if self.options.mip_presolve_clique_table {
            let result = self.remove_fixed_col_value_x1b(col, fixval);
            self.analysis_.logging_on_ = logging_on;
            if logging_on {
                self.analysis_stop_presolve_rule_log(kPresolveRuleFixedCol);
            }
            return result;
        }
        self.remove_fixed_col_value(col, fixval);
        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleFixedCol);
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::fixColToZero (lines 6884-6893, 10 lines)
    pub fn fix_col_to_zero(&mut self, postsolve_stack: &mut HighsPostsolveStack, col: HighsInt) {
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleFixedCol);
        }
        postsolve_stack.fixed_col_at_zero(
            col,
            self.model.col_cost_[col as usize],
            self.get_column_vector(col).iter(self),
        );
        self.remove_fixed_col_value(col, 0.0);
        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleFixedCol);
        }
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::fixColToZero (lines 8598-8608)
    /// C++: `Result fixColToZero(HighsPostsolveStack& postsolve_stack, HighsInt col)`
    ///
    /// HiGHS #3312 (mip_presolve_clique_table): the Result-returning form, called only from
    /// gated sites (1.15.1 `fix_col_to_zero` is unchanged).
    pub fn fix_col_to_zero_x1b(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleFixedCol);
        }
        postsolve_stack.fixed_col_at_zero(
            col,
            self.model.col_cost_[col as usize],
            self.get_column_vector(col).iter(self),
        );
        let result = self.remove_fixed_col_value_x1b(col, 0.0);
        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleFixedCol);
        }
        result
    }

    /// HiGHS latest (#3264, f2ae2b320): HPresolve.cpp::HPresolve::unlinkRow (lines 8562-8570)
    pub fn unlink_row(&mut self, row: HighsInt) {
        debug_assert!(row >= 0);
        debug_assert!((row as usize) < self.rowroot.len());
        self.store_row(row);
        let rowpositions = std::mem::take(&mut self.rowpositions);
        for &rowiter in rowpositions.iter() {
            debug_assert!(self.a_row[rowiter as usize] == row);
            self.unlink(rowiter);
        }
        self.rowpositions = rowpositions;
    }

    /// HiGHS: HPresolve.cpp::HPresolve::removeRow (lines 6895-6906, 12 lines)
    /// HiGHS latest (#3264, f2ae2b320): HPresolve.cpp::HPresolve::removeRow (lines 8572-8579)
    pub fn remove_row(&mut self, row: HighsInt) {
        debug_assert!(row >= 0);
        debug_assert!((row as usize) < self.rowroot.len());
        // first mark the row as logically deleted, so that it is not register as
        // singleton row upon removing its nonzeros
        self.mark_row_deleted(row);
        self.unlink_row(row);
    }

    /// HiGHS: HPresolve.cpp::HPresolve::removeFixedCol (lines 6908-6914, 7 lines)
    ///
    /// C++ overload `removeFixedCol(col)`.
    pub fn remove_fixed_col(&mut self, col: HighsInt) {
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleFixedCol);
        }
        let fixval = self.model.col_lower_[col as usize];
        self.remove_fixed_col_value(col, fixval);
        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleFixedCol);
        }
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::removeFixedCol (lines 8629-8636)
    /// C++: `Result removeFixedCol(HighsInt col)`
    ///
    /// HiGHS #3312 (mip_presolve_clique_table): called only from gated sites.
    pub fn remove_fixed_col_x1b(&mut self, col: HighsInt) -> HPresolveResult {
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleFixedCol);
        }
        let fixval = self.model.col_lower_[col as usize];
        let result = self.remove_fixed_col_value_x1b(col, fixval);
        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleFixedCol);
        }
        result
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::removeFixedCol (lines 8638-8663)
    /// C++: `Result removeFixedCol(HighsInt col, double fixval)`
    ///
    /// HiGHS #3312 (mip_presolve_clique_table): the clique-table update (line 8639), then the
    /// 1.15.1 body ([`HPresolve::remove_fixed_col_value`]) unchanged.
    pub fn remove_fixed_col_value_x1b(&mut self, col: HighsInt, fixval: f64) -> HPresolveResult {
        hpresolve_checked_call!(self.update_clique_table_fixed_col(col, fixval));
        self.remove_fixed_col_value(col, fixval);
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::removeFixedCol (lines 6916-6943, 28 lines)
    ///
    /// C++ overload `removeFixedCol(col, fixval)`.
    pub fn remove_fixed_col_value(&mut self, col: HighsInt, fixval: f64) {
        // mark the column as deleted first so that it is not registered as singleton
        // column upon removing its non-zeros
        self.mark_col_deleted(col);

        let h16_row_duals = self.options.h16_row_duals;
        let mut coliter = self.colhead[col as usize];
        while coliter != -1 {
            let colrow = self.a_row[coliter as usize];
            let colval = self.a_value[coliter as usize];
            debug_assert!(self.a_col[coliter as usize] == col);

            let colpos = coliter;
            coliter = self.a_next[coliter as usize];

            if h16_row_duals {
                // HiGHS 1.16 #3270 (h16_row_duals): HPresolve.cpp lines 4854-4855 (1.16.0)
                self.h16_add_to_row_lower(colrow, -HighsCDouble::from(colval) * fixval);
                self.h16_add_to_row_upper(colrow, -HighsCDouble::from(colval) * fixval);
            } else {
                if self.model.row_lower_[colrow as usize] != -kHighsInf {
                    self.model.row_lower_[colrow as usize] -= colval * fixval;
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 8602 at f2ae2b320 (1.15.1 6930)
                    self.x1_row_modified(colrow);
                }

                if self.model.row_upper_[colrow as usize] != kHighsInf {
                    self.model.row_upper_[colrow as usize] -= colval * fixval;
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 8603 at f2ae2b320 (1.15.1 6933)
                    self.x1_row_modified(colrow);
                }
            }

            self.unlink(colpos);

            self.reinsert_equation(colrow);
        }

        self.model.offset_ += self.model.col_cost_[col as usize] * fixval;
        debug_assert!(self.model.offset_.is_finite());
        self.model.col_cost_[col as usize] = 0.0;
    }

    /// HiGHS: HPresolve.cpp::HPresolve::removeRowSingletons (lines 6945-6958, 14 lines)
    /// HiGHS 1.16 (#3139): HPresolve.cpp::HPresolve::removeRowSingletons (lines 4867-4880)
    pub fn remove_row_singletons(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        let h16_singleton_loops = self.options.h16_singleton_loops;
        let mut i: usize = 0;
        // HiGHS 1.16 #3139 (h16_singleton_loops): HPresolve.cpp line 4869 (1.16.0)
        while if h16_singleton_loops {
            i < self.singleton_rows.len()
        } else {
            i != self.singleton_rows.len()
        } {
            let row = self.singleton_rows[i];
            if self.row_deleted[row as usize] == 0 && self.rowsize[row as usize] <= 1 {
                // row presolve will delegate to rowSingleton() if the row size is 1
                // if the singleton row has become empty it will also remove the row
                hpresolve_checked_call!(self.row_presolve(postsolve_stack, row));
            }
            i += 1;
        }

        self.singleton_rows.clear();

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::presolveColSingletons (lines 6960-6974, 15 lines)
    /// HiGHS 1.16 (#3139): HPresolve.cpp::HPresolve::presolveColSingletons (lines 4882-4896)
    pub fn presolve_col_singletons(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        let h16_singleton_loops = self.options.h16_singleton_loops;
        let mut i: usize = 0;
        // HiGHS 1.16 #3139 (h16_singleton_loops): HPresolve.cpp line 4884 (1.16.0)
        while if h16_singleton_loops {
            i < self.singleton_columns.len()
        } else {
            i != self.singleton_columns.len()
        } {
            let col = self.singleton_columns[i];
            if self.col_deleted[col as usize] == 0 {
                hpresolve_checked_call!(self.col_presolve(postsolve_stack, col));
            }
            i += 1;
        }
        let col_deleted = &self.col_deleted;
        let colsize = &self.colsize;
        self.singleton_columns
            .retain(|&col| !(col_deleted[col as usize] != 0 || colsize[col as usize] > 1));

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::presolveChangedRows (lines 6976-6988, 13 lines)
    pub fn presolve_changed_rows(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        // `changedRows.swap(changedRowIndices)`: the list is taken, and new
        // entries pushed during the loop go to the (now empty) member.
        let changed_rows = std::mem::take(&mut self.changed_row_indices);
        for row in changed_rows {
            if self.row_deleted[row as usize] != 0 {
                continue;
            }
            hpresolve_checked_call!(self.row_presolve(postsolve_stack, row));
            self.changed_row_flag[row as usize] = self.row_deleted[row as usize];
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::presolveChangedCols (lines 6990-7002, 13 lines)
    pub fn presolve_changed_cols(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        // `changedCols.swap(changedColIndices)`: see `presolve_changed_rows`.
        let changed_cols = std::mem::take(&mut self.changed_col_indices);
        for col in changed_cols {
            if self.col_deleted[col as usize] != 0 {
                continue;
            }
            hpresolve_checked_call!(self.col_presolve(postsolve_stack, col));
            self.changed_col_flag[col as usize] = self.col_deleted[col as usize];
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::removeDoubletonEquations (lines 7004-7021, 18 lines)
    ///
    /// `auto eq = equations.begin(); ... ++eq` over the key (module doc).
    pub fn remove_doubleton_equations(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        let mut eq: Option<(HighsInt, HighsInt)> = self.equations.iter().next().copied();
        while let Some(key) = eq {
            let eqrow = key.1;
            debug_assert!(self.row_deleted[eqrow as usize] == 0);
            debug_assert!(key.0 == self.rowsize[eqrow as usize]);
            debug_assert!(self.is_equation(eqrow));
            if self.rowsize[eqrow as usize] > 2 {
                return HPresolveResult::kOk;
            }
            hpresolve_checked_call!(self.row_presolve(postsolve_stack, eqrow));
            if self.row_deleted[eqrow as usize] != 0 {
                eq = self.equations.iter().next().copied();
            } else {
                eq = self
                    .equations
                    .range((std::ops::Bound::Excluded(key), std::ops::Bound::Unbounded))
                    .next()
                    .copied();
            }
        }

        HPresolveResult::kOk
    }
}

/// Time limit for the dependent-equations factorization (`factor.set_time_limit`).
///
/// HiGHS 1.16 #3167/#3398 (h16_dep_eq_time_cap): HPresolveDependent.cpp lines 91-93 (1.16.0),
/// kMaxDependentEquationsTime = 100; 1.15.1 HPresolve.cpp:6528 used 1000
/// (ff138f533 introduced the constant at 10, settled at 100 in 31b27b7eb).
/// With `cap_100` false the value is bit-identical to the 1.15.1 inline expression.
fn dependent_equations_time_limit(time_limit: f64, cap_100: bool) -> f64 {
    let max_time = if cap_100 { 100.0 } else { 1000.0 };
    cpp_max(1.0, cpp_min(0.01 * time_limit, max_time))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two columns, two rows: row 0 = `x0` (a singleton row, in [2, 5]) and
    /// row 1 = `x0 + x1` in [0, 10]; `x1` is a singleton column of row 1.
    fn singleton_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 3];
        lp.a_matrix_.index_ = vec![0, 1, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0];
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 10.0];
        lp.integrality_ = vec![HighsVarType::kContinuous, HighsVarType::kContinuous];
        lp.row_lower_ = vec![2.0, 0.0];
        lp.row_upper_ = vec![5.0, 10.0];
        lp
    }

    /// Runs the row and column singleton loops with the option set to `on` and returns
    /// the deleted-row and deleted-column flags afterwards.
    fn run_singleton_loops(on: bool) -> (Vec<u8>, Vec<u8>) {
        let mut options = HighsOptions::default();
        options.h16_singleton_loops = on;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(
            singleton_lp(),
            &options,
            options.presolve_reduction_limit,
            None
        ));
        // ok_set_input does not queue these singletons; queue them as unlink does
        // (links.rs) so the loops have work to do.
        presolve.singleton_rows.push(0);
        presolve.singleton_columns.push(1);
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(2, 2);
        assert!(matches!(
            presolve.remove_row_singletons(&mut stack),
            HPresolveResult::kOk
        ));
        assert!(presolve.singleton_rows.is_empty());
        assert!(matches!(
            presolve.presolve_col_singletons(&mut stack),
            HPresolveResult::kOk
        ));
        (presolve.row_deleted.clone(), presolve.col_deleted.clone())
    }

    #[test]
    fn h16_singleton_loops_on_terminates() {
        // With the option on the loops run `i < len`: the singleton row is processed
        // and removed, and remove_row_singletons empties its queue.
        let (row_deleted, _) = run_singleton_loops(true);
        assert!(row_deleted[0] != 0);
    }

    #[test]
    fn h16_singleton_loops_off_matches_1_15_1() {
        // With the option off the loops run `i != len` (1.15.1). On this input they
        // terminate, and they delete the same rows and columns as the option-on path.
        let off = run_singleton_loops(false);
        assert!(off.0[0] != 0);
        assert_eq!(off, run_singleton_loops(true));
    }
}

#[cfg(test)]
mod h16_subst_int_tests {
    use super::*;

    /// Row 0: `x0 + x1 = 1`, x0 integer (substcol), x1 implied integer (staycol).
    fn two_col_presolve(h16_subst_int: bool) -> HPresolve {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0; 2];
        lp.col_lower_ = vec![0.0; 2];
        lp.col_upper_ = vec![1.0; 2];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![1.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.h16_subst_int = h16_subst_int;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, None));
        // LP mode makes every column continuous. Set the MIP integrality afterwards and
        // the row counters the MIP path would have built (x0: integer, x1: implied integer).
        presolve.model.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kImplicitInteger];
        presolve.rowsize_integer[0] = 1;
        presolve.rowsize_impl_int[0] = 1;
        presolve
    }

    /// Integer substcol with an implied-integer staycol: the staycol is upgraded to
    /// kInteger and the row counters move one entry from implied to explicit integer.
    #[test]
    fn h16_subst_int_upgrades_implied_stay_col() {
        let mut presolve = two_col_presolve(true);
        presolve.substitute_col(0, 1, 0.0, 1.0);
        assert_eq!(presolve.model.integrality_[1], HighsVarType::kInteger);
        // row 0 now holds only x1, an integer column
        assert_eq!(presolve.rowsize_integer[0], 1);
        assert_eq!(presolve.rowsize_impl_int[0], 0);
        assert_eq!(presolve.rowsize[0], 1);
    }

    /// Option off: the staycol keeps kImplicitInteger and the counters follow 1.15.1.
    #[test]
    fn h16_subst_int_off_keeps_implied_stay_col() {
        let mut presolve = two_col_presolve(false);
        presolve.substitute_col(0, 1, 0.0, 1.0);
        assert_eq!(
            presolve.model.integrality_[1],
            HighsVarType::kImplicitInteger
        );
        assert_eq!(presolve.rowsize_integer[0], 0);
        assert_eq!(presolve.rowsize_impl_int[0], 1);
        assert_eq!(presolve.rowsize[0], 1);
    }
}

#[cfg(test)]
mod x1b_clqp_c_tests {
    use super::*;

    /// Row 0: `x0 - x1 = 0`; row 1: `x0 + x2 <= 1`.
    fn three_col_presolve() -> HPresolve {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![1.0, 0.0, 0.0];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![1.0; 3];
        lp.row_lower_ = vec![0.0, -kHighsInf];
        lp.row_upper_ = vec![0.0, 1.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 2, 3, 4];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, -1.0, 1.0];
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_presolve_clique_table = true;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, None));
        presolve
    }

    /// On: substitute_x1b with row 0 deletes x0 and row 0 itself, and moves x0 of row 1 to
    /// x1: the same state as the 1.15.1 caller sequence markColDeleted, removeRow,
    /// substitute (no MIP solver, so the clique-table update is a no-op).
    #[test]
    fn x1b_clqp_c_substitute_x1b_marks_and_removes() {
        let mut on = three_col_presolve();
        assert_eq!(on.substitute_x1b(0, 1, 0.0, 1.0, 0), HPresolveResult::kOk);
        assert_ne!(on.col_deleted[0], 0);
        assert_ne!(on.row_deleted[0], 0);
        assert_eq!(on.row_deleted[1], 0);
        assert_eq!(on.rowsize[1], 2);
        assert_eq!(on.colsize[0], 0);
        assert_eq!(on.model.col_cost_[1], 1.0);

        let mut off = three_col_presolve();
        off.mark_col_deleted(0);
        off.remove_row(0);
        off.substitute_col(0, 1, 0.0, 1.0);
        assert_eq!(on.col_deleted, off.col_deleted);
        assert_eq!(on.row_deleted, off.row_deleted);
        assert_eq!(on.rowsize, off.rowsize);
        assert_eq!(on.colsize, off.colsize);
        assert_eq!(on.model.col_cost_, off.model.col_cost_);

        // row = -1 keeps the row
        let mut keep = three_col_presolve();
        assert_eq!(
            keep.substitute_x1b(0, 1, 0.0, 1.0, -1),
            HPresolveResult::kOk
        );
        assert_eq!(keep.row_deleted[0], 0);
    }

    /// remove_fixed_col_x1b without a MIP solver is the 1.15.1 removeFixedCol plus kOk.
    #[test]
    fn x1b_clqp_c_remove_fixed_col_x1b_lp() {
        let mut on = three_col_presolve();
        on.model.col_upper_[2] = 0.0;
        assert_eq!(on.remove_fixed_col_x1b(2), HPresolveResult::kOk);
        let mut off = three_col_presolve();
        off.model.col_upper_[2] = 0.0;
        off.remove_fixed_col(2);
        assert_eq!(on.col_deleted, off.col_deleted);
        assert_eq!(on.rowsize, off.rowsize);
        assert_eq!(on.model.row_upper_, off.model.row_upper_);
    }
}

#[cfg(test)]
mod h16_dep_eq_time_cap_tests {
    use super::*;

    /// Option on (1.16.0): cap at kMaxDependentEquationsTime = 100. time_limit 5000
    /// gives 0.01 * 5000 = 50 < 100, so the cap does not bite there; 50 gives 0.5,
    /// clamped up to the floor of 1.
    #[test]
    fn h16_dep_eq_time_cap_on_caps_at_100() {
        assert_eq!(dependent_equations_time_limit(f64::INFINITY, true), 100.0);
        assert_eq!(dependent_equations_time_limit(50_000.0, true), 100.0);
        assert_eq!(dependent_equations_time_limit(5_000.0, true), 50.0);
        assert_eq!(dependent_equations_time_limit(50.0, true), 1.0);
    }

    /// Option off (1.15.1): cap at the literal 1000, bit-identical to the inline
    /// expression `cpp_max(1.0, cpp_min(0.01 * time_limit, 1000.0))`.
    #[test]
    fn h16_dep_eq_time_cap_off_caps_at_1000() {
        for time_limit in [f64::INFINITY, 50_000.0, 5_000.0, 50.0] {
            assert_eq!(
                dependent_equations_time_limit(time_limit, false),
                cpp_max(1.0, cpp_min(0.01 * time_limit, 1000.0))
            );
        }
        assert_eq!(dependent_equations_time_limit(f64::INFINITY, false), 1000.0);
        assert_eq!(dependent_equations_time_limit(50_000.0, false), 500.0);
        assert_eq!(dependent_equations_time_limit(5_000.0, false), 50.0);
        assert_eq!(dependent_equations_time_limit(50.0, false), 1.0);
    }
}
