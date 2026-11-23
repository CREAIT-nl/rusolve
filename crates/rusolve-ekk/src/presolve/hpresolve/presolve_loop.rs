//! HiGHS source: presolve/HPresolve.cpp (lines 5698-6045)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Implied bounds helpers and the main presolve loop. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::computeImpliedLowerBound (lines 5698-5705, 8 lines)
    ///
    /// C++ default arguments: `boundCol = -1`, `boundColValue = kHighsInf`, `boundColCoeffPattern = 0`.
    /// `double lowerBound;` (uninitialised, always written) -> `let mut lower_bound = 0.0;`.
    pub fn compute_implied_lower_bound(
        &mut self,
        col: HighsInt,
        bound_col: HighsInt,
        bound_col_value: f64,
        bound_col_coeff_pattern: HighsInt,
    ) -> f64 {
        let mut lower_bound = 0.0;
        self.compute_col_bounds(
            col,
            bound_col,
            bound_col_value,
            bound_col_coeff_pattern,
            Some(&mut lower_bound),
            None,
            None,
            None,
        );
        lower_bound
    }

    /// HiGHS: HPresolve.cpp::HPresolve::computeImpliedUpperBound (lines 5707-5714, 8 lines)
    ///
    /// C++ default arguments: `boundCol = -1`, `boundColValue = kHighsInf`, `boundColCoeffPattern = 0`.
    pub fn compute_implied_upper_bound(
        &mut self,
        col: HighsInt,
        bound_col: HighsInt,
        bound_col_value: f64,
        bound_col_coeff_pattern: HighsInt,
    ) -> f64 {
        let mut upper_bound = 0.0;
        self.compute_col_bounds(
            col,
            bound_col,
            bound_col_value,
            bound_col_coeff_pattern,
            None,
            Some(&mut upper_bound),
            None,
            None,
        );
        upper_bound
    }

    /// HiGHS: HPresolve.cpp::HPresolve::computeWorstCaseLowerBound (lines 5716-5723, 8 lines)
    ///
    /// C++ default arguments: `boundCol = -1`, `boundColValue = kHighsInf`, `boundColCoeffPattern = 0`.
    pub fn compute_worst_case_lower_bound(
        &mut self,
        col: HighsInt,
        bound_col: HighsInt,
        bound_col_value: f64,
        bound_col_coeff_pattern: HighsInt,
    ) -> f64 {
        let mut lower_bound = 0.0;
        self.compute_col_bounds(
            col,
            bound_col,
            bound_col_value,
            bound_col_coeff_pattern,
            None,
            None,
            Some(&mut lower_bound),
            None,
        );
        lower_bound
    }

    /// HiGHS: HPresolve.cpp::HPresolve::computeWorstCaseUpperBound (lines 5725-5732, 8 lines)
    ///
    /// C++ default arguments: `boundCol = -1`, `boundColValue = kHighsInf`, `boundColCoeffPattern = 0`.
    pub fn compute_worst_case_upper_bound(
        &mut self,
        col: HighsInt,
        bound_col: HighsInt,
        bound_col_value: f64,
        bound_col_coeff_pattern: HighsInt,
    ) -> f64 {
        let mut upper_bound = 0.0;
        self.compute_col_bounds(
            col,
            bound_col,
            bound_col_value,
            bound_col_coeff_pattern,
            None,
            None,
            None,
            Some(&mut upper_bound),
        );
        upper_bound
    }

    /// HiGHS: HPresolve.cpp::HPresolve::initialRowAndColPresolve (lines 5734-5758, 25 lines)
    pub fn initial_row_and_col_presolve(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        // do a full scan over the rows as the singleton arrays and the changed row
        // arrays are not initialized, also unset changedRowFlag so that the row will
        // be added to the changed row vector when it is changed after it was
        // processed
        for row in 0..self.model.num_row_ {
            if self.row_deleted[row as usize] != 0 {
                continue;
            }
            hpresolve_checked_call!(self.row_presolve(postsolve_stack, row));
            self.changed_row_flag[row as usize] = 0;
        }

        // same for the columns
        for col in 0..self.model.num_col_ {
            if self.col_deleted[col as usize] != 0 {
                continue;
            }
            // round and update bounds
            if self.model.integrality_[col as usize] != HighsVarType::kContinuous {
                hpresolve_checked_call!(self.change_col_bounds(
                    col,
                    self.model.col_lower_[col as usize],
                    self.model.col_upper_[col as usize]
                ));
            }
            hpresolve_checked_call!(self.col_presolve(postsolve_stack, col));
            self.changed_col_flag[col as usize] = 0;
        }

        self.check_limits(postsolve_stack)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::fastPresolveLoop (lines 5760-5778, 19 lines)
    pub fn fast_presolve_loop(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::FAST_LOOP);
        loop {
            self.store_current_problem_size();

            hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));

            hpresolve_checked_call!(self.presolve_changed_rows(postsolve_stack));

            hpresolve_checked_call!(self.remove_doubleton_equations(postsolve_stack));

            hpresolve_checked_call!(self.presolve_col_singletons(postsolve_stack));

            hpresolve_checked_call!(self.presolve_changed_cols(postsolve_stack));

            // do { ... } while (problemSizeReduction() > 0.01);
            if !(self.problem_size_reduction() > 0.01) {
                break;
            }
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::presolve (lines 5780-6045, 266 lines)
    ///
    /// The rule-order loop: every `analysis_.allow_rule_` gate, `problemSizeReduction()` threshold
    /// and `continue` decides the reduction sequence. `analysis_.setup(this->model, this->options, ...)` ->
    /// `self.analysis_.setup(&self.model, &self.options, silent)`; `analysePresolveRuleLog()` ->
    /// `self.analysis_.analyse_presolve_rule_log(false, &self.model, &self.options, self.num_deleted_rows, self.num_deleted_cols)`.
    /// `toCSC(model->a_matrix_...)` before `shrinkProblem`: take the three vectors (module doc).
    /// `presolveReturn` and the other lambdas: see the module doc.
    pub fn presolve(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HPresolveResult {
        // HiGHS 1.16 #3167 (h16_presolve_rules): HPresolve.cpp::presolve lines 3920-4256 (1.16.0),
        // in presolve_h16.rs.
        let h16_presolve_rules = self.options.h16_presolve_rules;
        if h16_presolve_rules {
            return self.h16_presolve(postsolve_stack);
        }

        // C++ lambda `presolveReturn` (captures `this`).
        fn presolve_return(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
        ) -> HPresolveResult {
            if this.mipsolver.is_some() {
                hpresolve_checked_call!(this.scale_mip(postsolve_stack));
            }

            // analysePresolveRuleLog() should return true - no errors
            debug_assert!(this.analysis_.analyse_presolve_rule_log(
                false,
                &this.model,
                &this.options,
                this.num_deleted_rows,
                this.num_deleted_cols
            ));
            // Possibly report presolve log
            this.analysis_.analyse_presolve_rule_log(
                true,
                &this.model,
                &this.options,
                this.num_deleted_rows,
                this.num_deleted_cols,
            );
            HPresolveResult::kOk
        }

        // C++ lambda `report` (captures `this` and `silent`).
        fn report(this: &HPresolve, silent: bool) {
            if !silent {
                let num_col = this.model.num_col_ - this.num_deleted_cols;
                let num_row = this.model.num_row_ - this.num_deleted_rows;
                let num_nonz = (this.a_value.len() - this.freeslots.len()) as HighsInt;
                // Only read the run time if it's to be printed
                let run_time = if this.options.output_flag {
                    this.timer.as_ref().unwrap().read()
                } else {
                    0.0
                };
                // highsTimeSecondToString: std::to_string(double) ("%f") in debug builds,
                // std::to_string(static_cast<int>(time)) under NDEBUG, plus "s".
                let mut time_str = if cfg!(debug_assertions) {
                    format!("{:.6}s", run_time)
                } else {
                    format!("{}s", run_time as HighsInt)
                };
                if this.options.timeless_log {
                    time_str = String::new();
                }
                crate::highs_log_user!(
                    &this.options.log_options,
                    HighsLogType::kInfo,
                    "{} rows, {} cols, {} nonzeros {}\n",
                    num_row,
                    num_col,
                    num_nonz,
                    time_str
                );
            }
        }

        // The shrink-and-rebuild sequence appears twice in C++ (main loop and dependent
        // equations): shrinkProblem, then toCSC/okFromCSC on the model's matrix.
        fn shrink_and_rebuild(this: &mut HPresolve, postsolve_stack: &mut HighsPostsolveStack) {
            this.shrink_problem(postsolve_stack);

            let mut a_value = std::mem::take(&mut this.model.a_matrix_.value_);
            let mut a_index = std::mem::take(&mut this.model.a_matrix_.index_);
            let mut a_start = std::mem::take(&mut this.model.a_matrix_.start_);
            this.to_csc(&mut a_value, &mut a_index, &mut a_start);
            this.model.a_matrix_.value_ = a_value;
            this.model.a_matrix_.index_ = a_index;
            this.model.a_matrix_.start_ = a_start;

            let a_value = std::mem::take(&mut this.model.a_matrix_.value_);
            let a_index = std::mem::take(&mut this.model.a_matrix_.index_);
            let a_start = std::mem::take(&mut this.model.a_matrix_.start_);
            this.ok_from_csc(&a_value, &a_index, &a_start);
            this.model.a_matrix_.value_ = a_value;
            this.model.a_matrix_.index_ = a_index;
            this.model.a_matrix_.start_ = a_start;
        }

        // for the inner most loop we take the order roughly from the old presolve
        // but we nest the rounds with a new outer loop which layers the newer
        // presolvers
        //    fast presolve loop
        //        - empty, forcing and dominated rows and row singletons immediately
        //        after each forcing row
        //        - doubleton equations and row singletons immediately after each
        //        successful substitution
        //        - col singletons (can this introduce row singletons? If yes then
        //        immediately remove)
        //        - empty, dominated and weakly dominated columns
        //        - row singletons
        //        - if( !has enough changes ) stop
        // main loop
        //    - fast presolve loop
        //    - parallel rows and columns
        //    - if (changes found) fast presolve loop
        //    - aggregator // add limit that catches many substitutions but stops when
        //    many failures, do not run exhaustively as now
        //    - if (changes found) start main loop from beginning
        //    - primal and dual matrix sparsification
        //    - if (changes found) fast presolve loop
        //    - stop

        // convert model to minimization problem
        if self.model.sense_ == ObjSense::kMaximize {
            for i in 0..self.model.num_col_ {
                self.model.col_cost_[i as usize] = -self.model.col_cost_[i as usize];
            }

            self.model.offset_ = -self.model.offset_;
            debug_assert!(self.model.offset_.is_finite());
            self.model.sense_ = ObjSense::kMinimize;
        }

        let silent = self.silent_log();
        if self.options.presolve != kHighsOffString {
            if !silent {
                crate::highs_log_user!(
                    &self.options.log_options,
                    HighsLogType::kInfo,
                    "Presolving model\n"
                );
            }
        }
        // Set up the logic to allow presolve rules, and logging for their
        // effectiveness
        self.analysis_.setup(&self.model, &self.options, silent);

        // HiGHS 1.16 #3342/#3348 (h16_dual_postsolve_rules): HPresolve.cpp::chooseRules (1.16.0
        // lines 503-534, the mipsolver == nullptr branch) switches off the rules that are
        // incompatible with dual or basis postsolve. The Sparsify basis rule-off belongs to
        // the deferred chooseRules rework (#3167/#3334; try_sparsify already covers it) and
        // ZeroCostSingleton is gated at its call site in singletons.rs (#2962).
        let h16_dual_postsolve_rules = self.options.h16_dual_postsolve_rules;
        if h16_dual_postsolve_rules && self.mipsolver.is_none() {
            let may_require_dual_postsolve = h16_may_require_dual_postsolve(&self.options);
            let may_require_basis_postsolve = h16_may_require_basis_postsolve(&self.options);
            if may_require_dual_postsolve {
                // Cannot use kPresolveRuleDualFixing rule (see
                // test-dual-fixing-primal-dual-postsolve)
                self.analysis_.allow_rule_[kPresolveRuleDualFixing as usize] = false;
            }
            if may_require_basis_postsolve || may_require_dual_postsolve {
                // Cannot use dominated column rule for either dual or basis
                // postsolve (see issue-3342)
                self.analysis_.allow_rule_[kPresolveRuleDominatedCol as usize] = false;
            }
        }

        if self.options.presolve != kHighsOffString {
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 6606 at 7b890c255
            // (removed: the flag is scoped to runProbing / enumerateSolutions)
            if self.mipsolver.is_some() && !self.options.mip_presolve_clique_table {
                self.mip_mut().mipdata_.cliquetable.set_presolve_flag(true);
            }

            // Need to check for time-out in checkLimits, so make sure that
            // the timer is well defined, and that its total time clock is
            // running
            debug_assert!(self.timer.is_some());
            debug_assert!(self.timer.as_ref().unwrap().running());

            if self.options.presolve_rule_test != 0 {
                hpresolve_checked_call!(self.presolve_rule_test(postsolve_stack));
                return presolve_return(self, postsolve_stack);
            }
            hpresolve_checked_call!(self.initial_row_and_col_presolve(postsolve_stack));
            // rusolve mip_presolve_xor_gf2
            hpresolve_checked_call!(self.xor_gf2_presolve(postsolve_stack));
            // rusolve mip_presolve_clique_varbound
            hpresolve_checked_call!(self.clique_varbound_presolve(postsolve_stack));

            let mut num_parallel_row_col_calls: HighsInt = 0;
            // ReductionType::kEqualityRowAddition(s) has no basis postsolve,
            // so can only be used when basis postsolve is not required - when
            // the problem is a MIP, IPM is run without crossover, or when
            // PDLP is used. However, if the LP is reduced to empty, the basis
            // must not be formed in the case of IPM without crossover or PDLP
            //
            // HighsOptions::lp_presolve_requires_basis_postsolve is true by
            // default, and only switched to false if the solver is IPM
            // without crossover or PDLP
            // ENABLE_SPARSIFY_FOR_LP is 0: the #else branch is the only one ported.
            let mut try_sparsify =
                self.mipsolver.is_some() || !self.options.lp_presolve_requires_basis_postsolve;
            let mut try_probing = self.mipsolver.is_some();
            let mut num_cliques_before_probing: HighsInt = -1;
            let mut domcol_after_probing_called = false;
            let mut dependent_equations_called = self.mipsolver.is_some();
            let mut last_print_size: HighsInt = kHighsIInf;

            // Start of main presolve loop
            //
            loop {
                let curr_size = self.model.num_col_ - self.num_deleted_cols + self.model.num_row_
                    - self.num_deleted_rows;
                if (curr_size as f64) < 0.85 * (last_print_size as f64) {
                    last_print_size = curr_size;
                    report(self, silent);
                }

                hpresolve_checked_call!(self.fast_presolve_loop(postsolve_stack));

                self.store_current_problem_size();

                // when presolving after a restart the clique table and implication
                // structure may contain substitutions which we apply directly before
                // running the aggregator as they might lose validity otherwise
                if self.mipsolver.is_some() {
                    let mut num_del_col: HighsInt = 0;
                    hpresolve_checked_call!(
                        self.apply_conflict_graph_substitutions(postsolve_stack, &mut num_del_col)
                    );
                }

                if self.analysis_.allow_rule_[kPresolveRuleAggregator as usize] {
                    hpresolve_checked_call!(self.aggregator(postsolve_stack));
                }

                // rusolve mip_presolve_eq_subst
                if self.options.mip_presolve_eq_subst {
                    self.store_current_problem_size();
                    hpresolve_checked_call!(self.eq_subst_bound_rows(postsolve_stack));
                    if self.problem_size_reduction() > 0.0 {
                        hpresolve_checked_call!(self.fast_presolve_loop(postsolve_stack));
                    }
                }

                if self.problem_size_reduction() > 0.05 {
                    continue;
                }

                if try_sparsify
                    && !self.ti_structure_protect()
                    && self.analysis_.allow_rule_[kPresolveRuleSparsify as usize]
                {
                    let num_nz = self.num_nonzeros();
                    hpresolve_checked_call!(self.sparsify(postsolve_stack));
                    let nz_reduction = 100.0 * (1.0 - (self.num_nonzeros() as f64 / num_nz as f64));

                    if nz_reduction > 0.0 {
                        crate::highs_log_dev!(
                            &self.options.log_options,
                            HighsLogType::kInfo,
                            "Sparsify removed {:.1}% of nonzeros\n",
                            nz_reduction
                        );

                        // #1710 exposes that this should not be
                        //
                        // fastPresolveLoop(postsolve_stack);
                        //
                        // but
                        hpresolve_checked_call!(self.fast_presolve_loop(postsolve_stack));
                    }
                    try_sparsify = false;
                }

                if self.analysis_.allow_rule_[kPresolveRuleParallelRowsAndCols as usize]
                    && num_parallel_row_col_calls < 5
                {
                    if self.shrink_problem_enabled
                        && (self.num_deleted_cols >= self.model.num_col_ / 2
                            || self.num_deleted_rows >= self.model.num_row_ / 2)
                    {
                        shrink_and_rebuild(self, postsolve_stack);
                    }
                    self.store_current_problem_size();
                    hpresolve_checked_call!(self.detect_parallel_rows_and_cols(postsolve_stack));
                    num_parallel_row_col_calls += 1;
                    if self.problem_size_reduction() > 0.05 {
                        continue;
                    }
                }

                hpresolve_checked_call!(self.fast_presolve_loop(postsolve_stack));

                if self.mipsolver.is_some() {
                    let mut num_strengthened: HighsInt = -1;
                    hpresolve_checked_call!(
                        self.strengthen_inequalities(postsolve_stack, &mut num_strengthened)
                    );
                    debug_assert!(num_strengthened >= 0);
                    if num_strengthened > 0 {
                        crate::highs_log_dev!(
                            &self.options.log_options,
                            HighsLogType::kInfo,
                            "Strengthened {} coefficients\n",
                            num_strengthened
                        );
                    }
                }

                hpresolve_checked_call!(self.fast_presolve_loop(postsolve_stack));

                if self.mipsolver.is_some() && num_cliques_before_probing == -1 {
                    num_cliques_before_probing = self.mip().mipdata_.cliquetable.num_cliques();
                    self.store_current_problem_size();
                    hpresolve_checked_call!(self.dominated_columns(postsolve_stack));
                    if self.problem_size_reduction() > 0.0 {
                        hpresolve_checked_call!(self.fast_presolve_loop(postsolve_stack));
                    }
                    if self.problem_size_reduction() > 0.05 {
                        continue;
                    }
                }

                // enumerate solutions
                if self.mipsolver.is_some()
                    && self.analysis_.allow_rule_[kPresolveRuleEnumeration as usize]
                {
                    self.store_current_problem_size();
                    hpresolve_checked_call!(self.enumerate_solutions(postsolve_stack));
                    if self.problem_size_reduction() > 0.05 {
                        continue;
                    }
                }

                if try_probing && self.analysis_.allow_rule_[kPresolveRuleProbing as usize] {
                    hpresolve_checked_call!(self.detect_implied_integers());
                    self.store_current_problem_size();
                    hpresolve_checked_call!(self.run_probing(postsolve_stack));
                    try_probing = self.probing_contingent > self.num_probed as i64
                        && (self.problem_size_reduction() > 1.0 || self.probing_early_abort);
                    try_sparsify = true;
                    if self.problem_size_reduction() > 0.05 || try_probing {
                        continue;
                    }
                    hpresolve_checked_call!(self.fast_presolve_loop(postsolve_stack));
                }

                if !dependent_equations_called {
                    if self.shrink_problem_enabled
                        && (self.num_deleted_cols >= self.model.num_col_ / 2
                            || self.num_deleted_rows >= self.model.num_row_ / 2)
                    {
                        shrink_and_rebuild(self, postsolve_stack);
                    }
                    self.store_current_problem_size();
                    if self.analysis_.allow_rule_[kPresolveRuleDependentEquations as usize] {
                        hpresolve_checked_call!(self.remove_dependent_equations(postsolve_stack));
                        dependent_equations_called = true;
                    }
                    if self.analysis_.allow_rule_[kPresolveRuleDependentFreeCols as usize] {
                        hpresolve_checked_call!(self.remove_dependent_free_cols(postsolve_stack));
                    }
                    if self.problem_size_reduction() > 0.05 {
                        continue;
                    }
                }

                if self.mipsolver.is_some()
                    && self.mip().mipdata_.cliquetable.num_cliques() > num_cliques_before_probing
                    && !domcol_after_probing_called
                {
                    domcol_after_probing_called = true;
                    self.store_current_problem_size();
                    hpresolve_checked_call!(self.dominated_columns(postsolve_stack));
                    if self.problem_size_reduction() > 0.0 {
                        hpresolve_checked_call!(self.fast_presolve_loop(postsolve_stack));
                    }
                    if self.problem_size_reduction() > 0.05 {
                        continue;
                    }
                }

                break;
            }

            // Now consider removing slacks
            if self.options.presolve_remove_slacks {
                hpresolve_checked_call!(self.remove_slacks(postsolve_stack));
            }

            report(self, silent);
        } else {
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "\nPresolve is switched off\n"
            );
        }

        presolve_return(self, postsolve_stack)
    }
}

/// HiGHS 1.16 (#3342): HighsOptions.cpp::mayRequireDualPostsolve (lines 1160-1179)
/// (1.16.0 lines 1160-1179). Dual postsolve is needed for PDLP, and for IPM-family
/// solvers when crossover is off or chosen.
pub(super) fn h16_may_require_dual_postsolve(options: &HighsOptions) -> bool {
    if use_pdlp(&options.solver) {
        return true;
    }
    if use_ipm(&options.solver)
        && (options.run_crossover == kHighsOffString || options.run_crossover == kHighsChooseString)
    {
        return true;
    }
    false
}

/// HiGHS 1.16 (#3342): HighsOptions.cpp::mayRequireBasisPostsolve (lines 1160-1179)
/// (1.16.0 lines 1160-1179). Basis postsolve is needed unless PDLP, or IPM-family
/// with crossover off.
pub(super) fn h16_may_require_basis_postsolve(options: &HighsOptions) -> bool {
    if use_pdlp(&options.solver) {
        return false;
    }
    if use_ipm(&options.solver) && options.run_crossover == kHighsOffString {
        return false;
    }
    true
}

#[cfg(test)]
mod h16_dual_postsolve_tests {
    use super::*;
    use crate::prelude::*;

    /// issue-3342 (3348.diff, TestPresolve.cpp): min -6w + 2x - y + 2z;
    /// 3 <= 2x + 2y - 2z; -4 <= w - 2z <= 0; w in [-6, 2], -4 <= x, y <= 2, z <= 1.
    /// Optimum (C++ comment): w(2; -4), x(1/2; 0), y(2; -3), z(1; 0), r0(3; 1),
    /// r1(0; -2), objective -11.
    fn issue_3342_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 4,
            num_row_: 2,
            col_cost_: vec![-6.0, 2.0, -1.0, 2.0],
            col_lower_: vec![-6.0, -4.0, -kHighsInf, -kHighsInf],
            col_upper_: vec![2.0, kHighsInf, 2.0, 1.0],
            row_lower_: vec![3.0, -4.0],
            row_upper_: vec![kHighsInf, 0.0],
            sense_: ObjSense::kMinimize,
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 1, 2, 3, 5];
        lp.a_matrix_.index_ = vec![1, 0, 0, 0, 1];
        lp.a_matrix_.value_ = vec![1.0, 2.0, 2.0, -2.0, -2.0];
        lp
    }

    /// test-dual-fixing-primal-dual-postsolve (3342.diff, TestPresolve.cpp):
    /// min -8w + z; v - 2z <= 8; w - x - y - z = -5; 0 <= v - w;
    /// v free, w <= 7, x in [0, 9], y <= 3, z in [0, 1]. Optimum objective -56 with
    /// w = 7, x = 9, y = 3, z = 0 (unique); dual fixing of z gives the dual
    /// infeasible z(0; -2.2) after IPM without crossover.
    fn dual_fixing_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 5,
            num_row_: 3,
            col_cost_: vec![0.0, -8.0, 0.0, 0.0, 1.0],
            col_lower_: vec![-kHighsInf, -kHighsInf, 0.0, -kHighsInf, 0.0],
            col_upper_: vec![kHighsInf, 7.0, 9.0, 3.0, 1.0],
            row_lower_: vec![-kHighsInf, -5.0, 0.0],
            row_upper_: vec![8.0, -5.0, kHighsInf],
            sense_: ObjSense::kMinimize,
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 5;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 2, 4, 5, 7, 8];
        lp.a_matrix_.index_ = vec![0, 2, 1, 2, 1, 0, 1, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, -1.0, -1.0, -2.0, -1.0, -1.0];
        lp
    }

    fn solve(lp: HighsLp, solver: &str, crossover: &str, h16_on: bool) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", solver),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("run_crossover", crossover),
            HighsStatus::kOk
        );
        // presolve_rule_logging is not a registered option in this port
        highs.options_.presolve_rule_logging = true;
        assert_eq!(
            highs.set_option_value_bool("h16_dual_postsolve_rules", h16_on),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        let _ = highs.run();
        highs
    }

    fn dominated_col_calls(highs: &Highs) -> HighsInt {
        highs.presolve_log_.rule[kPresolveRuleDominatedCol as usize].call
    }

    /// HiGHS 1.16 #3342 (h16_dual_postsolve_rules): HighsOptions.cpp::mayRequireDualPostsolve /
    /// mayRequireBasisPostsolve truth table (1.16.0 lines 1160-1179).
    #[test]
    fn h16_dual_postsolve_may_require_truth_table() {
        let mut options = HighsOptions::default();
        let cases: [(&str, &str, bool, bool); 5] = [
            ("pdlp", "on", true, false),
            ("ipm", "off", true, false),
            ("ipm", "choose", true, true),
            ("ipm", "on", false, true),
            ("simplex", "off", false, true),
        ];
        for (solver, crossover, dual, basis) in cases {
            options.solver = solver.to_string();
            options.run_crossover = crossover.to_string();
            assert_eq!(
                h16_may_require_dual_postsolve(&options),
                dual,
                "{solver}/{crossover} dual"
            );
            assert_eq!(
                h16_may_require_basis_postsolve(&options),
                basis,
                "{solver}/{crossover} basis"
            );
        }
    }

    /// HiGHS 1.16 #3348 (h16_dual_postsolve_rules), issue-3342: with IPM and no crossover the
    /// dominated column rule is never called on the on path (1.15.1 uses it to fix y
    /// and z), and the postsolved duals are the dual feasible 1.16 ones.
    #[test]
    fn h16_dual_postsolve_issue_3342_ipm_no_crossover() {
        let off = solve(issue_3342_lp(), "ipm", "off", false);
        assert!(dominated_col_calls(&off) > 0);

        let on = solve(issue_3342_lp(), "ipm", "off", true);
        assert_eq!(dominated_col_calls(&on), 0);
        assert_eq!(on.get_model_status(), HighsModelStatus::kOptimal);
        assert!((on.get_info().objective_function_value + 11.0).abs() < 1e-6);
        assert_eq!(on.get_info().num_dual_infeasibilities, 0);
        let sol = on.get_solution();
        let x = &sol.col_value;
        for (v, e) in x.iter().zip([2.0, 0.5, 2.0, 1.0]) {
            assert!((v - e).abs() < 1e-6, "{x:?}");
        }
        // The optimum is degenerate (w and r1 both at their upper bound), so the duals
        // are not unique: check stationarity d = c - A^T y and the 1.16 sign
        // conditions instead (w, y, z at upper bound: d <= 0; x interior: d = 0;
        // r0 at lower bound: y0 >= 0; r1 at upper bound: y1 <= 0).
        let (d, y) = (&sol.col_dual, &sol.row_dual);
        let cost: [f64; 4] = [-6.0, 2.0, -1.0, 2.0];
        let aty: [f64; 4] = [y[1], 2.0 * y[0], 2.0 * y[0], -2.0 * y[0] - 2.0 * y[1]];
        for j in 0..4 {
            assert!((cost[j] - aty[j] - d[j]).abs() < 1e-6, "{d:?} {y:?}");
        }
        assert!(
            d[0] <= 1e-7 && d[1].abs() <= 1e-7 && d[2] <= 1e-7 && d[3] <= 1e-7,
            "{d:?}"
        );
        assert!(y[0] >= -1e-7 && y[1] <= 1e-7, "{y:?}");
    }

    /// HiGHS 1.16 #3348: simplex needs basis postsolve, so on the on path the dominated column
    /// rule is off as well; the optimum is unchanged.
    #[test]
    fn h16_dual_postsolve_issue_3342_simplex() {
        let on = solve(issue_3342_lp(), "simplex", "on", true);
        assert_eq!(dominated_col_calls(&on), 0);
        assert_eq!(on.get_model_status(), HighsModelStatus::kOptimal);
        assert!((on.get_info().objective_function_value + 11.0).abs() < 1e-6);
        assert_eq!(on.get_info().num_dual_infeasibilities, 0);
    }

    /// HiGHS 1.16 #3342 (h16_dual_postsolve_rules), test-dual-fixing-primal-dual-postsolve: with
    /// IPM and no crossover dual fixing is off, so the postsolved duals are feasible
    /// (z keeps a nonnegative reduced cost at its lower bound).
    #[test]
    fn h16_dual_postsolve_dual_fixing_ipm_no_crossover() {
        let on = solve(dual_fixing_lp(), "ipm", "off", true);
        assert_eq!(on.get_model_status(), HighsModelStatus::kOptimal);
        assert!((on.get_info().objective_function_value + 56.0).abs() < 1e-6);
        assert_eq!(on.get_info().num_dual_infeasibilities, 0);
        let sol = on.get_solution();
        let x = &sol.col_value;
        for (j, e) in [(1, 7.0), (2, 9.0), (3, 3.0), (4, 0.0)] {
            assert!((x[j] - e).abs() < 1e-6, "{x:?}");
        }
        assert!(sol.col_dual[4] >= -1e-7, "{:?}", sol.col_dual);
    }
}
