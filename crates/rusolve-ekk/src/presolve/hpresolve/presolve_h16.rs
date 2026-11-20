//! HiGHS source: presolve/HPresolve.cpp (1.16.0 lines 65-210, 367-535, 3837-3893, 3920-4256)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiGHS 1.16 presolve driver, behind `h16_presolve_rules` (HiGHS #3167 #3210
//! #3272 #3300 #3334 #3342 #3348 #3498 #3536). With the option
//! off none of this runs and the 1.15.1 `ok_set_input` / `presolve` / `run` are used.
//!
//! What 1.16 moved, relative to the 1.15.1 code in `setup.rs` / `presolve_loop.rs`:
//! * `setInput` only chooses the rules (`chooseRules`, which now owns `allow_rule_`) and
//!   sets up the analysis; the presolve data structures are built inside `presolve`
//!   (`okSetupPresolveDataStructures`, then `setupSubstitutionOpportunities`).
//! * `presolve` logs the original size, then `checkOriginalModelBounds` rounds integer
//!   bounds (`ceil(l - feastol)`, `floor(u + feastol)`) and checks trivial
//!   infeasibility; `initialRowAndColPresolve` no longer rounds through
//!   `changeColBounds`.
//! * The progress report is gated on `output_flag` and also fires every 10 s; the final
//!   report is skipped when presolve reduced the problem to empty; `removeSlacks` is gone.
//!
//! The rules 1.16 added live in extra `allow_rule_` slots (hconst.rs,
//! `kH16PresolveRule*`); `presolve_rule_off` bits use the 1.16 numbering and are mapped
//! with `h16_presolve_rule_from_116`.
//!
//! The LP initial sweep (`HPresolveInitialSweep`, #3167 #3210 #3442) is in
//! `initial_sweep_h16.rs`. Not ported: the presolve timer clocks; the try/catch and
//! out-of-memory statuses of `run`.

use super::*;
use crate::prelude::*;

/// HiGHS 1.16 (#3498): HighsModelUtils.cpp::utilPresolveRuleTypeToString (1.16.0 lines
/// 1519-1572), indexed by the 1.16 rule number.
fn h16_rule_name(rule: HighsInt) -> &'static str {
    const NAMES: [&str; 24] = [
        "Empty row",
        "Singleton row",
        "Redundant row",
        "Empty column",
        "Fixed column",
        "Dominated column",
        "Forcing row",
        "Forcing column",
        "Free column substitution",
        "Doubleton equation",
        "Dependent equations",
        "Dependent free columns",
        "Aggregator",
        "Parallel rows and columns",
        "Sparsify",
        "Probing",
        "Enumeration",
        "Dual fixing",
        "Zero cost singleton",
        "Column stuffing",
        "Dual-fixing probing",
        "Initial sweep",
        "Fourier-Motzkin",
        "Weakly dominated column",
    ];
    if (0..24).contains(&rule) {
        NAMES[rule as usize]
    } else {
        "????"
    }
}

/// HiGHS: HighsIO.cpp::highsTimeSecondToString (`std::to_string(double)` in debug
/// builds, `std::to_string(static_cast<int>(time))` under NDEBUG, plus "s").
fn h16_time_second_to_string(time: f64) -> String {
    if cfg!(debug_assertions) {
        format!("{:.6}s", time)
    } else {
        format!("{}s", time as HighsInt)
    }
}

/// Port `allow_rule_` slot of a 1.16 rule number (`h16_presolve_rule_from_116`).
fn h16_rule_slot(rule116: HighsInt) -> usize {
    h16_presolve_rule_from_116(rule116) as usize
}

impl HPresolve {
    /// HiGHS 1.16 (#3167 #3300): HPresolve.cpp::setInput (lines 65-97), from the point
    /// after `model` / `options` / `timer` are assigned (done by the caller,
    /// `ok_set_input`). No data structures are set up here any more.
    pub(super) fn h16_set_input(&mut self, presolve_reduction_limit: HighsInt) -> bool {
        // Set up the logic to allow presolve rules
        self.h16_choose_rules();
        // Set up profiling for presolve rules and logging for their effectiveness
        self.analysis_.h16_setup(&self.model, &self.options);

        if self.mipsolver.is_none() {
            self.primal_feastol = self.options.primal_feasibility_tolerance;
            ok_assign(
                &mut self.model.integrality_,
                self.model.num_col_,
                HighsVarType::kContinuous,
            );
        } else {
            self.primal_feastol = self.options.mip_feasibility_tolerance;
        }

        // Take value passed in as reduction limit, allowing different
        // values to be used for initial presolve, and after restart
        self.reduction_limit = if presolve_reduction_limit < 0 {
            kHighsSize_tInf
        } else {
            presolve_reduction_limit as usize
        };
        if self.options.presolve != kHighsOffString && self.reduction_limit < kHighsSize_tInf {
            crate::highs_log_dev!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "HPresolve::setInput reductionLimit = {}\n",
                self.reduction_limit as i32
            );
        }
        true
    }

    /// HiGHS 1.16 (#3167 #3334 #3342 #3498 #3536): HPresolve.cpp::chooseRules (lines
    /// 367-535). `allow_rule_` has `kH16PresolveRuleCount` slots (1.15.1 numbering plus
    /// the four 1.16 rules appended, hconst.rs); the loops run over the 1.16 numbering.
    pub fn h16_choose_rules(&mut self) {
        let silent = self.silent_log();
        let count = kH16PresolveRuleCount as usize;
        let num_rules_116: HighsInt = 24;

        // All rules can be switched off, although the most fundamental can
        // only be switched off in initial sweep
        let allowed_off_in_initial_sweep = |rule: HighsInt| -> bool {
            rule == kPresolveRuleEmptyRow
                || rule == kPresolveRuleSingletonRow
                || rule == kPresolveRuleRedundantRow
                || rule == kPresolveRuleEmptyCol
                || rule == kPresolveRuleFixedCol
        };

        let log_rule = |options: &HighsOptions, rule: HighsInt| {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                " {:>1} Rule {:2} (set bit {:2} = {:7}): {}\n",
                if allowed_off_in_initial_sweep(rule) {
                    "*"
                } else {
                    " "
                },
                rule,
                rule,
                1i32 << rule,
                h16_rule_name(rule)
            );
        };

        if !silent && self.options.log_dev_level != 0 {
            // State which rules can be off, and what bit to set
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Permitted suppression of presolve rules via presolve_rule_off option:\n"
            );
            for rule in kPresolveRuleMin..num_rules_116 {
                log_rule(&self.options, rule);
            }
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                " * Only in initial sweep\n"
            );
        }

        // setRulesForFullMode: all rules on in full mode
        let full_mode = |allow: &mut Vec<bool>| {
            allow.clear();
            allow.resize(count, true);
        };
        // setRulesForLightMode: all rules on, then switch off the rules not used in
        // presolve light
        let light_mode = |allow: &mut Vec<bool>| {
            full_mode(allow);
            for rule in [
                kPresolveRuleDominatedCol,
                kPresolveRuleDependentEquations,
                kPresolveRuleDependentFreeCols,
                kPresolveRuleAggregator,
                kPresolveRuleParallelRowsAndCols,
                kPresolveRuleSparsify,
                kPresolveRuleProbing,
                kPresolveRuleEnumeration,
                kPresolveRuleDualFixing,
                kH16PresolveRuleZeroCostSingleton,
                kPresolveRuleColStuffing,
                kH16PresolveRuleFourierMotzkin,
            ] {
                allow[rule as usize] = false;
            }
        };
        // setRulesForInitialSweepMode: all rules off, then switch on the rules used in
        // initial sweep
        let initial_sweep_mode = |allow: &mut Vec<bool>| {
            allow.clear();
            allow.resize(count, false);
            for rule in [
                kPresolveRuleEmptyRow,
                kPresolveRuleSingletonRow,
                kPresolveRuleRedundantRow,
                kPresolveRuleEmptyCol,
                kPresolveRuleFixedCol,
                kPresolveRuleInitialSweep,
            ] {
                allow[rule as usize] = true;
            }
        };

        let mut allow = std::mem::take(&mut self.analysis_.allow_rule_);
        // Set the pre-defined rules on/off according to the mode
        let mode = self.options.presolve_mode;
        if mode == kPresolveModeDefault {
            // By default all presolve rules are allowed, but subsets for
            // specific solvers (esp PDLP) can be set
            if use_pdlp(&self.options.solver) {
                light_mode(&mut allow);
            } else {
                full_mode(&mut allow);
            }
        } else if mode == kPresolveModeInitialSweep {
            initial_sweep_mode(&mut allow);
        } else if mode == kPresolveModeLight {
            light_mode(&mut allow);
        } else if mode == kPresolveModeFull {
            full_mode(&mut allow);
        } else {
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kError,
                "Unrecognised presolve mode {}\n",
                mode
            );
            debug_assert!(false);
            // C++ leaves allow_rule_ unassigned (release); keep every rule on.
            full_mode(&mut allow);
        }
        // Possibly switch off specific rules defined by the user via
        // presolve_rule_off
        if self.options.presolve_rule_off != 0 {
            // Transform presolve_rule_off into logical settings in
            // allow_rule_[*], commenting on the rules switched off
            if !silent {
                crate::highs_log_user!(
                    &self.options.log_options,
                    HighsLogType::kInfo,
                    "Presolve rules switched off:\n"
                );
            }
            // Keep track of whether to add the qualifier "* Only in initial
            // sweep" to the logging!
            let mut off_in_initial_sweep = false;
            for rule in kPresolveRuleMin..num_rules_116 {
                let bit: HighsInt = 1 << rule;
                if self.options.presolve_rule_off & bit == 0 {
                    continue;
                }
                off_in_initial_sweep = allowed_off_in_initial_sweep(rule) || off_in_initial_sweep;
                allow[h16_rule_slot(rule)] = false;
                if !silent {
                    log_rule(&self.options, rule);
                }
            }
            if !silent && off_in_initial_sweep {
                crate::highs_log_user!(
                    &self.options.log_options,
                    HighsLogType::kInfo,
                    " * Only in initial sweep\n"
                );
            }
        }

        // Finally switch off rules that are incompatible with dual or basis
        // postsolve
        self.h16_may_require_basis_postsolve_ = false;
        self.h16_may_require_dual_postsolve_ = false;

        let rule_off_for_postsolve = |options: &HighsOptions,
                                      allow: &mut Vec<bool>,
                                      postsolve_type: &str,
                                      rule: HighsInt| {
            let slot = rule as usize;
            // If the rule is already switched off (by user action etc) then
            // no action is needed
            if !allow[slot] {
                return;
            }
            allow[slot] = false;
            if !silent && options.log_dev_level != 0 {
                let rule116 = (0..num_rules_116)
                    .find(|&r| h16_presolve_rule_from_116(r) == rule)
                    .unwrap_or(kPresolveRuleIllegal);
                crate::highs_log_user!(
                    &options.log_options,
                    HighsLogType::kWarning,
                    "May require {} postsolve, so switched off \"{}\" rule\n",
                    postsolve_type,
                    h16_rule_name(rule116)
                );
            }
        };

        if self.mipsolver.is_none() {
            // Record whether one or both of dual and basis postsolve may be
            // required, and switch off rules that are correspondingly
            // incompatible
            if super::presolve_loop::h16_may_require_dual_postsolve(&self.options) {
                self.h16_may_require_dual_postsolve_ = true;
                // Cannot use weakly dominated column rule (see
                // test-weakly-dominated-column-dual-postsolve)
                rule_off_for_postsolve(
                    &self.options,
                    &mut allow,
                    "dual",
                    kH16PresolveRuleWeaklyDominatedCol,
                );
                // Cannot use kPresolveRuleDualFixing rule (see
                // test-dual-fixing-dual-postsolve)
                rule_off_for_postsolve(&self.options, &mut allow, "dual", kPresolveRuleDualFixing);
                // Cannot use kPresolveRuleZeroCostSingleton (see
                // hipo-ipx-pdlp-no-crossover)
                rule_off_for_postsolve(
                    &self.options,
                    &mut allow,
                    "dual",
                    kH16PresolveRuleZeroCostSingleton,
                );
            }
            if super::presolve_loop::h16_may_require_basis_postsolve(&self.options) {
                self.h16_may_require_basis_postsolve_ = true;
                // Cannot use sparsify or any parallel row rule using
                // EqualityRowAddition(s). Latter are avoided logically in
                // HPresolve::detectParallelRowsAndCols
                rule_off_for_postsolve(&self.options, &mut allow, "basis", kPresolveRuleSparsify);
            }
            if self.h16_may_require_dual_postsolve_ || self.h16_may_require_basis_postsolve_ {
                // Cannot use dominated column rule for either dual or basis
                // postsolve (see issue-3342)
                rule_off_for_postsolve(
                    &self.options,
                    &mut allow,
                    "dual or basis",
                    kPresolveRuleDominatedCol,
                );
            }
        }
        // If dominated column rule is off, weakly dominated column rule
        // isn't currently considered, because checks for weakly dominated
        // columns are currently only in HPresolve::detectDominatedCol, but
        // this makes it explicit
        if !allow[kPresolveRuleDominatedCol as usize] {
            allow[kH16PresolveRuleWeaklyDominatedCol as usize] = false;
        }
        self.analysis_.allow_rule_ = allow;
    }

    /// HiGHS 1.16: HPresolve.h::hasPresolveDataStructures (line 189)
    pub fn h16_has_presolve_data_structures(&self) -> bool {
        !self.col_deleted.is_empty()
    }

    /// HiGHS 1.16 (#3167): HPresolve.cpp::okSetupPresolveDataStructures (lines 127-200).
    /// The resize / CSC / change-data steps are the 1.15.1 `okSetInput` ones.
    pub fn h16_ok_setup_presolve_data_structures(&mut self) -> bool {
        // Lambda for resizing the implied bound data
        let nc = self.model.num_col_;
        let nr = self.model.num_row_;
        if !ok_resize(&mut self.col_lower_source, nc, -1) {
            return false;
        }
        if !ok_resize(&mut self.col_upper_source, nc, -1) {
            return false;
        }
        if !ok_resize(&mut self.impl_col_lower, nc, -kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.impl_col_upper, nc, kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.col_impl_source_by_row, nr, BTreeSet::new()) {
            return false;
        }
        if !ok_resize(&mut self.impl_row_dual_source_by_col, nc, BTreeSet::new()) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_lower, nr, -kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_upper, nr, kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.impl_row_dual_lower, nr, -kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.impl_row_dual_upper, nr, kHighsInf) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_upper_source, nr, -1) {
            return false;
        }
        if !ok_resize(&mut self.row_dual_lower_source, nr, -1) {
            return false;
        }
        for i in 0..nr as usize {
            if self.model.row_lower_[i] == -kHighsInf {
                self.row_dual_upper[i] = 0.0;
            }
            if self.model.row_upper_[i] == kHighsInf {
                self.row_dual_lower[i] = 0.0;
            }
        }

        // Lambda for setting up the matrix
        let value = std::mem::take(&mut self.model.a_matrix_.value_);
        let index = std::mem::take(&mut self.model.a_matrix_.index_);
        let start = std::mem::take(&mut self.model.a_matrix_.start_);
        let ok = if self.model.a_matrix_.is_rowwise() {
            // Does this even happen?
            debug_assert!(self.model.a_matrix_.is_colwise());
            self.ok_from_csr(&value, &index, &start)
        } else {
            self.ok_from_csc(&value, &index, &start)
        };
        self.model.a_matrix_.value_ = value;
        self.model.a_matrix_.index_ = index;
        self.model.a_matrix_.start_ = start;
        if !ok {
            return false;
        }

        // numDeletedCols and numDeletedRows are not cumulative through the
        // whole of presolve, but "since the last time the model had no
        // deleted columns or rows" - ie from initialisation here, or from a
        // call to shrinkProblem
        self.num_deleted_cols = 0;
        self.num_deleted_rows = 0;
        // Need to reset current number of deleted rows and columns in logging
        self.analysis_.reset_num_deleted();

        // Lambda for resizing the change tracking data
        //
        // initialize everything as changed, but do not add all indices
        // since the first thing presolve will do is a scan for easy reductions
        // of each row and column and set the flag of processed columns to false
        // from then on they are added to the vector whenever there are changes
        if !ok_resize(&mut self.changed_row_flag, nr, 1u8) {
            return false;
        }
        if !ok_resize(&mut self.row_deleted, nr, 0u8) {
            return false;
        }
        if !ok_reserve(&mut self.changed_row_indices, nr) {
            return false;
        }
        if !ok_resize(&mut self.changed_col_flag, nc, 1u8) {
            return false;
        }
        if !ok_resize(&mut self.col_deleted, nc, 0u8) {
            return false;
        }
        if !ok_reserve(&mut self.changed_col_indices, nc) {
            return false;
        }
        // okReserve(liftingOpportunities, ...): a BTreeMap has no capacity (setup.rs).
        ok_resize(&mut self.single_equation_checked, nr, 0u8)
    }

    /// HiGHS 1.16 (#3167): HPresolve.cpp::setupSubstitutionOpportunities (lines 202-210)
    pub fn h16_setup_substitution_opportunities(&mut self) {
        for row in 0..self.model.num_row_ {
            if !self.is_dual_implied_free(row) {
                continue;
            }
            let mut it = self.get_row_vector(row);
            while let Some(nonzero) = it.next_nz(self) {
                if self.is_implied_free(nonzero.index()) {
                    self.substitution_opportunities.push((row, nonzero.index()));
                }
            }
        }
    }

    /// HiGHS 1.16 (#3272 #3300 #3334): HPresolve.cpp::checkOriginalModelBounds (lines
    /// 3837-3862). Integer bounds are rounded on the model arrays before the presolve
    /// data structures exist; `ceil(0 - feastol)` is `-0.0`, as in C++.
    pub fn h16_check_original_model_bounds(&mut self) -> HPresolveResult {
        let is_mip = self.mipsolver.is_some();
        debug_assert!(!is_mip || !self.model.integrality_.is_empty());
        for i_col in 0..self.model.num_col_ as usize {
            if is_mip && self.model.integrality_[i_col] != HighsVarType::kContinuous {
                // Perform integer rounding of bounds on integer variables
                self.model.col_lower_[i_col] =
                    (self.model.col_lower_[i_col] - self.primal_feastol).ceil();
                self.model.col_upper_[i_col] =
                    (self.model.col_upper_[i_col] + self.primal_feastol).floor();
            }
            // Check for trivial primal infeasibility
            if self.model.col_lower_[i_col] > self.model.col_upper_[i_col] {
                return HPresolveResult::kPrimalInfeasible;
            }
        }
        // Check for trivial primal infeasibility in rows
        for i_row in 0..self.model.num_row_ as usize {
            if self.model.row_lower_[i_row] > self.model.row_upper_[i_row] {
                return HPresolveResult::kPrimalInfeasible;
            }
        }
        HPresolveResult::kOk
    }

    /// HiGHS 1.16 (#3167): HPresolve.cpp::initialRowAndColPresolve (lines 3864-3893).
    /// The integer rounding moved to `checkOriginalModelBounds`.
    pub fn h16_initial_row_and_col_presolve(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        // Lambda for a full scan over the rows
        for row in 0..self.model.num_row_ {
            if self.row_deleted[row as usize] != 0 {
                continue;
            }
            hpresolve_checked_call!(self.row_presolve(postsolve_stack, row));
            self.changed_row_flag[row as usize] = 0;
        }
        // Lambda for a full scan over the cols
        for col in 0..self.model.num_col_ {
            if self.col_deleted[col as usize] != 0 {
                continue;
            }
            hpresolve_checked_call!(self.col_presolve(postsolve_stack, col));
            self.changed_col_flag[col as usize] = 0;
        }
        self.check_limits(postsolve_stack)
    }

    /// HiGHS 1.16 (#3167 #3300 #3334): HPresolve.cpp::reducedToEmpty (HPresolve.h lines
    /// 184-187)
    pub fn h16_reduced_to_empty(&self) -> bool {
        self.num_deleted_cols == self.model.num_col_ && self.num_deleted_rows == self.model.num_row_
    }

    /// HiGHS 1.16 (#3167 #3210 #3300 #3334): HPresolve.cpp::presolve (lines 3920-4256)
    ///
    /// The rule-order loop is the 1.15.1 one (presolve_loop.rs) except where marked. The
    /// LP initial sweep block (lines 3994-4019) is not ported (module doc).
    pub fn h16_presolve(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HPresolveResult {
        // C++ lambda `presolveReturn`.
        fn presolve_return(
            this: &mut HPresolve,
            postsolve_stack: &mut HighsPostsolveStack,
        ) -> HPresolveResult {
            if this.mipsolver.is_some() {
                hpresolve_checked_call!(this.scale_mip(postsolve_stack));
            }
            // rusolve RUSOLVE_PRESOLVE_PROF (debug output only)
            if presolve_prof::on() {
                let submip = this.mipsolver.is_some() && this.mip().submip;
                presolve_prof::report(if submip { "submip" } else { "main" });
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

        // C++ lambda `report` (captures `silent` and `current_time` by reference).
        fn report(this: &HPresolve, silent: bool, current_time: f64) {
            if !silent {
                let num_col = this.model.num_col_ - this.num_deleted_cols;
                let num_row = this.model.num_row_ - this.num_deleted_rows;
                let num_nonz = (this.a_value.len() - this.freeslots.len()) as HighsInt;
                let mut time_str = h16_time_second_to_string(current_time);
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

        // shrinkProblem followed by toCSC / okFromCSC (as presolve_loop.rs).
        fn shrink_and_rebuild(this: &mut HPresolve, postsolve_stack: &mut HighsPostsolveStack) {
            this.shrink_problem(postsolve_stack);
            let mut a_value = std::mem::take(&mut this.model.a_matrix_.value_);
            let mut a_index = std::mem::take(&mut this.model.a_matrix_.index_);
            let mut a_start = std::mem::take(&mut this.model.a_matrix_.start_);
            this.to_csc(&mut a_value, &mut a_index, &mut a_start);
            this.ok_from_csc(&a_value, &a_index, &a_start);
            this.model.a_matrix_.value_ = a_value;
            this.model.a_matrix_.index_ = a_index;
            this.model.a_matrix_.start_ = a_start;
        }

        // convert model to minimization problem
        if self.model.sense_ == ObjSense::kMaximize {
            for i in 0..self.model.num_col_ as usize {
                self.model.col_cost_[i] = -self.model.col_cost_[i];
            }
            self.model.offset_ = -self.model.offset_;
            debug_assert!(self.model.offset_.is_finite());
            self.model.sense_ = ObjSense::kMinimize;
        }

        // Need to check for time-out in checkLimits, so make sure that
        // the timer is well defined, and that its total time clock is
        // running
        debug_assert!(self.timer.is_some());

        let silent = self.silent_log();
        let report_frequency = 10.0;
        let mut current_time = self.timer.as_ref().map_or(0.0, |t| t.read());
        let mut last_report_size: HighsInt = self.model.num_col_ + self.model.num_row_;
        let mut last_report_time = current_time;
        if self.options.presolve != kHighsOffString && !silent {
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Presolving model\n"
            );
            let mut time_str = h16_time_second_to_string(current_time);
            if self.options.timeless_log {
                time_str = String::new();
            }
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "{} rows, {} cols, {} nonzeros {}\n",
                self.model.num_row_,
                self.model.num_col_,
                self.model.a_matrix_.num_nz(),
                time_str
            );
        }

        // Perform integer rounding of bounds on integer variables and check
        // for trivial bound violations - which yield
        // Result::kPrimalInfeasible
        hpresolve_checked_call!(self.h16_check_original_model_bounds());

        // HPresolve.cpp lines 3994-4019: the LP initial sweep (initial_sweep_h16.rs)
        hpresolve_checked_call!(self.h16_initial_sweep(postsolve_stack));

        if !self.h16_ok_setup_presolve_data_structures() {
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kError,
                "Insufficient memory for presolve data structures\n"
            );
            // C++ Result::kOutOfMemory; the port's Result has no such status.
            return HPresolveResult::kStopped;
        }

        // initialize substitution opportunities
        self.h16_setup_substitution_opportunities();

        if self.options.presolve != kHighsOffString {
            // 1.15.1 line kept: clique table presolve flag.
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 6606 at 7b890c255
            // (removed: the flag is scoped to runProbing / enumerateSolutions)
            if self.mipsolver.is_some() && !self.options.mip_presolve_clique_table {
                self.mip_mut().mipdata_.cliquetable.set_presolve_flag(true);
            }

            debug_assert!(self.timer.is_some());
            debug_assert!(self.timer.as_ref().unwrap().running());

            // Possibly just perform a bespoke presolve rule test
            if self.options.presolve_rule_test != 0 {
                hpresolve_checked_call!(self.presolve_rule_test(postsolve_stack));
                return presolve_return(self, postsolve_stack);
            }

            {
                let _prof = presolve_prof::scope(presolve_prof::INITIAL);
                hpresolve_checked_call!(self.h16_initial_row_and_col_presolve(postsolve_stack));
            }
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
            let mut try_sparsify = !self.h16_may_require_basis_postsolve_;
            let mut try_probing = self.mipsolver.is_some();
            let mut num_cliques_before_probing: HighsInt = -1;
            let mut domcol_after_probing_called = false;
            let mut dependent_equations_called = self.mipsolver.is_some();
            // HiGHS #3255 (mip_fm_presolve): HPresolve.cpp line 4070 (1.16.0)
            let mut try_fourier_motzkin = self.mipsolver.is_some();

            // Start of main presolve loop
            //
            loop {
                debug_assert!(self.model.num_col_ >= self.num_deleted_cols);
                debug_assert!(self.model.num_row_ >= self.num_deleted_rows);
                let current_size = self.model.num_col_ - self.num_deleted_cols
                    + self.model.num_row_
                    - self.num_deleted_rows;
                if self.options.output_flag {
                    current_time = self.timer.as_ref().unwrap().read();
                    if (current_size as f64) < 0.85 * (last_report_size as f64)
                        || current_time > last_report_time + report_frequency
                    {
                        last_report_size = current_size;
                        last_report_time = current_time;
                        report(self, silent, current_time);
                    }
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

                // HiGHS #3255 (mip_fm_presolve): HPresolve.cpp lines 4106-4111 (1.16.0)
                let mut num_cols_eliminated_fourier_motzkin: HighsInt = 0;
                if try_fourier_motzkin
                    && self.x1b_fm_presolve_enabled()
                    && self.analysis_.allow_rule_[kH16PresolveRuleFourierMotzkin as usize]
                {
                    hpresolve_checked_call!(self.fourier_motzkin(
                        postsolve_stack,
                        &mut num_cols_eliminated_fourier_motzkin
                    ));
                    try_fourier_motzkin = false;
                }

                if self.h16_reduced_to_empty() {
                    break;
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

            if !self.h16_reduced_to_empty() {
                report(self, silent, current_time);
            }
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
