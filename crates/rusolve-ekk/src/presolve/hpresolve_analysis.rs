//! HiGHS source: presolve/HPresolveAnalysis.h, presolve/HPresolveAnalysis.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Per-rule logging of LP presolve reductions.
//!
//! # Design: no stored pointers
//!
//! C++ keeps `const HighsLp* model`, `const HighsOptions* options` and
//! `const HighsInt* numDeletedRows/numDeletedCols` (pointing into the owning
//! HPresolve) after `setup`. Rust stores none of them: each method takes
//! what it reads as arguments.
//! * `*numDeletedRows` -> parameter `num_deleted_rows`, `*numDeletedCols`
//!   -> `num_deleted_cols`;
//! * `model->x` -> `model.x`, `options->x` / `options_->x` -> `options.x`.
//!
//! HPresolve calls these through the wrappers
//! `HPresolve::analysis_start_presolve_rule_log(rule)` /
//! `HPresolve::analysis_stop_presolve_rule_log(rule)` (see
//! `presolve/hpresolve/mod.rs`), so HPresolve translators write
//! `analysis_.startPresolveRuleLog(r)` as
//! `self.analysis_start_presolve_rule_log(r)`.
//!
//! `allow_rule` (a `const bool*`) is declared in C++ but never set or read:
//! omitted. `presolveReductionTypeToString` is declared but never defined:
//! omitted.

use crate::prelude::*;

/// HiGHS: HPresolveAnalysis.h::HPresolveAnalysis
///
/// C++ leaves every scalar uninitialised until `setup` /
/// `resetNumDeleted`; Rust defaults them to 0 / false (unobservable: HPresolve
/// calls `setup` in `run` before any use, and `presolve` assigns
/// `allow_logging_` / `logging_on_`).
#[derive(Clone, Debug, Default)]
pub struct HPresolveAnalysis {
    original_num_col_: HighsInt,
    original_num_row_: HighsInt,

    pub allow_rule_: Vec<bool>,

    pub allow_logging_: bool,
    pub logging_on_: bool,

    pub log_rule_type_: i32,
    pub num_deleted_rows0_: HighsInt,
    pub num_deleted_cols0_: HighsInt,
    pub presolve_log_: HighsPresolveLog,
}

impl HPresolveAnalysis {
    /// HiGHS: HPresolveAnalysis.cpp::HPresolveAnalysis::setup (lines 11-74, 64 lines)
    ///
    /// The `numDeletedRows_` / `numDeletedCols_` reference parameters only
    /// set the (dropped) pointers, so they are not passed.
    pub fn setup(&mut self, model: &HighsLp, options: &HighsOptions, silent: bool) {
        // C++ 16-19 assign the (dropped) pointer members model/options/
        // numDeletedRows/numDeletedCols; they are arguments here (module doc).
        self.allow_rule_ = vec![true; kPresolveRuleCount as usize];

        if options.presolve_rule_off != 0 || options.log_dev_level != 0 {
            // Some presolve rules are off
            //
            // Transform options->presolve_rule_off into logical settings in
            // allow_rule_[*], commenting on the rules switched off
            if !silent {
                if options.presolve_rule_off != 0 {
                    crate::highs_log_user!(
                        &options.log_options,
                        HighsLogType::kInfo,
                        "Presolve rules not allowed:\n"
                    );
                } else {
                    crate::highs_log_user!(
                        &options.log_options,
                        HighsLogType::kInfo,
                        "Permitted suppression of presolve rules via presolve_rule_off option:\n"
                    );
                }
            }
            let mut bit: HighsInt = 1;
            for rule_type in kPresolveRuleMin..kPresolveRuleCount {
                // Identify whether this rule is allowed
                let allow = (options.presolve_rule_off & bit) == 0;
                if rule_type >= kPresolveRuleFirstAllowOff {
                    // This is a rule that can be switched off, so comment
                    // positively if it is off
                    self.allow_rule_[rule_type as usize] = allow;
                    if !silent {
                        if !allow || (options.presolve_rule_off == 0 && options.log_dev_level != 0)
                        {
                            crate::highs_log_user!(
                                &options.log_options,
                                HighsLogType::kInfo,
                                "   Rule {:2} (set bit {:2} = {:5}): {}\n",
                                rule_type,
                                rule_type,
                                bit,
                                util_presolve_rule_type_to_string(rule_type)
                            );
                        }
                    }
                } else if !allow && !silent {
                    // This is a rule that cannot be switched off so, if an
                    // attempt is made, don't allow it to be off and comment
                    // negatively
                    crate::highs_log_user!(
                        &options.log_options,
                        HighsLogType::kWarning,
                        "Cannot disallow rule {:2} (bit {:2} = {:5}): {}\n",
                        rule_type,
                        rule_type,
                        bit,
                        util_presolve_rule_type_to_string(rule_type)
                    );
                }
                bit *= 2;
            }
        }
        // Allow logging if option is set and model is not a MIP
        self.allow_logging_ = options.presolve_rule_logging && !model.is_mip();
        self.logging_on_ = self.allow_logging_;
        self.log_rule_type_ = kPresolveRuleIllegal;
        self.reset_num_deleted();
        self.presolve_log_.clear();
        self.original_num_col_ = model.num_col_;
        self.original_num_row_ = model.num_row_;
    }

    /// HiGHS 1.16 (#3167 #3300): HPresolveAnalysis.cpp::HPresolveAnalysis::setup (lines 11-60)
    ///
    /// The 1.16 setup no longer decodes `presolve_rule_off` (allow_rule_ moved to
    /// HPresolve::chooseRules, see `presolve_h16.rs`). The presolve timer clocks
    /// (`analyse_presolve_time_`) are not ported. Called only under h16_presolve_rules.
    pub fn h16_setup(&mut self, model: &HighsLp, options: &HighsOptions) {
        let lp_presolve = !model.is_mip() || options.solve_relaxation;
        // Allow logging if option is set and LP presolve is being used
        self.allow_logging_ = options.presolve_rule_logging && lp_presolve;
        self.logging_on_ = self.allow_logging_;
        self.log_rule_type_ = kPresolveRuleIllegal;
        self.reset_num_deleted();
        self.presolve_log_.clear();
        self.original_num_col_ = model.num_col_;
        self.original_num_row_ = model.num_row_;
    }

    /// HiGHS: HPresolveAnalysis.cpp::HPresolveAnalysis::resetNumDeleted (lines 85-88, 4 lines)
    pub fn reset_num_deleted(&mut self) {
        self.num_deleted_rows0_ = 0;
        self.num_deleted_cols0_ = 0;
    }

    /// HiGHS: HPresolveAnalysis.cpp::HPresolveAnalysis::startPresolveRuleLog (lines 90-133, 44 lines)
    ///
    /// `printf` -> `print!`; `assert` -> `debug_assert!`.
    pub fn start_presolve_rule_log(
        &mut self,
        rule_type: HighsInt,
        num_deleted_rows: HighsInt,
        num_deleted_cols: HighsInt,
        model: &HighsLp,
    ) {
        let debug_print = false;
        debug_assert!(self.logging_on_);
        debug_assert!(rule_type >= kPresolveRuleMin && rule_type <= kPresolveRuleMax);
        debug_assert!(self.allow_rule_[rule_type as usize]);
        // Prevent any future calls to "start" until logging is on again
        self.logging_on_ = false;
        let check_rule: HighsInt = kPresolveRuleIllegal;
        if debug_print {
            println!(
                "   startPresolveRuleLog [{:6}, {:6}] for ({:2}) {}",
                num_deleted_rows,
                num_deleted_cols,
                rule_type,
                util_presolve_rule_type_to_string(rule_type)
            );
        }
        if rule_type == check_rule {
            println!(
                ">> startPresolveRuleLog [{:6}, {:6}] for ({:2}) {}",
                check_rule,
                num_deleted_rows,
                num_deleted_cols,
                util_presolve_rule_type_to_string(check_rule)
            );
        }
        self.presolve_log_.rule[rule_type as usize].call += 1;
        // Check that stop has been called since the last start
        debug_assert!(self.log_rule_type_ == kPresolveRuleIllegal);
        self.log_rule_type_ = rule_type;
        // Check that no un-logged reductions have been performed
        if self.num_deleted_rows0_ != num_deleted_rows
            || self.num_deleted_cols0_ != num_deleted_cols
        {
            println!(
                "ERROR: Model {}: {} = num_deleted_rows0_ != *numDeletedRows = {} ||{} = num_deleted_cols0_ != *numDeletedCols = {}",
                model.model_name_,
                self.num_deleted_rows0_,
                num_deleted_rows,
                self.num_deleted_cols0_,
                num_deleted_cols
            );
            // fflush(stdout) is omitted (stdout flushing is not a solver effect).
        }
        debug_assert!(self.num_deleted_rows0_ == num_deleted_rows);
        debug_assert!(self.num_deleted_cols0_ == num_deleted_cols);
        self.num_deleted_rows0_ = num_deleted_rows;
        self.num_deleted_cols0_ = num_deleted_cols;
        let check_num_deleted_rows0_: HighsInt = -255;
        let check_num_deleted_cols0_: HighsInt = -688;
        if self.num_deleted_rows0_ == check_num_deleted_rows0_
            && self.num_deleted_cols0_ == check_num_deleted_cols0_
        {
            println!(
                "num_deleted ({}, {})",
                self.num_deleted_rows0_, self.num_deleted_cols0_
            );
        }
    }

    /// HiGHS: HPresolveAnalysis.cpp::HPresolveAnalysis::stopPresolveRuleLog (lines 135-176, 42 lines)
    pub fn stop_presolve_rule_log(
        &mut self,
        rule_type: HighsInt,
        num_deleted_rows: HighsInt,
        num_deleted_cols: HighsInt,
        model: &HighsLp,
    ) {
        let _ = model; // C++ does not read model in this method
        let debug_print = false;
        debug_assert!(self.logging_on_);
        debug_assert!(rule_type == self.log_rule_type_);
        if debug_print {
            println!(
                "    stopPresolveRuleLog [{:6}, {:6}] for ({:2}) {}",
                num_deleted_rows,
                num_deleted_cols,
                rule_type,
                util_presolve_rule_type_to_string(rule_type)
            );
        }
        let check_rule: HighsInt = kPresolveRuleIllegal;
        if rule_type == check_rule {
            println!(
                ">>  stopPresolveRuleLog [{:6}, {:6}] for ({:2}) {}",
                check_rule,
                num_deleted_rows,
                num_deleted_cols,
                util_presolve_rule_type_to_string(check_rule)
            );
        }
        let num_removed_row: HighsInt = num_deleted_rows - self.num_deleted_rows0_;
        let num_removed_col: HighsInt = num_deleted_cols - self.num_deleted_cols0_;
        debug_assert!(num_removed_row >= 0);
        debug_assert!(num_removed_col >= 0);
        self.presolve_log_.rule[rule_type as usize].col_removed += num_removed_col;
        self.presolve_log_.rule[rule_type as usize].row_removed += num_removed_row;

        // Set the rule type to be illegal to indicate that stop has been
        // called, and update the record of num_deleted_rows/cols
        self.log_rule_type_ = kPresolveRuleIllegal;
        self.num_deleted_rows0_ = num_deleted_rows;
        self.num_deleted_cols0_ = num_deleted_cols;

        let report = false;
        if report {
            println!(
                "{:<25} Call {:9}: ({:3}, {:3}) ({:3}, {:3})",
                util_presolve_rule_type_to_string(rule_type),
                self.presolve_log_.rule[rule_type as usize].call,
                num_removed_col,
                num_removed_row,
                self.presolve_log_.rule[rule_type as usize].col_removed,
                self.presolve_log_.rule[rule_type as usize].row_removed
            );
        }
        let check_num_deleted_rows0_: HighsInt = -212;
        let check_num_deleted_cols0_: HighsInt = -637;
        if self.num_deleted_rows0_ == check_num_deleted_rows0_
            && self.num_deleted_cols0_ == check_num_deleted_cols0_
        {
            println!(
                "num_deleted ({}, {})",
                self.num_deleted_rows0_, self.num_deleted_cols0_
            );
        }
    }

    /// HiGHS: HPresolveAnalysis.cpp::HPresolveAnalysis::analysePresolveRuleLog (lines 178-239, 62 lines)
    ///
    /// C++ default arguments: `report = false`.
    pub fn analyse_presolve_rule_log(
        &mut self,
        report: bool,
        model: &HighsLp,
        options: &HighsOptions,
        num_deleted_rows: HighsInt,
        num_deleted_cols: HighsInt,
    ) -> bool {
        if !self.allow_logging_ {
            return true;
        }
        let log_options = &options.log_options;
        let mut sum_removed_row: HighsInt = 0;
        let mut sum_removed_col: HighsInt = 0;
        for rule_type in kPresolveRuleMin..kPresolveRuleCount {
            sum_removed_row += self.presolve_log_.rule[rule_type as usize].row_removed;
            sum_removed_col += self.presolve_log_.rule[rule_type as usize].col_removed;
        }
        if report && sum_removed_row + sum_removed_col != 0 {
            let rule = "-------------------------------------------------------";
            crate::highs_log_user!(log_options, HighsLogType::kInfo, "{}\n", rule);
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "{:<25}      Rows      Cols     Calls\n",
                "Presolve rule removed"
            );
            crate::highs_log_user!(log_options, HighsLogType::kInfo, "{}\n", rule);
            for rule_type in kPresolveRuleMin..kPresolveRuleCount {
                if self.presolve_log_.rule[rule_type as usize].call != 0
                    || self.presolve_log_.rule[rule_type as usize].row_removed != 0
                    || self.presolve_log_.rule[rule_type as usize].col_removed != 0
                {
                    crate::highs_log_user!(
                        log_options,
                        HighsLogType::kInfo,
                        "{:<25} {:9} {:9} {:9}\n",
                        util_presolve_rule_type_to_string(rule_type),
                        self.presolve_log_.rule[rule_type as usize].row_removed,
                        self.presolve_log_.rule[rule_type as usize].col_removed,
                        self.presolve_log_.rule[rule_type as usize].call
                    );
                }
            }
            crate::highs_log_user!(log_options, HighsLogType::kInfo, "{}\n", rule);
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "{:<25} {:9} {:9}\n",
                "Total reductions",
                sum_removed_row,
                sum_removed_col
            );
            crate::highs_log_user!(log_options, HighsLogType::kInfo, "{}\n", rule);
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "{:<25} {:9} {:9}\n",
                "Original  model",
                self.original_num_row_,
                self.original_num_col_
            );
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "{:<25} {:9} {:9}\n",
                "Presolved model",
                self.original_num_row_ - sum_removed_row,
                self.original_num_col_ - sum_removed_col
            );
            crate::highs_log_user!(log_options, HighsLogType::kInfo, "{}\n", rule);
        }
        if self.original_num_row_ == model.num_row_ && self.original_num_col_ == model.num_col_ {
            if sum_removed_row != num_deleted_rows {
                crate::highs_log_dev!(
                    log_options,
                    HighsLogType::kError,
                    "{} = sum_removed_row != numDeletedRows = {}\n",
                    sum_removed_row,
                    num_deleted_rows
                );
                // fflush(stdout) is omitted (stdout flushing is not a solver effect).
                debug_assert!(sum_removed_row == num_deleted_rows);
                return false;
            }
            if sum_removed_col != num_deleted_cols {
                crate::highs_log_dev!(
                    log_options,
                    HighsLogType::kError,
                    "{} = sum_removed_col != numDeletedCols = {}\n",
                    sum_removed_col,
                    num_deleted_cols
                );
                // fflush(stdout) is omitted (stdout flushing is not a solver effect).
                debug_assert!(sum_removed_col == num_deleted_cols);
                return false;
            }
        }
        true
    }
}
