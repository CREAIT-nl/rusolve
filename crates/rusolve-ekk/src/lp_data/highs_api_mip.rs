//! HiGHS source: lp_data/Highs.cpp, lp_data/HighsOptions.cpp,
//! lp_data/HighsInterface.cpp, Highs.h (the `Highs` API the MIP solver calls
//! that the LP slice did not need)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port. Callers in `crate::mip`: `HighsLpRelaxation`
//! (`lpsolver`), the repair LP of `HighsMipSolverData::trySolution`
//! (`tmpSolver`), and the heuristics' local LPs.
//!
//! * `setOptionValue(name, value)` is typed by overload: `_bool`, `_int`,
//!   `_double`, `_str`. Only the option names `crate::mip` passes are
//!   recognised (each C++ call site's literal name). Port
//!   `setLocalOptionValue`: an out-of-range or invalid value logs and
//!   returns `kError` without changing the option. On success sync
//!   `options_.log_options` (`sync_log_options`) and return
//!   `optionChangeAction()` (`option_change_action`).
//! * `setProfiling` is a no-op (timers are stubbed). `writeModel` is only
//!   reached by debug branches that are compiled out (`const bool ... =
//!   false`), so it is not ported.
//! * The option records are plain `HighsOptions` fields here, so the record
//!   type and bounds of each recognised name are tabled below
//!   (`option_kind`, `int_bounds`, `double_bounds`).

use crate::lp_data::highs_option_records::{self as registry, OptionValue};
use crate::lp_data::highs_run::Highs;
use crate::prelude::*;

/// HiGHS: HighsOptions.h::HighsOptionType.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OptionKind {
    Bool,
    Int,
    Double,
    Str,
}

/// HiGHS: HighsOptions.cpp::getOptionIndex (the name lookup) and the record's type: the
/// registry (`highs_option_records`) is the only name table. `None` is C++
/// `OptionStatus::kUnknownOption`.
fn option_kind(name: &str) -> Option<OptionKind> {
    match registry::option_type_name(name) {
        Some("bool") => Some(OptionKind::Bool),
        Some("HighsInt") => Some(OptionKind::Int),
        Some("double") => Some(OptionKind::Double),
        Some("string") => Some(OptionKind::Str),
        _ => None,
    }
}

/// HiGHS: HighsOptions.cpp::optionEntryTypeToString (lines 47-57)
fn option_entry_type_to_string(kind: OptionKind) -> &'static str {
    match kind {
        OptionKind::Bool => "bool",
        OptionKind::Int => "HighsInt",
        OptionKind::Double => "double",
        OptionKind::Str => "string",
    }
}

/// HiGHS: HighsOptions.cpp::getOptionIndex (lines 225-235), unknown-name log.
fn log_unknown_option(log_options: &HighsLogOptions, name: &str) {
    crate::highs_log_user!(
        log_options,
        HighsLogType::kError,
        "getOptionIndex: Option \"{}\" is unknown\n",
        name
    );
}

/// HiGHS: HighsOptions.cpp::boolFromString (lines 210-223)
fn bool_from_string(value: &str) -> Option<bool> {
    match value.to_ascii_lowercase().as_str() {
        "t" | "true" | "1" | kHighsOnString => Some(true),
        "f" | "false" | "0" | kHighsOffString => Some(false),
        _ => None,
    }
}

/// C `strtol(s, _, 10)` on the leading `[+-]?digits` of `s`, saturating as strtol
/// does. Returns the value and the number of characters consumed (0 when there
/// are no digits, i.e. no conversion).
fn c_strtol_prefix(s: &str) -> (i64, usize) {
    let b = s.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let digits_start = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == digits_start {
        return (0, 0);
    }
    let value = match s[..i].parse::<i64>() {
        Ok(v) => v,
        Err(_) if b[0] == b'-' => i64::MIN,
        Err(_) => i64::MAX,
    };
    (value, i)
}

/// C `atoi(s)`: `(int) strtol(s, _, 10)`, so the saturated value is truncated to int.
fn c_atoi(s: &str) -> HighsInt {
    c_strtol_prefix(s).0 as HighsInt
}

/// C `atof(s)`: the longest prefix of `[+-]?digits[.digits][e[+-]digits]` (0.0 when
/// there is no digit). Callers have already restricted `s` to ASCII.
fn c_atof(s: &str) -> f64 {
    let b = s.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let mut digits = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
        digits += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return 0.0;
    }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let exp_start = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > exp_start {
            i = j;
        }
    }
    s[..i].parse::<f64>().unwrap_or(0.0)
}

/// HiGHS: HighsOptions.h::OptionRecordInt lower/upper bounds, read from the registry. The
/// fallback is the (0, kHighsIInf) that most int records declare, and the exceptions are the
/// records whose 1.16.0 bounds the registry does not carry because the option is not wired.
fn int_bounds(name: &str) -> (HighsInt, HighsInt) {
    if let Some(bounds) = registry::option_int_bounds(name) {
        return bounds;
    }
    match name {
        "highs_analysis_level" => (kHighsAnalysisLevelMin, kHighsAnalysisLevelMax),
        "log_dev_level" => (kHighsLogDevLevelMin, kHighsLogDevLevelMax),
        "simplex_scale_strategy" => (kSimplexScaleStrategyMin, kSimplexScaleStrategyMax),
        "simplex_strategy" => (kSimplexStrategyMin, kSimplexStrategyMax),
        "iis_strategy" => (kIisStrategyMin, kIisStrategyMax),
        "presolve_mode" => (kPresolveModeMin, kPresolveModeMax),
        "icrash_iterations" => (0, 200),
        "icrash_approx_iter" => (0, 100),
        "pdlp_scaling_mode" => (kPdlpScalingMin, kPdlpScalingMax),
        "pdlp_restart_strategy" => (kPdlpRestartStrategyMin, kPdlpRestartStrategyMax),
        "pdlp_cupdlpc_restart_method" => (0, 2),
        "pdlp_step_size_strategy" => (kPdlpStepSizeStrategyMin, kPdlpStepSizeStrategyMax),
        "user_objective_scale" | "user_bound_scale" => (-kHighsIInf, kHighsIInf),
        "hipo_block_size" => (0, kHighsIInf),
        _ => (0, kHighsIInf),
    }
}

/// HiGHS: HighsOptions.h::OptionRecordDouble lower/upper bounds, read from the registry, with
/// the same fallback rule as [`int_bounds`].
fn double_bounds(name: &str) -> (f64, f64) {
    if let Some(bounds) = registry::option_double_bounds(name) {
        return bounds;
    }
    match name {
        "dual_feasibility_tolerance" | "primal_feasibility_tolerance" => {
            (kMinimumKktTolerance, kHighsInf)
        }
        "objective_bound" => (-kHighsInf, kHighsInf),
        "icrash_starting_weight" => (1e-10, 1e50),
        "pdlp_optimality_tolerance" => (kMinimumKktTolerance, kHighsInf),
        _ => (0.0, kHighsInf),
    }
}

/// HiGHS: HighsOptions.cpp::checkOptionValue(OptionRecordInt&, HighsInt) (lines 423-441)
fn check_int_value(
    log_options: &HighsLogOptions,
    name: &str,
    value: HighsInt,
    lower: HighsInt,
    upper: HighsInt,
) -> bool {
    if value < lower {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "checkOptionValue: Value {} for option \"{}\" is below lower bound of {}\n",
            value,
            name,
            lower
        );
        return false;
    } else if value > upper {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "checkOptionValue: Value {} for option \"{}\" is above upper bound of {}\n",
            value,
            name,
            upper
        );
        return false;
    }
    true
}

/// HiGHS: HighsOptions.cpp::checkOptionValue(OptionRecordDouble&, double) (lines 443-459)
fn check_double_value(
    log_options: &HighsLogOptions,
    name: &str,
    value: f64,
    lower: f64,
    upper: f64,
) -> bool {
    if value < lower {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "checkOptionValue: Value {} for option \"{}\" is below lower bound of {}\n",
            value,
            name,
            lower
        );
        return false;
    } else if value > upper {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "checkOptionValue: Value {} for option \"{}\" is above upper bound of {}\n",
            value,
            name,
            upper
        );
        return false;
    }
    true
}

/// HiGHS: HighsOptions.cpp::optionOffChooseOnOk (lines 59-70)
fn option_off_choose_on_ok(log_options: &HighsLogOptions, name: &str, value: &str) -> bool {
    if value == kHighsOffString || value == kHighsChooseString || value == kHighsOnString {
        return true;
    }
    crate::highs_log_user!(
        log_options,
        HighsLogType::kError,
        "Value \"{}\" for {} option is not one of \"{}\", \"{}\" or \"{}\"\n",
        value,
        name,
        kHighsOffString,
        kHighsChooseString,
        kHighsOnString
    );
    false
}

/// HiGHS: HighsOptions.cpp::optionSolverOk (lines 82-111).
fn option_solver_ok(log_options: &HighsLogOptions, value: &str, hipo_available: bool) -> bool {
    if value == kHighsChooseString
        || value == kSimplexString
        || value == kIpmString
        || (value == kHipoString && hipo_available)
        || value == kIpxString
        || value == kPdlpString
        || value == "qpasm"
        || value == kHiPdlpString
    {
        return true;
    }
    if value == kHipoString && !hipo_available {
        // HighsExternalApi::logUnavailable<HighsExtras::hipo> (HighsExternalApi.h:178-194):
        // "%s%sThe following features are unavailable: %s\n" with the formatted
        // message, " " (len > 0) and the std::set-ordered missing names of
        // hipo = require<amd, blas, metis, rcm> (HighsExtrasExternalDeps.h:220),
        // all missing when HiPO is unavailable.
        crate::highs_log_user!(
            log_options,
            HighsLogType::kError,
            "{}{}The following features are unavailable: {}\n",
            format!(
                "The HiPO solver was requested via the \"{}\" option.",
                kSolverString
            ),
            " ",
            "amd, blas, metis, rcm"
        );
        return false;
    }
    let prefix = if hipo_available { "hipo\", \"" } else { "" };
    crate::highs_log_user!(
        log_options,
        HighsLogType::kWarning,
        "Value \"{}\" for LP solver option (\"solver\") is not one of {}\"{}\", \"{}\", \"{}\", \"{}\" or \"{}\"\n",
        value,
        prefix,
        kHighsChooseString,
        kSimplexString,
        kIpmString,
        kIpxString,
        kPdlpString
    );
    false
}

/// HiGHS 1.16.0: HighsOptions.cpp::optionMipLpSolverOk (lines 133-142). HiPO is only a
/// valid choice when it is available, exactly as in C++.
fn option_mip_lp_solver_ok(
    log_options: &HighsLogOptions,
    value: &str,
    hipo_available: bool,
) -> bool {
    if value == kHighsChooseString
        || value == kSimplexString
        || value == kIpmString
        || value == kIpxString
        || (hipo_available && value == kHipoString)
    {
        return true;
    }
    crate::highs_log_user!(
        log_options,
        HighsLogType::kError,
        "Value \"{}\" for {} option (\"{}\") is not one of {}\"{}\", \"{}\", \"{}\", \"{}\"\n",
        value,
        "MIP LP solver",
        kMipLpSolverString,
        if hipo_available { "\"hipo\", " } else { "" },
        kHighsChooseString,
        kSimplexString,
        kIpmString,
        kIpxString
    );
    false
}

/// HiGHS 1.16.0: HighsOptions.cpp::optionMipIpmSolverOk (lines 144-152)
fn option_mip_ipm_solver_ok(
    log_options: &HighsLogOptions,
    value: &str,
    hipo_available: bool,
) -> bool {
    if value == kHighsChooseString
        || value == kIpmString
        || value == kIpxString
        || (hipo_available && value == kHipoString)
    {
        return true;
    }
    crate::highs_log_user!(
        log_options,
        HighsLogType::kError,
        "Value \"{}\" for {} option (\"{}\") is not one of {}\"{}\", \"{}\", \"{}\"\n",
        value,
        "MIP IPM solver",
        kMipIpmSolverString,
        if hipo_available { "\"hipo\", " } else { "" },
        kHighsChooseString,
        kIpmString,
        kIpxString
    );
    false
}

/// HiGHS: HighsOptions.cpp::optionHipoParallelTypeOk (lines 170-181)
fn option_hipo_parallel_type_ok(log_options: &HighsLogOptions, value: &str) -> bool {
    if value == kHipoNodeString || value == kHipoTreeString || value == kHipoBothString {
        return true;
    }
    crate::highs_log_user!(
        log_options,
        HighsLogType::kError,
        "Value \"{}\" for {} option is not one of \"{}\", \"{}\" or \"{}\"\n",
        value,
        kHipoParallelString,
        kHipoTreeString,
        kHipoNodeString,
        kHipoBothString
    );
    false
}

/// HiGHS: HighsOptions.cpp::optionHipoSystemOk (lines 183-194)
fn option_hipo_system_ok(log_options: &HighsLogOptions, value: &str) -> bool {
    if value == kHipoNormalEqString || value == kHipoAugmentedString || value == kHighsChooseString
    {
        return true;
    }
    crate::highs_log_user!(
        log_options,
        HighsLogType::kError,
        "Value \"{}\" for {} option is not one of \"{}\", \"{}\" or \"{}\"\n",
        value,
        kHipoSystemString,
        kHipoNormalEqString,
        kHipoAugmentedString,
        kHighsChooseString
    );
    false
}

/// HiGHS: HighsOptions.cpp::optionHipoOrderingOk (lines 196-208)
fn option_hipo_ordering_ok(log_options: &HighsLogOptions, value: &str) -> bool {
    if value == kHipoAmdString
        || value == kHipoMetisString
        || value == kHipoRcmString
        || value == kHighsChooseString
    {
        return true;
    }
    crate::highs_log_user!(
        log_options,
        HighsLogType::kError,
        "Value \"{}\" for {} option is not one of \"{}\", \"{}\", \"{}\" or \"{}\"\n",
        value,
        kHipoOrderingString,
        kHipoAmdString,
        kHipoMetisString,
        kHipoRcmString,
        kHighsChooseString
    );
    false
}

impl Highs {
    /// HiGHS: Highs.cpp::Highs::setOptionValue(const std::string&, bool) (lines 140-145)
    pub fn set_option_value_bool(&mut self, option: &str, value: bool) -> HighsStatus {
        let Some(kind) = option_kind(option) else {
            log_unknown_option(&self.options_.log_options, option);
            return HighsStatus::kError;
        };
        if kind != OptionKind::Bool {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "setLocalOptionValue: Option \"{}\" cannot be assigned a bool\n",
                option
            );
            return HighsStatus::kError;
        }
        match option {
            "output_flag" => self.options_.output_flag = value,
            "log_to_console" => self.options_.log_to_console = value,
            "run_centring" => self.options_.run_centring = value,
            "icrash" => self.options_.icrash = value,
            "icrash_dualize" => self.options_.icrash_dualize = value,
            "icrash_exact" => self.options_.icrash_exact = value,
            "icrash_breakpoints" => self.options_.icrash_breakpoints = value,
            "qp_allow_hot_start" => self.options_.qp_allow_hot_start = value,
            "solve_relaxation" => self.options_.solve_relaxation = value,
            "use_warm_start" => self.options_.use_warm_start = value,
            "blend_multi_objectives" => self.options_.blend_multi_objectives = value,
            "mip_improving_solution_save" => self.options_.mip_improving_solution_save = value,
            "mip_improving_solution_report_sparse" => {
                self.options_.mip_improving_solution_report_sparse = value
            }
            "mip_sb_budget_pre_restart" => self.options_.mip_sb_budget_pre_restart = value,
            "mip_dfs_node_limit" => self.options_.mip_dfs_node_limit = value,
            "mip_reduced_stack_fast" => self.options_.mip_reduced_stack_fast = value,
            "mip_cdouble_act_contrib" => self.options_.mip_cdouble_act_contrib = value,
            "mip_keep_branch_reason" => self.options_.mip_keep_branch_reason = value,
            "mip_redundant_vbd" => self.options_.mip_redundant_vbd = value,
            "mip_path_same_bds" => self.options_.mip_path_same_bds = value,
            "mip_strengthen_small_contrib" => self.options_.mip_strengthen_small_contrib = value,
            "mip_lurking_bounds_robust" => self.options_.mip_lurking_bounds_robust = value,
            "mip_dual_subst_limits" => self.options_.mip_dual_subst_limits = value,
            "mip_ranged_dual_implied" => self.options_.mip_ranged_dual_implied = value,
            "mip_machine_sched_sepa" => self.options_.mip_machine_sched_sepa = value,
            "mip_dual_fix_probing" => self.options_.mip_dual_fix_probing = value,
            "mip_aggregate_implications" => self.options_.mip_aggregate_implications = value,
            "mip_clique_normalise" => self.options_.mip_clique_normalise = value,
            "mip_coef_change_impl_bounds" => self.options_.mip_coef_change_impl_bounds = value,
            "mip_dual_side_lost" => self.options_.mip_dual_side_lost = value,
            "mip_cdouble_proofs" => self.options_.mip_cdouble_proofs = value,
            "mip_fix_last_active" => self.options_.mip_fix_last_active = value,
            "mip_keep_implications" => self.options_.mip_keep_implications = value,
            "mip_presolve_clique_table" => self.options_.mip_presolve_clique_table = value,
            "mip_presolve_row_append" => self.options_.mip_presolve_row_append = value,
            "mip_fm_presolve" => self.options_.mip_fm_presolve = value,
            "mip_sym_hash64" => self.options_.mip_sym_hash64 = value,
            "h16_ac_centring" => self.options_.h16_ac_centring = value,
            "h16_ac_fix" => self.options_.h16_ac_fix = value,
            "h16_add_no_modify_on_error" => self.options_.h16_add_no_modify_on_error = value,
            "h16_agg_term_zero" => self.options_.h16_agg_term_zero = value,
            "h16_assess_small" => self.options_.h16_assess_small = value,
            "h16_big_act" => self.options_.h16_big_act = value,
            "h16_bound_chain" => self.options_.h16_bound_chain = value,
            "h16_bound_range_inf" => self.options_.h16_bound_range_inf = value,
            "h16_cdouble" => self.options_.h16_cdouble = value,
            "h16_clique_links" => self.options_.h16_clique_links = value,
            "h16_col_stuffing" => self.options_.h16_col_stuffing = value,
            "h16_cost_perturb" => self.options_.h16_cost_perturb = value,
            "h16_crossover_x" => self.options_.h16_crossover_x = value,
            "h16_cupdlp_time" => self.options_.h16_cupdlp_time = value,
            "h16_cutpool_age" => self.options_.h16_cutpool_age = value,
            "h16_dep_eq_time_cap" => self.options_.h16_dep_eq_time_cap = value,
            "h16_dominated_cols" => self.options_.h16_dominated_cols = value,
            "h16_dual_free_list" => self.options_.h16_dual_free_list = value,
            "h16_dual_postsolve_rules" => self.options_.h16_dual_postsolve_rules = value,
            "h16_dup_row_undo" => self.options_.h16_dup_row_undo = value,
            "h16_dupcol_warn" => self.options_.h16_dupcol_warn = value,
            "h16_enum_sig64" => self.options_.h16_enum_sig64 = value,
            "h16_eq_row_adds" => self.options_.h16_eq_row_adds = value,
            "h16_feasrelax_sense" => self.options_.h16_feasrelax_sense = value,
            "h16_fj" => self.options_.h16_fj = value,
            "h16_frontier_merge" => self.options_.h16_frontier_merge = value,
            "h16_hash_common" => self.options_.h16_hash_common = value,
            "h16_hfactor_pivot" => self.options_.h16_hfactor_pivot = value,
            "h16_hipdlp_max" => self.options_.h16_hipdlp_max = value,
            "h16_hipo_display" => self.options_.h16_hipo_display = value,
            "h16_hipo_empty_rows" => self.options_.h16_hipo_empty_rows = value,
            "h16_impl_eq" => self.options_.h16_impl_eq = value,
            "h16_implint_bounds" => self.options_.h16_implint_bounds = value,
            "h16_int_bound_round" => self.options_.h16_int_bound_round = value,
            "h16_int_bound_rounding" => self.options_.h16_int_bound_rounding = value,
            "h16_int_round_assert_tol" => self.options_.h16_int_round_assert_tol = value,
            "h16_ipm_gap_offset" => self.options_.h16_ipm_gap_offset = value,
            "h16_ipx_free_rows" => self.options_.h16_ipx_free_rows = value,
            "h16_ipx_scaled" => self.options_.h16_ipx_scaled = value,
            "h16_ipx_stall" => self.options_.h16_ipx_stall = value,
            "h16_kkt_effective_costs" => self.options_.h16_kkt_effective_costs = value,
            "h16_max_duals" => self.options_.h16_max_duals = value,
            "h16_mip_start" => self.options_.h16_mip_start = value,
            "h16_mip_workers" => self.options_.h16_mip_workers = value,
            "h16_ninfmin" => self.options_.h16_ninfmin = value,
            "h16_nla_debug_gate" => self.options_.h16_nla_debug_gate = value,
            "h16_obj_cdouble" => self.options_.h16_obj_cdouble = value,
            "h16_parallel_rows" => self.options_.h16_parallel_rows = value,
            "h16_postsolve_basis_check" => self.options_.h16_postsolve_basis_check = value,
            "h16_presolve_rules" => self.options_.h16_presolve_rules = value,
            "presolve_fm_obj_reformulation" => self.options_.presolve_fm_obj_reformulation = value,
            "h16_row_duals" => self.options_.h16_row_duals = value,
            "h16_fourier_motzkin" => self.options_.h16_fourier_motzkin = value,
            "h16_qp_bounded_start" => self.options_.h16_qp_bounded_start = value,
            "h16_qp_convexity" => self.options_.h16_qp_convexity = value,
            "h16_ray_simplex" => self.options_.h16_ray_simplex = value,
            "h16_rebuild_reason" => self.options_.h16_rebuild_reason = value,
            "h16_repair_lp" => self.options_.h16_repair_lp = value,
            "h16_row_presolve" => self.options_.h16_row_presolve = value,
            "h16_row_redundant" => self.options_.h16_row_redundant = value,
            "h16_scaling_report" => self.options_.h16_scaling_report = value,
            "h16_singleton_loops" => self.options_.h16_singleton_loops = value,
            "h16_subst_int" => self.options_.h16_subst_int = value,
            "h16_sym_refine" => self.options_.h16_sym_refine = value,
            "h16_unresolved_unbounded_or_infeasible" => {
                self.options_.h16_unresolved_unbounded_or_infeasible = value
            }
            "h16_user_scale_offset" => self.options_.h16_user_scale_offset = value,
            "h16_user_sol_obj" => self.options_.h16_user_sol_obj = value,
            "h16_user_solution_assess" => self.options_.h16_user_solution_assess = value,
            "h16_watch_link" => self.options_.h16_watch_link = value,
            "h16_zero_cost_singleton" => self.options_.h16_zero_cost_singleton = value,
            "h16_sub_solver_flags" => self.options_.h16_sub_solver_flags = value,
            "mip_conflict_resolution" => self.options_.mip_conflict_resolution = value,
            "mip_conflict_resolution_tree_only" => {
                self.options_.mip_conflict_resolution_tree_only = value
            }
            "mip_symmetry_dynamic" => self.options_.mip_symmetry_dynamic = value,
            "mip_cut_selection_hybrid" => self.options_.mip_cut_selection_hybrid = value,
            "mip_cutoff_bigm" => self.options_.mip_cutoff_bigm = value,
            "mip_node_bigm" => self.options_.mip_node_bigm = value,
            "mip_node_bigm_rows" => self.options_.mip_node_bigm_rows = value,
            "lp_fold" => self.options_.lp_fold = value,
            "lp_crossover_ipm_crash" => self.options_.lp_crossover_ipm_crash = value,
            "lp_crossover_batch_push" => self.options_.lp_crossover_batch_push = value,
            "lp_dual_devex_start" => self.options_.lp_dual_devex_start = value,
            "lp_postsolve_snap_basis" => self.options_.lp_postsolve_snap_basis = value,
            "lp_ipm_support_reduction" => self.options_.lp_ipm_support_reduction = value,
            "lp_concurrent" => self.options_.lp_concurrent = value,
            "lp_solver_auto" => self.options_.lp_solver_auto = value,
            "lp_postsolve_fix_tiny" => self.options_.lp_postsolve_fix_tiny = value,
            "lp_fold_basis" => self.options_.lp_fold_basis = value,
            "lp_crossover_ipx_refine" => self.options_.lp_crossover_ipx_refine = value,
            "lp_ipx_no_drop" => self.options_.lp_ipx_no_drop = value,
            "lp_ipx_soft_drop" => self.options_.lp_ipx_soft_drop = value,
            "presolve_dep_eq_fixed_budget" => self.options_.presolve_dep_eq_fixed_budget = value,
            "lp_dualize_auto" => self.options_.lp_dualize_auto = value,
            "ipm_hipo_fast_blas" => self.options_.ipm_hipo_fast_blas = value,
            "ipm_hipo_cache_perm" => self.options_.ipm_hipo_cache_perm = value,
            "ipm_hipo_par_ne" => self.options_.ipm_hipo_par_ne = value,
            "ipm_hipo_stall_detect" => self.options_.ipm_hipo_stall_detect = value,
            "ipm_hipo_amd_first" => self.options_.ipm_hipo_amd_first = value,
            "ipm_hipo_dense_cols" => self.options_.ipm_hipo_dense_cols = value,
            "ipm_hipo_dense_schur" => self.options_.ipm_hipo_dense_schur = value,
            "ipm_hipo_dualize" => self.options_.ipm_hipo_dualize = value,
            "lp_ipm_failure_policy" => self.options_.lp_ipm_failure_policy = value,
            "mip_heuristic_run_oneopt" => self.options_.mip_heuristic_run_oneopt = value,
            "mip_heuristic_run_crossover" => self.options_.mip_heuristic_run_crossover = value,
            "mip_heuristic_run_fpr" => self.options_.mip_heuristic_run_fpr = value,
            "mip_heuristic_run_local_mip" => self.options_.mip_heuristic_run_local_mip = value,
            "mip_strengthen_cover_guard" => self.options_.mip_strengthen_cover_guard = value,
            "mip_ls_gub_swap" => self.options_.mip_ls_gub_swap = value,
            "mip_ls_int_lp_repair" => self.options_.mip_ls_int_lp_repair = value,
            "mip_alns" => self.options_.mip_alns = value,
            "mip_presolve_eq_subst" => self.options_.mip_presolve_eq_subst = value,
            "mip_presolve_clique_varbound" => self.options_.mip_presolve_clique_varbound = value,
            "mip_presolve_xor_gf2" => self.options_.mip_presolve_xor_gf2 = value,
            "mip_orbitope_root_fix" => self.options_.mip_orbitope_root_fix = value,
            "mip_orbitope_full_fix" => self.options_.mip_orbitope_full_fix = value,
            "mip_shift_and_propagate" => self.options_.mip_shift_and_propagate = value,
            "mip_proximity" => self.options_.mip_proximity = value,
            "mip_improve_early" => self.options_.mip_improve_early = value,
            "mip_primal_stall" => self.options_.mip_primal_stall = value,
            "mip_guided_dive" => self.options_.mip_guided_dive = value,
            "mip_local_branching" => self.options_.mip_local_branching = value,
            "mip_cdcl_dive" => self.options_.mip_cdcl_dive = value,
            "mip_root_bound_stall" => self.options_.mip_root_bound_stall = value,
            "mip_cdcl_uip" => self.options_.mip_cdcl_uip = value,
            "mip_zero_obj_cdcl" => self.options_.mip_zero_obj_cdcl = value,
            "mip_decompose" => self.options_.mip_decompose = value,
            "mip_lattice" => self.options_.mip_lattice = value,
            "mip_activation_dff" => self.options_.mip_activation_dff = value,
            "mip_fj_native" => self.options_.mip_fj_native = value,
            "mip_first_inc" => self.options_.mip_first_inc = value,
            "mip_first_inc_root" => self.options_.mip_first_inc_root = value,
            "mip_lattice_kernel" => self.options_.mip_lattice_kernel = value,
            "mip_root_lp_concurrent" => self.options_.mip_root_lp_concurrent = value,
            "mip_presolve_probing_cap" => self.options_.mip_presolve_probing_cap = value,
            "mip_presolve_domcol_full" => self.options_.mip_presolve_domcol_full = value,
            "mip_redcost_proof_uncapped" => self.options_.mip_redcost_proof_uncapped = value,
            "mip_root_stall_scip" => self.options_.mip_root_stall_scip = value,
            "mip_heur_lp_quota" => self.options_.mip_heur_lp_quota = value,
            "mip_heuristic_run_feasibility_jump" => {
                self.options_.mip_heuristic_run_feasibility_jump = value
            }
            // Registered bool option without an arm above: the registry maps it to its field
            // so no registered option is silently ignored.
            _ => {
                if let Err(err) =
                    registry::set_option_field(&mut self.options_, option, OptionValue::Bool(value))
                {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "setLocalOptionValue: {}\n",
                        err
                    );
                    return HighsStatus::kError;
                }
            }
        }
        self.options_.sync_log_options();
        self.option_change_action()
    }

    /// HiGHS: Highs.cpp::Highs::setOptionValue(const std::string&, HighsInt) (lines 147-153)
    pub fn set_option_value_int(&mut self, option: &str, value: HighsInt) -> HighsStatus {
        let Some(kind) = option_kind(option) else {
            log_unknown_option(&self.options_.log_options, option);
            return HighsStatus::kError;
        };
        match kind {
            OptionKind::Int => {}
            // C++ interprets an integer as a double for a double option
            OptionKind::Double => return self.set_option_value_double(option, value as f64),
            _ => {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "setLocalOptionValue: Option \"{}\" cannot be assigned an int\n",
                    option
                );
                return HighsStatus::kError;
            }
        }
        let (lower, upper) = int_bounds(option);
        if !check_int_value(&self.options_.log_options, option, value, lower, upper) {
            return HighsStatus::kError;
        }
        match option {
            "highs_analysis_level" => self.options_.highs_analysis_level = value,
            "ipm_iteration_limit" => self.options_.ipm_iteration_limit = value,
            "log_dev_level" => self.options_.log_dev_level = value,
            "threads" => self.options_.threads = value,
            "random_seed" => self.options_.random_seed = value,
            "simplex_iteration_limit" => self.options_.simplex_iteration_limit = value,
            "simplex_scale_strategy" => self.options_.simplex_scale_strategy = value,
            "simplex_strategy" => self.options_.simplex_strategy = value,
            "qp_iteration_limit" => self.options_.qp_iteration_limit = value,
            "qp_nullspace_limit" => self.options_.qp_nullspace_limit = value,
            "iis_strategy" => self.options_.iis_strategy = value,
            "presolve_mode" => self.options_.presolve_mode = value,
            "icrash_iterations" => self.options_.icrash_iterations = value,
            "icrash_approx_iter" => self.options_.icrash_approx_iter = value,
            "pdlp_iteration_limit" => self.options_.pdlp_iteration_limit = value,
            "pdlp_scaling_mode" => self.options_.pdlp_scaling_mode = value,
            "pdlp_ruiz_iterations" => self.options_.pdlp_ruiz_iterations = value,
            "pdlp_restart_strategy" => self.options_.pdlp_restart_strategy = value,
            "pdlp_cupdlpc_restart_method" => self.options_.pdlp_cupdlpc_restart_method = value,
            "pdlp_step_size_strategy" => self.options_.pdlp_step_size_strategy = value,
            "user_objective_scale" => self.options_.user_objective_scale = value,
            "user_bound_scale" => self.options_.user_bound_scale = value,
            "mip_max_start_nodes" => self.options_.mip_max_start_nodes = value,
            "mip_max_nodes" => self.options_.mip_max_nodes = value,
            "mip_max_leaves" => self.options_.mip_max_leaves = value,
            "hipo_block_size" => self.options_.hipo_block_size = value,
            _ => {
                if let Err(err) =
                    registry::set_option_field(&mut self.options_, option, OptionValue::Int(value))
                {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "setLocalOptionValue: {}\n",
                        err
                    );
                    return HighsStatus::kError;
                }
            }
        }
        self.options_.sync_log_options();
        self.option_change_action()
    }

    /// HiGHS: HighsInterface.cpp::Highs::optionChangeAction (lines 3047-3051)
    fn option_change_action(&mut self) -> HighsStatus {
        if self.iis_.valid_ && self.options_.iis_strategy != self.iis_.strategy_ {
            self.iis_.clear();
        }
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::setOptionValue(const std::string&, double) (lines 155-161)
    pub fn set_option_value_double(&mut self, option: &str, value: f64) -> HighsStatus {
        let Some(kind) = option_kind(option) else {
            log_unknown_option(&self.options_.log_options, option);
            return HighsStatus::kError;
        };
        if kind != OptionKind::Double {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "setLocalOptionValue: Option \"{}\" cannot be assigned a double\n",
                option
            );
            return HighsStatus::kError;
        }
        let (lower, upper) = double_bounds(option);
        if !check_double_value(&self.options_.log_options, option, value, lower, upper) {
            return HighsStatus::kError;
        }
        match option {
            "dual_feasibility_tolerance" => self.options_.dual_feasibility_tolerance = value,
            "objective_bound" => self.options_.objective_bound = value,
            "primal_feasibility_tolerance" => self.options_.primal_feasibility_tolerance = value,
            "primal_simplex_bound_perturbation_multiplier" => {
                self.options_.primal_simplex_bound_perturbation_multiplier = value
            }
            "time_limit" => self.options_.time_limit = value,
            "qp_regularization_value" => self.options_.qp_regularization_value = value,
            "iis_time_limit" => self.options_.iis_time_limit = value,
            "icrash_starting_weight" => self.options_.icrash_starting_weight = value,
            "pdlp_optimality_tolerance" => self.options_.pdlp_optimality_tolerance = value,
            "mip_rel_gap" => self.options_.mip_rel_gap = value,
            "mip_abs_gap" => self.options_.mip_abs_gap = value,
            _ => {
                if let Err(err) = registry::set_option_field(
                    &mut self.options_,
                    option,
                    OptionValue::Double(value),
                ) {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "setLocalOptionValue: {}\n",
                        err
                    );
                    return HighsStatus::kError;
                }
            }
        }
        self.option_change_action()
    }

    /// HiGHS: HighsOptions.cpp::setLocalOptionValue(..., const std::string& value_passed)
    /// (lines 575-660): the text parse for bool, int and double options. Called by
    /// `set_option_value_str` for every non-string kind.
    fn set_option_value_text(
        &mut self,
        option: &str,
        kind: OptionKind,
        value: &str,
    ) -> HighsStatus {
        // Trim any leading and trailing spaces
        let value_trim = value.trim_matches(' ');
        match kind {
            OptionKind::Bool => {
                let Some(value_bool) = bool_from_string(value_trim) else {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "setLocalOptionValue: Value \"{}\" cannot be interpreted as a bool\n",
                        value_trim
                    );
                    return HighsStatus::kError;
                };
                self.set_option_value_bool(option, value_bool)
            }
            OptionKind::Int => {
                // Check that the string only contains legitimate characters
                if value_trim.bytes().any(|b| !b"+-0123456789eE".contains(&b)) {
                    return HighsStatus::kError;
                }
                // sscanf(value, "%d%n", &value_int, &scanned_num_char). HiGHS quirk: when
                // no digit is scanned, C++ reads uninitialised variables; here that is
                // treated as "not converted" with value 0.
                let (scanned_value, scanned_num_char) = c_strtol_prefix(value_trim);
                let value_int = scanned_value as HighsInt;
                // Check that all characters in value_trim are scanned by the format
                if scanned_num_char == 0 || scanned_num_char != value_trim.len() {
                    crate::highs_log_dev!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "setLocalOptionValue: Value = \"{}\" converts via sscanf as {} by scanning {} of {} characters\n",
                        value_trim,
                        value_int,
                        scanned_num_char,
                        value_trim.len()
                    );
                    return HighsStatus::kError;
                }
                self.set_option_value_int(option, value_int)
            }
            OptionKind::Double => {
                // Check that the string only contains legitimate characters - after
                // handling +/- inf. C++ lower-cases the whole value here (tolower).
                let lower = value_trim.to_ascii_lowercase();
                let value_double = if lower == "inf" || lower == "+inf" {
                    kHighsInf
                } else if lower == "-inf" {
                    -kHighsInf
                } else {
                    if lower.bytes().any(|b| !b"+-.0123456789eE".contains(&b)) {
                        return HighsStatus::kError;
                    }
                    let value_int = c_atoi(&lower);
                    let value_double = c_atof(&lower);
                    let value_int_double = f64::from(value_int);
                    if value_double == value_int_double {
                        crate::highs_log_dev!(
                            &self.options_.log_options,
                            HighsLogType::kInfo,
                            "setLocalOptionValue: Value = \"{}\" converts via atoi as {} so is {} as double, and {} via atof\n",
                            lower,
                            value_int,
                            value_int_double,
                            value_double
                        );
                    }
                    value_double
                };
                self.set_option_value_double(option, value_double)
            }
            // Callers route string options to `set_option_value_str`'s own branch.
            OptionKind::Str => HighsStatus::kError,
        }
    }

    /// HiGHS: Highs.cpp::Highs::setOptionValue(const std::string&, const std::string&) (lines 163-170)
    pub fn set_option_value_str(&mut self, option: &str, value: &str) -> HighsStatus {
        let Some(kind) = option_kind(option) else {
            log_unknown_option(&self.options_.log_options, option);
            return HighsStatus::kError;
        };
        if kind != OptionKind::Str {
            return self.set_option_value_text(option, kind, value);
        }
        // Trim leading and trailing spaces; lower-case values except file names
        let mut value_trim = value.trim_matches(' ').to_string();
        // possibleLowerCaseOptionValue (HighsOptions.cpp:722-738) exempts file names
        if option != "log_file"
            && option != "mip_improving_solution_file"
            && option != "mip_debug_solution_file"
        {
            value_trim = value_trim.to_ascii_lowercase();
        }
        let log_options = &self.options_.log_options;
        let ok = match option {
            // C++: `!optionOffChooseOnOk(...) && value != "mip"`; the log is
            // written before the "mip" test, as in C++.
            "presolve" => {
                option_off_choose_on_ok(log_options, option, &value_trim) || value_trim == "mip"
            }
            "solver" => option_solver_ok(log_options, &value_trim, self.options_.hipo_available),
            "parallel" | "run_crossover" => {
                option_off_choose_on_ok(log_options, option, &value_trim)
            }
            "mip_lp_solver" => {
                option_mip_lp_solver_ok(log_options, &value_trim, self.options_.hipo_available)
            }
            "mip_ipm_solver" => {
                option_mip_ipm_solver_ok(log_options, &value_trim, self.options_.hipo_available)
            }
            "hipo_parallel_type" => option_hipo_parallel_type_ok(log_options, &value_trim),
            "hipo_system" => option_hipo_system_ok(log_options, &value_trim),
            "hipo_ordering" => option_hipo_ordering_ok(log_options, &value_trim),
            // log_file: no value check (C++ opens the file; not ported)
            // icrash_strategy: no value check in C++ (HighsOptions.h:1033-1036)
            _ => true,
        };
        if !ok {
            return HighsStatus::kError;
        }
        match option {
            "log_file" => self.options_.log_file = value_trim,
            "parallel" => self.options_.parallel = value_trim,
            "presolve" => self.options_.presolve = value_trim,
            "run_crossover" => self.options_.run_crossover = value_trim,
            "solver" => self.options_.solver = value_trim,
            "icrash_strategy" => self.options_.icrash_strategy = value_trim,
            "mip_improving_solution_file" => self.options_.mip_improving_solution_file = value_trim,
            "mip_debug_solution_file" => self.options_.mip_debug_solution_file = value_trim,
            "mip_lp_solver" => self.options_.mip_lp_solver = value_trim,
            "mip_ipm_solver" => self.options_.mip_ipm_solver = value_trim,
            "hipo_parallel_type" => self.options_.hipo_parallel_type = value_trim,
            "hipo_system" => self.options_.hipo_system = value_trim,
            "hipo_ordering" => self.options_.hipo_ordering = value_trim,
            _ => {}
        }
        self.option_change_action()
    }

    /// HiGHS: Highs.cpp::Highs::setCallback(HighsCallbackFunctionType, void*) (lines 2627-2637).
    /// The C++ `void* user_callback_data` is not a parameter: the closure owns its state.
    pub fn set_callback(&mut self, user_callback: HighsCallbackFunctionType) -> HighsStatus {
        let num_col = self.model_lp_.num_col_ as usize;
        self.callback_.clear(num_col);
        self.callback_.user_callback = Some(user_callback);

        self.options_.log_options.user_callback = self.callback_.user_callback.clone();
        self.options_.log_options.user_callback_active = false;
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::startCallback(const int callback_type) (lines 2660-2676)
    pub fn start_callback(&mut self, callback_type: HighsInt) -> HighsStatus {
        let callback_type_ok = callback_type >= kCallbackMin && callback_type <= kCallbackMax;
        debug_assert!(callback_type_ok);
        if !callback_type_ok {
            return HighsStatus::kError;
        }
        if self.callback_.user_callback.is_none() {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Cannot start callback when user_callback not defined\n"
            );
            return HighsStatus::kError;
        }
        debug_assert!(self.callback_.active.len() == kNumCallbackType as usize);
        self.callback_.active[callback_type as usize] = true;
        // Possibly modify the logging callback activity
        if callback_type == kCallbackLogging {
            self.options_.log_options.user_callback_active = true;
        }
        HighsStatus::kOk
    }

    /// HiGHS: Highs.cpp::Highs::stopCallback(const int callback_type) (lines 2696-2712).
    pub fn stop_callback(&mut self, callback_type: HighsInt) -> HighsStatus {
        let callback_type_ok = callback_type >= kCallbackMin && callback_type <= kCallbackMax;
        debug_assert!(callback_type_ok);
        if !callback_type_ok {
            return HighsStatus::kError;
        }
        if self.callback_.user_callback.is_none() {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "Cannot stop callback when user_callback not defined\n"
            );
            return HighsStatus::kWarning;
        }
        debug_assert!(self.callback_.active.len() == kNumCallbackType as usize);
        self.callback_.active[callback_type as usize] = false;
        // Possibly modify the logging callback activity
        if callback_type == kCallbackLogging {
            self.options_.log_options.user_callback_active = false;
        }
        HighsStatus::kOk
    }

    /// HiGHS: Highs.h::Highs::getOptionValue (bool overload; `output_flag`)
    pub fn get_option_value_bool(&self, option: &str) -> Option<bool> {
        let Some(kind) = option_kind(option) else {
            log_unknown_option(&self.options_.log_options, option);
            return None;
        };
        if kind != OptionKind::Bool {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "getLocalOptionValue: Option \"{}\" requires value of type {}, not bool\n",
                option,
                option_entry_type_to_string(kind)
            );
            return None;
        }
        match option {
            "log_to_console" => Some(self.options_.log_to_console),
            "output_flag" => Some(self.options_.output_flag),
            "run_centring" => Some(self.options_.run_centring),
            "icrash" => Some(self.options_.icrash),
            "icrash_dualize" => Some(self.options_.icrash_dualize),
            "icrash_exact" => Some(self.options_.icrash_exact),
            "icrash_breakpoints" => Some(self.options_.icrash_breakpoints),
            "qp_allow_hot_start" => Some(self.options_.qp_allow_hot_start),
            "solve_relaxation" => Some(self.options_.solve_relaxation),
            "use_warm_start" => Some(self.options_.use_warm_start),
            "blend_multi_objectives" => Some(self.options_.blend_multi_objectives),
            "mip_improving_solution_save" => Some(self.options_.mip_improving_solution_save),
            "mip_improving_solution_report_sparse" => {
                Some(self.options_.mip_improving_solution_report_sparse)
            }
            "mip_sb_budget_pre_restart" => Some(self.options_.mip_sb_budget_pre_restart),
            "mip_dfs_node_limit" => Some(self.options_.mip_dfs_node_limit),
            "mip_reduced_stack_fast" => Some(self.options_.mip_reduced_stack_fast),
            "mip_cdouble_act_contrib" => Some(self.options_.mip_cdouble_act_contrib),
            "mip_keep_branch_reason" => Some(self.options_.mip_keep_branch_reason),
            "mip_redundant_vbd" => Some(self.options_.mip_redundant_vbd),
            "mip_path_same_bds" => Some(self.options_.mip_path_same_bds),
            "mip_strengthen_small_contrib" => Some(self.options_.mip_strengthen_small_contrib),
            "mip_lurking_bounds_robust" => Some(self.options_.mip_lurking_bounds_robust),
            "mip_dual_subst_limits" => Some(self.options_.mip_dual_subst_limits),
            "mip_ranged_dual_implied" => Some(self.options_.mip_ranged_dual_implied),
            "mip_machine_sched_sepa" => Some(self.options_.mip_machine_sched_sepa),
            "mip_dual_fix_probing" => Some(self.options_.mip_dual_fix_probing),
            "mip_aggregate_implications" => Some(self.options_.mip_aggregate_implications),
            "mip_clique_normalise" => Some(self.options_.mip_clique_normalise),
            "mip_coef_change_impl_bounds" => Some(self.options_.mip_coef_change_impl_bounds),
            "mip_dual_side_lost" => Some(self.options_.mip_dual_side_lost),
            "mip_cdouble_proofs" => Some(self.options_.mip_cdouble_proofs),
            "mip_fix_last_active" => Some(self.options_.mip_fix_last_active),
            "mip_keep_implications" => Some(self.options_.mip_keep_implications),
            "mip_presolve_clique_table" => Some(self.options_.mip_presolve_clique_table),
            "mip_presolve_row_append" => Some(self.options_.mip_presolve_row_append),
            "mip_fm_presolve" => Some(self.options_.mip_fm_presolve),
            "mip_sym_hash64" => Some(self.options_.mip_sym_hash64),
            "h16_ac_centring" => Some(self.options_.h16_ac_centring),
            "h16_ac_fix" => Some(self.options_.h16_ac_fix),
            "h16_add_no_modify_on_error" => Some(self.options_.h16_add_no_modify_on_error),
            "h16_agg_term_zero" => Some(self.options_.h16_agg_term_zero),
            "h16_assess_small" => Some(self.options_.h16_assess_small),
            "h16_big_act" => Some(self.options_.h16_big_act),
            "h16_bound_chain" => Some(self.options_.h16_bound_chain),
            "h16_bound_range_inf" => Some(self.options_.h16_bound_range_inf),
            "h16_cdouble" => Some(self.options_.h16_cdouble),
            "h16_clique_links" => Some(self.options_.h16_clique_links),
            "h16_col_stuffing" => Some(self.options_.h16_col_stuffing),
            "h16_cost_perturb" => Some(self.options_.h16_cost_perturb),
            "h16_crossover_x" => Some(self.options_.h16_crossover_x),
            "h16_cupdlp_time" => Some(self.options_.h16_cupdlp_time),
            "h16_cutpool_age" => Some(self.options_.h16_cutpool_age),
            "h16_dep_eq_time_cap" => Some(self.options_.h16_dep_eq_time_cap),
            "h16_dominated_cols" => Some(self.options_.h16_dominated_cols),
            "h16_dual_free_list" => Some(self.options_.h16_dual_free_list),
            "h16_dual_postsolve_rules" => Some(self.options_.h16_dual_postsolve_rules),
            "h16_dup_row_undo" => Some(self.options_.h16_dup_row_undo),
            "h16_dupcol_warn" => Some(self.options_.h16_dupcol_warn),
            "h16_enum_sig64" => Some(self.options_.h16_enum_sig64),
            "h16_eq_row_adds" => Some(self.options_.h16_eq_row_adds),
            "h16_feasrelax_sense" => Some(self.options_.h16_feasrelax_sense),
            "h16_fj" => Some(self.options_.h16_fj),
            "h16_frontier_merge" => Some(self.options_.h16_frontier_merge),
            "h16_hash_common" => Some(self.options_.h16_hash_common),
            "h16_hfactor_pivot" => Some(self.options_.h16_hfactor_pivot),
            "h16_hipdlp_max" => Some(self.options_.h16_hipdlp_max),
            "h16_hipo_display" => Some(self.options_.h16_hipo_display),
            "h16_hipo_empty_rows" => Some(self.options_.h16_hipo_empty_rows),
            "h16_impl_eq" => Some(self.options_.h16_impl_eq),
            "h16_implint_bounds" => Some(self.options_.h16_implint_bounds),
            "h16_int_bound_round" => Some(self.options_.h16_int_bound_round),
            "h16_int_bound_rounding" => Some(self.options_.h16_int_bound_rounding),
            "h16_int_round_assert_tol" => Some(self.options_.h16_int_round_assert_tol),
            "h16_ipm_gap_offset" => Some(self.options_.h16_ipm_gap_offset),
            "h16_ipx_free_rows" => Some(self.options_.h16_ipx_free_rows),
            "h16_ipx_scaled" => Some(self.options_.h16_ipx_scaled),
            "h16_ipx_stall" => Some(self.options_.h16_ipx_stall),
            "h16_kkt_effective_costs" => Some(self.options_.h16_kkt_effective_costs),
            "h16_max_duals" => Some(self.options_.h16_max_duals),
            "h16_mip_start" => Some(self.options_.h16_mip_start),
            "h16_mip_workers" => Some(self.options_.h16_mip_workers),
            "h16_ninfmin" => Some(self.options_.h16_ninfmin),
            "h16_nla_debug_gate" => Some(self.options_.h16_nla_debug_gate),
            "h16_obj_cdouble" => Some(self.options_.h16_obj_cdouble),
            "h16_parallel_rows" => Some(self.options_.h16_parallel_rows),
            "h16_postsolve_basis_check" => Some(self.options_.h16_postsolve_basis_check),
            "h16_presolve_rules" => Some(self.options_.h16_presolve_rules),
            "presolve_fm_obj_reformulation" => Some(self.options_.presolve_fm_obj_reformulation),
            "h16_row_duals" => Some(self.options_.h16_row_duals),
            "h16_fourier_motzkin" => Some(self.options_.h16_fourier_motzkin),
            "h16_qp_bounded_start" => Some(self.options_.h16_qp_bounded_start),
            "h16_qp_convexity" => Some(self.options_.h16_qp_convexity),
            "h16_ray_simplex" => Some(self.options_.h16_ray_simplex),
            "h16_rebuild_reason" => Some(self.options_.h16_rebuild_reason),
            "h16_repair_lp" => Some(self.options_.h16_repair_lp),
            "h16_row_presolve" => Some(self.options_.h16_row_presolve),
            "h16_row_redundant" => Some(self.options_.h16_row_redundant),
            "h16_scaling_report" => Some(self.options_.h16_scaling_report),
            "h16_singleton_loops" => Some(self.options_.h16_singleton_loops),
            "h16_subst_int" => Some(self.options_.h16_subst_int),
            "h16_sym_refine" => Some(self.options_.h16_sym_refine),
            "h16_unresolved_unbounded_or_infeasible" => {
                Some(self.options_.h16_unresolved_unbounded_or_infeasible)
            }
            "h16_user_scale_offset" => Some(self.options_.h16_user_scale_offset),
            "h16_user_sol_obj" => Some(self.options_.h16_user_sol_obj),
            "h16_user_solution_assess" => Some(self.options_.h16_user_solution_assess),
            "h16_watch_link" => Some(self.options_.h16_watch_link),
            "h16_zero_cost_singleton" => Some(self.options_.h16_zero_cost_singleton),
            "h16_sub_solver_flags" => Some(self.options_.h16_sub_solver_flags),
            "mip_conflict_resolution" => Some(self.options_.mip_conflict_resolution),
            "mip_conflict_resolution_tree_only" => {
                Some(self.options_.mip_conflict_resolution_tree_only)
            }
            "mip_symmetry_dynamic" => Some(self.options_.mip_symmetry_dynamic),
            "mip_cut_selection_hybrid" => Some(self.options_.mip_cut_selection_hybrid),
            "mip_cutoff_bigm" => Some(self.options_.mip_cutoff_bigm),
            "mip_node_bigm" => Some(self.options_.mip_node_bigm),
            "mip_node_bigm_rows" => Some(self.options_.mip_node_bigm_rows),
            "lp_fold" => Some(self.options_.lp_fold),
            "lp_crossover_ipm_crash" => Some(self.options_.lp_crossover_ipm_crash),
            "lp_crossover_batch_push" => Some(self.options_.lp_crossover_batch_push),
            "lp_dual_devex_start" => Some(self.options_.lp_dual_devex_start),
            "lp_postsolve_snap_basis" => Some(self.options_.lp_postsolve_snap_basis),
            "lp_ipm_support_reduction" => Some(self.options_.lp_ipm_support_reduction),
            "lp_concurrent" => Some(self.options_.lp_concurrent),
            "lp_solver_auto" => Some(self.options_.lp_solver_auto),
            "lp_postsolve_fix_tiny" => Some(self.options_.lp_postsolve_fix_tiny),
            "lp_fold_basis" => Some(self.options_.lp_fold_basis),
            "lp_crossover_ipx_refine" => Some(self.options_.lp_crossover_ipx_refine),
            "lp_ipx_no_drop" => Some(self.options_.lp_ipx_no_drop),
            "lp_ipx_soft_drop" => Some(self.options_.lp_ipx_soft_drop),
            "presolve_dep_eq_fixed_budget" => Some(self.options_.presolve_dep_eq_fixed_budget),
            "lp_dualize_auto" => Some(self.options_.lp_dualize_auto),
            "ipm_hipo_fast_blas" => Some(self.options_.ipm_hipo_fast_blas),
            "ipm_hipo_cache_perm" => Some(self.options_.ipm_hipo_cache_perm),
            "ipm_hipo_par_ne" => Some(self.options_.ipm_hipo_par_ne),
            "ipm_hipo_stall_detect" => Some(self.options_.ipm_hipo_stall_detect),
            "ipm_hipo_amd_first" => Some(self.options_.ipm_hipo_amd_first),
            "ipm_hipo_dense_cols" => Some(self.options_.ipm_hipo_dense_cols),
            "ipm_hipo_dense_schur" => Some(self.options_.ipm_hipo_dense_schur),
            "ipm_hipo_dualize" => Some(self.options_.ipm_hipo_dualize),
            "lp_ipm_failure_policy" => Some(self.options_.lp_ipm_failure_policy),
            "mip_heuristic_run_oneopt" => Some(self.options_.mip_heuristic_run_oneopt),
            "mip_heuristic_run_crossover" => Some(self.options_.mip_heuristic_run_crossover),
            "mip_heuristic_run_fpr" => Some(self.options_.mip_heuristic_run_fpr),
            "mip_heuristic_run_local_mip" => Some(self.options_.mip_heuristic_run_local_mip),
            "mip_strengthen_cover_guard" => Some(self.options_.mip_strengthen_cover_guard),
            "mip_ls_gub_swap" => Some(self.options_.mip_ls_gub_swap),
            "mip_ls_int_lp_repair" => Some(self.options_.mip_ls_int_lp_repair),
            "mip_alns" => Some(self.options_.mip_alns),
            "mip_presolve_eq_subst" => Some(self.options_.mip_presolve_eq_subst),
            "mip_presolve_clique_varbound" => Some(self.options_.mip_presolve_clique_varbound),
            "mip_presolve_xor_gf2" => Some(self.options_.mip_presolve_xor_gf2),
            "mip_orbitope_root_fix" => Some(self.options_.mip_orbitope_root_fix),
            "mip_orbitope_full_fix" => Some(self.options_.mip_orbitope_full_fix),
            "mip_shift_and_propagate" => Some(self.options_.mip_shift_and_propagate),
            "mip_proximity" => Some(self.options_.mip_proximity),
            "mip_improve_early" => Some(self.options_.mip_improve_early),
            "mip_primal_stall" => Some(self.options_.mip_primal_stall),
            "mip_guided_dive" => Some(self.options_.mip_guided_dive),
            "mip_local_branching" => Some(self.options_.mip_local_branching),
            "mip_cdcl_dive" => Some(self.options_.mip_cdcl_dive),
            "mip_root_bound_stall" => Some(self.options_.mip_root_bound_stall),
            "mip_cdcl_uip" => Some(self.options_.mip_cdcl_uip),
            "mip_zero_obj_cdcl" => Some(self.options_.mip_zero_obj_cdcl),
            "mip_decompose" => Some(self.options_.mip_decompose),
            "mip_lattice" => Some(self.options_.mip_lattice),
            "mip_activation_dff" => Some(self.options_.mip_activation_dff),
            "mip_fj_native" => Some(self.options_.mip_fj_native),
            "mip_first_inc" => Some(self.options_.mip_first_inc),
            "mip_first_inc_root" => Some(self.options_.mip_first_inc_root),
            "mip_lattice_kernel" => Some(self.options_.mip_lattice_kernel),
            "mip_root_lp_concurrent" => Some(self.options_.mip_root_lp_concurrent),
            "mip_presolve_probing_cap" => Some(self.options_.mip_presolve_probing_cap),
            "mip_presolve_domcol_full" => Some(self.options_.mip_presolve_domcol_full),
            "mip_redcost_proof_uncapped" => Some(self.options_.mip_redcost_proof_uncapped),
            "mip_root_stall_scip" => Some(self.options_.mip_root_stall_scip),
            "mip_heur_lp_quota" => Some(self.options_.mip_heur_lp_quota),
            "mip_heuristic_run_feasibility_jump" => {
                Some(self.options_.mip_heuristic_run_feasibility_jump)
            }
            _ => match registry::get_option_field(&self.options_, option) {
                Some(OptionValue::Bool(value)) => Some(value),
                _ => None,
            },
        }
    }

    /// HiGHS: Highs.h::Highs::getOptionValue (HighsInt overload) -> HighsOptions.cpp::getLocalOptionValue
    ///
    pub fn get_option_value_int(&self, option: &str) -> Option<HighsInt> {
        let Some(kind) = option_kind(option) else {
            log_unknown_option(&self.options_.log_options, option);
            return None;
        };
        if kind != OptionKind::Int {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "getLocalOptionValue: Option \"{}\" requires value of type {}, not HighsInt\n",
                option,
                option_entry_type_to_string(kind)
            );
            return None;
        }
        match option {
            "highs_analysis_level" => Some(self.options_.highs_analysis_level),
            "ipm_iteration_limit" => Some(self.options_.ipm_iteration_limit),
            "log_dev_level" => Some(self.options_.log_dev_level),
            "threads" => Some(self.options_.threads),
            "random_seed" => Some(self.options_.random_seed),
            "simplex_iteration_limit" => Some(self.options_.simplex_iteration_limit),
            "simplex_scale_strategy" => Some(self.options_.simplex_scale_strategy),
            "simplex_strategy" => Some(self.options_.simplex_strategy),
            "qp_iteration_limit" => Some(self.options_.qp_iteration_limit),
            "qp_nullspace_limit" => Some(self.options_.qp_nullspace_limit),
            "iis_strategy" => Some(self.options_.iis_strategy),
            "presolve_mode" => Some(self.options_.presolve_mode),
            "icrash_iterations" => Some(self.options_.icrash_iterations),
            "icrash_approx_iter" => Some(self.options_.icrash_approx_iter),
            "pdlp_iteration_limit" => Some(self.options_.pdlp_iteration_limit),
            "pdlp_scaling_mode" => Some(self.options_.pdlp_scaling_mode),
            "pdlp_ruiz_iterations" => Some(self.options_.pdlp_ruiz_iterations),
            "pdlp_restart_strategy" => Some(self.options_.pdlp_restart_strategy),
            "pdlp_cupdlpc_restart_method" => Some(self.options_.pdlp_cupdlpc_restart_method),
            "pdlp_step_size_strategy" => Some(self.options_.pdlp_step_size_strategy),
            "user_objective_scale" => Some(self.options_.user_objective_scale),
            "user_bound_scale" => Some(self.options_.user_bound_scale),
            "mip_max_start_nodes" => Some(self.options_.mip_max_start_nodes),
            "mip_max_nodes" => Some(self.options_.mip_max_nodes),
            "mip_max_leaves" => Some(self.options_.mip_max_leaves),
            "hipo_block_size" => Some(self.options_.hipo_block_size),
            _ => match registry::get_option_field(&self.options_, option) {
                Some(OptionValue::Int(value)) => Some(value),
                _ => None,
            },
        }
    }

    /// HiGHS: Highs.h::Highs::getOptionValue (double overload). The port had no double
    /// getter; the registry supplies it for every registered double option.
    pub fn get_option_value_double(&self, option: &str) -> Option<f64> {
        let Some(kind) = option_kind(option) else {
            log_unknown_option(&self.options_.log_options, option);
            return None;
        };
        if kind != OptionKind::Double {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "getLocalOptionValue: Option \"{}\" requires value of type {}, not double\n",
                option,
                option_entry_type_to_string(kind)
            );
            return None;
        }
        match registry::get_option_field(&self.options_, option) {
            Some(OptionValue::Double(value)) => Some(value),
            _ => None,
        }
    }

    /// HiGHS: Highs.h::Highs::getOptionValue (string overload; `presolve`, `solver`)
    pub fn get_option_value_str(&self, option: &str) -> Option<String> {
        let Some(kind) = option_kind(option) else {
            log_unknown_option(&self.options_.log_options, option);
            return None;
        };
        if kind != OptionKind::Str {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "getLocalOptionValue: Option \"{}\" requires value of type {}, not string\n",
                option,
                option_entry_type_to_string(kind)
            );
            return None;
        }
        match option {
            "log_file" => Some(self.options_.log_file.clone()),
            "parallel" => Some(self.options_.parallel.clone()),
            "presolve" => Some(self.options_.presolve.clone()),
            "run_crossover" => Some(self.options_.run_crossover.clone()),
            "solver" => Some(self.options_.solver.clone()),
            "icrash_strategy" => Some(self.options_.icrash_strategy.clone()),
            "mip_improving_solution_file" => {
                Some(self.options_.mip_improving_solution_file.clone())
            }
            "mip_debug_solution_file" => Some(self.options_.mip_debug_solution_file.clone()),
            "hipo_system" => Some(self.options_.hipo_system.clone()),
            "hipo_parallel_type" => Some(self.options_.hipo_parallel_type.clone()),
            "hipo_ordering" => Some(self.options_.hipo_ordering.clone()),
            _ => match registry::get_option_field(&self.options_, option) {
                Some(OptionValue::Str(value)) => Some(value),
                _ => None,
            },
        }
    }

    /// HiGHS: Highs.cpp::Highs::passOptions (lines 198-203) -> HighsOptions.cpp::passLocalOptions (lines 740-833)
    ///
    /// Every field of `options` is copied into this object. This object keeps
    /// its own `log_options`, which are re-synchronised as C++ does through the
    /// shared pointers. The per-record bound checks of `passLocalOptions` are
    /// not repeated, because the `HighsOptions` fields here carry no bounds
    /// table; a source object built with valid defaults passes them in C++.
    pub fn pass_options(&mut self, options: &HighsOptions) -> HighsStatus {
        let log_options = self.options_.log_options.clone();
        self.options_ = options.clone();
        self.options_.log_options = log_options;
        self.options_.sync_log_options();
        self.option_change_action()
    }

    /// HiGHS: Highs.h::Highs::getOptions (line 314)
    pub fn get_options(&self) -> &HighsOptions {
        &self.options_
    }

    /// HiGHS: Highs.cpp::Highs::changeColsCost (mask overload, lines 3092-3104)
    ///
    /// Reached from the feasibility pump (HighsPrimalHeuristics.cpp:1681).
    /// C++ asserts that `create` succeeds; so does this port.
    pub fn change_cols_cost_mask(&mut self, mask: &[HighsInt], cost: &[f64]) -> HighsStatus {
        self.clear_derived_model_properties();
        let mut index_collection = HighsIndexCollection::default();
        let create_error = create_mask(&mut index_collection, mask, self.model_lp_.num_col_);
        debug_assert!(create_error == 0);
        let _ = create_error;
        let call_status = self.change_costs_interface(&index_collection, cost);
        let mut return_status = HighsStatus::kOk;
        return_status = crate::lp_data::highs_status::interpret_call_status(
            &self.options_.log_options,
            call_status,
            return_status,
            "changeCosts",
        );
        if return_status == HighsStatus::kError {
            return HighsStatus::kError;
        }
        self.return_from_highs(return_status)
    }

    /// HiGHS: HighsInterface.cpp::Highs::changeCostsInterface (lines 908-938, 31 lines)
    ///
    /// Dropped: `doubleUserDataNotNull` (slices are never NULL). The final
    /// `assert(lp.has_infinite_cost_ == lp.hasInfiniteCost(..))` is a debug_assert.
    pub(crate) fn change_costs_interface(
        &mut self,
        index_collection: &HighsIndexCollection,
        cost: &[f64],
    ) -> HighsStatus {
        let num_cost: HighsInt = data_size(index_collection);
        // If a non-positive number of costs (may) need changing nothing needs to be
        // done
        if num_cost <= 0 {
            return HighsStatus::kOk;
        }
        // Take a copy of the cost that can be normalised
        let mut local_col_cost: Vec<f64> = cost[..num_cost as usize].to_vec();
        let mut return_status = HighsStatus::kOk;
        let mut local_has_infinite_cost = false;
        return_status = crate::lp_data::highs_status::interpret_call_status(
            &self.options_.log_options,
            crate::lp_data::lp_assess::assess_costs(
                &self.options_,
                0,
                index_collection,
                &mut local_col_cost,
                &mut local_has_infinite_cost,
                self.options_.infinite_cost,
            ),
            return_status,
            "assessCosts",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        let infinite_cost = self.options_.infinite_cost;
        let lp: &mut HighsLp = &mut self.model_lp_;
        crate::lp_data::highs_lp_modify::change_lp_costs(
            lp,
            index_collection,
            &local_col_cost,
            infinite_cost,
        );

        // Interpret possible introduction of infinite costs
        lp.has_infinite_cost_ = lp.has_infinite_cost_ || local_has_infinite_cost;
        debug_assert!(lp.has_infinite_cost_ == lp.has_infinite_cost(infinite_cost));

        // Deduce the consequences of new costs
        self.invalidate_model_status_solution_and_info();
        // Determine any implications for simplex data
        self.ekk_instance_.update_status(LpAction::kNewCosts);
        HighsStatus::kOk
    }

    /// HiGHS: Highs.h::Highs::getDualEdgeWeights (lines 1400-1404); `None` = C++ `nullptr`.
    pub fn get_dual_edge_weights(&self) -> Option<&[f64]> {
        if self.ekk_instance_.status_.has_dual_steepest_edge_weights {
            Some(self.ekk_instance_.dual_edge_weight_.as_slice())
        } else {
            None
        }
    }

    /// HiGHS: Highs.h::Highs::getNumNz (line 768)
    pub fn get_num_nz(&self) -> HighsInt {
        self.model_lp_.a_matrix_.num_nz()
    }

    /// HiGHS: Highs.h::Highs::setProfiling (line 1362). No-op: timers are stubbed.
    pub fn set_profiling(&mut self) {}

    /// HiGHS: Highs.h::Highs::modelStatusToString
    pub fn model_status_to_string(&self, model_status: HighsModelStatus) -> String {
        util_model_status_to_string(model_status).to_string()
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

    #[test]
    fn s4_opt_default_off() {
        assert!(!crate::lp_data::highs_options::HighsOptions::default().mip_cut_selection_hybrid);
        let highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_cut_selection_hybrid"),
            Some(false)
        );
    }

    /// rusolve: `mip_heuristic_run_feasibility_jump` (HighsOptions.h:1229) is settable by
    /// name, default true, through the bool and string setters.
    #[test]
    fn s5_data_fj_option_set_get() {
        assert_eq!(
            option_kind("mip_heuristic_run_feasibility_jump"),
            Some(OptionKind::Bool)
        );
        let mut highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_heuristic_run_feasibility_jump"),
            Some(true)
        );
        assert_eq!(
            highs.set_option_value_str("mip_heuristic_run_feasibility_jump", "false"),
            HighsStatus::kOk
        );
        assert!(!highs.options_.mip_heuristic_run_feasibility_jump);
        assert_eq!(
            highs.set_option_value_bool("mip_heuristic_run_feasibility_jump", true),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_heuristic_run_feasibility_jump"),
            Some(true)
        );
    }

    #[test]
    fn s4_opt_set_get() {
        assert_eq!(
            option_kind("mip_cut_selection_hybrid"),
            Some(OptionKind::Bool)
        );
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("mip_cut_selection_hybrid", true),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_cut_selection_hybrid"),
            Some(true)
        );
        assert_eq!(
            highs.set_option_value_bool("mip_cut_selection_hybrid", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_cut_selection_hybrid"),
            Some(false)
        );
        assert_eq!(
            highs.set_option_value_str("mip_cut_selection_hybrid", "true"),
            HighsStatus::kOk
        );
    }

    /// rusolve: `mip_strengthen_cover_guard` is off in
    /// the 1.15.1 defaults and on in the recommended profile (via `RECOMMENDED_MIP_ON`, not
    /// `set_all_h16`: it is a rusolve guard, not a HiGHS 1.16.0 option), and settable and
    /// gettable by name.
    #[test]
    fn mip_strengthen_cover_guard_default_off() {
        use crate::lp_data::highs_options::RECOMMENDED_MIP_ON;
        assert!(!crate::lp_data::highs_options::HighsOptions::default().mip_strengthen_cover_guard);
        assert!(
            crate::lp_data::highs_options::HighsOptions::recommended().mip_strengthen_cover_guard
        );
        assert!(RECOMMENDED_MIP_ON.contains(&"mip_strengthen_cover_guard"));
        let highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_strengthen_cover_guard"),
            Some(false)
        );
    }

    #[test]
    fn mip_strengthen_cover_guard_set_get() {
        assert_eq!(
            option_kind("mip_strengthen_cover_guard"),
            Some(OptionKind::Bool)
        );
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("mip_strengthen_cover_guard", true),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_strengthen_cover_guard"),
            Some(true)
        );
        assert_eq!(
            highs.set_option_value_bool("mip_strengthen_cover_guard", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_strengthen_cover_guard"),
            Some(false)
        );
        assert_eq!(
            highs.set_option_value_str("mip_strengthen_cover_guard", "true"),
            HighsStatus::kOk
        );
        assert!(highs.options_.mip_strengthen_cover_guard);
    }

    #[test]
    fn s6_opt_default_off() {
        let defaults = crate::lp_data::highs_options::HighsOptions::default();
        assert!(!defaults.mip_heuristic_run_oneopt);
        assert!(!defaults.mip_heuristic_run_crossover);
        let highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_heuristic_run_oneopt"),
            Some(false)
        );
        assert_eq!(
            highs.get_option_value_bool("mip_heuristic_run_crossover"),
            Some(false)
        );
    }

    #[test]
    fn s5_opt_default_off() {
        assert!(!crate::lp_data::highs_options::HighsOptions::default().mip_heuristic_run_fpr);
        assert!(
            !crate::lp_data::highs_options::HighsOptions::default().mip_heuristic_run_local_mip
        );
        let highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_heuristic_run_fpr"),
            Some(false)
        );
        assert_eq!(
            highs.get_option_value_bool("mip_heuristic_run_local_mip"),
            Some(false)
        );
    }

    #[test]
    fn ls_gub_swap_opt_default_off_and_settable() {
        assert!(!crate::lp_data::highs_options::HighsOptions::default().mip_ls_gub_swap);
        assert_eq!(option_kind("mip_ls_gub_swap"), Some(OptionKind::Bool));
        let mut highs = Highs::default();
        assert_eq!(highs.get_option_value_bool("mip_ls_gub_swap"), Some(false));
        assert_eq!(
            highs.set_option_value_bool("mip_ls_gub_swap", true),
            HighsStatus::kOk
        );
        assert_eq!(highs.get_option_value_bool("mip_ls_gub_swap"), Some(true));
        assert_eq!(
            highs.set_option_value_str("mip_ls_gub_swap", "false"),
            HighsStatus::kOk
        );
        assert_eq!(highs.get_option_value_bool("mip_ls_gub_swap"), Some(false));
    }

    #[test]
    fn s6_opt_set_get() {
        for name in ["mip_heuristic_run_oneopt", "mip_heuristic_run_crossover"] {
            assert_eq!(option_kind(name), Some(OptionKind::Bool));
            let mut highs = Highs::default();
            assert_eq!(highs.set_option_value_bool(name, true), HighsStatus::kOk);
            assert_eq!(highs.get_option_value_bool(name), Some(true));
            assert_eq!(highs.set_option_value_bool(name, false), HighsStatus::kOk);
            assert_eq!(highs.get_option_value_bool(name), Some(false));
            assert_eq!(highs.set_option_value_str(name, "true"), HighsStatus::kOk);
        }
    }

    #[test]
    fn alns_opt_set_get() {
        assert_eq!(option_kind("mip_alns"), Some(OptionKind::Bool));
        let mut highs = Highs::default();
        assert_eq!(highs.get_option_value_bool("mip_alns"), Some(false));
        assert_eq!(
            highs.set_option_value_bool("mip_alns", true),
            HighsStatus::kOk
        );
        assert_eq!(highs.get_option_value_bool("mip_alns"), Some(true));
        assert_eq!(
            highs.set_option_value_str("mip_alns", "false"),
            HighsStatus::kOk
        );
        assert_eq!(highs.get_option_value_bool("mip_alns"), Some(false));
    }

    #[test]
    fn s5_opt_set_get() {
        assert_eq!(option_kind("mip_heuristic_run_fpr"), Some(OptionKind::Bool));
        assert_eq!(
            option_kind("mip_heuristic_run_local_mip"),
            Some(OptionKind::Bool)
        );
        let mut highs = Highs::default();
        for name in ["mip_heuristic_run_fpr", "mip_heuristic_run_local_mip"] {
            assert_eq!(highs.set_option_value_bool(name, true), HighsStatus::kOk);
            assert_eq!(highs.get_option_value_bool(name), Some(true));
            assert_eq!(highs.set_option_value_bool(name, false), HighsStatus::kOk);
            assert_eq!(highs.get_option_value_bool(name), Some(false));
            assert_eq!(highs.set_option_value_str(name, "true"), HighsStatus::kOk);
        }
    }

    /// The `pdlp_*` option records (HighsOptions.h 1339-1385)
    /// are settable by name with their bounds.
    #[test]
    fn w3_pdlp_options_round_trip() {
        let mut highs = Highs::default();
        assert_eq!(
            highs.get_option_value_int("pdlp_iteration_limit"),
            Some(kHighsIInf)
        );
        assert_eq!(
            highs.set_option_value_int("pdlp_iteration_limit", 2000),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.pdlp_iteration_limit, 2000);
        assert_eq!(
            highs.set_option_value_int("pdlp_scaling_mode", 8),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_int("pdlp_restart_strategy", 3),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_int("pdlp_cupdlpc_restart_method", 3),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_int("pdlp_step_size_strategy", 2),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_int("pdlp_ruiz_iterations", -1),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_double("pdlp_optimality_tolerance", 1e-11),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_double("pdlp_optimality_tolerance", 1e-4),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.pdlp_optimality_tolerance, 1e-4);
        assert_eq!(
            highs.set_option_value_str("solver", "hipdlp"),
            HighsStatus::kOk
        );
    }

    #[test]
    fn lp_ipm_failure_policy_round_trip() {
        let mut highs = Highs::default();
        assert!(highs.options_.lp_ipm_failure_policy);
        assert_eq!(
            highs.set_option_value_bool("lp_ipm_failure_policy", false),
            HighsStatus::kOk
        );
        assert!(!highs.options_.lp_ipm_failure_policy);
    }

    #[test]
    fn w4_hipo_options_round_trip() {
        let mut highs = Highs::default();
        assert_eq!(highs.get_option_value_int("hipo_block_size"), Some(128));
        assert_eq!(
            highs.set_option_value_int("hipo_block_size", -1),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_int("hipo_block_size", 64),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.hipo_block_size, 64);
        assert_eq!(
            highs.set_option_value_str("hipo_ordering", "metis"),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.hipo_ordering, "metis");
        assert_eq!(
            highs.set_option_value_str("hipo_ordering", "nd"),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_str("hipo_system", "augmented"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("hipo_system", "both"),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_str("hipo_parallel_type", "tree"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("hipo_parallel_type", "choose"),
            HighsStatus::kError
        );
        // optionSolverOk accepts "hipo" only when HiPO is available.
        assert_eq!(
            highs.set_option_value_str("solver", "hipo"),
            HighsStatus::kError
        );
        highs.options_.hipo_available = true;
        assert_eq!(
            highs.set_option_value_str("solver", "hipo"),
            HighsStatus::kOk
        );
    }

    #[test]
    fn w1_options_round_trip() {
        let mut highs = Highs::default();
        assert_eq!(
            highs.get_option_value_int("iis_strategy"),
            Some(kIisStrategyDefault)
        );
        assert_eq!(
            highs.set_option_value_int("iis_strategy", 6),
            HighsStatus::kOk
        );
        assert_eq!(highs.get_option_value_int("iis_strategy"), Some(6));
        assert_eq!(
            highs.set_option_value_int("iis_strategy", 32),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_int("qp_iteration_limit", -1),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_int("qp_nullspace_limit", 10),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.qp_nullspace_limit, 10);
        assert_eq!(
            highs.set_option_value_double("qp_regularization_value", 1e-6),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_double("icrash_starting_weight", 1e-11),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_bool("icrash", true),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("icrash_strategy", "penalty"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_str("icrash_strategy").as_deref(),
            Some("penalty")
        );
        assert_eq!(
            highs.set_option_value_int("icrash_iterations", 201),
            HighsStatus::kError
        );
    }

    #[test]
    fn changing_iis_strategy_clears_a_valid_iis() {
        let mut highs = Highs::default();
        highs.iis_.valid_ = true;
        highs.iis_.strategy_ = kIisStrategyDefault;
        highs.iis_.row_index_ = vec![0];
        assert_eq!(
            highs.set_option_value_int("iis_strategy", 2),
            HighsStatus::kOk
        );
        assert!(!highs.iis_.valid_ && highs.iis_.row_index_.is_empty());
    }

    /// Gap-fill API options (HighsOptions.h records): defaults, bounds and the file-name
    /// exemption from lower-casing.
    #[test]
    fn g1_options_round_trip() {
        let mut highs = Highs::default();
        assert_eq!(highs.get_option_value_int("user_objective_scale"), Some(0));
        assert_eq!(highs.get_option_value_int("user_bound_scale"), Some(0));
        assert_eq!(highs.get_option_value_bool("use_warm_start"), Some(true));
        assert_eq!(
            highs.get_option_value_bool("blend_multi_objectives"),
            Some(true)
        );
        assert_eq!(highs.get_option_value_bool("solve_relaxation"), Some(false));
        assert_eq!(highs.get_option_value_int("mip_max_start_nodes"), Some(500));
        assert_eq!(
            highs.set_option_value_int("user_objective_scale", -3),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.user_objective_scale, -3);
        assert_eq!(
            highs.set_option_value_int("user_bound_scale", 4),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.user_bound_scale, 4);
        assert_eq!(
            highs.set_option_value_int("mip_max_start_nodes", -1),
            HighsStatus::kError
        );
        assert_eq!(
            highs.set_option_value_str("use_warm_start", "false"),
            HighsStatus::kOk
        );
        assert!(!highs.options_.use_warm_start);
        assert_eq!(
            highs.set_option_value_str("mip_improving_solution_file", "Sols.TXT"),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.mip_improving_solution_file, "Sols.TXT");
    }

    #[test]
    fn x1_opt_defaults_off() {
        let opts = crate::lp_data::highs_options::HighsOptions::default();
        assert!(!opts.mip_sb_budget_pre_restart);
        assert!(!opts.mip_dfs_node_limit);
        assert!(!opts.mip_reduced_stack_fast);
        assert!(!opts.mip_cdouble_act_contrib);
        assert!(!opts.mip_keep_branch_reason);
        assert!(!opts.mip_redundant_vbd);
        assert!(!opts.mip_path_same_bds);
        assert!(!opts.mip_strengthen_small_contrib);
        assert!(!opts.mip_lurking_bounds_robust);
        assert!(!opts.mip_dual_subst_limits);
        assert!(!opts.mip_ranged_dual_implied);
        assert!(!opts.mip_machine_sched_sepa);
        assert!(!opts.mip_dual_fix_probing);
        assert!(!opts.mip_aggregate_implications);
        assert!(!opts.mip_clique_normalise);
        assert!(!opts.mip_coef_change_impl_bounds);
        assert!(!opts.mip_dual_side_lost);
        let highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_clique_normalise"),
            Some(false)
        );
    }

    #[test]
    fn x1_opt_set_by_name() {
        let mut highs = Highs::default();
        let names = [
            "mip_sb_budget_pre_restart",
            "mip_dfs_node_limit",
            "mip_reduced_stack_fast",
            "mip_cdouble_act_contrib",
            "mip_keep_branch_reason",
            "mip_redundant_vbd",
            "mip_path_same_bds",
            "mip_strengthen_small_contrib",
            "mip_lurking_bounds_robust",
            "mip_dual_subst_limits",
            "mip_ranged_dual_implied",
            "mip_machine_sched_sepa",
            "mip_dual_fix_probing",
            "mip_aggregate_implications",
            "mip_clique_normalise",
            "mip_coef_change_impl_bounds",
            "mip_dual_side_lost",
        ];
        for name in &names {
            assert_eq!(
                highs.set_option_value_str(name, "true"),
                HighsStatus::kOk,
                "set_option_value_str({name}, \"true\")"
            );
            assert_eq!(
                highs.get_option_value_bool(name),
                Some(true),
                "get_option_value_bool({name})"
            );
        }
        assert_eq!(
            highs.set_option_value_str("mip_no_such_option", "true"),
            HighsStatus::kError
        );
    }

    #[test]
    fn x1_opt_gaps_registered() {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_str("mip_rel_gap", "0"),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.mip_rel_gap, 0.0);
        assert_eq!(
            highs.set_option_value_str("mip_abs_gap", "1e-9"),
            HighsStatus::kOk
        );
        assert_eq!(highs.options_.mip_abs_gap, 1e-9);
        assert_eq!(
            highs.set_option_value_str("mip_rel_gap", "-1"),
            HighsStatus::kError
        );
    }

    #[test]
    fn x1b_opt_defaults_off() {
        let opts = crate::lp_data::highs_options::HighsOptions::default();
        assert!(!opts.mip_cdouble_proofs);
        assert!(!opts.mip_fix_last_active);
        assert!(!opts.mip_keep_implications);
        assert!(!opts.mip_presolve_clique_table);
        assert!(!opts.mip_presolve_row_append);
        assert!(!opts.mip_fm_presolve);
        assert!(!opts.mip_sym_hash64);
        let highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_cdouble_proofs"),
            Some(false)
        );
        assert_eq!(
            highs.get_option_value_bool("mip_fix_last_active"),
            Some(false)
        );
        assert_eq!(
            highs.get_option_value_bool("mip_keep_implications"),
            Some(false)
        );
        assert_eq!(
            highs.get_option_value_bool("mip_presolve_clique_table"),
            Some(false)
        );
    }

    #[test]
    fn x1b_opt_set_get_roundtrip() {
        let mut highs = Highs::default();
        let names = [
            "mip_cdouble_proofs",
            "mip_fix_last_active",
            "mip_keep_implications",
            "mip_presolve_clique_table",
            "mip_presolve_row_append",
            "mip_fm_presolve",
            "mip_sym_hash64",
        ];
        for name in &names {
            assert_eq!(
                highs.set_option_value_str(name, "true"),
                HighsStatus::kOk,
                "set_option_value_str({name}, \"true\")"
            );
            assert_eq!(
                highs.get_option_value_bool(name),
                Some(true),
                "get_option_value_bool({name})"
            );
        }
        assert_eq!(
            highs.set_option_value_str("mip_no_such_option", "true"),
            HighsStatus::kError
        );
    }

    /// Every `h16_*` option (HiGHS 1.16.0 backports; all bool, default false).
    const H16_OPT_NAMES: [&str; 70] = [
        "h16_ac_centring",
        "h16_ac_fix",
        "h16_add_no_modify_on_error",
        "h16_agg_term_zero",
        "h16_assess_small",
        "h16_big_act",
        "h16_bound_chain",
        "h16_bound_range_inf",
        "h16_cdouble",
        "h16_clique_links",
        "h16_col_stuffing",
        "h16_cost_perturb",
        "h16_crossover_x",
        "h16_cupdlp_time",
        "h16_cutpool_age",
        "h16_dep_eq_time_cap",
        "h16_dominated_cols",
        "h16_dual_free_list",
        "h16_dual_postsolve_rules",
        "h16_dup_row_undo",
        "h16_dupcol_warn",
        "h16_enum_sig64",
        "h16_eq_row_adds",
        "h16_feasrelax_sense",
        "h16_fj",
        "h16_frontier_merge",
        "h16_hash_common",
        "h16_hfactor_pivot",
        "h16_hipdlp_max",
        "h16_hipo_display",
        "h16_hipo_empty_rows",
        "h16_impl_eq",
        "h16_implint_bounds",
        "h16_int_bound_round",
        "h16_int_bound_rounding",
        "h16_int_round_assert_tol",
        "h16_ipm_gap_offset",
        "h16_ipx_free_rows",
        "h16_ipx_scaled",
        "h16_ipx_stall",
        "h16_kkt_effective_costs",
        "h16_max_duals",
        "h16_mip_start",
        "h16_mip_workers",
        "h16_ninfmin",
        "h16_nla_debug_gate",
        "h16_obj_cdouble",
        "h16_parallel_rows",
        "h16_postsolve_basis_check",
        "h16_presolve_rules",
        "h16_row_duals",
        "h16_fourier_motzkin",
        "h16_qp_bounded_start",
        "h16_qp_convexity",
        "h16_ray_simplex",
        "h16_rebuild_reason",
        "h16_repair_lp",
        "h16_row_presolve",
        "h16_row_redundant",
        "h16_scaling_report",
        "h16_singleton_loops",
        "h16_subst_int",
        "h16_sym_refine",
        "h16_unresolved_unbounded_or_infeasible",
        "h16_user_scale_offset",
        "h16_user_sol_obj",
        "h16_user_solution_assess",
        "h16_watch_link",
        "h16_zero_cost_singleton",
        "h16_sub_solver_flags",
    ];

    /// `h16_*` option test LP (the ipx_wrapper fixture): optimum 6 - 8/7 + 0.5 = 5.357142857...
    fn h16_opt_tiny_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 3,
            col_cost_: vec![3.0, -4.0],
            col_lower_: vec![0.0, -kHighsInf],
            col_upper_: vec![10.0, kHighsInf],
            row_lower_: vec![1.0, -kHighsInf, 2.0],
            row_upper_: vec![4.0, kHighsInf, 2.0],
            sense_: ObjSense::kMaximize,
            offset_: 0.5,
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 1, 2];
        lp.a_matrix_.value_ = vec![2.0, 3.0, 5.0, 7.0];
        lp
    }

    #[test]
    fn h16_opt_defaults_off() {
        let highs = Highs::default();
        for name in H16_OPT_NAMES {
            assert_eq!(highs.get_option_value_bool(name), Some(false), "{name}");
        }
        assert!(!crate::lp_data::highs_options::HighsOptions::default().h16_hipo_display);
    }

    /// `recommended()`: every `h16_*` option on except `H16_RECOMMENDED_OFF`, the
    /// `RECOMMENDED_MIP_ON` options on, the rest equal to the 1.15.1 defaults.
    #[test]
    fn h16_opt_recommended_profile() {
        use crate::lp_data::highs_options::{
            H16_RECOMMENDED_OFF, HighsOptions, RECOMMENDED_MIP_ON,
        };
        let highs = Highs::with_options(HighsOptions::recommended());
        for name in H16_OPT_NAMES {
            let want = !H16_RECOMMENDED_OFF.contains(&name);
            assert_eq!(highs.get_option_value_bool(name), Some(want), "{name}");
        }
        for name in RECOMMENDED_MIP_ON {
            assert_eq!(highs.get_option_value_bool(name), Some(true), "{name}");
        }
        let mut back = HighsOptions::recommended();
        back.set_all_h16(false);
        back.set_recommended_mip(false);
        assert_eq!(
            format!("{back:?}"),
            format!("{:?}", HighsOptions::default())
        );
    }

    #[test]
    fn h16_opt_set_get_roundtrip() {
        let mut highs = Highs::default();
        for name in H16_OPT_NAMES {
            assert_eq!(
                highs.set_option_value_str(name, "true"),
                HighsStatus::kOk,
                "set_option_value_str({name}, \"true\")"
            );
            assert_eq!(highs.get_option_value_bool(name), Some(true), "{name}");
            assert_eq!(
                highs.set_option_value_bool(name, false),
                HighsStatus::kOk,
                "{name}"
            );
            assert_eq!(highs.get_option_value_bool(name), Some(false), "{name}");
        }
        assert_eq!(
            highs.set_option_value_str("h16_no_such_option", "true"),
            HighsStatus::kError
        );
    }

    /// `h16_sub_solver_flags`: `copy_version_flags_from` carries
    /// every `h16_*` / `mip_*` / `mip_*` flag to a MIP-internal solver and nothing else.
    #[test]
    fn h16_sub_solver_flags_copy_version_flags() {
        use crate::lp_data::highs_options::HighsOptions;
        const X1_NAMES: [&str; 24] = [
            "mip_sb_budget_pre_restart",
            "mip_dfs_node_limit",
            "mip_reduced_stack_fast",
            "mip_cdouble_act_contrib",
            "mip_keep_branch_reason",
            "mip_redundant_vbd",
            "mip_path_same_bds",
            "mip_strengthen_small_contrib",
            "mip_lurking_bounds_robust",
            "mip_dual_subst_limits",
            "mip_ranged_dual_implied",
            "mip_machine_sched_sepa",
            "mip_dual_fix_probing",
            "mip_aggregate_implications",
            "mip_clique_normalise",
            "mip_coef_change_impl_bounds",
            "mip_dual_side_lost",
            "mip_cdouble_proofs",
            "mip_fix_last_active",
            "mip_keep_implications",
            "mip_presolve_clique_table",
            "mip_presolve_row_append",
            "mip_fm_presolve",
            "mip_sym_hash64",
        ];
        let mut src = Highs::default();
        for name in H16_OPT_NAMES.iter().chain(X1_NAMES.iter()) {
            assert_eq!(
                src.set_option_value_bool(name, true),
                HighsStatus::kOk,
                "{name}"
            );
        }
        // Non-version options must not travel with the flags.
        src.options_.mip_rel_gap = 0.5;
        src.options_.random_seed = 7;
        src.options_.mip_conflict_resolution = true;

        let mut dst = Highs::default();
        dst.options_.copy_version_flags_from(&src.options_);
        for name in H16_OPT_NAMES.iter().chain(X1_NAMES.iter()) {
            assert_eq!(dst.get_option_value_bool(name), Some(true), "{name}");
        }
        let mut back = dst.options_.clone();
        back.copy_version_flags_from(&HighsOptions::default());
        assert_eq!(
            format!("{back:?}"),
            format!("{:?}", HighsOptions::default())
        );
    }

    /// HiGHS 1.16 (#3592): ipm/hipo/ipm/Solver.cpp lines 69-70. Flag off keeps the
    /// 1.15.1 display rule; flag on sets display from output_flag alone.
    #[test]
    fn h16_opt_hipo_display_gate() {
        let mut opts = crate::lp_data::highs_options::HighsOptions::default();
        opts.output_flag = true;
        opts.log_to_console = false;
        let mut solver = crate::ipm::hipo::ipm::solver::Solver::default();
        solver.setOptions(&opts);
        assert!(!solver.options_.display);
        opts.h16_hipo_display = true;
        solver.setOptions(&opts);
        assert!(solver.options_.display);
        opts.output_flag = false;
        solver.setOptions(&opts);
        assert!(!solver.options_.display);
    }

    /// The `h16_*` HiPO options reach the HiPO Options (setOptions) and Model (load).
    #[test]
    fn h16_lp_plumb_hipo_fields() {
        let lp = h16_opt_tiny_lp();
        let mut opts = crate::lp_data::highs_options::HighsOptions::default();
        opts.h16_ipx_stall = true;
        opts.h16_ipx_scaled = true;
        opts.h16_ipm_gap_offset = true;
        opts.h16_ipx_free_rows = true;
        opts.h16_hipo_empty_rows = true;
        let mut solver = crate::ipm::hipo::ipm::solver::Solver::default();
        solver.setOptions(&opts);
        assert!(solver.options_.h16_ipx_stall);
        assert!(solver.options_.h16_ipx_scaled);
        assert!(solver.options_.h16_ipm_gap_offset);
        assert!(solver.options_.h16_ipx_free_rows);
        assert!(solver.options_.h16_hipo_empty_rows);
        assert_eq!(
            solver.load(&lp, &crate::model::highs_hessian::HighsHessian::default()),
            crate::ipm::hipo::ipm::status::kStatusOk
        );
        assert!(solver.model_.h16_ipx_free_rows_);
        assert!(solver.model_.h16_ipm_gap_offset_);
        assert!(solver.model_.h16_hipo_empty_rows_);
    }

    /// The `h16_hfactor_pivot` flag reaches HFactor after HSimplexNla::setup.
    #[test]
    fn h16_lp_plumb_hfactor_flag() {
        let lp = h16_opt_tiny_lp();
        let basic_index = vec![2, 3, 4];
        let mut opts = crate::lp_data::highs_options::HighsOptions::default();
        let mut nla = crate::simplex::nla::HSimplexNla::default();
        nla.setup(&lp, &basic_index, &opts, &lp.a_matrix_, 0.1);
        assert!(!nla.factor_.h16_hfactor_pivot);
        opts.h16_hfactor_pivot = true;
        nla.setup(&lp, &basic_index, &opts, &lp.a_matrix_, 0.1);
        assert!(nla.factor_.h16_hfactor_pivot);
    }

    /// The `h16_hipdlp_max` flag reaches PrimalDualParams after PDLPSolver::setup.
    #[test]
    fn h16_lp_plumb_pdlp_flag() {
        let timer = crate::util::highs_timer::HighsTimer::new();
        let mut opts = crate::lp_data::highs_options::HighsOptions::default();
        let mut pdlp = crate::pdlp::hipdlp::pdhg::PDLPSolver::default();
        pdlp.setup(&opts, &timer);
        assert!(!pdlp.params_.h16_hipdlp_max);
        opts.h16_hipdlp_max = true;
        pdlp.setup(&opts, &timer);
        assert!(pdlp.params_.h16_hipdlp_max);
    }

    /// Every solver path on the tiny LP with all `h16_*` flags at their defaults gives the
    /// optimum. (No all-on comparison: the flags stop being inert as their jobs land.)
    #[test]
    fn h16_lp_plumb_default_solve() {
        let expected = 6.0 - 8.0 / 7.0 + 0.5;
        for solver_name in ["simplex", "ipm", "hipo"] {
            let mut highs = Highs::default();
            if solver_name == "hipo" {
                highs.options_.hipo_available = true;
            }
            assert_eq!(
                highs.set_option_value_bool("log_to_console", false),
                HighsStatus::kOk
            );
            assert_eq!(
                highs.set_option_value_str("presolve", "off"),
                HighsStatus::kOk
            );
            assert_eq!(
                highs.set_option_value_str("solver", solver_name),
                HighsStatus::kOk
            );
            assert_eq!(highs.pass_model(h16_opt_tiny_lp()), HighsStatus::kOk);
            assert_eq!(highs.run(), HighsStatus::kOk, "{solver_name}");
            assert_eq!(
                highs.get_model_status(),
                HighsModelStatus::kOptimal,
                "{solver_name}"
            );
            let objective = highs.get_info().objective_function_value;
            assert!(
                (objective - expected).abs() < 1e-6,
                "{solver_name}: {objective}"
            );
        }
    }

    #[test]
    fn s2_opt_default_false() {
        let opts = crate::lp_data::highs_options::HighsOptions::default();
        assert!(!opts.mip_conflict_resolution);
        let highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_conflict_resolution"),
            Some(false)
        );
    }

    #[test]
    fn s2_opt_set_get() {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("mip_conflict_resolution", true),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_conflict_resolution"),
            Some(true)
        );
        assert_eq!(
            highs.set_option_value_bool("mip_conflict_resolution", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_conflict_resolution"),
            Some(false)
        );
    }

    #[test]
    fn s2_opt_kind_bool() {
        assert_eq!(
            option_kind("mip_conflict_resolution"),
            Some(OptionKind::Bool)
        );
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_str("mip_conflict_resolution", "true"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_conflict_resolution"),
            Some(true)
        );
    }

    #[test]
    fn s2_tree_only_opt_default_off() {
        assert!(
            !crate::lp_data::highs_options::HighsOptions::default()
                .mip_conflict_resolution_tree_only
        );
        assert_eq!(
            option_kind("mip_conflict_resolution_tree_only"),
            Some(OptionKind::Bool)
        );
        let mut highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_conflict_resolution_tree_only"),
            Some(false)
        );
        assert_eq!(
            highs.set_option_value_str("mip_conflict_resolution_tree_only", "true"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_conflict_resolution_tree_only"),
            Some(true)
        );
    }

    #[test]
    fn s3_opt_default_off() {
        assert!(!crate::lp_data::highs_options::HighsOptions::default().mip_symmetry_dynamic);
        let highs = Highs::default();
        assert_eq!(
            highs.get_option_value_bool("mip_symmetry_dynamic"),
            Some(false)
        );
    }

    #[test]
    fn s3_opt_set_get() {
        assert_eq!(option_kind("mip_symmetry_dynamic"), Some(OptionKind::Bool));
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("mip_symmetry_dynamic", true),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_symmetry_dynamic"),
            Some(true)
        );
        assert_eq!(
            highs.set_option_value_bool("mip_symmetry_dynamic", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.get_option_value_bool("mip_symmetry_dynamic"),
            Some(false)
        );
        assert_eq!(
            highs.set_option_value_str("mip_symmetry_dynamic", "true"),
            HighsStatus::kOk
        );
    }
}
