//! Table-driven option registry — port of HiGHS 1.16.0 `lp_data/HighsOptions.h`
//! (`HighsOptionsRecord` + `HighsOptions::initRecords`) plus the listing and file pieces of
//! `lp_data/HighsOptions.cpp` (`reportOptions`, `writeOptionsToFile`) and `io/LoadOptions.cpp`
//! (`loadOptionsFromFile`), MIT (see THIRD_PARTY_NOTICES.md).
//!
//! One table describing
//! every option — name, help text, type, default, range, advanced flag and source — is the
//! single source of truth for by-name access, `rusolve options` (table or JSON), options-file
//! read/write, and actionable unknown-option errors.
//!
//! Scope decisions:
//! * A record exists here exactly when the option is wired into this port's
//!   [`HighsOptions`]: the 27 of HiGHS' 163 1.16.0 records whose value no ported code reads
//!   (ranging, file I/O, crash weights, HiPO strategy) are omitted rather than exposed as
//!   silently-dead knobs — settable ⟺ wired. rusolve's own levers (HiGHS 1.16 / post-1.15.1 backports, SCIP ports,
//!   original techniques, trace and dump switches) join the same table, so one interface
//!   covers them all and each one's `source` names where it comes from.
//! * Validation stays in the HiGHS-faithful typed setters (`highs_api_mip.rs`), which read
//!   type and range from this table. One path validates, so `--opt`, an options file and
//!   `Highs::set_option_value_*` cannot drift apart.
//! * The table never feeds [`HighsOptions::default`] or [`HighsOptions::recommended`]: the
//!   default search is untouched (identity 20/20).
//!
//! Bounds and defaults are verbatim 1.16.0 record values, verified equal to the ranges the
//! port's previous hand-written tables enforced for every option that was settable before.

use serde::Serialize;

use super::highs_options::HighsOptions;

/// Registry entry kind: the HiGHS record type plus the record's declared default and range
/// (`setHighsOptionsRecord` argument order is name, advanced, current, lower, upper).
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum RecordKind {
    Bool {
        default: bool,
    },
    Int {
        lower: i32,
        default: i32,
        upper: i32,
    },
    Double {
        lower: f64,
        default: f64,
        upper: f64,
    },
    Str {
        default: &'static str,
    },
}

/// What an option can influence, so an agent can tell a search knob from a tolerance from an
/// output switch before spending time on a sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Affects {
    /// Changes the branch-and-bound / simplex / barrier search: sweeping it can move the
    /// answer or the time to it.
    Search,
    /// A tolerance or convergence threshold: numeric precision, not strategy.
    Tolerance,
    /// Changes what the solve means or produces rather than the path it takes.
    Semantics,
    /// Output, files, logging, profiling.
    Output,
    /// Threads, timing, memory, build capabilities.
    Runtime,
}

/// One row of the registry.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct OptionRecord {
    pub name: &'static str,
    /// Help text: verbatim HiGHS 1.16.0 for ported records; for rusolve levers, what the
    /// lever controls, citing the HiGHS PR, paper or solver it follows.
    pub description: &'static str,
    /// HiGHS `advanced` flag. rusolve levers are all advanced: they are not HiGHS behaviour.
    pub advanced: bool,
    pub kind: RecordKind,
    /// Where the option comes from: "HiGHS 1.16.0", a HiGHS backport, another solver (SCIP,
    /// mip-heuristics), or "rusolve" for an original lever.
    pub source: &'static str,
    pub affects: Affects,
    /// Worth including in an agent tuning sweep.
    pub tunable: bool,
}

/// A current option value, by the four HiGHS option types.
#[derive(Debug, Clone, PartialEq)]
pub enum OptionValue {
    Bool(bool),
    Int(i32),
    Double(f64),
    Str(String),
}

impl OptionValue {
    /// The short rendering used by options files and the plain listing. Doubles print with
    /// Rust's shortest round-trip form, so write → read → write is stable.
    pub fn to_token(&self) -> String {
        match self {
            OptionValue::Bool(value) => value.to_string(),
            OptionValue::Int(value) => value.to_string(),
            OptionValue::Double(value) => value.to_string(),
            OptionValue::Str(value) => value.clone(),
        }
    }

    /// HiGHS `optionEntryTypeToString` spellings, kept so tools see the same words.
    pub fn type_name(&self) -> &'static str {
        match self {
            OptionValue::Bool(_) => "bool",
            OptionValue::Int(_) => "HighsInt",
            OptionValue::Double(_) => "double",
            OptionValue::Str(_) => "string",
        }
    }
}

/// JSON form of a record plus its current value. Doubles that are not finite (range bounds
/// like `inf`) are emitted as strings: JSON numbers cannot hold them.
#[derive(Debug, Serialize)]
struct JsonRecord {
    name: &'static str,
    type_name: &'static str,
    description: &'static str,
    advanced: bool,
    source: &'static str,
    affects: &'static str,
    tunable: bool,
    default: serde_json::Value,
    lower: Option<serde_json::Value>,
    upper: Option<serde_json::Value>,
    value: serde_json::Value,
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    changed: bool,
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    affects_search: bool,
}

fn json_number(value: f64) -> serde_json::Value {
    if value.is_finite() {
        serde_json::Value::from(value)
    } else if value.is_infinite() {
        serde_json::Value::from(if value > 0.0 { "inf" } else { "-inf" })
    } else {
        serde_json::Value::from("nan")
    }
}

/// The registry: one record per option wired into this port, sorted by name (HiGHS keeps its
/// option index sorted too).
pub static OPTION_RECORDS: &[OptionRecord] = &[
    OptionRecord {
        name: "allow_unbounded_or_infeasible",
        description: "Allow ModelStatus::kUnboundedOrInfeasible",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "allowed_cost_scale_factor",
        description: "Largest power-of-two factor permitted when scaling the costs",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 20,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "allowed_matrix_scale_factor",
        description: "Largest power-of-two factor permitted when scaling the constraint matrix",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 20,
            upper: 30,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "blend_multi_objectives",
        description: "Blend multiple objectives or apply lexicographically",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "centring_ratio_tolerance",
        description: "Centring stops when the ratio max(x_js_j) / min(x_js_j) is below this tolerance (default = 100)",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 100.0,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "cost_scale_factor",
        description: "Scaling factor for costs",
        advanced: true,
        kind: RecordKind::Int {
            lower: -20,
            default: 0,
            upper: 20,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "dual_feasibility_tolerance",
        description: "Dual feasibility tolerance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 1e-07,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "dual_residual_tolerance",
        description: "Dual residual tolerance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 1e-07,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "dual_simplex_cost_perturbation_multiplier",
        description: "Dual simplex cost perturbation multiplier: 0 => no perturbation",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 1.0,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "dual_simplex_pivot_growth_tolerance",
        description: "Dual simplex pivot growth tolerance",
        advanced: true,
        kind: RecordKind::Double {
            lower: 1e-12,
            default: 1e-09,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "dual_steepest_edge_weight_error_tolerance",
        description: "Tolerance on dual steepest edge weight errors",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: f64::INFINITY,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "dual_steepest_edge_weight_log_error_threshold",
        description: "Threshold on dual steepest edge weight errors for Devex switch",
        advanced: true,
        kind: RecordKind::Double {
            lower: 1.0,
            default: 10.0,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "factor_pivot_threshold",
        description: "Matrix factorization pivot threshold",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0008,
            default: 0.1,
            upper: 0.5,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "factor_pivot_tolerance",
        description: "Matrix factorization pivot tolerance",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 1e-10,
            upper: 1.0,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "h16_ac_centring",
        description: "Ported behind this default-off lever: HiGHS #3456 #3522",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_ac_fix",
        description: "Ported behind this default-off lever: HiGHS #3408",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_add_no_modify_on_error",
        description: "Ported behind this default-off lever: HiGHS #3386",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_agg_term_zero",
        description: "Ported behind this default-off lever: HiGHS #3507",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_assess_small",
        description: "Ported behind this default-off lever: HiGHS #3169 #3382",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_big_act",
        description: "Ported behind this default-off lever: HiGHS #3267",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_bound_chain",
        description: "Ported behind this default-off lever: HiGHS #3196",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_bound_range_inf",
        description: "Ported behind this default-off lever: HiGHS #3525",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_cdouble",
        description: "Ported behind this default-off lever: HiGHS #3149",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_clique_links",
        description: "Ported behind this default-off lever: HiGHS #3540",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_col_stuffing",
        description: "Ported behind this default-off lever: HiGHS #3554",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_cost_perturb",
        description: "Ported behind this default-off lever: HiGHS #3226",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_crossover_x",
        description: "Ported behind this default-off lever: HiGHS #3466 #3185",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_cupdlp_time",
        description: "Ported behind this default-off lever: HiGHS #3450",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_cutpool_age",
        description: "Ported behind this default-off lever: HiGHS #3574",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_dep_eq_time_cap",
        description: "Ported behind this default-off lever: HiGHS #3167 #3398",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_dominated_cols",
        description: "Ported behind this default-off lever: HiGHS #3411 #3285 #3290 #3472",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_dual_free_list",
        description: "Ported behind this default-off lever: HiGHS #3353",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_dual_postsolve_rules",
        description: "Ported behind this default-off lever: HiGHS #3342 #3348 #2962 #3306 #3455",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_dup_row_undo",
        description: "Ported behind this default-off lever: HiGHS #3313",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_dupcol_warn",
        description: "Ported behind this default-off lever: HiGHS #3580",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_enum_sig64",
        description: "Ported behind this default-off lever: HiGHS #3529",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_eq_row_adds",
        description: "Ported behind this default-off lever: HiGHS #3324",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_feasrelax_sense",
        description: "Ported behind this default-off lever: HiGHS #3370",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_fj",
        description: "Ported behind this default-off lever: HiGHS #3310 #3269",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_fourier_motzkin",
        description: "Ported behind this default-off lever: HiGHS #3218",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_frontier_merge",
        description: "Ported behind this default-off lever: HiGHS #3443",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_hash_common",
        description: "Ported behind this default-off lever: HiGHS #3327",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_hfactor_pivot",
        description: "Ported behind this default-off lever: HiGHS #3352 #3354",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_hipdlp_max",
        description: "Ported behind this default-off lever: HiGHS #3469",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_hipo_display",
        description: "Ported behind this default-off lever: HiGHS #3592",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_hipo_empty_rows",
        description: "Ported behind this default-off lever: HiGHS #3376",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_impl_eq",
        description: "Ported behind this default-off lever: HiGHS #3385",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_implint_bounds",
        description: "Ported behind this default-off lever: HiGHS #3276",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_int_bound_round",
        description: "Ported behind this default-off lever: HiGHS #3320",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_int_bound_rounding",
        description: "Ported behind this default-off lever: HiGHS #3411 #3285 #3290 #3472",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_int_round_assert_tol",
        description: "Ported behind this default-off lever: HiGHS #3576",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Tolerance,
        tunable: false,
    },
    OptionRecord {
        name: "h16_ipm_gap_offset",
        description: "Ported behind this default-off lever: HiGHS #3427",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_ipx_free_rows",
        description: "Ported behind this default-off lever: HiGHS #3319",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_ipx_scaled",
        description: "Ported behind this default-off lever: HiGHS #3513",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_ipx_stall",
        description: "Ported behind this default-off lever: HiGHS #3401",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_kkt_effective_costs",
        description: "Ported behind this default-off lever: HiGHS #3299",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_max_duals",
        description: "Ported behind this default-off lever: HiGHS #3248",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_mip_start",
        description: "Ported behind this default-off lever: HiGHS #3135 #3523 #3272",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_mip_workers",
        description: "Ported behind this default-off lever: HiGHS #3124 #3192 #3237 #3336a #3406 #3521 #3578",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_ninfmin",
        description: "Ported behind this default-off lever: HiGHS #3275",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_nla_debug_gate",
        description: "Ported behind this default-off lever: HiGHS #3515",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "h16_obj_cdouble",
        description: "Ported behind this default-off lever: HiGHS #3513",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_parallel_rows",
        description: "Ported behind this default-off lever: HiGHS #3358 #3520 #3261",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_postsolve_basis_check",
        description: "Ported behind this default-off lever: HiGHS #3350 #3421",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_presolve_rules",
        description: "Ported behind this default-off lever: HiGHS #3167 #3210 #3300 #3536 #3334 #3342-family #3442 #3308 #3498",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_qp_bounded_start",
        description: "Ported behind this default-off lever: HiGHS #3153",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_qp_convexity",
        description: "Ported behind this default-off lever: HiGHS #3325 #3458",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_ray_simplex",
        description: "Ported behind this default-off lever: HiGHS #3435",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_rebuild_reason",
        description: "Ported behind this default-off lever: HiGHS #3351",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_repair_lp",
        description: "Ported behind this default-off lever: HiGHS #3181",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_row_duals",
        description: "Ported behind this default-off lever: HiGHS #3270-family #3367",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_row_presolve",
        description: "Ported behind this default-off lever: HiGHS #3260 #3392 #3224",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_row_redundant",
        description: "Ported behind this default-off lever: HiGHS #3575",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_scaling_report",
        description: "Ported behind this default-off lever: HiGHS #3350",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_singleton_loops",
        description: "Ported behind this default-off lever: HiGHS #3139",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_sub_solver_flags",
        description: "rusolve port plumbing behind this default-off lever: MIP sub-solvers copy the parent's 1.16.0 version flags",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_subst_int",
        description: "Ported behind this default-off lever: HiGHS #3268",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_sym_refine",
        description: "Ported behind this default-off lever: HiGHS #3182",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_unresolved_unbounded_or_infeasible",
        description: "h16 unresolved unbounded or infeasible: rusolve default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_user_scale_offset",
        description: "Ported behind this default-off lever: HiGHS #3186",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_user_sol_obj",
        description: "Ported behind this default-off lever: HiGHS #3492",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_user_solution_assess",
        description: "Ported behind this default-off lever: HiGHS #3556",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_watch_link",
        description: "Ported behind this default-off lever: HiGHS #3571",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "h16_zero_cost_singleton",
        description: "Ported behind this default-off lever: HiGHS #3342 #3348 #2962 #3306 #3455",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "highs_analysis_level",
        description: "Analysis level in HiGHS",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 511,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "highs_debug_level",
        description: "Debugging level in HiGHS",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 3,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "hipo_block_size",
        description: "Block size for dense linear algebra within HiPO",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 128,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "hipo_ordering",
        description: "HiPO matrix reordering: \"choose\", \"metis\", \"amd\" or \"rcm\"",
        advanced: false,
        kind: RecordKind::Str { default: "choose" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "hipo_parallel_type",
        description: "HiPO parallelism: \"tree\", \"node\", \"both\" or \"choose\"",
        advanced: false,
        kind: RecordKind::Str { default: "both" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "hipo_system",
        description: "HiPO Newton system: \"choose\", \"augmented\" or \"normaleq\"",
        advanced: false,
        kind: RecordKind::Str { default: "choose" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "icrash",
        description: "Run iCrash",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "icrash_approx_iter",
        description: "iCrash approximate minimization iterations",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 50,
            upper: 100,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "icrash_breakpoints",
        description: "Exact subproblem solution for iCrash",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "icrash_dualize",
        description: "Dualize strategy for iCrash",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "icrash_exact",
        description: "Exact subproblem solution for iCrash",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "icrash_iterations",
        description: "iCrash iterations",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 30,
            upper: 200,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "icrash_starting_weight",
        description: "iCrash starting weight",
        advanced: true,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 0.001,
            upper: 1e+50,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "icrash_strategy",
        description: "Strategy for iCrash",
        advanced: true,
        kind: RecordKind::Str { default: "ICA" },
        source: "HiGHS 1.16.0",
        affects: Affects::Output,
        tunable: false,
    },
    OptionRecord {
        name: "iis_strategy",
        description: "Strategy for IIS calculation: 0 => Light test; 1 => Try dual ray; 2 => Try elastic LP; 4 => Find true IIS; 8 => Prioritise columns; 16 => Find relaxation IIS for MIP",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 31,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "iis_time_limit",
        description: "Time limit for computing IIS (seconds)",
        advanced: false,
        kind: RecordKind::Double {
            lower: 0.0,
            default: f64::INFINITY,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "infinite_bound",
        description: "Limit on variable/constraint bound: values greater than or equal to this will be treated as infinite",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1000000000000000.0,
            default: 1e+20,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "infinite_cost",
        description: "Limit on cost coefficient: values greater than or equal to this will be treated as infinite",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1000000000000000.0,
            default: 1e+20,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "ipm_iteration_limit",
        description: "Iteration limit for IPM solver",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2147483647,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "ipm_optimality_tolerance",
        description: "IPM optimality tolerance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-12,
            default: 1e-08,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "ipx_dualize_strategy",
        description: "Strategy for dualizing before IPX",
        advanced: true,
        kind: RecordKind::Int {
            lower: -1,
            default: 2,
            upper: 3,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "kkt_tolerance",
        description: "If changed from its default value, this tolerance is used for all feasibility and optimality (KKT) measures",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 1e-07,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "large_matrix_value",
        description: "Upper limit on matrix entries: values greater than or equal to this will be treated as infinite",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1.0,
            default: 1000000000000000.0,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "less_infeasible_DSE_check",
        description: "Check whether LP is candidate for LiDSE",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "less_infeasible_DSE_choose_row",
        description: "Use LiDSE if LP has right properties",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "log_dev_level",
        description: "Output development messages: 0 => none; 1 => info; 2 => detailed; 3 => verbose",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 3,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Output,
        tunable: false,
    },
    OptionRecord {
        name: "log_file",
        description: "Log file",
        advanced: false,
        kind: RecordKind::Str { default: "" },
        source: "HiGHS 1.16.0",
        affects: Affects::Output,
        tunable: false,
    },
    OptionRecord {
        name: "log_to_console",
        description: "Enables or disables console logging",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Output,
        tunable: false,
    },
    OptionRecord {
        name: "lp_fold",
        description: "Flag. Fold the LP by colour refinement, solve the quotient, lift and cross over",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_fold_basis",
        description: "Default-off lever. With lp_fold: lift the quotient basis and finish with simplex instead of crossover",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_dualize_auto",
        description: "Flag. Solve the presolved LP through its dual in simplex when it has at least 3x more rows than columns (Andersen-Andersen 1995, Gondzio 1997)",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_crossover_ipm_crash",
        description: "Default-off lever. IPX crossover from a supplied interior point (the HiPO hand-off): crash basis ranks candidates by their distance to the bounds instead of by sparsity only",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_crossover_batch_push",
        description: "Default-off lever. IPX crossover push phase: superbasics are pushed to their bounds in groups with one FTRAN per group when no basic variable blocks, falling back to single pushes",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_dual_devex_start",
        description: "Default-off lever. Dual simplex from a non-logical basis uses Devex instead of exact steepest edge weights",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_concurrent",
        description: "Flag. With threads >= 2, race dual simplex against the HiPO interior point plus crossover on the LP; the first optimal wins",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_solver_auto",
        description: "Flag. With one thread and solver=choose, pick the HiPO interior point over simplex when the LP's structure predicts it wins; default on",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_ipm_support_reduction",
        description: "Default-off lever. After HiPO converges and before crossover, extra barrier iterations on a randomly cost-perturbed problem move the point to a low-support vertex of the optimal face",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_postsolve_snap_basis",
        description: "Default-off lever. After postsolve, a nonbasic boxed variable gets the status of the bound its recovered value is at",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_postsolve_fix_tiny",
        description: "Flag. The clean-up solve after postsolve fixes columns whose bound range is within the primal feasibility tolerance at their lower bound",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_crossover_ipx_refine",
        description: "Flag. HiPO then IPX crossover: IPX iterates on HiPO's point to a 100x tighter gap (1e-10) before the push phase",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "ipm_hipo_fast_blas",
        description: "Flag. HiPO uses a packed dgemm and grouped dtrsm instead of reference BLAS: AVX2/FMA or NEON 8x6 tile, 1.5-2x faster HiPO on fome13/nug08-3rd, not bit-identical",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "ipm_hipo_cache_perm",
        description: "Flag. HiPO replays a recorded permutation and column sort of the matrix instead of redoing it each factorisation; bit-identical",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: false,
    },
    OptionRecord {
        name: "ipm_hipo_par_ne",
        description: "Flag. HiPO builds the normal-equations matrix on several threads, same summation order as serial, bit-identical",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: false,
    },
    OptionRecord {
        name: "ipm_hipo_dualize",
        description: "Flag. HiPO solves the dual of an LP with at least 3x more rows than columns, then maps the dual crossover basis back to the primal (Andersen-Andersen 1995, Gondzio 1997)",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "ipm_hipo_dense_cols",
        description: "Default-off lever. HiPO keeps dense columns out of the normal-equations factor and recovers them by preconditioned CG (Andersen 1996; negative result)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "ipm_hipo_dense_schur",
        description: "Default-off lever. HiPO keeps up to 128 dense columns out of the normal-equations factor and adds them back by a Woodbury update inside CG (Andersen 1996; negative result)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "ipm_hipo_stall_detect",
        description: "Default-off lever. HiPO stops when dual infeasibility stays flat while mu falls by 1e4 over 5 iterations, and the LP goes to IPX",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "ipm_hipo_amd_first",
        description: "Flag. HiPO analyses NE first and skips AS when NE is cheap; Metis only if AMD flops >= 1e9; default on",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_ipx_soft_drop",
        description: "Flag. IPX drops a nonbasic only if z*x over its finite bounds is also within the drop tolerance (not on MIP analytic-centre solves)",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_dep_eq_fixed_budget",
        description: "Default-off lever. The dependent-equations search of presolve gets its fixed 100 s cap instead of 1 % of time_limit, so presolve does not change with the time limit",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_ipx_no_drop",
        description: "Default-off lever. IPX keeps every basic variable in the interior point solve (never drops one as implied), so the postprocessed point stays complementary for crossover",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_ipm_failure_policy",
        description: "Flag. Cold IPX restart after a HiPO failure, crossover on a near-optimal IPX stall, simplex on IPM error",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "lp_presolve_requires_basis_postsolve",
        description: "Prevents LP presolve steps for which postsolve cannot maintain a basis: redundant option!",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "max_centring_steps",
        description: "Maximum number of steps to use (default = 5) when computing the analytic centre",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 5,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "max_dual_simplex_cleanup_level",
        description: "Max level of dual simplex cleanup",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 1,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "max_dual_simplex_phase1_cleanup_level",
        description: "Max level of dual simplex phase 1 cleanup",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 2,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_abs_gap",
        description: "Tolerance on absolute gap of MIP, ub-lb, to determine whether optimality has been reached for a MIP instance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 1e-06,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "mip_allow_cut_separation_at_nodes",
        description: "Whether cut separation at nodes other than the root node is permitted",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_allow_restart",
        description: "Whether MIP restart is permitted",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_alns",
        description: "Flag. Bandit-scheduled LNS over RENS, RINS, mutation and crossover sub-MIPs on a work budget (Hendel 2022; SCIP heur_alns.c); default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "SCIP 11.0.0 / Hendel 2022",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_presolve_clique_varbound",
        description: "Default-off lever. Presolve merges variable-bound rows z + a_j x_j >= L whose binaries lie in one clique into a single row (Achterberg et al. 2020)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_presolve_xor_gf2",
        description: "Flag. Presolve solves parity (XOR) rows over GF(2) by Gaussian elimination, fixing forced binaries",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_orbitope_full_fix",
        description: "Default-off lever. At the root, apply lexicographic orbitopal fixing (Kaibel-Peinhardt-Pfetsch 2011) over all rows of every non-packing full orbitope",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_orbitope_root_fix",
        description: "Flag. Fix the staircase of packing orbitopes at the root after the first fractional LP (when large enough); disjoint dynamic symmetries stay",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_presolve_eq_subst",
        description: "Default-off lever. Presolve substitutes a unit-coefficient integer column out of a short integer equation and adds its non-implied bounds as one row",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_shift_and_propagate",
        description: "Default-off lever. Pre-root Shift-and-Propagate primal heuristic, only without an incumbent",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_fj_native",
        description: "Default-off lever. Run the native engine's Feasibility Jump (LP-free local search) in place of the HiGHS FJ",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_first_inc",
        description: "Default-off lever. Without an incumbent after the pre-root Feasibility Jump, run the native FJ and then the FPR and Local-MIP presolve heuristics before the root LP",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_first_inc_root",
        description: "Default-off lever. Without an incumbent after the first root LP and its rounding heuristics (before the cut loop), run the native FJ and then the FPR and Local-MIP presolve heuristics",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_activation_dff",
        description: "Default-off lever. Before presolve add up to two dual-feasible-function rows sum f(a_i) x_i <= y per activation knapsack sum a_i x_i <= b y of a model with an objective",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_decompose",
        description: "Default-off lever. Solve the connected components of a MIP as separate MIPs",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_lattice",
        description: "Flag. Small pure general-integer MIPs whose columns LLL shortens a lot: solve in lattice coordinates x = x0 + U y",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_lattice_kernel",
        description: "Flag. Small equality-knapsack MIPs with integral data (market split): solve in the Aardal-Hurkens-Lenstra kernel-lattice coordinates x = x0 + Q y",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_zero_obj_cdcl",
        description: "Flag. All-zero objective: CDCL feasibility search at the root",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_cdcl_uip",
        description: "Default-off lever. CDCL dive on a first-UIP, two-watched-literal engine, pure-binary models",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_cdcl_dive",
        description: "Default-off lever. LP-free conflict-driven DFS with Luby restarts as a primal heuristic (pre-root and after the root) on pure-integer models; for feasibility-type trees",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_root_bound_stall",
        description: "Flag. Default on. Stop root cut separation after 5 rounds with an unchanged root bound once an incumbent exists within 2% of it, so RENS and sub-MIP run sooner",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heur_lp_quota",
        description: "Flag. Root heuristics obey the heuristic LP-iteration quota also before the first incumbent (SCIP maxlpiterquot)",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_root_stall_scip",
        description: "Default-off lever. Root cut loop stops after 10 stalled separation rounds (SCIP separating/maxstallroundsroot) instead of 3",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_redcost_proof_uncapped",
        description: "Default-off lever. Reduced-cost fixing builds a dual proof for any number of bound changes instead of at most 100",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_presolve_probing_cap",
        description: "Default-off lever. MIP presolve stops probing after 10000 probed binaries instead of the HiGHS contingent",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_presolve_domcol_full",
        description: "Default-off lever. Dominated-column presolve never switches off its predictive bound analysis (HiGHS stops it after a check budget)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_presolve_probing_useless",
        description: "Presolve probing stops after this many successive probes with no deletion and no new clique (SCIP probing maxuseless); 0 = off; default-off lever",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 2147483647,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_conflict_off_mask",
        description: "Ablation lever. Bits: 1 tree and untagged (strong branching, sub-MIPs), 2 heur, 4/8 inf/bound-LP pool part only, 16 degen, 32 cut-pool cuts incl. heur calls",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 63,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: false,
    },
    OptionRecord {
        name: "mip_racing_threads",
        description: "Number of diversified copies of the MIP search raced concurrently, sharing incumbents (Shinano et al. 2018, Carvajal et al. 2014); 0 or 1 = off; default-off lever",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 64,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_root_lp_concurrent",
        description: "Default-off lever. Race the root LP between dual simplex and IPM with crossover on two threads (needs threads >= 2); the first optimal basis wins (Bixby 2002)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_conflict_resolution",
        description: "Generalized-resolution conflict analysis (SCIP 11.0.0 conflict_resolution.c) behind this default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "SCIP 11.0.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_conflict_resolution_tree_only",
        description: "Flag. With mip_conflict_resolution: skip the generalized-resolution step before the first node is counted (root dives, sub-MIP roots); default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_cut_selection_hybrid",
        description: "SCIP 11.0.0 hybrid cut selection (cutsel_hybrid.c + sepastore apply) behind this default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "SCIP 11.0.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_knapsack_cover_sepa",
        description: "SCIP linear-row knapsack cuts: lifted minimal covers and extended weight inequalities; default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_dff_separator",
        description: "Separate Fekete-Schepers dual-feasible-function rows sum f(a/b) x <= y of activation knapsacks sum a x <= b y at the violated LP point; default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_tableau_loosen",
        description: "Tableau separator without HiGHS weight-range 1e4, density and row-length caps, latched off after 2 stalled calls; 30n20b8 root 130 to 151, hurts neos-911970/abz6; default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_tableau_loosen_auto",
        description: "mip_tableau_loosen only when the first tableau round had >= 90% of its aggregations rejected by the density cap and made no cut; default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_ti_precedence_sepa",
        description: "Disaggregated (step) precedence cuts from aggregated sum-t-x precedence rows over exactly-one start-time sets (Sousa and Wolsey 1992); time-indexed job shop root bound",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_netdesign_sepa",
        description: "Strong linking cuts x <= min(d,u) y on capacity rows; each commodity is checked (cost >= 0, lower 0, network columns, monotone side rows) and skipped when the check fails",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_netdesign_heur",
        description: "Slope scaling primal heuristic for fixed-charge network design (Crainic, Gendron, Hernu 2004), once at the root",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_netdesign_cutset",
        description: "Cut-set + MIR cuts for network design, greedy node sets (Bienstock 1998, Atamturk 2002)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_skip_analytic_centre",
        description: "Skip the analytic centre IPM solve at the MIP root (Gondzio 1996 centre; HiGHS runs it on a spare thread); frees the single-thread root of huge models; default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_ti_list_sched",
        description: "Time-indexed scheduling primal heuristic: LP alpha-point list scheduling with hill climb (Kolisch 1996; Hall et al. 1997); first incumbent for time-indexed job shop",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_cont_prop_improve",
        description: "Relative range gain propagation needs to accept a continuous bound (HiGHS 0.3, SCIP 0.05); default 0.05, lower helps big-M scheduling",
        advanced: true,
        kind: RecordKind::Double {
            lower: 1e-6,
            default: 0.05,
            upper: 1.0,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_pscost_variance_tol",
        description: "Pseudocost reliability from sample variance (Hendel 2015): 95% CI half-width within this fraction of the mean; 0 = HiGHS count rule",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 0.0,
            upper: 10.0,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_cutoff_bigm",
        description: "Default-off lever. Restart when the incumbent cutoff lets presolve tighten big-M coefficients (Savelsbergh coefficient tightening on cutoff-propagated bounds)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_node_bigm_rows",
        description: "Default-off lever. With mip_node_bigm, add the tightened big-M row as a node-local LP row instead of changing the coefficient in place, so the basis survives",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_rens_min_fixing_rate",
        description: "RENS solves no sub-MIP when it fixed less than this share of the integer columns (SCIP RENS minfixingrate 0.5, Berthold 2014); 0 = off",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 0.0,
            upper: 1.0,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_node_bigm_min_ratio",
        description: "With mip_node_bigm, smallest coefficient reduction d relative to |c_z| worth applying",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 0.001,
            upper: 1.0,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_node_bigm_active",
        description: "With mip_node_bigm, only tighten rows whose slack in the previous LP solution is at most this multiple of |c_z|; 0 = all rows",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 0.0,
            upper: 1e9,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_node_bigm",
        description: "Default-off lever. Tighten big-M row coefficients of binary variables against the local node bounds in the node LP (Savelsbergh 1994)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_debug_solution_file",
        description: "Solution file for debug solution of the MIP solver",
        advanced: false,
        kind: RecordKind::Str { default: "" },
        source: "HiGHS 1.16.0",
        affects: Affects::Output,
        tunable: false,
    },
    OptionRecord {
        name: "mip_detect_symmetry",
        description: "Whether MIP symmetry should be detected",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_feasibility_tolerance",
        description: "MIP integrality tolerance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 1e-06,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_effort",
        description: "Effort spent for MIP heuristics",
        advanced: false,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 0.05,
            upper: 1.0,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_crossover",
        description: "SCIP 11.0.0 crossover primal heuristic (heur_crossover.c) behind this default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "SCIP 11.0.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_feasibility_jump",
        description: "Use the feasibility jump heuristic",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_fpr",
        description: "FPR presolve heuristic from mip-heuristics (Spoorendonk, MIT) behind this default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "mip-heuristics (MIT)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_local_mip",
        description: "Local-MIP presolve heuristic from mip-heuristics (Spoorendonk, MIT) behind this default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "mip-heuristics (MIT)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_oneopt",
        description: "SCIP 11.0.0 oneopt primal heuristic (heur_oneopt.c) behind this default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "SCIP 11.0.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_rens",
        description: "Use the RENS heuristic",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_rins",
        description: "Use the RINS heuristic",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_root_reduced_cost",
        description: "Use the rootReducedCost heuristic",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_shifting",
        description: "Use the Shifting heuristic",
        advanced: false,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_heuristic_run_zi_round",
        description: "Use the ZI Round heuristic",
        advanced: false,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_improving_solution_file",
        description: "File for reporting improving MIP solutions: not reported for an empty string \\\\\"\\\\\"",
        advanced: false,
        kind: RecordKind::Str { default: "" },
        source: "HiGHS 1.16.0",
        affects: Affects::Output,
        tunable: false,
    },
    OptionRecord {
        name: "mip_improving_solution_report_sparse",
        description: "Whether improving MIP solutions should be reported in sparse format",
        advanced: false,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_improving_solution_save",
        description: "Whether improving MIP solutions should be saved",
        advanced: false,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_ipm_solver",
        description: "MIP IPM solver: \"choose\", \"ipx\" or \"hipo\"",
        advanced: false,
        kind: RecordKind::Str { default: "choose" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_lifting_for_probing",
        description: "Level of lifting for probing that is used",
        advanced: false,
        kind: RecordKind::Int {
            lower: -1,
            default: -1,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_local_branching",
        description: "Default-off lever. Local branching (Fischetti-Lodi 2003): sub-MIP with the original objective inside a Hamming ball around the incumbent, radius adapted on the outcome",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_guided_dive",
        description: "Default-off lever. Guided diving (SCIP heur_guideddiving.c): LP dive that rounds fractional integers toward the incumbent, one backtrack, LP-iteration quota 5% + 1000",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_lp_age_limit",
        description: "Maximal age of dynamic LP rows before they are removed from the LP relaxation in the MIP solver",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 10,
            upper: 32767,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_lp_solver",
        description: "MIP LP solver: \"choose\", \"simplex\", \"ipm\", \"ipx\" or \"hipo\"",
        advanced: false,
        kind: RecordKind::Str { default: "choose" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_ls_gub_swap",
        description: "Flag. GUB swap / ejection-chain moves in Local-MIP (Benoist et al. 2011); acts only with mip_heuristic_run_local_mip; default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_ls_int_lp_repair",
        description: "Flag. Fix the integers of local-search end points and solve an LP over the rest to repair them; default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_max_improving_sols",
        description: "Limit on the number of improving solutions found to stop the MIP solver prematurely",
        advanced: false,
        kind: RecordKind::Int {
            lower: 1,
            default: 2147483647,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_max_leaves",
        description: "MIP solver max number of leaf nodes",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2147483647,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "mip_max_nodes",
        description: "MIP solver max number of nodes",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2147483647,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "mip_max_stall_nodes",
        description: "MIP solver max number of nodes where estimate is above cutoff bound",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2147483647,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "mip_max_start_nodes",
        description: "MIP solver max number of nodes when completing a partial MIP start",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 500,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_min_cliquetable_entries_for_parallelism",
        description: "Minimal number of entries in the MIP solver cliquetable before neighbourhood queries of the conflict graph use parallel processing",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 100000,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_min_logging_interval",
        description: "MIP minimum logging interval",
        advanced: false,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 5.0,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_pool_age_limit",
        description: "Maximal age of rows in the MIP solver cutpool before they are deleted",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 30,
            upper: 1000,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_pool_soft_limit",
        description: "Soft limit on the number of rows in the MIP solver cutpool for dynamic age adjustment",
        advanced: false,
        kind: RecordKind::Int {
            lower: 1,
            default: 10000,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_improve_early_work",
        description: "Work budget of the mip_improve_early chain as a multiple of the pre-root work; each step gets a fifth of it",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 20.0,
            upper: 1e6,
        },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_improve_early",
        description: "Default-off lever. Pre-root improvement chain: proximity-search sub-MIPs on the first (feasibility jump) incumbent before the root LP, work-capped",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_primal_stall",
        description: "Default-off lever. Stall-triggered proximity search: a proximity sub-MIP around the incumbent once neither the incumbent nor the dual bound has moved for a share of the search work",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_proximity",
        description: "Default-off lever. Proximity search (Fischetti-Monaci 2014): sub-MIP minimising the Hamming distance to the incumbent under an objective cutoff row, after a new incumbent has aged",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_pscost_minreliable",
        description: "Minimal number of observations before MIP solver pseudo costs are considered reliable",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 8,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_rel_gap",
        description: "Tolerance on relative gap, ub-lb/ub, to determine whether optimality has been reached for a MIP instance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 0.0001,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "mip_report_level",
        description: "MIP solver reporting level",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 1,
            upper: 2,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_root_presolve_only",
        description: "Whether MIP presolve is only applied at the root node",
        advanced: false,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_search_simulate_concurrency",
        description: "Simulate MIP search concurrency on a single thread",
        advanced: false,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_strengthen_cover_guard",
        description: "Flag. Exit the 1.16.0 strengthenInequalities cover loop when slackupper < 0 (HiGHS asserts it), where it can loop forever; on in recommended; lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "rusolve",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_symmetry_dynamic",
        description: "SCIP 11.0.0 dynamic symmetry handling (orbitopal reduction + LexRed, propagating/symmetry/usedynamicprop) behind this default-off lever",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "SCIP 11.0.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_aggregate_implications",
        description: "Ported behind this default-off lever: HiGHS #3562 #3596",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_cdouble_act_contrib",
        description: "Ported behind this default-off lever: HiGHS #3184, Upgrade activity contributions to HighsCDouble",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_clique_normalise",
        description: "Ported behind this default-off lever: HiGHS #3286, Enhance clique extraction (clique-row normalisation)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_coef_change_impl_bounds",
        description: "Ported behind this default-off lever: HiGHS #3318, implied-bound bookkeeping on coefficient change (matrixNonZeroChanged, link/unlink/addToMatrix)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_dfs_node_limit",
        description: "Ported behind this default-off lever: HiGHS #3445, Limit nodes in the RENS/RINS depth-first search",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_dual_fix_probing",
        description: "Ported behind this default-off lever: HiGHS #3244, Dual fixing inside probing",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_dual_side_lost",
        description: "Ported behind this default-off lever: HiGHS #3332, Reset row-dual implied bounds when the dual side is lost",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_dual_subst_limits",
        description: "Ported behind this default-off lever: HiGHS #3404, Early stop for dual substitution, row size limit",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_keep_branch_reason",
        description: "Ported behind this default-off lever: HiGHS #3357, Always retain branch as reason",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_lurking_bounds_robust",
        description: "Ported behind this default-off lever: HiGHS #3481 (+ part of #3474), robust lurking bounds",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_machine_sched_sepa",
        description: "Ported behind this default-off lever: HiGHS #3591",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0 (backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_path_same_bds",
        description: "Ported behind this default-off lever: HiGHS #3179, Same bound types in path-mixing transform",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_ranged_dual_implied",
        description: "Ported behind this default-off lever: HiGHS #3437, Substitution candidates when a ranged row becomes dual implied free",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_reduced_stack_fast",
        description: "Ported behind this default-off lever: HiGHS #3336, dive limit 10000",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_redundant_vbd",
        description: "Ported behind this default-off lever: HiGHS #3177, Handling of redundant variable bounds",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_sb_budget_pre_restart",
        description: "Ported behind this default-off lever: HiGHS #3213, Reduce strong branching budget before two restarts",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_strengthen_small_contrib",
        description: "Ported behind this default-off lever: HiGHS #3468, Small contributions in strengthenInequalities",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_cdouble_proofs",
        description: "Ported behind this default-off lever: HiGHS #3223, Fix traininstance6 bug (closes #2874)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_fix_last_active",
        description: "Ported behind this default-off lever: HiGHS #3220, fixLastActiveAndRemove",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_fm_presolve",
        description: "Ported behind this default-off lever: HiGHS #3255, Fourier-Motzkin presolve",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_keep_implications",
        description: "Ported behind this default-off lever: HiGHS #3209, Keep implications between restarts",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_presolve_clique_table",
        description: "Ported behind this default-off lever: HiGHS #3312, Clique table through presolve",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_presolve_row_append",
        description: "Ported behind this default-off lever: HiGHS #3218, FMPresolve5-peripheral (row/col append infrastructure)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "mip_sym_hash64",
        description: "Ported behind this default-off lever: HiGHS #3535, add-symmetry-timelimit (u64 vertex hashes, detection time limit)",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS (post-1.15.1 PR backport)",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "no_unnecessary_rebuild_refactor",
        description: "No unnecessary refactorization on simplex rebuild",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "objective_bound",
        description: "Objective bound for termination of the dual simplex solver",
        advanced: false,
        kind: RecordKind::Double {
            lower: -f64::INFINITY,
            default: f64::INFINITY,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "objective_target",
        description: "Objective target for termination of the MIP solver",
        advanced: false,
        kind: RecordKind::Double {
            lower: -f64::INFINITY,
            default: -f64::INFINITY,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "optimality_tolerance",
        description: "Optimality tolerance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 1e-07,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "output_flag",
        description: "Enables or disables solver output",
        advanced: false,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Output,
        tunable: false,
    },
    OptionRecord {
        name: "parallel",
        description: "Parallel: \"off\", \"choose\" or \"on\"",
        advanced: false,
        kind: RecordKind::Str { default: "choose" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "pdlp_cupdlpc_restart_method",
        description: "Restart mode for cuPDLP-C solver: 0 => none; 1 => GPU (default); 2 => CPU",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 1,
            upper: 2,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "pdlp_iteration_limit",
        description: "Iteration limit for PDLP solver",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2147483647,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "pdlp_optimality_tolerance",
        description: "PDLP optimality tolerance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 1e-07,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "pdlp_restart_strategy",
        description: "Restart strategy for PDLP solver: 0 => off; 1 => fixed; 2 => adaptive; 3 => Halpern",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2,
            upper: 3,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "pdlp_ruiz_iterations",
        description: "Number of Ruiz scaling iteraitons for PDLP solver",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 10,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "pdlp_scaling_mode",
        description: "Scaling mode for PDLP solver (default = 5): 1 => Ruiz; 2 => L2; 4 => PC",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 5,
            upper: 7,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "pdlp_step_size_strategy",
        description: "Step size strategy for PDLP solver: 0 => fixed; 1 => adaptive; 2 => Malitsky-Pock; 3 => PID",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 1,
            upper: 3,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve",
        description: "Presolve: \"off\", \"choose\" or \"on\"",
        advanced: false,
        kind: RecordKind::Str { default: "choose" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_fm_obj_reformulation",
        description: "Fourier-Motzkin objective reformulation",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_mode",
        description: "Presolve mode 0 => default; 1 => initial sweep; 2 => light; 3 => full",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 3,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_pivot_threshold",
        description: "Matrix factorization pivot threshold for substitutions in presolve",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0008,
            default: 0.01,
            upper: 0.5,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_reduction_limit",
        description: "Limit on number of presolve reductions -1 => no limit",
        advanced: true,
        kind: RecordKind::Int {
            lower: -1,
            default: -1,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_remove_slacks",
        description: "Remove slacks after presolve: redundant option!",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_rule_logging",
        description: "Log effectiveness of presolve rules for LP",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_rule_off",
        description: "Bit mask of presolve rules that are not allowed",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_rule_test",
        description: "Presolve rule to test - DEV only!",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 19,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "presolve_substitution_maxfillin",
        description: "Maximal fillin allowed for substitutions in presolve",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 10,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "primal_feasibility_tolerance",
        description: "Primal feasibility tolerance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 1e-07,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "primal_residual_tolerance",
        description: "Primal residual tolerance",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-10,
            default: 1e-07,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "primal_simplex_bound_perturbation_multiplier",
        description: "Primal simplex bound perturbation multiplier: 0 => no perturbation",
        advanced: true,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 1.0,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "qp_allow_hot_start",
        description: "Allow the active set QP solver to hot start",
        advanced: false,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "qp_iteration_limit",
        description: "Iteration limit for the active set QP solver",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2147483647,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "qp_nullspace_limit",
        description: "Nullspace limit for the active set QP solver",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 4000,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "qp_regularization_value",
        description: "Regularization value added to the Hessian in the active set QP solver",
        advanced: false,
        kind: RecordKind::Double {
            lower: 0.0,
            default: 1e-07,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "random_seed",
        description: "Random seed used in HiGHS",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "rebuild_refactor_solution_error_tolerance",
        description: "Tolerance on solution error when considering refactorization on simplex rebuild",
        advanced: true,
        kind: RecordKind::Double {
            lower: -f64::INFINITY,
            default: 1e-08,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "restart_presolve_reduction_limit",
        description: "Limit on number of further presolve reductions on restart in MIP solver -1 => no limit, otherwise, must be positive",
        advanced: true,
        kind: RecordKind::Int {
            lower: -1,
            default: -1,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "run_centring",
        description: "Perform centring steps or not",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "run_crossover",
        description: "Run IPM crossover: \"off\", \"choose\" or \"on\"",
        advanced: false,
        kind: RecordKind::Str { default: "on" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_dual_edge_weight_strategy",
        description: "Strategy for simplex dual edge weights: Choose / Dantzig / Devex / Steepest Edge (-1/0/1/2)",
        advanced: false,
        kind: RecordKind::Int {
            lower: -1,
            default: -1,
            upper: 2,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_dualize_strategy",
        description: "Strategy for dualizing before simplex",
        advanced: true,
        kind: RecordKind::Int {
            lower: -1,
            default: -1,
            upper: 1,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_iteration_limit",
        description: "Iteration limit for simplex solver when solving LPs, but not subproblems in the MIP solver",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2147483647,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_max_concurrency",
        description: "Maximum level of concurrency in parallel simplex",
        advanced: false,
        kind: RecordKind::Int {
            lower: 1,
            default: 8,
            upper: 8,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_min_concurrency",
        description: "Minimum level of concurrency in parallel simplex",
        advanced: true,
        kind: RecordKind::Int {
            lower: 1,
            default: 1,
            upper: 8,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_permute_strategy",
        description: "Strategy for permuting before simplex",
        advanced: true,
        kind: RecordKind::Int {
            lower: -1,
            default: -1,
            upper: 1,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_price_strategy",
        description: "Strategy for PRICE in simplex",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 3,
            upper: 3,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_primal_edge_weight_strategy",
        description: "Strategy for simplex primal edge weights: Choose / Dantzig / Devex / Steepest Edge (-1/0/1/2)",
        advanced: false,
        kind: RecordKind::Int {
            lower: -1,
            default: -1,
            upper: 2,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_scale_strategy",
        description: "Simplex scaling strategy: off / choose / equilibration (default) / forced equilibration / max value (0/1/2/3/4)",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 2,
            upper: 4,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_strategy",
        description: "Strategy for simplex solver 0 => Choose; 1 => Dual (serial); 2 => Dual (SIP); 3 => Dual (PAMI); 4 => Primal",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 1,
            upper: 4,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_unscaled_solution_strategy",
        description: "Strategy for solving unscaled LP in simplex",
        advanced: true,
        kind: RecordKind::Int {
            lower: 0,
            default: 1,
            upper: 2,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "simplex_update_limit",
        description: "Limit on the number of simplex UPDATE operations",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 5000,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "small_matrix_value",
        description: "Lower limit on matrix entries: values less than or equal to this will be treated as zero",
        advanced: false,
        kind: RecordKind::Double {
            lower: 1e-12,
            default: 1e-09,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "solve_relaxation",
        description: "Solve the relaxation of discrete model components",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "solver",
        description: "LP/QP solver: \"choose\", \"simplex\", \"ipm\", \"ipx\", \"hipo\", \"pdlp\", \"qpasm\" or \"hipdlp\",",
        advanced: false,
        kind: RecordKind::Str { default: "choose" },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "start_crossover_tolerance",
        description: "Tolerance to be satisfied before IPM crossover will start",
        advanced: true,
        kind: RecordKind::Double {
            lower: 1e-12,
            default: 1e-08,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Tolerance,
        tunable: true,
    },
    OptionRecord {
        name: "threads",
        description: "Maximum number of threads used by HiGHS (0: automatic)",
        advanced: false,
        kind: RecordKind::Int {
            lower: 0,
            default: 0,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Runtime,
        tunable: false,
    },
    OptionRecord {
        name: "time_limit",
        description: "Time limit (seconds)",
        advanced: false,
        kind: RecordKind::Double {
            lower: 0.0,
            default: f64::INFINITY,
            upper: f64::INFINITY,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "timeless_log",
        description: "Suppression of time-based data in logging",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "use_implied_bounds_from_presolve",
        description: "Use relaxed implied bounds from presolve: redundant option!",
        advanced: true,
        kind: RecordKind::Bool { default: false },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "use_warm_start",
        description: "Use any warm start that is available",
        advanced: true,
        kind: RecordKind::Bool { default: true },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "user_bound_scale",
        description: "Exponent of power-of-two bound scaling for model",
        advanced: false,
        kind: RecordKind::Int {
            lower: -2147483647,
            default: 0,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
    OptionRecord {
        name: "user_objective_scale",
        description: "Exponent of power-of-two objective scaling for model",
        advanced: false,
        kind: RecordKind::Int {
            lower: -2147483647,
            default: 0,
            upper: 2147483647,
        },
        source: "HiGHS 1.16.0",
        affects: Affects::Search,
        tunable: true,
    },
];

/// Fields with no record because nothing wires them by name: `hipo_available` (a build
/// capability HiGHS reads internally), `log_options` (derived from the log options) and
/// `pdlp_features_off` (a build-time feature switch). Kept explicit so the field-coverage
/// test stays exhaustive.
pub const INTERNAL_FIELDS: &[&str] = &["hipo_available", "log_options", "pdlp_features_off"];

fn affects_name(affects: Affects) -> &'static str {
    match affects {
        Affects::Search => "search",
        Affects::Tolerance => "tolerance",
        Affects::Semantics => "semantics",
        Affects::Output => "output",
        Affects::Runtime => "runtime",
    }
}

/// Whether changing this option can move the search (basis of the JSON `affects_search`).
pub fn affects_search(affects: Affects) -> bool {
    matches!(affects, Affects::Search)
}

fn record_for(name: &str) -> Option<&'static OptionRecord> {
    OPTION_RECORDS.iter().find(|record| record.name == name)
}

/// Number of registered options (used by the CLI's listing header).
pub fn option_record_count() -> usize {
    OPTION_RECORDS.len()
}

/// Registry type name for a registered option (HiGHS `getOptionType` analogue); `None` is
/// HiGHS `OptionStatus::kUnknownOption`.
pub fn option_type_name(name: &str) -> Option<&'static str> {
    let kind = &record_for(name)?.kind;
    Some(match kind {
        RecordKind::Bool { .. } => "bool",
        RecordKind::Int { .. } => "HighsInt",
        RecordKind::Double { .. } => "double",
        RecordKind::Str { .. } => "string",
    })
}

/// The registry's declared range for an int option.
pub fn option_int_bounds(name: &str) -> Option<(i32, i32)> {
    match &record_for(name)?.kind {
        RecordKind::Int { lower, upper, .. } => Some((*lower, *upper)),
        _ => None,
    }
}

/// The registry's declared range for a double option.
pub fn option_double_bounds(name: &str) -> Option<(f64, f64)> {
    match &record_for(name)?.kind {
        RecordKind::Double { lower, upper, .. } => Some((*lower, *upper)),
        _ => None,
    }
}

/// The record's declared default, used by `--non-default` listings and the drift test.
pub fn option_default(name: &str) -> Option<OptionValue> {
    match &record_for(name)?.kind {
        RecordKind::Bool { default } => Some(OptionValue::Bool(*default)),
        RecordKind::Int { default, .. } => Some(OptionValue::Int(*default)),
        RecordKind::Double { default, .. } => Some(OptionValue::Double(*default)),
        RecordKind::Str { default } => Some(OptionValue::Str((*default).to_string())),
    }
}

/// Generated by `option_field_accessors!`: variant-preserving get/set of a struct field, so
/// the typed setters in `highs_api_mip.rs` and the listing share one name→field map.
macro_rules! option_field_accessors {
    ($($name:literal => $field:ident : $variant:path;)*) => {
        pub(crate) fn get_option_field(opts: &HighsOptions, name: &str) -> Option<OptionValue> {
            match name {
                $($name => Some($variant(opts.$field.clone())),)*
                _ => None,
            }
        }

        pub(crate) fn set_option_field(
            opts: &mut HighsOptions,
            name: &str,
            value: OptionValue,
        ) -> Result<(), String> {
            match name {
                $($name => match value {
                    $variant(v) => {
                        opts.$field = v;
                        Ok(())
                    }
                    _ => Err(format!("internal: option '{}' changed type", name)),
                },)*
                _ => Err(format!("internal: unregistered option '{}'", name)),
            }
        }
    };
}

option_field_accessors! {
    "allow_unbounded_or_infeasible" => allow_unbounded_or_infeasible : OptionValue::Bool;
    "allowed_cost_scale_factor" => allowed_cost_scale_factor : OptionValue::Int;
    "allowed_matrix_scale_factor" => allowed_matrix_scale_factor : OptionValue::Int;
    "blend_multi_objectives" => blend_multi_objectives : OptionValue::Bool;
    "centring_ratio_tolerance" => centring_ratio_tolerance : OptionValue::Double;
    "cost_scale_factor" => cost_scale_factor : OptionValue::Int;
    "dual_feasibility_tolerance" => dual_feasibility_tolerance : OptionValue::Double;
    "dual_residual_tolerance" => dual_residual_tolerance : OptionValue::Double;
    "dual_simplex_cost_perturbation_multiplier" => dual_simplex_cost_perturbation_multiplier : OptionValue::Double;
    "dual_simplex_pivot_growth_tolerance" => dual_simplex_pivot_growth_tolerance : OptionValue::Double;
    "dual_steepest_edge_weight_error_tolerance" => dual_steepest_edge_weight_error_tolerance : OptionValue::Double;
    "dual_steepest_edge_weight_log_error_threshold" => dual_steepest_edge_weight_log_error_threshold : OptionValue::Double;
    "factor_pivot_threshold" => factor_pivot_threshold : OptionValue::Double;
    "factor_pivot_tolerance" => factor_pivot_tolerance : OptionValue::Double;
    "h16_ac_centring" => h16_ac_centring : OptionValue::Bool;
    "h16_ac_fix" => h16_ac_fix : OptionValue::Bool;
    "h16_add_no_modify_on_error" => h16_add_no_modify_on_error : OptionValue::Bool;
    "h16_agg_term_zero" => h16_agg_term_zero : OptionValue::Bool;
    "h16_assess_small" => h16_assess_small : OptionValue::Bool;
    "h16_big_act" => h16_big_act : OptionValue::Bool;
    "h16_bound_chain" => h16_bound_chain : OptionValue::Bool;
    "h16_bound_range_inf" => h16_bound_range_inf : OptionValue::Bool;
    "h16_cdouble" => h16_cdouble : OptionValue::Bool;
    "h16_clique_links" => h16_clique_links : OptionValue::Bool;
    "h16_col_stuffing" => h16_col_stuffing : OptionValue::Bool;
    "h16_cost_perturb" => h16_cost_perturb : OptionValue::Bool;
    "h16_crossover_x" => h16_crossover_x : OptionValue::Bool;
    "h16_cupdlp_time" => h16_cupdlp_time : OptionValue::Bool;
    "h16_cutpool_age" => h16_cutpool_age : OptionValue::Bool;
    "h16_dep_eq_time_cap" => h16_dep_eq_time_cap : OptionValue::Bool;
    "h16_dominated_cols" => h16_dominated_cols : OptionValue::Bool;
    "h16_dual_free_list" => h16_dual_free_list : OptionValue::Bool;
    "h16_dual_postsolve_rules" => h16_dual_postsolve_rules : OptionValue::Bool;
    "h16_dup_row_undo" => h16_dup_row_undo : OptionValue::Bool;
    "h16_dupcol_warn" => h16_dupcol_warn : OptionValue::Bool;
    "h16_enum_sig64" => h16_enum_sig64 : OptionValue::Bool;
    "h16_eq_row_adds" => h16_eq_row_adds : OptionValue::Bool;
    "h16_feasrelax_sense" => h16_feasrelax_sense : OptionValue::Bool;
    "h16_fj" => h16_fj : OptionValue::Bool;
    "h16_fourier_motzkin" => h16_fourier_motzkin : OptionValue::Bool;
    "h16_frontier_merge" => h16_frontier_merge : OptionValue::Bool;
    "h16_hash_common" => h16_hash_common : OptionValue::Bool;
    "h16_hfactor_pivot" => h16_hfactor_pivot : OptionValue::Bool;
    "h16_hipdlp_max" => h16_hipdlp_max : OptionValue::Bool;
    "h16_hipo_display" => h16_hipo_display : OptionValue::Bool;
    "h16_hipo_empty_rows" => h16_hipo_empty_rows : OptionValue::Bool;
    "h16_impl_eq" => h16_impl_eq : OptionValue::Bool;
    "h16_implint_bounds" => h16_implint_bounds : OptionValue::Bool;
    "h16_int_bound_round" => h16_int_bound_round : OptionValue::Bool;
    "h16_int_bound_rounding" => h16_int_bound_rounding : OptionValue::Bool;
    "h16_int_round_assert_tol" => h16_int_round_assert_tol : OptionValue::Bool;
    "h16_ipm_gap_offset" => h16_ipm_gap_offset : OptionValue::Bool;
    "h16_ipx_free_rows" => h16_ipx_free_rows : OptionValue::Bool;
    "h16_ipx_scaled" => h16_ipx_scaled : OptionValue::Bool;
    "h16_ipx_stall" => h16_ipx_stall : OptionValue::Bool;
    "h16_kkt_effective_costs" => h16_kkt_effective_costs : OptionValue::Bool;
    "h16_max_duals" => h16_max_duals : OptionValue::Bool;
    "h16_mip_start" => h16_mip_start : OptionValue::Bool;
    "h16_mip_workers" => h16_mip_workers : OptionValue::Bool;
    "h16_ninfmin" => h16_ninfmin : OptionValue::Bool;
    "h16_nla_debug_gate" => h16_nla_debug_gate : OptionValue::Bool;
    "h16_obj_cdouble" => h16_obj_cdouble : OptionValue::Bool;
    "h16_parallel_rows" => h16_parallel_rows : OptionValue::Bool;
    "h16_postsolve_basis_check" => h16_postsolve_basis_check : OptionValue::Bool;
    "h16_presolve_rules" => h16_presolve_rules : OptionValue::Bool;
    "h16_qp_bounded_start" => h16_qp_bounded_start : OptionValue::Bool;
    "h16_qp_convexity" => h16_qp_convexity : OptionValue::Bool;
    "h16_ray_simplex" => h16_ray_simplex : OptionValue::Bool;
    "h16_rebuild_reason" => h16_rebuild_reason : OptionValue::Bool;
    "h16_repair_lp" => h16_repair_lp : OptionValue::Bool;
    "h16_row_duals" => h16_row_duals : OptionValue::Bool;
    "h16_row_presolve" => h16_row_presolve : OptionValue::Bool;
    "h16_row_redundant" => h16_row_redundant : OptionValue::Bool;
    "h16_scaling_report" => h16_scaling_report : OptionValue::Bool;
    "h16_singleton_loops" => h16_singleton_loops : OptionValue::Bool;
    "h16_sub_solver_flags" => h16_sub_solver_flags : OptionValue::Bool;
    "h16_subst_int" => h16_subst_int : OptionValue::Bool;
    "h16_sym_refine" => h16_sym_refine : OptionValue::Bool;
    "h16_unresolved_unbounded_or_infeasible" => h16_unresolved_unbounded_or_infeasible : OptionValue::Bool;
    "h16_user_scale_offset" => h16_user_scale_offset : OptionValue::Bool;
    "h16_user_sol_obj" => h16_user_sol_obj : OptionValue::Bool;
    "h16_user_solution_assess" => h16_user_solution_assess : OptionValue::Bool;
    "h16_watch_link" => h16_watch_link : OptionValue::Bool;
    "h16_zero_cost_singleton" => h16_zero_cost_singleton : OptionValue::Bool;
    "highs_analysis_level" => highs_analysis_level : OptionValue::Int;
    "highs_debug_level" => highs_debug_level : OptionValue::Int;
    "hipo_block_size" => hipo_block_size : OptionValue::Int;
    "hipo_ordering" => hipo_ordering : OptionValue::Str;
    "hipo_parallel_type" => hipo_parallel_type : OptionValue::Str;
    "hipo_system" => hipo_system : OptionValue::Str;
    "icrash" => icrash : OptionValue::Bool;
    "icrash_approx_iter" => icrash_approx_iter : OptionValue::Int;
    "icrash_breakpoints" => icrash_breakpoints : OptionValue::Bool;
    "icrash_dualize" => icrash_dualize : OptionValue::Bool;
    "icrash_exact" => icrash_exact : OptionValue::Bool;
    "icrash_iterations" => icrash_iterations : OptionValue::Int;
    "icrash_starting_weight" => icrash_starting_weight : OptionValue::Double;
    "icrash_strategy" => icrash_strategy : OptionValue::Str;
    "iis_strategy" => iis_strategy : OptionValue::Int;
    "iis_time_limit" => iis_time_limit : OptionValue::Double;
    "infinite_bound" => infinite_bound : OptionValue::Double;
    "infinite_cost" => infinite_cost : OptionValue::Double;
    "ipm_iteration_limit" => ipm_iteration_limit : OptionValue::Int;
    "ipm_optimality_tolerance" => ipm_optimality_tolerance : OptionValue::Double;
    "ipx_dualize_strategy" => ipx_dualize_strategy : OptionValue::Int;
    "kkt_tolerance" => kkt_tolerance : OptionValue::Double;
    "large_matrix_value" => large_matrix_value : OptionValue::Double;
    "less_infeasible_DSE_check" => less_infeasible_DSE_check : OptionValue::Bool;
    "less_infeasible_DSE_choose_row" => less_infeasible_DSE_choose_row : OptionValue::Bool;
    "log_dev_level" => log_dev_level : OptionValue::Int;
    "log_file" => log_file : OptionValue::Str;
    "log_to_console" => log_to_console : OptionValue::Bool;
    "lp_fold" => lp_fold : OptionValue::Bool;
    "lp_crossover_ipm_crash" => lp_crossover_ipm_crash : OptionValue::Bool;
    "lp_crossover_batch_push" => lp_crossover_batch_push : OptionValue::Bool;
    "lp_dual_devex_start" => lp_dual_devex_start : OptionValue::Bool;
    "lp_postsolve_snap_basis" => lp_postsolve_snap_basis : OptionValue::Bool;
    "lp_ipm_support_reduction" => lp_ipm_support_reduction : OptionValue::Bool;
    "lp_concurrent" => lp_concurrent : OptionValue::Bool;
    "lp_solver_auto" => lp_solver_auto : OptionValue::Bool;
    "lp_postsolve_fix_tiny" => lp_postsolve_fix_tiny : OptionValue::Bool;
    "lp_fold_basis" => lp_fold_basis : OptionValue::Bool;
    "lp_crossover_ipx_refine" => lp_crossover_ipx_refine : OptionValue::Bool;
    "lp_ipx_no_drop" => lp_ipx_no_drop : OptionValue::Bool;
    "lp_ipx_soft_drop" => lp_ipx_soft_drop : OptionValue::Bool;
    "presolve_dep_eq_fixed_budget" => presolve_dep_eq_fixed_budget : OptionValue::Bool;
    "lp_dualize_auto" => lp_dualize_auto : OptionValue::Bool;
    "ipm_hipo_fast_blas" => ipm_hipo_fast_blas : OptionValue::Bool;
    "ipm_hipo_cache_perm" => ipm_hipo_cache_perm : OptionValue::Bool;
    "ipm_hipo_par_ne" => ipm_hipo_par_ne : OptionValue::Bool;
    "ipm_hipo_stall_detect" => ipm_hipo_stall_detect : OptionValue::Bool;
    "ipm_hipo_amd_first" => ipm_hipo_amd_first : OptionValue::Bool;
    "ipm_hipo_dense_cols" => ipm_hipo_dense_cols : OptionValue::Bool;
    "ipm_hipo_dense_schur" => ipm_hipo_dense_schur : OptionValue::Bool;
    "ipm_hipo_dualize" => ipm_hipo_dualize : OptionValue::Bool;
    "lp_ipm_failure_policy" => lp_ipm_failure_policy : OptionValue::Bool;
    "lp_presolve_requires_basis_postsolve" => lp_presolve_requires_basis_postsolve : OptionValue::Bool;
    "max_centring_steps" => max_centring_steps : OptionValue::Int;
    "max_dual_simplex_cleanup_level" => max_dual_simplex_cleanup_level : OptionValue::Int;
    "max_dual_simplex_phase1_cleanup_level" => max_dual_simplex_phase1_cleanup_level : OptionValue::Int;
    "mip_abs_gap" => mip_abs_gap : OptionValue::Double;
    "mip_allow_cut_separation_at_nodes" => mip_allow_cut_separation_at_nodes : OptionValue::Bool;
    "mip_allow_restart" => mip_allow_restart : OptionValue::Bool;
    "mip_alns" => mip_alns : OptionValue::Bool;
    "mip_presolve_eq_subst" => mip_presolve_eq_subst : OptionValue::Bool;
    "mip_presolve_clique_varbound" => mip_presolve_clique_varbound : OptionValue::Bool;
    "mip_presolve_xor_gf2" => mip_presolve_xor_gf2 : OptionValue::Bool;
    "mip_orbitope_root_fix" => mip_orbitope_root_fix : OptionValue::Bool;
    "mip_orbitope_full_fix" => mip_orbitope_full_fix : OptionValue::Bool;
    "mip_shift_and_propagate" => mip_shift_and_propagate : OptionValue::Bool;
    "mip_cdcl_dive" => mip_cdcl_dive : OptionValue::Bool;
    "mip_root_bound_stall" => mip_root_bound_stall : OptionValue::Bool;
    "mip_cdcl_uip" => mip_cdcl_uip : OptionValue::Bool;
    "mip_zero_obj_cdcl" => mip_zero_obj_cdcl : OptionValue::Bool;
    "mip_decompose" => mip_decompose : OptionValue::Bool;
    "mip_lattice" => mip_lattice : OptionValue::Bool;
    "mip_activation_dff" => mip_activation_dff : OptionValue::Bool;
    "mip_fj_native" => mip_fj_native : OptionValue::Bool;
    "mip_first_inc" => mip_first_inc : OptionValue::Bool;
    "mip_first_inc_root" => mip_first_inc_root : OptionValue::Bool;
    "mip_lattice_kernel" => mip_lattice_kernel : OptionValue::Bool;
    "mip_presolve_probing_cap" => mip_presolve_probing_cap : OptionValue::Bool;
    "mip_presolve_domcol_full" => mip_presolve_domcol_full : OptionValue::Bool;
    "mip_redcost_proof_uncapped" => mip_redcost_proof_uncapped : OptionValue::Bool;
    "mip_root_stall_scip" => mip_root_stall_scip : OptionValue::Bool;
    "mip_heur_lp_quota" => mip_heur_lp_quota : OptionValue::Bool;
    "mip_racing_threads" => mip_racing_threads : OptionValue::Int;
    "mip_conflict_off_mask" => mip_conflict_off_mask : OptionValue::Int;
    "mip_presolve_probing_useless" => mip_presolve_probing_useless : OptionValue::Int;
    "mip_root_lp_concurrent" => mip_root_lp_concurrent : OptionValue::Bool;
    "mip_conflict_resolution" => mip_conflict_resolution : OptionValue::Bool;
    "mip_conflict_resolution_tree_only" => mip_conflict_resolution_tree_only : OptionValue::Bool;
    "mip_cut_selection_hybrid" => mip_cut_selection_hybrid : OptionValue::Bool;
    "mip_knapsack_cover_sepa" => mip_knapsack_cover_sepa : OptionValue::Bool;
    "mip_dff_separator" => mip_dff_separator : OptionValue::Bool;
    "mip_tableau_loosen" => mip_tableau_loosen : OptionValue::Bool;
    "mip_tableau_loosen_auto" => mip_tableau_loosen_auto : OptionValue::Bool;
    "mip_ti_precedence_sepa" => mip_ti_precedence_sepa : OptionValue::Bool;
    "mip_skip_analytic_centre" => mip_skip_analytic_centre : OptionValue::Bool;
    "mip_netdesign_sepa" => mip_netdesign_sepa : OptionValue::Bool;
    "mip_netdesign_heur" => mip_netdesign_heur : OptionValue::Bool;
    "mip_netdesign_cutset" => mip_netdesign_cutset : OptionValue::Bool;
    "mip_ti_list_sched" => mip_ti_list_sched : OptionValue::Bool;
    "mip_cont_prop_improve" => mip_cont_prop_improve : OptionValue::Double;
    "mip_pscost_variance_tol" => mip_pscost_variance_tol : OptionValue::Double;
    "mip_cutoff_bigm" => mip_cutoff_bigm : OptionValue::Bool;
    "mip_node_bigm" => mip_node_bigm : OptionValue::Bool;
    "mip_node_bigm_rows" => mip_node_bigm_rows : OptionValue::Bool;
    "mip_node_bigm_min_ratio" => mip_node_bigm_min_ratio : OptionValue::Double;
    "mip_rens_min_fixing_rate" => mip_rens_min_fixing_rate : OptionValue::Double;
    "mip_node_bigm_active" => mip_node_bigm_active : OptionValue::Double;
    "mip_debug_solution_file" => mip_debug_solution_file : OptionValue::Str;
    "mip_detect_symmetry" => mip_detect_symmetry : OptionValue::Bool;
    "mip_feasibility_tolerance" => mip_feasibility_tolerance : OptionValue::Double;
    "mip_heuristic_effort" => mip_heuristic_effort : OptionValue::Double;
    "mip_heuristic_run_crossover" => mip_heuristic_run_crossover : OptionValue::Bool;
    "mip_heuristic_run_feasibility_jump" => mip_heuristic_run_feasibility_jump : OptionValue::Bool;
    "mip_heuristic_run_fpr" => mip_heuristic_run_fpr : OptionValue::Bool;
    "mip_heuristic_run_local_mip" => mip_heuristic_run_local_mip : OptionValue::Bool;
    "mip_heuristic_run_oneopt" => mip_heuristic_run_oneopt : OptionValue::Bool;
    "mip_heuristic_run_rens" => mip_heuristic_run_rens : OptionValue::Bool;
    "mip_heuristic_run_rins" => mip_heuristic_run_rins : OptionValue::Bool;
    "mip_heuristic_run_root_reduced_cost" => mip_heuristic_run_root_reduced_cost : OptionValue::Bool;
    "mip_heuristic_run_shifting" => mip_heuristic_run_shifting : OptionValue::Bool;
    "mip_heuristic_run_zi_round" => mip_heuristic_run_zi_round : OptionValue::Bool;
    "mip_improving_solution_file" => mip_improving_solution_file : OptionValue::Str;
    "mip_improving_solution_report_sparse" => mip_improving_solution_report_sparse : OptionValue::Bool;
    "mip_improving_solution_save" => mip_improving_solution_save : OptionValue::Bool;
    "mip_ipm_solver" => mip_ipm_solver : OptionValue::Str;
    "mip_lifting_for_probing" => mip_lifting_for_probing : OptionValue::Int;
    "mip_guided_dive" => mip_guided_dive : OptionValue::Bool;
    "mip_local_branching" => mip_local_branching : OptionValue::Bool;
    "mip_lp_age_limit" => mip_lp_age_limit : OptionValue::Int;
    "mip_lp_solver" => mip_lp_solver : OptionValue::Str;
    "mip_ls_gub_swap" => mip_ls_gub_swap : OptionValue::Bool;
    "mip_ls_int_lp_repair" => mip_ls_int_lp_repair : OptionValue::Bool;
    "mip_max_improving_sols" => mip_max_improving_sols : OptionValue::Int;
    "mip_max_leaves" => mip_max_leaves : OptionValue::Int;
    "mip_max_nodes" => mip_max_nodes : OptionValue::Int;
    "mip_max_stall_nodes" => mip_max_stall_nodes : OptionValue::Int;
    "mip_max_start_nodes" => mip_max_start_nodes : OptionValue::Int;
    "mip_min_cliquetable_entries_for_parallelism" => mip_min_cliquetable_entries_for_parallelism : OptionValue::Int;
    "mip_min_logging_interval" => mip_min_logging_interval : OptionValue::Double;
    "mip_pool_age_limit" => mip_pool_age_limit : OptionValue::Int;
    "mip_pool_soft_limit" => mip_pool_soft_limit : OptionValue::Int;
    "mip_proximity" => mip_proximity : OptionValue::Bool;
    "mip_improve_early" => mip_improve_early : OptionValue::Bool;
    "mip_primal_stall" => mip_primal_stall : OptionValue::Bool;
    "mip_improve_early_work" => mip_improve_early_work : OptionValue::Double;
    "mip_pscost_minreliable" => mip_pscost_minreliable : OptionValue::Int;
    "mip_rel_gap" => mip_rel_gap : OptionValue::Double;
    "mip_report_level" => mip_report_level : OptionValue::Int;
    "mip_root_presolve_only" => mip_root_presolve_only : OptionValue::Bool;
    "mip_search_simulate_concurrency" => mip_search_simulate_concurrency : OptionValue::Bool;
    "mip_strengthen_cover_guard" => mip_strengthen_cover_guard : OptionValue::Bool;
    "mip_symmetry_dynamic" => mip_symmetry_dynamic : OptionValue::Bool;
    "mip_aggregate_implications" => mip_aggregate_implications : OptionValue::Bool;
    "mip_cdouble_act_contrib" => mip_cdouble_act_contrib : OptionValue::Bool;
    "mip_clique_normalise" => mip_clique_normalise : OptionValue::Bool;
    "mip_coef_change_impl_bounds" => mip_coef_change_impl_bounds : OptionValue::Bool;
    "mip_dfs_node_limit" => mip_dfs_node_limit : OptionValue::Bool;
    "mip_dual_fix_probing" => mip_dual_fix_probing : OptionValue::Bool;
    "mip_dual_side_lost" => mip_dual_side_lost : OptionValue::Bool;
    "mip_dual_subst_limits" => mip_dual_subst_limits : OptionValue::Bool;
    "mip_keep_branch_reason" => mip_keep_branch_reason : OptionValue::Bool;
    "mip_lurking_bounds_robust" => mip_lurking_bounds_robust : OptionValue::Bool;
    "mip_machine_sched_sepa" => mip_machine_sched_sepa : OptionValue::Bool;
    "mip_path_same_bds" => mip_path_same_bds : OptionValue::Bool;
    "mip_ranged_dual_implied" => mip_ranged_dual_implied : OptionValue::Bool;
    "mip_reduced_stack_fast" => mip_reduced_stack_fast : OptionValue::Bool;
    "mip_redundant_vbd" => mip_redundant_vbd : OptionValue::Bool;
    "mip_sb_budget_pre_restart" => mip_sb_budget_pre_restart : OptionValue::Bool;
    "mip_strengthen_small_contrib" => mip_strengthen_small_contrib : OptionValue::Bool;
    "mip_cdouble_proofs" => mip_cdouble_proofs : OptionValue::Bool;
    "mip_fix_last_active" => mip_fix_last_active : OptionValue::Bool;
    "mip_fm_presolve" => mip_fm_presolve : OptionValue::Bool;
    "mip_keep_implications" => mip_keep_implications : OptionValue::Bool;
    "mip_presolve_clique_table" => mip_presolve_clique_table : OptionValue::Bool;
    "mip_presolve_row_append" => mip_presolve_row_append : OptionValue::Bool;
    "mip_sym_hash64" => mip_sym_hash64 : OptionValue::Bool;
    "no_unnecessary_rebuild_refactor" => no_unnecessary_rebuild_refactor : OptionValue::Bool;
    "objective_bound" => objective_bound : OptionValue::Double;
    "objective_target" => objective_target : OptionValue::Double;
    "optimality_tolerance" => optimality_tolerance : OptionValue::Double;
    "output_flag" => output_flag : OptionValue::Bool;
    "parallel" => parallel : OptionValue::Str;
    "pdlp_cupdlpc_restart_method" => pdlp_cupdlpc_restart_method : OptionValue::Int;
    "pdlp_iteration_limit" => pdlp_iteration_limit : OptionValue::Int;
    "pdlp_optimality_tolerance" => pdlp_optimality_tolerance : OptionValue::Double;
    "pdlp_restart_strategy" => pdlp_restart_strategy : OptionValue::Int;
    "pdlp_ruiz_iterations" => pdlp_ruiz_iterations : OptionValue::Int;
    "pdlp_scaling_mode" => pdlp_scaling_mode : OptionValue::Int;
    "pdlp_step_size_strategy" => pdlp_step_size_strategy : OptionValue::Int;
    "presolve" => presolve : OptionValue::Str;
    "presolve_fm_obj_reformulation" => presolve_fm_obj_reformulation : OptionValue::Bool;
    "presolve_mode" => presolve_mode : OptionValue::Int;
    "presolve_pivot_threshold" => presolve_pivot_threshold : OptionValue::Double;
    "presolve_reduction_limit" => presolve_reduction_limit : OptionValue::Int;
    "presolve_remove_slacks" => presolve_remove_slacks : OptionValue::Bool;
    "presolve_rule_logging" => presolve_rule_logging : OptionValue::Bool;
    "presolve_rule_off" => presolve_rule_off : OptionValue::Int;
    "presolve_rule_test" => presolve_rule_test : OptionValue::Int;
    "presolve_substitution_maxfillin" => presolve_substitution_maxfillin : OptionValue::Int;
    "primal_feasibility_tolerance" => primal_feasibility_tolerance : OptionValue::Double;
    "primal_residual_tolerance" => primal_residual_tolerance : OptionValue::Double;
    "primal_simplex_bound_perturbation_multiplier" => primal_simplex_bound_perturbation_multiplier : OptionValue::Double;
    "qp_allow_hot_start" => qp_allow_hot_start : OptionValue::Bool;
    "qp_iteration_limit" => qp_iteration_limit : OptionValue::Int;
    "qp_nullspace_limit" => qp_nullspace_limit : OptionValue::Int;
    "qp_regularization_value" => qp_regularization_value : OptionValue::Double;
    "random_seed" => random_seed : OptionValue::Int;
    "rebuild_refactor_solution_error_tolerance" => rebuild_refactor_solution_error_tolerance : OptionValue::Double;
    "restart_presolve_reduction_limit" => restart_presolve_reduction_limit : OptionValue::Int;
    "run_centring" => run_centring : OptionValue::Bool;
    "run_crossover" => run_crossover : OptionValue::Str;
    "simplex_dual_edge_weight_strategy" => simplex_dual_edge_weight_strategy : OptionValue::Int;
    "simplex_dualize_strategy" => simplex_dualize_strategy : OptionValue::Int;
    "simplex_iteration_limit" => simplex_iteration_limit : OptionValue::Int;
    "simplex_max_concurrency" => simplex_max_concurrency : OptionValue::Int;
    "simplex_min_concurrency" => simplex_min_concurrency : OptionValue::Int;
    "simplex_permute_strategy" => simplex_permute_strategy : OptionValue::Int;
    "simplex_price_strategy" => simplex_price_strategy : OptionValue::Int;
    "simplex_primal_edge_weight_strategy" => simplex_primal_edge_weight_strategy : OptionValue::Int;
    "simplex_scale_strategy" => simplex_scale_strategy : OptionValue::Int;
    "simplex_strategy" => simplex_strategy : OptionValue::Int;
    "simplex_unscaled_solution_strategy" => simplex_unscaled_solution_strategy : OptionValue::Int;
    "simplex_update_limit" => simplex_update_limit : OptionValue::Int;
    "small_matrix_value" => small_matrix_value : OptionValue::Double;
    "solve_relaxation" => solve_relaxation : OptionValue::Bool;
    "solver" => solver : OptionValue::Str;
    "start_crossover_tolerance" => start_crossover_tolerance : OptionValue::Double;
    "threads" => threads : OptionValue::Int;
    "time_limit" => time_limit : OptionValue::Double;
    "timeless_log" => timeless_log : OptionValue::Bool;
    "use_implied_bounds_from_presolve" => use_implied_bounds_from_presolve : OptionValue::Bool;
    "use_warm_start" => use_warm_start : OptionValue::Bool;
    "user_bound_scale" => user_bound_scale : OptionValue::Int;
    "user_objective_scale" => user_objective_scale : OptionValue::Int;
}

/// Current value of a registered option by name.
pub fn get_option_value(opts: &HighsOptions, name: &str) -> Option<OptionValue> {
    get_option_field(opts, name)
}

// ---------------------------------------------------------------------------
// Actionable unknown-option errors: closest registry names.
// ---------------------------------------------------------------------------

fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0usize; b.len() + 1];
    for (i, ca) in a.chars().enumerate() {
        cur[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let substitution = prev[j] + usize::from(ca != *cb);
            cur[j + 1] = substitution.min(prev[j + 1] + 1).min(cur[j] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

/// Up to three closest registry names, closest first. A candidate only qualifies when it is
/// genuinely close (within a quarter of the longer name) and not the name itself, so an
/// unrelated word gets no bogus suggestion and a real name is never "corrected".
pub fn nearest_option_names(name: &str) -> Vec<&'static str> {
    let mut scored: Vec<(usize, usize, &str)> = OPTION_RECORDS
        .iter()
        .map(|record| {
            (
                edit_distance(name, record.name),
                record.name.len(),
                record.name,
            )
        })
        .collect();
    scored.sort();
    let candidates: Vec<&str> = scored
        .into_iter()
        .filter(|(distance, _len, _candidate)| *distance > 0)
        .filter(|(distance, _len, candidate)| *distance <= name.len().max(candidate.len()) / 4)
        .map(|(_, _, candidate)| candidate)
        .collect();
    candidates.iter().take(3).copied().collect()
}

/// The A7 hint for an unrecognised name: `Some("did you mean X?")` when one name is clearly
/// closest, `Some("closest options are X, Y, Z")` when several are equally near.
pub fn unknown_option_hint(name: &str) -> Option<String> {
    let candidates = nearest_option_names(name);
    match candidates.as_slice() {
        [] => None,
        [single] => Some(format!("did you mean {}?", single)),
        many => {
            let distance = edit_distance(name, many[0]);
            let tied = many[1..]
                .iter()
                .any(|other| edit_distance(name, other) == distance);
            if tied {
                Some(format!("closest options are {}", many.join(", ")))
            } else {
                Some(format!("did you mean {}?", many[0]))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Options files (HiGHS 1.16.0 io/LoadOptions.cpp loadOptionsFromFile +
// lp_data/HighsOptions.cpp reportOptions/writeOptionsToFile, kFull style).
// ---------------------------------------------------------------------------

/// One `name = value` line of an options file, with the line it came from so errors can point
/// at it.
#[derive(Debug, Clone, PartialEq)]
pub struct OptionsFileEntry {
    pub line: usize,
    pub name: String,
    pub value: String,
}

/// Read an options file with HiGHS `loadOptionsFromFile` semantics: blank lines and lines
/// starting with '#' in column 0 are skipped; every other line must have a usable '=' with a
/// non-blank name and value; names and values are trimmed of whitespace and quote characters.
///
/// Two documented deviations from HiGHS, both for actionable errors: an unusable line is
/// reported by number, and an unknown name is reported with its closest registry names, so an
/// agent can fix the file instead of guessing. Applying the entries is the caller's job,
/// through the same typed setter `--opt` uses, so validation cannot differ between the two.
pub fn read_options_file(path: &str) -> Result<Vec<OptionsFileEntry>, String> {
    let contents = std::fs::read_to_string(path)
        .map_err(|err| format!("options file {} cannot be opened: {}", path, err))?;
    parse_options_file(&contents)
}

/// The parsed lines of an options file (see [`read_options_file`]).
pub fn parse_options_file(contents: &str) -> Result<Vec<OptionsFileEntry>, String> {
    let mut entries = Vec::new();
    for (index, line) in contents.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let separator = match line.find('=') {
            Some(position) if !line[..position].trim().is_empty() => position,
            _ => {
                return Err(format!(
                    "line {} has no usable 'name = value' pair: '{}'",
                    index + 1,
                    line
                ));
            }
        };
        let name = trim_quotes(line[..separator].trim());
        let value = trim_quotes(line[separator + 1..].trim());
        if name.is_empty() {
            return Err(format!("line {} gives no option name", index + 1));
        }
        if value.is_empty() && option_type_name(&name) != Some("string") {
            // Only a string option may legitimately hold the empty value ("" = no file);
            // everywhere else an empty right-hand side is the mistake HiGHS' reader rejects.
            return Err(format!(
                "line {}: gives no value for option '{}'",
                index + 1,
                name
            ));
        }
        entries.push(OptionsFileEntry {
            line: index + 1,
            name,
            value,
        });
    }
    Ok(entries)
}

fn trim_quotes(text: &str) -> String {
    text.trim_matches(|c: char| c.is_whitespace() || "\"'".contains(c))
        .to_string()
}

/// HiGHS kFull-style text for the given options: per option a description comment, a
/// `[type, advanced, range, default]` comment, then `name = value`. Doubles use Rust's
/// shortest round-trip form so a written file reads back bit-identically.
pub fn report_options(opts: &HighsOptions, non_default_only: Option<bool>) -> String {
    let non_default_only = non_default_only.unwrap_or(false);
    let mut out = String::new();
    for record in OPTION_RECORDS {
        let Some(value) = get_option_field(opts, record.name) else {
            continue;
        };
        let changed = option_default(record.name).is_some_and(|default| value != default);
        if non_default_only && !changed {
            continue;
        }
        let (type_name, range, default) = match &record.kind {
            RecordKind::Bool { default } => ("bool", String::new(), default.to_string()),
            RecordKind::Int {
                lower,
                default,
                upper,
            } => (
                "HighsInt",
                format!("range: [{}, {}], ", lower, upper),
                default.to_string(),
            ),
            RecordKind::Double {
                lower,
                default,
                upper,
            } => (
                "double",
                format!("range: [{}, {}], ", lower, upper),
                default.to_string(),
            ),
            RecordKind::Str { default } => ("string", String::new(), (*default).to_string()),
        };
        out.push_str(&format!("# {}\n", record.description));
        out.push_str(&format!(
            "# [type: {}, advanced: {}, {}default: {}]\n",
            type_name, record.advanced, range, default
        ));
        out.push_str(&format!("{} = {}\n", record.name, value.to_token()));
    }
    out
}

/// Write the options to a file in the format [`report_options`] produces.
pub fn write_options_file(opts: &HighsOptions, path: &str) -> Result<(), String> {
    std::fs::write(path, report_options(opts, None))
        .map_err(|err| format!("options file {}: {}", path, err))
}

/// JSON listing of the registry with current values: one entry per record,
/// `changed` against the record default, `affects_search` for sweep selection.
pub fn option_records_json(
    opts: &HighsOptions,
    non_default_only: Option<bool>,
) -> Result<String, String> {
    let non_default_only = non_default_only.unwrap_or(false);
    let mut entries = Vec::new();
    for record in OPTION_RECORDS {
        let Some(value) = get_option_field(opts, record.name) else {
            continue;
        };
        let default = option_default(record.name);
        let changed = default.as_ref().is_some_and(|default| value != *default);
        if non_default_only && !changed {
            continue;
        }
        let (default_json, lower, upper) = match &record.kind {
            RecordKind::Bool { default } => (serde_json::Value::from(*default), None, None),
            RecordKind::Int {
                lower,
                default,
                upper,
            } => (
                serde_json::Value::from(*default),
                Some(serde_json::Value::from(*lower)),
                Some(serde_json::Value::from(*upper)),
            ),
            RecordKind::Double {
                lower,
                default,
                upper,
            } => (
                json_number(*default),
                Some(json_number(*lower)),
                Some(json_number(*upper)),
            ),
            RecordKind::Str { default } => (serde_json::Value::from(*default), None, None),
        };
        let value_json = match &value {
            OptionValue::Bool(v) => serde_json::Value::from(*v),
            OptionValue::Int(v) => serde_json::Value::from(*v),
            OptionValue::Double(v) => json_number(*v),
            OptionValue::Str(v) => serde_json::Value::from(v.clone()),
        };
        entries.push(JsonRecord {
            name: record.name,
            type_name: value.type_name(),
            description: record.description,
            advanced: record.advanced,
            source: record.source,
            affects: affects_name(record.affects),
            tunable: record.tunable,
            default: default_json,
            lower,
            upper,
            value: value_json,
            changed,
            affects_search: affects_search(record.affects),
        });
    }
    serde_json::to_string_pretty(&entries).map_err(|err| err.to_string())
}

/// The plain (non-JSON) listing: one line per option with its current value, the registry's
/// type/range, and the source, so an agent can read the whole surface at once.
pub fn report_options_table(opts: &HighsOptions, non_default_only: Option<bool>) -> String {
    let non_default_only = non_default_only.unwrap_or(false);
    let mut out = String::new();
    for record in OPTION_RECORDS {
        let Some(value) = get_option_field(opts, record.name) else {
            continue;
        };
        let changed = option_default(record.name).is_some_and(|default| value != default);
        if non_default_only && !changed {
            continue;
        }
        let range = match &record.kind {
            RecordKind::Int { lower, upper, .. } => format!(", range [{}, {}]", lower, upper),
            RecordKind::Double { lower, upper, .. } => {
                format!(", range [{}, {}]", lower, upper)
            }
            _ => String::new(),
        };
        out.push_str(&format!(
            "{} = {} [{}{}; {}{}; {}]\n",
            record.name,
            value.to_token(),
            value.type_name(),
            range,
            affects_name(record.affects),
            if record.affects_search() {
                ", affects search"
            } else {
                ""
            },
            record.source
        ));
    }
    out
}

impl OptionRecord {
    /// Whether this option can move the search.
    pub fn affects_search(&self) -> bool {
        affects_search(self.affects)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::highs_run::Highs;
    use crate::prelude::HighsStatus;

    fn all_field_names() -> Vec<&'static str> {
        let mut names = vec![
            "allow_unbounded_or_infeasible",
            "allowed_cost_scale_factor",
            "allowed_matrix_scale_factor",
            "blend_multi_objectives",
            "centring_ratio_tolerance",
            "cost_scale_factor",
            "dual_feasibility_tolerance",
            "dual_residual_tolerance",
            "dual_simplex_cost_perturbation_multiplier",
            "dual_simplex_pivot_growth_tolerance",
            "dual_steepest_edge_weight_error_tolerance",
            "dual_steepest_edge_weight_log_error_threshold",
            "factor_pivot_threshold",
            "factor_pivot_tolerance",
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
            "h16_fourier_motzkin",
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
            "h16_qp_bounded_start",
            "h16_qp_convexity",
            "h16_ray_simplex",
            "h16_rebuild_reason",
            "h16_repair_lp",
            "h16_row_duals",
            "h16_row_presolve",
            "h16_row_redundant",
            "h16_scaling_report",
            "h16_singleton_loops",
            "h16_sub_solver_flags",
            "h16_subst_int",
            "h16_sym_refine",
            "h16_unresolved_unbounded_or_infeasible",
            "h16_user_scale_offset",
            "h16_user_sol_obj",
            "h16_user_solution_assess",
            "h16_watch_link",
            "h16_zero_cost_singleton",
            "highs_analysis_level",
            "highs_debug_level",
            "hipo_available",
            "hipo_block_size",
            "hipo_ordering",
            "hipo_parallel_type",
            "hipo_system",
            "icrash",
            "icrash_approx_iter",
            "icrash_breakpoints",
            "icrash_dualize",
            "icrash_exact",
            "icrash_iterations",
            "icrash_starting_weight",
            "icrash_strategy",
            "iis_strategy",
            "iis_time_limit",
            "infinite_bound",
            "infinite_cost",
            "ipm_iteration_limit",
            "ipm_optimality_tolerance",
            "ipx_dualize_strategy",
            "kkt_tolerance",
            "large_matrix_value",
            "less_infeasible_DSE_check",
            "less_infeasible_DSE_choose_row",
            "log_dev_level",
            "log_file",
            "log_options",
            "log_to_console",
            "lp_fold",
            "lp_crossover_ipm_crash",
            "lp_crossover_batch_push",
            "lp_dual_devex_start",
            "lp_postsolve_snap_basis",
            "lp_ipm_support_reduction",
            "lp_concurrent",
            "lp_solver_auto",
            "lp_postsolve_fix_tiny",
            "lp_fold_basis",
            "lp_crossover_ipx_refine",
            "lp_ipx_no_drop",
            "lp_ipx_soft_drop",
            "presolve_dep_eq_fixed_budget",
            "lp_dualize_auto",
            "ipm_hipo_fast_blas",
            "ipm_hipo_cache_perm",
            "ipm_hipo_par_ne",
            "ipm_hipo_stall_detect",
            "ipm_hipo_amd_first",
            "ipm_hipo_dense_cols",
            "ipm_hipo_dense_schur",
            "ipm_hipo_dualize",
            "lp_ipm_failure_policy",
            "lp_presolve_requires_basis_postsolve",
            "max_centring_steps",
            "max_dual_simplex_cleanup_level",
            "max_dual_simplex_phase1_cleanup_level",
            "mip_abs_gap",
            "mip_allow_cut_separation_at_nodes",
            "mip_allow_restart",
            "mip_alns",
            "mip_shift_and_propagate",
            "mip_cdcl_dive",
            "mip_root_bound_stall",
            "mip_cdcl_uip",
            "mip_zero_obj_cdcl",
            "mip_decompose",
            "mip_lattice",
            "mip_activation_dff",
            "mip_fj_native",
            "mip_first_inc",
            "mip_first_inc_root",
            "mip_lattice_kernel",
            "mip_racing_threads",
            "mip_conflict_off_mask",
            "mip_presolve_probing_useless",
            "mip_root_lp_concurrent",
            "mip_conflict_resolution",
            "mip_conflict_resolution_tree_only",
        ];
        names.extend_from_slice(&[
            "mip_cut_selection_hybrid",
            "mip_knapsack_cover_sepa",
            "mip_dff_separator",
            "mip_tableau_loosen",
            "mip_tableau_loosen_auto",
            "mip_ti_precedence_sepa",
            "mip_skip_analytic_centre",
            "mip_netdesign_sepa",
            "mip_netdesign_heur",
            "mip_netdesign_cutset",
            "mip_ti_list_sched",
            "mip_cont_prop_improve",
            "mip_pscost_variance_tol",
            "mip_cutoff_bigm",
            "mip_node_bigm",
            "mip_node_bigm_rows",
            "mip_node_bigm_min_ratio",
            "mip_rens_min_fixing_rate",
            "mip_node_bigm_active",
            "mip_debug_solution_file",
            "mip_detect_symmetry",
            "mip_feasibility_tolerance",
            "mip_heuristic_effort",
            "mip_heuristic_run_crossover",
            "mip_heuristic_run_feasibility_jump",
            "mip_heuristic_run_fpr",
            "mip_heuristic_run_local_mip",
            "mip_heuristic_run_oneopt",
            "mip_heuristic_run_rens",
            "mip_heuristic_run_rins",
            "mip_heuristic_run_root_reduced_cost",
            "mip_heuristic_run_shifting",
            "mip_heuristic_run_zi_round",
            "mip_improving_solution_file",
            "mip_improving_solution_report_sparse",
            "mip_improving_solution_save",
            "mip_ipm_solver",
            "mip_lifting_for_probing",
            "mip_guided_dive",
            "mip_local_branching",
            "mip_lp_age_limit",
            "mip_lp_solver",
            "mip_ls_gub_swap",
            "mip_ls_int_lp_repair",
            "mip_max_improving_sols",
            "mip_max_leaves",
            "mip_max_nodes",
            "mip_max_stall_nodes",
            "mip_max_start_nodes",
            "mip_min_cliquetable_entries_for_parallelism",
            "mip_min_logging_interval",
            "mip_pool_age_limit",
            "mip_pool_soft_limit",
            "mip_presolve_eq_subst",
            "mip_presolve_clique_varbound",
            "mip_presolve_xor_gf2",
            "mip_orbitope_root_fix",
            "mip_orbitope_full_fix",
            "mip_proximity",
            "mip_improve_early",
            "mip_primal_stall",
            "mip_improve_early_work",
            "mip_presolve_probing_cap",
            "mip_presolve_domcol_full",
            "mip_redcost_proof_uncapped",
            "mip_root_stall_scip",
            "mip_heur_lp_quota",
            "mip_pscost_minreliable",
            "mip_rel_gap",
            "mip_report_level",
            "mip_root_presolve_only",
            "mip_search_simulate_concurrency",
            "mip_strengthen_cover_guard",
            "mip_symmetry_dynamic",
            "mip_aggregate_implications",
            "mip_cdouble_act_contrib",
            "mip_clique_normalise",
            "mip_coef_change_impl_bounds",
            "mip_dfs_node_limit",
            "mip_dual_fix_probing",
            "mip_dual_side_lost",
            "mip_dual_subst_limits",
            "mip_keep_branch_reason",
            "mip_lurking_bounds_robust",
            "mip_machine_sched_sepa",
            "mip_path_same_bds",
            "mip_ranged_dual_implied",
            "mip_reduced_stack_fast",
            "mip_redundant_vbd",
            "mip_sb_budget_pre_restart",
            "mip_strengthen_small_contrib",
            "mip_cdouble_proofs",
            "mip_fix_last_active",
            "mip_fm_presolve",
            "mip_keep_implications",
            "mip_presolve_clique_table",
            "mip_presolve_row_append",
            "mip_sym_hash64",
            "no_unnecessary_rebuild_refactor",
            "objective_bound",
            "objective_target",
            "optimality_tolerance",
            "output_flag",
            "parallel",
            "pdlp_cupdlpc_restart_method",
            "pdlp_features_off",
            "pdlp_iteration_limit",
            "pdlp_optimality_tolerance",
            "pdlp_restart_strategy",
            "pdlp_ruiz_iterations",
            "pdlp_scaling_mode",
            "pdlp_step_size_strategy",
            "presolve",
            "presolve_fm_obj_reformulation",
            "presolve_mode",
            "presolve_pivot_threshold",
            "presolve_reduction_limit",
            "presolve_remove_slacks",
            "presolve_rule_logging",
            "presolve_rule_off",
            "presolve_rule_test",
            "presolve_substitution_maxfillin",
            "primal_feasibility_tolerance",
            "primal_residual_tolerance",
            "primal_simplex_bound_perturbation_multiplier",
            "qp_allow_hot_start",
            "qp_iteration_limit",
            "qp_nullspace_limit",
            "qp_regularization_value",
            "random_seed",
            "rebuild_refactor_solution_error_tolerance",
            "restart_presolve_reduction_limit",
            "run_centring",
            "run_crossover",
            "simplex_dual_edge_weight_strategy",
            "simplex_dualize_strategy",
            "simplex_iteration_limit",
            "simplex_max_concurrency",
            "simplex_min_concurrency",
            "simplex_permute_strategy",
            "simplex_price_strategy",
            "simplex_primal_edge_weight_strategy",
            "simplex_scale_strategy",
            "simplex_strategy",
            "simplex_unscaled_solution_strategy",
            "simplex_update_limit",
            "small_matrix_value",
            "solve_relaxation",
            "solver",
            "start_crossover_tolerance",
            "threads",
            "time_limit",
            "timeless_log",
            "use_implied_bounds_from_presolve",
            "use_warm_start",
            "user_bound_scale",
            "user_objective_scale",
        ]);
        names.sort_unstable();
        names
    }

    #[test]
    fn opt_registry_field_coverage() {
        // Every option field is either a settable record or an explicitly internal field, so
        // no option can exist that `rusolve options`, the options file and `--opt` do not
        // know about. The binding below is exhaustive: a new field fails to compile until it
        // is registered or declared internal.
        let opts = HighsOptions::default();
        let HighsOptions {
            allow_unbounded_or_infeasible: _,
            allowed_cost_scale_factor: _,
            allowed_matrix_scale_factor: _,
            blend_multi_objectives: _,
            centring_ratio_tolerance: _,
            cost_scale_factor: _,
            dual_feasibility_tolerance: _,
            dual_residual_tolerance: _,
            dual_simplex_cost_perturbation_multiplier: _,
            dual_simplex_pivot_growth_tolerance: _,
            dual_steepest_edge_weight_error_tolerance: _,
            dual_steepest_edge_weight_log_error_threshold: _,
            factor_pivot_threshold: _,
            factor_pivot_tolerance: _,
            h16_ac_centring: _,
            h16_ac_fix: _,
            h16_add_no_modify_on_error: _,
            h16_agg_term_zero: _,
            h16_assess_small: _,
            h16_big_act: _,
            h16_bound_chain: _,
            h16_bound_range_inf: _,
            h16_cdouble: _,
            h16_clique_links: _,
            h16_col_stuffing: _,
            h16_cost_perturb: _,
            h16_crossover_x: _,
            h16_cupdlp_time: _,
            h16_cutpool_age: _,
            h16_dep_eq_time_cap: _,
            h16_dominated_cols: _,
            h16_dual_free_list: _,
            h16_dual_postsolve_rules: _,
            h16_dup_row_undo: _,
            h16_dupcol_warn: _,
            h16_enum_sig64: _,
            h16_eq_row_adds: _,
            h16_feasrelax_sense: _,
            h16_fj: _,
            h16_fourier_motzkin: _,
            h16_frontier_merge: _,
            h16_hash_common: _,
            h16_hfactor_pivot: _,
            h16_hipdlp_max: _,
            h16_hipo_display: _,
            h16_hipo_empty_rows: _,
            h16_impl_eq: _,
            h16_implint_bounds: _,
            h16_int_bound_round: _,
            h16_int_bound_rounding: _,
            h16_int_round_assert_tol: _,
            h16_ipm_gap_offset: _,
            h16_ipx_free_rows: _,
            h16_ipx_scaled: _,
            h16_ipx_stall: _,
            h16_kkt_effective_costs: _,
            h16_max_duals: _,
            h16_mip_start: _,
            h16_mip_workers: _,
            h16_ninfmin: _,
            h16_nla_debug_gate: _,
            h16_obj_cdouble: _,
            h16_parallel_rows: _,
            h16_postsolve_basis_check: _,
            h16_presolve_rules: _,
            h16_qp_bounded_start: _,
            h16_qp_convexity: _,
            h16_ray_simplex: _,
            h16_rebuild_reason: _,
            h16_repair_lp: _,
            h16_row_duals: _,
            h16_row_presolve: _,
            h16_row_redundant: _,
            h16_scaling_report: _,
            h16_singleton_loops: _,
            h16_sub_solver_flags: _,
            h16_subst_int: _,
            h16_sym_refine: _,
            h16_unresolved_unbounded_or_infeasible: _,
            h16_user_scale_offset: _,
            h16_user_sol_obj: _,
            h16_user_solution_assess: _,
            h16_watch_link: _,
            h16_zero_cost_singleton: _,
            highs_analysis_level: _,
            highs_debug_level: _,
            hipo_available: _,
            hipo_block_size: _,
            hipo_ordering: _,
            hipo_parallel_type: _,
            hipo_system: _,
            icrash: _,
            icrash_approx_iter: _,
            icrash_breakpoints: _,
            icrash_dualize: _,
            icrash_exact: _,
            icrash_iterations: _,
            icrash_starting_weight: _,
            icrash_strategy: _,
            iis_strategy: _,
            iis_time_limit: _,
            infinite_bound: _,
            infinite_cost: _,
            ipm_iteration_limit: _,
            ipm_optimality_tolerance: _,
            ipx_dualize_strategy: _,
            kkt_tolerance: _,
            large_matrix_value: _,
            less_infeasible_DSE_check: _,
            less_infeasible_DSE_choose_row: _,
            log_dev_level: _,
            log_file: _,
            log_options: _,
            log_to_console: _,
            lp_fold: _,
            lp_crossover_ipm_crash: _,
            lp_crossover_batch_push: _,
            lp_dual_devex_start: _,
            lp_postsolve_snap_basis: _,
            lp_ipm_support_reduction: _,
            lp_concurrent: _,
            lp_solver_auto: _,
            lp_postsolve_fix_tiny: _,
            lp_fold_basis: _,
            lp_crossover_ipx_refine: _,
            lp_ipx_no_drop: _,
            lp_ipx_soft_drop: _,
            presolve_dep_eq_fixed_budget: _,
            lp_dualize_auto: _,
            ipm_hipo_fast_blas: _,
            ipm_hipo_cache_perm: _,
            ipm_hipo_par_ne: _,
            ipm_hipo_stall_detect: _,
            ipm_hipo_amd_first: _,
            ipm_hipo_dense_cols: _,
            ipm_hipo_dense_schur: _,
            ipm_hipo_dualize: _,
            lp_ipm_failure_policy: _,
            lp_presolve_requires_basis_postsolve: _,
            max_centring_steps: _,
            max_dual_simplex_cleanup_level: _,
            max_dual_simplex_phase1_cleanup_level: _,
            mip_abs_gap: _,
            mip_alns: _,
            mip_presolve_eq_subst: _,
            mip_presolve_clique_varbound: _,
            mip_presolve_xor_gf2: _,
            mip_orbitope_root_fix: _,
            mip_orbitope_full_fix: _,
            mip_shift_and_propagate: _,
            mip_cdcl_dive: _,
            mip_root_bound_stall: _,
            mip_cdcl_uip: _,
            mip_zero_obj_cdcl: _,
            mip_decompose: _,
            mip_lattice: _,
            mip_activation_dff: _,
            mip_fj_native: _,
            mip_first_inc: _,
            mip_first_inc_root: _,
            mip_lattice_kernel: _,
            mip_presolve_probing_cap: _,
            mip_presolve_domcol_full: _,
            mip_redcost_proof_uncapped: _,
            mip_root_stall_scip: _,
            mip_heur_lp_quota: _,
            mip_racing_threads: _,
            mip_conflict_off_mask: _,
            mip_presolve_probing_useless: _,
            mip_root_lp_concurrent: _,
            mip_allow_cut_separation_at_nodes: _,
            mip_allow_restart: _,
            mip_conflict_resolution: _,
            mip_conflict_resolution_tree_only: _,
            mip_cut_selection_hybrid: _,
            mip_knapsack_cover_sepa: _,
            mip_dff_separator: _,
            mip_tableau_loosen: _,
            mip_tableau_loosen_auto: _,
            mip_ti_precedence_sepa: _,
            mip_skip_analytic_centre: _,
            mip_netdesign_sepa: _,
            mip_netdesign_heur: _,
            mip_netdesign_cutset: _,
            mip_ti_list_sched: _,
            mip_cont_prop_improve: _,
            mip_pscost_variance_tol: _,
            mip_cutoff_bigm: _,
            mip_node_bigm: _,
            mip_node_bigm_rows: _,
            mip_node_bigm_min_ratio: _,
            mip_rens_min_fixing_rate: _,
            mip_node_bigm_active: _,
            mip_debug_solution_file: _,
            mip_detect_symmetry: _,
            mip_feasibility_tolerance: _,
            mip_heuristic_effort: _,
            mip_heuristic_run_crossover: _,
            mip_heuristic_run_feasibility_jump: _,
            mip_heuristic_run_fpr: _,
            mip_heuristic_run_local_mip: _,
            mip_heuristic_run_oneopt: _,
            mip_heuristic_run_rens: _,
            mip_heuristic_run_rins: _,
            mip_heuristic_run_root_reduced_cost: _,
            mip_heuristic_run_shifting: _,
            mip_heuristic_run_zi_round: _,
            mip_improving_solution_file: _,
            mip_improving_solution_report_sparse: _,
            mip_improving_solution_save: _,
            mip_ipm_solver: _,
            mip_lifting_for_probing: _,
            mip_guided_dive: _,
            mip_local_branching: _,
            mip_lp_age_limit: _,
            mip_lp_solver: _,
            mip_ls_gub_swap: _,
            mip_ls_int_lp_repair: _,
            mip_max_improving_sols: _,
            mip_max_leaves: _,
            mip_max_nodes: _,
            mip_max_stall_nodes: _,
            mip_max_start_nodes: _,
            mip_min_cliquetable_entries_for_parallelism: _,
            mip_min_logging_interval: _,
            mip_pool_age_limit: _,
            mip_pool_soft_limit: _,
            mip_proximity: _,
            mip_improve_early: _,
            mip_primal_stall: _,
            mip_improve_early_work: _,
            mip_pscost_minreliable: _,
            mip_rel_gap: _,
            mip_report_level: _,
            mip_root_presolve_only: _,
            mip_search_simulate_concurrency: _,
            mip_strengthen_cover_guard: _,
            mip_symmetry_dynamic: _,
            mip_aggregate_implications: _,
            mip_cdouble_act_contrib: _,
            mip_clique_normalise: _,
            mip_coef_change_impl_bounds: _,
            mip_dfs_node_limit: _,
            mip_dual_fix_probing: _,
            mip_dual_side_lost: _,
            mip_dual_subst_limits: _,
            mip_keep_branch_reason: _,
            mip_lurking_bounds_robust: _,
            mip_machine_sched_sepa: _,
            mip_path_same_bds: _,
            mip_ranged_dual_implied: _,
            mip_reduced_stack_fast: _,
            mip_redundant_vbd: _,
            mip_sb_budget_pre_restart: _,
            mip_strengthen_small_contrib: _,
            mip_cdouble_proofs: _,
            mip_fix_last_active: _,
            mip_fm_presolve: _,
            mip_keep_implications: _,
            mip_presolve_clique_table: _,
            mip_presolve_row_append: _,
            mip_sym_hash64: _,
            no_unnecessary_rebuild_refactor: _,
            objective_bound: _,
            objective_target: _,
            optimality_tolerance: _,
            output_flag: _,
            parallel: _,
            pdlp_cupdlpc_restart_method: _,
            pdlp_features_off: _,
            pdlp_iteration_limit: _,
            pdlp_optimality_tolerance: _,
            pdlp_restart_strategy: _,
            pdlp_ruiz_iterations: _,
            pdlp_scaling_mode: _,
            pdlp_step_size_strategy: _,
            presolve: _,
            presolve_fm_obj_reformulation: _,
            presolve_mode: _,
            presolve_pivot_threshold: _,
            presolve_reduction_limit: _,
            presolve_remove_slacks: _,
            presolve_rule_logging: _,
            presolve_rule_off: _,
            presolve_rule_test: _,
            presolve_substitution_maxfillin: _,
            primal_feasibility_tolerance: _,
            primal_residual_tolerance: _,
            primal_simplex_bound_perturbation_multiplier: _,
            qp_allow_hot_start: _,
            qp_iteration_limit: _,
            qp_nullspace_limit: _,
            qp_regularization_value: _,
            random_seed: _,
            rebuild_refactor_solution_error_tolerance: _,
            restart_presolve_reduction_limit: _,
            run_centring: _,
            run_crossover: _,
            simplex_dual_edge_weight_strategy: _,
            simplex_dualize_strategy: _,
            simplex_iteration_limit: _,
            simplex_max_concurrency: _,
            simplex_min_concurrency: _,
            simplex_permute_strategy: _,
            simplex_price_strategy: _,
            simplex_primal_edge_weight_strategy: _,
            simplex_scale_strategy: _,
            simplex_strategy: _,
            simplex_unscaled_solution_strategy: _,
            simplex_update_limit: _,
            small_matrix_value: _,
            solve_relaxation: _,
            solver: _,
            start_crossover_tolerance: _,
            threads: _,
            time_limit: _,
            timeless_log: _,
            use_implied_bounds_from_presolve: _,
            use_warm_start: _,
            user_bound_scale: _,
            user_objective_scale: _,
        } = &opts;
        let mut covered: Vec<&str> = OPTION_RECORDS.iter().map(|record| record.name).collect();
        covered.extend_from_slice(INTERNAL_FIELDS);
        covered.sort_unstable();
        assert_eq!(covered, all_field_names());
        assert_eq!(OPTION_RECORDS.len() + INTERNAL_FIELDS.len(), 321);
    }

    #[test]
    fn opt_registry_defaults_match_options_default() {
        // The registry describes HighsOptions::default(), it never feeds it. A record whose
        // declared default drifts from the struct default would mislead every agent that
        // reads `rusolve options` before deciding what to change.
        let opts = HighsOptions::default();
        for record in OPTION_RECORDS {
            let current = get_option_field(&opts, record.name).unwrap();
            assert_eq!(
                option_default(record.name).unwrap(),
                current,
                "registry default drifts for {}",
                record.name
            );
        }
    }

    /// The typed setter lower-cases string values except for the file-name options (as the
    /// port's 1.15.1 `possibleLowerCaseOptionValue` does), so `icrash_strategy`'s "ICA"
    /// default is stored lower-cased. The round-trip expectations account for that.
    fn stored(record: &OptionRecord, value: &OptionValue) -> OptionValue {
        const CASE_PRESERVED: &[&str] = &[
            "log_file",
            "mip_improving_solution_file",
            "mip_debug_solution_file",
        ];
        match (value, record.name) {
            (OptionValue::Str(text), name) if !CASE_PRESERVED.contains(&name) => {
                OptionValue::Str(text.to_lowercase())
            }
            _ => value.clone(),
        }
    }

    fn typed_get(h: &Highs, record: &OptionRecord) -> OptionValue {
        match &record.kind {
            RecordKind::Bool { .. } => {
                OptionValue::Bool(h.get_option_value_bool(record.name).unwrap())
            }
            RecordKind::Int { .. } => {
                OptionValue::Int(h.get_option_value_int(record.name).unwrap())
            }
            RecordKind::Double { .. } => {
                OptionValue::Double(h.get_option_value_double(record.name).unwrap())
            }
            RecordKind::Str { .. } => {
                OptionValue::Str(h.get_option_value_str(record.name).unwrap())
            }
        }
    }

    #[test]
    fn opt_registry_round_trip_per_record() {
        // Every registered option is reachable by name through the typed setter and getter,
        // and re-setting its own default value is accepted.
        let mut highs = Highs::default();
        for record in OPTION_RECORDS {
            let default = option_default(record.name).unwrap();
            let token = default.to_token();
            highs.set_option_value_str(record.name, &token);
            assert_eq!(
                typed_get(&highs, record),
                stored(record, &default),
                "re-setting the default of {} to {} did not stick",
                record.name,
                token
            );
            assert_eq!(
                option_type_name(record.name).unwrap(),
                default.type_name(),
                "type drift for {}",
                record.name
            );
        }
    }

    #[test]
    fn opt_registry_previously_failing_names_take_effect() {
        // These names used to be rejected by the CLI path (unknown option), so a sweep that
        // set them silently did nothing. The registry now knows them and the typed setter
        // applies them.
        let mut highs = Highs::default();
        for (name, token) in [
            ("mip_heuristic_effort", "0.1"),
            ("mip_detect_symmetry", "false"),
            ("mip_feasibility_tolerance", "0.00000001"),
            ("optimality_tolerance", "0.000000001"),
            ("mip_lp_age_limit", "17"),
            ("mip_lp_solver", "ipm"),
            ("mip_ipm_solver", "ipm"),
            ("mip_allow_restart", "true"),
            ("mip_pool_soft_limit", "3"),
            ("mip_pscost_minreliable", "0"),
        ] {
            highs.set_option_value_str(name, token);
            let kind = record_for(name).unwrap_or_else(|| panic!("{name} is not registered"));
            let got = match &kind.kind {
                RecordKind::Bool { .. } => {
                    OptionValue::Bool(highs.get_option_value_bool(name).unwrap())
                }
                RecordKind::Int { .. } => {
                    OptionValue::Int(highs.get_option_value_int(name).unwrap())
                }
                RecordKind::Double { .. } => {
                    OptionValue::Double(highs.get_option_value_double(name).unwrap())
                }
                RecordKind::Str { .. } => {
                    OptionValue::Str(highs.get_option_value_str(name).unwrap())
                }
            };
            assert_eq!(got.to_token(), token, "{name} did not take the set value");
        }
    }

    #[test]
    fn opt_registry_range_errors() {
        // Declared ranges are enforced by the setter: out-of-range is an error, the boundary
        // values are not.
        let mut highs = Highs::default();
        for record in OPTION_RECORDS {
            let (low, high) = match &record.kind {
                RecordKind::Int { lower, upper, .. } => (*lower as i64, *upper as i64),
                _ => continue,
            };
            if low > i64::from(i32::MIN) {
                assert_eq!(
                    highs.set_option_value_int(record.name, (low - 1) as i32),
                    HighsStatus::kError,
                    "{} accepted a value below its range",
                    record.name
                );
            }
            if high < i64::from(i32::MAX) {
                assert_eq!(
                    highs.set_option_value_int(record.name, (high + 1) as i32),
                    HighsStatus::kError,
                    "{} accepted a value above its range",
                    record.name
                );
            }
            if low > 0 {
                assert_ne!(
                    highs.set_option_value_int(record.name, low as i32),
                    HighsStatus::kError,
                    "{} rejected its own lower bound",
                    record.name
                );
            }
            if high < i64::from(i32::MAX) {
                assert_ne!(
                    highs.set_option_value_int(record.name, high as i32),
                    HighsStatus::kError,
                    "{} rejected its own upper bound",
                    record.name
                );
            }
        }
        for (name, too_big) in [
            ("mip_rel_gap", -1.0),
            ("mip_heuristic_effort", 1.5),
            ("icrash_starting_weight", 1e60),
        ] {
            assert_eq!(
                highs.set_option_value_double(name, too_big),
                HighsStatus::kError,
                "{} accepted an out-of-range double",
                name
            );
        }
    }

    #[test]
    fn opt_registry_unknown_name_suggests_closest() {
        // Requirement A7: a typo names the option the caller almost certainly meant.
        assert_eq!(
            unknown_option_hint("mip_heuristic_effor").as_deref(),
            Some("did you mean mip_heuristic_effort?")
        );
        assert_eq!(
            unknown_option_hint("mip_detect_symetry").as_deref(),
            Some("did you mean mip_detect_symmetry?")
        );
        assert_eq!(
            unknown_option_hint("optimality_tolerence").as_deref(),
            Some("did you mean optimality_tolerance?")
        );
        // An unrelated word gets no bogus suggestion.
        assert_eq!(unknown_option_hint("banana"), None);
        assert!(
            unknown_option_hint("threads").is_none(),
            "a real name is not a typo"
        );
    }

    #[test]
    fn opt_registry_options_file_round_trip() {
        // What the writer emits, the reader reads back: line-for-line, and every value
        // survives the trip.
        let mut highs = Highs::default();
        highs.set_option_value_str("mip_rel_gap", "0.02");
        highs.set_option_value_str("mip_detect_symmetry", "false");
        highs.set_option_value_str("threads", "3");
        highs.set_option_value_str("mip_lp_age_limit", "9");
        highs.set_option_value_str("mip_lp_solver", "ipm");
        let text = report_options(highs.get_options(), None);
        assert!(
            text.contains("# [type: double, advanced: false, range: [0, inf], default: 0.0001]"),
            "{text}"
        );
        let entries = parse_options_file(&text).unwrap();
        assert_eq!(entries.len(), OPTION_RECORDS.len());
        let mut applied = Highs::default();
        for entry in &entries {
            applied.set_option_value_str(&entry.name, &entry.value);
        }
        for record in OPTION_RECORDS {
            assert_eq!(
                typed_get(&applied, record),
                stored(record, &typed_get(&highs, record)),
                "options file round trip lost {}",
                record.name
            );
        }
    }

    #[test]
    fn opt_registry_options_file_bad_lines_report_their_line() {
        // A comment must start with '#' in column 0, exactly as HiGHS' reader requires.
        let err = parse_options_file("mip_rel_gap 0.01\n").unwrap_err();
        assert!(err.starts_with("line 1 "), "{err}");
        assert!(err.contains("no usable"), "{err}");
        let err = parse_options_file("# comment\n # indented comment\n").unwrap_err();
        assert!(err.starts_with("line 2 "), "{err}");
        let err = parse_options_file("\n\nthreads =\n").unwrap_err();
        assert!(
            err.contains("line 3: gives no value for option 'threads'"),
            "{err}"
        );
        let err = parse_options_file(" = 3\n").unwrap_err();
        assert!(err.starts_with("line 1 "), "{err}");
        assert!(err.contains("no usable"), "{err}");
        // Comments and blank lines are skipped, quotes are trimmed.
        let entries = parse_options_file("# header\n\nmip_rel_gap = '0.05'\n").unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].name, "mip_rel_gap");
        assert_eq!(entries[0].value, "0.05");
    }

    #[test]
    fn opt_registry_json_listing_parses() {
        let opts = HighsOptions::default();
        let text = option_records_json(&opts, None).unwrap();
        let parsed: Vec<serde_json::Map<String, serde_json::Value>> =
            serde_json::from_str(&text).unwrap();
        assert_eq!(parsed.len(), OPTION_RECORDS.len());
        for entry in &parsed {
            for key in [
                "name",
                "type_name",
                "description",
                "advanced",
                "source",
                "affects",
                "tunable",
                "default",
                "value",
            ] {
                assert!(
                    entry.contains_key(key),
                    "{key} missing from {}",
                    entry["name"]
                );
            }
            assert!(
                !entry["description"].as_str().unwrap().is_empty(),
                "empty description for {}",
                entry["name"]
            );
            assert!(
                !entry["source"].as_str().unwrap().is_empty(),
                "empty source for {}",
                entry["name"]
            );
        }
        // Nothing is `changed` at the default profile; one option moved shows alone in the
        // non-default listing, with its range and search flag.
        let mut moved = Highs::default();
        moved.set_option_value_str("mip_rel_gap", "0.05");
        let only = option_records_json(moved.get_options(), Some(true)).unwrap();
        let parsed: Vec<serde_json::Map<String, serde_json::Value>> =
            serde_json::from_str(&only).unwrap();
        assert_eq!(parsed.len(), 1, "{only}");
        assert_eq!(parsed[0]["name"], "mip_rel_gap");
        assert_eq!(parsed[0]["changed"], true);
        assert_eq!(parsed[0]["affects"], "tolerance");
        // `false` flags are skipped, so an agent sees only what is worth reading.
        assert!(!parsed[0].contains_key("affects_search"));
        let all: Vec<serde_json::Map<String, serde_json::Value>> =
            serde_json::from_str(&option_records_json(moved.get_options(), None).unwrap()).unwrap();
        let search = all
            .iter()
            .find(|entry| entry["affects"] == "search")
            .expect("no search option in the listing");
        assert_eq!(search["affects_search"], true);
        assert_eq!(search["tunable"], true);
        assert_eq!(parsed[0]["upper"], "inf");
    }

    #[test]
    fn opt_registry_record_metadata_is_sane() {
        let mut seen = std::collections::HashSet::new();
        for record in OPTION_RECORDS {
            assert!(seen.insert(record.name), "duplicate record {}", record.name);
            assert!(
                !record.description.is_empty(),
                "{} has no description",
                record.name
            );
            assert!(!record.source.is_empty(), "{} has no source", record.name);
            assert!(
                record.description.len() < 200,
                "{} description too long",
                record.name
            );
            match &record.kind {
                RecordKind::Int {
                    lower,
                    default,
                    upper,
                } => {
                    assert!(
                        lower <= default && default <= upper,
                        "{} default out of range",
                        record.name
                    );
                }
                RecordKind::Double {
                    lower,
                    default,
                    upper,
                } => {
                    assert!(
                        lower <= default && default <= upper,
                        "{} default out of range",
                        record.name
                    );
                }
                _ => {}
            }
            assert_eq!(
                option_type_name(record.name).unwrap(),
                option_default(record.name).unwrap().type_name(),
                "{} kind disagrees with its default",
                record.name
            );
        }
        // rusolve's own levers are all default-off bools, so the default search is unchanged. The exceptions are
        // `mip_cont_prop_improve`, a numeric lever (default 0.05, SCIP's boundstreps; HiGHS hard-codes 0.3),
        // `lp_fold`, `lp_dualize_auto`, `ipm_hipo_fast_blas`, `lp_postsolve_fix_tiny`, `mip_lattice`, `ipm_hipo_dualize`, `mip_lattice_kernel` and `mip_zero_obj_cdcl` (default-on levers).
        for record in OPTION_RECORDS
            .iter()
            .filter(|record| record.source != "HiGHS 1.16.0")
            .filter(|record| {
                !matches!(
                    record.name,
                    "mip_cont_prop_improve"
                        | "mip_improve_early_work"
                        | "mip_node_bigm_min_ratio"
                        | "mip_rens_min_fixing_rate"
                        | "mip_node_bigm_active"
                        | "mip_pscost_variance_tol"
                        | "mip_presolve_xor_gf2"
                        | "mip_orbitope_root_fix"
                        | "lp_fold"
                        | "lp_concurrent"
                        | "lp_solver_auto"
                        | "lp_crossover_ipx_refine"
                        | "lp_ipm_failure_policy"
                        | "mip_ti_list_sched"
                        | "mip_ti_precedence_sepa"
                        | "mip_netdesign_sepa"
                        | "mip_netdesign_heur"
                        | "mip_netdesign_cutset"
                        | "mip_heur_lp_quota"
                        | "lp_dualize_auto"
                        | "lp_ipx_soft_drop"
                        | "ipm_hipo_fast_blas"
                        | "mip_root_bound_stall"
                        | "mip_lattice"
                        | "mip_lattice_kernel"
                        | "mip_zero_obj_cdcl"
                        | "ipm_hipo_cache_perm"
                        | "ipm_hipo_par_ne"
                        | "lp_postsolve_fix_tiny"
                        | "ipm_hipo_dualize"
                        | "ipm_hipo_amd_first"
                )
            })
        {
            assert!(record.advanced, "{} should be advanced", record.name);
            // the one count-valued lever: 0 = off
            let off = if record.name == "mip_racing_threads" {
                RecordKind::Int {
                    lower: 0,
                    default: 0,
                    upper: 64,
                }
            } else if record.name == "mip_conflict_off_mask" {
                RecordKind::Int {
                    lower: 0,
                    default: 0,
                    upper: 63,
                }
            } else if record.name == "mip_presolve_probing_useless" {
                RecordKind::Int {
                    lower: 0,
                    default: 0,
                    upper: 2147483647,
                }
            } else {
                RecordKind::Bool { default: false }
            };
            assert_eq!(
                record.kind, off,
                "{} is not a default-off lever",
                record.name
            );
            assert!(
                record.description.contains("MANIFEST")
                    || record.description.contains("docs/")
                    || record.description.contains("HiGHS #")
                    || record.description.contains("job ")
                    || record.description.contains("lever"),
                "{} cites nothing: {}",
                record.name,
                record.description
            );
        }
    }
}
