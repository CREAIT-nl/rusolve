//! Typed simplex-trace records (FORMAT.md, `simplex/<lp>.<dual|primal>.jsonl`),
//! mirroring the builders in `tools/highs-oracle/src/simplex.cpp`.
//!
//! Each struct has `from_rec` (decode a fixture record) and `to_obj` (emit
//! the record in the oracle's encoding). A Rust port produces its trace by
//! filling these structs from its own state and calling `to_obj().build()`;
//! [`super::trace::compare_traces`] then diffs the two streams. The
//! round-trip `from_rec(r).to_obj().build() == *r.v` holds for every fixture
//! record (checked by `tests/fixture_selftest.rs`), so the builders encode
//! field names and types exactly as the oracle does.
//!
//! Hash fields (`h_*`, `hash`) hold the u64 FNV-1a value: compute them with
//! [`super::hash`] from the port's arrays.

use serde_json::Value;

use super::emit::Obj;
use super::fixture::Rec;
use super::records::{BasisRec, LpRec, MatrixRec, WorkArraysRec};

/// `HighsSimplexStatus` (`statusObj`).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct SimplexStatus {
    pub initialised_for_new_lp: bool,
    pub is_dualized: bool,
    pub is_permuted: bool,
    pub initialised_for_solve: bool,
    pub has_basis: bool,
    pub has_ar_matrix: bool,
    pub has_nla: bool,
    pub has_dual_steepest_edge_weights: bool,
    pub has_invert: bool,
    pub has_fresh_invert: bool,
    pub has_fresh_rebuild: bool,
    pub has_dual_objective_value: bool,
    pub has_primal_objective_value: bool,
}

impl SimplexStatus {
    pub fn from_rec(r: &Rec) -> Self {
        SimplexStatus {
            initialised_for_new_lp: r.bool("initialised_for_new_lp"),
            is_dualized: r.bool("is_dualized"),
            is_permuted: r.bool("is_permuted"),
            initialised_for_solve: r.bool("initialised_for_solve"),
            has_basis: r.bool("has_basis"),
            has_ar_matrix: r.bool("has_ar_matrix"),
            has_nla: r.bool("has_nla"),
            has_dual_steepest_edge_weights: r.bool("has_dual_steepest_edge_weights"),
            has_invert: r.bool("has_invert"),
            has_fresh_invert: r.bool("has_fresh_invert"),
            has_fresh_rebuild: r.bool("has_fresh_rebuild"),
            has_dual_objective_value: r.bool("has_dual_objective_value"),
            has_primal_objective_value: r.bool("has_primal_objective_value"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::new()
            .b("initialised_for_new_lp", self.initialised_for_new_lp)
            .b("is_dualized", self.is_dualized)
            .b("is_permuted", self.is_permuted)
            .b("initialised_for_solve", self.initialised_for_solve)
            .b("has_basis", self.has_basis)
            .b("has_ar_matrix", self.has_ar_matrix)
            .b("has_nla", self.has_nla)
            .b(
                "has_dual_steepest_edge_weights",
                self.has_dual_steepest_edge_weights,
            )
            .b("has_invert", self.has_invert)
            .b("has_fresh_invert", self.has_fresh_invert)
            .b("has_fresh_rebuild", self.has_fresh_rebuild)
            .b("has_dual_objective_value", self.has_dual_objective_value)
            .b(
                "has_primal_objective_value",
                self.has_primal_objective_value,
            )
    }
}

/// `scalars` (`ekkScalars`): HEkk / HighsSimplexInfo scalars.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct EkkScalars {
    pub iteration_count: i32,
    pub update_count: i32,
    pub updated_dual_objective_value: f64,
    pub updated_primal_objective_value: f64,
    pub dual_objective_value: f64,
    pub primal_objective_value: f64,
    pub num_primal_infeasibilities: i32,
    pub sum_primal_infeasibilities: f64,
    pub num_dual_infeasibilities: i32,
    pub sum_dual_infeasibilities: f64,
    pub costs_perturbed: bool,
    pub costs_shifted: bool,
    pub bounds_perturbed: bool,
    pub bounds_shifted: bool,
    /// `HEkk::cost_perturbation_base_`
    pub cost_perturbation_base: f64,
    /// `HEkk::cost_perturbation_max_abs_cost_`
    pub cost_perturbation_max_abs_cost: f64,
    /// `HEkk::build_synthetic_tick_`
    pub build_synthetic_tick: f64,
    /// `HEkk::total_synthetic_tick_`
    pub total_synthetic_tick: f64,
    pub col_aq_density: f64,
    pub row_ep_density: f64,
    pub row_ap_density: f64,
    pub row_dse_density: f64,
    pub col_bfrt_density: f64,
    pub primal_col_density: f64,
    pub dual_col_density: f64,
    /// `HEkk::dual_simplex_cleanup_level_`
    pub dual_simplex_cleanup_level: i32,
    /// `HEkk::simplex_in_scaled_space_`
    pub simplex_in_scaled_space: bool,
}

impl EkkScalars {
    pub fn from_rec(r: &Rec) -> Self {
        EkkScalars {
            iteration_count: r.int("iteration_count"),
            update_count: r.int("update_count"),
            updated_dual_objective_value: r.f64("updated_dual_objective_value"),
            updated_primal_objective_value: r.f64("updated_primal_objective_value"),
            dual_objective_value: r.f64("dual_objective_value"),
            primal_objective_value: r.f64("primal_objective_value"),
            num_primal_infeasibilities: r.int("num_primal_infeasibilities"),
            sum_primal_infeasibilities: r.f64("sum_primal_infeasibilities"),
            num_dual_infeasibilities: r.int("num_dual_infeasibilities"),
            sum_dual_infeasibilities: r.f64("sum_dual_infeasibilities"),
            costs_perturbed: r.bool("costs_perturbed"),
            costs_shifted: r.bool("costs_shifted"),
            bounds_perturbed: r.bool("bounds_perturbed"),
            bounds_shifted: r.bool("bounds_shifted"),
            cost_perturbation_base: r.f64("cost_perturbation_base"),
            cost_perturbation_max_abs_cost: r.f64("cost_perturbation_max_abs_cost"),
            build_synthetic_tick: r.f64("build_synthetic_tick"),
            total_synthetic_tick: r.f64("total_synthetic_tick"),
            col_aq_density: r.f64("col_aq_density"),
            row_ep_density: r.f64("row_ep_density"),
            row_ap_density: r.f64("row_ap_density"),
            row_dse_density: r.f64("row_DSE_density"),
            col_bfrt_density: r.f64("col_BFRT_density"),
            primal_col_density: r.f64("primal_col_density"),
            dual_col_density: r.f64("dual_col_density"),
            dual_simplex_cleanup_level: r.int("dual_simplex_cleanup_level"),
            simplex_in_scaled_space: r.bool("simplex_in_scaled_space"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::new()
            .i("iteration_count", self.iteration_count)
            .i("update_count", self.update_count)
            .f(
                "updated_dual_objective_value",
                self.updated_dual_objective_value,
            )
            .f(
                "updated_primal_objective_value",
                self.updated_primal_objective_value,
            )
            .f("dual_objective_value", self.dual_objective_value)
            .f("primal_objective_value", self.primal_objective_value)
            .i(
                "num_primal_infeasibilities",
                self.num_primal_infeasibilities,
            )
            .f(
                "sum_primal_infeasibilities",
                self.sum_primal_infeasibilities,
            )
            .i("num_dual_infeasibilities", self.num_dual_infeasibilities)
            .f("sum_dual_infeasibilities", self.sum_dual_infeasibilities)
            .b("costs_perturbed", self.costs_perturbed)
            .b("costs_shifted", self.costs_shifted)
            .b("bounds_perturbed", self.bounds_perturbed)
            .b("bounds_shifted", self.bounds_shifted)
            .f("cost_perturbation_base", self.cost_perturbation_base)
            .f(
                "cost_perturbation_max_abs_cost",
                self.cost_perturbation_max_abs_cost,
            )
            .f("build_synthetic_tick", self.build_synthetic_tick)
            .f("total_synthetic_tick", self.total_synthetic_tick)
            .f("col_aq_density", self.col_aq_density)
            .f("row_ep_density", self.row_ep_density)
            .f("row_ap_density", self.row_ap_density)
            .f("row_DSE_density", self.row_dse_density)
            .f("col_BFRT_density", self.col_bfrt_density)
            .f("primal_col_density", self.primal_col_density)
            .f("dual_col_density", self.dual_col_density)
            .i(
                "dual_simplex_cleanup_level",
                self.dual_simplex_cleanup_level,
            )
            .b("simplex_in_scaled_space", self.simplex_in_scaled_space)
    }
}

/// `hash` object (`hashObj`): FNV-1a of `basis_.basicIndex_` (i32),
/// `basis_.nonbasicMove_` (i8) and the work arrays except `workRange`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct EkkHash {
    pub basic_index: u64,
    pub nonbasic_move: u64,
    pub work_cost: u64,
    pub work_dual: u64,
    pub work_shift: u64,
    pub work_lower: u64,
    pub work_upper: u64,
    pub work_value: u64,
    pub base_value: u64,
    pub base_lower: u64,
    pub base_upper: u64,
    pub dual_edge_weight: u64,
}

impl EkkHash {
    /// Hash the port's arrays exactly as `hashObj` does.
    pub fn of(basic_index: &[i32], nonbasic_move: &[i8], a: &WorkArraysRec) -> Self {
        use super::hash::{fnv_f64, fnv_i8, fnv_i32};
        EkkHash {
            basic_index: fnv_i32(basic_index),
            nonbasic_move: fnv_i8(nonbasic_move),
            work_cost: fnv_f64(&a.work_cost),
            work_dual: fnv_f64(&a.work_dual),
            work_shift: fnv_f64(&a.work_shift),
            work_lower: fnv_f64(&a.work_lower),
            work_upper: fnv_f64(&a.work_upper),
            work_value: fnv_f64(&a.work_value),
            base_value: fnv_f64(&a.base_value),
            base_lower: fnv_f64(&a.base_lower),
            base_upper: fnv_f64(&a.base_upper),
            dual_edge_weight: fnv_f64(&a.dual_edge_weight),
        }
    }

    pub fn from_rec(r: &Rec) -> Self {
        EkkHash {
            basic_index: r.hash("basicIndex"),
            nonbasic_move: r.hash("nonbasicMove"),
            work_cost: r.hash("workCost"),
            work_dual: r.hash("workDual"),
            work_shift: r.hash("workShift"),
            work_lower: r.hash("workLower"),
            work_upper: r.hash("workUpper"),
            work_value: r.hash("workValue"),
            base_value: r.hash("baseValue"),
            base_lower: r.hash("baseLower"),
            base_upper: r.hash("baseUpper"),
            dual_edge_weight: r.hash("dual_edge_weight"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::new()
            .h("basicIndex", self.basic_index)
            .h("nonbasicMove", self.nonbasic_move)
            .h("workCost", self.work_cost)
            .h("workDual", self.work_dual)
            .h("workShift", self.work_shift)
            .h("workLower", self.work_lower)
            .h("workUpper", self.work_upper)
            .h("workValue", self.work_value)
            .h("baseValue", self.base_value)
            .h("baseLower", self.base_lower)
            .h("baseUpper", self.base_upper)
            .h("dual_edge_weight", self.dual_edge_weight)
    }
}

pub fn basis_obj(b: &BasisRec) -> Obj {
    Obj::new()
        .vi("basicIndex", &b.basic_index)
        .vi8("nonbasicFlag", &b.nonbasic_flag)
        .vi8("nonbasicMove", &b.nonbasic_move)
}

pub fn arrays_obj(a: &WorkArraysRec) -> Obj {
    a.named().iter().fold(Obj::new(), |o, (k, v)| o.vf(k, v))
}

pub fn lp_obj(lp: &LpRec) -> Obj {
    let s = &lp.scale;
    Obj::new()
        .i("num_col", lp.num_col)
        .i("num_row", lp.num_row)
        .i("sense", lp.sense)
        .f("offset", lp.offset)
        .vf("col_cost", &lp.col_cost)
        .vf("col_lower", &lp.col_lower)
        .vf("col_upper", &lp.col_upper)
        .vf("row_lower", &lp.row_lower)
        .vf("row_upper", &lp.row_upper)
        .mat("a", &lp.a.view())
        .b("is_scaled", lp.is_scaled)
        .b("is_moved", lp.is_moved)
        .o(
            "scale",
            Obj::new()
                .i("strategy", s.strategy)
                .b("has_scaling", s.has_scaling)
                .i("num_col", s.num_col)
                .i("num_row", s.num_row)
                .f("cost", s.cost)
                .vf("col", &s.col)
                .vf("row", &s.row),
        )
}

/// `primal state` (`HighsOracleAccess::primalState`).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct PrimalState {
    pub solve_phase: i32,
    pub rebuild_reason: i32,
    pub edge_weight_mode: i32,
    pub primal_feasibility_tolerance: f64,
    pub dual_feasibility_tolerance: f64,
    pub num_free_col: i32,
    pub num_flip_since_rebuild: i32,
    /// FNV-1a of `HEkkPrimal::edge_weight_`.
    pub h_edge_weight: u64,
}

impl PrimalState {
    pub fn from_rec(r: &Rec) -> Self {
        PrimalState {
            solve_phase: r.int("solve_phase"),
            rebuild_reason: r.int("rebuild_reason"),
            edge_weight_mode: r.int("edge_weight_mode"),
            primal_feasibility_tolerance: r.f64("primal_feasibility_tolerance"),
            dual_feasibility_tolerance: r.f64("dual_feasibility_tolerance"),
            num_free_col: r.int("num_free_col"),
            num_flip_since_rebuild: r.int("num_flip_since_rebuild"),
            h_edge_weight: r.hash("h_edge_weight"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::new()
            .i("solve_phase", self.solve_phase)
            .i("rebuild_reason", self.rebuild_reason)
            .i("edge_weight_mode", self.edge_weight_mode)
            .f(
                "primal_feasibility_tolerance",
                self.primal_feasibility_tolerance,
            )
            .f(
                "dual_feasibility_tolerance",
                self.dual_feasibility_tolerance,
            )
            .i("num_free_col", self.num_free_col)
            .i("num_flip_since_rebuild", self.num_flip_since_rebuild)
            .h("h_edge_weight", self.h_edge_weight)
    }
}

/// `trace_header`.
#[derive(Clone, Debug, PartialEq)]
pub struct TraceHeader {
    pub name: String,
    pub simplex_strategy: i32,
    pub random_seed: i32,
    pub highs_version: String,
    pub highs_githash: String,
}

impl TraceHeader {
    pub fn from_rec(r: &Rec) -> Self {
        TraceHeader {
            name: r.str("name").into(),
            simplex_strategy: r.int("simplex_strategy"),
            random_seed: r.int("random_seed"),
            highs_version: r.str("highs_version").into(),
            highs_githash: r.str("highs_githash").into(),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::typed("trace_header")
            .s("name", &self.name)
            .i("simplex_strategy", self.simplex_strategy)
            .i("random_seed", self.random_seed)
            .s("highs_version", &self.highs_version)
            .s("highs_githash", &self.highs_githash)
    }
}

/// `ekk_solve_enter`: `HEkk::solve` before `initialiseForSolve`.
#[derive(Clone, Debug, PartialEq)]
pub struct EkkSolveEnter {
    pub solve_call: i32,
    pub iteration_count: i32,
    pub status: SimplexStatus,
    /// Every HighsOptions record (except file paths), kept as raw JSON:
    /// doubles are hex strings, ints/bools numbers, strings strings.
    pub options: Value,
    pub cost_scale: f64,
    pub lp: LpRec,
    /// Only when `status.has_basis`.
    pub basis: Option<BasisRec>,
}

impl EkkSolveEnter {
    pub fn from_rec(r: &Rec) -> Self {
        EkkSolveEnter {
            solve_call: r.int("solve_call"),
            iteration_count: r.int("iteration_count"),
            status: SimplexStatus::from_rec(&r.obj("status")),
            options: r.value("options").clone(),
            cost_scale: r.f64("cost_scale"),
            lp: r.lp("lp"),
            basis: r.has("basis").then(|| r.basis("basis")),
        }
    }

    pub fn to_obj(&self) -> Obj {
        let mut o = Obj::typed("ekk_solve_enter")
            .i("solve_call", self.solve_call)
            .i("iteration_count", self.iteration_count)
            .o("status", self.status.to_obj())
            .raw("options", self.options.clone())
            .f("cost_scale", self.cost_scale)
            .o("lp", lp_obj(&self.lp));
        if let Some(b) = &self.basis {
            o = o.o("basis", basis_obj(b));
        }
        o
    }
}

/// `invert`: end of every `HEkk::computeFactor` that ran INVERT.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct InvertEvent {
    pub solve_call: i32,
    pub iteration_count: i32,
    pub rank_deficiency: i32,
    /// `HFactor::build_synthetic_tick`
    pub factor_build_synthetic_tick: f64,
    /// `HEkk::build_synthetic_tick_`
    pub build_synthetic_tick: f64,
    /// `HEkk::total_synthetic_tick_`
    pub total_synthetic_tick: f64,
    pub basis_matrix_num_el: i32,
    pub invert_num_el: i32,
    pub kernel_dim: i32,
    pub kernel_num_el: i32,
    pub has_invert: bool,
    /// FNV-1a of `basicIndex_` after INVERT's permutation.
    pub h_basic_index: u64,
}

impl InvertEvent {
    pub fn from_rec(r: &Rec) -> Self {
        InvertEvent {
            solve_call: r.int("solve_call"),
            iteration_count: r.int("iteration_count"),
            rank_deficiency: r.int("rank_deficiency"),
            factor_build_synthetic_tick: r.f64("factor_build_synthetic_tick"),
            build_synthetic_tick: r.f64("build_synthetic_tick"),
            total_synthetic_tick: r.f64("total_synthetic_tick"),
            basis_matrix_num_el: r.int("basis_matrix_num_el"),
            invert_num_el: r.int("invert_num_el"),
            kernel_dim: r.int("kernel_dim"),
            kernel_num_el: r.int("kernel_num_el"),
            has_invert: r.bool("has_invert"),
            h_basic_index: r.hash("h_basicIndex"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::typed("invert")
            .i("solve_call", self.solve_call)
            .i("iteration_count", self.iteration_count)
            .i("rank_deficiency", self.rank_deficiency)
            .f(
                "factor_build_synthetic_tick",
                self.factor_build_synthetic_tick,
            )
            .f("build_synthetic_tick", self.build_synthetic_tick)
            .f("total_synthetic_tick", self.total_synthetic_tick)
            .i("basis_matrix_num_el", self.basis_matrix_num_el)
            .i("invert_num_el", self.invert_num_el)
            .i("kernel_dim", self.kernel_dim)
            .i("kernel_num_el", self.kernel_num_el)
            .b("has_invert", self.has_invert)
            .h("h_basicIndex", self.h_basic_index)
    }
}

/// `ekk_solve_initialised`: after `initialiseForSolve`.
#[derive(Clone, Debug, PartialEq)]
pub struct EkkSolveInitialised {
    pub solve_call: i32,
    pub status: SimplexStatus,
    pub basis: BasisRec,
    pub scalars: EkkScalars,
    /// `info_.numTotRandomValue_` (from the SECOND
    /// `initialiseSimplexLpRandomVectors` draw)
    pub num_tot_random_value: Vec<f64>,
    /// `info_.numTotPermutation_`
    pub num_tot_permutation: Vec<i32>,
    /// `info_.numColPermutation_`
    pub num_col_permutation: Vec<i32>,
    pub simplex_strategy: i32,
    pub dual_edge_weight_strategy: i32,
    pub primal_edge_weight_strategy: i32,
    pub price_strategy: i32,
    pub dual_simplex_cost_perturbation_multiplier: f64,
    pub primal_simplex_bound_perturbation_multiplier: f64,
    pub factor_pivot_threshold: f64,
    pub update_limit: i32,
    pub arrays: WorkArraysRec,
    /// `HEkk::ar_matrix_` (row-wise partitioned, format 3)
    pub ar_matrix: MatrixRec,
}

impl EkkSolveInitialised {
    pub fn from_rec(r: &Rec) -> Self {
        EkkSolveInitialised {
            solve_call: r.int("solve_call"),
            status: SimplexStatus::from_rec(&r.obj("status")),
            basis: r.basis("basis"),
            scalars: EkkScalars::from_rec(&r.obj("scalars")),
            num_tot_random_value: r.f64s("numTotRandomValue"),
            num_tot_permutation: r.ints("numTotPermutation"),
            num_col_permutation: r.ints("numColPermutation"),
            simplex_strategy: r.int("simplex_strategy"),
            dual_edge_weight_strategy: r.int("dual_edge_weight_strategy"),
            primal_edge_weight_strategy: r.int("primal_edge_weight_strategy"),
            price_strategy: r.int("price_strategy"),
            dual_simplex_cost_perturbation_multiplier: r
                .f64("dual_simplex_cost_perturbation_multiplier"),
            primal_simplex_bound_perturbation_multiplier: r
                .f64("primal_simplex_bound_perturbation_multiplier"),
            factor_pivot_threshold: r.f64("factor_pivot_threshold"),
            update_limit: r.int("update_limit"),
            arrays: r.arrays("arrays"),
            ar_matrix: r.mat("ar_matrix"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::typed("ekk_solve_initialised")
            .i("solve_call", self.solve_call)
            .o("status", self.status.to_obj())
            .o("basis", basis_obj(&self.basis))
            .o("scalars", self.scalars.to_obj())
            .vf("numTotRandomValue", &self.num_tot_random_value)
            .vi("numTotPermutation", &self.num_tot_permutation)
            .vi("numColPermutation", &self.num_col_permutation)
            .i("simplex_strategy", self.simplex_strategy)
            .i("dual_edge_weight_strategy", self.dual_edge_weight_strategy)
            .i(
                "primal_edge_weight_strategy",
                self.primal_edge_weight_strategy,
            )
            .i("price_strategy", self.price_strategy)
            .f(
                "dual_simplex_cost_perturbation_multiplier",
                self.dual_simplex_cost_perturbation_multiplier,
            )
            .f(
                "primal_simplex_bound_perturbation_multiplier",
                self.primal_simplex_bound_perturbation_multiplier,
            )
            .f("factor_pivot_threshold", self.factor_pivot_threshold)
            .i("update_limit", self.update_limit)
            .o("arrays", arrays_obj(&self.arrays))
            .mat("ar_matrix", &self.ar_matrix.view())
    }
}

/// `dual_start`: `HEkkDual::solve` before the major loop.
#[derive(Clone, Debug, PartialEq)]
pub struct DualStart {
    pub solve_call: i32,
    pub dual_call: i32,
    pub solve_phase: i32,
    /// 0 Dantzig, 1 Devex, 2 steepest edge
    pub edge_weight_mode: i32,
    pub initial_basis_is_logical: bool,
    /// Absent unless `edge_weight_mode` = 2 (uninitialised in HEkkDual).
    pub allow_dual_steepest_edge_to_devex_switch: Option<bool>,
    pub tp: f64,
    pub td: f64,
    pub objective_bound: f64,
    pub force_phase2: bool,
    pub scalars: EkkScalars,
    pub arrays: WorkArraysRec,
}

impl DualStart {
    pub fn from_rec(r: &Rec) -> Self {
        DualStart {
            solve_call: r.int("solve_call"),
            dual_call: r.int("dual_call"),
            solve_phase: r.int("solve_phase"),
            edge_weight_mode: r.int("edge_weight_mode"),
            initial_basis_is_logical: r.bool("initial_basis_is_logical"),
            allow_dual_steepest_edge_to_devex_switch: r
                .has("allow_dual_steepest_edge_to_devex_switch")
                .then(|| r.bool("allow_dual_steepest_edge_to_devex_switch")),
            tp: r.f64("Tp"),
            td: r.f64("Td"),
            objective_bound: r.f64("objective_bound"),
            force_phase2: r.bool("force_phase2"),
            scalars: EkkScalars::from_rec(&r.obj("scalars")),
            arrays: r.arrays("arrays"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::typed("dual_start")
            .i("solve_call", self.solve_call)
            .i("dual_call", self.dual_call)
            .i("solve_phase", self.solve_phase)
            .i("edge_weight_mode", self.edge_weight_mode)
            .b("initial_basis_is_logical", self.initial_basis_is_logical)
            .opt_b(
                "allow_dual_steepest_edge_to_devex_switch",
                self.allow_dual_steepest_edge_to_devex_switch,
            )
            .f("Tp", self.tp)
            .f("Td", self.td)
            .f("objective_bound", self.objective_bound)
            .b("force_phase2", self.force_phase2)
            .o("scalars", self.scalars.to_obj())
            .o("arrays", arrays_obj(&self.arrays))
    }
}

/// `primal_start`: `HEkkPrimal::solve` before the major loop.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimalStart {
    pub solve_call: i32,
    pub primal_call: i32,
    pub state: PrimalState,
    pub scalars: EkkScalars,
    pub arrays: WorkArraysRec,
}

impl PrimalStart {
    pub fn from_rec(r: &Rec) -> Self {
        PrimalStart {
            solve_call: r.int("solve_call"),
            primal_call: r.int("primal_call"),
            state: PrimalState::from_rec(&r.obj("state")),
            scalars: EkkScalars::from_rec(&r.obj("scalars")),
            arrays: r.arrays("arrays"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::typed("primal_start")
            .i("solve_call", self.solve_call)
            .i("primal_call", self.primal_call)
            .o("state", self.state.to_obj())
            .o("scalars", self.scalars.to_obj())
            .o("arrays", arrays_obj(&self.arrays))
    }
}

/// Algorithm-specific part of `rebuild` / `cleanup` / `phase_end`.
#[derive(Clone, Debug, PartialEq)]
pub enum AlgoPart {
    /// dual: the fields named in each record's doc
    Dual,
    /// primal: `state`
    Primal(PrimalState),
}

/// `rebuild`: after every `rebuild()`.
#[derive(Clone, Debug, PartialEq)]
pub struct Rebuild {
    pub solve_call: i32,
    /// `rebuild_reason` on entry (0 first rebuild of a phase, 1 update limit,
    /// 2 synthetic clock, 3 possibly optimal, 4 possibly phase-1 feasible,
    /// 5 possibly primal unbounded, 6 possibly dual unbounded, 7 possibly
    /// singular, 8 primal infeasible in primal, 9 CHUZC fail, 10 force
    /// refactor, 11 excessive primal value, -1 cleanup)
    pub reason: i32,
    /// dual only: `solve_phase`, `rebuild_reason_after`, `dualInfeasCount`,
    /// `edge_weight_mode`
    pub dual: Option<DualRebuildPart>,
    /// primal only: `state`
    pub primal_state: Option<PrimalState>,
    pub scalars: EkkScalars,
    pub hash: EkkHash,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct DualRebuildPart {
    pub solve_phase: i32,
    pub rebuild_reason_after: i32,
    pub dual_infeas_count: i32,
    pub edge_weight_mode: i32,
}

impl Rebuild {
    pub fn from_rec(r: &Rec) -> Self {
        let dual = r.str("algorithm") == "dual";
        Rebuild {
            solve_call: r.int("solve_call"),
            reason: r.int("reason"),
            dual: dual.then(|| DualRebuildPart {
                solve_phase: r.int("solve_phase"),
                rebuild_reason_after: r.int("rebuild_reason_after"),
                dual_infeas_count: r.int("dualInfeasCount"),
                edge_weight_mode: r.int("edge_weight_mode"),
            }),
            primal_state: (!dual).then(|| PrimalState::from_rec(&r.obj("state"))),
            scalars: EkkScalars::from_rec(&r.obj("scalars")),
            hash: EkkHash::from_rec(&r.obj("hash")),
        }
    }

    pub fn to_obj(&self) -> Obj {
        let mut o = Obj::typed("rebuild")
            .s(
                "algorithm",
                if self.dual.is_some() {
                    "dual"
                } else {
                    "primal"
                },
            )
            .i("solve_call", self.solve_call)
            .i("reason", self.reason);
        if let Some(d) = &self.dual {
            o = o
                .i("solve_phase", d.solve_phase)
                .i("rebuild_reason_after", d.rebuild_reason_after)
                .i("dualInfeasCount", d.dual_infeas_count)
                .i("edge_weight_mode", d.edge_weight_mode);
        }
        if let Some(s) = &self.primal_state {
            o = o.o("state", s.to_obj());
        }
        o.o("scalars", self.scalars.to_obj())
            .o("hash", self.hash.to_obj())
    }
}

/// Dual `iter` (an iteration happened) or `iter_none` (it did not).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct DualIter {
    /// `iteration_count_` after the iteration
    pub iteration: i32,
    pub solve_phase: i32,
    pub rebuild_reason: i32,
    pub row_out: i32,
    /// `None` for `iter_none`.
    pub done: Option<DualIterDone>,
}

/// The fields of a dual `iter` beyond those of `iter_none`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct DualIterDone {
    pub variable_out: i32,
    pub variable_in: i32,
    pub move_out: i32,
    pub delta_primal: f64,
    pub theta_dual: f64,
    pub theta_primal: f64,
    pub alpha_col: f64,
    pub alpha_row: f64,
    pub numerical_trouble: f64,
    /// DSE weight of `row_out` computed in CHUZR (the weight actually used)
    pub computed_edge_weight: f64,
    pub edge_weight_mode: i32,
    /// `dualRow.workData[0..workCount).first`: BFRT-flipped variables in
    /// order; `num_flips` is its length.
    pub flips: Vec<i32>,
    pub update_count: i32,
    pub updated_dual_objective_value: f64,
    pub total_synthetic_tick: f64,
    pub row_ep_density: f64,
    pub row_ap_density: f64,
    pub col_aq_density: f64,
    pub row_dse_density: f64,
    /// `dual_edge_weight_[row_out]` after the update; the oracle emits it only
    /// when `0 <= row_out < dual_edge_weight_.len()`.
    pub edge_weight_row_out_after: Option<f64>,
    pub h_basic_index: u64,
    pub h_work_dual: u64,
    pub h_base_value: u64,
    pub h_dual_edge_weight: u64,
}

impl DualIter {
    pub fn from_rec(r: &Rec) -> Self {
        let done = r.t() == "iter";
        DualIter {
            iteration: r.int("iteration"),
            solve_phase: r.int("solve_phase"),
            rebuild_reason: r.int("rebuild_reason"),
            row_out: r.int("row_out"),
            done: done.then(|| {
                let flips = r.ints("flips");
                assert_eq!(
                    flips.len() as i32,
                    r.int("num_flips"),
                    "{}: num_flips != flips.len()",
                    r.ctx
                );
                DualIterDone {
                    variable_out: r.int("variable_out"),
                    variable_in: r.int("variable_in"),
                    move_out: r.int("move_out"),
                    delta_primal: r.f64("delta_primal"),
                    theta_dual: r.f64("theta_dual"),
                    theta_primal: r.f64("theta_primal"),
                    alpha_col: r.f64("alpha_col"),
                    alpha_row: r.f64("alpha_row"),
                    numerical_trouble: r.f64("numericalTrouble"),
                    computed_edge_weight: r.f64("computed_edge_weight"),
                    edge_weight_mode: r.int("edge_weight_mode"),
                    flips,
                    update_count: r.int("update_count"),
                    updated_dual_objective_value: r.f64("updated_dual_objective_value"),
                    total_synthetic_tick: r.f64("total_synthetic_tick"),
                    row_ep_density: r.f64("row_ep_density"),
                    row_ap_density: r.f64("row_ap_density"),
                    col_aq_density: r.f64("col_aq_density"),
                    row_dse_density: r.f64("row_DSE_density"),
                    edge_weight_row_out_after: r.opt_f64("edge_weight_row_out_after"),
                    h_basic_index: r.hash("h_basicIndex"),
                    h_work_dual: r.hash("h_workDual"),
                    h_base_value: r.hash("h_baseValue"),
                    h_dual_edge_weight: r.hash("h_dual_edge_weight"),
                }
            }),
        }
    }

    pub fn to_obj(&self) -> Obj {
        let t = if self.done.is_some() {
            "iter"
        } else {
            "iter_none"
        };
        let o = Obj::typed(t)
            .s("algorithm", "dual")
            .i("iteration", self.iteration)
            .i("solve_phase", self.solve_phase)
            .i("rebuild_reason", self.rebuild_reason)
            .i("row_out", self.row_out);
        let Some(d) = &self.done else { return o };
        o.i("variable_out", d.variable_out)
            .i("variable_in", d.variable_in)
            .i("move_out", d.move_out)
            .f("delta_primal", d.delta_primal)
            .f("theta_dual", d.theta_dual)
            .f("theta_primal", d.theta_primal)
            .f("alpha_col", d.alpha_col)
            .f("alpha_row", d.alpha_row)
            .f("numericalTrouble", d.numerical_trouble)
            .f("computed_edge_weight", d.computed_edge_weight)
            .i("edge_weight_mode", d.edge_weight_mode)
            .i("num_flips", d.flips.len() as i64)
            .vi("flips", &d.flips)
            .i("update_count", d.update_count)
            .f(
                "updated_dual_objective_value",
                d.updated_dual_objective_value,
            )
            .f("total_synthetic_tick", d.total_synthetic_tick)
            .f("row_ep_density", d.row_ep_density)
            .f("row_ap_density", d.row_ap_density)
            .f("col_aq_density", d.col_aq_density)
            .f("row_DSE_density", d.row_dse_density)
            .opt_f("edge_weight_row_out_after", d.edge_weight_row_out_after)
            .h("h_basicIndex", d.h_basic_index)
            .h("h_workDual", d.h_work_dual)
            .h("h_baseValue", d.h_base_value)
            .h("h_dual_edge_weight", d.h_dual_edge_weight)
    }
}

/// The `p` object of a primal `iter` / `iter_none` (`primalIter`).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct PrimalIterP {
    pub solve_phase: i32,
    pub rebuild_reason: i32,
    pub variable_in: i32,
    /// -1 = bound flip
    pub row_out: i32,
    /// The pivot fields, `move_in` .. `numericalTrouble`: `None` in an
    /// `iter_none` before the first iteration of the `HEkkPrimal::solve`
    /// call (uninitialised members, not emitted).
    pub move_in: Option<i32>,
    pub variable_out: Option<i32>,
    pub move_out: Option<i32>,
    pub theta_dual: Option<f64>,
    pub theta_primal: Option<f64>,
    pub value_in: Option<f64>,
    pub alpha_col: Option<f64>,
    pub alpha_row: Option<f64>,
    pub numerical_trouble: Option<f64>,
    pub num_flip_since_rebuild: i32,
    pub edge_weight_mode: i32,
    /// FNV-1a of `HEkkPrimal::edge_weight_`
    pub h_edge_weight: u64,
    /// `edge_weight_[variable_in]`; the oracle emits it only when
    /// `0 <= variable_in < edge_weight_.len()`.
    pub edge_weight_in: Option<f64>,
}

impl PrimalIterP {
    pub fn from_rec(r: &Rec) -> Self {
        PrimalIterP {
            solve_phase: r.int("solve_phase"),
            rebuild_reason: r.int("rebuild_reason"),
            variable_in: r.int("variable_in"),
            row_out: r.int("row_out"),
            move_in: r.opt_int("move_in"),
            variable_out: r.opt_int("variable_out"),
            move_out: r.opt_int("move_out"),
            theta_dual: r.opt_f64("theta_dual"),
            theta_primal: r.opt_f64("theta_primal"),
            value_in: r.opt_f64("value_in"),
            alpha_col: r.opt_f64("alpha_col"),
            alpha_row: r.opt_f64("alpha_row"),
            numerical_trouble: r.opt_f64("numericalTrouble"),
            num_flip_since_rebuild: r.int("num_flip_since_rebuild"),
            edge_weight_mode: r.int("edge_weight_mode"),
            h_edge_weight: r.hash("h_edge_weight"),
            edge_weight_in: r.opt_f64("edge_weight_in"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::new()
            .i("solve_phase", self.solve_phase)
            .i("rebuild_reason", self.rebuild_reason)
            .i("variable_in", self.variable_in)
            .i("row_out", self.row_out)
            .opt_i("move_in", self.move_in)
            .opt_i("variable_out", self.variable_out)
            .opt_i("move_out", self.move_out)
            .opt_f("theta_dual", self.theta_dual)
            .opt_f("theta_primal", self.theta_primal)
            .opt_f("value_in", self.value_in)
            .opt_f("alpha_col", self.alpha_col)
            .opt_f("alpha_row", self.alpha_row)
            .opt_f("numericalTrouble", self.numerical_trouble)
            .i("num_flip_since_rebuild", self.num_flip_since_rebuild)
            .i("edge_weight_mode", self.edge_weight_mode)
            .h("h_edge_weight", self.h_edge_weight)
            .opt_f("edge_weight_in", self.edge_weight_in)
    }
}

/// Primal `iter` / `iter_none`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct PrimalIter {
    pub iteration: i32,
    pub p: PrimalIterP,
    /// `None` for `iter_none`.
    pub done: Option<PrimalIterDone>,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct PrimalIterDone {
    pub update_count: i32,
    pub updated_primal_objective_value: f64,
    pub total_synthetic_tick: f64,
    pub h_basic_index: u64,
    pub h_work_dual: u64,
    pub h_base_value: u64,
    pub h_work_value: u64,
}

impl PrimalIter {
    pub fn from_rec(r: &Rec) -> Self {
        let done = r.t() == "iter";
        PrimalIter {
            iteration: r.int("iteration"),
            p: PrimalIterP::from_rec(&r.obj("p")),
            done: done.then(|| PrimalIterDone {
                update_count: r.int("update_count"),
                updated_primal_objective_value: r.f64("updated_primal_objective_value"),
                total_synthetic_tick: r.f64("total_synthetic_tick"),
                h_basic_index: r.hash("h_basicIndex"),
                h_work_dual: r.hash("h_workDual"),
                h_base_value: r.hash("h_baseValue"),
                h_work_value: r.hash("h_workValue"),
            }),
        }
    }

    pub fn to_obj(&self) -> Obj {
        let t = if self.done.is_some() {
            "iter"
        } else {
            "iter_none"
        };
        let o = Obj::typed(t)
            .s("algorithm", "primal")
            .i("iteration", self.iteration)
            .o("p", self.p.to_obj());
        let Some(d) = &self.done else { return o };
        o.i("update_count", d.update_count)
            .f(
                "updated_primal_objective_value",
                d.updated_primal_objective_value,
            )
            .f("total_synthetic_tick", d.total_synthetic_tick)
            .h("h_basicIndex", d.h_basic_index)
            .h("h_workDual", d.h_work_dual)
            .h("h_baseValue", d.h_base_value)
            .h("h_workValue", d.h_work_value)
    }
}

/// `cleanup`: end of `cleanup()`.
#[derive(Clone, Debug, PartialEq)]
pub struct Cleanup {
    pub solve_call: i32,
    /// dual only
    pub solve_phase: Option<i32>,
    /// primal only
    pub primal_state: Option<PrimalState>,
    pub scalars: EkkScalars,
    pub hash: EkkHash,
}

impl Cleanup {
    pub fn from_rec(r: &Rec) -> Self {
        let dual = r.str("algorithm") == "dual";
        Cleanup {
            solve_call: r.int("solve_call"),
            solve_phase: dual.then(|| r.int("solve_phase")),
            primal_state: (!dual).then(|| PrimalState::from_rec(&r.obj("state"))),
            scalars: EkkScalars::from_rec(&r.obj("scalars")),
            hash: EkkHash::from_rec(&r.obj("hash")),
        }
    }

    pub fn to_obj(&self) -> Obj {
        let mut o = Obj::typed("cleanup")
            .s(
                "algorithm",
                if self.primal_state.is_none() {
                    "dual"
                } else {
                    "primal"
                },
            )
            .i("solve_call", self.solve_call);
        if let Some(p) = self.solve_phase {
            o = o.i("solve_phase", p);
        }
        if let Some(s) = &self.primal_state {
            o = o.o("state", s.to_obj());
        }
        o.o("scalars", self.scalars.to_obj())
            .o("hash", self.hash.to_obj())
    }
}

/// `phase_end`: after `solvePhase1()` / `solvePhase2()` returns.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseEnd {
    pub solve_call: i32,
    /// dual only: (`solve_phase` decided next, 0 = exit; `rebuild_reason`)
    pub dual: Option<(i32, i32)>,
    /// primal only
    pub primal_state: Option<PrimalState>,
    pub model_status: i32,
    pub iteration: i32,
}

impl PhaseEnd {
    pub fn from_rec(r: &Rec) -> Self {
        let dual = r.str("algorithm") == "dual";
        PhaseEnd {
            solve_call: r.int("solve_call"),
            dual: dual.then(|| (r.int("solve_phase"), r.int("rebuild_reason"))),
            primal_state: (!dual).then(|| PrimalState::from_rec(&r.obj("state"))),
            model_status: r.int("model_status"),
            iteration: r.int("iteration"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        let mut o = Obj::typed("phase_end")
            .s(
                "algorithm",
                if self.dual.is_some() {
                    "dual"
                } else {
                    "primal"
                },
            )
            .i("solve_call", self.solve_call);
        if let Some((phase, reason)) = self.dual {
            o = o.i("solve_phase", phase).i("rebuild_reason", reason);
        }
        if let Some(s) = &self.primal_state {
            o = o.o("state", s.to_obj());
        }
        o.i("model_status", self.model_status)
            .i("iteration", self.iteration)
    }
}

/// `ekk_solve_exit`: start of `HEkk::returnFromEkkSolve` (after
/// `HEkk::returnFromSolve`, which calls `zeroBasicDuals()`, so `workDual` of
/// basic variables is exactly 0 here but not in the last `cleanup` hash).
#[derive(Clone, Debug, PartialEq)]
pub struct EkkSolveExit {
    pub solve_call: i32,
    /// HighsModelStatus as int (7 = optimal)
    pub model_status: i32,
    pub exit_algorithm: i32,
    pub solve_bailout: bool,
    pub scalars: EkkScalars,
    pub status: SimplexStatus,
    pub basis: BasisRec,
    pub arrays: WorkArraysRec,
}

impl EkkSolveExit {
    pub fn from_rec(r: &Rec) -> Self {
        EkkSolveExit {
            solve_call: r.int("solve_call"),
            model_status: r.int("model_status"),
            exit_algorithm: r.int("exit_algorithm"),
            solve_bailout: r.bool("solve_bailout"),
            scalars: EkkScalars::from_rec(&r.obj("scalars")),
            status: SimplexStatus::from_rec(&r.obj("status")),
            basis: r.basis("basis"),
            arrays: r.arrays("arrays"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::typed("ekk_solve_exit")
            .i("solve_call", self.solve_call)
            .i("model_status", self.model_status)
            .i("exit_algorithm", self.exit_algorithm)
            .b("solve_bailout", self.solve_bailout)
            .o("scalars", self.scalars.to_obj())
            .o("status", self.status.to_obj())
            .o("basis", basis_obj(&self.basis))
            .o("arrays", arrays_obj(&self.arrays))
    }
}

/// `result`: the `Highs::run` outcome (unscaled solution).
#[derive(Clone, Debug, PartialEq)]
pub struct RunResult {
    pub run_status: i32,
    pub model_status: i32,
    pub model_status_string: String,
    pub objective_function_value: f64,
    pub simplex_iteration_count: i32,
    pub primal_solution_status: i32,
    pub dual_solution_status: i32,
    pub num_primal_infeasibilities: i32,
    pub num_dual_infeasibilities: i32,
    pub basis_valid: bool,
    /// HighsBasisStatus as int
    pub col_status: Vec<i32>,
    pub row_status: Vec<i32>,
    pub col_value: Vec<f64>,
    pub col_dual: Vec<f64>,
    pub row_value: Vec<f64>,
    pub row_dual: Vec<f64>,
}

impl RunResult {
    pub fn from_rec(r: &Rec) -> Self {
        RunResult {
            run_status: r.int("run_status"),
            model_status: r.int("model_status"),
            model_status_string: r.str("model_status_string").into(),
            objective_function_value: r.f64("objective_function_value"),
            simplex_iteration_count: r.int("simplex_iteration_count"),
            primal_solution_status: r.int("primal_solution_status"),
            dual_solution_status: r.int("dual_solution_status"),
            num_primal_infeasibilities: r.int("num_primal_infeasibilities"),
            num_dual_infeasibilities: r.int("num_dual_infeasibilities"),
            basis_valid: r.bool("basis_valid"),
            col_status: r.ints("col_status"),
            row_status: r.ints("row_status"),
            col_value: r.f64s("col_value"),
            col_dual: r.f64s("col_dual"),
            row_value: r.f64s("row_value"),
            row_dual: r.f64s("row_dual"),
        }
    }

    pub fn to_obj(&self) -> Obj {
        Obj::typed("result")
            .i("run_status", self.run_status)
            .i("model_status", self.model_status)
            .s("model_status_string", &self.model_status_string)
            .f("objective_function_value", self.objective_function_value)
            .i("simplex_iteration_count", self.simplex_iteration_count)
            .i("primal_solution_status", self.primal_solution_status)
            .i("dual_solution_status", self.dual_solution_status)
            .i(
                "num_primal_infeasibilities",
                self.num_primal_infeasibilities,
            )
            .i("num_dual_infeasibilities", self.num_dual_infeasibilities)
            .b("basis_valid", self.basis_valid)
            .vi("col_status", &self.col_status)
            .vi("row_status", &self.row_status)
            .vf("col_value", &self.col_value)
            .vf("col_dual", &self.col_dual)
            .vf("row_value", &self.row_value)
            .vf("row_dual", &self.row_dual)
    }
}

/// Any simplex-trace record.
#[derive(Clone, Debug, PartialEq)]
pub enum Event {
    Header(TraceHeader),
    SolveEnter(Box<EkkSolveEnter>),
    Invert(InvertEvent),
    Initialised(Box<EkkSolveInitialised>),
    DualStart(Box<DualStart>),
    PrimalStart(Box<PrimalStart>),
    Rebuild(Box<Rebuild>),
    DualIter(Box<DualIter>),
    PrimalIter(Box<PrimalIter>),
    Cleanup(Box<Cleanup>),
    PhaseEnd(PhaseEnd),
    SolveExit(Box<EkkSolveExit>),
    Result(Box<RunResult>),
}

impl Event {
    pub fn from_rec(r: &Rec) -> Self {
        match r.t() {
            "trace_header" => Event::Header(TraceHeader::from_rec(r)),
            "ekk_solve_enter" => Event::SolveEnter(Box::new(EkkSolveEnter::from_rec(r))),
            "invert" => Event::Invert(InvertEvent::from_rec(r)),
            "ekk_solve_initialised" => {
                Event::Initialised(Box::new(EkkSolveInitialised::from_rec(r)))
            }
            "dual_start" => Event::DualStart(Box::new(DualStart::from_rec(r))),
            "primal_start" => Event::PrimalStart(Box::new(PrimalStart::from_rec(r))),
            "rebuild" => Event::Rebuild(Box::new(Rebuild::from_rec(r))),
            "iter" | "iter_none" => match r.str("algorithm") {
                "dual" => Event::DualIter(Box::new(DualIter::from_rec(r))),
                "primal" => Event::PrimalIter(Box::new(PrimalIter::from_rec(r))),
                a => panic!("{}: unknown algorithm {a:?}", r.ctx),
            },
            "cleanup" => Event::Cleanup(Box::new(Cleanup::from_rec(r))),
            "phase_end" => Event::PhaseEnd(PhaseEnd::from_rec(r)),
            "ekk_solve_exit" => Event::SolveExit(Box::new(EkkSolveExit::from_rec(r))),
            "result" => Event::Result(Box::new(RunResult::from_rec(r))),
            t => panic!("{}: unknown simplex trace record type {t:?}", r.ctx),
        }
    }

    pub fn to_obj(&self) -> Obj {
        match self {
            Event::Header(e) => e.to_obj(),
            Event::SolveEnter(e) => e.to_obj(),
            Event::Invert(e) => e.to_obj(),
            Event::Initialised(e) => e.to_obj(),
            Event::DualStart(e) => e.to_obj(),
            Event::PrimalStart(e) => e.to_obj(),
            Event::Rebuild(e) => e.to_obj(),
            Event::DualIter(e) => e.to_obj(),
            Event::PrimalIter(e) => e.to_obj(),
            Event::Cleanup(e) => e.to_obj(),
            Event::PhaseEnd(e) => e.to_obj(),
            Event::SolveExit(e) => e.to_obj(),
            Event::Result(e) => e.to_obj(),
        }
    }

    pub fn to_value(&self) -> Value {
        self.to_obj().build()
    }
}
