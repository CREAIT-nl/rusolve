//! HiGHS source: simplex/HEkkDual.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The dual simplex solver struct. Bodies: `solve.rs`, `phases.rs`, `iterate.rs`,
//! `update.rs`, `assess.rs` (HEkkDual.cpp by line range) and `multi/` (HEkkDualMulti.cpp by
//! line range).

use crate::prelude::*;

mod assess;
mod iterate;
mod multi;
mod phases;
mod solve;
mod update;

/// HiGHS: HEkkDual.h::HEkkDual
/// # Removed C++ aliases (pattern (a))
///
/// | C++ member          | Rust expression                         |
/// |---------------------|-----------------------------------------|
/// | `ekk_instance_`     | `ekk` (method parameter)                |
/// | `a_matrix`          | `ekk.lp_.a_matrix_`                     |
/// | `simplex_nla`       | `ekk.simplex_nla_`                      |
/// | `analysis`          | `ekk.analysis_`                         |
/// | `jMove`             | `ekk.basis_.nonbasic_move_`             |
/// | `workDual`          | `ekk.info_.work_dual_`                  |
/// | `workValue`         | `ekk.info_.work_value_`                 |
/// | `workRange`         | `ekk.info_.work_range_`                 |
/// | `baseLower`         | `ekk.info_.base_lower_`                 |
/// | `baseUpper`         | `ekk.info_.base_upper_`                 |
/// | `baseValue`         | `ekk.info_.base_value_`                 |
/// | `colLower`, `colUpper`, `rowLower`, `rowUpper`, `nonbasicFlag` | never assigned in HiGHS; dropped |
///
/// `dualRow.f(...)` -> `self.dual_row.f(ekk, ...)`, `dualRHS.f(...)` -> `self.dual_rhs.f(ekk, ...)`.
/// `simplex_nla->ftran(col_aq, d, clock)` -> build `NlaRefs::new(&ekk.lp_, &ekk.basis_.basic_index_, &ekk.options_)`
/// and call `ekk.simplex_nla_.ftran(&refs, &mut self.col_aq, d, None)` (borrow `ekk` fields disjointly).
/// PAMI members `multi_choice` / `multi_finish` are the fixed arrays
/// `[MChoice; kSimplexConcurrencyLimit]` / `[MFinish; kSimplexConcurrencyLimit]`.
#[derive(Clone, Debug)]
pub struct HEkkDual {
    pub algorithm: SimplexAlgorithm,
    pub num_devex_iterations: HighsInt,
    pub new_devex_framework: bool,
    pub minor_new_devex_framework: bool,
    pub solver_num_row: HighsInt,
    pub solver_num_col: HighsInt,
    pub solver_num_tot: HighsInt,
    pub inv_solver_num_row: f64,
    pub initial_basis_is_logical_: bool,
    pub edge_weight_mode: EdgeWeightMode,
    pub allow_dual_steepest_edge_to_devex_switch: bool,
    pub tp: f64,
    pub primal_feasibility_tolerance: f64,
    pub td: f64,
    pub dual_feasibility_tolerance: f64,
    pub objective_bound: f64,
    pub force_phase2: bool,
    pub solve_phase: HighsInt,
    pub rebuild_reason: HighsInt,
    pub row_ep: HVector,
    pub row_ap: HVector,
    pub col_aq: HVector,
    pub col_bfrt: HVector,
    pub col_dse: HVector,
    pub dev_row_ep: HVector,
    pub dev_col_dse: HVector,
    pub dual_row: HEkkDualRow,
    pub dual_infeas_count: HighsInt,
    pub dual_rhs: HEkkDualRHS,
    pub row_out: HighsInt,
    pub variable_out: HighsInt,
    pub move_out: HighsInt,
    pub variable_in: HighsInt,
    pub delta_primal: f64,
    pub theta_dual: f64,
    pub theta_primal: f64,
    pub alpha_col: f64,
    pub alpha_row: f64,
    pub numerical_trouble: f64,
    pub computed_edge_weight: f64,
    pub check_invert_condition: bool,
    pub slice_num: HighsInt,
    pub slice_price: HighsInt,
    pub slice_start: Vec<HighsInt>,
    pub slice_a_matrix: Vec<HighsSparseMatrix>,
    pub slice_ar_matrix: Vec<HighsSparseMatrix>,
    pub slice_row_ap: Vec<HVector>,
    pub slice_dual_row: Vec<HEkkDualRow>,
    pub multi_num: HighsInt,
    pub multi_chosen: HighsInt,
    pub multi_i_choice: HighsInt,
    pub multi_n_finish: HighsInt,
    pub multi_iteration: HighsInt,
    pub multi_choose_again: HighsInt,
    /// HiGHS: HEkkDual.h `MChoice multi_choice[kSimplexConcurrencyLimit]` (line 509)
    pub multi_choice: [MChoice; kSimplexConcurrencyLimit as usize],
    /// HiGHS: HEkkDual.h `MFinish multi_finish[kSimplexConcurrencyLimit]` (line 510)
    pub multi_finish: [MFinish; kSimplexConcurrencyLimit as usize],
}

/// HiGHS: HEkkDual.h::HEkkDual::MChoice (lines 469-480): multiple CHUZR data.
#[derive(Clone, Debug, Default)]
pub struct MChoice {
    pub row_out: HighsInt,
    pub base_value: f64,
    pub base_lower: f64,
    pub base_upper: f64,
    pub infeas_value: f64,
    pub infeas_ed_wt: f64,
    pub infeas_limit: f64,
    pub row_ep: HVector,
    pub col_aq: HVector,
    pub col_bfrt: HVector,
}

/// HiGHS: HEkkDual.h::HEkkDual::MFinish (lines 485-501): multiple minor iteration data.
///
/// C++ `HVector_ptr row_ep, col_aq, col_BFRT` are only ever assigned together, to
/// `&multi_choice[iChoice].{row_ep, col_aq, col_BFRT}` (HEkkDualMulti.cpp 288-290), so the
/// three pointers are one index: `choice = Some(iChoice)`; `None` is the C++ indeterminate
/// value before the first minor iteration (never read).
#[derive(Clone, Debug, Default)]
pub struct MFinish {
    pub move_in: i8,
    pub shift_out: f64,
    pub flip_list: Vec<HighsInt>,
    pub row_out: HighsInt,
    pub variable_out: HighsInt,
    pub variable_in: HighsInt,
    pub alpha_row: f64,
    pub theta_primal: f64,
    pub basic_bound: f64,
    pub basic_value: f64,
    pub ed_wt: f64,
    /// C++ `row_ep` / `col_aq` / `col_BFRT` = `&multi_choice[choice].{..}`
    pub choice: Option<usize>,
}

impl Default for HEkkDual {
    fn default() -> Self {
        Self {
            algorithm: SimplexAlgorithm::kDual,
            num_devex_iterations: 0,
            new_devex_framework: false,
            minor_new_devex_framework: false,
            solver_num_row: Default::default(),
            solver_num_col: Default::default(),
            solver_num_tot: Default::default(),
            inv_solver_num_row: Default::default(),
            initial_basis_is_logical_: Default::default(),
            edge_weight_mode: Default::default(),
            allow_dual_steepest_edge_to_devex_switch: Default::default(),
            tp: Default::default(),
            primal_feasibility_tolerance: Default::default(),
            td: Default::default(),
            dual_feasibility_tolerance: Default::default(),
            objective_bound: Default::default(),
            force_phase2: Default::default(),
            solve_phase: Default::default(),
            rebuild_reason: Default::default(),
            row_ep: Default::default(),
            row_ap: Default::default(),
            col_aq: Default::default(),
            col_bfrt: Default::default(),
            col_dse: Default::default(),
            dev_row_ep: Default::default(),
            dev_col_dse: Default::default(),
            dual_row: Default::default(),
            dual_infeas_count: Default::default(),
            dual_rhs: Default::default(),
            row_out: Default::default(),
            variable_out: Default::default(),
            move_out: Default::default(),
            variable_in: Default::default(),
            delta_primal: Default::default(),
            theta_dual: Default::default(),
            theta_primal: Default::default(),
            alpha_col: Default::default(),
            alpha_row: Default::default(),
            numerical_trouble: Default::default(),
            computed_edge_weight: Default::default(),
            check_invert_condition: false,
            slice_num: Default::default(),
            slice_price: Default::default(),
            slice_start: Default::default(),
            slice_a_matrix: Default::default(),
            slice_ar_matrix: Default::default(),
            slice_row_ap: Default::default(),
            slice_dual_row: Default::default(),
            multi_num: Default::default(),
            multi_chosen: Default::default(),
            multi_i_choice: Default::default(),
            multi_n_finish: Default::default(),
            multi_iteration: Default::default(),
            multi_choose_again: Default::default(),
            multi_choice: Default::default(),
            multi_finish: Default::default(),
        }
    }
}

impl HEkkDual {
    /// HiGHS: HEkkDual.h::HEkkDual::HEkkDual(HEkk& simplex) (header constructor)
    ///
    /// C++: `ekk_instance_(simplex), dualRow(simplex), dualRHS(simplex) { initialiseInstance();
    /// dualRow.setup(); dualRHS.setup(); if (simplex_strategy != kSimplexStrategyDualPlain)
    /// initialiseInstanceParallel(simplex); }`
    pub fn new(ekk: &mut HEkk) -> Self {
        let mut dual = HEkkDual::default();
        // PORT NOTE (perf): the seven local vectors take recycled storage (HVectorPool):
        // each equals `HVector::default()`, so `initialise_instance` sets them up to exactly
        // the values fresh members get. `recycle` hands them back after the solve.
        dual.row_ep = ekk.hvector_pool_.take();
        dual.row_ap = ekk.hvector_pool_.take();
        dual.col_aq = ekk.hvector_pool_.take();
        dual.col_bfrt = ekk.hvector_pool_.take();
        dual.col_dse = ekk.hvector_pool_.take();
        dual.dev_row_ep = ekk.hvector_pool_.take();
        dual.dev_col_dse = ekk.hvector_pool_.take();
        dual.initialise_instance(ekk);
        dual.dual_row.setup(ekk);
        dual.dual_rhs.setup(ekk);
        if ekk.info_.simplex_strategy != kSimplexStrategyDualPlain {
            dual.initialise_instance_parallel(ekk);
        }
        dual
    }

    /// Not in HiGHS (perf only): returns the local vectors taken in `new` to
    /// `ekk.hvector_pool_` once the solver is finished with them. Dropping the solver
    /// instead (as C++ does at end of scope) is equally correct, only slower.
    pub fn recycle(mut self, ekk: &mut HEkk) {
        let pool = &mut ekk.hvector_pool_;
        pool.give(std::mem::take(&mut self.row_ep));
        pool.give(std::mem::take(&mut self.row_ap));
        pool.give(std::mem::take(&mut self.col_aq));
        pool.give(std::mem::take(&mut self.col_bfrt));
        pool.give(std::mem::take(&mut self.col_dse));
        pool.give(std::mem::take(&mut self.dev_row_ep));
        pool.give(std::mem::take(&mut self.dev_col_dse));
    }
}
