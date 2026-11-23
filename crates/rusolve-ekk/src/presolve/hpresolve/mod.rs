//! HiGHS source: presolve/HPresolve.h, presolve/HPresolve.cpp (split by line range),
//! presolve/HPresolveTest.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The LP presolve engine. One struct, [`HPresolve`]; its methods live in the
//! sibling files, one contiguous C++ line range each (named in each
//! file header). This file holds the struct, its `Default`, the
//! result types, the shared macros and the trivial inline members of
//! `HPresolve.h`.
//!
//! # Ownership of the C++ pointer members
//!
//! | C++ member                        | Rust field                      | C++ expression -> Rust |
//! |-----------------------------------|---------------------------------|------------------------|
//! | `HighsLp* model`                  | `model: HighsLp` (owned)        | `model->x` -> `self.model.x` |
//! | `const HighsOptions* options`     | `options: HighsOptions` (clone) | `options->x` -> `self.options.x` |
//! | `HighsTimer* timer`               | `timer: Option<HighsTimer>`     | `this->timer->read()` -> `self.timer.as_ref().unwrap().read()`; `assert(this->timer)` -> `debug_assert!(self.timer.is_some())` |
//! | `HighsMipSolver* mipsolver = nullptr` | `mipsolver: Option<Box<HighsMipSolver>>` (`None` on the LP path) | see "MIP mode" below |
//!
//! `ok_set_input` takes the LP **by value** and stores it in `self.model`;
//! `PresolveComponent::run` moves it in (`std::mem::take`) and back out
//! with [`HPresolve::take_model`] after `run`, so the caller sees every
//! in-place modification exactly as through the C++ pointer.
//!
//! # MIP mode
//!
//! `HighsMipSolverData::runMipPresolve` moves the whole solver in with
//! [`HPresolve::ok_set_input_mip`] (`Box<HighsMipSolver>`) and takes it back with
//! [`HPresolve::take_mip_solver`] after `run`. While it is in:
//!
//! * the presolved model (C++ `*mipsolver->model_` = `mipdata_->presolvedModel` =
//!   `*model`, one object) lives in `self.model`; `ms.mipdata_.presolved_model` is empty.
//!   So `mipsolver->model_->x`, `mipsolver->numCol()`, `mipsolver->variableType(c)`
//!   -> `self.model.x` / `self.model.num_col_` / `self.model.integrality_[c]`, never
//!   through `ms` (they would read the empty model).
//! * `mipsolver->orig_model_->x` -> `self.mip().orig_model_.x`; other members
//!   (`mipdata_->epsilon`, `options_mip_`, `submip`) -> `self.mip().x` /
//!   `self.mip_mut().x`. `mipsolver->timer_` is `self.timer` (a clone taken by
//!   `ok_set_input_mip`; the stub timer only reads the wall clock).
//! * A call into `crate::mip` that takes `ms` (`domain.propagate(ms)`,
//!   `cliquetable.x(ms, ..)`, `implications.add_vlb(ms, ..)`, `mipdata_->checkLimits()`)
//!   runs inside `with_mip_model!(self, ms => body)`: it takes the solver out of
//!   `self.mipsolver` and swaps `self.model` into `ms.mipdata_.presolved_model` for the
//!   body, so the callee sees the same model the C++ shares. In the body use `ms` (and
//!   `mip_model!(ms)` for the model), plain fields of `self` other than `model`, and
//!   no `self` method (`self.mipsolver` is `None` there, so shared code would take the
//!   LP branches) and no early `return` / `?`.
//!
//! `model->integrality_` is assigned `kContinuous` for every column in
//! `ok_set_input` (C++ line 99), so on the LP path integrality tests in shared code
//! take the continuous branches.
//!
//! # Result types and the checked-call macro
//!
//! * `HPresolve::Result` -> [`HPresolveResult`] (`Result::kOk` ->
//!   `HPresolveResult::kOk`). The name `Result` is not used because the
//!   prelude glob would shadow `std::result::Result`.
//! * `StatusResult` -> [`StatusResult`]: `StatusResult(false)` ->
//!   `StatusResult::from(false)`, `StatusResult(Result::kX)` ->
//!   `StatusResult::from(HPresolveResult::kX)`; `if (sr)` / `!sr` ->
//!   `sr.flag()` / `!sr.flag()`; `static_cast<Result>(sr)` -> `sr.result()`.
//! * `HPRESOLVE_CHECKED_CALL(x)` -> `hpresolve_checked_call!(x)`
//!   (C++ lines 39-43: returns the result from the enclosing function
//!   unless it is `kOk`). The macro is in textual scope in every sibling.
//!
//! # Field names (C++ -> Rust)
//!
//! `Avalue` -> `a_value`, `Arow` -> `a_row`, `Acol` -> `a_col`,
//! `Anext` -> `a_next`, `Aprev` -> `a_prev`, `ARleft` -> `a_r_left`,
//! `ARright` -> `a_r_right`; all-lowercase names stay (`colhead`,
//! `rowroot`, `rowsize`, `colsize`, `rowpositions`, `freeslots`); camelCase
//! is snake_cased (`rowsizeInteger` -> `rowsize_integer`, `implColLower` ->
//! `impl_col_lower`, `colLowerSource` -> `col_lower_source`,
//! `impliedRowBounds` -> `implied_row_bounds`, `eqiters` -> `eqiters`,
//! `numDeletedRows` -> `num_deleted_rows`, `presolve_status_` ->
//! `presolve_status_`, `analysis_` -> `analysis_`).
//!
//! # Containers (all keep the C++ iteration order)
//!
//! | C++                                         | Rust |
//! |---------------------------------------------|------|
//! | `std::set<HighsInt>`                        | `BTreeSet<HighsInt>` |
//! | `std::set<std::pair<HighsInt,HighsInt>> equations` | `BTreeSet<(HighsInt, HighsInt)>` (tuple order = pair order) |
//! | `std::vector<std::set<...>::iterator> eqiters` | `Vec<Option<(HighsInt, HighsInt)>>`: the *key* the iterator points at; `equations.end()` -> `None` |
//! | `std::set<double>` / `std::map<double, T>`  | `BTreeSet<OrderedF64>` / `BTreeMap<OrderedF64, T>` ([`OrderedF64`]) |
//! | `map.lower_bound(x)`                        | `map.range(OrderedF64(x)..).next()` |
//! | `std::unordered_multimap<u64, HighsInt>`    | [`HighsMultiBuckets`] (libc++ order; see its doc) |
//! | `HighsHashTable<K, V>`                      | [`HighsHashTable`] (faithful layout and iteration) |
//! | `std::function<bool(HighsInt,bool,bool)>`   | `&mut dyn FnMut(HighsInt, bool, bool) -> bool` |
//! | `pdqsort(b, e, comp)`                       | `pdqsort_by(&mut v[..], comp)` ([`pdqsort_by`], faithful branchy pdqsort; `comp` returns `bool` = C++ `<`) |
//! | `pdqsort(b, e)` (no comparator)             | `pdqsort(&mut v[..])` (`util::highs_sort::pdqsort`, default `<`) |
//! | `HighsHashHelpers::hash(HighsInt)`          | `hash_i32(x)` |
//! | `HighsHashHelpers::hash(std::make_pair(u32 a, u32 b))` | `hash_u32_pair(a, b)` |
//! | `HighsHashHelpers::sparse_combine(h, i, v)` | `sparse_combine_value(&mut h, i, v)`; `double_hash_code(x)` returns `u32`, pass `as u64` |
//! | `std::vector<T>::iterator` loops            | indices |
//!
//! `equations` iterator rules: `eqiters[row] = equations.emplace(a, b).first`
//! -> `self.equations.insert((a, b)); self.eqiters[row as usize] = Some((a, b));`;
//! `equations.erase(eqiters[row])` -> `self.equations.remove(&key)`;
//! `eqiters[row]->first` -> `key.0`; `eqiters.assign(n, equations.end())`
//! -> `self.eqiters = vec![None; n]`. A loop `auto eq = equations.begin();
//! ... ++eq` holds the current key: `begin()` ->
//! `self.equations.iter().next().copied()`, `++eq` ->
//! `self.equations.range((Excluded(k), Unbounded)).next().copied()`
//! (`std::ops::Bound::{Excluded, Unbounded}`).
//!
//! # Slices
//!
//! `HPresolve` implements [`HighsTripletStorage`]; the getters below return
//! cursors. `for (const HighsSliceNonzero& nz : getColumnVector(col))`
//! with a read-only body -> `for nz in self.get_column_vector(col).iter(self)`;
//! with a body that mutates `self` -> `let mut it =
//! self.get_column_vector(col); while let Some(nz) = it.next_nz(self) {..}`.
//! `nz.index()` / `nz.value()` are copies taken when the cursor advanced
//! (see `util/highs_matrix_slice.rs` for the one caveat). A
//! `const HighsMatrixSlice<F>&` parameter (`equalityRowAddition`,
//! `computeDynamism`) is a generic `S: HighsTripletSlice` cursor taken by
//! value; pass it on to the postsolve stack as `rowvector.iter(self)`.
//!
//! # Linear-sum bound arrays
//!
//! `impliedRowBounds.f(args)` where `f` reads the bound arrays (see
//! `util/highs_linear_sum_bounds.rs`) ->
//! `self.implied_row_bounds.f(&col_bound_arrays!(self), args)`, and
//! `impliedDualRowBounds.f(args)` ->
//! `self.implied_dual_row_bounds.f(&row_dual_bound_arrays!(self), args)`.
//! These are the arrays of the C++ `setBoundArrays` calls (lines 2576-2583):
//! `model->col_lower_/col_upper_`, `implColLower/Upper`,
//! `colLower/UpperSource`; and `rowDualLower/Upper`, `implRowDualLower/Upper`,
//! `rowDualLower/UpperSource`. Build them at the call (they then see the
//! current values, as the C++ pointers do).
//!
//! # Members passed to methods of `self`; lambdas
//!
//! * A member passed by reference to another `self` method (e.g.
//!   `toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
//!   model->a_matrix_.start_)`, `getRowPositions(row, rowpositions)`):
//!   take each vector out with `std::mem::take`, call, put it back (crate
//!   root rule (f)). Take only the vectors passed, not their parent struct.
//! * A C++ lambda capturing `[&]` that calls `self` methods: write it as a
//!   local closure with an explicit `this: &mut HPresolve` (or `&HPresolve`)
//!   first parameter plus explicit parameters for captured locals it
//!   mutates, or as a nested `fn` inside the body. That is a body-local
//!   change; signatures stay fixed.
//!
//! # Analysis, logging, timers
//!
//! `analysis_.startPresolveRuleLog(r)` -> `self.analysis_start_presolve_rule_log(r)`;
//! `analysis_.stopPresolveRuleLog(r)` -> `self.analysis_stop_presolve_rule_log(r)`;
//! other `analysis_.x` members are accessed directly (`self.analysis_.allow_rule_[r as usize]`).
//! `highsLogUser/highsLogDev(options->log_options, ...)` ->
//! `crate::highs_log_user!(&self.options.log_options, ...)` /
//! `crate::highs_log_dev!(...)`. `printf` -> `print!`.

#![allow(clippy::diverging_sub_expression)]

use crate::prelude::*;
use crate::presolve::hpresolve_clique_table::HPresolveCliqueTable;

/// HiGHS: HPresolve.cpp::HPRESOLVE_CHECKED_CALL (lines 39-43)
#[allow(unused_macros)]
macro_rules! hpresolve_checked_call {
    ($call:expr) => {{
        let __result: HPresolveResult = $call;
        if __result != HPresolveResult::kOk {
            return __result;
        }
    }};
}

/// Bound arrays of `impliedRowBounds` (HPresolve.cpp 2576-2578):
/// `self.implied_row_bounds.add(&col_bound_arrays!(self), row, col, val)`.
#[allow(unused_macros)]
macro_rules! col_bound_arrays {
    ($s:expr) => {
        HighsLinearSumBoundArrays {
            var_lower: &$s.model.col_lower_,
            var_upper: &$s.model.col_upper_,
            impl_var_lower: &$s.impl_col_lower,
            impl_var_upper: &$s.impl_col_upper,
            impl_var_lower_source: &$s.col_lower_source,
            impl_var_upper_source: &$s.col_upper_source,
        }
    };
}

/// Bound arrays of `impliedDualRowBounds` (HPresolve.cpp 2580-2583):
/// `self.implied_dual_row_bounds.add(&row_dual_bound_arrays!(self), col, row, val)`.
#[allow(unused_macros)]
macro_rules! row_dual_bound_arrays {
    ($s:expr) => {
        HighsLinearSumBoundArrays {
            var_lower: &$s.row_dual_lower,
            var_upper: &$s.row_dual_upper,
            impl_var_lower: &$s.impl_row_dual_lower,
            impl_var_upper: &$s.impl_row_dual_upper,
            impl_var_lower_source: &$s.row_dual_lower_source,
            impl_var_upper_source: &$s.row_dual_upper_source,
        }
    };
}

/// MIP mode: run `$body` with the solver out of `$self.mipsolver` as `$ms`
/// (`Box<HighsMipSolver>`) and the presolved model swapped from `$self.model` into
/// `$ms.mipdata_.presolved_model`, then undo both. The value of `$body` is returned.
/// See the module doc ("MIP mode") for what the body may touch.
///
/// `with_mip_model!(self, ms => detach!(ms.mipdata_.domains[0], |d| d.propagate(&mut ms)))`
#[allow(unused_macros)] // used by the MIP-only presolve methods (mip_probing.rs, mip_dominated.rs, mip_stubs.rs)
macro_rules! with_mip_model {
    ($self:ident, $ms:ident => $body:expr) => {{
        let mut $ms = $self
            .mipsolver
            .take()
            .expect("with_mip_model!: not in MIP mode");
        std::mem::swap(&mut $self.model, &mut $ms.mipdata_.presolved_model);
        let __result = $body;
        std::mem::swap(&mut $self.model, &mut $ms.mipdata_.presolved_model);
        $self.mipsolver = Some($ms);
        __result
    }};
}

/// As `with_mip_model!`, but without swapping the presolved model: `$ms` is a
/// `&mut HighsMipSolver` borrowed from `$self.mipsolver`. Only for bodies that never
/// read `ms.mipdata_.presolved_model` (nor `ms.model_`), e.g. the
/// `HighsImplications::addVLB/addVUB` calls of `extractVarBounds`, which read only
/// their own tables, `feastol` and the debug solution.
///
/// PORT NOTE (perf): the two `HighsLp` swaps of
/// `with_mip_model!` move ~1.7 KB per call; C++ HPresolve's `model` simply *is*
/// `&mipdata_->presolvedModel`, so it has no such cost.
macro_rules! with_mip_solver {
    ($self:ident, $ms:ident => $body:expr) => {{
        let $ms: &mut HighsMipSolver = $self
            .mipsolver
            .as_deref_mut()
            .expect("with_mip_solver!: not in MIP mode");
        $body
    }};
}

/// As `with_mip_model!`, and also swap the live `$pss: &mut HighsPostsolveStack` into
/// `$ms.mipdata_.post_solve_stack` for the body. In C++ the stack HPresolve runs on *is*
/// `mipsolver->mipdata_->postSolveStack` (HighsMipSolverData.cpp:819), and callees read it
/// there (`HighsImplications::rebuild`, `HighsCliqueTable::extractCliques`). The port moves
/// the stack out of `mipdata_` while presolve runs, so the bodies that reach it through `ms`
/// see it only inside this macro.
macro_rules! with_mip_model_pss {
    ($self:ident, $pss:ident, $ms:ident => $body:expr) => {{
        with_mip_model!($self, $ms => {
            std::mem::swap(&mut *$pss, &mut $ms.mipdata_.post_solve_stack);
            let __r = $body;
            std::mem::swap(&mut *$pss, &mut $ms.mipdata_.post_solve_stack);
            __r
        })
    }};
}

mod aggregator;
mod clique_update;
mod clique_varbound;
mod col_presolve;
mod col_stuffing;
mod csx;
mod dual_fixing;
mod eq_subst;
mod fourier_motzkin;
mod initial_sweep_h16;
mod links;
mod matrix;
mod mip_dominated;
mod mip_probing;
mod mip_stubs;
mod parallel;
mod presolve_h16;
mod presolve_loop;
pub mod presolve_prof;
mod row_bounds_h16;
mod row_presolve;
mod run;
mod setup;
mod shrink;
mod singletons;
mod sparsify;
mod xor_gf2;
mod zero_cost_singleton;

/// HiGHS: HPresolve.h::HPresolve::Result
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HPresolveResult {
    kOk,
    kPrimalInfeasible,
    kDualInfeasible,
    kStopped,
}

/// HiGHS: HPresolve.h::HPresolve::StatusResult
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StatusResult {
    my_flag: bool,
    my_result: HPresolveResult,
}

impl From<bool> for StatusResult {
    /// `StatusResult(bool flag) : my_flag(flag), my_result(Result::kOk)`
    fn from(flag: bool) -> Self {
        Self {
            my_flag: flag,
            my_result: HPresolveResult::kOk,
        }
    }
}

impl From<HPresolveResult> for StatusResult {
    /// `StatusResult(Result result) : my_result(result) { my_flag = (result == Result::kOk); }`
    fn from(result: HPresolveResult) -> Self {
        Self {
            my_flag: result == HPresolveResult::kOk,
            my_result: result,
        }
    }
}

impl StatusResult {
    /// `explicit operator bool()`
    pub fn flag(&self) -> bool {
        self.my_flag
    }

    /// `explicit operator Result()` (`static_cast<Result>(sr)`)
    pub fn result(&self) -> HPresolveResult {
        self.my_result
    }
}

/// HiGHS: HPresolve.h::HPresolve
///
/// C++ declares every member without an initialiser except
/// `mipsolver = nullptr`; the scalars are assigned by `okSetInput` / `run`
/// before use. Rust defaults: empty containers, scalars 0 / false,
/// `presolve_status_ = kNotSet` (overwritten first thing in `run`),
/// `reduction_limit = 0` (assigned in `ok_set_input`).
#[derive(Clone, Debug)]
pub struct HPresolve {
    // pointer to model and options that where presolved
    pub model: HighsLp,
    pub options: HighsOptions,
    pub timer: Option<HighsTimer>,
    pub mipsolver: Option<Box<HighsMipSolver>>,
    pub primal_feastol: f64,

    // triplet storage
    pub a_value: Vec<f64>,
    pub a_row: Vec<HighsInt>,
    pub a_col: Vec<HighsInt>,

    // linked list links for column based links for each nonzero
    pub colhead: Vec<HighsInt>,
    pub a_next: Vec<HighsInt>,
    pub a_prev: Vec<HighsInt>,

    // splay tree links for row based iteration and nonzero lookup
    pub rowroot: Vec<HighsInt>,
    pub a_r_left: Vec<HighsInt>,
    pub a_r_right: Vec<HighsInt>,

    // length of rows and columns
    pub rowsize: Vec<HighsInt>,
    pub rowsize_integer: Vec<HighsInt>,
    pub rowsize_impl_int: Vec<HighsInt>,
    /// HiGHS 1.16 (#3367): HPresolve.h `rowsizeSingleton` (line 71), the number of
    /// singleton columns in each row. Maintained by link/unlink in every mode (a pure
    /// counter); only the `h16_row_duals` parallel-rows scan reads it.
    pub rowsize_singleton: Vec<HighsInt>,
    pub colsize: Vec<HighsInt>,

    // vector to store the nonzero positions of a row
    pub rowpositions: Vec<HighsInt>,
    /// Rust-only: the C++ `rowpositions` buffer as `getStoredRow()`'s raw pointer sees it.
    /// `clear()` + `push_back` keep the old tail, so a slice taken before a nested
    /// `storeRow` of a shorter row reads the new prefix then the stale tail (this
    /// never shrinks). A nested row longer than the capacity reallocates in C++ (the
    /// outer slice then reads freed memory, UB); here it reads the new contents.
    pub rowpositions_mem: Vec<HighsInt>,

    // stack to reuse free slots
    pub freeslots: Vec<HighsInt>,

    // implied bounds on primal and dual variables and their origins
    pub impl_col_lower: Vec<f64>,
    pub impl_col_upper: Vec<f64>,
    pub col_lower_source: Vec<HighsInt>,
    pub col_upper_source: Vec<HighsInt>,
    pub row_dual_lower: Vec<f64>,
    pub row_dual_upper: Vec<f64>,
    pub impl_row_dual_lower: Vec<f64>,
    pub impl_row_dual_upper: Vec<f64>,
    pub row_dual_lower_source: Vec<HighsInt>,
    pub row_dual_upper_source: Vec<HighsInt>,
    pub col_impl_source_by_row: Vec<BTreeSet<HighsInt>>,
    pub impl_row_dual_source_by_col: Vec<BTreeSet<HighsInt>>,

    // implied bounds on values of primal and dual rows
    pub implied_row_bounds: HighsLinearSumBounds,
    pub implied_dual_row_bounds: HighsLinearSumBounds,

    pub changed_row_indices: Vec<HighsInt>,
    pub changed_row_flag: Vec<u8>,
    pub changed_col_indices: Vec<HighsInt>,
    pub changed_col_flag: Vec<u8>,

    pub substitution_opportunities: Vec<(HighsInt, HighsInt)>,

    /// C++ `std::unordered_map<HighsInt, HighsHashTree<std::pair<HighsInt,
    /// HighsInt>, double>> liftingOpportunities`. Only MIP probing
    /// (`liftingForProbing`) inserts; the LP path only reserves it and
    /// clears entries that are never present.
    ///
    /// PORT NOTE: `BTreeMap` (row order, no `RandomState`) instead of the libc++
    /// `unordered_map` bucket order. The value is `util::highs_hash_tree`, keyed by
    /// `(col, val)` as `HighsUInt` pairs (`HighsHashKey` is implemented for `(u32, u32)`).
    /// `clear()` = C++ `clear()`.
    pub lifting_opportunities: BTreeMap<HighsInt, HighsHashTree<(HighsUInt, HighsUInt), f64>>,

    // equation rows sorted by (size, index), and the key of each row's entry
    pub equations: BTreeSet<(HighsInt, HighsInt)>,
    pub eqiters: Vec<Option<(HighsInt, HighsInt)>>,

    pub shrink_problem_enabled: bool,
    pub reduction_limit: usize,

    // vectors storing singleton rows and columns
    pub singleton_rows: Vec<HighsInt>,
    pub singleton_columns: Vec<HighsInt>,

    // flags to mark rows/columns as deleted
    pub row_deleted: Vec<u8>,
    pub col_deleted: Vec<u8>,

    // flags to skip repeated single-equation handling on unchanged rows
    pub single_equation_checked: Vec<u8>,

    pub num_probes: Vec<u16>,

    pub probing_contingent: i64,
    pub probing_num_del_col: HighsInt,
    pub num_probed: HighsInt,

    // counters for number of deleted rows and columns
    pub num_deleted_rows: HighsInt,
    pub num_deleted_cols: HighsInt,

    // old problem sizes to compute percentage reductions
    pub old_num_col: HighsInt,
    pub old_num_row: HighsInt,
    pub probing_early_abort: bool,

    pub presolve_status_: HighsPresolveStatus,
    pub analysis_: HPresolveAnalysis,

    /// HiGHS 1.16 (#3342): HPresolve.h::may_require_basis_postsolve_ /
    /// may_require_dual_postsolve_, set by chooseRules (h16_presolve_rules only).
    pub h16_may_require_basis_postsolve_: bool,
    pub h16_may_require_dual_postsolve_: bool,

    /// HiGHS latest (#3312, 7b890c255): HPresolve.h::HPresolve::presolveCliqueTable (line 44).
    /// HiGHS #3312 (mip_presolve_clique_table): inert (empty) unless the option is on.
    pub presolve_clique_table: HPresolveCliqueTable,
}

impl Default for HPresolve {
    fn default() -> Self {
        Self {
            model: HighsLp::default(),
            options: HighsOptions::default(),
            timer: None,
            mipsolver: None,
            primal_feastol: 0.0,
            a_value: Vec::new(),
            a_row: Vec::new(),
            a_col: Vec::new(),
            colhead: Vec::new(),
            a_next: Vec::new(),
            a_prev: Vec::new(),
            rowroot: Vec::new(),
            a_r_left: Vec::new(),
            a_r_right: Vec::new(),
            rowsize: Vec::new(),
            rowsize_integer: Vec::new(),
            rowsize_impl_int: Vec::new(),
            rowsize_singleton: Vec::new(),
            colsize: Vec::new(),
            rowpositions: Vec::new(),
            rowpositions_mem: Vec::new(),
            freeslots: Vec::new(),
            impl_col_lower: Vec::new(),
            impl_col_upper: Vec::new(),
            col_lower_source: Vec::new(),
            col_upper_source: Vec::new(),
            row_dual_lower: Vec::new(),
            row_dual_upper: Vec::new(),
            impl_row_dual_lower: Vec::new(),
            impl_row_dual_upper: Vec::new(),
            row_dual_lower_source: Vec::new(),
            row_dual_upper_source: Vec::new(),
            col_impl_source_by_row: Vec::new(),
            impl_row_dual_source_by_col: Vec::new(),
            implied_row_bounds: HighsLinearSumBounds::default(),
            implied_dual_row_bounds: HighsLinearSumBounds::default(),
            changed_row_indices: Vec::new(),
            changed_row_flag: Vec::new(),
            changed_col_indices: Vec::new(),
            changed_col_flag: Vec::new(),
            substitution_opportunities: Vec::new(),
            lifting_opportunities: BTreeMap::new(),
            equations: BTreeSet::new(),
            eqiters: Vec::new(),
            shrink_problem_enabled: false,
            reduction_limit: 0,
            singleton_rows: Vec::new(),
            singleton_columns: Vec::new(),
            row_deleted: Vec::new(),
            col_deleted: Vec::new(),
            single_equation_checked: Vec::new(),
            num_probes: Vec::new(),
            probing_contingent: 0,
            probing_num_del_col: 0,
            num_probed: 0,
            num_deleted_rows: 0,
            num_deleted_cols: 0,
            old_num_col: 0,
            old_num_row: 0,
            probing_early_abort: false,
            presolve_status_: HighsPresolveStatus::kNotSet,
            analysis_: HPresolveAnalysis::default(),
            h16_may_require_basis_postsolve_: false,
            h16_may_require_dual_postsolve_: false,
            presolve_clique_table: HPresolveCliqueTable::default(),
        }
    }
}

impl HighsTripletStorage for HPresolve {
    fn triplet_arrays(&self) -> HighsTripletArrays<'_> {
        HighsTripletArrays {
            a_value: &self.a_value,
            a_row: &self.a_row,
            a_col: &self.a_col,
            a_next: &self.a_next,
            a_r_left: &self.a_r_left,
            a_r_right: &self.a_r_right,
            row_positions: &self.rowpositions_mem,
        }
    }
}

impl HPresolve {
    /// Rust-only: move the presolved model out (the C++ caller keeps
    /// observing it through the `HighsLp*` it passed to `okSetInput`).
    pub fn take_model(&mut self) -> HighsLp {
        std::mem::take(&mut self.model)
    }

    /// HiGHS: HPresolve.h::HPresolve::clearLiftingOpportunities (lines 250-253, 4 lines)
    pub fn clear_lifting_opportunities(&mut self, row: HighsInt) {
        if let Some(search) = self.lifting_opportunities.get_mut(&row) {
            search.clear();
        }
    }

    /// HiGHS: HPresolve.h::HPresolve::setReductionLimit (lines 371-373, 3 lines)
    pub fn set_reduction_limit(&mut self, reduction_limit: usize) {
        self.reduction_limit = reduction_limit;
    }

    /// HiGHS: HPresolve.h::HPresolve::numNonzeros (line 375, 1 line)
    pub fn num_nonzeros(&self) -> HighsInt {
        (self.a_value.len() - self.freeslots.len()) as HighsInt
    }

    /// HiGHS: HPresolve.h::HPresolve::getPresolveLog (lines 493-495, 3 lines)
    pub fn get_presolve_log(&self) -> &HighsPresolveLog {
        &self.analysis_.presolve_log_
    }

    /// HiGHS: HPresolve.h::HPresolve::getPresolveStatus (line 497, 1 line)
    pub fn get_presolve_status(&self) -> HighsPresolveStatus {
        self.presolve_status_
    }

    /// HiGHS: HPresolve.cpp::HPresolve::getColumnVector (lines 2134-2137, 4 lines)
    pub fn get_column_vector(&self, col: HighsInt) -> HighsTripletListSlice {
        HighsTripletListSlice::new(self.colhead[col as usize])
    }

    /// HiGHS: HPresolve.cpp::HPresolve::getRowVector (lines 2139-2142, 4 lines)
    pub fn get_row_vector(&self, row: HighsInt) -> HighsTripletTreeSlicePreOrder {
        HighsTripletTreeSlicePreOrder::new(self.rowroot[row as usize])
    }

    /// HiGHS: HPresolve.cpp::HPresolve::getSortedRowVector (lines 2144-2147, 4 lines)
    pub fn get_sorted_row_vector(&self, row: HighsInt) -> HighsTripletTreeSliceInOrder {
        HighsTripletTreeSliceInOrder::new(self.rowroot[row as usize])
    }

    /// HiGHS: HPresolve.cpp::HPresolve::storeRow (line 2552, 1 line)
    pub fn store_row(&mut self, row: HighsInt) {
        let mut rowpositions = std::mem::take(&mut self.rowpositions);
        self.get_row_positions(row, &mut rowpositions);
        self.rowpositions = rowpositions;
        self.sync_rowpositions_mem();
    }

    /// Rust-only: write `rowpositions` over the front of [`HPresolve::rowpositions_mem`].
    pub fn sync_rowpositions_mem(&mut self) {
        let n = self.rowpositions.len();
        if self.rowpositions_mem.len() < n {
            self.rowpositions_mem.resize(n, 0);
        }
        self.rowpositions_mem[..n].copy_from_slice(&self.rowpositions);
    }

    /// HiGHS: HPresolve.cpp::HPresolve::getStoredRow (lines 2554-2557, 4 lines)
    pub fn get_stored_row(&self) -> HighsTripletPositionSlice {
        HighsTripletPositionSlice::new(self.rowpositions.len())
    }

    /// `analysis_.startPresolveRuleLog(rule_type)`: passes the counters and
    /// model that C++ `HPresolveAnalysis` reads through its stored pointers.
    pub fn analysis_start_presolve_rule_log(&mut self, rule_type: HighsInt) {
        self.analysis_.start_presolve_rule_log(
            rule_type,
            self.num_deleted_rows,
            self.num_deleted_cols,
            &self.model,
        );
    }

    /// `analysis_.stopPresolveRuleLog(rule_type)`; see
    /// [`HPresolve::analysis_start_presolve_rule_log`].
    pub fn analysis_stop_presolve_rule_log(&mut self, rule_type: HighsInt) {
        self.analysis_.stop_presolve_rule_log(
            rule_type,
            self.num_deleted_rows,
            self.num_deleted_cols,
            &self.model,
        );
    }
}

#[cfg(test)]
mod stored_row_tests {
    use super::*;
    use crate::util::highs_matrix_slice::HighsTripletSlice;

    /// rows 0 = x0 + x1 + x2 (3 nonzeros) and 1 = x0 (1 nonzero), colwise
    fn two_row_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![1.0; 3];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![1.0; 3];
        lp.row_lower_ = vec![-kHighsInf; 2];
        lp.row_upper_ = vec![2.0, 1.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 2, 3, 4];
        lp.a_matrix_.index_ = vec![0, 1, 0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 1.0];
        lp
    }

    /// C++ `getStoredRow()` holds `rowpositions.data()` and the size at creation: a
    /// nested `storeRow` of a shorter row leaves the outer slice reading the new
    /// prefix and then the stale tail (air05 seed 3 hit this in dominatedColumns).
    #[test]
    fn stored_row_slice_reads_stale_tail_after_nested_store() {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(
            two_row_lp(),
            &options,
            options.presolve_reduction_limit,
            None
        ));
        presolve.store_row(0);
        let row0 = presolve.rowpositions.clone();
        assert_eq!(row0.len(), 3);
        let mut it = presolve.get_stored_row();
        let first = it.next_nz(&presolve).unwrap();
        assert_eq!(first.position(), row0[0]);
        presolve.store_row(1);
        let row1 = presolve.rowpositions.clone();
        assert_eq!(row1.len(), 1);
        let rest: Vec<HighsInt> =
            std::iter::from_fn(|| it.next_nz(&presolve).map(|nz| nz.position())).collect();
        assert_eq!(rest, vec![row0[1], row0[2]]);
        // a fresh slice sees only the stored row
        let mut it = presolve.get_stored_row();
        assert_eq!(it.next_nz(&presolve).unwrap().position(), row1[0]);
        assert!(it.next_nz(&presolve).is_none());
    }
}
