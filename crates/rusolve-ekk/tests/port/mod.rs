//! Port-side glue for the golden tests: builds `rusolve-ekk` values from
//! fixture records and views port values in the oracle's encoding.
//!
//! Unlike `tests/common/` (which is port-independent and checked by
//! `fixture_selftest`), this module compiles against the real port API. Use
//! it from an integration test with `mod common; mod port;`.
//!
//! Every builder mirrors what the oracle did in C++ to create its inputs
//! (`tools/highs-oracle/src/*.cpp`), so a replayed call sees the same state
//! HiGHS saw. Where the port API forced a choice, the comment says so.
#![allow(dead_code, unused_imports)]

use crate::common::events::{EkkScalars, EkkSolveEnter, SimplexStatus};
use crate::common::records::ScaleRec;
use crate::common::*;
use rusolve_ekk::prelude::*;
use serde_json::Value;

pub mod g1;
pub mod ipx;
pub mod sim;
pub mod w1;
pub mod w3;
pub mod w4;
pub mod w5;
pub use sim::*;

// ---------------------------------------------------------------------------
// Option profile
// ---------------------------------------------------------------------------

/// Env var: the option profile the replay builders start from. Unset or
/// `default`: [`HighsOptions::default`] (HiGHS 1.15.1 semantics, every `h16_*`
/// option off). `recommended`: [`HighsOptions::recommended`]. `all-h16`:
/// the default with every `h16_*` option on. The oracle records then override
/// what they carry, as before.
pub const TEST_PROFILE_ENV: &str = "RUSOLVE_TEST_PROFILE";

/// The base [`HighsOptions`] of [`TEST_PROFILE_ENV`].
pub fn profile_options() -> HighsOptions {
    match std::env::var(TEST_PROFILE_ENV).as_deref() {
        Err(_) | Ok("") | Ok("default") => HighsOptions::default(),
        Ok("recommended") => HighsOptions::recommended(),
        Ok("all-h16") => {
            let mut o = HighsOptions::default();
            o.set_all_h16(true);
            o
        }
        Ok(other) => panic!("{TEST_PROFILE_ENV}={other}: expected default, recommended or all-h16"),
    }
}

// ---------------------------------------------------------------------------
// Fixture sets
// ---------------------------------------------------------------------------

/// Every manifest fixture `<dir>/*.jsonl` that is present: committed ones
/// always, target-only ones when generated (otherwise a SKIP banner is
/// printed by `load`, or a panic with RUSOLVE_REQUIRE_TARGET_FIXTURES=1).
pub fn fixtures_in(dir: &str) -> Vec<Fixture> {
    let mut v: Vec<_> = manifest().select(&format!("{dir}/"), ".jsonl");
    v.sort_by_key(|e| (e.location != Location::Committed, e.file.clone()));
    v.into_iter().filter_map(|e| load(&e.file)).collect()
}

/// The records of type `ty`; panics if there are none, so a renamed record
/// type cannot turn a golden test into a vacuous pass.
pub fn recs_of<'a>(fx: &'a Fixture, ty: &str) -> Vec<Rec<'a>> {
    let v: Vec<_> = fx.of_type(ty).collect();
    assert!(!v.is_empty(), "{}: no {ty:?} records", fx.file);
    v
}

/// `HighsSparseMatrix::priceBy*` `debug_report` argument the oracle used
/// (the C++ default).
pub const DEBUG_REPORT_OFF: HighsInt = kDebugReportOff;

// ---------------------------------------------------------------------------
// Matrices and vectors
// ---------------------------------------------------------------------------

/// `MatrixFormat` from the oracle's integer encoding.
pub fn matrix_format(format: i32) -> MatrixFormat {
    match format {
        1 => MatrixFormat::kColwise,
        2 => MatrixFormat::kRowwise,
        3 => MatrixFormat::kRowwisePartitioned,
        f => panic!("unknown MatrixFormat {f}"),
    }
}

/// A `HighsSparseMatrix` holding exactly the fixture's arrays.
pub fn matrix_from_rec(m: &MatrixRec) -> HighsSparseMatrix {
    HighsSparseMatrix {
        format_: matrix_format(m.format),
        num_col_: m.num_col,
        num_row_: m.num_row,
        start_: m.start.clone(),
        p_end_: m.p_end.clone().unwrap_or_default(),
        index_: m.index.clone(),
        value_: m.value.clone(),
    }
}

/// View a port matrix in the oracle's encoding (`Obj::mat`).
pub fn matrix_view(m: &HighsSparseMatrix) -> MatrixView<'_> {
    MatrixView {
        format: m.format_ as i32,
        num_col: m.num_col_,
        num_row: m.num_row_,
        start: &m.start_,
        p_end: Some(&m.p_end_),
        index: &m.index_,
        value: &m.value_,
    }
}

/// An `HVector` as the oracle built its inputs: `setup(size)` (via
/// `makeVec`: setup + clear), then `count`, `index[0..count)`, `array`,
/// `synthetic_tick` and, when dumped, the `pack*` fields. `index` past
/// `count` is 0 (HiGHS leaves stale entries there; no HiGHS routine reads
/// them). A dense-mode record (`count < 0`) gets `count = -1`.
pub fn hvector_from_rec(r: &HVectorRec) -> HVector {
    let mut v = HVector::default();
    v.setup(r.size);
    v.clear();
    v.count = r.count;
    v.index = r.index_padded();
    v.array = r.array();
    v.synthetic_tick = r.synthetic_tick;
    if let Some(p) = &r.pack {
        v.pack_flag = p.pack_flag;
        v.pack_count = p.pack_count;
        v.pack_index[..p.pack_index.len()].copy_from_slice(&p.pack_index);
        v.pack_value[..p.pack_value.len()].copy_from_slice(&p.pack_value);
    }
    v
}

/// View a port HVector in the oracle's encoding (`Obj::hv`).
pub fn hv_view(v: &HVector, with_pack: bool) -> HVectorView<'_> {
    HVectorView {
        size: v.size,
        count: v.count,
        index: &v.index,
        array: &v.array,
        synthetic_tick: v.synthetic_tick,
        pack: with_pack.then_some(PackView {
            pack_flag: v.pack_flag,
            pack_count: v.pack_count,
            pack_index: &v.pack_index,
            pack_value: &v.pack_value,
        }),
    }
}

/// Compare an HVector with a fixture record (pack fields only if the record
/// has them).
pub fn check_hv(chk: &mut Checker, field: &str, expected: &HVectorRec, actual: &HVector) {
    chk.hv(field, expected, &hv_view(actual, expected.pack.is_some()));
}

/// Compare `HFactor::getInvert()` with a fixture `invert` object, array by
/// array (full vector length, as dumped).
pub fn check_invert(chk: &mut Checker, field: &str, e: &InvertRec, a: &InvertibleRepresentation) {
    macro_rules! ints {
        ($($k:ident),*) => {$( chk.ints(&format!("{field}.{}", stringify!($k)), &e.$k, &a.$k); )*};
    }
    macro_rules! f64s {
        ($($k:ident),*) => {$( chk.f64s(&format!("{field}.{}", stringify!($k)), &e.$k, &a.$k); )*};
    }
    ints!(
        l_pivot_index,
        l_pivot_lookup,
        l_start,
        l_index,
        lr_start,
        lr_index,
        u_pivot_lookup,
        u_pivot_index,
        u_start,
        u_last_p,
        u_index,
        ur_start,
        ur_lastp,
        ur_space,
        ur_index,
        pf_start,
        pf_index,
        pf_pivot_index
    );
    f64s!(
        l_value,
        lr_value,
        u_pivot_value,
        u_value,
        ur_value,
        pf_value,
        pf_pivot_value
    );
}

// ---------------------------------------------------------------------------
// HFactor exactly as the factor oracle sets it up
// ---------------------------------------------------------------------------

/// `setupFactor` in factor.cpp: `HFactor::setupGeneral(num_col, num_row,
/// num_row, a, basic_index, kDefaultPivotThreshold, kDefaultPivotTolerance,
/// kHighsDebugLevelMin, NULL, true, update_method)`.
///
/// Signature note: `setup_general_raw` takes `&HighsLogOptions` where the
/// C++ takes a nullable pointer (the oracle passes NULL). With NULL, HiGHS
/// sets output_flag = false, log_to_console = true, log_dev_level = 0; we pass
/// exactly those values. Logging never affects numerics.
pub fn setup_factor(f: &mut HFactor, a: &HighsSparseMatrix, update_method: HighsInt) {
    let log = null_log_options();
    f.setup_general_raw(
        a.num_col_,
        a.num_row_,
        a.num_row_,
        &a.start_,
        &a.index_,
        &a.value_,
        kDefaultPivotThreshold,
        kDefaultPivotTolerance,
        kHighsDebugLevelMin,
        &log,
        true,
        update_method,
    );
}

/// The log options HFactor::setupGeneral installs when given NULL.
pub fn null_log_options() -> HighsLogOptions {
    HighsLogOptions {
        log_stream: None,
        output_flag: false,
        log_to_console: true,
        log_dev_level: 0,
        user_callback: None,
        user_callback_active: false,
    }
}

/// `methodName` in factor.cpp, inverted.
pub fn update_method_from_name(name: &str) -> HighsInt {
    match name {
        "FT" => kUpdateMethodFt,
        "PF" => kUpdateMethodPf,
        "MPF" => kUpdateMethodMpf,
        "APF" => kUpdateMethodApf,
        m => panic!("unknown update method {m}"),
    }
}

/// Something that can replay the oracle's `solves` entries: `HFactor`
/// (`ftranCall` / `btranCall`, both overloads) or `HSimplexNla`.
pub trait Solver {
    /// The HiGHS function(s) under test, for messages.
    fn name(&self) -> &'static str;
    fn hv(&mut self, ftran: bool, v: &mut HVector, expected_density: f64);
    fn vec(&mut self, ftran: bool, v: &mut Vec<f64>);
}

impl Solver for HFactor {
    fn name(&self) -> &'static str {
        "HFactor::ftranCall/btranCall"
    }
    fn hv(&mut self, ftran: bool, v: &mut HVector, ed: f64) {
        if ftran {
            self.ftran_call(v, ed, None);
        } else {
            self.btran_call(v, ed, None);
        }
    }
    fn vec(&mut self, ftran: bool, v: &mut Vec<f64>) {
        if ftran {
            self.ftran_call_vec(v, None);
        } else {
            self.btran_call_vec(v, None);
        }
    }
}

/// Replay a `solves` list in order on `s` and compare every output.
/// Field names are `solves[i]<op rhs_kind ed=...>`.
pub fn check_solves(chk: &mut Checker, field: &str, solves: &[SolveRec], s: &mut dyn Solver) {
    for (i, sr) in solves.iter().enumerate() {
        let ftran = sr.is_ftran();
        match (&sr.input, &sr.out) {
            (SolveIo::Hv(input), SolveIo::Hv(out)) => {
                let ed = sr.expected_density.unwrap_or_else(|| {
                    panic!("{}: HVector solve without expected_density", sr.ctx)
                });
                let f = format!(
                    "{field}[{i}]<{} rhs_kind={} expected_density={ed} [{}]>",
                    sr.op,
                    sr.rhs_kind,
                    s.name()
                );
                let mut v = hvector_from_rec(input);
                s.hv(ftran, &mut v, ed);
                check_hv(chk, &f, out, &v);
            }
            (SolveIo::Vec(input), SolveIo::Vec(out)) => {
                let f = format!(
                    "{field}[{i}]<{} rhs_kind={} (std::vector overload) [{}]>",
                    sr.op,
                    sr.rhs_kind,
                    s.name()
                );
                let mut v = input.clone();
                s.vec(ftran, &mut v);
                chk.f64s(&f, out, &v);
            }
            _ => panic!("{}: solve in/out kinds differ", sr.ctx),
        }
    }
}

/// Compare the post-`build()` state of an HFactor with a `result` object
/// (`buildObj` in factor.cpp).
pub fn check_build_result(
    chk: &mut Checker,
    e: &BuildResultRec,
    f: &HFactor,
    rank_deficiency: HighsInt,
    basic_index_after: &[HighsInt],
) {
    chk.int("rank_deficiency", e.rank_deficiency, rank_deficiency);
    chk.ints("basic_index_after", &e.basic_index_after, basic_index_after);
    chk.ints(
        "row_with_no_pivot",
        &e.row_with_no_pivot,
        &f.row_with_no_pivot,
    );
    chk.ints(
        "col_with_no_pivot",
        &e.col_with_no_pivot,
        &f.col_with_no_pivot,
    );
    chk.ints(
        "var_with_no_pivot",
        &e.var_with_no_pivot,
        &f.var_with_no_pivot,
    );
    chk.f64(
        "build_synthetic_tick",
        e.build_synthetic_tick,
        f.build_synthetic_tick,
    );
    chk.int(
        "basis_matrix_num_el",
        e.basis_matrix_num_el,
        f.basis_matrix_num_el,
    );
    chk.int("invert_num_el", e.invert_num_el, f.invert_num_el);
    chk.int("kernel_dim", e.kernel_dim, f.kernel_dim);
    chk.int("kernel_num_el", e.kernel_num_el, f.kernel_num_el);
    chk.ints(
        "refactor_pivot_row",
        &e.refactor_pivot_row,
        &f.refactor_info_.pivot_row,
    );
    chk.ints(
        "refactor_pivot_var",
        &e.refactor_pivot_var,
        &f.refactor_info_.pivot_var,
    );
    chk.i8s(
        "refactor_pivot_type",
        &e.refactor_pivot_type,
        &f.refactor_info_.pivot_type,
    );
    chk.f64(
        "refactor_build_synthetic_tick",
        e.refactor_build_synthetic_tick,
        f.refactor_info_.build_synthetic_tick,
    );
    check_invert(chk, "invert", &e.invert, &f.get_invert());
}

/// The scaled constraint matrix of a factor fixture (`factor_lp.a`).
pub fn factor_matrix(fx: &Fixture) -> (HighsSparseMatrix, String) {
    let r = fx.first("factor_lp");
    (matrix_from_rec(&r.mat("a")), r.str("name").to_string())
}

/// `buildCase` in factor.cpp, first half: a fresh FT factor built on
/// `basic_index`. Returns (factor, rank_deficiency, basic_index after).
pub fn build_factor(
    a: &HighsSparseMatrix,
    basic_index: &[HighsInt],
) -> (HFactor, HighsInt, Vec<HighsInt>) {
    let mut f = HFactor::default();
    let mut bi = basic_index.to_vec();
    setup_factor(&mut f, a, kUpdateMethodFt);
    let rd = f.build(&mut bi, None);
    (f, rd, bi)
}

thread_local! {
    static AT: std::cell::RefCell<String> = const { std::cell::RefCell::new(String::new()) };
}

/// Record where a case is (e.g. "step k=3, HFactor::update"), so a panic
/// inside the port caught by [`Failures::run`] names the step.
pub fn at(progress: String) {
    AT.with(|a| *a.borrow_mut() = progress);
}

/// Collects the failures of independent cases (e.g. one per fixture record)
/// so one run reports every failing case, each with its own first
/// divergence, then panics once in [`Failures::finish`].
#[derive(Default)]
pub struct Failures {
    pub cases: usize,
    pub msgs: Vec<String>,
}

impl Failures {
    /// Record a finished checker; returns true when it passed.
    pub fn take(&mut self, chk: &Checker) -> bool {
        self.cases += 1;
        if chk.ok() {
            return true;
        }
        let mut s = format!("{}: {} mismatch(es)", chk.ctx, chk.failures.len());
        for m in chk.failures.iter().take(12) {
            s.push_str(&format!("\n  {m}"));
        }
        if chk.failures.len() > 12 {
            s.push_str(&format!("\n  ... {} more", chk.failures.len() - 12));
        }
        self.msgs.push(s);
        false
    }

    /// Run one case, turning a panic inside the port into a reported
    /// failure (with `ctx`) so the other cases still run.
    pub fn run(&mut self, ctx: &str, case: impl FnOnce() -> Option<Checker>) {
        at(String::new());
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(case));
        let progress = AT.with(|a| std::mem::take(&mut *a.borrow_mut()));
        let ctx = if progress.is_empty() {
            ctx.to_string()
        } else {
            format!("{ctx}, at {progress}")
        };
        let ctx = ctx.as_str();
        match r {
            Ok(None) => self.cases += 1,
            Ok(Some(chk)) => {
                self.take(&chk);
            }
            Err(e) => {
                let msg = e
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                    .unwrap_or_else(|| "<non-string panic>".into());
                self.fail(format!("{ctx}: PANIC inside the port: {msg}"));
            }
        }
    }

    pub fn fail(&mut self, msg: String) {
        self.cases += 1;
        self.msgs.push(msg);
    }

    #[track_caller]
    pub fn finish(&self, what: &str) {
        assert!(self.cases > 0, "{what}: no cases ran");
        if self.msgs.is_empty() {
            eprintln!("{what}: all {} case(s) passed", self.cases);
            return;
        }
        panic!(
            "{what}: {} of {} case(s) failed\n\n{}",
            self.msgs.len(),
            self.cases,
            self.msgs.join("\n\n")
        );
    }
}

// ---------------------------------------------------------------------------
// LP, options and the HApp::solveLpSimplex prefix
// ---------------------------------------------------------------------------

pub fn scale_from_rec(s: &ScaleRec) -> HighsScale {
    HighsScale {
        strategy: s.strategy,
        has_scaling: s.has_scaling,
        num_col: s.num_col,
        num_row: s.num_row,
        cost: s.cost,
        col: s.col.clone(),
        row: s.row.clone(),
    }
}

/// The HighsLp HEkk receives (`lp` of `ekk_solve_enter`): the scaled LP,
/// with names and integrality empty (the oracle drops integrality).
pub fn lp_from_rec(r: &LpRec) -> HighsLp {
    HighsLp {
        num_col_: r.num_col,
        num_row_: r.num_row,
        col_cost_: r.col_cost.clone(),
        col_lower_: r.col_lower.clone(),
        col_upper_: r.col_upper.clone(),
        row_lower_: r.row_lower.clone(),
        row_upper_: r.row_upper.clone(),
        a_matrix_: matrix_from_rec(&r.a),
        sense_: match r.sense {
            1 => ObjSense::kMinimize,
            -1 => ObjSense::kMaximize,
            s => panic!("unknown ObjSense {s}"),
        },
        offset_: r.offset,
        scale_: scale_from_rec(&r.scale),
        is_scaled_: r.is_scaled,
        is_moved_: r.is_moved,
        ..HighsLp::default()
    }
}

/// `LpRec` view of a port LP, for comparison with `ekk_solve_enter.lp`.
pub fn lp_rec_of(lp: &HighsLp) -> LpRec {
    LpRec {
        num_col: lp.num_col_,
        num_row: lp.num_row_,
        sense: lp.sense_ as i32,
        offset: lp.offset_,
        col_cost: lp.col_cost_.clone(),
        col_lower: lp.col_lower_.clone(),
        col_upper: lp.col_upper_.clone(),
        row_lower: lp.row_lower_.clone(),
        row_upper: lp.row_upper_.clone(),
        a: MatrixRec {
            format: lp.a_matrix_.format_ as i32,
            num_col: lp.a_matrix_.num_col_,
            num_row: lp.a_matrix_.num_row_,
            start: lp.a_matrix_.start_.clone(),
            p_end: (lp.a_matrix_.format_ == MatrixFormat::kRowwisePartitioned)
                .then(|| lp.a_matrix_.p_end_.clone()),
            index: lp.a_matrix_.index_.clone(),
            value: lp.a_matrix_.value_.clone(),
        },
        is_scaled: lp.is_scaled_,
        is_moved: lp.is_moved_,
        scale: ScaleRec {
            strategy: lp.scale_.strategy,
            has_scaling: lp.scale_.has_scaling,
            num_col: lp.scale_.num_col,
            num_row: lp.scale_.num_row,
            cost: lp.scale_.cost,
            col: lp.scale_.col.clone(),
            row: lp.scale_.row.clone(),
        },
    }
}

/// `HighsOptions` from the `options` object of `ekk_solve_enter` (every
/// HiGHS option record by name). Each port field is read from its HiGHS
/// option name, and a missing name panics, so the port cannot silently run
/// with a default the oracle did not use.
///
/// Signature note: the port's `HighsOptions` has no setter by name, so
/// the mapping is spelled out here; a field added to the port must be added
/// here too. `log_options` is derived as HiGHS does (pointers to
/// output_flag / log_to_console / log_dev_level) via `sync_log_options`.
pub fn options_from_json(v: &Value, ctx: &str) -> HighsOptions {
    let r = Rec::new(v, format!("{ctx}.options"));
    let mut o = profile_options();
    macro_rules! f {
        ($($k:ident),*) => {$( o.$k = r.f64(stringify!($k)); )*};
    }
    macro_rules! i {
        ($($k:ident),*) => {$( o.$k = r.int(stringify!($k)); )*};
    }
    macro_rules! b {
        ($($k:ident),*) => {$( o.$k = r.bool(stringify!($k)); )*};
    }
    f!(
        time_limit,
        objective_bound,
        objective_target,
        primal_feasibility_tolerance,
        dual_feasibility_tolerance,
        infinite_cost,
        infinite_bound,
        small_matrix_value,
        large_matrix_value,
        dual_simplex_cost_perturbation_multiplier,
        primal_simplex_bound_perturbation_multiplier,
        dual_simplex_pivot_growth_tolerance,
        dual_steepest_edge_weight_error_tolerance,
        dual_steepest_edge_weight_log_error_threshold,
        rebuild_refactor_solution_error_tolerance,
        factor_pivot_threshold,
        factor_pivot_tolerance
    );
    i!(
        threads,
        random_seed,
        highs_debug_level,
        highs_analysis_level,
        simplex_strategy,
        simplex_scale_strategy,
        simplex_dual_edge_weight_strategy,
        simplex_primal_edge_weight_strategy,
        simplex_iteration_limit,
        simplex_update_limit,
        simplex_min_concurrency,
        simplex_max_concurrency,
        simplex_price_strategy,
        simplex_unscaled_solution_strategy,
        max_dual_simplex_cleanup_level,
        max_dual_simplex_phase1_cleanup_level,
        allowed_cost_scale_factor,
        cost_scale_factor,
        log_dev_level
    );
    b!(
        allow_unbounded_or_infeasible,
        no_unnecessary_rebuild_refactor,
        less_infeasible_DSE_check,
        less_infeasible_DSE_choose_row,
        output_flag,
        log_to_console,
        timeless_log
    );
    o.parallel = r.str("parallel").to_string();
    o.presolve = r.str("presolve").to_string();
    o.solver = r.str("solver").to_string();
    o.sync_log_options();
    o
}

/// The HApp-level options the port does not model; the tests replay the
/// path these values select, so assert them.
pub fn assert_happ_options(v: &Value, ctx: &str) {
    let r = Rec::new(v, format!("{ctx}.options"));
    // kHighsOptionOff = -1: HApp::solveLpSimplex neither dualizes nor permutes.
    assert_eq!(
        r.int("simplex_dualize_strategy"),
        -1,
        "{ctx}: test replays the no-dualize path"
    );
    assert_eq!(
        r.int("simplex_permute_strategy"),
        -1,
        "{ctx}: test replays the no-permute path"
    );
    assert_eq!(
        r.int("simplex_crash_strategy"),
        0,
        "{ctx}: test replays the no-crash path"
    );
    assert!(r.bool("use_original_HFactor_logic"), "{ctx}");
}

/// The objects a `HighsLpSolverObject` borrows (members of `Highs`).
pub struct HighsParts {
    pub lp: HighsLp,
    pub basis: HighsBasis,
    pub solution: HighsSolution,
    pub info: HighsInfo,
    pub callback: HighsCallback,
    pub options: HighsOptions,
    pub timer: HighsTimer,
}

/// Fresh `Highs` members holding the fixture LP and options.
pub fn highs_parts(enter: &EkkSolveEnter, ctx: &str) -> HighsParts {
    HighsParts {
        lp: lp_from_rec(&enter.lp),
        basis: HighsBasis::default(),
        solution: HighsSolution::default(),
        info: HighsInfo::default(),
        callback: HighsCallback::default(),
        options: options_from_json(&enter.options, ctx),
        timer: HighsTimer::default(),
    }
}

/// HApp::solveLpSimplex up to (not including) `ekk_instance.solve()`, for a
/// fresh `Highs` with no basis, as in the traces (one solve call, no
/// dualize/permute, LP already scaled and stored in `ekk_solve_enter.lp`):
///
/// ```text
/// ekk_instance.iteration_count_ = highs_info.simplex_iteration_count; // 0
/// ekk_instance.initialiseSimplexStats();
/// ekk_instance.moveLp(solver_object);
/// ```
///
/// `iteration_limit`, when given, overrides `simplex_iteration_limit` (used
/// to stop HEkk::solve after exactly k iterations).
pub fn ekk_after_move_lp(
    enter: &EkkSolveEnter,
    ctx: &str,
    iteration_limit: Option<HighsInt>,
) -> (HEkk, HighsParts) {
    assert_happ_options(&enter.options, ctx);
    assert!(
        enter.basis.is_none() && !enter.status.has_basis,
        "{ctx}: the tests replay the logical-basis path (no basis at solve entry)"
    );
    let mut parts = highs_parts(enter, ctx);
    if let Some(k) = iteration_limit {
        parts.options.simplex_iteration_limit = k;
    }
    // HApp: ekk_instance.iteration_count_ = highs_info.simplex_iteration_count;
    let mut ekk = HEkk {
        iteration_count_: parts.info.simplex_iteration_count,
        ..HEkk::default()
    };
    ekk.initialise_simplex_stats();
    {
        let mut so = HighsLpSolverObject::new(
            &mut parts.lp,
            &mut parts.basis,
            &mut parts.solution,
            &mut parts.info,
            &mut parts.callback,
            &mut parts.options,
            &mut parts.timer,
        );
        ekk.move_lp(&mut so);
    }
    (ekk, parts)
}

/// The statements of HEkk::solve before the `kOracleEkkSolveEnter` hook
/// (HEkk.cpp lines 1002-1015): the state the `ekk_solve_enter` record shows.
pub fn ekk_solve_prefix(ekk: &mut HEkk) {
    ekk.debug_initialise();
    ekk.initialise_analysis();
    ekk.initialise_control();
    ekk.dual_simplex_cleanup_level_ = 0;
    ekk.dual_simplex_phase1_cleanup_level_ = 0;
    ekk.previous_iteration_cycling_detected = -kHighsIInf;
}

// ---------------------------------------------------------------------------
// HEkk state in the trace encoding
// ---------------------------------------------------------------------------

pub fn status_of(s: &HighsSimplexStatus) -> SimplexStatus {
    SimplexStatus {
        initialised_for_new_lp: s.initialised_for_new_lp,
        is_dualized: s.is_dualized,
        is_permuted: s.is_permuted,
        initialised_for_solve: s.initialised_for_solve,
        has_basis: s.has_basis,
        has_ar_matrix: s.has_ar_matrix,
        has_nla: s.has_nla,
        has_dual_steepest_edge_weights: s.has_dual_steepest_edge_weights,
        has_invert: s.has_invert,
        has_fresh_invert: s.has_fresh_invert,
        has_fresh_rebuild: s.has_fresh_rebuild,
        has_dual_objective_value: s.has_dual_objective_value,
        has_primal_objective_value: s.has_primal_objective_value,
    }
}

/// `ekkScalars` in simplex.cpp.
pub fn scalars_of(e: &HEkk) -> EkkScalars {
    let i = &e.info_;
    EkkScalars {
        iteration_count: e.iteration_count_,
        update_count: i.update_count,
        updated_dual_objective_value: i.updated_dual_objective_value,
        updated_primal_objective_value: i.updated_primal_objective_value,
        dual_objective_value: i.dual_objective_value,
        primal_objective_value: i.primal_objective_value,
        num_primal_infeasibilities: i.num_primal_infeasibilities,
        sum_primal_infeasibilities: i.sum_primal_infeasibilities,
        num_dual_infeasibilities: i.num_dual_infeasibilities,
        sum_dual_infeasibilities: i.sum_dual_infeasibilities,
        costs_perturbed: i.costs_perturbed,
        costs_shifted: i.costs_shifted,
        bounds_perturbed: i.bounds_perturbed,
        bounds_shifted: i.bounds_shifted,
        cost_perturbation_base: e.cost_perturbation_base_,
        cost_perturbation_max_abs_cost: e.cost_perturbation_max_abs_cost_,
        build_synthetic_tick: e.build_synthetic_tick_,
        total_synthetic_tick: e.total_synthetic_tick_,
        col_aq_density: i.col_aq_density,
        row_ep_density: i.row_ep_density,
        row_ap_density: i.row_ap_density,
        row_dse_density: i.row_dse_density,
        col_bfrt_density: i.col_bfrt_density,
        primal_col_density: i.primal_col_density,
        dual_col_density: i.dual_col_density,
        dual_simplex_cleanup_level: e.dual_simplex_cleanup_level_,
        simplex_in_scaled_space: e.simplex_in_scaled_space_,
    }
}

/// `workArrays` in simplex.cpp.
pub fn arrays_of(e: &HEkk) -> WorkArraysRec {
    let i = &e.info_;
    WorkArraysRec {
        work_cost: i.work_cost_.clone(),
        work_dual: i.work_dual_.clone(),
        work_shift: i.work_shift_.clone(),
        work_lower: i.work_lower_.clone(),
        work_upper: i.work_upper_.clone(),
        work_range: i.work_range_.clone(),
        work_value: i.work_value_.clone(),
        base_lower: i.base_lower_.clone(),
        base_upper: i.base_upper_.clone(),
        base_value: i.base_value_.clone(),
        dual_edge_weight: e.dual_edge_weight_.clone(),
    }
}

/// `basisObj` in simplex.cpp.
pub fn basis_of(b: &SimplexBasis) -> BasisRec {
    BasisRec {
        basic_index: b.basic_index_.clone(),
        nonbasic_flag: b.nonbasic_flag_.clone(),
        nonbasic_move: b.nonbasic_move_.clone(),
    }
}
