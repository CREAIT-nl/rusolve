//! IPX glue for the golden tests `ipx_model` and `ipx_linalg`: builds
//! `rusolve_ekk::ipm::ipx` values from the records written by
//! `tools/highs-oracle/src/ipx.cpp` and compares port state against them
//! (FORMAT.md "ipx-model / ipx-linalg").

use std::rc::Rc;

use crate::common::*;
use rusolve_ekk::ipm::ipx::basis::Basis;
use rusolve_ekk::ipm::ipx::control::Control;
use rusolve_ekk::ipm::ipx::iterate::Iterate;
use rusolve_ekk::ipm::ipx::lu_update::LuOracleStats;
use rusolve_ekk::ipm::ipx::model::Model;
use rusolve_ekk::ipm::ipx::sparse_matrix::SparseMatrix;
use rusolve_ekk::ipm::ipx::{Info, Int, Parameters, Vector};
use rusolve_ekk::lp_data::highs_callback::HighsCallback;

/// `ipx::Parameters` of a `params` record (ipx.cpp paramsObj); the C++
/// `logfile` is null and `log_options` unset (silent, as the oracle).
pub fn ipx_params(r: &Rec) -> Parameters {
    let profile = super::profile_options();
    Parameters {
        display: r.int("display"),
        logfile: None,
        print_interval: r.f64("print_interval"),
        time_limit: r.f64("time_limit"),
        analyse_basis_data: r.int("analyse_basis_data") != 0,
        dualize: r.int("dualize"),
        scale: r.int("scale"),
        ipm_maxiter: r.int("ipm_maxiter"),
        ipm_feasibility_tol: r.f64("ipm_feasibility_tol"),
        ipm_optimality_tol: r.f64("ipm_optimality_tol"),
        ipm_drop_primal: r.f64("ipm_drop_primal"),
        ipm_drop_dual: r.f64("ipm_drop_dual"),
        ipm_drop_dual_soft: false,
        kkt_tol: r.f64("kkt_tol"),
        crash_basis: r.int("crash_basis"),
        dependency_tol: r.f64("dependency_tol"),
        volume_tol: r.f64("volume_tol"),
        rows_per_slice: r.int("rows_per_slice"),
        maxskip_updates: r.int("maxskip_updates"),
        lu_kernel: r.int("lu_kernel"),
        lu_pivottol: r.f64("lu_pivottol"),
        run_crossover: r.int("run_crossover"),
        start_crossover_tol: r.f64("start_crossover_tol"),
        pfeasibility_tol: r.f64("pfeasibility_tol"),
        dfeasibility_tol: r.f64("dfeasibility_tol"),
        debug: r.int("debug"),
        switchiter: r.int("switchiter"),
        stop_at_switch: r.int("stop_at_switch"),
        update_heuristic: r.int("update_heuristic"),
        maxpasses: r.int("maxpasses"),
        run_centring: r.int("run_centring"),
        max_centring_steps: r.int("max_centring_steps"),
        centring_ratio_tolerance: r.f64("centring_ratio_tolerance"),
        centring_ratio_reduction: r.f64("centring_ratio_reduction"),
        centring_alpha_scaling: r.f64("centring_alpha_scaling"),
        bad_products_tolerance: r.int("bad_products_tolerance"),
        highs_logging: r.int("highs_logging") != 0,
        timeless_log: r.int("timeless_log") != 0,
        log_options: None,
        // `h16_*` flags: the oracle records carry none; they follow the test
        // profile (off by default, the base behaviour).
        h16_ipx_stall: profile.h16_ipx_stall,
        h16_ipx_scaled: profile.h16_ipx_scaled,
        h16_ipm_gap_offset: profile.h16_ipm_gap_offset,
        stall_to_crossover: false,
        crossover_ipm_crash: false,
        crossover_batch_push: false,
    }
}

/// The Control of LpSolver as IpxWrapper.cpp sets it up: the parameters
/// and a HighsCallback without user callback (`Control::InterruptCheck`
/// asserts that one is set; the NDEBUG oracle drives IPX without one, which
/// takes the same path).
pub fn ipx_control(fx: &Fixture) -> Control {
    let control = Control::new();
    control.set_parameters(&ipx_params(&fx.first("params")));
    control.callback(Some(HighsCallback::default()));
    control
}

/// `Model::Load` of the `input` record (the arrays of IpxWrapper.cpp
/// fillInIpxData); returns the errflag.
pub fn ipx_load_model(control: &Control, input: &Rec, model: &mut Model) -> Int {
    let (ap, ai, ax) = (input.ints("Ap"), input.ints("Ai"), input.f64s("Ax"));
    let (rhs, obj, lb, ub) = (
        input.f64s("rhs"),
        input.f64s("obj"),
        input.f64s("lb"),
        input.f64s("ub"),
    );
    let ct: Vec<u8> = input.ints("constr_type").iter().map(|&c| c as u8).collect();
    fn or_null<T>(v: &[T]) -> Option<&[T]> {
        if v.is_empty() { None } else { Some(v) }
    }
    let (opt_i, opt_f) = (or_null::<i32>, or_null::<f64>);
    model.load(
        control,
        input.int("num_constr"),
        input.int("num_var"),
        opt_i(&ap),
        opt_i(&ai),
        opt_f(&ax),
        opt_f(&rhs),
        or_null(&ct),
        input.f64("offset"),
        opt_f(&obj),
        opt_f(&lb),
        opt_f(&ub),
    )
}

/// `ipx::Info` of an info object; the wall-clock `time_*` fields the
/// oracle omits are 0 (they never feed back into the algorithm).
pub fn ipx_info(r: &Rec) -> Info {
    Info {
        status: r.int("status"),
        status_ipm: r.int("status_ipm"),
        status_crossover: r.int("status_crossover"),
        errflag: r.int("errflag"),
        num_var: r.int("num_var"),
        num_constr: r.int("num_constr"),
        num_entries: r.int("num_entries"),
        num_rows_solver: r.int("num_rows_solver"),
        num_cols_solver: r.int("num_cols_solver"),
        num_entries_solver: r.int("num_entries_solver"),
        dualized: r.int("dualized"),
        dense_cols: r.int("dense_cols"),
        centring_tried: r.int("centring_tried"),
        centring_success: r.int("centring_success"),
        dependent_rows: r.int("dependent_rows"),
        dependent_cols: r.int("dependent_cols"),
        rows_inconsistent: r.int("rows_inconsistent"),
        cols_inconsistent: r.int("cols_inconsistent"),
        primal_dropped: r.int("primal_dropped"),
        dual_dropped: r.int("dual_dropped"),
        abs_presidual: r.f64("abs_presidual"),
        abs_dresidual: r.f64("abs_dresidual"),
        rel_presidual: r.f64("rel_presidual"),
        rel_dresidual: r.f64("rel_dresidual"),
        pobjval: r.f64("pobjval"),
        dobjval: r.f64("dobjval"),
        rel_objgap: r.f64("rel_objgap"),
        complementarity: r.f64("complementarity"),
        normx: r.f64("normx"),
        normy: r.f64("normy"),
        normz: r.f64("normz"),
        objval: r.f64("objval"),
        primal_infeas: r.f64("primal_infeas"),
        dual_infeas: r.f64("dual_infeas"),
        iter: r.int("iter"),
        kktiter1: r.int("kktiter1"),
        kktiter2: r.int("kktiter2"),
        basis_repairs: r.int("basis_repairs"),
        updates_start: r.int("updates_start"),
        updates_ipm: r.int("updates_ipm"),
        updates_crossover: r.int("updates_crossover"),
        ftran_sparse: r.f64("ftran_sparse"),
        btran_sparse: r.f64("btran_sparse"),
        mean_fill: r.f64("mean_fill"),
        max_fill: r.f64("max_fill"),
        maxvol_updates: r.int("maxvol_updates"),
        maxvol_skipped: r.int("maxvol_skipped"),
        maxvol_passes: r.int("maxvol_passes"),
        tbl_nnz: r.int("tbl_nnz"),
        tbl_max: r.f64("tbl_max"),
        frobnorm_squared: r.f64("frobnorm_squared"),
        lambdamax: r.f64("lambdamax"),
        volume_increase: r.f64("volume_increase"),
        ..Info::default()
    }
}

/// Every non-time field of `ipx::Info`.
pub fn check_ipx_info(chk: &mut Checker, field: &str, e: &Rec, a: &Info) {
    let ei = ipx_info(e);
    macro_rules! ints {
        ($($f:ident)*) => { $( chk.int(&format!("{field}.{}", stringify!($f)), ei.$f, a.$f); )* };
    }
    macro_rules! dbls {
        ($($f:ident)*) => { $( chk.f64(&format!("{field}.{}", stringify!($f)), ei.$f, a.$f); )* };
    }
    ints!(status status_ipm status_crossover errflag num_var num_constr num_entries
        num_rows_solver num_cols_solver num_entries_solver dualized dense_cols centring_tried
        centring_success dependent_rows dependent_cols rows_inconsistent cols_inconsistent
        primal_dropped dual_dropped iter kktiter1 kktiter2 basis_repairs updates_start
        updates_ipm updates_crossover maxvol_updates maxvol_skipped maxvol_passes tbl_nnz);
    dbls!(abs_presidual abs_dresidual rel_presidual rel_dresidual pobjval dobjval rel_objgap
        complementarity normx normy normz objval primal_infeas dual_infeas ftran_sparse
        btran_sparse mean_fill max_fill tbl_max frobnorm_squared lambdamax volume_increase);
}

/// An `ipx::SparseMatrix` against matObj (dims, colptr, the `entries()`
/// leading part of rowidx/values).
pub fn check_ipx_matrix(chk: &mut Checker, field: &str, e: &Rec, a: &SparseMatrix) {
    chk.int(&format!("{field}.rows"), e.int("rows"), a.rows());
    chk.int(&format!("{field}.cols"), e.int("cols"), a.cols());
    let nz = a.entries() as usize;
    let cp = &a.colptr()[..(a.cols() as usize + 1).min(a.colptr().len())];
    chk.ints(&format!("{field}.colptr"), &e.ints("colptr"), cp);
    chk.ints(
        &format!("{field}.rowidx"),
        &e.ints("rowidx"),
        &a.rowidx()[..nz],
    );
    chk.f64s(
        &format!("{field}.values"),
        &e.f64s("values"),
        &a.values()[..nz],
    );
}

pub fn check_vec(chk: &mut Checker, field: &str, e: &Rec, key: &str, a: &[f64]) {
    chk.f64s(&format!("{field}.{key}"), &e.f64s(key), a);
}

/// A fresh iterate of `model` with the LpSolver tolerances of
/// lp_solver.cc 368-373.
pub fn ipx_new_iterate(control: &Control, model: &Rc<Model>) -> Iterate {
    let mut it = Iterate::new(model.clone());
    it.set_feasibility_tol(control.ipm_feasibility_tol());
    it.set_optimality_tol(control.ipm_optimality_tol());
    if control.run_crossover() != 0 {
        it.set_start_crossover_tol(control.start_crossover_tol());
    }
    it
}

/// An iterate in the state of an `iterate` record.
pub fn ipx_load_iterate(control: &Control, model: &Rc<Model>, e: &Rec) -> Iterate {
    let mut it = ipx_new_iterate(control, model);
    it.oracle_load(
        &e.f64s("x"),
        &e.f64s("xl"),
        &e.f64s("xu"),
        &e.f64s("y"),
        &e.f64s("zl"),
        &e.f64s("zu"),
        &e.ints("state"),
        e.bool("postprocessed"),
    );
    it
}

/// Everything iterateObj dumps: the state, then (through the evaluating
/// getters) the objectives, residuals and complementarity measures.
pub fn check_ipx_iterate(chk: &mut Checker, field: &str, e: &Rec, it: &Iterate) {
    check_vec(chk, field, e, "x", it.x());
    check_vec(chk, field, e, "xl", it.xl());
    check_vec(chk, field, e, "xu", it.xu());
    check_vec(chk, field, e, "y", it.y());
    check_vec(chk, field, e, "zl", it.zl());
    check_vec(chk, field, e, "zu", it.zu());
    chk.ints(
        &format!("{field}.state"),
        &e.ints("state"),
        &it.oracle_state(),
    );
    chk.f64(
        &format!("{field}.feasibility_tol"),
        e.f64("feasibility_tol"),
        it.feasibility_tol(),
    );
    chk.f64(
        &format!("{field}.optimality_tol"),
        e.f64("optimality_tol"),
        it.optimality_tol(),
    );
    chk.f64(
        &format!("{field}.start_crossover_tol"),
        e.f64("start_crossover_tol"),
        it.start_crossover_tol(),
    );
    chk.f64(
        &format!("{field}.bounds_measure"),
        e.f64("bounds_measure"),
        it.bounds_measure_,
    );
    chk.f64(
        &format!("{field}.costs_measure"),
        e.f64("costs_measure"),
        it.costs_measure_,
    );
    chk.bool(
        &format!("{field}.postprocessed"),
        e.bool("postprocessed"),
        it.oracle_postprocessed(),
    );
    let scalars: [(&str, f64); 10] = [
        ("pobjective", it.pobjective()),
        ("dobjective", it.dobjective()),
        ("pobjective_after_postproc", it.pobjective_after_postproc()),
        ("dobjective_after_postproc", it.dobjective_after_postproc()),
        ("presidual", it.presidual()),
        ("dresidual", it.dresidual()),
        ("complementarity", it.complementarity()),
        ("mu", it.mu()),
        ("mu_min", it.mu_min()),
        ("mu_max", it.mu_max()),
    ];
    for (k, v) in scalars {
        chk.f64(&format!("{field}.{k}"), e.f64(k), v);
    }
    chk.bool(
        &format!("{field}.feasible"),
        e.bool("feasible"),
        it.feasible(),
    );
    chk.bool(&format!("{field}.optimal"), e.bool("optimal"), it.optimal());
    chk.bool(
        &format!("{field}.term_crit_reached"),
        e.bool("term_crit_reached"),
        it.term_crit_reached(),
    );
    check_vec(chk, field, e, "rb", &it.rb());
    check_vec(chk, field, e, "rl", &it.rl());
    check_vec(chk, field, e, "ru", &it.ru());
    check_vec(chk, field, e, "rc", &it.rc());
    let n = it.x().len() as Int;
    let scaling: Vec<f64> = (0..n).map(|j| it.scaling_factor(j)).collect();
    chk.f64s(&format!("{field}.scaling"), &e.f64s("scaling"), &scaling);
}

/// Everything basisObj dumps, incl. the LU memory state (BasicLu: the
/// sizes after the Reallocate loop) and, when fresh, the LU factors.
pub fn check_ipx_basis(chk: &mut Checker, field: &str, e: &Rec, basis: &Basis) {
    {
        let v = basis.oracle_view();
        chk.ints(&format!("{field}.basis"), &e.ints("basis"), v.basis);
        chk.ints(
            &format!("{field}.map2basis"),
            &e.ints("map2basis"),
            v.map2basis,
        );
        chk.bool(&format!("{field}.fresh"), e.bool("fresh"), v.fresh);
        for (k, a) in [
            ("factorizations", v.factorizations),
            ("updates", v.updates),
            ("ftran", v.ftran),
            ("btran", v.btran),
            ("ftran_sparse", v.ftran_sparse),
            ("btran_sparse", v.btran_sparse),
        ] {
            chk.int(&format!("{field}.{k}"), e.int(k), a);
        }
        chk.f64s(
            &format!("{field}.fill_factors"),
            &e.f64s("fill_factors"),
            v.fill_factors,
        );
        chk.f64(
            &format!("{field}.sum_ftran_density"),
            e.f64("sum_ftran_density"),
            v.sum_ftran_density,
        );
        chk.f64(
            &format!("{field}.sum_btran_density"),
            e.f64("sum_btran_density"),
            v.sum_btran_density,
        );
        chk.f64(&format!("{field}.pivottol"), e.f64("pivottol"), v.pivottol);
        match (e.str("lu"), &v.lu) {
            (
                "basiclu",
                LuOracleStats::BasicLu {
                    lmem,
                    umem,
                    wmem,
                    istore_size,
                    xstore_size,
                    li_size,
                    ui_size,
                    wi_size,
                    fill_factor,
                },
            ) => {
                for (k, a) in [
                    ("lmem", *lmem as i64),
                    ("umem", *umem as i64),
                    ("wmem", *wmem as i64),
                    ("istore_size", *istore_size as i64),
                    ("xstore_size", *xstore_size as i64),
                    ("Li_size", *li_size as i64),
                    ("Ui_size", *ui_size as i64),
                    ("Wi_size", *wi_size as i64),
                ] {
                    chk.int(&format!("{field}.{k}"), e.i64(k), a);
                }
                chk.f64(
                    &format!("{field}.fill_factor"),
                    e.f64("fill_factor"),
                    *fill_factor,
                );
            }
            (
                "forrest_tomlin",
                LuOracleStats::ForrestTomlin {
                    ft_updates,
                    fill_factor,
                },
            ) => {
                chk.int(
                    &format!("{field}.ft_updates"),
                    e.i64("ft_updates"),
                    *ft_updates as i64,
                );
                chk.f64(
                    &format!("{field}.fill_factor"),
                    e.f64("fill_factor"),
                    *fill_factor,
                );
            }
            (k, _) => chk.push(
                &format!("{field}.lu"),
                k.to_string(),
                "other kernel".into(),
                String::new(),
            ),
        }
    }
    chk.f64(
        &format!("{field}.mean_fill"),
        e.f64("mean_fill"),
        basis.mean_fill(),
    );
    chk.f64(
        &format!("{field}.max_fill"),
        e.f64("max_fill"),
        basis.max_fill(),
    );
    if let Some(f) = e.opt_obj("factors") {
        if !basis.factorization_is_fresh() {
            return;
        }
        let m = basis.model().rows() as usize;
        let (mut l, mut u) = (SparseMatrix::new(), SparseMatrix::new());
        let (mut rowperm, mut colperm) = (vec![0; m], vec![0; m]);
        basis.get_lu_factors(&mut l, &mut u, &mut rowperm, &mut colperm);
        check_ipx_matrix(chk, &format!("{field}.factors.L"), &f.obj("L"), &l);
        check_ipx_matrix(chk, &format!("{field}.factors.U"), &f.obj("U"), &u);
        chk.ints(
            &format!("{field}.factors.rowperm"),
            &f.ints("rowperm"),
            &rowperm,
        );
        chk.ints(
            &format!("{field}.factors.colperm"),
            &f.ints("colperm"),
            &colperm,
        );
    }
}

pub fn vector_of(e: &Rec, key: &str) -> Vector {
    Vector::from_slice(&e.f64s(key))
}
