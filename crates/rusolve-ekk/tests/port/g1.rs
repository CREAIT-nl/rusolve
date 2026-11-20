//! Oracle glue for the gap-fill API around `Highs::run`): fixture
//! discovery, the base options with the `overrides`, the `g1_header` model
//! and the replay of a fixture's `g1_call` sequence, each call checked
//! against the oracle's return status and outputs, each `g1_result` against
//! the port state. Used by `tests/t_g1_*.rs`; record layout in
//! tools/highs-oracle/FORMAT.md ("G1 records").

use crate::common::records::LpRec;
use crate::common::*;
use crate::port::w1::wave_fixtures;
use crate::port::*;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;

/// Every `g1/<dir>/*.jsonl` fixture available. Empty after a SKIP banner
/// ("SKIPPED, NOT PASSED") when none was generated, unless
/// RUSOLVE_REQUIRE_TARGET_FIXTURES=1.
pub fn g1_fixtures(dir: &str) -> Vec<Fixture> {
    wave_fixtures(&format!("g1/{dir}/"), "gap-fill oracle commands")
}

/// The oracle's base options (output off, serial, `random_seed=0`; presolve
/// is left at its default unless overridden, unlike the PDLP glue) plus the
/// `overrides` of the `g1_header` (`k=v,k=v`).
pub fn g1_options(overrides: &str) -> HighsOptions {
    let mut o = HighsOptions {
        output_flag: false,
        threads: 1,
        parallel: "off".to_string(),
        random_seed: 0,
        ..crate::port::profile_options()
    };
    for kv in overrides.split(',') {
        if kv == "default" || kv.is_empty() {
            continue;
        }
        let (k, v) = kv.split_once('=').expect("key=value override");
        let i = || v.parse::<HighsInt>().unwrap();
        let f = || v.parse::<f64>().unwrap();
        let b = || match v {
            "true" => true,
            "false" => false,
            _ => panic!("bool override {k}={v}"),
        };
        match k {
            "solver" => o.solver = v.to_string(),
            "presolve" => o.presolve = v.to_string(),
            "run_crossover" => o.run_crossover = v.to_string(),
            "user_objective_scale" => o.user_objective_scale = i(),
            "user_bound_scale" => o.user_bound_scale = i(),
            "infinite_cost" => o.infinite_cost = f(),
            "solve_relaxation" => o.solve_relaxation = b(),
            "blend_multi_objectives" => o.blend_multi_objectives = b(),
            "use_warm_start" => o.use_warm_start = b(),
            "mip_improving_solution_save" => o.mip_improving_solution_save = b(),
            "mip_rel_gap" => o.mip_rel_gap = f(),
            "mip_max_nodes" => o.mip_max_nodes = i(),
            "simplex_strategy" => o.simplex_strategy = i(),
            "mip_improving_solution_report_sparse" => o.mip_improving_solution_report_sparse = b(),
            // The oracle's file name, written under the system temp dir here
            "mip_improving_solution_file" => o.mip_improving_solution_file = g1_temp_path(v),
            other => panic!("unknown gap-fill override {other}"),
        }
    }
    o.sync_log_options();
    o
}

/// Where the port writes a file the oracle wrote as `name` (the
/// `mip_improving_solution_file` override).
pub fn g1_temp_path(name: &str) -> String {
    let dir = std::env::temp_dir().join(format!("rusolve-g1-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name).to_string_lossy().to_string()
}

/// The `g1_header` record (first line of every gap-fill fixture).
pub fn g1_header(fx: &Fixture) -> Rec<'_> {
    let h = fx.rec(0);
    assert_eq!(h.t(), "g1_header", "{}: first record", fx.file);
    h
}

/// `HighsVarType` from its C++ `uint8_t` value.
pub fn var_type(v: i32) -> HighsVarType {
    match v {
        0 => HighsVarType::kContinuous,
        1 => HighsVarType::kInteger,
        2 => HighsVarType::kSemiContinuous,
        3 => HighsVarType::kSemiInteger,
        4 => HighsVarType::kImplicitInteger,
        t => panic!("unknown HighsVarType {t}"),
    }
}

/// An `lp` object (`LpRec`, unscaled) plus its optional `integrality`.
pub fn g1_lp_from(r: &Rec) -> HighsLp {
    let mut lp = lp_from_rec(&LpRec::from_rec(&r.obj("lp")));
    if r.has("integrality") {
        lp.integrality_ = r.ints("integrality").into_iter().map(var_type).collect();
    }
    lp
}

/// The model the oracle passed to `Highs::passModel` (as read from the
/// `.mps` / `.lp` file, then any header edits).
pub fn g1_lp(fx: &Fixture) -> HighsLp {
    g1_lp_from(&g1_header(fx))
}

/// A `HighsSolution` from the optional `col_value`, `row_value`,
/// `col_dual`, `row_dual` of a `g1_call` (what the oracle built in C++).
pub fn solution_from(r: &Rec) -> HighsSolution {
    let mut s = HighsSolution::default();
    if r.has("col_value") {
        s.col_value = r.f64s("col_value");
        s.value_valid = true;
    }
    if r.has("row_value") {
        s.row_value = r.f64s("row_value");
    }
    if r.has("col_dual") {
        s.col_dual = r.f64s("col_dual");
        s.row_dual = r.f64s("row_dual");
        s.dual_valid = true;
    }
    s
}

/// The `objectives` list of a `pass_linear_objectives` call.
pub fn linear_objectives_from(r: &Rec) -> Vec<HighsLinearObjective> {
    r.list("objectives")
        .iter()
        .map(|o| HighsLinearObjective {
            weight: o.f64("weight"),
            offset: o.f64("offset"),
            coefficients: o.f64s("coefficients"),
            abs_tolerance: o.f64("abs_tolerance"),
            rel_tolerance: o.f64("rel_tolerance"),
            priority: o.int("priority"),
        })
        .collect()
}

/// Compares two LPs field by field (`lp` objects: the presolved LP, the
/// fixed LP).
pub fn check_lp(chk: &mut Checker, what: &str, want: &HighsLp, got: &HighsLp) {
    let mut got = got.clone();
    got.ensure_colwise();
    chk.int(&format!("{what}.num_col"), want.num_col_, got.num_col_)
        .int(&format!("{what}.num_row"), want.num_row_, got.num_row_)
        .int(
            &format!("{what}.sense"),
            want.sense_ as i32,
            got.sense_ as i32,
        )
        .f64(&format!("{what}.offset"), want.offset_, got.offset_)
        .f64s(&format!("{what}.col_cost"), &want.col_cost_, &got.col_cost_)
        .f64s(
            &format!("{what}.col_lower"),
            &want.col_lower_,
            &got.col_lower_,
        )
        .f64s(
            &format!("{what}.col_upper"),
            &want.col_upper_,
            &got.col_upper_,
        )
        .f64s(
            &format!("{what}.row_lower"),
            &want.row_lower_,
            &got.row_lower_,
        )
        .f64s(
            &format!("{what}.row_upper"),
            &want.row_upper_,
            &got.row_upper_,
        )
        .ints(
            &format!("{what}.a.start"),
            &want.a_matrix_.start_,
            &got.a_matrix_.start_,
        )
        .ints(
            &format!("{what}.a.index"),
            &want.a_matrix_.index_,
            &got.a_matrix_.index_,
        )
        .f64s(
            &format!("{what}.a.value"),
            &want.a_matrix_.value_,
            &got.a_matrix_.value_,
        );
    let int = |lp: &HighsLp| {
        lp.integrality_
            .iter()
            .map(|t| *t as i32)
            .collect::<Vec<i32>>()
    };
    chk.ints(&format!("{what}.integrality"), &int(want), &int(&got));
}

/// The `g1_result` fields (FORMAT.md, "G1 records") against the port after
/// the preceding solve call (`run`, `postsolve`, `crossover`). Optional
/// fields are checked when present.
pub fn check_g1_result(chk: &mut Checker, r: &Rec, highs: &Highs) {
    let info = highs.get_info();
    let sol = highs.get_solution();
    chk.int(
        "model_status",
        r.int("model_status"),
        highs.get_model_status() as i32,
    )
    .f64(
        "objective_function_value",
        r.f64("objective_function_value"),
        info.objective_function_value,
    )
    .int(
        "primal_solution_status",
        r.int("primal_solution_status"),
        info.primal_solution_status,
    )
    .int(
        "dual_solution_status",
        r.int("dual_solution_status"),
        info.dual_solution_status,
    )
    .int(
        "simplex_iteration_count",
        r.int("simplex_iteration_count"),
        info.simplex_iteration_count,
    )
    .f64(
        "max_primal_infeasibility",
        r.f64("max_primal_infeasibility"),
        info.max_primal_infeasibility,
    )
    .f64(
        "max_dual_infeasibility",
        r.f64("max_dual_infeasibility"),
        info.max_dual_infeasibility,
    )
    .bool("value_valid", r.bool("value_valid"), sol.value_valid)
    .bool("dual_valid", r.bool("dual_valid"), sol.dual_valid)
    .bool(
        "basis_valid",
        r.bool("basis_valid"),
        highs.get_basis().valid,
    )
    .f64s("col_value", &r.f64s("col_value"), &sol.col_value)
    .f64s("row_value", &r.f64s("row_value"), &sol.row_value)
    .f64s("col_dual", &r.f64s("col_dual"), &sol.col_dual)
    .f64s("row_dual", &r.f64s("row_dual"), &sol.row_dual);
    if r.has("mip_node_count") {
        chk.int(
            "mip_node_count",
            r.i64("mip_node_count"),
            info.mip_node_count,
        )
        .f64(
            "mip_dual_bound",
            r.f64("mip_dual_bound"),
            info.mip_dual_bound,
        )
        .f64("mip_gap", r.f64("mip_gap"), info.mip_gap)
        .f64(
            "max_integrality_violation",
            r.f64("max_integrality_violation"),
            info.max_integrality_violation,
        );
    }
    if r.has("col_status") {
        let basis = highs.get_basis();
        let st = |v: &[HighsBasisStatus]| v.iter().map(|s| *s as i32).collect::<Vec<i32>>();
        chk.ints("col_status", &r.ints("col_status"), &st(&basis.col_status))
            .ints("row_status", &r.ints("row_status"), &st(&basis.row_status));
    }
    if r.has("lp") {
        check_lp(chk, "model", &g1_lp_from(r), highs.get_lp());
    }
}

/// Executes one `g1_call` (`op` plus arguments) on the port, checks the
/// returned status against `status` and any op outputs against the record.
pub fn apply_g1_call(chk: &mut Checker, highs: &mut Highs, rc: &Rec) {
    let op = rc.str("op");
    let f = |k: &str| rc.f64(k);
    let i = |k: &str| rc.int(k);
    let status = match op {
        "run" => highs.run(),
        "set_solution" => highs.set_solution(&solution_from(rc)),
        "set_solution_sparse" => {
            let index = rc.ints("index");
            highs.set_solution_sparse(index.len() as HighsInt, &index, &rc.f64s("value"))
        }
        "pass_linear_objectives" => {
            let objs = linear_objectives_from(rc);
            highs.pass_linear_objectives(objs.len() as HighsInt, &objs)
        }
        "clear_linear_objectives" => highs.clear_linear_objectives(),
        "presolve" => {
            let s = highs.presolve();
            chk.int(
                "model_presolve_status",
                rc.int("model_presolve_status"),
                highs.get_model_presolve_status() as i32,
            );
            if rc.has("presolved") {
                check_lp(
                    chk,
                    "presolved",
                    &g1_lp_from(&rc.obj("presolved")),
                    highs.get_presolved_lp(),
                );
            }
            s
        }
        "postsolve" => highs.postsolve_solution(&solution_from(rc)),
        "crossover" => highs.crossover(&solution_from(rc)),
        "get_standard_form_lp" => {
            let (mut nc, mut nr, mut nz, mut off) = (0, 0, 0, 0.0);
            let s0 = highs.get_standard_form_lp(
                &mut nc, &mut nr, &mut nz, &mut off, None, None, None, None, None,
            );
            chk.int("num_col", rc.int("num_col"), nc)
                .int("num_row", rc.int("num_row"), nr)
                .int("num_nz", rc.int("num_nz"), nz);
            if s0 == HighsStatus::kError {
                s0
            } else {
                let mut cost = vec![0.0; nc as usize];
                let mut rhs = vec![0.0; nr as usize];
                let mut start = vec![0; nc as usize + 1];
                let mut index = vec![0; nz as usize];
                let mut value = vec![0.0; nz as usize];
                let s = highs.get_standard_form_lp(
                    &mut nc,
                    &mut nr,
                    &mut nz,
                    &mut off,
                    Some(&mut cost),
                    Some(&mut rhs),
                    Some(&mut start),
                    Some(&mut index),
                    Some(&mut value),
                );
                chk.f64("offset", rc.f64("offset"), off)
                    .f64s("cost", &rc.f64s("cost"), &cost)
                    .f64s("rhs", &rc.f64s("rhs"), &rhs)
                    .ints("start", &rc.ints("start"), &start)
                    .ints("index", &rc.ints("index"), &index)
                    .f64s("value", &rc.f64s("value"), &value);
                s
            }
        }
        "get_primal_ray" | "get_dual_ray" | "get_dual_unboundedness_direction" => {
            let n = rc.int("dim") as usize;
            let mut has = false;
            let mut v = vec![0.0; n];
            let s = match op {
                "get_primal_ray" => highs.get_primal_ray(&mut has, Some(&mut v)),
                "get_dual_ray" => highs.get_dual_ray(&mut has, Some(&mut v)),
                _ => highs.get_dual_unboundedness_direction(&mut has, Some(&mut v)),
            };
            chk.bool("has", rc.bool("has"), has);
            if has {
                chk.f64s("value", &rc.f64s("value"), &v);
            }
            s
        }
        "get_fixed_lp" => {
            let mut lp = HighsLp::default();
            let s = highs.get_fixed_lp(&mut lp);
            if s != HighsStatus::kError {
                check_lp(chk, "fixed", &g1_lp_from(&rc.obj("fixed")), &lp);
            }
            s
        }
        "get_ill_conditioning" => {
            let mut ic = HighsIllConditioning::default();
            let s = highs.get_ill_conditioning(
                &mut ic,
                rc.bool("constraint"),
                i("method"),
                f("ill_conditioning_bound"),
            );
            let index: Vec<i32> = ic.record.iter().map(|r| r.index).collect();
            let mult: Vec<f64> = ic.record.iter().map(|r| r.multiplier).collect();
            chk.ints("index", &rc.ints("index"), &index).f64s(
                "multiplier",
                &rc.f64s("multiplier"),
                &mult,
            );
            s
        }
        "get_saved_mip_solutions" => {
            let saved = highs.get_saved_mip_solutions();
            let obj: Vec<f64> = saved.iter().map(|s| s.objective).collect();
            chk.f64s("objectives", &rc.f64s("objectives"), &obj);
            HighsStatus::kOk
        }
        "change_coeff" => highs.change_coeff(i("row"), i("col"), f("value")),
        "scale_col" => highs.scale_col(i("col"), f("value")),
        "scale_row" => highs.scale_row(i("row"), f("value")),
        "add_col" => {
            let index = rc.ints("index");
            highs.add_col(
                f("cost"),
                f("lower"),
                f("upper"),
                index.len() as HighsInt,
                &index,
                &rc.f64s("value"),
            )
        }
        "add_row" => {
            let index = rc.ints("index");
            highs.add_row(
                f("lower"),
                f("upper"),
                index.len() as HighsInt,
                &index,
                &rc.f64s("value"),
            )
        }
        other => panic!("unknown gap-fill op {other}"),
    };
    chk.int(&format!("{op}.status"), rc.int("status"), status as i32);
}

/// Replays one gap-fill fixture: `passModel(lp)` with the header overrides, then
/// every `g1_call` in order, checking each `g1_result` against the port
/// state at that point.
pub fn replay_g1(fx: &Fixture) -> Option<Checker> {
    let header = g1_header(fx);
    let mut chk = Checker::new(fx.file.clone());
    let mut highs = Highs::with_options(g1_options(header.str("overrides")));
    chk.int(
        "pass_model",
        HighsStatus::kOk as i32,
        highs.pass_model(g1_lp(fx)) as i32,
    );
    for (k, r) in fx.recs().enumerate().skip(1) {
        match r.t() {
            "g1_call" => apply_g1_call(&mut chk, &mut highs, &r),
            "g1_result" => check_g1_result(&mut chk, &r, &highs),
            "g1_file" => {
                // The text the oracle's `run` wrote to mip_improving_solution_file
                let got = std::fs::read_to_string(&highs.options_.mip_improving_solution_file)
                    .unwrap_or_default();
                chk.str("improving_solution_file", r.str("text"), &got);
            }
            other => panic!("{}: record {k}: unexpected type {other}", fx.file),
        }
    }
    Some(chk)
}

/// The `t_g1_*` body: replays every fixture of `g1/<dir>/`.
pub fn run_g1(dir: &str, test: &str) {
    let fxs = g1_fixtures(dir);
    if fxs.is_empty() {
        return;
    }
    let mut fails = Failures::default();
    for fx in &fxs {
        fails.run(&fx.file, || replay_g1(fx));
    }
    fails.finish(test);
}
