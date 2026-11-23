//! Oracle glue for QP, ranging, IIS and iCrash: fixture discovery, the base
//! options, the `w1_header` builders and the common `w1_result` check.
//! Used by `tests/t_w1_{qp,ranging,iis,icrash}.rs`; record layout in
//! tools/highs-oracle/FORMAT.md ("W1 records").
//!
//! Fixture discovery differs from `fixtures_in` in one respect: until the
//! oracle regenerates MANIFEST.tsv, `w1/...` files present in a target
//! directory but not yet listed are read without the sha256 check (a note
//! says so). Once listed, `load` (with the check) is used.

use crate::common::fixture::{FIXTURE_DIR_ENV, REQUIRE_ENV, target_dirs};
use crate::common::records::LpRec;
use crate::common::*;
use crate::port::*;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use serde_json::Value;
use std::path::PathBuf;

/// Every `w1/<dir>/*.jsonl` fixture available. Empty (after a SKIP banner)
/// when none was generated, unless RUSOLVE_REQUIRE_TARGET_FIXTURES=1.
pub fn w1_fixtures(dir: &str) -> Vec<Fixture> {
    wave_fixtures(&format!("w1/{dir}/"), "QP/IIS/ranging oracle commands")
}

/// Every `<prefix>*.jsonl` fixture available (`prefix` like `w3/cupdlp/`);
/// `producer` names the oracle commands for the SKIP banner. Shared by the
/// QP and PDLP glue.
pub fn wave_fixtures(prefix: &str, producer: &str) -> Vec<Fixture> {
    let prefix = prefix.to_string();
    let listed = manifest().select(&prefix, ".jsonl");
    if !listed.is_empty() {
        let mut files: Vec<String> = listed.iter().map(|e| e.file.clone()).collect();
        files.sort();
        return files.iter().filter_map(|f| load(f)).collect();
    }
    for root in target_dirs() {
        let d = root.join(&prefix);
        let Ok(rd) = std::fs::read_dir(&d) else {
            continue;
        };
        let mut paths: Vec<PathBuf> = rd
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "jsonl"))
            .collect();
        if paths.is_empty() {
            continue;
        }
        paths.sort();
        eprintln!(
            "note: {prefix}*.jsonl not in MANIFEST.tsv yet; read from {} without the sha256 check",
            d.display()
        );
        return paths
            .into_iter()
            .map(|p| read_unlisted(&prefix, p))
            .collect();
    }
    let reason =
        format!("no {prefix}*.jsonl fixtures (not in MANIFEST.tsv, not in any target dir)");
    if std::env::var(REQUIRE_ENV).is_ok_and(|v| v == "1") {
        panic!("{reason} ({REQUIRE_ENV}=1)");
    }
    eprintln!(
        "\n######## SKIPPED, NOT PASSED ########\n\
         # test:    {}\n\
         # reason:  {reason}\n\
         # fix:     run tools/highs-oracle/gen_fixtures.sh (needs reference/HiGHS and the\n\
         #          {producer}), or set {FIXTURE_DIR_ENV}=<dir>.\n\
         #          Set {REQUIRE_ENV}=1 to turn this skip into a failure.\n\
         ####################################\n",
        std::thread::current().name().unwrap_or("<unnamed>")
    );
    Vec::new()
}

fn read_unlisted(prefix: &str, path: PathBuf) -> Fixture {
    let name = path.file_name().unwrap().to_string_lossy().to_string();
    let file = format!("{prefix}{name}");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    let records = text
        .lines()
        .enumerate()
        .map(|(n, l)| {
            serde_json::from_str::<Value>(l)
                .unwrap_or_else(|e| panic!("fixture {file} line {}: invalid JSON: {e}", n + 1))
        })
        .collect();
    Fixture {
        file,
        path,
        location: Location::Target,
        records,
    }
}

/// The oracle's base options (`oracleBaseOptions`: output off, serial,
/// presolve off, `random_seed=0`) plus the `overrides` of a record
/// (`k=v,k=v` or `default`).
pub fn w1_options(overrides: &str) -> HighsOptions {
    let mut o = HighsOptions {
        output_flag: false,
        threads: 1,
        parallel: "off".to_string(),
        presolve: "off".to_string(),
        random_seed: 0,
        ..crate::port::profile_options()
    };
    for kv in overrides.split(',') {
        if kv == "default" || kv.is_empty() {
            continue;
        }
        let (k, v) = kv.split_once('=').expect("key=value override");
        let b = || v == "true" || v == "on" || v == "1";
        match k {
            "solver" => o.solver = v.to_string(),
            "simplex_strategy" => o.simplex_strategy = v.parse().unwrap(),
            "iis_strategy" => o.iis_strategy = v.parse().unwrap(),
            "qp_iteration_limit" => o.qp_iteration_limit = v.parse().unwrap(),
            "qp_nullspace_limit" => o.qp_nullspace_limit = v.parse().unwrap(),
            "icrash" => o.icrash = b(),
            "icrash_strategy" => o.icrash_strategy = v.to_string(),
            "icrash_iterations" => o.icrash_iterations = v.parse().unwrap(),
            other => panic!("unknown QP/IIS/ranging override {other}"),
        }
    }
    o.sync_log_options();
    o
}

/// The `w1_header` record (first line of every QP/IIS/ranging fixture).
pub fn w1_header(fx: &Fixture) -> Rec<'_> {
    let h = fx.rec(0);
    assert_eq!(h.t(), "w1_header", "{}: first record", fx.file);
    h
}

/// The LP the oracle passed to `Highs::passModel` (as read from the MPS).
pub fn header_lp(fx: &Fixture) -> HighsLp {
    lp_from_rec(&LpRec::from_rec(&w1_header(fx).obj("lp")))
}

/// The Hessian the oracle passed (`hessian` of `w1_header`, w1/qp only).
pub fn header_hessian(fx: &Fixture) -> HighsHessian {
    let h = w1_header(fx).obj("hessian");
    HighsHessian {
        dim_: h.int("dim"),
        format_: match h.int("format") {
            1 => HessianFormat::kTriangular,
            2 => HessianFormat::kSquare,
            f => panic!("{}: unknown HessianFormat {f}", h.ctx),
        },
        start_: h.ints("start"),
        index_: h.ints("index"),
        value_: h.f64s("value"),
    }
}

/// The hook lines (`t` in `types`) of `fx` in file order, verbatim.
pub fn hook_lines<'a>(fx: &Fixture, text: &'a str, types: &[&str]) -> Vec<&'a str> {
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), fx.len(), "{}", fx.file);
    (0..fx.len())
        .filter(|&i| types.contains(&fx.rec(i).t()))
        .map(|i| lines[i])
        .collect()
}

/// Byte-for-byte hook-line comparison; reports the first differing line
/// field by field.
pub fn diff_trace(chk: &mut Checker, want: &[&str], got: &str) {
    let got: Vec<&str> = got.lines().collect();
    for (k, (w, g)) in want.iter().zip(got.iter()).enumerate() {
        if w != g {
            let wv: Value = serde_json::from_str(w).unwrap();
            let gv: Value = serde_json::from_str(g).unwrap();
            for m in diff_values(&format!("trace[{k}]"), &wv, &gv) {
                chk.push(&m.path, m.expected, m.actual, m.note);
            }
            if chk.ok() {
                chk.push(
                    &format!("trace[{k}]"),
                    w.to_string(),
                    g.to_string(),
                    String::new(),
                );
            }
            return;
        }
    }
    chk.int("trace.len", want.len() as i32, got.len() as i32);
}

/// The common `w1_result` fields (FORMAT.md, "W1 records") against the port
/// after `run`.
pub fn check_w1_result(chk: &mut Checker, r: &Rec, run_status: HighsStatus, highs: &Highs) {
    let info = highs.get_info();
    let sol = highs.get_solution();
    let basis = highs.get_basis();
    chk.int("run_status", r.int("run_status"), run_status as i32)
        .int(
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
            "simplex_iteration_count",
            r.int("simplex_iteration_count"),
            info.simplex_iteration_count,
        )
        .int(
            "qp_iteration_count",
            r.int("qp_iteration_count"),
            info.qp_iteration_count,
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
            "num_primal_infeasibilities",
            r.int("num_primal_infeasibilities"),
            info.num_primal_infeasibilities,
        )
        .f64(
            "max_primal_infeasibility",
            r.f64("max_primal_infeasibility"),
            info.max_primal_infeasibility,
        )
        .f64(
            "sum_primal_infeasibilities",
            r.f64("sum_primal_infeasibilities"),
            info.sum_primal_infeasibilities,
        )
        .int(
            "num_dual_infeasibilities",
            r.int("num_dual_infeasibilities"),
            info.num_dual_infeasibilities,
        )
        .f64(
            "max_dual_infeasibility",
            r.f64("max_dual_infeasibility"),
            info.max_dual_infeasibility,
        )
        .f64(
            "sum_dual_infeasibilities",
            r.f64("sum_dual_infeasibilities"),
            info.sum_dual_infeasibilities,
        )
        .bool("basis_valid", r.bool("basis_valid"), basis.valid);
    let col_status: Vec<i32> = basis.col_status.iter().map(|s| *s as i32).collect();
    let row_status: Vec<i32> = basis.row_status.iter().map(|s| *s as i32).collect();
    chk.ints("col_status", &r.ints("col_status"), &col_status)
        .ints("row_status", &r.ints("row_status"), &row_status)
        .f64s("col_value", &r.f64s("col_value"), &sol.col_value)
        .f64s("col_dual", &r.f64s("col_dual"), &sol.col_dual)
        .f64s("row_value", &r.f64s("row_value"), &sol.row_value)
        .f64s("row_dual", &r.f64s("row_dual"), &sol.row_dual);
}
