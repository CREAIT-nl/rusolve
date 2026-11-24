//! `mip_engine=cpsat` (docs/design/cpsat-engine.md). The OR-Tools library is loaded at run time:
//! the refusal and layout tests run everywhere, the solve tests skip (and say so) when the library
//! is not found. Point `RUSOLVE_ORTOOLS_LIB` at it to run them
//! (`RUSOLVE_ORTOOLS_LIB=~/opt/or-tools/lib/libortools.so.9 cargo test --release -p rusolve-cli --test cpsat`).

use serde_json::Value;
use std::path::PathBuf;
use std::process::{Command, Output};

fn write_tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cpsat-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn mip(path: &PathBuf, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("mip")
        .arg(path)
        .args(["--format", "json", "--time-limit", "60"])
        .args(extra)
        .output()
        .unwrap()
}

fn json(out: &Output) -> Value {
    let text = String::from_utf8_lossy(&out.stdout);
    serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("stdout is not one JSON object ({e}): {text}"))
}

/// `rusolve-cli version --format json` of `exe` with `env` set; `clean` first removes
/// `RUSOLVE_ORTOOLS_LIB` and `LD_LIBRARY_PATH`, so only the layout around `exe` can find the library.
fn version_of(exe: &std::path::Path, env: &[(&str, &str)], clean: bool) -> Value {
    let mut c = Command::new(exe);
    c.args(["version", "--format", "json"]);
    if clean {
        c.env_remove("RUSOLVE_ORTOOLS_LIB")
            .env_remove("LD_LIBRARY_PATH");
    }
    for (k, v) in env {
        c.env(k, v);
    }
    let out = c.output().unwrap();
    assert_eq!(out.status.code(), Some(0), "{out:?}");
    json(&out)
}

/// The library this machine has, if any, as the loader of the test binary finds it.
fn found_lib() -> Option<PathBuf> {
    let v = version_of(
        std::path::Path::new(env!("CARGO_BIN_EXE_rusolve-cli")),
        &[],
        false,
    );
    v["cpsat"]["library"].as_str().map(PathBuf::from)
}

/// Skips (returns from) a test that needs the OR-Tools library when there is none.
macro_rules! need_lib {
    () => {
        if found_lib().is_none() {
            eprintln!("skipped: the OR-Tools library was not found (set RUSOLVE_ORTOOLS_LIB)");
            return;
        }
    };
}

/// max 5x + 4y  s.t. 6x + 4y <= 24, x + 2y <= 6, x, y integer in [0, 10]: optimum 20.
const SMALL_MIP: &str = "\
NAME SMALL
OBJSENSE
  MAX
ROWS
 N  obj
 L  r0
 L  r1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  5  r0  6
    x  r1  1
    y  obj  4  r0  4
    y  r1  2
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  r0  24  r1  6
BOUNDS
 UP bnd x 10
 UP bnd y 10
ENDATA
";

/// The same model with a continuous column.
const CONTINUOUS_MIP: &str = "\
NAME CONT
ROWS
 N  obj
 G  r0
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  r0  1
    MARKER  'MARKER'  'INTEND'
    y  obj  1  r0  1
RHS
    rhs  r0  2.5
BOUNDS
 UP bnd x 10
 UP bnd y 10
ENDATA
";

#[test]
fn a_missing_library_is_a_typed_error_that_says_where_it_looked() {
    let p = write_tmp("nolib.mps", SMALL_MIP);
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .arg("mip")
        .arg(&p)
        .args(["--format", "json", "--opt", "mip_engine=cpsat"])
        .env("RUSOLVE_ORTOOLS_LIB", "/nonexistent/libortools.so.9")
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(3), "{out:?}");
    let text = json(&out).to_string();
    assert!(text.contains("RS0012"), "{text}");
    assert!(text.contains("/nonexistent/libortools.so.9"), "{text}");
    assert!(text.contains("install.sh --with-cpsat"), "{text}");
    assert!(text.contains("pip install"), "{text}");
}

#[test]
fn version_says_what_is_available_without_solving() {
    let exe = std::path::Path::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    let v = version_of(
        exe,
        &[("RUSOLVE_ORTOOLS_LIB", "/nonexistent/libortools.so.9")],
        false,
    );
    assert_eq!(v["name"], "rusolve");
    assert_eq!(v["cpsat"]["available"], false, "{v}");
    assert!(v["cpsat"]["library"].is_null());
    let searched = v["cpsat"]["searched"].to_string();
    assert!(searched.contains("/nonexistent/libortools.so.9"), "{v}");
    assert!(
        v["cpsat"]["error"]
            .as_str()
            .unwrap()
            .contains("--with-cpsat")
    );
}

#[test]
fn search_order_install_prefix_then_the_ortools_package() {
    need_lib!();
    let real = found_lib().unwrap();
    let name = rusolve_cpsat::runtime::LIB_NAME;
    let root = std::env::temp_dir().join(format!("rusolve-layout-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    let exe_src = std::path::Path::new(env!("CARGO_BIN_EXE_rusolve-cli"));
    // install.sh and the wheel both put the library next to the files it depends on (its run path
    // is $ORIGIN, taken from the path it was opened by), so link every file of the real directory
    let link = |to: &PathBuf| {
        let dir = to.parent().unwrap();
        std::fs::create_dir_all(dir).unwrap();
        for e in std::fs::read_dir(real.parent().unwrap()).unwrap() {
            let e = e.unwrap().path();
            let target = dir.join(e.file_name().unwrap());
            #[cfg(unix)]
            std::os::unix::fs::symlink(&e, &target).unwrap();
            #[cfg(not(unix))]
            std::fs::copy(&e, &target).unwrap();
        }
        // the wheel's vendored Abseil libraries: `$ORIGIN/../../ortools.libs` of the library's
        // directory (install.sh: <prefix>/lib/ortools.libs; wheel: site-packages/ortools.libs)
        let vendored = real
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("ortools.libs");
        if vendored.is_dir() {
            let dest = dir.parent().unwrap().parent().unwrap().join("ortools.libs");
            std::fs::create_dir_all(&dest).unwrap();
            for e in std::fs::read_dir(&vendored).unwrap() {
                let e = e.unwrap().path();
                #[cfg(unix)]
                std::os::unix::fs::symlink(&e, dest.join(e.file_name().unwrap())).unwrap();
                #[cfg(not(unix))]
                std::fs::copy(&e, dest.join(e.file_name().unwrap())).unwrap();
            }
        }
    };
    // (a) install.sh layout: <prefix>/bin/rusolve-cli, <prefix>/lib/rusolve/ortools/<lib>
    let prefix = root.join("prefix");
    std::fs::create_dir_all(prefix.join("bin")).unwrap();
    let exe_a = prefix.join("bin/rusolve-cli");
    std::fs::copy(exe_src, &exe_a).unwrap();
    let lib_a = prefix.join("lib/rusolve/ortools").join(name);
    link(&lib_a);
    std::fs::write(lib_a.with_file_name("VERSION"), "9.15.6755\n").unwrap();
    let v = version_of(&exe_a, &[], true);
    assert_eq!(v["cpsat"]["available"], true, "{v}");
    assert_eq!(v["cpsat"]["library"], lib_a.to_str().unwrap(), "{v}");
    assert_eq!(v["cpsat"]["version"], "9.15.6755", "{v}");
    // (b) wheel layout: site-packages/rusolve/_bin/rusolve-cli, site-packages/ortools/.libs/<lib>
    let site = root.join("site-packages");
    std::fs::create_dir_all(site.join("rusolve/_bin")).unwrap();
    let exe_b = site.join("rusolve/_bin/rusolve-cli");
    std::fs::copy(exe_src, &exe_b).unwrap();
    let lib_b = site.join("ortools/.libs").join(name);
    link(&lib_b);
    std::fs::create_dir_all(site.join("ortools-9.15.6755.dist-info")).unwrap();
    let v = version_of(&exe_b, &[], true);
    assert_eq!(v["cpsat"]["library"], lib_b.to_str().unwrap(), "{v}");
    assert_eq!(v["cpsat"]["version"], "9.15.6755", "{v}");
    // the environment variable wins over both
    let v = version_of(
        &exe_b,
        &[("RUSOLVE_ORTOOLS_LIB", real.to_str().unwrap())],
        true,
    );
    assert_eq!(v["cpsat"]["library"], real.to_str().unwrap(), "{v}");
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn an_unknown_engine_is_a_usage_error() {
    let p = write_tmp("badengine.mps", SMALL_MIP);
    let out = mip(&p, &["--opt", "mip_engine=nosuch"]);
    assert_eq!(out.status.code(), Some(2), "{out:?}");
    assert!(json(&out).to_string().contains("mip_engine=nosuch"));
}

#[test]
fn the_default_engine_is_untouched() {
    let p = write_tmp("default.mps", SMALL_MIP);
    let v = json(&mip(&p, &[]));
    assert_eq!(v["status"], "Optimal");
    assert_eq!(v["engine"], "highs-port");
    assert!(v.get("cpsat").is_none());
}

mod with_cpsat {
    use super::*;

    /// Two jobs on two machines, disjunctive big-M: J1 = M1(3) then M2(2), J2 = M2(4) then M1(1).
    const JOBSHOP: &str = "\
NAME JS22
ROWS
 N  cmax
 G  j1
 G  j2
 G  m1a
 G  m1b
 G  m2a
 G  m2b
 G  c1
 G  c2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    s11  j1  -1  m1a  -1
    s11  m1b  1
    s12  j1  1  m2a  -1
    s12  m2b  1  c1  -1
    s21  j2  -1  m2a  1
    s21  m2b  -1
    s22  j2  1  m1a  1
    s22  m1b  -1  c2  -1
    y1  m1a  -20  m1b  20
    y2  m2a  -20  m2b  20
    cmax  cmax  1  c1  1
    cmax  c2  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  j1  3  j2  4
    rhs  m1a  -17  m1b  1
    rhs  m2a  -18  m2b  4
    rhs  c1  -2  c2  -1
BOUNDS
 UP bnd s11 20
 UP bnd s12 20
 UP bnd s21 20
 UP bnd s22 20
 UP bnd cmax 40
 UP bnd y1 1
 UP bnd y2 1
ENDATA
";

    fn solved_alike(name: &str, text: &str, cp_extra: &[&str]) -> (Value, Value) {
        let p = write_tmp(name, text);
        let port = json(&mip(&p, &[]));
        let mut args = vec!["--opt", "mip_engine=cpsat"];
        args.extend_from_slice(cp_extra);
        let out = mip(&p, &args);
        assert_eq!(out.status.code(), Some(0), "{out:?}");
        (port, json(&out))
    }

    #[test]
    fn a_small_integer_program_matches_the_port() {
        need_lib!();
        let (port, cp) = solved_alike("small.mps", SMALL_MIP, &["--solution"]);
        assert_eq!(port["status"], "Optimal");
        assert_eq!(cp["status"], "Optimal");
        assert_eq!(cp["objective"], port["objective"]);
        assert_eq!(cp["objective"], 20.0);
        assert_eq!(cp["engine"], "cp-sat");
        assert_eq!(cp["cpsat"]["status"], "OPTIMAL");
        assert_eq!(cp["cpsat"]["checked"]["verdict"], "feasible");
        assert_eq!(cp["cpsat"]["objective_kind"], "integer");
        assert!(
            cp["cpsat"]["ortools_version"]
                .as_str()
                .unwrap()
                .starts_with("9.")
        );
        assert!(cp["cpsat"]["wall_time_s"].as_f64().unwrap() >= 0.0);
        assert_eq!(cp["cpsat"]["workers"], 1);
        assert_eq!(cp["solution"]["columns"].as_array().unwrap().len(), 2);
    }

    #[test]
    fn a_job_shop_matches_the_port() {
        need_lib!();
        let (port, cp) = solved_alike("js22.mps", JOBSHOP, &["--threads", "2"]);
        assert_eq!(port["status"], "Optimal");
        assert_eq!(cp["status"], "Optimal");
        assert_eq!(cp["objective"], port["objective"]);
        assert_eq!(cp["cpsat"]["workers"], 2);
        assert_eq!(cp["dual_bound"], cp["objective"]);
    }

    #[test]
    fn a_continuous_column_is_refused_by_name() {
        need_lib!();
        let p = write_tmp("continuous.mps", CONTINUOUS_MIP);
        let out = mip(&p, &["--opt", "mip_engine=cpsat"]);
        assert_eq!(out.status.code(), Some(3), "{out:?}");
        let text = json(&out).to_string();
        assert!(
            text.contains("RS0012") && text.contains("column y"),
            "{text}"
        );
        assert!(text.contains("continuous"), "{text}");
    }

    #[test]
    fn an_infeasible_model_is_proved() {
        need_lib!();
        let text = SMALL_MIP.replace("x  obj  5  r0  6", "x  obj  5  r0  6\n    x  r2  1");
        let text = text
            .replace(" L  r1\n", " L  r1\n G  r2\n")
            .replace("rhs  r0  24  r1  6", "rhs  r0  24  r1  6  r2  50");
        let p = write_tmp("infeasible.mps", &text);
        let out = mip(&p, &["--opt", "mip_engine=cpsat"]);
        assert_eq!(out.status.code(), Some(0), "{out:?}");
        let v = json(&out);
        assert_eq!(v["status"], "Infeasible", "{v}");
        assert_eq!(v["cpsat"]["status"], "INFEASIBLE");
        assert_eq!(v["has_solution"], false);
    }

    #[test]
    fn unsupported_flags_are_refused() {
        need_lib!();
        let p = write_tmp("flags.mps", SMALL_MIP);
        let out = mip(&p, &["--opt", "mip_engine=cpsat", "--watchdog"]);
        assert_eq!(out.status.code(), Some(2), "{out:?}");
    }
}
