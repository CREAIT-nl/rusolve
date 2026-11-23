//! End-to-end checks of the `rusolve mip`
//! CSV line for two reader/reporting defects found by a gap audit —
//! the MPS objective constant (`Model::obj_offset`) and `LI`/`UI` bounds.

use std::path::PathBuf;
use std::process::Command;

fn write_tmp(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusolve-cli-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

/// `(status, objective, bound)` of `rusolve mip <path> --csv <extra>`.
fn mip_csv(path: &PathBuf, extra: &[&str]) -> (String, f64, f64) {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "0") // native search (`--cutoff` is native-only)
        .arg("mip")
        .arg(path)
        .arg("--csv")
        .args(extra)
        .output()
        .unwrap();
    let stdout = String::from_utf8(out.stdout).unwrap();
    let line = stdout.lines().nth(1).unwrap_or_else(|| panic!("{stdout}"));
    let f: Vec<&str> = line.split(',').collect();
    (
        f[0].to_string(),
        f[1].parse().unwrap(),
        f[2].parse().unwrap(),
    )
}

/// max x + 2y s.t. x + y ≤ 3.5, x, y integer in [0, 2]: optimum (1, 2),
/// c·x = 5. The N-row RHS entry −10 is the constant +10 (HiGHS
/// `obj_offset = -val`), so the reported objective is 15.
const OFFSET_MIP: &str = "\
NAME OFFMIP
OBJSENSE
  MAX
ROWS
 N  obj
 L  r0
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  r0  1
    y  obj  2  r0  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  obj  -10  r0  3.5
BOUNDS
 UP bnd x 2
 UP bnd y 2
ENDATA
";

#[test]
fn mip_reports_the_objective_constant() {
    let p = write_tmp("offset.mps", OFFSET_MIP);
    let (status, obj, bound) = mip_csv(&p, &[]);
    assert_eq!(status, "Optimal");
    assert_eq!(obj, 15.0, "c·x = 5 plus the constant 10");
    assert_eq!(bound, 15.0);
    // `--cutoff` lives in the same reported space: 15 is achievable, so a
    // cutoff at it keeps the optimum (a raw-space cutoff of 15 would
    // declare the problem infeasible, max c·x = 5 < 15).
    let (status, obj, _) = mip_csv(&p, &["--cutoff", "14.5"]);
    assert_eq!((status.as_str(), obj), ("Optimal", 15.0));
}

/// min −x − y − z with general integers bounded only by LI/UI lines (HiGHS
/// writer shape): x ∈ [2, 7] (LI then UI), y ≤ 4 (UI first in BOUNDS),
/// z ∈ [−3, ∞) capped by a row z ≤ 1. Optimum (7, 4, 1) → −12; the old
/// reader errored on the leading UI line.
#[test]
fn mip_honours_li_and_ui_bounds() {
    let p = write_tmp(
        "liui.mps",
        "\
NAME LIUI
ROWS
 N  obj
 L  r0
 L  r1
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  -1  r0  1
    y  obj  -1  r0  1
    z  obj  -1  r1  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  r0  100  r1  1
BOUNDS
 UI bnd y 4
 LI bnd x 2
 UI bnd x 7
 LI bnd z -3
ENDATA
",
    );
    let (status, obj, _) = mip_csv(&p, &[]);
    assert_eq!((status.as_str(), obj), ("Optimal", -12.0));
    // x ≥ 2 is binding when the objective pushes x down.
    let p = write_tmp(
        "li_low.mps",
        "\
NAME LILOW
ROWS
 N  obj
 L  r0
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  r0  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  r0  10
BOUNDS
 LI bnd x 2
ENDATA
",
    );
    let (status, obj, _) = mip_csv(&p, &[]);
    assert_eq!((status.as_str(), obj), ("Optimal", 2.0));
}

/// max 8a + 11b + 6c + 4d s.t. 5a + 7b + 4c + 3d <= 14, binaries: optimum b, c, d = 21.
const KNAPSACK_MIP: &str = "\
NAME KNAP
OBJSENSE
  MAX
ROWS
 N  obj
 L  cap
COLUMNS
    MARKER  'MARKER'  'INTORG'
    a  obj  8  cap  5
    b  obj  11  cap  7
    c  obj  6  cap  4
    d  obj  4  cap  3
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  cap  14
BOUNDS
 BV bnd a
 BV bnd b
 BV bnd c
 BV bnd d
ENDATA
";

/// The conflict ablation switch with every source off still reaches the optimum
/// the default search reaches (the switch only drops conflicts, which are redundant cuts).
#[test]
fn conflict_off_mask_63_reaches_the_same_optimum_as_mask_0() {
    let p = write_tmp("knap_conflict_mask.mps", KNAPSACK_MIP);
    for mask in ["0", "63"] {
        let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
            .env("RUSOLVE_MIP_HIGHS", "1")
            .arg("mip")
            .arg(&p)
            .arg("--csv")
            .args(["--opt", &format!("mip_conflict_off_mask={mask}")])
            .output()
            .unwrap();
        let stdout = String::from_utf8(out.stdout).unwrap();
        let line = stdout.lines().nth(1).unwrap_or_else(|| panic!("{stdout}"));
        let f: Vec<&str> = line.split(',').collect();
        assert!(f[0].contains("ptimal"), "mask {mask}: status {}", f[0]);
        let obj: f64 = f[1].parse().unwrap();
        assert!((obj - 21.0).abs() < 1e-6, "mask {mask}: objective {obj}");
    }
}
