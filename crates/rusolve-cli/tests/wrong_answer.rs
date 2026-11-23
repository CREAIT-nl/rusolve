//! Wrong-answer regressions: instances where a HiGHS release (and the port, faithfully) reported
//! a non-optimal objective as `Optimal`. Each case pins the true optimum (HiGHS 1.15.1 and
//! presolve-off runs agree on it) and runs the default options, which include the 1.16.0 levers.
//!
//! int_production s22-s24 (a wrong optimum in HiGHS 1.16.0). Implied-bound cuts on
//! implicit-integer columns were flagged integral with a fractional rhs, and the tableau CMIR
//! then cut off the optimum (`HighsImplications::implied_cut_is_integral`).

use std::path::Path;
use std::process::Command;

use serde_json::Value;

const CRATE_DIR: &str = env!("CARGO_MANIFEST_DIR");

/// (seed, true optimum, objective 1.16.0 reported before the fix)
const INT_PRODUCTION: [(u32, f64, f64); 3] = [
    (22, 85020.868, 85260.868),
    (23, 63288.834, 65208.323),
    (24, 67668.734, 75723.734),
];

fn solve(model: &Path, extra: &[&str]) -> Value {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(model)
        .args(["--opt", "mip_rel_gap=0", "--opt", "mip_abs_gap=1e-9"])
        .args(extra)
        .args(["--format", "json"])
        .output()
        .unwrap();
    assert_eq!(
        out.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout).unwrap()
}

fn check_int_production(extra: &[&str]) {
    for (seed, optimum, wrong) in INT_PRODUCTION {
        let model = Path::new(CRATE_DIR).join(format!(
            "../../benches/data/wrong_answer/int_production_s{seed}.mps"
        ));
        let r = solve(&model, extra);
        assert_eq!(r["status"], "Optimal", "s{seed} {extra:?}");
        let obj = r["objective"].as_f64().unwrap();
        assert!(
            (obj - optimum).abs() <= 1e-3,
            "s{seed} {extra:?}: objective {obj}, true optimum {optimum} (pre-fix value {wrong})"
        );
    }
}

#[test]
fn int_production_default_options_reach_the_true_optimum() {
    check_int_production(&[]);
}

#[test]
fn int_production_with_116_bound_rounding_reaches_the_true_optimum() {
    // #3320 rounds implied bounds on integral columns, which is what exposed the bad cut flag.
    check_int_production(&["--opt", "h16_int_bound_round=true"]);
}

#[test]
fn int_production_with_1151_behaviour_reaches_the_true_optimum() {
    check_int_production(&["--opt", "h16_int_bound_round=false"]);
}
