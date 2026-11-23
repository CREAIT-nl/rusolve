//! `mip_shift_and_propagate` (Berthold and Hendel 2015; SCIP
//! `heur_shiftandpropagate.c`) finds a feasible point of a small
//! zero-objective bin-packing model before the root, with feasibility jump and the other
//! root heuristics off, and credits it to its own solution source.

use std::process::Command;

fn bin_packing_lp() -> String {
    let sizes = [4, 5, 6, 7, 8, 9, 3, 2];
    let (bins, cap) = (4, 15);
    let mut s = String::from("Minimize\n obj: 0 x_0_0\nSubject To\n");
    for i in 0..sizes.len() {
        let terms: Vec<String> = (0..bins).map(|b| format!("x_{i}_{b}")).collect();
        s += &format!(" a{i}: {} = 1\n", terms.join(" + "));
    }
    for b in 0..bins {
        let terms: Vec<String> = sizes
            .iter()
            .enumerate()
            .map(|(i, w)| format!("{w} x_{i}_{b}"))
            .collect();
        s += &format!(" c{b}: {} <= {cap}\n", terms.join(" + "));
    }
    s += "Binary\n";
    for i in 0..sizes.len() {
        for b in 0..bins {
            s += &format!(" x_{i}_{b}");
        }
    }
    s + "\nEnd\n"
}

fn solve(path: &std::path::Path, on: bool) -> serde_json::Value {
    let out = Command::new(env!("CARGO_BIN_EXE_rusolve-cli"))
        .env("RUSOLVE_MIP_HIGHS", "1")
        .arg("mip")
        .arg(path)
        .args(["--threads", "1", "--format", "json"])
        .args(["--opt", &format!("mip_shift_and_propagate={on}")])
        .args(["--opt", "mip_heuristic_run_feasibility_jump=false"])
        .args(["--opt", "mip_alns=false", "--opt", "presolve=off"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout).unwrap()
}

#[test]
fn shift_and_propagate_finds_a_bin_packing_point() {
    let dir = std::env::temp_dir().join(format!("rusolve-shiftprop-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("bp.lp");
    std::fs::write(&path, bin_packing_lp()).unwrap();

    let on = solve(&path, true);
    assert_eq!(on["status"], "Optimal", "{on}");
    let by_source = &on["stats"]["phases"]["heuristics"]["improving_incumbents_by_source"];
    assert_eq!(by_source[0]["source"], "Shift and propagate", "{on}");

    let off = solve(&path, false);
    assert_eq!(off["status"], "Optimal");
    let by_source =
        off["stats"]["phases"]["heuristics"]["improving_incumbents_by_source"].to_string();
    assert!(!by_source.contains("Shift and propagate"), "{by_source}");
    let _ = std::fs::remove_dir_all(&dir);
}
