//! Wasm user-separator overhead: time per `separate` call, native Rust against the
//! `cover_cut.wat` component, and end to end on the 4-item cover model of
//! `tests/wasm_separator.rs` (median of 201 solves). Run:
//! `cargo run --release -p rusolve-wasm --example separator_overhead`; numbers go to
//! docs/design/wasm-hooks.md.

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::hooks::*;
use rusolve_ekk::prelude::*;
use rusolve_wasm::{WasmHooks, WasmLimits};
use std::time::Instant;

const WAT: &str = include_str!("../components/cover_cut.wat");

struct Native;
impl MipHooks for Native {
    fn separate(&mut self, i: &SeparationInfo<'_>) -> Vec<UserCut> {
        if i.solution[0] + i.solution[1] + i.solution[2] > 2.0001 {
            vec![UserCut {
                indices: vec![0, 1, 2],
                values: vec![1.0; 3],
                rhs: 2.0,
            }]
        } else {
            Vec::new()
        }
    }
}

#[allow(clippy::field_reassign_with_default)]
fn model() -> HighsLp {
    let mut lp = HighsLp::default();
    lp.num_col_ = 4;
    lp.num_row_ = 1;
    lp.sense_ = ObjSense::kMaximize;
    lp.col_cost_ = vec![10.0, 11.0, 12.0, 3.0];
    lp.col_lower_ = vec![0.0; 4];
    lp.col_upper_ = vec![1.0; 4];
    lp.row_lower_ = vec![-kHighsInf];
    lp.row_upper_ = vec![14.0];
    lp.integrality_ = vec![HighsVarType::kInteger; 4];
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 4;
    lp.a_matrix_.num_row_ = 1;
    lp.a_matrix_.start_ = vec![0, 1, 2, 3, 4];
    lp.a_matrix_.index_ = vec![0; 4];
    lp.a_matrix_.value_ = vec![5.0, 5.0, 5.0, 4.0];
    lp
}

fn run(hooks: Option<Box<dyn MipHooks>>) -> (f64, f64) {
    let mut h = Highs::default();
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_str("presolve", "off");
    h.set_option_value_bool("mip_detect_symmetry", false);
    h.set_option_value_double("mip_heuristic_effort", 0.0);
    h.set_option_value_bool("mip_heuristic_run_feasibility_jump", false);
    h.set_option_value_int("threads", 1);
    h.pass_model(model());
    if let Some(k) = hooks {
        h.set_mip_hooks(k);
    }
    let t = Instant::now();
    h.run();
    (
        h.get_info().objective_function_value,
        t.elapsed().as_secs_f64(),
    )
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() {
    let lo = vec![0.0; 14];
    let up = vec![1.0; 14];
    let sol: Vec<f64> = (0..14).map(|j| if j < 3 { 0.95 } else { 0.1 }).collect();
    let info = SeparationInfo {
        lp_objective: 0.0,
        upper_limit: 0.0,
        feastol: 1e-6,
        num_nodes: 0,
        solution: &sol,
        col_lower: &lo,
        col_upper: &up,
    };
    let reps = 200_000u32;
    let mut w = WasmHooks::new("w", WAT.as_bytes(), WasmLimits::default()).unwrap();
    let mut n = Native;
    let t = Instant::now();
    for _ in 0..reps {
        std::hint::black_box(w.separate(&info));
    }
    let wns = t.elapsed().as_secs_f64() / reps as f64 * 1e9;
    let t = Instant::now();
    for _ in 0..reps {
        std::hint::black_box(n.separate(&info));
    }
    let nns = t.elapsed().as_secs_f64() / reps as f64 * 1e9;
    let st = w.stats().lock().unwrap().clone();
    println!(
        "separate per call (ns), 14 columns, 1 cut returned: wasm {wns:.0} native {nns:.0}; fuel {} per call",
        st.fuel_used / st.calls
    );
    let mut t = [vec![], vec![], vec![]];
    let mut obj = [0.0; 3];
    for _ in 0..201 {
        let (o, s) = run(None);
        t[0].push(s);
        obj[0] = o;
        let (o, s) = run(Some(Box::new(Native)));
        t[1].push(s);
        obj[1] = o;
        let h = WasmHooks::new("w", WAT.as_bytes(), WasmLimits::default()).unwrap();
        let (o, s) = run(Some(Box::new(h)));
        t[2].push(s);
        obj[2] = o;
    }
    for (name, i) in [
        ("no hooks", 0),
        ("native separator", 1),
        ("wasm separator", 2),
    ] {
        println!(
            "{name:17} obj {} median {:.1} us",
            obj[i],
            median(t[i].clone()) * 1e6
        );
    }
}
