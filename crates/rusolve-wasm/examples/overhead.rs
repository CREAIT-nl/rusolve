//! Wasm hook overhead: time per call and end-to-end on the knapsack MIP of the tests.
//! Run: `cargo run --release -p rusolve-wasm --example overhead`. Numbers go to
//! docs/design/wasm-hooks.md.

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::hooks::*;
use rusolve_ekk::prelude::*;
use rusolve_wasm::{WasmHooks, WasmLimits};
use std::time::Instant;

const WAT: &str = include_str!("../components/most_fractional.wat");
const N: usize = 28;
const M: usize = 4;

fn lcg(s: &mut u64) -> u64 {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *s >> 33
}

#[allow(clippy::field_reassign_with_default)]
fn model() -> HighsLp {
    let mut s = 12345u64;
    let values: Vec<f64> = (0..N).map(|_| (20 + lcg(&mut s) % 60) as f64).collect();
    let mut w = vec![vec![0.0; N]; M];
    for r in w.iter_mut() {
        for x in r.iter_mut() {
            *x = (10 + lcg(&mut s) % 50) as f64;
        }
    }
    let mut lp = HighsLp::default();
    lp.num_col_ = N as HighsInt;
    lp.num_row_ = M as HighsInt;
    lp.sense_ = ObjSense::kMaximize;
    lp.col_cost_ = values;
    lp.col_lower_ = vec![0.0; N];
    lp.col_upper_ = vec![1.0; N];
    lp.row_lower_ = vec![-kHighsInf; M];
    lp.row_upper_ = w
        .iter()
        .map(|r| (r.iter().sum::<f64>() * 0.35).floor())
        .collect();
    lp.integrality_ = vec![HighsVarType::kInteger; N];
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = N as HighsInt;
    lp.a_matrix_.num_row_ = M as HighsInt;
    let mut start = vec![0 as HighsInt];
    let (mut index, mut value) = (Vec::new(), Vec::new());
    for j in 0..N {
        for (i, r) in w.iter().enumerate() {
            index.push(i as HighsInt);
            value.push(r[j]);
        }
        start.push(index.len() as HighsInt);
    }
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    lp
}

fn run(hooks: Option<Box<dyn MipHooks>>) -> (f64, i64, f64) {
    let mut h = Highs::default();
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_str("presolve", "off");
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
        h.get_info().mip_node_count,
        t.elapsed().as_secs_f64(),
    )
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

/// The same two rules as the WAT example, native.
struct Native;
impl MipHooks for Native {
    fn branching_score(&mut self, i: &BranchingInfo<'_>) -> Option<Vec<f64>> {
        Some(
            i.candidates
                .iter()
                .map(|cd| 0.5 - ((cd.value - cd.value.floor()) - 0.5).abs())
                .collect(),
        )
    }
    fn cut_score(&mut self, i: &CutScoreInfo) -> Option<f64> {
        Some(i.violation)
    }
}

fn main() {
    // per-call cost through the trait, outside the solver
    let cands: Vec<BranchingCandidate> = (0..10)
        .map(|j| BranchingCandidate {
            col: j,
            value: 0.37 + j as f64,
        })
        .collect();
    let zeros = vec![0.0; 10];
    let ones = vec![1.0; 10];
    let sol = vec![0.5; 10];
    let binfo = BranchingInfo {
        depth: 3,
        node_lower_bound: 0.0,
        upper_limit: 0.0,
        candidates: &cands,
        solution: &sol,
        col_lower: &zeros,
        col_upper: &ones,
    };
    let cinfo = CutScoreInfo {
        pool_row: 0,
        violation: 0.3,
        num_nonzeros: 5,
        num_active_nonzeros: 5,
        active_norm: 2.0,
        default_score: 0.1,
    };
    let reps = 200_000u32;
    let mut w = WasmHooks::new("w", WAT.as_bytes(), WasmLimits::default()).unwrap();
    let mut n = Native;
    let time = |f: &mut dyn FnMut()| {
        let t = Instant::now();
        for _ in 0..reps {
            f();
        }
        t.elapsed().as_secs_f64() / reps as f64 * 1e9
    };
    let wc = time(&mut || {
        std::hint::black_box(w.cut_score(&cinfo));
    });
    let wb = time(&mut || {
        std::hint::black_box(w.branching_score(&binfo));
    });
    let nc = time(&mut || {
        std::hint::black_box(n.cut_score(&cinfo));
    });
    let nb = time(&mut || {
        std::hint::black_box(n.branching_score(&binfo));
    });
    println!(
        "per call (ns): cut_score wasm {wc:.0} native {nc:.0}; branch_score(10 cand) wasm {wb:.0} native {nb:.0}"
    );
    let st = w.stats().lock().unwrap().clone();
    println!(
        "fuel: {} per call on average over {} calls",
        st.fuel_used / st.calls,
        st.calls
    );

    // end to end, median of 9
    let mut t = [vec![], vec![], vec![]];
    let mut info = [(0.0, 0i64); 3];
    for _ in 0..9 {
        let (o, nd, s) = run(None);
        t[0].push(s);
        info[0] = (o, nd);
        let (o, nd, s) = run(Some(Box::new(Native)));
        t[1].push(s);
        info[1] = (o, nd);
        let h = WasmHooks::new("w", WAT.as_bytes(), WasmLimits::default()).unwrap();
        let (o, nd, s) = run(Some(Box::new(h)));
        t[2].push(s);
        info[2] = (o, nd);
    }
    for (name, i) in [("no hooks", 0), ("native hooks", 1), ("wasm hooks", 2)] {
        println!(
            "{name:13} obj {} nodes {} median {:.4} s",
            info[i].0,
            info[i].1,
            median(t[i].clone())
        );
    }
}
