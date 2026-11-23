//! The Wasm component written in Rust (`components/rust`, built to `components/hooks_rust.wasm`
//! by `tools/build_wasm_components.sh` and committed, so these tests need no wasm32 target) behaves
//! like the hand-written WAT examples: same optimum, same tree, same cut, falls back when out of
//! fuel. `cargo test --release -p rusolve-wasm --test wasm_rust_component -- --ignored --nocapture`
//! prints the per-call overhead of the Rust build against the WAT.

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::hooks::*;
use rusolve_ekk::prelude::*;
use rusolve_wasm::{WasmHooks, WasmLimits, WasmStats};
use std::time::Instant;

const RUST: &[u8] = include_bytes!("../components/hooks_rust.wasm");
const MOST_FRACTIONAL: &str = include_str!("../components/most_fractional.wat");
const DEEPEST_FIRST: &str = include_str!("../components/deepest_first.wat");
const COVER_CUT: &str = include_str!("../components/cover_cut.wat");

/// The body of a one-memory example module (everything after its `memory` line, without the
/// closing parenthesis), so that several examples merge into the one module the Rust build is.
fn body(wat: &str) -> String {
    let at = wat.find("(memory (export \"memory\") 1)").unwrap() + 31;
    let b = &wat[at..];
    b[..b.rfind(')').unwrap()].to_string()
}

/// `branch_score`, `cut_score` and `select_node` of the two WAT examples in one module.
fn wat_score_hooks() -> String {
    format!(
        "(module (memory (export \"memory\") 1) {} {})",
        body(MOST_FRACTIONAL),
        body(DEEPEST_FIRST)
    )
}

const N: usize = 28;
const M: usize = 4;

fn lcg(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

#[allow(clippy::field_reassign_with_default)]
fn model() -> HighsLp {
    let mut s = 12345u64;
    let values: Vec<f64> = (0..N).map(|_| (20 + lcg(&mut s) % 60) as f64).collect();
    let mut weights = vec![vec![0.0; N]; M];
    for row in weights.iter_mut() {
        for w in row.iter_mut() {
            *w = (10 + lcg(&mut s) % 50) as f64;
        }
    }
    let caps: Vec<f64> = weights
        .iter()
        .map(|r| (r.iter().sum::<f64>() * 0.35).floor())
        .collect();
    let mut lp = HighsLp::default();
    lp.num_col_ = N as HighsInt;
    lp.num_row_ = M as HighsInt;
    lp.sense_ = ObjSense::kMaximize;
    lp.col_cost_ = values;
    lp.col_lower_ = vec![0.0; N];
    lp.col_upper_ = vec![1.0; N];
    lp.row_lower_ = vec![-kHighsInf; M];
    lp.row_upper_ = caps;
    lp.integrality_ = vec![HighsVarType::kInteger; N];
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = N as HighsInt;
    lp.a_matrix_.num_row_ = M as HighsInt;
    let mut start = vec![0 as HighsInt];
    let (mut index, mut value) = (Vec::new(), Vec::new());
    for j in 0..N {
        for (i, row) in weights.iter().enumerate() {
            index.push(i as HighsInt);
            value.push(row[j]);
        }
        start.push(index.len() as HighsInt);
    }
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    lp
}

fn new_highs(lp: HighsLp) -> Highs {
    let mut h = Highs::default();
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_str("presolve", "off");
    h.set_option_value_double("mip_heuristic_effort", 0.0);
    h.set_option_value_bool("mip_heuristic_run_feasibility_jump", false);
    h.set_option_value_bool("mip_detect_symmetry", false);
    h.set_option_value_int("threads", 1);
    h.pass_model(lp);
    h
}

/// (objective, nodes, component stats, hook stats) of a solve of `lp` with the optional component.
fn solve(
    lp: HighsLp,
    src: Option<(&[u8], WasmLimits)>,
) -> (f64, i64, Option<WasmStats>, HookStats) {
    let mut h = new_highs(lp);
    let stats = src.map(|(s, limits)| {
        let hooks = WasmHooks::new("wasm", s, limits).expect("component");
        let st = hooks.stats();
        h.set_mip_hooks(Box::new(hooks));
        st
    });
    assert_eq!(h.run(), HighsStatus::kOk);
    assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
    let info = h.get_info();
    let snap = stats.map(|s| s.lock().unwrap().clone());
    (
        info.objective_function_value,
        info.mip_node_count,
        snap,
        h.mip_hook_stats(),
    )
}

/// max 10x0+11x1+12x2+3x3, 5x0+5x1+5x2+4x3 <= 14: x0+x1+x2 <= 2 is a cover cut; optimum 26.
#[allow(clippy::field_reassign_with_default)]
fn cover_model() -> HighsLp {
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

/// The Rust build and the equivalent WAT module take the same decisions: same optimum, same tree,
/// same number of calls and of fuel-independent outcomes, and the tree differs from the default
/// (the rules do take decisions).
#[test]
fn rust_component_matches_wat_and_the_default_on_the_knapsack() {
    let (obj, nodes, _, _) = solve(model(), None);
    let wat = wat_score_hooks();
    let (wobj, wnodes, wst, _) = solve(model(), Some((wat.as_bytes(), WasmLimits::default())));
    let (robj, rnodes, rst, rhs) = solve(model(), Some((RUST, WasmLimits::default())));
    let (wst, rst) = (wst.unwrap(), rst.unwrap());
    assert!((robj - obj).abs() < 1e-6, "{robj} vs default {obj}");
    assert!((robj - wobj).abs() < 1e-6, "{robj} vs wat {wobj}");
    assert_eq!(rnodes, wnodes, "the Rust build explores the WAT tree");
    assert_ne!(rnodes, nodes, "the rules change the tree: they are in use");
    assert!(rst.calls > 0 && rst.ok == rst.calls, "{rst:?}");
    assert_eq!(rst.traps + rst.out_of_fuel + rst.invalid, 0, "{rst:?}");
    // the Rust module also exports `separate` (the cover cut of the other model), which the WAT
    // score module does not: its extra calls find no violated cut on this knapsack
    assert!(rst.calls >= wst.calls, "{rst:?} {wst:?}");
    assert_eq!(rhs.user_cuts_offered, 0, "no cover cut here: {rhs:?}");
}

/// Out of fuel with a tiny budget: every call falls back, nothing crashes, the tree is the default.
#[test]
fn rust_component_out_of_fuel_falls_back() {
    let (obj, nodes, _, _) = solve(model(), None);
    let tiny = WasmLimits {
        fuel_per_call: 3,
        ..WasmLimits::default()
    };
    let (robj, rnodes, rst, _) = solve(model(), Some((RUST, tiny)));
    let rst = rst.unwrap();
    assert!(rst.calls > 0 && rst.out_of_fuel == rst.calls, "{rst:?}");
    assert_eq!(rst.ok, 0);
    assert!((robj - obj).abs() < 1e-6);
    assert_eq!(rnodes, nodes, "every call fell back: same tree as default");
}

/// The Rust `separate` returns the same cover cut as `cover_cut.wat`.
#[test]
fn rust_component_separates_the_cover_cut_like_the_wat() {
    let (_, nodes, _, _) = solve(cover_model(), None);
    let (wobj, wnodes, wst, whs) = solve(
        cover_model(),
        Some((COVER_CUT.as_bytes(), WasmLimits::default())),
    );
    let (robj, rnodes, rst, rhs) = solve(cover_model(), Some((RUST, WasmLimits::default())));
    let rst = rst.unwrap();
    assert!(
        (robj - 26.0).abs() < 1e-6 && (wobj - 26.0).abs() < 1e-6,
        "{robj} {wobj}"
    );
    assert!(rst.calls > 0 && rst.ok == rst.calls, "{rst:?}");
    assert!(rhs.user_cuts_offered >= 1, "{rhs:?}");
    assert_eq!(
        rhs.user_cuts_rejected_malformed + rhs.user_cuts_rejected_debug_solution,
        0
    );
    assert_eq!(
        rhs.user_cuts_offered, whs.user_cuts_offered,
        "same cuts as the WAT"
    );
    assert_eq!(rnodes, wnodes);
    let _ = (nodes, wst);
}

/// Fuel use is a property of the module, not of the machine or the run.
#[test]
fn rust_component_fuel_is_deterministic() {
    let a = solve(model(), Some((RUST, WasmLimits::default())))
        .2
        .unwrap();
    let b = solve(model(), Some((RUST, WasmLimits::default())))
        .2
        .unwrap();
    assert_eq!(a, b);
}

fn per_call_ns(reps: u32, f: &mut dyn FnMut()) -> f64 {
    let t = Instant::now();
    for _ in 0..reps {
        f();
    }
    t.elapsed().as_secs_f64() / reps as f64 * 1e9
}

/// Overhead of the Rust build against the WAT, per call through the `MipHooks` trait (release
/// build only means anything). Numbers are recorded in docs/design/wasm-hooks.md.
#[test]
#[ignore]
fn rust_component_overhead_against_wat() {
    let reps = 200_000u32;
    let lim = WasmLimits::default;
    let mut wat = WasmHooks::new("w", wat_score_hooks().as_bytes(), lim()).unwrap();
    let mut rs = WasmHooks::new("r", RUST, lim()).unwrap();
    let mut wsep = WasmHooks::new("ws", COVER_CUT.as_bytes(), lim()).unwrap();
    let cands: Vec<BranchingCandidate> = (0..10)
        .map(|j| BranchingCandidate {
            col: j,
            value: 0.37 + j as f64,
        })
        .collect();
    let (zeros, ones, sol) = (vec![0.0; 10], vec![1.0; 10], vec![0.5; 10]);
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
    let open: Vec<OpenNodeInfo> = (0..32)
        .map(|i| OpenNodeInfo {
            lower_bound: i as f64,
            estimate: 0.0,
            depth: (i % 7) as HighsInt,
            num_domain_changes: 3,
        })
        .collect();
    let lo = vec![0.0; 14];
    let up = vec![1.0; 14];
    let xs: Vec<f64> = (0..14).map(|j| if j < 3 { 0.95 } else { 0.1 }).collect();
    let sinfo = SeparationInfo {
        lp_objective: 0.0,
        upper_limit: 0.0,
        feastol: 1e-6,
        num_nodes: 0,
        solution: &xs,
        col_lower: &lo,
        col_upper: &up,
    };
    let mut rows: Vec<(&str, f64, f64, u64, u64)> = Vec::new();
    macro_rules! cmp {
        ($name:expr, $wat:ident, $call:expr) => {{
            let a = per_call_ns(reps, &mut || {
                std::hint::black_box($call(&mut $wat));
            });
            let b = per_call_ns(reps, &mut || {
                std::hint::black_box($call(&mut rs));
            });
            let fa = {
                let s = $wat.stats().lock().unwrap().clone();
                s.fuel_used / s.calls.max(1)
            };
            let fb = {
                let s = rs.stats().lock().unwrap().clone();
                s.fuel_used / s.calls.max(1)
            };
            rows.push(($name, a, b, fa, fb));
            for h in [&mut $wat, &mut rs] {
                *h.stats().lock().unwrap() = WasmStats::default();
            }
        }};
    }
    cmp!("cut_score", wat, |h: &mut WasmHooks| h.cut_score(&cinfo));
    cmp!("branch_score, 10 candidates", wat, |h: &mut WasmHooks| h
        .branching_score(&binfo));
    cmp!("select_node, 32 open nodes", wat, |h: &mut WasmHooks| h
        .select_node(&open));
    cmp!("separate, 14 columns, 1 cut", wsep, |h: &mut WasmHooks| h
        .separate(&sinfo));
    println!("per call through MipHooks, ns (fuel per call): wat | rust");
    for (n, a, b, fa, fb) in rows {
        println!("{n:32} {a:7.0} ({fa}) | {b:7.0} ({fb})");
    }
}
