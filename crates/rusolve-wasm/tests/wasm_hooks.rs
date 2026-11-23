//! Fuel-metered Wasm components plugged into `MipHooks`. A working component, one that traps
//! and one that never terminates must all leave the solve finishing with the right optimum.
//! Model: the multi-dimensional knapsack of the ekk `hooks_mip` tests (needs a real tree).

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use rusolve_wasm::{WasmHooks, WasmLimits, WasmStats};

const MOST_FRACTIONAL: &str = include_str!("../components/most_fractional.wat");
const DEEPEST_FIRST: &str = include_str!("../components/deepest_first.wat");

const TRAPS: &str = r#"(module (memory (export "memory") 1)
  (func (export "branch_score") (param i32 i32 i32 i32) (result i32) unreachable)
  (func (export "cut_score") (param f64 i32 i32 f64 f64) (result f64) unreachable))"#;

const LOOPS: &str = r#"(module (memory (export "memory") 1)
  (func (export "branch_score") (param i32 i32 i32 i32) (result i32) (loop $l (br $l)) (i32.const 0))
  (func (export "cut_score") (param f64 i32 i32 f64 f64) (result f64) (loop $l (br $l)) (f64.const 0)))"#;

const NAN_SCORES: &str = r#"(module (memory (export "memory") 1)
  (func (export "branch_score") (param $d i32) (param $n i32) (param $in i32) (param $out i32) (result i32)
    (f64.store (local.get $out) (f64.const nan)) (i32.const 0))
  (func (export "cut_score") (param f64 i32 i32 f64 f64) (result f64) (f64.const inf)))"#;

/// select_node that always answers index 999: out of range, invalid on every call.
const NODE_OUT_OF_RANGE: &str = r#"(module (memory (export "memory") 1)
  (func (export "select_node") (param i32 i32) (result i32) (i32.const 999)))"#;

/// select_node that never terminates: out of fuel on every call.
const NODE_LOOPS: &str = r#"(module (memory (export "memory") 1)
  (func (export "select_node") (param i32 i32) (result i32) (loop $l (br $l)) (i32.const -1)))"#;

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

fn new_highs() -> Highs {
    let mut h = Highs::default();
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_str("presolve", "off");
    h.set_option_value_double("mip_heuristic_effort", 0.0);
    h.set_option_value_bool("mip_heuristic_run_feasibility_jump", false);
    h.set_option_value_int("threads", 1);
    h.pass_model(model());
    h
}

/// (objective, nodes, component stats) of a solve with the optional component installed.
fn solve(src: Option<(&str, WasmLimits)>) -> (f64, i64, Option<WasmStats>) {
    let mut h = new_highs();
    let stats = src.map(|(s, limits)| {
        let hooks = WasmHooks::new("wasm", s.as_bytes(), limits).expect("component");
        let st = hooks.stats();
        h.set_mip_hooks(Box::new(hooks));
        st
    });
    assert_eq!(h.run(), HighsStatus::kOk);
    assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
    let info = h.get_info();
    let snap = stats.map(|s| s.lock().unwrap().clone());
    (info.objective_function_value, info.mip_node_count, snap)
}

fn small() -> WasmLimits {
    WasmLimits {
        fuel_per_call: 20_000,
        ..WasmLimits::default()
    }
}

#[test]
fn wasm_working_component_runs_and_reaches_the_optimum() {
    let (obj, _, _) = solve(None);
    let (wobj, _, st) = solve(Some((MOST_FRACTIONAL, WasmLimits::default())));
    let st = st.unwrap();
    assert!((wobj - obj).abs() < 1e-6, "{wobj} vs {obj}");
    assert!(st.calls > 0 && st.ok == st.calls, "{st:?}");
    assert_eq!(st.traps + st.out_of_fuel + st.invalid, 0);
    assert!(st.fuel_used > 0);
}

#[test]
fn wasm_trapping_component_falls_back_and_search_is_unchanged() {
    let (obj, nodes, _) = solve(None);
    let (wobj, wnodes, st) = solve(Some((TRAPS, small())));
    let st = st.unwrap();
    assert!(st.traps > 0 && st.traps == st.calls, "{st:?}");
    assert_eq!(st.ok, 0);
    assert!((wobj - obj).abs() < 1e-6);
    assert_eq!(wnodes, nodes, "every call fell back: same tree as default");
}

#[test]
fn wasm_infinite_loop_runs_out_of_fuel_and_search_is_unchanged() {
    let (obj, nodes, _) = solve(None);
    let (wobj, wnodes, st) = solve(Some((LOOPS, small())));
    let st = st.unwrap();
    assert!(st.out_of_fuel > 0 && st.out_of_fuel == st.calls, "{st:?}");
    assert!(st.fuel_used >= st.out_of_fuel * 19_000, "{st:?}");
    assert!((wobj - obj).abs() < 1e-6);
    assert_eq!(wnodes, nodes);
}

#[test]
fn wasm_invalid_answers_fall_back() {
    let (obj, nodes, _) = solve(None);
    let (wobj, wnodes, st) = solve(Some((NAN_SCORES, small())));
    let st = st.unwrap();
    assert!(st.invalid > 0 && st.invalid == st.calls, "{st:?}");
    assert!((wobj - obj).abs() < 1e-6);
    assert_eq!(wnodes, nodes);
}

/// Node selection: the deepest-first component explores a different tree (it dives) and
/// still reaches the same optimum.
#[test]
fn wasm_select_node_deepest_first_reaches_the_optimum() {
    let (obj, _, _) = solve(None);
    let (wobj, _, st) = solve(Some((DEEPEST_FIRST, WasmLimits::default())));
    let st = st.unwrap();
    assert!((wobj - obj).abs() < 1e-6, "{wobj} vs {obj}");
    assert!(st.calls > 0 && st.ok == st.calls, "{st:?}");
    assert_eq!(st.traps + st.out_of_fuel + st.invalid, 0);
    assert!(st.fuel_used > 0);
}

/// An out-of-range pick is invalid and falls back per call: the tree is the default one.
#[test]
fn wasm_select_node_out_of_range_falls_back_and_search_is_unchanged() {
    let (obj, nodes, _) = solve(None);
    let (wobj, wnodes, st) = solve(Some((NODE_OUT_OF_RANGE, small())));
    let st = st.unwrap();
    assert!(st.invalid > 0 && st.invalid == st.calls, "{st:?}");
    assert!((wobj - obj).abs() < 1e-6);
    assert_eq!(wnodes, nodes, "every call fell back: same tree as default");
}

/// A select_node that never terminates runs out of fuel per call: the tree is the default one.
#[test]
fn wasm_select_node_infinite_loop_runs_out_of_fuel_and_search_is_unchanged() {
    let (obj, nodes, _) = solve(None);
    let (wobj, wnodes, st) = solve(Some((NODE_LOOPS, small())));
    let st = st.unwrap();
    assert!(st.out_of_fuel > 0 && st.out_of_fuel == st.calls, "{st:?}");
    assert!(st.fuel_used >= st.out_of_fuel * 19_000, "{st:?}");
    assert!((wobj - obj).abs() < 1e-6);
    assert_eq!(wnodes, nodes);
}

#[test]
fn wasm_fuel_use_is_deterministic() {
    let a = solve(Some((MOST_FRACTIONAL, WasmLimits::default())))
        .2
        .unwrap();
    let b = solve(Some((MOST_FRACTIONAL, WasmLimits::default())))
        .2
        .unwrap();
    assert_eq!(a, b);
}

#[test]
fn wasm_modules_with_imports_or_bad_signatures_are_rejected() {
    let imp = r#"(module (import "env" "now" (func (result i64))))"#;
    let e = WasmHooks::new("x", imp.as_bytes(), WasmLimits::default())
        .err()
        .unwrap();
    assert!(e.contains("no host access"), "{e}");
    let bad = r#"(module (memory (export "memory") 1) (func (export "cut_score") (param i32)))"#;
    assert!(WasmHooks::new("x", bad.as_bytes(), WasmLimits::default()).is_err());
    let bad_select =
        r#"(module (memory (export "memory") 1) (func (export "select_node") (param i32)))"#;
    assert!(WasmHooks::new("x", bad_select.as_bytes(), WasmLimits::default()).is_err());
    let nomem = r#"(module (func (export "branch_score") (param i32 i32 i32 i32) (result i32) (i32.const 1)))"#;
    assert!(WasmHooks::new("x", nomem.as_bytes(), WasmLimits::default()).is_err());
    let nomem_select =
        r#"(module (func (export "select_node") (param i32 i32) (result i32) (i32.const -1)))"#;
    assert!(WasmHooks::new("x", nomem_select.as_bytes(), WasmLimits::default()).is_err());
}

#[test]
fn wasm_memory_cap_is_enforced() {
    let big = r#"(module (memory (export "memory") 64))"#; // 4 MiB declared, cap 1 MiB
    assert!(WasmHooks::new("x", big.as_bytes(), WasmLimits::default()).is_err());
}
