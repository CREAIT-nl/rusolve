//! The user separator in Wasm (variable-size cut result). A working component adds a valid
//! cover cut; a trap, an infinite loop, an oversized or truncated answer and an invalid cut all
//! leave the solve at the right optimum (the host cut gate drops what is malformed).

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::hooks::*;
use rusolve_ekk::prelude::*;
use rusolve_wasm::{WasmHooks, WasmLimits, WasmStats};

const COVER_CUT: &str = include_str!("../components/cover_cut.wat");

const TRAPS: &str = r#"(module (memory (export "memory") 1)
  (func (export "separate") (param i32 i32 i32 i32) (result i32) unreachable))"#;
const LOOPS: &str = r#"(module (memory (export "memory") 1)
  (func (export "separate") (param i32 i32 i32 i32) (result i32) (loop $l (br $l)) (i32.const 0)))"#;
/// two records in one answer (variable size): x0+x1+x2 <= 2 and the weaker x0+x1+x2+x3 <= 3
const TWO_CUTS: &str = r#"(module (memory (export "memory") 1)
  (func (export "separate") (param i32 i32 i32 i32) (result i32)
    (f64.store (local.get 2) (f64.const 3))
    (f64.store offset=8 (local.get 2) (f64.const 2))
    (f64.store offset=16 (local.get 2) (f64.const 0))
    (f64.store offset=24 (local.get 2) (f64.const 1))
    (f64.store offset=32 (local.get 2) (f64.const 1))
    (f64.store offset=40 (local.get 2) (f64.const 1))
    (f64.store offset=48 (local.get 2) (f64.const 2))
    (f64.store offset=56 (local.get 2) (f64.const 1))
    (f64.store offset=64 (local.get 2) (f64.const 4))
    (f64.store offset=72 (local.get 2) (f64.const 3))
    (f64.store offset=80 (local.get 2) (f64.const 0))
    (f64.store offset=88 (local.get 2) (f64.const 1))
    (f64.store offset=96 (local.get 2) (f64.const 1))
    (f64.store offset=104 (local.get 2) (f64.const 1))
    (f64.store offset=112 (local.get 2) (f64.const 2))
    (f64.store offset=120 (local.get 2) (f64.const 1))
    (f64.store offset=128 (local.get 2) (f64.const 3))
    (f64.store offset=136 (local.get 2) (f64.const 1))
    (i32.const 18)))"#;

/// claims more f64 slots than the host offered
const TOO_MANY: &str = r#"(module (memory (export "memory") 1)
  (func (export "separate") (param i32 i32 i32 i32) (result i32) (i32.const 1000000)))"#;
/// a record that says nnz = 5 but only 4 slots were returned
const TRUNCATED: &str = r#"(module (memory (export "memory") 1)
  (func (export "separate") (param i32 i32 i32 i32) (result i32)
    (f64.store (local.get 2) (f64.const 5)) (i32.const 4)))"#;
/// a well-formed record whose column does not exist: rejected by the host cut gate, not by the ABI
const BAD_COLUMN: &str = r#"(module (memory (export "memory") 1)
  (func (export "separate") (param i32 i32 i32 i32) (result i32)
    (f64.store (local.get 2) (f64.const 1))
    (f64.store offset=8 (local.get 2) (f64.const 1))
    (f64.store offset=16 (local.get 2) (f64.const 77))
    (f64.store offset=24 (local.get 2) (f64.const 1))
    (i32.const 4)))"#;

/// max 10x0+11x1+12x2+3x3, 5x0+5x1+5x2+4x3 <= 14: x0+x1+x2 <= 2 is a cover cut; optimum 26.
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

fn solve(src: Option<(&str, WasmLimits)>) -> (f64, Option<WasmStats>, HookStats) {
    let mut h = Highs::default();
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_str("presolve", "off");
    h.set_option_value_double("mip_heuristic_effort", 0.0);
    h.set_option_value_bool("mip_heuristic_run_feasibility_jump", false);
    h.set_option_value_bool("mip_detect_symmetry", false);
    h.set_option_value_int("threads", 1);
    h.pass_model(model());
    let stats = src.map(|(s, limits)| {
        let hooks = WasmHooks::new("wasm-sep", s.as_bytes(), limits).expect("component");
        let st = hooks.stats();
        h.set_mip_hooks(Box::new(hooks));
        st
    });
    assert_eq!(h.run(), HighsStatus::kOk);
    assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
    let snap = stats.map(|s| s.lock().unwrap().clone());
    (
        h.get_info().objective_function_value,
        snap,
        h.mip_hook_stats(),
    )
}

fn small() -> WasmLimits {
    WasmLimits {
        fuel_per_call: 20_000,
        ..WasmLimits::default()
    }
}

#[test]
fn wasm_separator_adds_a_valid_cut_and_reaches_the_optimum() {
    let (plain, _, hs0) = solve(None);
    assert!((plain - 26.0).abs() < 1e-6, "{plain}");
    assert_eq!(hs0.user_cuts_offered, 0);
    let (obj, st, hs) = solve(Some((COVER_CUT, WasmLimits::default())));
    let st = st.unwrap();
    assert!((obj - 26.0).abs() < 1e-6, "{obj}");
    assert!(st.calls > 0 && st.ok == st.calls, "{st:?}");
    assert!(st.fuel_used > 0);
    // the built-in knapsack-cover separator may already have found the same cut (a duplicate)
    assert!(hs.user_cuts_offered >= 1, "{hs:?}");
    assert_eq!(
        hs.user_cuts_accepted + hs.user_cuts_duplicate,
        hs.user_cuts_offered,
        "{hs:?}"
    );
    assert_eq!(
        hs.user_cuts_rejected_malformed + hs.user_cuts_rejected_debug_solution,
        0
    );
}

#[test]
fn wasm_separator_trap_loop_and_bad_answers_fall_back() {
    for (name, src, field) in [
        ("traps", TRAPS, 0),
        ("loops", LOOPS, 1),
        ("too many", TOO_MANY, 2),
        ("truncated", TRUNCATED, 2),
    ] {
        let (obj, st, hs) = solve(Some((src, small())));
        let st = st.unwrap();
        assert!((obj - 26.0).abs() < 1e-6, "{name}: {obj}");
        assert!(st.calls > 0, "{name}: {st:?}");
        let bad = [st.traps, st.out_of_fuel, st.invalid][field];
        assert_eq!(bad, st.calls, "{name}: every call must fall back: {st:?}");
        assert_eq!(hs.user_cuts_offered, 0, "{name}: no cut may reach the gate");
    }
}

#[test]
fn wasm_separator_invalid_cut_is_rejected_by_the_gate() {
    let (obj, st, hs) = solve(Some((BAD_COLUMN, small())));
    let st = st.unwrap();
    assert!((obj - 26.0).abs() < 1e-6, "{obj}");
    assert!(st.ok > 0 && st.ok == st.calls, "{st:?}");
    assert!(hs.user_cuts_offered > 0, "{hs:?}");
    assert_eq!(hs.user_cuts_rejected_malformed, hs.user_cuts_offered);
    assert_eq!(hs.user_cuts_accepted, 0);
}

#[test]
fn wasm_separator_fuel_is_deterministic_and_loading_checks_the_signature() {
    let (_, a, _) = solve(Some((COVER_CUT, WasmLimits::default())));
    let (_, b, _) = solve(Some((COVER_CUT, WasmLimits::default())));
    assert_eq!(a, b);
    let bad = r#"(module (memory (export "memory") 1) (func (export "separate") (param i32) (result i32) (i32.const 0)))"#;
    let e = WasmHooks::new("x", bad.as_bytes(), WasmLimits::default())
        .err()
        .expect("must reject");
    assert!(e.contains("separate signature"), "{e}");
    let nomem =
        r#"(module (func (export "separate") (param i32 i32 i32 i32) (result i32) (i32.const 0)))"#;
    assert!(WasmHooks::new("x", nomem.as_bytes(), WasmLimits::default()).is_err());
}

#[test]
fn wasm_separator_returns_several_cuts_in_one_call() {
    let (obj, st, hs) = solve(Some((TWO_CUTS, small())));
    let st = st.unwrap();
    assert!((obj - 26.0).abs() < 1e-6, "{obj}");
    assert!(st.ok > 0 && st.ok == st.calls, "{st:?}");
    assert_eq!(hs.user_cuts_offered, 2 * st.calls, "{hs:?}");
    assert_eq!(hs.user_cuts_rejected_malformed, 0);
}
