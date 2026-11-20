//! rusolve MIP hook points, the ported HiGHS
//! callback API and the validity gates. See docs/design/hooks.md.
//!
//! Model: a multi-dimensional 0/1 knapsack (maximise) that needs a real tree with presolve and
//! heuristics off, so every hook point is reached.

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::hooks::*;
use rusolve_ekk::prelude::*;
use std::sync::{Arc, Mutex};

const N: usize = 28;
const M: usize = 4;

fn lcg(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

/// weights[i][j] for row i, item j, and the capacity of each row.
fn data() -> (Vec<f64>, Vec<Vec<f64>>, Vec<f64>) {
    let mut s = 12345u64;
    let values: Vec<f64> = (0..N).map(|_| (20 + lcg(&mut s) % 60) as f64).collect();
    let mut weights = vec![vec![0.0; N]; M];
    for row in weights.iter_mut() {
        for w in row.iter_mut() {
            *w = (10 + lcg(&mut s) % 50) as f64;
        }
    }
    let caps = weights
        .iter()
        .map(|r| (r.iter().sum::<f64>() * 0.35).floor())
        .collect();
    (values, weights, caps)
}

#[allow(clippy::field_reassign_with_default)]
fn model() -> HighsLp {
    let (values, weights, caps) = data();
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

#[derive(Debug, Clone, PartialEq)]
struct Outcome {
    status: HighsModelStatus,
    objective: f64,
    nodes: i64,
    lp_iters: i64,
}

fn outcome(h: &Highs) -> Outcome {
    let info = h.get_info();
    Outcome {
        status: h.get_model_status(),
        objective: info.objective_function_value,
        nodes: info.mip_node_count,
        lp_iters: info.simplex_iteration_count as i64,
    }
}

fn baseline() -> (Outcome, Vec<f64>) {
    let mut h = new_highs();
    assert_eq!(h.run(), HighsStatus::kOk);
    (outcome(&h), h.get_solution().col_value.clone())
}

/// A debug-solution file in the HiGHS solution format (read by column position).
fn write_debug_solution(sol: &[f64], tag: &str) -> String {
    let path =
        std::env::temp_dir().join(format!("hooks_debug_sol_{tag}_{}.sol", std::process::id()));
    let mut text = String::from(
        "Model status\nOptimal\n\n# Primal solution values\nFeasible\n\n# Columns 28\n",
    );
    for (j, v) in sol.iter().enumerate() {
        text.push_str(&format!("C{j} {v}\n"));
    }
    text.push_str("\n# Rows 0\n");
    std::fs::write(&path, text).unwrap();
    path.to_string_lossy().into_owned()
}

#[test]
fn hooks_none_equals_baseline() {
    let (base, _) = baseline();
    assert_eq!(base.status, HighsModelStatus::kOptimal);
    assert!(
        base.nodes > 10,
        "model must need a tree, nodes {}",
        base.nodes
    );
    // never installed == installed then cleared == installed with all-default (no-op) methods
    let mut h = new_highs();
    h.set_mip_hooks(Box::new(Noop));
    h.clear_mip_hooks();
    assert_eq!(h.run(), HighsStatus::kOk);
    assert_eq!(outcome(&h), base);
    let mut h = new_highs();
    h.set_mip_hooks(Box::new(Noop));
    assert_eq!(h.run(), HighsStatus::kOk);
    assert_eq!(
        outcome(&h),
        base,
        "default no-op methods must not move the search"
    );
    let st = h.mip_hook_stats();
    assert!(st.separator_calls > 0 && st.node_select_calls > 0 && st.cut_score_calls > 0);
    assert_eq!(st.branching_overrides + st.node_select_overrides, 0);
}

struct Noop;
impl MipHooks for Noop {}

#[test]
fn hooks_interrupt_callback_stops_with_interrupted_status() {
    let mut h = new_highs();
    let calls = Arc::new(Mutex::new(0usize));
    let cnt = Arc::clone(&calls);
    h.set_callback(HighsCallbackFunctionType::new(
        move |ty, _msg, _out, data_in| {
            if ty == kCallbackMipInterrupt {
                *cnt.lock().unwrap() += 1;
                data_in.user_interrupt = true;
            }
        },
    ));
    h.start_callback(kCallbackMipInterrupt);
    assert_eq!(h.run(), HighsStatus::kWarning);
    assert_eq!(h.get_model_status(), HighsModelStatus::kInterrupt);
    assert!(*calls.lock().unwrap() >= 1);
}

#[test]
fn hooks_user_solution_callback_reduces_nodes() {
    let (base, opt) = baseline();
    let mut h = new_highs();
    let offered = Arc::new(Mutex::new(false));
    let o = Arc::clone(&offered);
    h.set_callback(HighsCallbackFunctionType::new(
        move |ty, _m, _out, data_in| {
            if ty == kCallbackMipUserSolution {
                let mut done = o.lock().unwrap();
                if !*done {
                    *done = true;
                    data_in.user_has_solution = true;
                    data_in.user_solution = opt.clone();
                }
            }
        },
    ));
    h.start_callback(kCallbackMipUserSolution);
    assert_eq!(h.run(), HighsStatus::kOk);
    let with = outcome(&h);
    assert_eq!(with.status, HighsModelStatus::kOptimal);
    assert!((with.objective - base.objective).abs() < 1e-6);
    assert!(
        with.nodes < base.nodes,
        "nodes {} vs baseline {}",
        with.nodes,
        base.nodes
    );
    let st = h.mip_hook_stats();
    assert_eq!(st.user_solutions_checked_feasible, 1);
    assert_eq!(st.user_solutions_rejected, 0);
}

#[test]
fn hooks_infeasible_user_solution_is_rejected() {
    let (base, _) = baseline();
    let mut h = new_highs();
    h.set_callback(HighsCallbackFunctionType::new(
        move |ty, _m, _out, data_in| {
            if ty == kCallbackMipUserSolution {
                data_in.user_has_solution = true;
                data_in.user_solution = vec![1.0; N]; // violates every knapsack row
            }
        },
    ));
    h.start_callback(kCallbackMipUserSolution);
    assert_eq!(h.run(), HighsStatus::kOk);
    let with = outcome(&h);
    assert_eq!(with.status, HighsModelStatus::kOptimal);
    assert!(
        (with.objective - base.objective).abs() < 1e-6,
        "an infeasible solution became the incumbent"
    );
    let st = h.mip_hook_stats();
    assert!(st.user_solutions_rejected >= 1);
    assert_eq!(st.user_solutions_checked_feasible, 0);
}

/// Scores columns by index: reverses the HiGHS priorities.
struct ByIndex {
    high_first: bool,
}
impl MipHooks for ByIndex {
    fn name(&self) -> &str {
        "by_index"
    }
    fn branching_score(&mut self, info: &BranchingInfo<'_>) -> Option<Vec<f64>> {
        let sign = if self.high_first { 1.0 } else { -1.0 };
        Some(
            info.candidates
                .iter()
                .map(|cand| sign * cand.col as f64)
                .collect(),
        )
    }
}

#[test]
fn hooks_branching_score_changes_nodes_not_optimum() {
    let (base, _) = baseline();
    let mut counts = vec![];
    for high_first in [true, false] {
        let mut h = new_highs();
        h.set_mip_hooks(Box::new(ByIndex { high_first }));
        assert_eq!(h.run(), HighsStatus::kOk);
        let o = outcome(&h);
        assert_eq!(o.status, HighsModelStatus::kOptimal);
        assert!((o.objective - base.objective).abs() < 1e-6);
        assert!(h.mip_hook_stats().branching_overrides > 0);
        counts.push(o.nodes);
    }
    assert!(
        counts.iter().any(|&n| n != base.nodes),
        "branching hook did not change the node count: {counts:?} vs {}",
        base.nodes
    );
}

/// Always dives into the deepest open node (a depth-first order).
struct DeepestFirst;
impl MipHooks for DeepestFirst {
    fn select_node(&mut self, open: &[OpenNodeInfo]) -> Option<usize> {
        (0..open.len()).max_by_key(|&i| open[i].depth)
    }
}

#[test]
fn hooks_node_selection_changes_nodes_not_optimum() {
    let (base, _) = baseline();
    let mut h = new_highs();
    h.set_mip_hooks(Box::new(DeepestFirst));
    assert_eq!(h.run(), HighsStatus::kOk);
    let o = outcome(&h);
    assert_eq!(o.status, HighsModelStatus::kOptimal);
    assert!((o.objective - base.objective).abs() < 1e-6);
    assert!(h.mip_hook_stats().node_select_overrides > 0);
    assert_ne!(
        o.nodes, base.nodes,
        "node-selection hook did not change the search"
    );
}

/// Prefers the cut with the smallest default score.
struct WorstCutFirst;
impl MipHooks for WorstCutFirst {
    fn cut_score(&mut self, info: &CutScoreInfo) -> Option<f64> {
        Some(-info.default_score)
    }
}

#[test]
fn hooks_cut_score_runs_and_keeps_optimum() {
    let (base, _) = baseline();
    let mut h = new_highs();
    h.set_mip_hooks(Box::new(WorstCutFirst));
    assert_eq!(h.run(), HighsStatus::kOk);
    let o = outcome(&h);
    assert_eq!(o.status, HighsModelStatus::kOptimal);
    assert!((o.objective - base.objective).abs() < 1e-6);
    assert!(h.mip_hook_stats().cut_score_calls > 0);
}

/// Returns one fixed cut on every call.
struct OneCut(UserCut);
impl MipHooks for OneCut {
    fn name(&self) -> &str {
        "one_cut"
    }
    fn separate(&mut self, _info: &SeparationInfo<'_>) -> Vec<UserCut> {
        vec![self.0.clone()]
    }
}

/// A knapsack cover cut of row 0: a greedy cover C (sum of weights > capacity) gives
/// sum_{j in C} x_j <= |C| - 1, valid for every feasible point.
fn cover_cut() -> UserCut {
    let (_, weights, caps) = data();
    let mut order: Vec<usize> = (0..N).collect();
    order.sort_by(|&a, &b| weights[0][b].partial_cmp(&weights[0][a]).unwrap());
    let (mut cover, mut sum) = (Vec::new(), 0.0);
    for j in order {
        cover.push(j);
        sum += weights[0][j];
        if sum > caps[0] {
            break;
        }
    }
    assert!(sum > caps[0]);
    UserCut {
        indices: cover.iter().map(|&j| j as HighsInt).collect(),
        values: vec![1.0; cover.len()],
        rhs: cover.len() as f64 - 1.0,
    }
}

#[test]
fn hooks_user_separator_valid_cut_is_accepted() {
    let (base, _) = baseline();
    let mut h = new_highs();
    h.set_mip_hooks(Box::new(OneCut(cover_cut())));
    assert_eq!(h.run(), HighsStatus::kOk);
    let o = outcome(&h);
    assert_eq!(o.status, HighsModelStatus::kOptimal);
    assert!((o.objective - base.objective).abs() < 1e-6);
    let st = h.mip_hook_stats();
    assert!(st.user_cuts_accepted >= 1, "{st:?}");
    assert_eq!(
        st.user_cuts_rejected_debug_solution + st.user_cuts_rejected_malformed,
        0
    );
}

#[test]
fn hooks_user_separator_invalid_cut_is_flagged_by_debug_solution_gate() {
    let (base, opt) = baseline();
    // sum_j x_j <= |opt| - 1 cuts off the optimum (and only the optimum-sized sets)
    let card = opt.iter().filter(|&&v| v > 0.5).count();
    let bad = UserCut {
        indices: (0..N as HighsInt).collect(),
        values: vec![1.0; N],
        rhs: card as f64 - 1.0,
    };
    let path = write_debug_solution(&opt, "invalid_cut");
    let mut h = new_highs();
    h.set_option_value_str("mip_debug_solution_file", &path);
    h.set_mip_hooks(Box::new(OneCut(bad)));
    assert_eq!(h.run(), HighsStatus::kOk);
    let st = h.mip_hook_stats();
    assert!(st.user_cuts_rejected_debug_solution >= 1, "{st:?}");
    assert_eq!(st.user_cuts_accepted, 0);
    // rejected, so it never reached the LP: the optimum is intact
    assert!((outcome(&h).objective - base.objective).abs() < 1e-6);
    let _ = std::fs::remove_file(path);
}

#[test]
fn hooks_user_separator_valid_cut_passes_debug_solution_gate() {
    let (base, opt) = baseline();
    let path = write_debug_solution(&opt, "valid_cut");
    let mut h = new_highs();
    h.set_option_value_str("mip_debug_solution_file", &path);
    h.set_mip_hooks(Box::new(OneCut(cover_cut())));
    assert_eq!(h.run(), HighsStatus::kOk);
    let st = h.mip_hook_stats();
    assert!(
        st.user_cuts_accepted >= 1 && st.user_cuts_rejected_debug_solution == 0,
        "{st:?}"
    );
    assert!((outcome(&h).objective - base.objective).abs() < 1e-6);
    let _ = std::fs::remove_file(path);
}

#[test]
fn hooks_malformed_user_cut_is_rejected() {
    let mut h = new_highs();
    h.set_mip_hooks(Box::new(OneCut(UserCut {
        indices: vec![0, 9999],
        values: vec![1.0, f64::NAN],
        rhs: 1.0,
    })));
    assert_eq!(h.run(), HighsStatus::kOk);
    let st = h.mip_hook_stats();
    assert!(st.user_cuts_rejected_malformed >= 1);
    assert_eq!(st.user_cuts_accepted, 0);
}
