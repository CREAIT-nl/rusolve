//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! end-to-end tests for `mip_symmetry_dynamic` on graph-coloring and
//! bin-packing instances with full orbitopes.

use super::*;
use crate::parallel::highs_parallel::initialize_scheduler;
use crate::parallel::highs_task_executor::HighsTaskExecutor;

// ---------------------------------------------------------------------------
// Graph-coloring MIP builder: 5 nodes, 4 colors, 5-cycle + one chord.
// ---------------------------------------------------------------------------

const GC_N: usize = 5;
const GC_C: usize = 4;

/// Edges of the 5-cycle plus one chord (0-2).
const GC_EDGES: [(usize, usize); 6] = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 2)];

/// Build the graph-coloring MIP.
/// Variables: x[v][c] for v in 0..5, c in 0..4 (20 vars), then y[c] for c in 0..4 (4 vars).
/// Total 24 columns.
/// Rows: 5 (one per node: sum_c x[v][c] = 1) + 6*4 = 24 (edge-color: x[u][c]+x[v][c] <= y[c]).
/// Total 29 rows.
fn graph_coloring_lp() -> HighsLp {
    let n_vars = GC_N * GC_C + GC_C; // 24
    let n_rows = GC_N + GC_EDGES.len() * GC_C; // 5 + 24 = 29

    // Build CSR matrix.
    let mut start: Vec<HighsInt> = vec![0; n_vars + 1];
    let mut index: Vec<HighsInt> = Vec::new();
    let mut value: Vec<f64> = Vec::new();

    // For each column, collect (row, value) pairs.
    let mut col_entries: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_vars];

    // Node rows: row v (0..5): sum_c x[v][c] = 1
    for v in 0..GC_N {
        for col_c in 0..GC_C {
            let col = v * GC_C + col_c;
            col_entries[col].push((v, 1.0));
        }
    }

    // Edge-color rows: for edge e (index in GC_EDGES), color c:
    // row = GC_N + e * GC_C + c
    // x[u][c] + x[v][c] - y[c] <= 0
    for (e, &(u, v)) in GC_EDGES.iter().enumerate() {
        for col_c in 0..GC_C {
            let row = GC_N + e * GC_C + col_c;
            let col_u = u * GC_C + col_c;
            let col_v = v * GC_C + col_c;
            let col_y = GC_N * GC_C + col_c;
            col_entries[col_u].push((row, 1.0));
            col_entries[col_v].push((row, 1.0));
            col_entries[col_y].push((row, -1.0));
        }
    }

    // Build CSR from col_entries.
    let mut pos: HighsInt = 0;
    for j in 0..n_vars {
        start[j] = pos;
        for &(row, val) in &col_entries[j] {
            index.push(row as HighsInt);
            value.push(val);
            pos += 1;
        }
    }
    start[n_vars] = pos;

    // Objective: min sum_c y_c
    let mut col_cost = vec![0.0f64; n_vars];
    for col_c in 0..GC_C {
        col_cost[GC_N * GC_C + col_c] = 1.0;
    }

    // Bounds
    let col_lower = vec![0.0f64; n_vars];
    let col_upper = vec![1.0f64; n_vars];

    // Row bounds
    let mut row_lower = vec![0.0f64; n_rows];
    let mut row_upper = vec![0.0f64; n_rows];
    // Node rows: = 1
    for v in 0..GC_N {
        row_lower[v] = 1.0;
        row_upper[v] = 1.0;
    }
    // Edge-color rows: <= 0
    for e in 0..GC_EDGES.len() {
        for col_c in 0..GC_C {
            let row = GC_N + e * GC_C + col_c;
            row_lower[row] = -kHighsInf;
            row_upper[row] = 0.0;
        }
    }

    let integrality = vec![HighsVarType::kInteger; n_vars];

    let mut lp = HighsLp::default();
    lp.num_col_ = n_vars as HighsInt;
    lp.num_row_ = n_rows as HighsInt;
    lp.col_cost_ = col_cost;
    lp.col_lower_ = col_lower;
    lp.col_upper_ = col_upper;
    lp.row_lower_ = row_lower;
    lp.row_upper_ = row_upper;
    lp.integrality_ = integrality;
    lp.a_matrix_.num_col_ = n_vars as HighsInt;
    lp.a_matrix_.num_row_ = n_rows as HighsInt;
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    lp
}

/// Brute-force optimum for graph coloring: enumerate all 4^5 assignments.
/// Each assignment maps node v -> color c. Check edge constraints.
/// Objective = number of colors used.
fn gc_brute_force_optimum() -> f64 {
    let mut best = GC_C as f64;
    for mask in 0..GC_C.pow(GC_N as u32) {
        let mut colors = [0usize; GC_N];
        let mut tmp = mask;
        for v in 0..GC_N {
            colors[v] = tmp % GC_C;
            tmp /= GC_C;
        }
        // Check edges
        let mut feasible = true;
        for &(u, v) in GC_EDGES.iter() {
            if colors[u] == colors[v] {
                feasible = false;
                break;
            }
        }
        if !feasible {
            continue;
        }
        // Count distinct colors used
        let mut used = [false; GC_C];
        for v in 0..GC_N {
            used[colors[v]] = true;
        }
        let obj = used.iter().filter(|&&b| b).count() as f64;
        if obj < best {
            best = obj;
        }
    }
    best
}

/// Feasibility check for graph coloring solution.
fn gc_check_feasible(sol: &[f64]) -> bool {
    let n_vars = GC_N * GC_C + GC_C;
    if sol.len() < n_vars {
        return false;
    }
    // Check integrality
    for &v in &sol[..n_vars] {
        if (v - v.round()).abs() > 1e-6 {
            return false;
        }
    }
    // Check bounds [0,1]
    for &v in &sol[..n_vars] {
        if v < -1e-6 || v > 1.0 + 1e-6 {
            return false;
        }
    }
    // Check node rows: sum_c x[v][c] = 1
    for v in 0..GC_N {
        let s: f64 = (0..GC_C).map(|col_c| sol[v * GC_C + col_c]).sum();
        if (s - 1.0).abs() > 1e-6 {
            return false;
        }
    }
    // Check edge-color rows: x[u][c] + x[v][c] - y[c] <= 0
    for (e, &(u, v)) in GC_EDGES.iter().enumerate() {
        for col_c in 0..GC_C {
            let col_u = u * GC_C + col_c;
            let col_v = v * GC_C + col_c;
            let col_y = GC_N * GC_C + col_c;
            let lhs = sol[col_u] + sol[col_v] - sol[col_y];
            if lhs > 1e-6 {
                return false;
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Bin-packing MIP builder: 6 items, 4 bins.
// ---------------------------------------------------------------------------

const BP_N: usize = 6;
const BP_B: usize = 4;
const BP_CAP: f64 = 10.0;
const BP_WEIGHTS: [f64; BP_N] = [3.0, 4.0, 5.0, 2.0, 3.0, 4.0];

/// Build the bin-packing MIP.
/// Variables: x[i][b] for i in 0..6, b in 0..4 (24 vars), then y[b] for b in 0..4 (4 vars).
/// Total 28 columns.
/// Rows: 6 (item assignment: sum_b x[i][b] = 1) + 4 (bin capacity: sum_i w_i x[i][b] - C y[b] <= 0).
/// Total 10 rows.
fn bin_packing_lp() -> HighsLp {
    let n_vars = BP_N * BP_B + BP_B; // 28
    let n_rows = BP_N + BP_B; // 10

    let mut col_entries: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_vars];

    // Item rows: row i (0..6): sum_b x[i][b] = 1
    for i in 0..BP_N {
        for b in 0..BP_B {
            let col = i * BP_B + b;
            col_entries[col].push((i, 1.0));
        }
    }

    // Bin capacity rows: row BP_N + b: sum_i w_i x[i][b] - C y[b] <= 0
    for b in 0..BP_B {
        let row = BP_N + b;
        for i in 0..BP_N {
            let col = i * BP_B + b;
            col_entries[col].push((row, BP_WEIGHTS[i]));
        }
        let col_y = BP_N * BP_B + b;
        col_entries[col_y].push((row, -BP_CAP));
    }

    // Build CSR.
    let mut start: Vec<HighsInt> = vec![0; n_vars + 1];
    let mut index: Vec<HighsInt> = Vec::new();
    let mut value: Vec<f64> = Vec::new();

    let mut pos: HighsInt = 0;
    for j in 0..n_vars {
        start[j] = pos;
        for &(row, val) in &col_entries[j] {
            index.push(row as HighsInt);
            value.push(val);
            pos += 1;
        }
    }
    start[n_vars] = pos;

    // Objective: min sum_b y_b
    let mut col_cost = vec![0.0f64; n_vars];
    for b in 0..BP_B {
        col_cost[BP_N * BP_B + b] = 1.0;
    }

    let col_lower = vec![0.0f64; n_vars];
    let col_upper = vec![1.0f64; n_vars];

    let mut row_lower = vec![0.0f64; n_rows];
    let mut row_upper = vec![0.0f64; n_rows];
    // Item rows: = 1
    for i in 0..BP_N {
        row_lower[i] = 1.0;
        row_upper[i] = 1.0;
    }
    // Bin capacity rows: <= 0
    for b in 0..BP_B {
        let row = BP_N + b;
        row_lower[row] = -kHighsInf;
        row_upper[row] = 0.0;
    }

    let integrality = vec![HighsVarType::kInteger; n_vars];

    let mut lp = HighsLp::default();
    lp.num_col_ = n_vars as HighsInt;
    lp.num_row_ = n_rows as HighsInt;
    lp.col_cost_ = col_cost;
    lp.col_lower_ = col_lower;
    lp.col_upper_ = col_upper;
    lp.row_lower_ = row_lower;
    lp.row_upper_ = row_upper;
    lp.integrality_ = integrality;
    lp.a_matrix_.num_col_ = n_vars as HighsInt;
    lp.a_matrix_.num_row_ = n_rows as HighsInt;
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    lp
}

/// Brute-force optimum for bin packing: enumerate all 4^6 assignments.
/// Each assignment maps item i -> bin b. Check capacity.
/// Objective = number of bins used.
fn bp_brute_force_optimum() -> f64 {
    let mut best = BP_B as f64;
    for mask in 0..BP_B.pow(BP_N as u32) {
        let mut bins = [0usize; BP_N];
        let mut tmp = mask;
        for i in 0..BP_N {
            bins[i] = tmp % BP_B;
            tmp /= BP_B;
        }
        // Check capacity
        let mut load = [0.0f64; BP_B];
        for i in 0..BP_N {
            load[bins[i]] += BP_WEIGHTS[i];
        }
        let mut feasible = true;
        for b in 0..BP_B {
            if load[b] > BP_CAP + 1e-9 {
                feasible = false;
                break;
            }
        }
        if !feasible {
            continue;
        }
        // Count bins used
        let mut used = [false; BP_B];
        for i in 0..BP_N {
            used[bins[i]] = true;
        }
        let obj = used.iter().filter(|&&b| b).count() as f64;
        if obj < best {
            best = obj;
        }
    }
    best
}

/// Feasibility check for bin-packing solution.
fn bp_check_feasible(sol: &[f64]) -> bool {
    let n_vars = BP_N * BP_B + BP_B;
    if sol.len() < n_vars {
        return false;
    }
    // Check integrality
    for &v in &sol[..n_vars] {
        if (v - v.round()).abs() > 1e-6 {
            return false;
        }
    }
    // Check bounds [0,1]
    for &v in &sol[..n_vars] {
        if v < -1e-6 || v > 1.0 + 1e-6 {
            return false;
        }
    }
    // Check item rows: sum_b x[i][b] = 1
    for i in 0..BP_N {
        let s: f64 = (0..BP_B).map(|b| sol[i * BP_B + b]).sum();
        if (s - 1.0).abs() > 1e-6 {
            return false;
        }
    }
    // Check bin capacity rows: sum_i w_i x[i][b] - C y[b] <= 0
    for b in 0..BP_B {
        let mut lhs = 0.0f64;
        for i in 0..BP_N {
            lhs += BP_WEIGHTS[i] * sol[i * BP_B + b];
        }
        let col_y = BP_N * BP_B + b;
        lhs -= BP_CAP * sol[col_y];
        if lhs > 1e-6 {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Solve helper
// ---------------------------------------------------------------------------

struct SolveResult {
    modelstatus: HighsModelStatus,
    objective: f64,
    solution: Vec<f64>,
    node_count: i64,
}

fn solve_with_flag(lp: &HighsLp, flag: bool) -> SolveResult {
    let lp = lp.clone();
    std::thread::spawn(move || {
        initialize_scheduler(1);
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.presolve = kHighsOffString.to_string();
        options.mip_heuristic_effort = 0.0;
        options.mip_heuristic_run_feasibility_jump = false;
        options.mip_heuristic_run_rens = false;
        options.mip_heuristic_run_rins = false;
        options.mip_heuristic_run_root_reduced_cost = false;
        options.mip_heuristic_run_zi_round = false;
        options.mip_heuristic_run_shifting = false;
        options.mip_symmetry_dynamic = flag;
        let solution = HighsSolution::default();
        let mut callback = HighsCallback::default();
        let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
        ms.run();
        let out = SolveResult {
            modelstatus: ms.modelstatus_,
            objective: ms.solution_objective_,
            solution: ms.solution_.clone(),
            node_count: ms.node_count_,
        };
        drop(ms);
        HighsTaskExecutor::shutdown(true);
        out
    })
    .join()
    .expect("solve thread panicked")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn s3_e2e_graph_coloring_flag_off_matches_brute_force() {
    let lp = graph_coloring_lp();
    let brute = gc_brute_force_optimum();
    let res = solve_with_flag(&lp, false);
    assert_eq!(res.modelstatus, HighsModelStatus::kOptimal);
    assert!(
        (res.objective - brute).abs() < 1e-6,
        "flag off: objective {} != brute force {}",
        res.objective,
        brute
    );
    assert!(
        gc_check_feasible(&res.solution),
        "flag off: infeasible solution"
    );
    eprintln!(
        "s3_e2e_graph_coloring_flag_off: obj={}, nodes={}",
        res.objective, res.node_count
    );
}

#[test]
fn s3_e2e_graph_coloring_flag_on_matches_brute_force() {
    let lp = graph_coloring_lp();
    let brute = gc_brute_force_optimum();
    let res_off = solve_with_flag(&lp, false);
    let res_on = solve_with_flag(&lp, true);
    assert_eq!(res_on.modelstatus, HighsModelStatus::kOptimal);
    assert!(
        (res_on.objective - brute).abs() < 1e-6,
        "flag on: objective {} != brute force {}",
        res_on.objective,
        brute
    );
    assert!(
        gc_check_feasible(&res_on.solution),
        "flag on: infeasible solution"
    );
    assert!(
        (res_on.objective - res_off.objective).abs() < 1e-6,
        "flag on objective {} != flag off objective {}",
        res_on.objective,
        res_off.objective
    );
    // Reporting only: flag-on node count vs flag off.
    eprintln!(
        "s3_e2e_graph_coloring_flag_on: obj={}, nodes_on={}, nodes_off={}, on_minus_off={}",
        res_on.objective,
        res_on.node_count,
        res_off.node_count,
        res_on.node_count - res_off.node_count
    );
}

#[test]
fn s3_e2e_bin_packing_flag_off_matches_brute_force() {
    let lp = bin_packing_lp();
    let brute = bp_brute_force_optimum();
    let res = solve_with_flag(&lp, false);
    assert_eq!(res.modelstatus, HighsModelStatus::kOptimal);
    assert!(
        (res.objective - brute).abs() < 1e-6,
        "flag off: objective {} != brute force {}",
        res.objective,
        brute
    );
    assert!(
        bp_check_feasible(&res.solution),
        "flag off: infeasible solution"
    );
    eprintln!(
        "s3_e2e_bin_packing_flag_off: obj={}, nodes={}",
        res.objective, res.node_count
    );
}

#[test]
fn s3_e2e_bin_packing_flag_on_matches_brute_force() {
    let lp = bin_packing_lp();
    let brute = bp_brute_force_optimum();
    let res_off = solve_with_flag(&lp, false);
    let res_on = solve_with_flag(&lp, true);
    assert_eq!(res_on.modelstatus, HighsModelStatus::kOptimal);
    assert!(
        (res_on.objective - brute).abs() < 1e-6,
        "flag on: objective {} != brute force {}",
        res_on.objective,
        brute
    );
    assert!(
        bp_check_feasible(&res_on.solution),
        "flag on: infeasible solution"
    );
    assert!(
        (res_on.objective - res_off.objective).abs() < 1e-6,
        "flag on objective {} != flag off objective {}",
        res_on.objective,
        res_off.objective
    );
    // Reporting only: flag-on node count vs flag off.
    eprintln!(
        "s3_e2e_bin_packing_flag_on: obj={}, nodes_on={}, nodes_off={}, on_minus_off={}",
        res_on.objective,
        res_on.node_count,
        res_off.node_count,
        res_on.node_count - res_off.node_count
    );
}
