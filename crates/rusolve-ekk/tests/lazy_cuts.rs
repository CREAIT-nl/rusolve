//! Lazy constraints (TSP subtour elimination) and
//! user cuts through `MipHooks::separate`, each checked against a model that lists the same rows.
//! The examples in `crates/rusolve-ekk/examples/` use the same code.

#[path = "../examples/common/mod.rs"]
mod common;

use common::*;
use rusolve_ekk::mip::hooks::*;
use rusolve_ekk::mip::lazy::LazyConstraints;
use rusolve_ekk::prelude::*;

#[test]
fn lazy_tsp_matches_brute_force_and_the_full_formulation() {
    let mut lazy_rows = 0;
    for (n, seed) in [
        (6usize, 1u64),
        (7, 2),
        (8, 3),
        (9, 4),
        (9, 5),
        (9, 6),
        (10, 7),
        (10, 8),
    ] {
        let pts = cities(n, seed);
        let brute = brute_force_tour(&pts);
        // lazy: degree rows only, subtour rows on demand
        let mut h = quiet(tsp_model(&pts, false));
        let mut oracle = SubtourOracle::new(n);
        let (st, rep) = h.run_with_lazy_constraints(&mut oracle, 50);
        assert_eq!(st, HighsStatus::kOk);
        assert!(rep.accepted, "{rep:?}");
        let obj = h.get_info().objective_function_value;
        assert!((obj - brute).abs() < 1e-6, "n={n} lazy {obj} brute {brute}");
        let tour = tour_components(&h.get_solution().col_value, n);
        assert_eq!(tour.len(), 1, "accepted solution must be one tour");
        // the full formulation lists every subtour row up front
        let mut f = quiet(tsp_model(&pts, true));
        assert_eq!(f.run(), HighsStatus::kOk);
        let full = f.get_info().objective_function_value;
        assert!((full - brute).abs() < 1e-6, "full {full} brute {brute}");
        lazy_rows += rep.rows_added;
        eprintln!(
            "n={n} seed={seed} rounds={} rows_added={} tour={obj:.4}",
            rep.rounds, rep.rows_added
        );
    }
    // the loop must have had real work to do: some instance needed subtour rows
    assert!(lazy_rows > 0);
}

#[test]
fn lazy_loop_drops_malformed_and_satisfied_rows_and_always_terminates() {
    let pts = cities(6, 5);
    let n = 6;
    let mut h = quiet(tsp_model(&pts, false));
    let mut inner = SubtourOracle::new(n);
    let mut calls = 0;
    let mut oracle = |sol: &[f64]| {
        calls += 1;
        let mut v = inner.lazy_constraints(sol);
        if v.is_empty() {
            return v; // an empty answer accepts; junk is only mixed into real rows
        }
        v.push(UserCut {
            indices: vec![0, 0],
            values: vec![1.0, 1.0],
            rhs: 1.0,
        }); // duplicate column
        v.push(UserCut {
            indices: vec![999],
            values: vec![1.0],
            rhs: 1.0,
        }); // out of range
        v.push(UserCut {
            indices: vec![0],
            values: vec![1.0],
            rhs: 5.0,
        }); // satisfied: x0 <= 5
        v
    };
    let (_, rep) = h.run_with_lazy_constraints(&mut oracle, 50);
    assert!(rep.accepted);
    assert_eq!(rep.rows_rejected_malformed, 2 * (rep.rounds - 1));
    assert_eq!(rep.rows_not_violated, rep.rounds - 1);
    assert_eq!(calls, rep.rounds);
    // round limit is reported, not hidden
    let mut h = quiet(tsp_model(&cities(8, 3), false));
    let mut never = |_: &[f64]| {
        vec![UserCut {
            indices: vec![0],
            values: vec![1.0],
            rhs: -1.0,
        }]
    };
    let (_, rep) = h.run_with_lazy_constraints(&mut never, 2);
    assert!(!rep.accepted && rep.rounds <= 2, "{rep:?}");
}

#[test]
fn lazy_oracle_returning_only_malformed_rows_is_not_accepted() {
    let n = 9;
    let pts = cities(n, 4);
    let mut h = quiet(tsp_model(&pts, false));
    let mut oracle = |_: &[f64]| {
        vec![UserCut {
            indices: vec![999],
            values: vec![1.0],
            rhs: 1.0,
        }]
    };
    let (_, rep) = h.run_with_lazy_constraints(&mut oracle, 50);
    assert!(!rep.accepted, "{rep:?}");
    assert_eq!(rep.rounds, 1);
    assert_eq!(rep.rows_rejected_malformed, 1);
    assert!(tour_components(&h.get_solution().col_value, n).len() > 1);
}

#[test]
fn user_cover_cuts_equal_the_same_cuts_in_the_model() {
    let (values, weights, cap) = knapsack_data();
    let covers = minimal_covers(&weights, cap);
    assert!(covers.len() > 20, "{} covers", covers.len());
    let plain = {
        let mut h = quiet_with(knapsack_model(&values, &weights, cap, &[]), false);
        assert_eq!(h.run(), HighsStatus::kOk);
        h.get_info().objective_function_value
    };
    // cuts as model rows
    let in_model = {
        let mut h = quiet_with(knapsack_model(&values, &weights, cap, &covers), false);
        assert_eq!(h.run(), HighsStatus::kOk);
        h.get_info().objective_function_value
    };
    // the same cuts from the user separator
    let mut h = quiet_with(knapsack_model(&values, &weights, cap, &[]), false);
    h.set_mip_hooks(Box::new(CoverSeparator {
        covers: covers.clone(),
    }));
    assert_eq!(h.run(), HighsStatus::kOk);
    let hooked = h.get_info().objective_function_value;
    assert!((plain - in_model).abs() < 1e-6, "{plain} {in_model}");
    assert!((plain - hooked).abs() < 1e-6, "{plain} {hooked}");
    let st = h.mip_hook_stats();
    assert!(
        st.separator_calls > 0 && st.user_cuts_accepted > 0,
        "{st:?}"
    );
    assert_eq!(
        st.user_cuts_rejected_malformed + st.user_cuts_rejected_debug_solution,
        0
    );
}

struct CoverSeparator {
    covers: Vec<Vec<usize>>,
}
impl MipHooks for CoverSeparator {
    fn name(&self) -> &str {
        "cover-separator"
    }
    fn separate(&mut self, _info: &SeparationInfo<'_>) -> Vec<UserCut> {
        cover_cuts(&self.covers)
    }
}
