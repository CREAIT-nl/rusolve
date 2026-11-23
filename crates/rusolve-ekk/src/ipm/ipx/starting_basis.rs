//! HiGHS source: ipm/ipx/starting_basis.h, starting_basis.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::basis::{BASIC, BASIC_FREE, Basis, NONBASIC, NONBASIC_FIXED};
use super::iterate::{Iterate, State};
use super::sparse_matrix::scatter_column;
use super::timer::Timer;
use super::*;

/// HiGHS: ipx/starting_basis.cc::AssertConsistency (static, lines 11-48)
///
/// Only called under `assert`; call it inside `if cfg!(debug_assertions)`.
fn assert_consistency(iterate: &Iterate, basis: &Basis) {
    let model = basis.model();
    let m = model.rows();
    let n = model.cols();
    let lb = model.lb();
    let ub = model.ub();

    for j in 0..n + m {
        let ju = j as usize;
        if lb[ju] == ub[ju] {
            match iterate.state_of(j) {
                State::fixed => debug_assert!(basis.status_of(j) == NONBASIC_FIXED),
                State::free => {
                    debug_assert!(j >= n); // must be slack variable
                    debug_assert!(basis.status_of(j) == BASIC_FREE);
                }
                _ => debug_assert!(false),
            }
        } else if lb[ju].is_infinite() && ub[ju].is_infinite() {
            match iterate.state_of(j) {
                State::fixed => debug_assert!(basis.status_of(j) == NONBASIC_FIXED),
                State::free => debug_assert!(basis.status_of(j) == BASIC_FREE),
                _ => debug_assert!(false),
            }
        } else {
            debug_assert!(matches!(iterate.state_of(j), State::barrier));
            debug_assert!(basis.status_of(j) == NONBASIC || basis.status_of(j) == BASIC);
        }
    }
}

/// HiGHS: ipx/starting_basis.cc::PostprocessDependencies (static, lines 50-125)
fn postprocess_dependencies(iterate: &mut Iterate, p_basis: &mut Basis, info: &mut Info) {
    // Rc clone: `iterate` is mutated by `update` below, so the model cannot
    // be borrowed out of it for the whole function.
    let model = iterate.model_rc();
    let m = model.rows();
    let n = model.cols();
    let ai = model.AI();
    let lb = model.lb();
    let ub = model.ub();
    let mut dependent_rows: Vec<Int> = Vec::new();
    let mut dependent_cols: Vec<Int> = Vec::new();
    let mut dx = Vector::new(n + m);
    let mut dy = Vector::new(m);

    // Postprocess linearly dependent columns. If a free variable could not be
    // made basic, then the variable is redundant and can be fixed at an
    // arbitrary value. We fix it at zero and adapt the basic variables to keep
    // AI*x unchanged.
    if info.dependent_cols > 0 {
        let mut dxbasic = Vector::new(m);
        for j in 0..n {
            let ju = j as usize;
            if lb[ju].is_infinite() && ub[ju].is_infinite() && p_basis.is_nonbasic(j) {
                debug_assert!(matches!(iterate.state_of(j), State::free));
                debug_assert!(p_basis.status_of(j) == NONBASIC_FIXED);
                dx[ju] = -iterate.x_at(j);
                scatter_column(ai, j, iterate.x_at(j), &mut dxbasic[..]);
                dependent_cols.push(j);
            }
        }
        // IPX quirk: C++ aliases rhs and lhs (`SolveDense(dxbasic, dxbasic, 'N')`);
        // BASICLU copies rhs before writing lhs, so cloning rhs is exact.
        let rhs = dxbasic.clone();
        p_basis.solve_dense(&rhs, &mut dxbasic, b'N');
        for p in 0..m {
            dx[p_basis.at(p) as usize] = dxbasic[p as usize];
        }
    }

    // Postprocess linearly dependent rows. If the slack variable to an equality
    // constraint cannot be made nonbasic, then the constraint is redundant. We
    // move its dual variable y[i] to zero without altering A'*y and make the
    // constraint "implied" (i.e. the slack variable is treated as free by the
    // interior point solver).
    if info.dependent_rows > 0 {
        for p in 0..m {
            let j = p_basis.at(p);
            let ju = j as usize;
            // The starting basis is constructed such that structural variables
            // with weight zero will always be nonbasic.
            if lb[ju] == ub[ju] {
                debug_assert!(j >= n);
            }
            if j >= n && lb[ju] == ub[ju] {
                debug_assert!(p_basis.status_of(j) == BASIC_FREE);
                debug_assert!(iterate.has_barrier_lb(j));
                debug_assert!(iterate.has_barrier_ub(j));
                dy[p as usize] = -iterate.y_at(j - n);
                dependent_rows.push(j - n);
            }
        }
        // IPX quirk: aliased SolveDense, rhs cloned as above.
        let rhs = dy.clone();
        p_basis.solve_dense(&rhs, &mut dy, b'T');
        for &i in &dependent_rows {
            dy[i as usize] = -iterate.y_at(i); // would be already in exact arithmetic
        }
    }

    // Now adjust the iterate.
    iterate.update(
        1.0,
        Some(&dx[..]),
        None,
        None,
        1.0,
        Some(&dy[..]),
        None,
        None,
    );

    for &j in &dependent_cols {
        debug_assert!(iterate.x_at(j) == 0.0);
        iterate.make_fixed_value(j, 0.0);
    }
    for &i in &dependent_rows {
        debug_assert!(iterate.y_at(i) == 0.0);
        iterate.make_implied_eq(n + i); // sets zl[n+i] = zu[n+i] = 0
        debug_assert!(iterate.zl_at(n + i) == 0.0);
        debug_assert!(iterate.zu_at(n + i) == 0.0);
    }
}

/// HiGHS: ipx/starting_basis.cc::StartingBasis (lines 127-180)
pub fn starting_basis(iterate: &mut Iterate, p_basis: &mut Basis, info: &mut Info) {
    // Rc clone: `iterate` is mutated below, so the model cannot be borrowed
    // out of it for the whole function.
    let model = iterate.model_rc();
    let m = model.rows();
    let n = model.cols();
    let lb = model.lb();
    let ub = model.ub();
    let mut colscale = Vector::new(n + m);
    info.errflag = 0;
    let timer = Timer::new(0.0);

    // Construct starting basis. The column weights for the crash procedure are
    // the interior point scaling factors from the current iterate, except that
    // fixed variables get weight zero.
    for j in 0..n + m {
        let ju = j as usize;
        colscale[ju] = iterate.scaling_factor(j);
        if lb[ju].is_infinite() && ub[ju].is_infinite() {
            debug_assert!(colscale[ju] == f64::INFINITY);
        } else {
            debug_assert!(colscale[ju].is_finite());
        }
        if lb[ju] == ub[ju] {
            colscale[ju] = 0.0;
        }
    }
    p_basis.construct_basis_from_weights(&colscale[..], info);
    if info.errflag != 0 {
        return;
    }

    // Change status of free variables either to BASIC_FREE or NONBASIC_FIXED.
    // Change status of fixed variables either to NONBASIC_FIXED or BASIC_FREE.
    for j in 0..n + m {
        let ju = j as usize;
        if colscale[ju] == 0.0 || colscale[ju].is_infinite() {
            if p_basis.is_basic(j) {
                p_basis.free_basic_variable(j);
            } else {
                p_basis.fix_nonbasic_variable(j);
            }
        }
    }

    // Variables with equal lower and upper bounds can be removed from the
    // interior point solve if their basic status is NONBASIC_FIXED. Setting
    // x[j] to its bound alters the residual in AI*x=b; this can be ignored
    // because the iterate has not been primal feasible if x[j] was
    // significantly away from its bound.
    for j in 0..n + m {
        let ju = j as usize;
        if lb[ju] == ub[ju] && p_basis.status_of(j) == NONBASIC_FIXED {
            iterate.make_fixed_value(j, lb[ju]);
        }
    }

    // Linearly dependent rows and columns require adjusting the iterate.
    postprocess_dependencies(iterate, p_basis, info);

    if cfg!(debug_assertions) {
        assert_consistency(iterate, p_basis);
    }
    info.time_starting_basis += timer.elapsed();
}
