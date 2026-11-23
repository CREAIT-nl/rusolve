//! HiGHS source: ipm/basiclu/basiclu_factorize.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU source)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_build_factors::lu_build_factors;
use super::lu_condest::lu_condest;
use super::lu_factorize_bump::lu_factorize_bump;
use super::lu_internal::{lu_load, lu_reset, lu_save};
use super::lu_residual_test::lu_residual_test;
use super::lu_setup_bump::lu_setup_bump;
use super::lu_singletons::lu_singletons;
use super::*;

/// HiGHS: basiclu/basiclu_factorize.c::basiclu_factorize (lines 10-132)
///
/// PORT NOTE: the C `NULL` checks of `Li..Bx` (`BASICLU_ERROR_argument_missing`)
/// are dropped: slices cannot be null and IPX never passes `NULL`. The
/// `goto` resume points become a `match this.task` fall-through ladder.
pub fn basiclu_factorize(
    istore: &mut [LuInt],
    xstore: &mut [f64],
    li: &mut [LuInt],
    lx: &mut [f64],
    ui: &mut [LuInt],
    ux: &mut [f64],
    wi: &mut [LuInt],
    wx: &mut [f64],
    bbegin: &[LuInt],
    bend: &[LuInt],
    bi: &[LuInt],
    bx: &[f64],
    c0ntinue: LuInt,
) -> LuInt {
    let mut this = match lu_load(istore, xstore, li, lx, ui, ux, wi, wx) {
        Ok(this) => this,
        Err(status) => return status,
    };

    if c0ntinue == 0 {
        lu_reset(&mut this);
        this.task = SINGLETONS;
    }

    // continue factorization: entry point of the goto ladder
    // (0 = singletons, 1 = setup_bump, 2 = factorize_bump, 3 = build_factors)
    let entry = match this.task {
        SINGLETONS => 0,
        SETUP_BUMP => 1,
        FACTORIZE_BUMP => 2,
        BUILD_FACTORS => 3,
        _ => return lu_save(&mut this, BASICLU_ERROR_invalid_call),
    };

    // Each of the following four parts of the factorization calls a routine
    // lu_do_something() which may request reallocation. In this case return
    // to the caller immediately, keeping the entry point in this.task.
    let mut status = BASICLU_OK;

    if entry <= 0 {
        // singletons:
        this.task = SINGLETONS;
        status = lu_singletons(&mut this, bbegin, bend, bi, bx);
        if status != BASICLU_OK {
            return lu_save(&mut this, status);
        }
    }

    if entry <= 1 {
        // setup_bump:
        this.task = SETUP_BUMP;
        status = lu_setup_bump(&mut this, bbegin, bend, bi, bx);
        if status != BASICLU_OK {
            return lu_save(&mut this, status);
        }
    }

    if entry <= 2 {
        // factorize_bump:
        this.task = FACTORIZE_BUMP;
        status = lu_factorize_bump(&mut this);
        if status != BASICLU_OK {
            return lu_save(&mut this, status);
        }
    }

    if entry <= 3 {
        // build_factors:
        this.task = BUILD_FACTORS;
        status = lu_build_factors(&mut this);
        if status != BASICLU_OK {
            return lu_save(&mut this, status);
        }
    }

    // factorization successfully finished
    this.task = NO_TASK;
    this.nupdate = 0; // make factorization valid
    this.ftran_for_update = -1;
    this.btran_for_update = -1;
    this.nfactorize += 1;

    let m = this.m;
    let mu = m as usize;
    // Lbegin = Wbegin + m+1, p = Wblink + m+1
    this.condest_l = lu_condest(
        m,
        &this.wbegin[mu + 1..],
        &this.lindex[..],
        &this.lvalue[..],
        None,
        &this.wblink[mu + 1..],
        0,
        &mut this.work1[..],
        Some(&mut this.norm_l),
        Some(&mut this.normest_linv),
    );
    this.condest_u = lu_condest(
        m,
        &this.ubegin[..],
        &this.uindex[..],
        &this.uvalue[..],
        Some(&this.row_pivot[..]),
        &this.wblink[mu + 1..],
        1,
        &mut this.work1[..],
        Some(&mut this.norm_u),
        Some(&mut this.normest_uinv),
    );

    // measure numerical stability of the factorization
    lu_residual_test(&mut this, bbegin, bend, bi, bx);

    // factor_cost is a deterministic measure of the factorization cost.
    // The parameters have been adjusted such that (on my computer)
    // 1e-6 * factor_cost =~ time_factorize.
    //
    // update_cost measures the accumulated cost of updates/solves compared
    // to the last factorization. It is computed from
    //
    //   update_cost = update_cost_numer / update_cost_denom.
    //
    // update_cost_denom is fixed here.
    // update_cost_numer is zero here and increased by solves/updates.
    let factor_cost = 0.04 * this.m as f64
        + 0.07 * this.matrix_nz as f64
        + 0.20 * this.bump_nz as f64
        + 0.20 * this.nsearch_pivot as f64
        + 0.008 * this.factor_flops as f64;

    this.update_cost_denom = factor_cost * 250.0;

    if this.rank < this.m {
        status = BASICLU_WARNING_singular_matrix;
    }

    // return_to_caller:
    lu_save(&mut this, status)
}
