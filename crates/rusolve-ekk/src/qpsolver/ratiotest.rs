//! HiGHS source: qpsolver/ratiotest.hpp, qpsolver/ratiotest.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `ratiotest_twopass` ends with `fmax(alpha, 0)`: translate with
//! the C `fmax` semantics (NaN-ignoring), not `cpp_max`.

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: ratiotest.hpp::RatiotestResult (lines 15-19)
#[derive(Clone, Copy, Debug, Default)]
pub struct RatiotestResult {
    pub alpha: f64,
    pub limitingconstraint: HighsInt,
    pub nowactiveatlower: bool,
}

/// HiGHS: ratiotest.cpp::step (lines 10-18)
fn step(x: f64, p: f64, l: f64, u: f64, t: f64) -> f64 {
    if p < -t && l > -f64::INFINITY {
        (l - x) / p
    } else if p > t && u < f64::INFINITY {
        (u - x) / p
    } else {
        f64::INFINITY
    }
}

/// HiGHS: ratiotest.cpp::ratiotest_textbook (lines 20-55)
fn ratiotest_textbook(
    rt: &Runtime,
    p: &QpVector,
    rowmove: &QpVector,
    instance: &Instance,
    alphastart: f64,
) -> RatiotestResult {
    let mut result = RatiotestResult {
        limitingconstraint: -1,
        alpha: alphastart,
        nowactiveatlower: false,
    };

    // check ratio towards variable bounds
    for j in 0..p.num_nz {
        let i = p.index[j as usize];
        let alpha_i = step(
            rt.primal.value[i as usize],
            p.value[i as usize],
            instance.var_lo[i as usize],
            instance.var_up[i as usize],
            rt.settings.ratiotest_t,
        );
        if alpha_i < result.alpha {
            result.alpha = alpha_i;
            result.limitingconstraint = instance.num_con + i;
            result.nowactiveatlower = p.value[i as usize] < 0.0;
        }
    }

    // check ratio towards constraint bounds
    for j in 0..rowmove.num_nz {
        let i = rowmove.index[j as usize];
        let alpha_i = step(
            rt.rowactivity.value[i as usize],
            rowmove.value[i as usize],
            instance.con_lo[i as usize],
            instance.con_up[i as usize],
            rt.settings.ratiotest_t,
        );
        if alpha_i < result.alpha {
            result.alpha = alpha_i;
            result.limitingconstraint = i;
            result.nowactiveatlower = rowmove.value[i as usize] < 0.0;
        }
    }

    result
}

/// HiGHS: ratiotest.cpp::ratiotest_twopass (lines 57-105)
fn ratiotest_twopass(
    runtime: &Runtime,
    p: &QpVector,
    rowmove: &QpVector,
    relaxed: &Instance,
    alphastart: f64,
) -> RatiotestResult {
    let res1 = ratiotest_textbook(runtime, p, rowmove, relaxed, alphastart);

    if res1.limitingconstraint == -1 {
        return res1;
    }

    let mut result = res1;

    let mut max_pivot: f64 = 0.0;
    if res1.limitingconstraint != -1 {
        if result.limitingconstraint < runtime.instance.num_con {
            max_pivot = rowmove.value[result.limitingconstraint as usize];
        } else {
            max_pivot = p.value[(result.limitingconstraint - runtime.instance.num_con) as usize];
        }
    }

    for i in 0..runtime.instance.num_con {
        let step_i = step(
            runtime.rowactivity.value[i as usize],
            rowmove.value[i as usize],
            runtime.instance.con_lo[i as usize],
            runtime.instance.con_up[i as usize],
            runtime.settings.ratiotest_t,
        );
        if rowmove.value[i as usize].abs() >= max_pivot.abs() && step_i <= res1.alpha {
            max_pivot = rowmove.value[i as usize];
            result.limitingconstraint = i;
            result.alpha = step_i;
            result.nowactiveatlower = rowmove.value[i as usize] < 0.0;
        }
    }

    for i in 0..runtime.instance.num_var {
        let step_i = step(
            runtime.primal.value[i as usize],
            p.value[i as usize],
            runtime.instance.var_lo[i as usize],
            runtime.instance.var_up[i as usize],
            runtime.settings.ratiotest_t,
        );
        if p.value[i as usize].abs() >= max_pivot.abs() && step_i <= res1.alpha {
            max_pivot = p.value[i as usize];
            result.limitingconstraint = runtime.instance.num_con + i;
            result.alpha = step_i;
            result.nowactiveatlower = p.value[i as usize] < 0.0;
        }
    }

    // PORT NOTE: C fmax
    result.alpha = result.alpha.max(0.0);
    result
}

/// HiGHS: ratiotest.cpp::ratiotest_relax_instance (lines 107-133)
pub fn ratiotest_relax_instance(runtime: &Runtime) -> Instance {
    let mut relaxed_instance = runtime.instance.clone();
    for bound in relaxed_instance.con_lo.iter_mut() {
        if *bound != -f64::INFINITY {
            *bound -= runtime.settings.ratiotest_d;
        }
    }

    for bound in relaxed_instance.con_up.iter_mut() {
        if *bound != f64::INFINITY {
            *bound += runtime.settings.ratiotest_d;
        }
    }

    for bound in relaxed_instance.var_lo.iter_mut() {
        if *bound != -f64::INFINITY {
            *bound -= runtime.settings.ratiotest_d;
        }
    }

    for bound in relaxed_instance.var_up.iter_mut() {
        if *bound != f64::INFINITY {
            *bound += runtime.settings.ratiotest_d;
        }
    }
    relaxed_instance
}

/// HiGHS: ratiotest.cpp::ratiotest (lines 135-146)
pub fn ratiotest(
    runtime: &Runtime,
    p: &QpVector,
    rowmove: &QpVector,
    alphastart: f64,
) -> RatiotestResult {
    match runtime.settings.ratiotest {
        RatiotestStrategy::Textbook => {
            ratiotest_textbook(runtime, p, rowmove, &runtime.instance, alphastart)
        }
        RatiotestStrategy::TwoPass => ratiotest_twopass(
            runtime,
            p,
            rowmove,
            &runtime.relaxed_for_ratiotest,
            alphastart,
        ),
    }
}
