//! HiGHS source: ipm/ipx/ipm.h, ipm.cc (struct Step, class IPM)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `ipm.cc` (895 lines) is split: [`part1`] lines 11-448 (Step, ctor,
//! StartingPoint, Driver, ComputeStartingPoint, StepToBoundary, Predictor,
//! AddCorrector), [`part2`] lines 450-895 (centring, step sizes, Newton
//! system, output).
//!
//! C++ `IPM::Driver`/`StartingPoint` store `kkt_`, `iterate_`, `info_` as
//! members for the duration of the call. Rust drops those members: the
//! public entry points build an [`IpmRefs`] from their arguments and pass it
//! to every private method (`iterate_->x()` -> `r.iterate.x()`).

use super::control::Control;
use super::info::*;
use super::iterate::Iterate;
use super::kkt_solver::KKTSolver;
use super::*;
use std::rc::Rc;

mod part1;
mod part2;

/// The C++ pointer members `kkt_`, `iterate_`, `info_` of `IPM`.
pub struct IpmRefs<'a> {
    pub kkt: &'a mut dyn KKTSolver,
    pub iterate: &'a mut Iterate,
    pub info: &'a mut Info,
}

/// HiGHS: ipx/ipm.cc::IPM::Step (lines 11-19)
#[derive(Debug, Clone)]
pub struct Step {
    pub x: Vector,
    pub xl: Vector,
    pub xu: Vector,
    pub y: Vector,
    pub zl: Vector,
    pub zu: Vector,
}

/// HiGHS: ipx/ipm.h::IPM
pub struct IPM {
    control_: Rc<Control>,
    step_primal_: f64,
    step_dual_: f64,
    num_bad_iter_: Int,
    best_complementarity_: f64,
    maxiter_: Int,
    centring_ratio: f64,
    bad_products: Int,
}

impl IPM {
    pub const kDivergeTol: f64 = 1e6;

    pub fn maxiter(&self) -> Int {
        self.maxiter_
    }
    /// HiGHS: ipx/ipm.h::IPM::maxiter(Int)
    pub fn set_maxiter(&mut self, i: Int) {
        self.maxiter_ = i;
    }
}
