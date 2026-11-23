//! HiGHS source: ipm/ipx/kkt_solver.h, kkt_solver.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::basis::Basis;
use super::iterate::Iterate;
use super::timer::Timer;
use super::*;

/// HiGHS: ipx/kkt_solver.h::KKTSolver (abstract)
///
/// Implementors: `KKTSolverDiag` (IPM before the switch) and
/// `KKTSolverBasis<'a>` (after). `Factorize(nullptr, info)` is called by
/// `IPM::ComputeStartingPoint` (ipm.cc:207), hence `Option<&mut Iterate>`.
pub trait KKTSolver {
    /// HiGHS: ipx/kkt_solver.cc::KKTSolver::Factorize (lines 6-10)
    fn factorize(&mut self, iterate: Option<&mut Iterate>, info: &mut Info) {
        let timer = Timer::default();
        self._factorize(iterate, info);
        info.time_kkt_factorize += timer.elapsed();
    }
    /// HiGHS: ipx/kkt_solver.cc::KKTSolver::Solve (lines 12-17)
    fn solve(
        &mut self,
        a: &Vector,
        b: &Vector,
        tol: f64,
        x: &mut Vector,
        y: &mut Vector,
        info: &mut Info,
    ) {
        let timer = Timer::default();
        self._solve(a, b, tol, x, y, info);
        info.time_kkt_solve += timer.elapsed();
    }
    /// HiGHS: ipx/kkt_solver.cc::KKTSolver::iter (line 19)
    fn iter(&self) -> Int {
        self._iter()
    }
    /// HiGHS: ipx/kkt_solver.cc::KKTSolver::basis_changes (line 20)
    fn basis_changes(&self) -> Int {
        self._basis_changes()
    }
    /// HiGHS: ipx/kkt_solver.cc::KKTSolver::basis (line 21)
    fn basis(&self) -> Option<&Basis> {
        self._basis()
    }

    // private virtuals (kkt_solver.h)
    fn _factorize(&mut self, iterate: Option<&mut Iterate>, info: &mut Info);
    fn _solve(
        &mut self,
        a: &Vector,
        b: &Vector,
        tol: f64,
        x: &mut Vector,
        y: &mut Vector,
        info: &mut Info,
    );
    fn _iter(&self) -> Int;
    /// HiGHS: ipx/kkt_solver.h::KKTSolver::_basis_changes (default 0)
    fn _basis_changes(&self) -> Int {
        0
    }
    /// HiGHS: ipx/kkt_solver.h::KKTSolver::_basis (default nullptr)
    fn _basis(&self) -> Option<&Basis> {
        None
    }
}
