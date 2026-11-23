//! HiGHS source: ipm/hipo/ipm/LinearSolver.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The `LinearSolver` interface (complete: it has no out-of-line code).
//! `valid_` is reached through `valid()` / `set_valid()`. The references a
//! C++ implementation caches at construction (`FactorHiGHSSolver::kkt_`,
//! `options_`, `model_`, `regul_`, `info_`, `data_`, `logger_`) are bundled
//! into [`LsCtx`] and passed to every method that reads them; build one
//! with `Solver::split_ls`.

use super::info::Info;
use super::ipm_data::IpmData;
use super::kkt_matrix::KktMatrix;
use super::model::Model;
use super::options::Options;
use super::parameters::Regularisation;
use crate::ipm::hipo::Int;
use crate::ipm::hipo::auxiliary::logger::Logger;

/// The C++ references of a `LinearSolver` implementation, per call.
pub struct LsCtx<'c> {
    pub kkt: &'c mut KktMatrix,
    pub options: &'c mut Options,
    pub model: &'c Model,
    pub regul: &'c Regularisation,
    pub info: &'c mut Info,
    pub data: &'c mut IpmData,
    pub logger: &'c Logger,
}

/// HiGHS: LinearSolver.h `class LinearSolver`
pub trait LinearSolver {
    /// HiGHS: LinearSolver.h `bool valid_ = false` (read)
    fn valid(&self) -> bool;
    /// HiGHS: LinearSolver.h `bool valid_` (write)
    fn set_valid(&mut self, valid: bool);

    /// HiGHS: LinearSolver.h::factorAS (pure virtual)
    fn factorAS(&mut self, ctx: &mut LsCtx<'_>, scaling: &[f64]) -> Int;
    /// HiGHS: LinearSolver.h::solveAS (pure virtual)
    fn solveAS(
        &mut self,
        ctx: &mut LsCtx<'_>,
        rhs_x: &[f64],
        rhs_y: &[f64],
        lhs_x: &mut Vec<f64>,
        lhs_y: &mut Vec<f64>,
    ) -> Int;
    /// HiGHS: LinearSolver.h::factorNE (pure virtual)
    fn factorNE(&mut self, ctx: &mut LsCtx<'_>, scaling: &[f64]) -> Int;
    /// HiGHS: LinearSolver.h::solveNE (pure virtual)
    fn solveNE(&mut self, ctx: &mut LsCtx<'_>, rhs: &[f64], lhs: &mut Vec<f64>) -> Int;
    /// HiGHS: LinearSolver.h::clear (pure virtual)
    fn clear(&mut self);

    /// HiGHS: LinearSolver.h::setup (virtual, `{ return 0; }`)
    fn setup(&mut self, ctx: &mut LsCtx<'_>) -> Int {
        0
    }
    /// HiGHS: LinearSolver.h::flops (virtual, `{ return 0; }`)
    fn flops(&self) -> f64 {
        0.0
    }
    /// HiGHS: LinearSolver.h::spops (virtual, `{ return 0; }`)
    fn spops(&self) -> f64 {
        0.0
    }
    /// HiGHS: LinearSolver.h::nz (virtual, `{ return 0; }`)
    fn nz(&self) -> f64 {
        0.0
    }
    /// HiGHS: LinearSolver.h::getReg (virtual, `{}`)
    fn getReg(&mut self, reg: &mut Vec<f64>) {}
}
