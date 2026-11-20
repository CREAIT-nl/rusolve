//! HiGHS source: ipm/hipo/factorhighs/SolveHandler.h, ipm/hipo/factorhighs/SolveHandler.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Base of the triangular-solve handlers (virtual-dispatch mapping in
//! `hipo/mod.rs`).

use super::data_collector::DataCollector;
use super::symbolic::Symbolic;

/// HiGHS: SolveHandler.h `class SolveHandler` (data members)
pub struct SolveHandler<'a> {
    pub S_: &'a Symbolic,
    pub sn_columns_: &'a [Vec<f64>],
    pub data_: &'a DataCollector,
}

impl<'a> SolveHandler<'a> {
    /// HiGHS: SolveHandler.cpp::SolveHandler::SolveHandler (lines 5-8)
    pub fn new(S: &'a Symbolic, sn_columns: &'a [Vec<f64>], data: &'a DataCollector) -> Self {
        SolveHandler {
            S_: S,
            sn_columns_: sn_columns,
            data_: data,
        }
    }
}

/// HiGHS: SolveHandler.h `class SolveHandler` (virtual methods)
pub trait SolveHandlerOps<'a> {
    fn base(&self) -> &SolveHandler<'a>;
    fn forwardSolve(&self, x: &mut Vec<f64>);
    fn backwardSolve(&self, x: &mut Vec<f64>);
    fn diagSolve(&self, x: &mut Vec<f64>);
}
