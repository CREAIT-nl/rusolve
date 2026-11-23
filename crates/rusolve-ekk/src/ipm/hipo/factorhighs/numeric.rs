//! HiGHS source: ipm/hipo/factorhighs/Numeric.h, ipm/hipo/factorhighs/Numeric.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The numeric factorisation. The C++ keeps pointers to the `Symbolic`, the
//! `sn_columns` and the `DataCollector` that `Factorise::run` set; the Rust
//! passes them to `solve` / `getReg` and keeps `factorised_` for the
//! `!sn_columns_ || !S_` test.

use super::super::Int;
use super::super::auxiliary::{permuteVector, permuteVectorInverse};
use super::data_collector::DataCollector;
use super::hybrid_solve_handler::HybridSolveHandler;
use super::solve_handler::SolveHandlerOps;
use super::symbolic::Symbolic;
use super::*;

/// HiGHS: Numeric.h `class Numeric`
#[derive(Clone, Debug, Default)]
pub struct Numeric {
    /// Rust only: true once `Factorise::run` has set the C++ pointers.
    pub factorised_: bool,
    pub swaps_: Vec<Vec<Int>>,
    pub pivot_2x2_: Vec<Vec<f64>>,
    pub total_reg_: Vec<f64>,
}

impl Numeric {
    /// HiGHS: Numeric.cpp::Numeric::solve (lines 15-44)
    /// Returns the number of solves performed (or `kRetInvalidPointer`).
    pub fn solve(
        &self,
        S: &Symbolic,
        sn_columns: &[Vec<f64>],
        data: &DataCollector,
        x: &mut Vec<f64>,
    ) -> Int {
        if !self.factorised_ {
            return kRetInvalidPointer;
        }

        // HIPO_CLOCK_CREATE; (line 20) — level-0 no-op, dropped.

        // initialise solve handler
        let SH = HybridSolveHandler::new(S, sn_columns, &self.swaps_, &self.pivot_2x2_, data);

        // permute rhs
        // HIPO_CLOCK_START(2); (line 26) — level-0 no-op, dropped.
        permuteVectorInverse(x, S.iperm());
        // HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare); (line 28) — dropped.

        // solve
        // HIPO_CLOCK_START(2); (line 31) — level-0 no-op, dropped.
        SH.forwardSolve(x);
        SH.diagSolve(x);
        SH.backwardSolve(x);
        // HIPO_CLOCK_STOP(2, *data_, kTimeSolveSolve); (line 35) — dropped.

        // unpermute solution
        // HIPO_CLOCK_START(2); (line 38) — level-0 no-op, dropped.
        permuteVector(x, S.iperm());
        // HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare); (line 40) — dropped.

        // HIPO_CLOCK_STOP(1, *data_, kTimeSolve); (line 42) — dropped.
        kRetOk
    }

    /// HiGHS: Numeric.cpp::Numeric::getReg (lines 46-51)
    pub fn getReg(&mut self, S: &Symbolic, reg: &mut Vec<f64>) {
        // unpermute regularisation
        permuteVector(&mut self.total_reg_, S.iperm());

        *reg = std::mem::take(&mut self.total_reg_);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn solve_without_factor_is_invalid_pointer() {
        let n = Numeric::default();
        let mut x = vec![1.0];
        assert_eq!(
            n.solve(&Symbolic::default(), &[], &DataCollector::new(), &mut x),
            kRetInvalidPointer
        );
    }
}
