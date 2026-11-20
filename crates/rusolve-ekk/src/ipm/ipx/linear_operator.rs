//! HiGHS source: ipm/ipx/linear_operator.h, linear_operator.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::*;

/// HiGHS: ipx/linear_operator.h::LinearOperator (abstract)
///
/// Implementors: `NormalMatrix`, `SplittedNormalMatrix`, `DiagonalPrecond`.
/// `rhs_dot_lhs` may be `nullptr`.
pub trait LinearOperator {
    /// HiGHS: ipx/linear_operator.cc::LinearOperator::Apply (lines 5-8)
    fn apply(&mut self, rhs: &Vector, lhs: &mut Vector, rhs_dot_lhs: Option<&mut f64>) {
        self._apply(rhs, lhs, rhs_dot_lhs);
    }

    /// HiGHS: ipx/linear_operator.h::LinearOperator::_Apply (pure virtual)
    fn _apply(&mut self, rhs: &Vector, lhs: &mut Vector, rhs_dot_lhs: Option<&mut f64>);
}
