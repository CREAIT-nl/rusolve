//! HiGHS source: qpsolver/dantzigpricing.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! See `pricing.rs` for the ownership note.

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: dantzigpricing.hpp::DantzigPricing (lines 18-76)
#[derive(Clone, Debug, Default)]
pub struct DantzigPricing {}

impl DantzigPricing {
    /// HiGHS: dantzigpricing.hpp::DantzigPricing::chooseconstrainttodrop (lines 24-56)
    fn chooseconstrainttodrop(
        &self,
        rt: &Runtime,
        basis: &mut Basis,
        lambda: &QpVector,
    ) -> HighsInt {
        // C++ `auto` copies both vectors; the copies also free `basis` for
        // `getstatus`, which takes `&mut self` (std::map operator[] inserts).
        let active_constraint_index = basis.getactive().clone();
        let constraintindexinbasisfactor = basis.getindexinfactor().clone();

        let mut minidx: HighsInt = -1;
        let mut maxabslambda: f64 = 0.0;
        for i in 0..active_constraint_index.len() {
            let iVar = active_constraint_index[i];
            let indexinbasis = constraintindexinbasisfactor[iVar as usize];
            if indexinbasis == -1 {
                println!("error");
            }
            debug_assert!(indexinbasis != -1);
            if rt.instance.is_equality(iVar) {
                continue;
            }
            if basis.getstatus(iVar) == BasisStatus::kActiveAtLower
                && -lambda.value[indexinbasis as usize] > maxabslambda
            {
                minidx = iVar;
                maxabslambda = -lambda.value[indexinbasis as usize];
            } else if basis.getstatus(iVar) == BasisStatus::kActiveAtUpper
                && lambda.value[indexinbasis as usize] > maxabslambda
            {
                minidx = iVar;
                maxabslambda = lambda.value[indexinbasis as usize];
            } else {
                // C++ lines 46-48: empty else branch (TODO), dropped.
            }
        }

        if maxabslambda <= rt.settings.lambda_zero_threshold {
            // C++ line 52: commented-out printf, dropped.
            return -1;
        }
        minidx
    }

    /// HiGHS: dantzigpricing.hpp::DantzigPricing::DantzigPricing (lines 59-61)
    pub fn new(_rt: &Runtime) -> Self {
        DantzigPricing {}
    }
}

impl Pricing for DantzigPricing {
    /// HiGHS: dantzigpricing.hpp::DantzigPricing::price (lines 63-66)
    fn price(
        &mut self,
        rt: &Runtime,
        basis: &mut Basis,
        redcosts: &mut ReducedCosts,
        gradient: &mut Gradient,
        _x: &QpVector,
    ) -> HighsInt {
        // C++ evaluates the reduced-costs argument before the call; the
        // borrow is split so `basis` is not borrowed twice.
        let lambda = redcosts.get_reduced_costs(rt, basis, gradient);
        self.chooseconstrainttodrop(rt, basis, lambda)
    }

    /// HiGHS: dantzigpricing.hpp::DantzigPricing::recompute (lines 68-70)
    fn recompute(&mut self, _rt: &Runtime, _basis: &mut Basis) {
        // C++ lines 68-70: do nothing.
    }

    /// HiGHS: dantzigpricing.hpp::DantzigPricing::update_weights (lines 72-75)
    fn update_weights(
        &mut self,
        _rt: &Runtime,
        _basis: &mut Basis,
        _aq: &QpVector,
        _ep: &QpVector,
        _p: HighsInt,
        _q: HighsInt,
    ) {
        // C++ lines 72-75: does nothing.
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use crate::qpsolver::*;

    /// min 0.5 (2 x0^2 + 4 x1^2) - 2 x0 - 4 x1, 0 <= x <= 10 and, when
    /// `with_row`, x0 + x1 <= 10. Optimum x = (1, 1), objective -3.
    fn box_qp(with_row: bool) -> Instance {
        let nc = if with_row { 1 } else { 0 };
        let mut inst = Instance::new(2, nc);
        inst.c = QpVector {
            num_nz: 2,
            dim: 2,
            index: vec![0, 1],
            value: vec![-2.0, -4.0],
        };
        inst.Q.mat = MatrixBase {
            num_row: 2,
            num_col: 2,
            start: vec![0, 1, 2],
            index: vec![0, 1],
            value: vec![2.0, 4.0],
        };
        inst.A.mat = if with_row {
            MatrixBase {
                num_row: 1,
                num_col: 2,
                start: vec![0, 1, 2],
                index: vec![0, 0],
                value: vec![1.0, 1.0],
            }
        } else {
            MatrixBase {
                num_row: 0,
                num_col: 2,
                start: vec![0, 0, 0],
                index: vec![],
                value: vec![],
            }
        };
        inst.con_lo = vec![-kHighsInf; nc as usize];
        inst.con_up = vec![10.0; nc as usize];
        inst.var_lo = vec![0.0, 0.0];
        inst.var_up = vec![10.0, 10.0];
        inst
    }

    /// Both variables active at their lower bound (constraint ids
    /// `num_con + i`); the row (id 0) is not in the basis.
    fn at_lower(rt: &mut Runtime) -> Basis {
        // basis.cpp lines 27-30: `status[i]` belongs to `active[i]`.
        let status = vec![BasisStatus::kActiveAtLower, BasisStatus::kActiveAtLower];
        Basis::new(rt, vec![1, 2], status, vec![])
    }

    /// Both variables in the nullspace (inactive in basis): Z = I.
    fn all_free(rt: &mut Runtime) -> Basis {
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kInactiveInBasis,
            BasisStatus::kInactiveInBasis,
        ];
        Basis::new(rt, vec![], status, vec![1, 2])
    }

    fn close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-9)
    }

    #[test]
    fn drops_most_negative_multiplier() {
        // At x = 0 with both bounds active the multipliers are the gradient
        // (-2, -4); Dantzig drops the most violated: constraint id 2 (x1).
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut basis = at_lower(&mut rt);
        let mut gradient = Gradient::new(&rt);
        let mut redcosts = ReducedCosts::new(&rt);
        let mut pricing = DantzigPricing::new(&rt);
        let x = rt.primal.clone();
        let minidx = pricing.price(&rt, &mut basis, &mut redcosts, &mut gradient, &x);
        assert_eq!(minidx, 2);
    }
}
