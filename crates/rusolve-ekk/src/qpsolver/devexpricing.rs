//! HiGHS source: qpsolver/devexpricing.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! See `pricing.rs` for the ownership note.

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: devexpricing.hpp::DevexPricing (lines 18-104)
#[derive(Clone, Debug, Default)]
pub struct DevexPricing {
    weights: Vec<f64>,
}

impl DevexPricing {
    /// HiGHS: devexpricing.hpp::DevexPricing::chooseconstrainttodrop (lines 26-59)
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
            let val = lambda.value[indexinbasis as usize] * lambda.value[indexinbasis as usize]
                / self.weights[indexinbasis as usize];
            if val > maxabslambda
                && lambda.value[indexinbasis as usize].abs() > rt.settings.lambda_zero_threshold
            {
                if basis.getstatus(iVar) == BasisStatus::kActiveAtLower
                    && -lambda.value[indexinbasis as usize] > 0.0
                {
                    minidx = iVar;
                    maxabslambda = val;
                } else if basis.getstatus(iVar) == BasisStatus::kActiveAtUpper
                    && lambda.value[indexinbasis as usize] > 0.0
                {
                    minidx = iVar;
                    maxabslambda = val;
                } else {
                    // C++ lines 52-54: empty else branch (TODO), dropped.
                }
            }
        }

        minidx
    }

    /// HiGHS: devexpricing.hpp::DevexPricing::DevexPricing (lines 62-67)
    pub fn new(rt: &Runtime) -> Self {
        DevexPricing {
            weights: vec![1.0; rt.instance.num_var as usize],
        }
    }
}

impl Pricing for DevexPricing {
    /// HiGHS: devexpricing.hpp::DevexPricing::price (lines 77-81)
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

    /// HiGHS: devexpricing.hpp::DevexPricing::recompute (lines 83-85)
    fn recompute(&mut self, _rt: &Runtime, _basis: &mut Basis) {
        // C++ lines 83-85: do nothing.
    }

    /// HiGHS: devexpricing.hpp::DevexPricing::update_weights (lines 87-103)
    fn update_weights(
        &mut self,
        rt: &Runtime,
        basis: &mut Basis,
        aq: &QpVector,
        _ep: &QpVector,
        p: HighsInt,
        _q: HighsInt,
    ) {
        let rowindex_p = basis.constraintindexinbasisfactor[p as usize];
        let weight_p = self.weights[rowindex_p as usize];
        for i in 0..rt.instance.num_var {
            if i == rowindex_p {
                self.weights[i as usize] =
                    weight_p / (aq.value[rowindex_p as usize] * aq.value[rowindex_p as usize]);
            } else {
                self.weights[i as usize] += (aq.value[i as usize] * aq.value[i as usize])
                    / (aq.value[rowindex_p as usize] * aq.value[rowindex_p as usize])
                    * weight_p
                    * weight_p;
            }
            if self.weights[i as usize] > 1e7 {
                self.weights[i as usize] = 1.0;
            }
        }
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

    /// min 0.5 (2 x0^2 + 4 x1^2) - 2 x0 - 4 x1, 0 <= x <= 10, x0 + x1 <= 10.
    fn box_qp() -> Instance {
        let mut inst = Instance::new(2, 1);
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
        inst.A.mat = MatrixBase {
            num_row: 1,
            num_col: 2,
            start: vec![0, 1, 2],
            index: vec![0, 0],
            value: vec![1.0, 1.0],
        };
        inst.con_lo = vec![-kHighsInf];
        inst.con_up = vec![10.0];
        inst.var_lo = vec![0.0, 0.0];
        inst.var_up = vec![10.0, 10.0];
        inst
    }

    /// Both variable bounds (ids 1, 2) active at lower: B = I.
    fn at_lower(rt: &mut Runtime) -> Basis {
        // basis.cpp lines 27-30: `status[i]` belongs to `active[i]`.
        let status = vec![BasisStatus::kActiveAtLower, BasisStatus::kActiveAtLower];
        Basis::new(rt, vec![1, 2], status, vec![])
    }

    fn dense(v: &[f64]) -> QpVector {
        QpVector {
            num_nz: v.len() as HighsInt,
            dim: v.len() as HighsInt,
            index: (0..v.len() as HighsInt).collect(),
            value: v.to_vec(),
        }
    }

    #[test]
    fn unit_weights_drop_most_negative_multiplier() {
        // Weights start at 1, so val = lambda^2 and lambda = (-2, -4) at x = 0:
        // constraint id 2 (x1) is chosen, as with Dantzig.
        let mut rt = Runtime::new(&box_qp(), Statistics::default());
        let mut basis = at_lower(&mut rt);
        let mut gradient = Gradient::new(&rt);
        let mut redcosts = ReducedCosts::new(&rt);
        let mut pricing = DevexPricing::new(&rt);
        assert_eq!(pricing.weights, vec![1.0, 1.0]);
        let x = rt.primal.clone();
        let minidx = pricing.price(&rt, &mut basis, &mut redcosts, &mut gradient, &x);
        assert_eq!(minidx, 2);
    }

    #[test]
    fn update_weights_then_price_prefers_light_row() {
        // C++ lines 89-102 with p = 1 (x0's bound), aq[r] = 2, aq[other] = 3,
        // weight_p = 1: weights[r] = 1 / 4, weights[other] = 1 + 9 / 4.
        let mut rt = Runtime::new(&box_qp(), Statistics::default());
        let mut basis = at_lower(&mut rt);
        let mut gradient = Gradient::new(&rt);
        let mut redcosts = ReducedCosts::new(&rt);
        let mut pricing = DevexPricing::new(&rt);
        let r = basis.getindexinfactor()[1] as usize;
        let other = 1 - r;
        let mut aq = vec![0.0; 2];
        aq[r] = 2.0;
        aq[other] = 3.0;
        let ep = QpVector::new(2);
        pricing.update_weights(&rt, &mut basis, &dense(&aq), &ep, 1, 0);
        assert_eq!(pricing.weights[r], 0.25);
        assert_eq!(pricing.weights[other], 3.25);
        // lambda_x0 = -2 now scores 4 / 0.25 = 16 against 16 / 3.25.
        let x = rt.primal.clone();
        let minidx = pricing.price(&rt, &mut basis, &mut redcosts, &mut gradient, &x);
        assert_eq!(minidx, 1);
    }

    #[test]
    fn update_weights_resets_huge_weights() {
        // C++ lines 99-101: any weight above 1e7 is reset to 1.
        let mut rt = Runtime::new(&box_qp(), Statistics::default());
        let mut basis = at_lower(&mut rt);
        let mut pricing = DevexPricing::new(&rt);
        let r = basis.getindexinfactor()[1] as usize;
        let mut aq = vec![1.0; 2];
        aq[r] = 1e-4;
        let ep = QpVector::new(2);
        pricing.update_weights(&rt, &mut basis, &dense(&aq), &ep, 1, 0);
        assert_eq!(pricing.weights, vec![1.0, 1.0]);
    }
}
