//! HiGHS source: qpsolver/steepestedgepricing.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! See `pricing.rs` for the ownership note. The constructor calls
//! `compute_exact_weights`, which BTRANs through the basis, so `new` takes it.

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing (lines 17-169)
#[derive(Clone, Debug, Default)]
pub struct SteepestEdgePricing {
    weights: Vec<f64>,
}

impl SteepestEdgePricing {
    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::chooseconstrainttodrop (lines 24-57)
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
        let mut maxval: f64 = 0.0;
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
            if val > maxval
                && lambda.value[indexinbasis as usize].abs() > rt.settings.lambda_zero_threshold
            {
                if basis.getstatus(iVar) == BasisStatus::kActiveAtLower
                    && -lambda.value[indexinbasis as usize] > 0.0
                {
                    minidx = iVar;
                    maxval = val;
                } else if basis.getstatus(iVar) == BasisStatus::kActiveAtUpper
                    && lambda.value[indexinbasis as usize] > 0.0
                {
                    minidx = iVar;
                    maxval = val;
                } else {
                    // C++ lines 50-52: empty else branch (TODO), dropped.
                }
            }
        }

        minidx
    }

    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::SteepestEdgePricing (lines 60-66)
    pub fn new(rt: &Runtime, bas: &mut Basis) -> Self {
        let mut sp = SteepestEdgePricing {
            weights: vec![1.0; rt.instance.num_var as usize],
        };
        sp.compute_exact_weights(rt, bas);
        sp
    }

    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::compute_exact_weight (lines 73-76)
    pub fn compute_exact_weight(&self, rt: &Runtime, basis: &mut Basis, i: HighsInt) -> f64 {
        let unit = QpVector::unit(rt.instance.num_var, i);
        let y_i = basis.btran(&unit, false, -1);
        y_i.dot(&y_i)
    }

    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::compute_exact_weights (lines 78-82)
    pub fn compute_exact_weights(&mut self, rt: &Runtime, basis: &mut Basis) {
        for i in 0..rt.instance.num_var {
            self.weights[i as usize] = self.compute_exact_weight(rt, basis, i);
        }
    }

    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::check_weight (lines 84-93)
    pub fn check_weight(&self, rt: &Runtime, basis: &mut Basis, i: HighsInt) -> bool {
        let weight_is = self.weights[i as usize];
        let weight_comp = self.compute_exact_weight(rt, basis, i);
        if (weight_comp - weight_is).abs() > 1e-2 {
            // C++ lines 88-89: commented-out printf, dropped.
            return false;
        }
        true
    }

    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::check_weights (lines 95-120)
    pub fn check_weights(&self, rt: &Runtime, basis: &mut Basis) -> bool {
        let mut correct_weights: Vec<i32> = Vec::new();
        let mut incorrect_weights: Vec<i32> = Vec::new();
        for i in 0..rt.instance.num_var {
            let correct = self.check_weight(rt, basis, i);
            if correct {
                correct_weights.push(i);
            } else {
                incorrect_weights.push(i);
            }
        }

        print!("correct weights: ");
        for i in &correct_weights {
            print!("{} ", i);
        }
        println!();

        print!("incorrect weights: ");
        for i in &incorrect_weights {
            print!("{} ", i);
        }
        println!();

        incorrect_weights.len() == 0
    }
}

impl Pricing for SteepestEdgePricing {
    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::price (lines 68-71)
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

    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::recompute (line 122)
    fn recompute(&mut self, rt: &Runtime, basis: &mut Basis) {
        self.compute_exact_weights(rt, basis);
    }

    /// HiGHS: steepestedgepricing.hpp::SteepestEdgePricing::update_weights (lines 124-168)
    fn update_weights(
        &mut self,
        rt: &Runtime,
        basis: &mut Basis,
        aq: &QpVector,
        ep: &QpVector,
        p: HighsInt,
        _q: HighsInt,
    ) {
        let rowindex_p = basis.constraintindexinbasisfactor[p as usize];
        // C++ lines 127-128: commented-out printf, dropped.

        // C++ lines 130-133: commented-out check_weights block, dropped.

        let delta = basis.ftran(aq, false, -1);

        // C++ line 137: commented-out old_weight_p_updated, dropped.
        // C++ line 138: comment, dropped.
        let old_weight_p_computed = ep.dot(ep);

        // C++ lines 141-144: commented-out discrepancy check, dropped.

        let weight_p = old_weight_p_computed;

        let t_p = aq.value[rowindex_p as usize];
        for i in 0..rt.instance.num_var {
            if i != rowindex_p {
                let t_i = aq.value[i as usize];
                self.weights[i as usize] = self.weights[i as usize]
                    - 2.0 * (t_i / t_p) * delta.value[i as usize]
                    + ((t_i * t_i) / (t_p * t_p)) * weight_p;
                // C++ line 154: commented-out printf, dropped.
            }
        }
        // C++ lines 157-158: commented-out new_ep / computed_weight, dropped.
        let new_weight_p_updated = weight_p / (t_p * t_p);

        // C++ lines 161-166: commented-out discrepancy check, dropped.
        self.weights[rowindex_p as usize] = new_weight_p_updated;
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
    fn exact_weights_of_identity_basis_are_one() {
        // B = I, so btran(e_i) = e_i and every weight is e_i . e_i = 1.
        let mut rt = Runtime::new(&box_qp(), Statistics::default());
        let mut basis = at_lower(&mut rt);
        let pricing = SteepestEdgePricing::new(&rt, &mut basis);
        assert_eq!(pricing.weights, vec![1.0, 1.0]);
        assert!(pricing.check_weights(&rt, &mut basis));
    }

    #[test]
    fn price_drops_most_negative_multiplier() {
        let mut rt = Runtime::new(&box_qp(), Statistics::default());
        let mut basis = at_lower(&mut rt);
        let mut gradient = Gradient::new(&rt);
        let mut redcosts = ReducedCosts::new(&rt);
        let mut pricing = SteepestEdgePricing::new(&rt, &mut basis);
        let x = rt.primal.clone();
        let minidx = pricing.price(&rt, &mut basis, &mut redcosts, &mut gradient, &x);
        assert_eq!(minidx, 2);
    }

    #[test]
    fn update_weights_uses_ep_norm_then_recompute_restores() {
        // C++ lines 135-167 with B = I (delta = aq), p = 1, aq[r] = 2,
        // aq[other] = 3, ep = (1, 1) so weight_p = 2:
        //   weights[other] = 1 - 2 * (3 / 2) * 3 + (9 / 4) * 2 = -3.5
        //   weights[r]     = 2 / 4 = 0.5
        let mut rt = Runtime::new(&box_qp(), Statistics::default());
        let mut basis = at_lower(&mut rt);
        let mut pricing = SteepestEdgePricing::new(&rt, &mut basis);
        let r = basis.getindexinfactor()[1] as usize;
        let other = 1 - r;
        let mut aq = vec![0.0; 2];
        aq[r] = 2.0;
        aq[other] = 3.0;
        pricing.update_weights(&rt, &mut basis, &dense(&aq), &dense(&[1.0, 1.0]), 1, 0);
        assert_eq!(pricing.weights[r], 0.5);
        assert_eq!(pricing.weights[other], -3.5);
        assert!(!pricing.check_weight(&rt, &mut basis, other as HighsInt));
        pricing.recompute(&rt, &mut basis);
        assert_eq!(pricing.weights, vec![1.0, 1.0]);
    }
}
