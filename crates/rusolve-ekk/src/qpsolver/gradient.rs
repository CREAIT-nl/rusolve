//! HiGHS source: qpsolver/gradient.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! C++ stores `Runtime& runtime`; here `rt` is a parameter of the
//! methods that read it (lib.rs rule (a)).

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: gradient.hpp::Gradient (lines 14-44)
#[derive(Clone, Debug, Default)]
pub struct Gradient {
    gradient: QpVector,
    uptodate: bool,
    numupdates: HighsInt,
}

impl Gradient {
    /// HiGHS: gradient.hpp::Gradient::Gradient (lines 22-23)
    pub fn new(rt: &Runtime) -> Self {
        Self {
            gradient: QpVector::new(rt.instance.num_var),
            uptodate: false,
            numupdates: 0,
        }
    }

    /// HiGHS: gradient.hpp::Gradient::recompute (lines 25-30)
    pub fn recompute(&mut self, rt: &Runtime) {
        rt.instance.Q.vec_mat_into(&rt.primal, &mut self.gradient);
        self.gradient.op_add_assign(&rt.instance.c);
        self.uptodate = true;
        self.numupdates = 0;
    }

    /// HiGHS: gradient.hpp::Gradient::getGradient (lines 32-38)
    pub fn get_gradient(&mut self, rt: &Runtime) -> &mut QpVector {
        if !self.uptodate || self.numupdates >= rt.settings.gradientrecomputefrequency {
            self.recompute(rt);
        }
        &mut self.gradient
    }

    /// HiGHS: gradient.hpp::Gradient::update (lines 40-43)
    pub fn update(&mut self, buffer_Qp: &QpVector, stepsize: f64) {
        self.gradient.saxpy(stepsize, buffer_Qp);
        self.numupdates += 1;
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
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kActiveAtLower,
            BasisStatus::kActiveAtLower,
        ];
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
    fn gradient_is_c_plus_qx() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut g = Gradient::new(&rt);
        assert!(close(&g.get_gradient(&rt).value, &[-2.0, -4.0]));
        rt.primal.value = vec![1.0, 1.0];
        rt.primal.resparsify();
        g.recompute(&rt);
        assert!(close(&g.get_gradient(&rt).value, &[0.0, 0.0]));
    }
}
