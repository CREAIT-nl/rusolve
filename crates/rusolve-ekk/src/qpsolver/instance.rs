//! HiGHS source: qpsolver/instance.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `Instance(nv = 0, nc = 0)`: `Default` is `Instance::new(0, 0)`.

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: instance.hpp::SumNum (lines 16-19)
#[derive(Clone, Copy, Debug, Default)]
pub struct SumNum {
    pub sum: f64,
    pub num: HighsInt,
}

/// HiGHS: instance.hpp::Instance (lines 21-73)
#[derive(Clone, Debug)]
pub struct Instance {
    pub sense: HighsInt,
    pub num_var: HighsInt,
    pub num_con: HighsInt,
    pub offset: f64,
    pub c: QpVector,
    pub Q: Matrix,
    pub con_lo: Vec<f64>,
    pub con_up: Vec<f64>,
    pub A: Matrix,
    pub var_lo: Vec<f64>,
    pub var_up: Vec<f64>,
}

impl Default for Instance {
    /// HiGHS: instance.hpp::Instance::Instance (default arguments nv = 0, nc = 0)
    fn default() -> Self {
        Self::new(0, 0)
    }
}

impl Instance {
    /// HiGHS: instance.hpp::Instance::Instance (lines 34-39)
    ///
    /// C++ default arguments: nv = 0, nc = 0
    pub fn new(nv: HighsInt, nc: HighsInt) -> Self {
        Instance {
            sense: 1,
            num_var: nv,
            num_con: nc,
            offset: 0.0,
            c: QpVector::new(nv),
            Q: Matrix::new(nv, nv),
            con_lo: Vec::new(),
            con_up: Vec::new(),
            A: Matrix::new(nc, nv),
            var_lo: Vec::new(),
            var_up: Vec::new(),
        }
    }

    /// HiGHS: instance.hpp::Instance::objval (lines 41-43)
    pub fn objval(&self, x: &QpVector) -> f64 {
        self.c.dot(x) + 0.5 * (self.Q.vec_mat(x).dot(x)) + self.offset
    }

    /// HiGHS: instance.hpp::Instance::isEquality (lines 45-49)
    pub fn is_equality(&self, i_var: HighsInt) -> bool {
        if i_var < self.num_con {
            self.con_lo[i_var as usize] == self.con_up[i_var as usize]
        } else {
            self.var_lo[(i_var - self.num_con) as usize]
                == self.var_up[(i_var - self.num_con) as usize]
        }
    }

    /// HiGHS: instance.hpp::Instance::sumnumprimalinfeasibilities (lines 50-72)
    pub fn sumnumprimalinfeasibilities(&self, x: &QpVector, rowactivity: &QpVector) -> SumNum {
        let mut res = SumNum::default();
        for row in 0..self.num_con {
            if rowactivity.value[row as usize] < self.con_lo[row as usize] {
                res.sum += self.con_lo[row as usize] - rowactivity.value[row as usize];
                res.num += 1;
            } else if rowactivity.value[row as usize] > self.con_up[row as usize] {
                res.sum += rowactivity.value[row as usize] - self.con_up[row as usize];
                res.num += 1;
            }
        }
        for var in 0..self.num_var {
            if x.value[var as usize] < self.var_lo[var as usize] {
                res.sum += self.var_lo[var as usize] - x.value[var as usize];
                res.num += 1;
            } else if x.value[var as usize] > self.var_up[var as usize] {
                res.sum += x.value[var as usize] - self.var_up[var as usize];
                res.num += 1;
            }
        }
        res
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

    #[test]
    fn objval_and_infeasibilities() {
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
        inst.con_up = vec![1.0];
        inst.var_lo = vec![0.0, 0.0];
        inst.var_up = vec![10.0, 10.0];
        let mut x = QpVector::new(2);
        x.value = vec![1.0, 1.0];
        x.resparsify();
        // -2 - 4 + 0.5 (2 + 4)
        assert_eq!(inst.objval(&x), -3.0);
        let ra = inst.A.mat_vec(&x);
        let s = inst.sumnumprimalinfeasibilities(&x, &ra);
        assert_eq!((s.num, s.sum), (1, 1.0));
    }

    #[test]
    fn is_equality_and_lower_violations() {
        let mut inst = Instance::new(2, 1);
        inst.con_lo = vec![2.0];
        inst.con_up = vec![2.0];
        inst.var_lo = vec![0.0, 0.0];
        inst.var_up = vec![2.0, 2.0];
        assert!(inst.is_equality(0));
        assert!(!inst.is_equality(1));
        let mut ra = QpVector::new(1);
        ra.value = vec![1.0];
        ra.resparsify();
        let mut x = QpVector::new(2);
        x.value = vec![3.0, -0.5];
        x.resparsify();
        // row below lo (1), var0 above up (1), var1 below lo (0.5)
        let s = inst.sumnumprimalinfeasibilities(&x, &ra);
        assert_eq!((s.num, s.sum), (3, 2.5));
    }
}
