//! HiGHS source: qpsolver/perturbation.hpp, qpsolver/perturbation.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Uses `std::minstd_rand0` and `std::uniform_real_distribution<double>(1e-5, 1e-4)`;
//! `Settings::perturbation` is false by default, so this is off the default
//! path. If ported bit-exactly it needs a libstdc++-faithful
//! `uniform_real_distribution` (generate_canonical) over minstd_rand0.

use crate::prelude::*;
use crate::qpsolver::*;

// PORT NOTE: std::default_random_engine = libstdc++ minstd_rand0
/// HiGHS: perturbation.cpp::perturb (lines 12-41)
struct MinstdRand0 {
    state: u64,
}

impl MinstdRand0 {
    fn new() -> Self {
        MinstdRand0 { state: 1 }
    }

    fn next(&mut self) -> u64 {
        self.state = (self.state * 16807) % 2147483647;
        self.state
    }
}

/// HiGHS: perturbation.cpp::perturb (lines 12-41)
/// PORT NOTE: std::uniform_real_distribution<double>(1e-5, 1e-4) via generate_canonical
fn uniform_real_distribution(re: &mut MinstdRand0) -> f64 {
    let a: f64 = 1e-5;
    let b: f64 = 1e-4;
    // generate_canonical<double, 53>(re)
    let r: f64 = 2147483646.0;
    let k: i32 = 2;
    let mut sum: f64 = 0.0;
    let mut tmp: f64 = 1.0;
    for _ in 0..k {
        let x = re.next();
        sum += ((x - 1) as f64) * tmp;
        tmp *= r;
    }
    let mut u: f64 = sum / tmp;
    if u >= 1.0 {
        u = 1.0 - 2.0f64.powi(-53);
    }
    u * (b - a) + a
}

/// HiGHS: perturbation.cpp::perturb (lines 12-41)
pub fn perturb(rt: &mut Runtime) {
    rt.perturbed = rt.instance.clone();
    if !rt.settings.perturbation {
        return;
    }

    let mut re = MinstdRand0::new();

    let mut i: HighsInt = 0;
    while i < rt.perturbed.num_con {
        if rt.perturbed.con_lo[i as usize] != rt.perturbed.con_up[i as usize] {
            if rt.perturbed.con_lo[i as usize] != -f64::INFINITY {
                rt.perturbed.con_lo[i as usize] -= uniform_real_distribution(&mut re);
            }
            if rt.perturbed.con_up[i as usize] != f64::INFINITY {
                rt.perturbed.con_up[i as usize] += uniform_real_distribution(&mut re);
            }
        }
        i += 1;
    }

    let mut i: HighsInt = 0;
    while i < rt.perturbed.num_var {
        if rt.perturbed.var_lo[i as usize] != rt.perturbed.var_up[i as usize] {
            if rt.perturbed.var_lo[i as usize] != -f64::INFINITY {
                rt.perturbed.var_lo[i as usize] -= uniform_real_distribution(&mut re);
            }
            if rt.perturbed.var_up[i as usize] != f64::INFINITY {
                rt.perturbed.var_up[i as usize] += uniform_real_distribution(&mut re);
            }
        }
        i += 1;
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
    fn perturbation_off_copies_instance() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        assert!(!rt.settings.perturbation);
        perturb(&mut rt);
        assert_eq!(rt.perturbed.con_up, rt.instance.con_up);
        assert_eq!(rt.perturbed.var_lo, rt.instance.var_lo);
    }
}
