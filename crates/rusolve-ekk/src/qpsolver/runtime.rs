//! HiGHS source: qpsolver/runtime.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//!
//!
//! PORT NOTE (ownership): C++ `Runtime` holds `Statistics& statistics`
//! (a reference to the caller's `Statistics`). Here `Runtime` OWNS it; the
//! caller (`solveqp_actual`, a_asm.cpp) moves its `Statistics` in with
//! `std::mem::take` and moves it back out after `Quass::solve`, which is
//! observably identical (nothing else touches the caller's object meanwhile).

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: runtime.hpp::Runtime (lines 17-43)
#[derive(Clone, Debug)]
pub struct Runtime {
    pub instance: Instance,
    pub relaxed_for_ratiotest: Instance,
    pub scaled: Instance,
    pub perturbed: Instance,
    pub settings: Settings,
    pub statistics: Statistics,

    pub primal: QpVector,
    pub rowactivity: QpVector,
    pub dualvar: QpVector,
    pub dualcon: QpVector,
    pub status: QpModelStatus,

    pub status_var: Vec<BasisStatus>,
    pub status_con: Vec<BasisStatus>,
}

impl Runtime {
    /// HiGHS: runtime.hpp::Runtime::Runtime (lines 34-42)
    ///
    /// `status` starts at `QpModelStatus::kUndetermined` (line 29).
    pub fn new(inst: &Instance, stats: Statistics) -> Self {
        Runtime {
            instance: inst.clone(),
            relaxed_for_ratiotest: Instance::default(),
            scaled: Instance::default(),
            perturbed: Instance::default(),
            settings: Settings::default(),
            statistics: stats,
            primal: QpVector::new(inst.num_var),
            rowactivity: QpVector::new(inst.num_con),
            dualvar: QpVector::new(inst.num_var),
            dualcon: QpVector::new(inst.num_con),
            status: QpModelStatus::kUndetermined,
            status_var: vec![BasisStatus::default(); inst.num_var as usize],
            status_con: vec![BasisStatus::default(); inst.num_con as usize],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_sizes_and_defaults() {
        let inst = Instance::new(3, 2);
        let rt = Runtime::new(&inst, Statistics::default());
        assert_eq!(rt.instance.num_var, 3);
        assert_eq!(rt.instance.num_con, 2);
        assert_eq!((rt.primal.dim, rt.dualvar.dim), (3, 3));
        assert_eq!((rt.rowactivity.dim, rt.dualcon.dim), (2, 2));
        assert_eq!(rt.status, QpModelStatus::kUndetermined);
        assert_eq!(rt.status_var, vec![BasisStatus::kInactive; 3]);
        assert_eq!(rt.status_con, vec![BasisStatus::kInactive; 2]);
        assert_eq!(
            (rt.scaled.num_var, rt.scaled.num_con, rt.scaled.sense),
            (0, 0, 1)
        );
    }
}
