//! HiGHS source: qpsolver/a_asm.hpp, qpsolver/a_asm.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//!
//!
//! PORT NOTE (ownership): `solveqp_actual` builds `Runtime rt(instance, stats)`
//! whose `statistics` is a C++ reference to `stats`. The Rust `Runtime` owns
//! its `Statistics` (see `runtime.rs`): move `stats` in with
//! `std::mem::take(stats)` and move `rt.statistics` back into `*stats` before
//! every return. Likewise `Runtime::instance` is a C++ `Instance&` to the
//! caller's `instance` but an owned clone here (`Runtime::new`): move it back
//! into `*instance` with `std::mem::take` before the return, so any change
//! `Quass::solve` makes to `runtime.instance` reaches the caller as in C++.

use crate::prelude::*;
use crate::qpsolver::*;
use crate::util::highs_cdouble::HighsCDouble;
use crate::util::highs_utils::cpp_max;

/// HiGHS: a_asm.hpp::QpAsmStatus (lines 17-23)
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QpAsmStatus {
    kOk = 0,
    kWarning,
    kError,
}

/// HiGHS: a_asm.hpp::QpSolution (lines 25-41)
#[derive(Clone, Debug)]
pub struct QpSolution {
    pub primal: QpVector,
    pub rowactivity: QpVector,
    pub dualvar: QpVector,
    pub dualcon: QpVector,

    pub status_var: Vec<BasisStatus>,
    pub status_con: Vec<BasisStatus>,
}

impl QpSolution {
    /// HiGHS: a_asm.hpp::QpSolution::QpSolution (lines 34-40)
    pub fn new(instance: &Instance) -> Self {
        QpSolution {
            primal: QpVector::new(instance.num_var),
            rowactivity: QpVector::new(instance.num_con),
            dualvar: QpVector::new(instance.num_var),
            dualcon: QpVector::new(instance.num_con),
            status_var: vec![BasisStatus::kInactive; instance.num_var as usize],
            status_con: vec![BasisStatus::kInactive; instance.num_con as usize],
        }
    }
}

/// HiGHS: a_asm.hpp::QpHotstartInformation (lines 43-52)
#[derive(Clone, Debug)]
pub struct QpHotstartInformation {
    pub active: Vec<HighsInt>,
    pub inactive: Vec<HighsInt>,
    pub status: Vec<BasisStatus>,
    pub primal: QpVector,
    pub rowact: QpVector,
}

impl QpHotstartInformation {
    /// HiGHS: a_asm.hpp::QpHotstartInformation::QpHotstartInformation (lines 50-51)
    pub fn new(num_var: HighsInt, num_row: HighsInt) -> Self {
        QpHotstartInformation {
            active: Vec::new(),
            inactive: Vec::new(),
            status: Vec::new(),
            primal: QpVector::new(num_var),
            rowact: QpVector::new(num_row),
        }
    }
}

/// HiGHS: a_asm.cpp::solveqp_actual (lines 13-35)
pub fn solveqp_actual(
    instance: &mut Instance,
    settings: &mut Settings,
    startinfo: &mut QpHotstartInformation,
    stats: &mut Statistics,
    status: &mut QpModelStatus,
    solution: &mut QpSolution,
    qp_timer: &mut HighsTimer,
) -> QpAsmStatus {
    let mut rt = Runtime::new(instance, std::mem::take(stats));
    rt.settings = settings.clone();
    let mut quass = Quass::new(&rt);

    let mut basis = Basis::new(
        &mut rt,
        startinfo.active.clone(),
        startinfo.status.clone(),
        startinfo.inactive.clone(),
    );

    quass.solve(
        &mut rt,
        &startinfo.primal,
        &startinfo.rowact,
        &mut basis,
        qp_timer,
    );

    *status = rt.status;

    solution.status_var = rt.status_var.clone();
    solution.status_con = rt.status_con.clone();
    solution.primal = rt.primal.clone();
    solution.dualvar = rt.dualvar.clone();
    solution.rowactivity = rt.rowactivity.clone();
    solution.dualcon = rt.dualcon.clone();

    // PORT NOTE: write back the C++ references (`Runtime::statistics`,
    // `Runtime::instance`); see module doc.
    *stats = std::mem::take(&mut rt.statistics);
    *instance = std::mem::take(&mut rt.instance);

    QpAsmStatus::kOk
}

/// HiGHS: a_asm.cpp::qpBasisStatusToString (lines 37-50)
pub fn qp_basis_status_to_string(qp_basis_status: BasisStatus) -> String {
    match qp_basis_status {
        BasisStatus::kInactive => "Inactive".to_string(),
        BasisStatus::kActiveAtLower => "Active at lower bound".to_string(),
        BasisStatus::kActiveAtUpper => "Active at upper bound".to_string(),
        BasisStatus::kInactiveInBasis => "Inactive in basis".to_string(),
    }
}

/// HiGHS: a_asm.cpp::qpModelStatusToString (lines 52-77)
pub fn qp_model_status_to_string(qp_model_status: QpModelStatus) -> String {
    match qp_model_status {
        QpModelStatus::kNotset => "Not set".to_string(),
        QpModelStatus::kUndetermined => "Undetermined".to_string(),
        QpModelStatus::kOptimal => "Optimal".to_string(),
        QpModelStatus::kUnbounded => "Unbounded".to_string(),
        QpModelStatus::kInfeasible => "Infeasible".to_string(),
        QpModelStatus::kIterationLimit => "Iteration limit".to_string(),
        QpModelStatus::kTimeLimit => "Time ;limit".to_string(),
        QpModelStatus::kLargeNullspace => "Large nullspace".to_string(),
        QpModelStatus::kNonConvex => "Non-convex".to_string(),
        QpModelStatus::kError => "Error".to_string(),
        _ => "Unidentified QP model status".to_string(),
    }
}

/// HiGHS: a_asm.cpp::assessQpPrimalFeasibility (lines 79-141)
pub fn assess_qp_primal_feasibility(
    instance: &Instance,
    primal_feasibility_tolerance: f64,
    var_value: &[f64],
    con_value: &[f64],
    num_var_infeasibilities: &mut HighsInt,
    max_var_infeasibility: &mut f64,
    sum_var_infeasibilities: &mut f64,
    num_con_infeasibilities: &mut HighsInt,
    max_con_infeasibility: &mut f64,
    sum_con_infeasibilities: &mut f64,
    max_con_residual: &mut f64,
    sum_con_residuals: &mut f64,
) {
    *num_var_infeasibilities = 0;
    *max_var_infeasibility = 0.0;
    *sum_var_infeasibilities = 0.0;
    *num_con_infeasibilities = 0;
    *max_con_infeasibility = 0.0;
    *sum_con_infeasibilities = 0.0;
    *max_con_residual = 0.0;
    *sum_con_residuals = 0.0;
    // Valid solution, but is it feasible?
    let mut con_value_quad: Vec<HighsCDouble> =
        vec![HighsCDouble::from(0.0); instance.num_con as usize];
    for i_var in 0..instance.num_var {
        let lower = instance.var_lo[i_var as usize];
        let upper = instance.var_up[i_var as usize];
        let primal = var_value[i_var as usize];
        let mut var_infeasibility = 0.0;
        if primal < lower - primal_feasibility_tolerance {
            var_infeasibility = lower - primal;
        } else if primal > upper + primal_feasibility_tolerance {
            var_infeasibility = primal - upper;
        }
        if var_infeasibility > 0.0 {
            if var_infeasibility > primal_feasibility_tolerance {
                *num_var_infeasibilities += 1;
            }
            *max_var_infeasibility = cpp_max(var_infeasibility, *max_var_infeasibility);
            *sum_var_infeasibilities += var_infeasibility;
        }
        for i_el in instance.A.mat.start[i_var as usize]..instance.A.mat.start[(i_var + 1) as usize]
        {
            con_value_quad[instance.A.mat.index[i_el as usize] as usize] +=
                primal * instance.A.mat.value[i_el as usize];
        }
    }
    for i_con in 0..instance.num_con {
        let lower = instance.con_lo[i_con as usize];
        let upper = instance.con_up[i_con as usize];
        let primal = con_value[i_con as usize];
        let mut con_infeasibility = 0.0;
        if primal < lower - primal_feasibility_tolerance {
            con_infeasibility = lower - primal;
        } else if primal > upper + primal_feasibility_tolerance {
            con_infeasibility = primal - upper;
        }
        if con_infeasibility > 0.0 {
            if con_infeasibility > primal_feasibility_tolerance {
                *num_con_infeasibilities += 1;
            }
            *max_con_infeasibility = cpp_max(con_infeasibility, *max_con_infeasibility);
            *sum_con_infeasibilities += con_infeasibility;
        }
        let con_residual = (primal - f64::from(con_value_quad[i_con as usize])).abs();
        *max_con_residual = cpp_max(con_residual, *max_con_residual);
        *sum_con_residuals += con_residual;
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
    /// `num_con + i`); the row (id 0) is not in the basis. `status[i]` is the
    /// status of `active[i]` (basis.cpp:27-30), not of constraint id `i`.
    fn at_lower(rt: &mut Runtime) -> Basis {
        let status = vec![BasisStatus::kActiveAtLower, BasisStatus::kActiveAtLower];
        Basis::new(rt, vec![1, 2], status, vec![])
    }

    /// Both variables in the nullspace (inactive in basis): Z = I.
    fn all_free(rt: &mut Runtime) -> Basis {
        Basis::new(rt, vec![], vec![], vec![1, 2])
    }

    fn close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-9)
    }

    #[test]
    fn status_strings() {
        assert_eq!(
            qp_model_status_to_string(QpModelStatus::kTimeLimit),
            "Time ;limit"
        );
        assert_eq!(
            qp_basis_status_to_string(BasisStatus::kActiveAtUpper),
            "Active at upper bound"
        );
    }

    fn box_start(status: Vec<BasisStatus>) -> QpHotstartInformation {
        let mut start = QpHotstartInformation::new(2, 1);
        start.active = vec![1, 2];
        start.inactive = vec![];
        start.status = status;
        start.primal = QpVector::new(2);
        start.rowact = QpVector::new(1);
        start
    }

    fn run_box(
        start: &mut QpHotstartInformation,
    ) -> (QpAsmStatus, QpModelStatus, Statistics, QpSolution) {
        let mut inst = box_qp(true);
        let mut settings = Settings::default();
        // C++ Settings default iteration_limit is numeric_limits<HighsInt>::infinity() == 0
        // (settings.hpp:76), so the default would stop at once; Highs::callSolveQp always
        // sets it from qp_iteration_limit. Set a generous limit here (test fix, see commit).
        settings.iteration_limit = 1000;
        let mut stats = Statistics::default();
        let mut status = QpModelStatus::kNotset;
        let mut solution = QpSolution::new(&inst);
        let mut timer = HighsTimer::new();
        let r = solveqp_actual(
            &mut inst,
            &mut settings,
            start,
            &mut stats,
            &mut status,
            &mut solution,
            &mut timer,
        );
        // The instance written back is the one the caller passed in.
        assert_eq!((inst.num_var, inst.num_con), (2, 1));
        assert_eq!(inst.Q.mat.value, vec![2.0, 4.0]);
        (r, status, stats, solution)
    }

    /// Expected values are the output of the real HiGHS 1.15.1
    /// `solveqp_actual` (oracle build, -ffp-contract=off) on this instance and
    /// start: kOptimal after 4 iterations, x = (1, 1), row activity 2, both
    /// variables kInactiveInBasis, all duals 0. Exact equality is intended.
    #[test]
    fn solveqp_actual_box_qp() {
        // QpHotstartInformation::status is parallel to `active`
        // (basis.cpp:27-30: basisstatus[active[i]] = status[i]). The skeleton
        // test indexed it by constraint id ([kInactive, kAtLower, kAtLower]),
        // which in C++ marks x0 kInactive while it sits in the active set.
        let mut start = box_start(vec![
            BasisStatus::kActiveAtLower,
            BasisStatus::kActiveAtLower,
        ]);
        let (r, status, stats, solution) = run_box(&mut start);
        assert_eq!(r, QpAsmStatus::kOk);
        assert_eq!(status, QpModelStatus::kOptimal);
        assert_eq!(stats.num_iterations, 4);
        assert_eq!(solution.primal.value, vec![1.0, 1.0]);
        assert!(close(&solution.primal.value, &[1.0, 1.0]));
        assert_eq!(solution.rowactivity.value, vec![2.0]);
        assert_eq!(
            solution.status_var,
            vec![BasisStatus::kInactiveInBasis, BasisStatus::kInactiveInBasis]
        );
        assert_eq!(solution.status_con, vec![BasisStatus::kInactive]);
        assert_eq!(solution.dualvar.value, vec![0.0, 0.0]);
        assert_eq!(solution.dualcon.value, vec![0.0]);
    }

    /// The skeleton's constraint-id-indexed status vector, kept as a
    /// regression: real HiGHS 1.15.1 returns kOptimal after 2 iterations with
    /// x = (0, 1) (x0 is active but has status kInactive, so it is never
    /// released), row activity 1, status_var [kInactive, kInactiveInBasis],
    /// dualvar (-2, 0).
    #[test]
    fn solveqp_actual_box_qp_status_by_constraint_id() {
        let mut start = box_start(vec![
            BasisStatus::kInactive,
            BasisStatus::kActiveAtLower,
            BasisStatus::kActiveAtLower,
        ]);
        let (r, status, stats, solution) = run_box(&mut start);
        assert_eq!(r, QpAsmStatus::kOk);
        assert_eq!(status, QpModelStatus::kOptimal);
        assert_eq!(stats.num_iterations, 2);
        assert_eq!(solution.primal.value, vec![0.0, 1.0]);
        assert_eq!(solution.rowactivity.value, vec![1.0]);
        assert_eq!(
            solution.status_var,
            vec![BasisStatus::kInactive, BasisStatus::kInactiveInBasis]
        );
        assert_eq!(solution.status_con, vec![BasisStatus::kInactive]);
        assert_eq!(solution.dualvar.value, vec![-2.0, 0.0]);
        assert_eq!(solution.dualcon.value, vec![0.0]);
    }

    #[allow(clippy::type_complexity)]
    fn assess(var: &[f64], con: &[f64]) -> (HighsInt, f64, f64, HighsInt, f64, f64, f64, f64) {
        let inst = box_qp(true);
        let (mut nv, mut mv, mut sv) = (-1, -1.0, -1.0);
        let (mut nc, mut mc, mut sc) = (-1, -1.0, -1.0);
        let (mut mr, mut sr) = (-1.0, -1.0);
        assess_qp_primal_feasibility(
            &inst, 1e-7, var, con, &mut nv, &mut mv, &mut sv, &mut nc, &mut mc, &mut sc, &mut mr,
            &mut sr,
        );
        (nv, mv, sv, nc, mc, sc, mr, sr)
    }

    /// Expected values are the output of the real HiGHS 1.15.1
    /// `assessQpPrimalFeasibility` (oracle build) on the same inputs.
    #[test]
    fn assess_primal_feasibility_matches_cpp() {
        // Both bound sides violated, row above its upper bound, residual 0.5.
        assert_eq!(
            assess(&[-1.0, 12.0], &[11.5]),
            (2, 2.0, 3.0, 1, 1.5, 1.5, 0.5, 0.5)
        );
        // Violations within tolerance are not counted; residual |3 - 10| = 7.
        assert_eq!(
            assess(&[-5e-8, 10.0 + 5e-8], &[3.0]),
            (0, 0.0, 0.0, 0, 0.0, 0.0, 7.0, 7.0)
        );
        // Row activity accumulated in HighsCDouble: 0.1 + 0.2 rounds to
        // 0.30000000000000004, so the residual is exactly 0.
        assert_eq!(
            assess(&[0.1, 0.2], &[0.30000000000000004]),
            (0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0)
        );
        // std::max(NaN, 0) returns its first argument: max_con_residual is NaN.
        let r = assess(&[0.5, 0.5], &[f64::NAN]);
        assert_eq!((r.0, r.1, r.2, r.3, r.4, r.5), (0, 0.0, 0.0, 0, 0.0, 0.0));
        assert!(r.6.is_nan() && r.7.is_nan());
    }
}
