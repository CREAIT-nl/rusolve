//! HiGHS source: qpsolver/a_quass.hpp, qpsolver/a_quass.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `Eventhandler::fire` takes `&mut T`: `quass2highs` fires
//! `settings.qp_model_status_log` on `qp_model_status` (C++ passes `T&`).

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: a_quass.cpp::quass2highs (lines 14-128)
/// HiGHS 1.16 (#3153): a_quass.cpp::quass2highs (lines 14-132), kNonConvex arm only
fn quass2highs(
    instance: &mut Instance,
    settings: &mut Settings,
    stats: &mut Statistics,
    qp_model_status: &mut QpModelStatus,
    qp_solution: &mut QpSolution,
    highs_model_status: &mut HighsModelStatus,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
) -> QpAsmStatus {
    settings.qp_model_status_log.fire(qp_model_status);
    let qp_asm_return_status = match *qp_model_status {
        QpModelStatus::kOptimal => {
            *highs_model_status = HighsModelStatus::kOptimal;
            QpAsmStatus::kOk
        }
        QpModelStatus::kUnbounded => {
            *highs_model_status = HighsModelStatus::kUnbounded;
            QpAsmStatus::kOk
        }
        QpModelStatus::kInfeasible => {
            *highs_model_status = HighsModelStatus::kInfeasible;
            QpAsmStatus::kOk
        }
        QpModelStatus::kIterationLimit => {
            *highs_model_status = HighsModelStatus::kIterationLimit;
            QpAsmStatus::kWarning
        }
        QpModelStatus::kTimeLimit => {
            *highs_model_status = HighsModelStatus::kTimeLimit;
            QpAsmStatus::kWarning
        }
        QpModelStatus::kInterrupt => {
            *highs_model_status = HighsModelStatus::kInterrupt;
            QpAsmStatus::kWarning
        }
        QpModelStatus::kUndetermined => {
            *highs_model_status = HighsModelStatus::kSolveError;
            return QpAsmStatus::kError;
        }
        QpModelStatus::kLargeNullspace => {
            *highs_model_status = HighsModelStatus::kSolveError;
            return QpAsmStatus::kError;
        }
        // HiGHS 1.16 #3153 (h16_qp_bounded_start): HiGHS 1.16 a_quass.cpp lines 55-57 (kNonConvex arm).
        // Only reachable with the flag on; with it off kNonConvex takes the `_` arm as before.
        QpModelStatus::kNonConvex if settings.h16_qp_bounded_start => {
            *highs_model_status = HighsModelStatus::kSolveError;
            return QpAsmStatus::kError;
        }
        QpModelStatus::kError => {
            *highs_model_status = HighsModelStatus::kSolveError;
            return QpAsmStatus::kError;
        }
        QpModelStatus::kNotset => {
            *highs_model_status = HighsModelStatus::kNotset;
            return QpAsmStatus::kError;
        }
        _ => {
            *highs_model_status = HighsModelStatus::kNotset;
            return QpAsmStatus::kError;
        }
    };

    debug_assert!(qp_asm_return_status != QpAsmStatus::kError);
    // extract variable values
    highs_solution
        .col_value
        .resize(instance.num_var as usize, 0.0);
    highs_solution
        .col_dual
        .resize(instance.num_var as usize, 0.0);
    for i_col in 0..instance.num_var {
        highs_solution.col_value[i_col as usize] = qp_solution.primal.value[i_col as usize];
        highs_solution.col_dual[i_col as usize] =
            instance.sense as f64 * qp_solution.dualvar.value[i_col as usize];
    }
    // extract constraint activity
    highs_solution
        .row_value
        .resize(instance.num_con as usize, 0.0);
    highs_solution
        .row_dual
        .resize(instance.num_con as usize, 0.0);
    // Negate the vector and Hessian
    for i_row in 0..instance.num_con {
        highs_solution.row_value[i_row as usize] = qp_solution.rowactivity.value[i_row as usize];
        highs_solution.row_dual[i_row as usize] =
            instance.sense as f64 * qp_solution.dualcon.value[i_row as usize];
    }
    highs_solution.value_valid = true;
    highs_solution.dual_valid = true;

    // extract basis status
    highs_basis
        .col_status
        .resize(instance.num_var as usize, HighsBasisStatus::kLower);
    highs_basis
        .row_status
        .resize(instance.num_con as usize, HighsBasisStatus::kLower);

    // lines 91-92, 94-96, 110-112: `debug_report` printf (compile-time false) dropped.
    for i in 0..instance.num_var {
        if qp_solution.status_var[i as usize] == BasisStatus::kActiveAtLower {
            highs_basis.col_status[i as usize] = HighsBasisStatus::kLower;
        } else if qp_solution.status_var[i as usize] == BasisStatus::kActiveAtUpper {
            highs_basis.col_status[i as usize] = HighsBasisStatus::kUpper;
        } else if qp_solution.status_var[i as usize] == BasisStatus::kInactiveInBasis {
            highs_basis.col_status[i as usize] = HighsBasisStatus::kNonbasic;
        } else {
            debug_assert!(qp_solution.status_var[i as usize] == BasisStatus::kInactive);
            highs_basis.col_status[i as usize] = HighsBasisStatus::kBasic;
        }
    }

    for i in 0..instance.num_con {
        if qp_solution.status_con[i as usize] == BasisStatus::kActiveAtLower {
            highs_basis.row_status[i as usize] = HighsBasisStatus::kLower;
        } else if qp_solution.status_con[i as usize] == BasisStatus::kActiveAtUpper {
            highs_basis.row_status[i as usize] = HighsBasisStatus::kUpper;
        } else if qp_solution.status_con[i as usize] == BasisStatus::kInactiveInBasis {
            highs_basis.row_status[i as usize] = HighsBasisStatus::kNonbasic;
        } else {
            debug_assert!(qp_solution.status_con[i as usize] == BasisStatus::kInactive);
            highs_basis.row_status[i as usize] = HighsBasisStatus::kBasic;
        }
    }
    highs_basis.valid = true;
    highs_basis.alien = false;
    highs_basis.useful = true;
    qp_asm_return_status
}

/// HiGHS: a_quass.cpp::solveqp (lines 130-194)
/// HiGHS 1.16 (#3153): a_quass.cpp::solveqp (lines 167-189, bounded-start gate and return block);
/// the 1.16 regularisation wrapper (lines 148-162) and callback block (200-219) are not ported
pub fn solveqp(
    instance: &mut Instance,
    settings: &mut Settings,
    stats: &mut Statistics,
    highs_model_status: &mut HighsModelStatus,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
    qp_timer: &mut HighsTimer,
) -> QpAsmStatus {
    let mut qp_model_status = QpModelStatus::kUndetermined;

    let mut qp_solution = QpSolution::new(instance);

    // presolve

    // scale instance, store scaling factors

    // perturb instance, store perturbance information

    // regularize
    for i in 0..instance.num_var {
        let mut index = instance.Q.mat.start[i as usize];
        while index < instance.Q.mat.start[(i + 1) as usize] {
            if instance.Q.mat.index[index as usize] == i {
                instance.Q.mat.value[index as usize] += settings.hessian_regularization_value;
            }
            index += 1;
        }
    }

    // compute initial feasible point
    let mut startinfo = QpHotstartInformation::new(instance.num_var, instance.num_con);
    // HiGHS 1.16 #3153 (h16_qp_bounded_start): flag read once at the gate
    let h16_qp_bounded_start: bool = settings.h16_qp_bounded_start;
    // HiGHS 1.16 (#3153): a_quass.cpp::solveqp (lines 167-168), bounded-start gate
    let bounded_start_gate: bool = if h16_qp_bounded_start {
        instance.num_con == 0 && (instance.num_var <= 5000 || instance.Q.mat.is_diagonal())
    } else {
        instance.num_con == 0 && instance.num_var <= 15000
    };
    if bounded_start_gate {
        compute_starting_point_bounded(
            instance,
            settings,
            stats,
            &mut qp_model_status,
            &mut startinfo,
            qp_timer,
        );
        if h16_qp_bounded_start {
            // HiGHS 1.16 (#3153): a_quass.cpp::solveqp (lines 175-189)
            if qp_model_status == QpModelStatus::kOptimal
                || qp_model_status == QpModelStatus::kUnbounded
                || qp_model_status == QpModelStatus::kNonConvex
            {
                if qp_model_status == QpModelStatus::kOptimal {
                    // Optimality occurs only if the unconstrained solution is
                    // feasible, so copy in the primal solution, zero the column
                    // duals and set the column status to inactive
                    qp_solution.primal = startinfo.primal.clone();
                    // C++ assign(num_var, 0) replaces the contents
                    qp_solution.dualvar.value = vec![0.0; instance.num_var as usize];
                    qp_solution.status_var =
                        vec![BasisStatus::kInactive; instance.num_var as usize];
                }
                return quass2highs(
                    instance,
                    settings,
                    stats,
                    &mut qp_model_status,
                    &mut qp_solution,
                    highs_model_status,
                    highs_basis,
                    highs_solution,
                );
            }
        } else {
            if qp_model_status == QpModelStatus::kOptimal {
                qp_solution.primal = startinfo.primal.clone();
                return quass2highs(
                    instance,
                    settings,
                    stats,
                    &mut qp_model_status,
                    &mut qp_solution,
                    highs_model_status,
                    highs_basis,
                    highs_solution,
                );
            }
            if qp_model_status == QpModelStatus::kUnbounded {
                return quass2highs(
                    instance,
                    settings,
                    stats,
                    &mut qp_model_status,
                    &mut qp_solution,
                    highs_model_status,
                    highs_basis,
                    highs_solution,
                );
            }
        }
    } else {
        compute_starting_point_highs(
            instance,
            settings,
            stats,
            &mut qp_model_status,
            &mut startinfo,
            highs_model_status,
            highs_basis,
            highs_solution,
            qp_timer,
        );
        if qp_model_status != QpModelStatus::kNotset {
            return quass2highs(
                instance,
                settings,
                stats,
                &mut qp_model_status,
                &mut qp_solution,
                highs_model_status,
                highs_basis,
                highs_solution,
            );
        }
    }

    // solve
    solveqp_actual(
        instance,
        settings,
        &mut startinfo,
        stats,
        &mut qp_model_status,
        &mut qp_solution,
        qp_timer,
    );

    // undo perturbation and resolve

    // undo scaling and resolve

    // postsolve

    // Transform QP status and qp_solution to HiGHS highs_basis and highs_solution
    quass2highs(
        instance,
        settings,
        stats,
        &mut qp_model_status,
        &mut qp_solution,
        highs_model_status,
        highs_basis,
        highs_solution,
    )
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
    fn solveqp_box_qp() {
        let mut inst = box_qp(true);
        let mut settings = Settings::default();
        // C++ settings.hpp:76 `iteration_limit = numeric_limits<HighsInt>::infinity()`
        // is 0 for an integer type, so a default Settings stops at iteration 0
        // (quass.cpp:353). Callers set the limit; the test does the same.
        settings.iteration_limit = 1000;
        // solveqp adds hessian_regularization_value (1e-7) to the Hessian diagonal;
        // the hand-derived optimum is for the unregularized problem, so turn it off.
        settings.hessian_regularization_value = 0.0;
        let mut stats = Statistics::default();
        let mut hms = HighsModelStatus::kNotset;
        let mut hb = HighsBasis::default();
        let mut hs = HighsSolution::default();
        let mut timer = HighsTimer::new();
        let r = solveqp(
            &mut inst,
            &mut settings,
            &mut stats,
            &mut hms,
            &mut hb,
            &mut hs,
            &mut timer,
        );
        assert_eq!(r, QpAsmStatus::kOk);
        assert_eq!(hms, HighsModelStatus::kOptimal);
        assert!(close(&hs.col_value, &[1.0, 1.0]));
        assert!(close(&hs.row_value, &[2.0]));
    }

    /// Runs solveqp and returns its return status, model status, basis and solution.
    fn run_solveqp(
        inst: &mut Instance,
        settings: &mut Settings,
    ) -> (QpAsmStatus, HighsModelStatus, HighsBasis, HighsSolution) {
        let mut stats = Statistics::default();
        let mut hms = HighsModelStatus::kNotset;
        let mut hb = HighsBasis::default();
        let mut hs = HighsSolution::default();
        let mut timer = HighsTimer::new();
        let r = solveqp(
            inst, settings, &mut stats, &mut hms, &mut hb, &mut hs, &mut timer,
        );
        (r, hms, hb, hs)
    }

    /// min 0.5 x'Qx + c'x, Q = [[4,1,1],[1,3,0],[1,0,2]] stored full symmetric (CSC, as
    /// triangular_to_square_hessian builds it), c = -Q (1,2,1) = -(7,7,3), so the
    /// unconstrained minimiser is x = (1,2,1); -10 <= x <= 10, no rows. Q(1,2) = 0 but the
    /// Cholesky factor has fill-in there, which the 1.15.1 column-wise factorisation skips
    /// (it only visits stored entries), so 1.15.1 returns a wrong x on this input.
    fn pd3_bounded_qp() -> Instance {
        let mut inst = Instance::new(3, 0);
        inst.c = QpVector {
            num_nz: 3,
            dim: 3,
            index: vec![0, 1, 2],
            value: vec![-7.0, -7.0, -3.0],
        };
        inst.Q.mat = MatrixBase {
            num_row: 3,
            num_col: 3,
            start: vec![0, 3, 5, 7],
            index: vec![0, 1, 2, 0, 1, 0, 2],
            value: vec![4.0, 1.0, 1.0, 1.0, 3.0, 1.0, 2.0],
        };
        inst.A.mat = MatrixBase {
            num_row: 0,
            num_col: 3,
            start: vec![0, 0, 0, 0],
            index: vec![],
            value: vec![],
        };
        inst.con_lo = vec![];
        inst.con_up = vec![];
        inst.var_lo = vec![-10.0; 3];
        inst.var_up = vec![10.0; 3];
        inst
    }

    /// min x (c = 1, Q diagonal entry 0, so the regularised diagonal is 1e-7), x free.
    /// Unbounded below: the 1.16 test `res < -0.5 / reg` with `var_lo == -inf` catches it.
    fn free_1d_qp() -> Instance {
        let mut inst = Instance::new(1, 0);
        inst.c = QpVector {
            num_nz: 1,
            dim: 1,
            index: vec![0],
            value: vec![1.0],
        };
        inst.Q.mat = MatrixBase {
            num_row: 1,
            num_col: 1,
            start: vec![0, 1],
            index: vec![0],
            value: vec![0.0],
        };
        inst.A.mat = MatrixBase {
            num_row: 0,
            num_col: 1,
            start: vec![0, 0],
            index: vec![],
            value: vec![],
        };
        inst.con_lo = vec![];
        inst.con_up = vec![];
        inst.var_lo = vec![-kHighsInf];
        inst.var_up = vec![kHighsInf];
        inst
    }

    #[test]
    fn h16_qpbnd_diagonal_interior_is_optimal_with_zero_duals() {
        let mut inst = box_qp(false);
        let mut settings = Settings::default();
        settings.iteration_limit = 1000;
        settings.hessian_regularization_value = 0.0;
        settings.h16_qp_bounded_start = true;
        let (r, hms, hb, hs) = run_solveqp(&mut inst, &mut settings);
        assert_eq!(r, QpAsmStatus::kOk);
        assert_eq!(hms, HighsModelStatus::kOptimal);
        assert!(close(&hs.col_value, &[1.0, 1.0]));
        assert!(close(&hs.col_dual, &[0.0, 0.0]));
        assert!(hb.col_status.iter().all(|s| *s == HighsBasisStatus::kBasic));
    }

    #[test]
    fn h16_qpbnd_off_diagonal_interior_unchanged() {
        let mut inst = box_qp(false);
        let mut settings = Settings::default();
        settings.iteration_limit = 1000;
        settings.hessian_regularization_value = 0.0;
        let (r, hms, _, hs) = run_solveqp(&mut inst, &mut settings);
        assert_eq!(r, QpAsmStatus::kOk);
        assert_eq!(hms, HighsModelStatus::kOptimal);
        assert!(close(&hs.col_value, &[1.0, 1.0]));
    }

    #[test]
    fn h16_qpbnd_cholesky_matches_hand_solution() {
        let mut inst = pd3_bounded_qp();
        let mut settings = Settings::default();
        settings.hessian_regularization_value = 0.0;
        settings.h16_qp_bounded_start = true;
        let (r, hms, _, hs) = run_solveqp(&mut inst, &mut settings);
        assert_eq!(r, QpAsmStatus::kOk);
        assert_eq!(hms, HighsModelStatus::kOptimal);
        assert!(close(&hs.col_value, &[1.0, 2.0, 1.0]));
    }

    #[test]
    fn h16_qpbnd_off_cholesky_is_1_15_1() {
        // Same input with the flag off: the 1.15.1 Cholesky path. Hand trace of the
        // 1.15.1 C++: U = [[2, .5, .5], [0, sqrt(2.75), 0 (fill-in skipped)],
        // [0, 0, sqrt(2 - .25)]], so y = (3.5, 5.25/sqrt(2.75), 1.25/sqrt(1.75)) and
        // x = (337/308, 21/11, 5/7), not the minimiser (1, 2, 1) (#3153 'wholly incorrect').
        let mut inst = pd3_bounded_qp();
        let mut settings = Settings::default();
        settings.hessian_regularization_value = 0.0;
        let (r, hms, _, hs) = run_solveqp(&mut inst, &mut settings);
        assert_eq!(r, QpAsmStatus::kOk);
        assert_eq!(hms, HighsModelStatus::kOptimal);
        assert!(close(
            &hs.col_value,
            &[337.0 / 308.0, 21.0 / 11.0, 5.0 / 7.0]
        ));
    }

    #[test]
    fn h16_qpbnd_nonconvex_diagonal_is_solve_error() {
        let mut inst = box_qp(false);
        inst.c = QpVector {
            num_nz: 0,
            dim: 2,
            index: vec![],
            value: vec![0.0, 0.0],
        };
        inst.Q.mat.value = vec![1.0, -1.0];
        let mut settings = Settings::default();
        settings.h16_qp_bounded_start = true;
        let (r, hms, _, _) = run_solveqp(&mut inst, &mut settings);
        assert_eq!(r, QpAsmStatus::kError);
        assert_eq!(hms, HighsModelStatus::kSolveError);
    }

    #[test]
    fn h16_qpbnd_unbounded_free_var_sign_fix() {
        let mut inst = free_1d_qp();
        let mut settings = Settings::default();
        settings.h16_qp_bounded_start = true;
        let (r, hms, _, _) = run_solveqp(&mut inst, &mut settings);
        assert_eq!(r, QpAsmStatus::kOk);
        assert_eq!(hms, HighsModelStatus::kUnbounded);
    }

    #[test]
    fn h16_qpbnd_off_free_var_is_1_15_1() {
        let mut inst = free_1d_qp();
        let mut settings = Settings::default();
        let (r, hms, _, _) = run_solveqp(&mut inst, &mut settings);
        assert_eq!(r, QpAsmStatus::kOk);
        // 1.15.1 misses the unboundedness (`var_lo == +inf` typo): the point is interior.
        assert_eq!(hms, HighsModelStatus::kOptimal);
    }
}
