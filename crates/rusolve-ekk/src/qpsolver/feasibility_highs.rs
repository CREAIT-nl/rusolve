//! HiGHS source: qpsolver/feasibility_highs.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Builds a nested `Highs` for the phase-1 LP (C++ `Highs highs;`
//! with `output_flag=false`, `presolve=kHighsOnString`, `time_limit`, and, on the
//! hot-start path, `setBasis`); use `Highs::default()` plus the Rust
//! `set_option_value_*` / `pass_model` / `run` / getters already in the crate.

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: feasibility_highs.hpp::computeStartingPointHighs (lines 14-302)
pub fn compute_starting_point_highs(
    instance: &mut Instance,
    settings: &mut Settings,
    stats: &mut Statistics,
    modelstatus: &mut QpModelStatus,
    result: &mut QpHotstartInformation,
    highs_model_status: &mut HighsModelStatus,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
    timer: &mut HighsTimer,
) {
    let mut have_starting_point = false;
    // C++ lines 41-50, 186-189, 264-272: printf under const bool debug_report = false; dropped.
    if settings.allow_hot_start && highs_solution.value_valid {
        // #1350 add primal_feasibility_tolerance to settings
        let primal_feasibility_tolerance = settings.lambda_zero_threshold;

        let mut num_var_infeasibilities: HighsInt = 0;
        let mut max_var_infeasibility: f64 = 0.0;
        let mut sum_var_infeasibilities: f64 = 0.0;
        let mut num_con_infeasibilities: HighsInt = 0;
        let mut max_con_infeasibility: f64 = 0.0;
        let mut sum_con_infeasibilities: f64 = 0.0;
        let mut max_con_residual: f64 = 0.0;
        let mut sum_con_residuals: f64 = 0.0;

        assess_qp_primal_feasibility(
            instance,
            primal_feasibility_tolerance,
            &highs_solution.col_value,
            &highs_solution.row_value,
            &mut num_var_infeasibilities,
            &mut max_var_infeasibility,
            &mut sum_var_infeasibilities,
            &mut num_con_infeasibilities,
            &mut max_con_infeasibility,
            &mut sum_con_infeasibilities,
            &mut max_con_residual,
            &mut sum_con_residuals,
        );

        have_starting_point =
            num_var_infeasibilities == 0 && num_con_infeasibilities == 0 && highs_basis.valid;
    }
    // compute initial feasible point
    let mut use_basis: HighsBasis;
    let mut use_solution: HighsSolution;
    let isfreevar = |i_var: HighsInt| -> bool {
        instance.var_lo[i_var as usize] == -kHighsInf
            && instance.var_up[i_var as usize] == kHighsInf
    };
    if have_starting_point {
        use_basis = highs_basis.clone();
        use_solution = highs_solution.clone();
        // Have to assume that the supplied basis is feasible
        *modelstatus = QpModelStatus::kNotset;
    } else {
        let mut highs = Highs::default();

        // set HiGHS to be silent
        highs.set_option_value_bool("output_flag", false);
        highs.set_option_value_str("presolve", kHighsOnString);
        // Set the residual time limit
        let use_time_limit = cpp_max(settings.time_limit - timer.read(), 0.001);
        highs.set_option_value_double("time_limit", use_time_limit);

        let mut lp = HighsLp::default();
        lp.a_matrix_.index_ = instance.A.mat.index.clone();
        lp.a_matrix_.start_ = instance.A.mat.start.clone();
        lp.a_matrix_.value_ = instance.A.mat.value.clone();
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.col_cost_ = vec![0.0; instance.num_var as usize];
        // C++ line 83: lp.col_cost_ = runtime.instance.c.value; commented out.
        lp.col_lower_ = instance.var_lo.clone();
        lp.col_upper_ = instance.var_up.clone();
        lp.row_lower_ = instance.con_lo.clone();
        lp.row_upper_ = instance.con_up.clone();
        lp.num_col_ = instance.num_var;
        lp.num_row_ = instance.num_con;

        // create artificial bounds for free variables: false by default
        debug_assert!(!settings.phase1boundfreevars);
        if settings.phase1boundfreevars {
            for i in 0..instance.num_var {
                if isfreevar(i) {
                    lp.col_lower_[i as usize] = -1E5;
                    lp.col_upper_[i as usize] = 1E5;
                }
            }
        }

        highs.pass_model(lp);
        // Make free variables basic: false by default
        debug_assert!(!settings.phase1movefreevarsbasic);
        if settings.phase1movefreevarsbasic {
            let mut basis = HighsBasis::default();
            basis.alien = true; // Set true when basis is instantiated
            for _i in 0..instance.num_con {
                basis.row_status.push(HighsBasisStatus::kNonbasic);
            }

            for i in 0..instance.num_var {
                // make free variables basic
                if instance.var_lo[i as usize] == -kHighsInf
                    && instance.var_up[i as usize] == kHighsInf
                {
                    // free variable
                    basis.col_status.push(HighsBasisStatus::kBasic);
                } else {
                    basis.col_status.push(HighsBasisStatus::kNonbasic);
                }
            }

            highs.set_basis(&basis, "");

            highs.set_option_value_int("simplex_strategy", kSimplexStrategyPrimal);
        }

        let status = highs.run();
        if status == HighsStatus::kError {
            *modelstatus = QpModelStatus::kError;
            return;
        }

        let phase1stat = highs.get_model_status();
        match phase1stat {
            HighsModelStatus::kOptimal => {
                *modelstatus = QpModelStatus::kNotset;
            }
            HighsModelStatus::kInfeasible => {
                *modelstatus = QpModelStatus::kInfeasible;
            }
            HighsModelStatus::kTimeLimit => {
                *modelstatus = QpModelStatus::kTimeLimit;
            }
            HighsModelStatus::kInterrupt => {
                *modelstatus = QpModelStatus::kInterrupt;
            }
            _ => {
                *modelstatus = QpModelStatus::kError;
            }
        }

        stats.phase1_iterations = highs.get_info().simplex_iteration_count;

        if *modelstatus != QpModelStatus::kNotset {
            return;
        }

        // Should only get here if feasibility problem is solved to
        // optimality - hence there is a feasible basis
        debug_assert!(phase1stat == HighsModelStatus::kOptimal);

        use_basis = highs.get_basis().clone();
        use_solution = highs.get_solution().clone();
    }

    let mut num_small_x0: HighsInt = 0;
    let mut num_small_ra: HighsInt = 0;
    let zero_activity_tolerance: f64 = if have_starting_point { 0.0 } else { 1e-4 };
    let mut x0 = QpVector::new(instance.num_var);
    let mut ra = QpVector::new(instance.num_con);
    for i in 0..x0.dim {
        if use_solution.col_value[i as usize].abs() > zero_activity_tolerance {
            x0.value[i as usize] = use_solution.col_value[i as usize];
            x0.index[x0.num_nz as usize] = i;
            x0.num_nz += 1;
        } else if use_solution.col_value[i as usize].abs() > 0.0 {
            num_small_x0 += 1;
        }
    }

    for i in 0..ra.dim {
        if use_solution.row_value[i as usize].abs() > zero_activity_tolerance {
            ra.value[i as usize] = use_solution.row_value[i as usize];
            ra.index[ra.num_nz as usize] = i;
            ra.num_nz += 1;
        } else if use_solution.row_value[i as usize].abs() > 0.0 {
            num_small_ra += 1;
        }
    }
    // C++ lines 186-189: printf under debug_report dropped; the counters are only read there.
    let _ = (num_small_x0, num_small_ra);
    let mut initial_active: Vec<HighsInt> = Vec::new();
    let mut initial_inactive: Vec<HighsInt> = Vec::new();
    let mut initial_status: Vec<BasisStatus> = Vec::new();

    let num_highs_basis_status: HighsInt = (HighsBasisStatus::kNonbasic as i32) + 1;
    let mut debug_row_status_count: Vec<HighsInt> = vec![0; num_highs_basis_status as usize];
    for i in 0..(use_basis.row_status.len() as HighsInt) {
        let status = use_basis.row_status[i as usize];
        debug_row_status_count[status as i32 as usize] += 1;
        if status == HighsBasisStatus::kLower {
            initial_active.push(i);
            initial_status.push(BasisStatus::kActiveAtLower);
        } else if status == HighsBasisStatus::kUpper {
            initial_active.push(i);
            initial_status.push(BasisStatus::kActiveAtUpper);
        } else if status == HighsBasisStatus::kZero {
            // Shouldn't happen, since free rows are basic in a logical
            // basis and remain basic, or are removed by presolve and
            // restored as basic in postsolve
            // C++ `assert(111 == 222)`: always false; debug_assert_eq! avoids clippy's assertions_on_constants.
            debug_assert_eq!(111, 222);
            // That said, a free row that is nonbasic in the Highs basis
            // must be counted as inactive in the QP basis for accounting
            // purposes
            initial_inactive.push(i);
        } else if status != HighsBasisStatus::kBasic {
            debug_assert!(status == HighsBasisStatus::kNonbasic);
            // Surely an error, but not a problem before, since simplex
            // solver cannot return a HighsBasisStatus::kNonbasic
            // variable. Does matter now, since a saved QP basis will
            // generally have such variables.
            //
            //      initial_inactive.push_back(instance.num_con + i);
            //
            // A HighsBasisStatus::kNonbasic variable corresponds one-to-one
            // with being inactive in the QP basis
            initial_inactive.push(i);
        } else {
            debug_assert!(status == HighsBasisStatus::kBasic);
        }
    }

    let mut debug_col_status_count: Vec<HighsInt> = vec![0; num_highs_basis_status as usize];
    for i in 0..(use_basis.col_status.len() as HighsInt) {
        let status = use_basis.col_status[i as usize];
        debug_col_status_count[status as i32 as usize] += 1;
        if status == HighsBasisStatus::kLower {
            if isfreevar(i) {
                initial_inactive.push(instance.num_con + i);
            } else {
                initial_active.push(instance.num_con + i);
                initial_status.push(BasisStatus::kActiveAtLower);
            }
        } else if status == HighsBasisStatus::kUpper {
            if isfreevar(i) {
                initial_inactive.push(instance.num_con + i);
            } else {
                initial_active.push(instance.num_con + i);
                initial_status.push(BasisStatus::kActiveAtUpper);
            }
        } else if status == HighsBasisStatus::kZero {
            initial_inactive.push(instance.num_con + i);
        } else if status != HighsBasisStatus::kBasic {
            debug_assert!(status == HighsBasisStatus::kNonbasic);
            initial_inactive.push(instance.num_con + i);
        } else {
            debug_assert!(status == HighsBasisStatus::kBasic);
        }
    }

    // C++ lines 264-272: printf under debug_report dropped.

    // This used to be an assert
    if ((initial_active.len() + initial_inactive.len()) as HighsInt) != instance.num_var {
        *modelstatus = QpModelStatus::kError;
        return;
    }

    if !have_starting_point {
        // When starting from a feasible basis, there will generally be
        // inactive variables in the basis that aren't free
        for ia in initial_inactive.iter() {
            if *ia < instance.num_con {
                // printf("free row %d\n", (int)ia);
                debug_assert!(instance.con_lo[*ia as usize] == -kHighsInf);
                debug_assert!(instance.con_up[*ia as usize] == kHighsInf);
            } else {
                // printf("free col %d\n", (int)ia);
                debug_assert!(instance.var_lo[(*ia - instance.num_con) as usize] == -kHighsInf);
                debug_assert!(instance.var_up[(*ia - instance.num_con) as usize] == kHighsInf);
            }
        }
    }

    result.status = initial_status;
    result.active = initial_active;
    result.inactive = initial_inactive;
    result.primal = x0;
    result.rowact = ra;
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
    fn highs_start_is_feasible() {
        let mut inst = box_qp(true);
        let mut settings = Settings::default();
        let mut stats = Statistics::default();
        let mut status = QpModelStatus::kNotset;
        let mut result = QpHotstartInformation::new(2, 1);
        let mut hms = HighsModelStatus::kNotset;
        let mut hb = HighsBasis::default();
        let mut hs = HighsSolution::default();
        let mut timer = HighsTimer::new();
        compute_starting_point_highs(
            &mut inst,
            &mut settings,
            &mut stats,
            &mut status,
            &mut result,
            &mut hms,
            &mut hb,
            &mut hs,
            &mut timer,
        );
        assert_ne!(status, QpModelStatus::kInfeasible);
        let (x, ra) = (&result.primal.value, &result.rowact.value);
        assert_eq!(ra[0], x[0] + x[1]);
        assert!(ra[0] <= 10.0 && x.iter().all(|v| (0.0..=10.0).contains(v)));
        assert_eq!(result.active.len() + result.inactive.len(), 2);
    }
}
