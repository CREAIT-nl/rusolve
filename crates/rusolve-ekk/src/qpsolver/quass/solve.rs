//! HiGHS source: qpsolver/quass.cpp (lines 309-569)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! See `mod.rs` for the ownership note. It also places the two
//! oracle trace calls (`qpsolver/oracle_trace.rs`, not upstream HiGHS):
//! `oracle_trace::qp_iter` at the end of every pass of the `while (true)`
//! loop and `oracle_trace::qp_end` after the post-loop `loginformation`,
//! at exactly the points of the OR1 C++ hook patch.

use super::*;
use crate::prelude::*;
use crate::qpsolver::snippets::indexof;
use crate::qpsolver::*;

impl Quass {
    /// HiGHS: quass.cpp::Quass::solve (lines 309-569)
    pub fn solve(
        &mut self,
        rt: &mut Runtime,
        x0: &QpVector,
        ra: &QpVector,
        b0: &mut Basis,
        timer: &mut HighsTimer,
    ) {
        // C++ line 311: feenableexcept dropped (debug-only).
        rt.statistics.time_start = Some(std::time::Instant::now());
        let basis = b0;
        rt.primal = x0.clone();

        // C++ lines 317-319: TODO comment dropped.

        let mut gradient = Gradient::new(rt);
        let mut redcosts = ReducedCosts::new(rt);
        let mut redgrad = ReducedGradient::new(rt);
        let mut factor = CholeskyFactor::new(rt, basis);
        rt.instance.A.mat_vec_into(&rt.primal, &mut rt.rowactivity);
        let mut pricing = get_pricing(rt, basis);

        let mut p = QpVector::new(rt.instance.num_var);
        let mut rowmove = QpVector::new(rt.instance.num_con);

        let mut buffer_yp = QpVector::new(rt.instance.num_var);
        let mut buffer_gyp = QpVector::new(rt.instance.num_var);
        let mut buffer_l = QpVector::new(rt.instance.num_var);
        let mut buffer_m = QpVector::new(rt.instance.num_var);

        let mut buffer_Qp = QpVector::new(rt.instance.num_var);

        // buffers for reduction
        let mut buffer_d = QpVector::new(rt.instance.num_var);

        regularize(rt);

        rt.relaxed_for_ratiotest = ratiotest_relax_instance(rt);

        let mut last_logging_iteration = rt.statistics.num_iterations - 1;
        let mut last_logging_time = 0.0;
        let mut logging_time_interval = 10.0;

        let current_num_active = basis.getnumactive();
        let mut atfsep = current_num_active == rt.instance.num_var;
        loop {
            // check iteration limit
            if rt.statistics.num_iterations >= rt.settings.iteration_limit {
                rt.status = QpModelStatus::kIterationLimit;
                break;
            }

            // check time limit
            if timer.read() >= rt.settings.time_limit {
                rt.status = QpModelStatus::kTimeLimit;
                break;
            }

            if basis.getnuminactive() > rt.settings.nullspace_limit {
                rt.settings
                    .nullspace_limit_log
                    .fire(&mut rt.settings.nullspace_limit);
                rt.status = QpModelStatus::kLargeNullspace;
                return;
            }

            // LOGGING
            let run_time = timer.read();
            if (rt.statistics.num_iterations % rt.settings.reportingfequency == 0
                || run_time - last_logging_time > logging_time_interval)
                && rt.statistics.num_iterations > last_logging_iteration
            {
                let mut log_report = true;
                if rt.statistics.num_iterations > 10 * rt.settings.reportingfequency {
                    rt.settings.reportingfequency *= 10;
                    log_report = false;
                }
                if run_time > 10.0 * logging_time_interval {
                    logging_time_interval *= 2.0;
                }
                if log_report {
                    last_logging_time = run_time;
                    last_logging_iteration = rt.statistics.num_iterations;
                    loginformation(rt, basis, &mut factor, timer);
                    rt.settings.iteration_log.fire(&mut rt.statistics);
                }
            }

            // REINVERSION
            let mut status = QpSolverStatus::OK;

            if check_reinvert_due(basis) {
                status = reinvert(
                    rt,
                    basis,
                    &mut factor,
                    &mut gradient,
                    &mut redcosts,
                    &mut redgrad,
                    &mut pricing,
                );
                if status != QpSolverStatus::OK {
                    // C++ lines 395-402: notOkReturn lambda inlined.
                    if status == QpSolverStatus::NOTPOSITIVDEFINITE {
                        rt.status = QpModelStatus::kNonConvex;
                    } else {
                        rt.status = QpModelStatus::kError;
                    }
                    return;
                }
            }

            let mut zero_curvature_direction = false;
            let mut maxsteplength = 1.0;
            // Oracle trace bookkeeping only (not in C++): the values the oracle hook
            // patch reports for this pass.
            let mut trace_minidx: HighsInt = -1;
            let mut trace_limiting: HighsInt = -2;
            let mut trace_alpha: f64 = 0.0;
            if atfsep {
                // Determine a variable to relax from being active. If there is
                // none, then basis is optimal
                // C++ evaluates the `gradient.getGradient()` argument (which may
                // recompute) before `price` runs (pricing.rs PORT NOTE).
                gradient.get_gradient(rt);
                let minidx = pricing.price(rt, basis, &mut redcosts, &mut gradient, &rt.primal);
                if minidx == -1 {
                    rt.status = QpModelStatus::kOptimal;
                    break;
                }
                trace_minidx = minidx;
                // Now perform a real iteration
                rt.statistics.num_iterations += 1;

                let unit = basis.getindexinfactor()[minidx as usize];
                QpVector::unit_into(rt.instance.num_var, unit, &mut buffer_yp);
                // C++ `basis.btran(buffer_yp, buffer_yp, ...)` aliases in/out; btran copies rhs
                // into its hvec first, so a clone of the input is the same computation.
                let buffer_yp_rhs = buffer_yp.clone();
                basis.btran_into(&buffer_yp_rhs, &mut buffer_yp, true, minidx);

                buffer_l.dim = basis.getnuminactive();
                buffer_m.dim = basis.getnuminactive();
                computesearchdirection_major(
                    rt,
                    basis,
                    &mut factor,
                    &buffer_yp,
                    &mut gradient,
                    &mut buffer_gyp,
                    &mut buffer_l,
                    &mut buffer_m,
                    &mut p,
                    &mut status,
                );
                if status != QpSolverStatus::OK {
                    if status == QpSolverStatus::NOTPOSITIVDEFINITE {
                        rt.status = QpModelStatus::kNonConvex;
                    } else {
                        rt.status = QpModelStatus::kError;
                    }
                    return;
                }
                basis.deactivate(minidx);
                computerowmove(rt, basis, &mut p, &mut rowmove);
                tidyup(&mut p, &mut rowmove, basis, rt);
                // C++ line 434 (`maxsteplength = infinity`) is overwritten by line 436 before any read; dropped.
                // C++ line 435: commented-out if dropped.
                maxsteplength = computemaxsteplength(
                    rt,
                    &p,
                    &mut gradient,
                    &mut buffer_Qp,
                    &mut zero_curvature_direction,
                );
                if !zero_curvature_direction {
                    // HiGHS 1.16 #3325 (h16_qp_convexity): quass.cpp lines 459-460 (1.16.0);
                    // the dot is computed only on the on path.
                    let curvature = if rt.settings.h16_qp_convexity {
                        p.dot(&buffer_Qp)
                    } else {
                        0.0
                    };
                    status = factor.expand(
                        &buffer_yp,
                        &mut buffer_gyp,
                        &mut buffer_l,
                        &mut buffer_m,
                        curvature,
                    );
                    if status != QpSolverStatus::OK {
                        if status == QpSolverStatus::NOTPOSITIVDEFINITE {
                            rt.status = QpModelStatus::kNonConvex;
                        } else {
                            rt.status = QpModelStatus::kError;
                        }
                        return;
                    }
                }
                redgrad.expand(rt, &mut gradient, &buffer_yp);
            } else {
                // Compute a search direction - which may be zero, in which case
                // atfsep is set true and the loop repeats with this
                // condition. In particular, this happens when the current basis
                // is optimal
                computesearchdirection_minor(
                    rt,
                    basis,
                    &mut factor,
                    &mut redgrad,
                    &mut gradient,
                    &mut p,
                    &mut status,
                );
                if status != QpSolverStatus::OK {
                    if status == QpSolverStatus::NOTPOSITIVDEFINITE {
                        rt.status = QpModelStatus::kNonConvex;
                    } else {
                        rt.status = QpModelStatus::kError;
                    }
                    return;
                }
                computerowmove(rt, basis, &mut p, &mut rowmove);
                tidyup(&mut p, &mut rowmove, basis, rt);
                rt.instance.Q.mat_vec_into(&p, &mut buffer_Qp);
            }
            if p.norm2() < rt.settings.pnorm_zero_threshold
                || maxsteplength == 0.0
                || (false
                    && gradient.get_gradient(rt).dot(&p).abs()
                        < rt.settings.improvement_zero_threshold)
            {
                atfsep = true;
            } else {
                // Now perform a real iteration
                rt.statistics.num_iterations += 1;

                let stepres = ratiotest(rt, &p, &rowmove, maxsteplength);
                trace_limiting = stepres.limitingconstraint;
                trace_alpha = stepres.alpha;
                if stepres.limitingconstraint != -1 {
                    let mut constrainttodrop: HighsInt = 0;
                    let mut maxabsd: HighsInt = 0;
                    // C++ lines 467-471: comment about #2734 dropped.
                    const kIllegalLogDValue: f64 = kHighsInf;
                    let mut log_d = kIllegalLogDValue;
                    status = reduce(
                        rt,
                        basis,
                        stepres.limitingconstraint,
                        &mut buffer_d,
                        &mut maxabsd,
                        &mut constrainttodrop,
                        &mut log_d,
                    );
                    if status != QpSolverStatus::OK {
                        // Possibly perform degeneracy failure logging
                        if status == QpSolverStatus::DEGENERATE {
                            // Ensure that log_d has been set
                            debug_assert!(log_d != kIllegalLogDValue);
                            let mut degeneracy_fail_data = (maxabsd, log_d);
                            rt.settings
                                .degeneracy_fail_log
                                .fire(&mut degeneracy_fail_data);
                        }
                        rt.status = QpModelStatus::kUndetermined;
                        return;
                    }
                    if !zero_curvature_direction {
                        factor.reduce(
                            &buffer_d,
                            maxabsd,
                            indexof(basis.getinactive(), &stepres.limitingconstraint) != -1,
                        );
                    } else if rt.settings.h16_qp_convexity {
                        // HiGHS 1.16 #3325 (h16_qp_convexity): quass.cpp lines 513-517 (1.16.0)
                        factor.invalidate();
                    }
                    redgrad.reduce(&buffer_d, maxabsd);
                    redgrad.update(stepres.alpha, false);

                    status = basis.activate(
                        rt,
                        stepres.limitingconstraint,
                        if stepres.nowactiveatlower {
                            BasisStatus::kActiveAtLower
                        } else {
                            BasisStatus::kActiveAtUpper
                        },
                        constrainttodrop,
                        pricing.as_mut(),
                    );
                    if status != QpSolverStatus::OK {
                        rt.status = QpModelStatus::kUndetermined;
                        return;
                    }
                    if basis.getnumactive() != rt.instance.num_var {
                        atfsep = false;
                    }
                } else {
                    if stepres.alpha == f64::INFINITY {
                        // unbounded
                        rt.status = QpModelStatus::kUnbounded;
                        return;
                    }
                    atfsep = false;
                    redgrad.update(stepres.alpha, false);
                }

                rt.primal.saxpy(stepres.alpha, &p);
                rt.rowactivity.saxpy(stepres.alpha, &rowmove);

                gradient.update(&buffer_Qp, stepres.alpha);
                redcosts.update();
            }

            // Oracle trace: qp_iter at the end of every pass (C++ hook patch: after
            // redcosts.update() on a step, after atfsep = true on a null step).
            if oracle_trace::enabled() {
                // C++ `gradient.oracleGradient()` reads the stored gradient without
                // the recompute `getGradient()` may trigger. gradient.rs has
                // no such accessor yet, so the read goes through a clone: the solver
                // state is untouched, but on a pass where the stored gradient is stale
                // (numupdates >= gradientrecomputefrequency) h_gradient is the
                // recomputed vector. FOLLOW-UP: add `Gradient::oracle_gradient`.
                let mut trace_gradient = gradient.clone();
                let num_active = basis.getnumactive() as i64;
                let num_inactive = basis.getnuminactive() as i64;
                oracle_trace::qp_iter(
                    rt.statistics.num_iterations as i64,
                    atfsep,
                    trace_minidx as i64,
                    trace_limiting as i64,
                    trace_alpha,
                    num_active,
                    num_inactive,
                    &rt.primal,
                    &rt.rowactivity,
                    trace_gradient.get_gradient(rt),
                    &p,
                );
            }
        }

        loginformation(rt, basis, &mut factor, timer);
        // Oracle trace: qp_end between loginformation and iteration_log.fire, as in
        // the oracle hook patch.
        oracle_trace::qp_end(
            rt.statistics.num_iterations as i64,
            rt.status as i64,
            &rt.primal,
            &rt.rowactivity,
        );
        rt.settings.iteration_log.fire(&mut rt.statistics);
        // C++ line 528: basis.report() commented out, dropped.

        let rowact = rt.instance.A.mat_vec(&rt.primal);
        rt.instance.sumnumprimalinfeasibilities(&rt.primal, &rowact);

        let lambda = redcosts.get_reduced_costs(rt, basis, &mut gradient);
        let active = basis.getactive().clone();
        let indexinfactor = basis.getindexinfactor().clone();
        for &e in &active {
            let indexinbasis = indexinfactor[e as usize];
            if e >= rt.instance.num_con {
                // active variable bound
                let var = e - rt.instance.num_con;
                rt.dualvar.value[var as usize] = lambda.value[indexinbasis as usize];
            } else {
                rt.dualcon.value[e as usize] = lambda.value[indexinbasis as usize];
            }
        }
        rt.dualcon.resparsify();
        rt.dualvar.resparsify();

        // C++ lines 547-553: commented-out debug code dropped.

        // extract basis status
        for i in 0..rt.instance.num_var {
            rt.status_var[i as usize] = basis.getstatus(rt.instance.num_con + i);
        }

        for i in 0..rt.instance.num_con {
            rt.status_con[i as usize] = basis.getstatus(i);
        }

        if basis.getnumactive() == rt.instance.num_var {
            rt.primal = basis.recomputex(&rt.instance);
        }
        // C++ line 567: x.report("x") commented out, dropped.
        rt.statistics.time_end = Some(std::time::Instant::now());
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
    /// `num_con + i`); the row (id 0) is not in the basis. C++ `Basis(rt,
    /// active, status, inactive)` reads `status[i]` for `active[i]` (basis.cpp
    /// 13-42), so `status` is parallel to `active`, not indexed by constraint id.
    fn at_lower(rt: &mut Runtime) -> Basis {
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
    fn solves_box_qp_from_lower_bounds() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        // C++ settings.hpp:76 `HighsInt iteration_limit = numeric_limits<HighsInt>::infinity()`
        // is 0 (integer numeric_limits::infinity() returns 0), so the C++ default stops at
        // once num_iterations >= 0. This test wants a real solve, so it sets a finite limit.
        rt.settings.iteration_limit = 1000;
        let mut basis = at_lower(&mut rt);
        let x0 = QpVector::new(2);
        let ra = QpVector::new(1);
        let mut timer = HighsTimer::new();
        let mut quass = Quass::new(&rt);
        super::super::oracle_trace::capture_begin();
        quass.solve(&mut rt, &x0, &ra, &mut basis, &mut timer);
        let trace = super::super::oracle_trace::capture_end();
        assert_eq!(rt.status, QpModelStatus::kOptimal);
        assert!(close(&rt.primal.value, &[1.0, 1.0]));
        assert!(rt.statistics.num_iterations > 0);
        assert!(trace.lines().any(|l| l.starts_with("{\"t\":\"qp_iter\"")));
        assert!(
            trace
                .lines()
                .last()
                .unwrap()
                .starts_with("{\"t\":\"qp_end\"")
        );
    }

    /// Column-wise LP builder for the PR #3325 model: `cols[j] = (cost, lower,
    /// upper, [(row, value)])`.
    #[allow(clippy::type_complexity)]
    fn lp_from_cols(cols: &[(f64, f64, f64, &[(HighsInt, f64)])], rows: &[(f64, f64)]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = cols.len() as HighsInt;
        lp.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = vec![0];
        for &(cost, lower, upper, entries) in cols {
            lp.col_cost_.push(cost);
            lp.col_lower_.push(lower);
            lp.col_upper_.push(upper);
            for &(r, v) in entries {
                lp.a_matrix_.index_.push(r);
                lp.a_matrix_.value_.push(v);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        for &(l, u) in rows {
            lp.row_lower_.push(l);
            lp.row_upper_.push(u);
        }
        lp
    }

    /// check/instances/3322.mps (PR #3325, TestQpSolver issue-3322), built by
    /// hand: MAX, 15 columns in file order, 3 E rows (norm, PolyF norm,
    /// coupling). The only quadratic entry is z,z = -11.390625. Rows 0 and 1
    /// have rhs 1 and 0; row 2 has rhs 0. Fixed columns are r, omega, gamma.
    fn issue_3322() -> (HighsLp, HighsHessian) {
        let inf = kHighsInf;
        let cols: &[(f64, f64, f64, &[(HighsInt, f64)])] = &[
            // MPB_lambda
            (0.0, 0.0, inf, &[(0, 1.0), (1, -1.0)]),
            // MPB_r (fixed at 0)
            (0.0, 0.0, 0.0, &[(0, 1.0)]),
            // MPB_omega (fixed at 0)
            (0.0, 0.0, 0.0, &[(0, -1.0)]),
            // MPB_z (free)
            (0.0, -inf, inf, &[(2, 1.0)]),
            // MPB_s_plus
            (-57.1072724509965, 0.0, inf, &[(2, 1.0)]),
            // MPB_s_minus
            (-999942.892727548, 0.0, inf, &[(2, -1.0)]),
            // PolyF_gamma (fixed at 0)
            (0.0, 0.0, 0.0, &[(1, 1.0)]),
            // c_ekk0..c_ekk7
            (
                -53905.2151686991,
                0.0,
                inf,
                &[(1, 2.05180190919353e-05), (2, 1797.35634311114)],
            ),
            (
                8.37729885461158e-12,
                0.0,
                inf,
                &[(1, 2.05181177843192e-05), (2, -0.00483858476063246)],
            ),
            (
                -26325.1652341927,
                0.0,
                inf,
                &[(1, 3.42995534290902e-05), (2, 1377.10573369115)],
            ),
            (
                -22199.0995786157,
                0.0,
                inf,
                &[(1, 0.000132194621708021), (2, 2412.49746327259)],
            ),
            (
                -10343.8913476677,
                0.0,
                inf,
                &[(1, 0.00030220365038542), (2, 3309.02687219242)],
            ),
            (
                -1159.73839015062,
                0.0,
                inf,
                &[(1, 0.000370124363750065), (2, 2701.79458025431)],
            ),
            (
                -0.0825272315620604,
                0.0,
                inf,
                &[(1, 0.000347604032739213), (2, 1268.66054205703)],
            ),
            (
                0.0,
                0.0,
                inf,
                &[(1, 0.000799020940383498), (2, -0.188425204796894)],
            ),
        ];
        let rows = [(1.0, 1.0), (0.0, 0.0), (0.0, 0.0)];
        let mut lp = lp_from_cols(cols, &rows);
        lp.sense_ = ObjSense::kMaximize;
        let hessian = HighsHessian {
            dim_: 15,
            format_: HessianFormat::kTriangular,
            start_: (0..=15).map(|j| if j <= 3 { 0 } else { 1 }).collect(),
            index_: vec![3],
            value_: vec![-11.390625],
        };
        (lp, hessian)
    }

    fn run_issue_3322(flag: bool, regularization: f64) -> (HighsStatus, HighsModelStatus, f64) {
        let (lp, hessian) = issue_3322();
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_qp_convexity", flag),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_double("qp_regularization_value", regularization),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        assert_eq!(highs.pass_hessian(hessian), HighsStatus::kOk);
        let r = highs.run();
        (
            r,
            highs.get_model_status(),
            highs.get_info().objective_function_value,
        )
    }

    /// HiGHS 1.16 #3325 on path: the issue-3322 QP is convex and solves to the PR's
    /// objective, for both regularisations of the PR test.
    #[test]
    fn h16_qpconv_issue_3322_on_is_optimal() {
        for regularization in [1e-7, 0.0] {
            let (r, hms, obj) = run_issue_3322(true, regularization);
            assert_eq!(r, HighsStatus::kOk);
            assert_eq!(hms, HighsModelStatus::kOptimal);
            let expected = -1.5339266953e-02;
            assert!((obj - expected).abs() / f64::max(1.0, expected.abs()) < 1e-5);
        }
    }

    /// HiGHS 1.16 #3325 off path: the same model on the 1.15.1 path, observed on this
    /// port. With regularisation 1e-7 the run returns kError with the model
    /// status kNotset. With regularisation 0 it returns kOptimal at objective
    /// -736.2225 instead of -0.0153, i.e. the 1.15.1 path is wrong here. The
    /// PR's fix is what the on-path test above checks.
    #[test]
    fn h16_qpconv_off_issue_3322_is_1_15_1() {
        let (r, hms, _) = run_issue_3322(false, 1e-7);
        assert_eq!(r, HighsStatus::kError);
        assert_eq!(hms, HighsModelStatus::kNotset);
        let (r, hms, obj) = run_issue_3322(false, 0.0);
        assert_eq!(r, HighsStatus::kOk);
        assert_eq!(hms, HighsModelStatus::kOptimal);
        assert!((obj - -736.2225041433807).abs() <= 1e-6 * 736.0);
    }
}
