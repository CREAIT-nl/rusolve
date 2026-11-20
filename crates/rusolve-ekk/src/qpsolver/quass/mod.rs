//! HiGHS source: qpsolver/quass.hpp, qpsolver/quass.cpp (lines 31-307; see `solve.rs` for 309-569)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//!
//!
//! PORT NOTE (ownership): C++ `Quass` stores only `Runtime& runtime`; here it
//! is a unit struct and `Quass::solve` takes `rt: &mut Runtime`. The static
//! helpers keep their C++ parameter lists, with `Runtime&` -> `&Runtime`
//! where only read and `&mut Runtime` where written; helpers that reach
//! `Gradient` / `ReducedCosts` / `CholeskyFactor` methods that need `rt` or
//! the basis get them as extra parameters (see each object's module doc).
//! `std::unique_ptr<Pricing>` -> `Box<dyn Pricing>`.

pub mod solve;

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: quass.hpp::Quass (lines 17-25)
#[derive(Clone, Copy, Debug, Default)]
pub struct Quass;

impl Quass {
    /// HiGHS: quass.cpp::Quass::Quass (line 31)
    pub fn new(rt: &Runtime) -> Self {
        Quass
    }
}

/// HiGHS: quass.cpp::loginformation (lines 33-46)
pub(crate) fn loginformation(
    rt: &mut Runtime,
    basis: &mut Basis,
    factor: &mut CholeskyFactor,
    timer: &mut HighsTimer,
) {
    rt.statistics.iteration.push(rt.statistics.num_iterations);
    rt.statistics
        .nullspacedimension
        .push(rt.instance.num_var - basis.getnumactive());
    rt.statistics.objval.push(rt.instance.objval(&rt.primal));
    rt.statistics.time.push(timer.read());
    let sm = rt
        .instance
        .sumnumprimalinfeasibilities(&rt.primal, &rt.rowactivity);
    rt.statistics.sum_primal_infeasibilities.push(sm.sum);
    rt.statistics.num_primal_infeasibilities.push(sm.num);
    rt.statistics.density_factor.push(factor.density());
    rt.statistics.density_nullspace.push(0.0);
}

/// HiGHS: quass.cpp::tidyup (lines 48-57)
///
/// NOTE: the C++ loop variable is `unsigned acon` over `basis.getactive()`.
pub(crate) fn tidyup(
    p: &mut QpVector,
    rowmove: &mut QpVector,
    basis: &mut Basis,
    runtime: &Runtime,
) {
    for &acon in basis.getactive() {
        if (acon as HighsInt) >= runtime.instance.num_con {
            p.value[(acon - runtime.instance.num_con) as usize] = 0.0;
        } else {
            rowmove.value[acon as usize] = 0.0;
        }
    }
}

/// HiGHS: quass.cpp::computerowmove (lines 59-84)
pub(crate) fn computerowmove(
    runtime: &Runtime,
    basis: &mut Basis,
    p: &mut QpVector,
    rowmove: &mut QpVector,
) {
    runtime.instance.A.mat_vec_into(p, rowmove);
    // C++ lines 63-83: unreachable code after return, dropped.
}

/// HiGHS: quass.cpp::computesearchdirection_minor (lines 87-101)
///
/// Returns nothing: C++ returns the `p` reference it was given.
pub(crate) fn computesearchdirection_minor(
    rt: &Runtime,
    bas: &mut Basis,
    cf: &mut CholeskyFactor,
    redgrad: &mut ReducedGradient,
    gradient: &mut Gradient,
    p: &mut QpVector,
    status: &mut QpSolverStatus,
) {
    *status = QpSolverStatus::OK;
    let mut g2 = redgrad.get(rt, bas, gradient).op_neg();
    g2.sanitize(1e-14);
    *status = cf.solve(rt, bas, &mut g2);
    if *status != QpSolverStatus::OK {
        bas.Zprod(&g2, p);
        return;
    }
    g2.sanitize(1e-14);

    bas.Zprod(&g2, p);
}

/// HiGHS 1.16 (#3325): quass.cpp::computesearchdirection_major (lines 103-138)
/// HiGHS: quass.cpp::computesearchdirection_major (lines 104-134)
///
/// Returns nothing: C++ returns the `p` reference it was given.
pub(crate) fn computesearchdirection_major(
    runtime: &Runtime,
    basis: &mut Basis,
    factor: &mut CholeskyFactor,
    yp: &QpVector,
    gradient: &mut Gradient,
    gyp: &mut QpVector,
    l: &mut QpVector,
    m: &mut QpVector,
    p: &mut QpVector,
    status: &mut QpSolverStatus,
) {
    *status = QpSolverStatus::OK;
    let mut yyp = yp.clone();
    // C++ lines 110-112: commented-out code, dropped.
    runtime.instance.Q.mat_vec_into(&yyp, gyp);
    if basis.getnumactive() < runtime.instance.num_var {
        basis.Ztprod(gyp, m, false, -1);
        *l = m.clone();
        *status = factor.solve_l(runtime, basis, l);
        if *status != QpSolverStatus::OK {
            p.saxpy_ab(-1.0, 1.0, &yyp);
            return;
        }
        let mut v = l.clone();
        factor.solve_lt(&mut v);
        // HiGHS 1.16 #3325 (h16_qp_convexity): quass.cpp line 124 (1.16.0)
        if runtime.settings.h16_qp_convexity {
            v.resparsify();
        }
        basis.Zprod(&v, p);
        if gradient.get_gradient(runtime).dot(&yyp) < 0.0 {
            p.saxpy_ab(-1.0, 1.0, &yyp);
        } else if runtime.settings.h16_qp_convexity {
            // HiGHS 1.16 #3325 (h16_qp_convexity): quass.cpp line 131 (1.16.0)
            p.saxpy_ab(1.0, -1.0, &yyp);
        } else {
            p.saxpy_ab(-1.0, -1.0, &yyp);
        }
    } else {
        let dot_val = gradient.get_gradient(runtime).dot(yp);
        p.repopulate(yp).scale(-dot_val);
    }
}

/// HiGHS: quass.cpp::computemaxsteplength (lines 136-151)
pub(crate) fn computemaxsteplength(
    runtime: &Runtime,
    p: &QpVector,
    gradient: &mut Gradient,
    buffer_Qp: &mut QpVector,
    zcd: &mut bool,
) -> f64 {
    runtime.instance.Q.mat_vec_into(p, buffer_Qp);
    let denominator = p.dot(buffer_Qp);
    if denominator.abs() > runtime.settings.pQp_zero_threshold {
        let numerator = -(p.dot(gradient.get_gradient(runtime)));
        if numerator < 0.0 {
            0.0
        } else {
            numerator / denominator
        }
    } else {
        *zcd = true;
        f64::INFINITY
    }
}

/// HiGHS: quass.cpp::reduce (lines 153-185)
///
/// `rt` is `&mut` because `rt.instance.A.t()` lazily caches the transpose.
pub(crate) fn reduce(
    rt: &mut Runtime,
    basis: &mut Basis,
    newactivecon: HighsInt,
    buffer_d: &mut QpVector,
    maxabsd: &mut HighsInt,
    constrainttodrop: &mut HighsInt,
    log_d: &mut f64,
) -> QpSolverStatus {
    let idx = crate::qpsolver::snippets::indexof(basis.getinactive(), &newactivecon);
    if idx != -1 {
        *maxabsd = idx;
        *constrainttodrop = newactivecon;
        QpVector::unit_into(basis.getinactive().len() as HighsInt, idx, buffer_d);
        return QpSolverStatus::OK;
    }

    let aq = rt.instance.A.t().extractcol(newactivecon);
    basis.Ztprod(&aq, buffer_d, true, newactivecon);

    *maxabsd = 0;
    for i in 0..buffer_d.num_nz {
        if buffer_d.value[buffer_d.index[i as usize] as usize].abs()
            > buffer_d.value[*maxabsd as usize].abs()
        {
            *maxabsd = buffer_d.index[i as usize];
        }
    }
    *constrainttodrop = basis.getinactive()[*maxabsd as usize];
    if buffer_d.value[*maxabsd as usize].abs() < rt.settings.d_zero_threshold {
        *log_d = buffer_d.value[*maxabsd as usize].abs().log10();
        return QpSolverStatus::DEGENERATE;
    }
    QpSolverStatus::OK
}

/// HiGHS: quass.cpp::getPricing (lines 187-201)
pub(crate) fn get_pricing(runtime: &Runtime, basis: &mut Basis) -> Box<dyn Pricing> {
    match runtime.settings.pricing {
        PricingStrategy::SteepestEdge => Box::new(SteepestEdgePricing::new(runtime, basis)),
        PricingStrategy::Devex => Box::new(DevexPricing::new(runtime)),
        PricingStrategy::DantzigWolfe => Box::new(DantzigPricing::new(runtime)),
    }
}

/// HiGHS: quass.cpp::regularize (lines 203-221)
pub(crate) fn regularize(rt: &mut Runtime) {
    debug_assert!(!rt.settings.hessian_regularization_on);
    if !rt.settings.hessian_regularization_on {
        return;
    }
    // C++ lines 208-209: comment about regularization not being called.
    for i in 0..rt.instance.num_var {
        let mut index = rt.instance.Q.mat.start[i as usize];
        while index < rt.instance.Q.mat.start[(i + 1) as usize] {
            if rt.instance.Q.mat.index[index as usize] == i {
                rt.instance.Q.mat.value[index as usize] += rt.settings.hessian_regularization_value;
            }
            index += 1;
        }
    }
}

/// HiGHS: quass.cpp::compute_actual_duals (lines 224-254)
/// C++ lines 223-287 are inside `#if 0`; translated as written in the skeleton.
pub(crate) fn compute_actual_duals(
    rt: &Runtime,
    basis: &mut Basis,
    lambda: &mut QpVector,
    dual_con: &mut QpVector,
    dual_var: &mut QpVector,
) {
    // Borrow split: `getstatus` needs `&mut basis`, so iterate over a copy of
    // the active list (C++ iterates the reference; nothing changes it here).
    let active = basis.getactive().clone();
    for &e in active.iter() {
        let indexinbasis = basis.getindexinfactor()[e as usize];
        let status = basis.getstatus(e);
        if e >= rt.instance.num_con {
            let var = e - rt.instance.num_con;

            if status == BasisStatus::kActiveAtUpper {
                dual_var.value[var as usize] = -lambda.value[indexinbasis as usize];
            } else if status == BasisStatus::kActiveAtLower {
                dual_var.value[var as usize] = lambda.value[indexinbasis as usize];
            } else {
                debug_assert!(lambda.value[indexinbasis as usize] == 0.0);
                let _ = dual_var.value[var as usize];
            }
        } else {
            if status == BasisStatus::kActiveAtUpper {
                dual_con.value[e as usize] = -lambda.value[indexinbasis as usize];
            } else if status == BasisStatus::kActiveAtLower {
                dual_con.value[e as usize] = lambda.value[indexinbasis as usize];
            } else {
                debug_assert!(lambda.value[indexinbasis as usize] == 0.0);
                dual_con.value[e as usize] = 0.0;
            }
        }
    }
    dual_con.resparsify();
    dual_var.resparsify();
}

/// HiGHS: quass.cpp::compute_primal_violation (lines 256-272)
/// C++ lines 223-287 are inside `#if 0`; translated as written in the skeleton.
pub(crate) fn compute_primal_violation(rt: &Runtime) -> f64 {
    let mut maxviolation = 0.0;
    let rowact = rt.instance.A.mat_vec(&rt.primal);
    for i in 0..rt.instance.num_con {
        let violation = rt.instance.con_lo[i as usize] - rowact.value[i as usize];
        maxviolation = cpp_max(violation, maxviolation);
        let violation = rowact.value[i as usize] - rt.instance.con_up[i as usize];
        maxviolation = cpp_max(violation, maxviolation);
    }
    for i in 0..rt.instance.num_var {
        let violation = rt.instance.var_lo[i as usize] - rt.primal.value[i as usize];
        maxviolation = cpp_max(violation, maxviolation);
        let violation = rt.primal.value[i as usize] - rt.instance.var_up[i as usize];
        maxviolation = cpp_max(violation, maxviolation);
    }
    maxviolation
}

/// HiGHS: quass.cpp::compute_dual_violation (lines 274-286)
///
/// `instance` is `&mut` because `instance.A.t()` lazily caches the transpose.
/// C++ lines 223-287 are inside `#if 0`; translated as written in the skeleton.
pub(crate) fn compute_dual_violation(
    instance: &mut Instance,
    primal: &QpVector,
    dual_con: &QpVector,
    dual_var: &QpVector,
) -> f64 {
    let mut maxviolation = 0.0;

    let q_res = instance.Q.mat_vec(primal);
    let c_res = q_res.op_add(&instance.c);
    let at_res = instance.A.t().mat_vec(dual_con);
    let ac_res = c_res.op_add(&at_res);
    let residuals = ac_res.op_add(dual_var);

    for i in 0..instance.num_var {
        let violation = residuals.value[i as usize];
        maxviolation = cpp_max(violation, maxviolation);
        let violation = -residuals.value[i as usize];
        maxviolation = cpp_max(violation, maxviolation);
    }
    maxviolation
}

/// HiGHS: quass.cpp::check_reinvert_due (lines 289-292)
pub(crate) fn check_reinvert_due(basis: &Basis) -> bool {
    basis.getreinversionhint()
}

/// HiGHS: quass.cpp::reinvert (lines 294-307)
pub(crate) fn reinvert(
    rt: &Runtime,
    basis: &mut Basis,
    factor: &mut CholeskyFactor,
    grad: &mut Gradient,
    rc: &mut ReducedCosts,
    rg: &mut ReducedGradient,
    pricing: &mut Box<dyn Pricing>,
) -> QpSolverStatus {
    basis.rebuild();
    let status = factor.recompute(rt, basis);
    if status != QpSolverStatus::OK {
        return status;
    }
    grad.recompute(rt);
    rc.recompute(rt, basis, grad);
    rg.recompute(rt, basis, grad);
    // C++ line 305: pricing->recompute() is commented out.
    QpSolverStatus::OK
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
    fn helpers_on_fresh_basis() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let basis = at_lower(&mut rt);
        assert!(!check_reinvert_due(&basis));
        let q = rt.instance.Q.mat.value.clone();
        regularize(&mut rt);
        assert_eq!(rt.instance.Q.mat.value, q);
        assert_eq!(compute_primal_violation(&rt), 0.0);
        let _quass = Quass::new(&rt);
    }
}
