//! HiGHS source: model/HighsModel.h, model/HighsModel.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `operator==` / `equalButForNames` (HighsModel.cpp lines 15-25)
//! are deferred: they need `HighsLp::equalButForNames` / `equalNames`, which
//! the crate does not port (no solver path compares models). `needsMods` /
//! `hasMods` are dropped with `HighsLp::mods_` (see `highs_lp.rs`).

use crate::model::highs_hessian::HighsHessian;
use crate::prelude::*;

/// HiGHS: HighsModel.h::HighsModel (lines 21-40)
#[derive(Clone, Debug, Default)]
pub struct HighsModel {
    pub lp_: HighsLp,
    pub hessian_: HighsHessian,
}

impl HighsModel {
    /// HiGHS: HighsModel.h::HighsModel::isQp (line 27)
    pub fn is_qp(&self) -> bool {
        self.hessian_.dim_ != 0
    }

    /// HiGHS: HighsModel.h::HighsModel::isMip (line 28)
    pub fn is_mip(&self) -> bool {
        self.lp_.is_mip()
    }

    /// HiGHS: HighsModel.h::HighsModel::isEmpty (lines 29-31)
    pub fn is_empty(&self) -> bool {
        self.lp_.num_col_ == 0 && self.lp_.num_row_ == 0
    }

    /// HiGHS: HighsModel.cpp::HighsModel::clear (lines 27-30)
    ///
    /// Architect-filled (reached by `HighsIis::clearData` on the LP path).
    pub fn clear(&mut self) {
        self.lp_.clear();
        self.hessian_.clear();
    }

    /// HiGHS: HighsModel.cpp::HighsModel::objectiveValue (lines 32-35)
    ///
    /// Delegates to [`model_objective_value`], the form `Highs` (which keeps
    /// `model_lp_` and `model_hessian_` apart) calls.
    ///
    /// HiGHS 1.16 #3513: `h16_obj_cdouble` is passed through to
    /// [`HighsLp::objective_value`] (no `HighsOptions` in reach).
    pub fn objective_value(&self, solution: &[f64], h16_obj_cdouble: bool) -> f64 {
        model_objective_value(&self.lp_, &self.hessian_, solution, h16_obj_cdouble)
    }

    /// HiGHS: HighsModel.cpp::HighsModel::objectiveGradient (lines 37-46)
    ///
    /// Delegates to [`model_objective_gradient`].
    pub fn objective_gradient(&self, solution: &[f64], gradient: &mut Vec<f64>) {
        model_objective_gradient(&self.lp_, &self.hessian_, solution, gradient)
    }
}

/// HiGHS: HighsModel.cpp::HighsModel::objectiveValue (lines 32-35)
///
/// PORT NOTE: `HighsModel` split into its `lp_` and `hessian_` so that
/// `Highs` can call it on `model_lp_` / `model_hessian_` without building a
/// `HighsModel`.
///
/// HiGHS 1.16 #3513 (h16_obj_cdouble): HighsModel.cpp lines 58-61 (1.16.0) reach
/// HighsLp::objectiveValue, whose CDouble sum the flag selects; off, the
/// 1.15.1 sum is used.
pub fn model_objective_value(
    lp: &HighsLp,
    hessian: &HighsHessian,
    solution: &[f64],
    h16_obj_cdouble: bool,
) -> f64 {
    hessian.objective_value(solution) + lp.objective_value(solution, h16_obj_cdouble)
}

/// HiGHS: HighsModel.cpp::HighsModel::objectiveGradient (lines 37-46)
///
/// PORT NOTE: split form, as [`model_objective_value`].
pub fn model_objective_gradient(
    lp: &HighsLp,
    hessian: &HighsHessian,
    solution: &[f64],
    gradient: &mut Vec<f64>,
) {
    if hessian.dim_ > 0 {
        hessian.product(solution, gradient);
    } else {
        gradient.clear();
        gradient.resize(lp.num_col_ as usize, 0.0);
    }
    let mut i_col: HighsInt = 0;
    while i_col < lp.num_col_ {
        gradient[i_col as usize] += lp.col_cost_[i_col as usize];
        i_col += 1;
    }
}

/// HiGHS: HighsSolution.cpp::getKktFailures(options, model, solution, basis,
/// highs_info) (lines 30-36) -> getKktFailures(.., primal_dual_errors,
/// get_residuals = false) (lines 38-51)
///
/// PORT NOTE: split-model form (see [`model_objective_value`]). The two
/// overloads are folded into one: `getPrimalDualBasisErrors` /
/// `getPrimalDualGlpsolErrors` only write the discarded
/// `primal_dual_errors`, so they are not ported (as for the LP form in
/// `presolve_support.rs`). Body: `model_objective_gradient`, then
/// `get_kkt_failures(options, hessian.dim_ != 0, lp, &gradient, solution,
/// highs_info, false)`.
pub fn get_kkt_failures_model(
    options: &HighsOptions,
    lp: &HighsLp,
    hessian: &HighsHessian,
    solution: &HighsSolution,
    basis: &HighsBasis,
    highs_info: &mut HighsInfo,
) {
    let mut gradient: Vec<f64> = Vec::new();
    model_objective_gradient(lp, hessian, &solution.col_value, &mut gradient);
    let is_qp: bool = hessian.dim_ != 0;
    crate::lp_data::presolve_support::get_kkt_failures(
        options, is_qp, lp, &gradient, solution, highs_info, false,
    );
    // C++ lines 48-50: getPrimalDualBasisErrors / getPrimalDualGlpsolErrors
    // write only into the discarded primal_dual_errors; not ported.
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

    /// [[2, 1], [1, 4]] in triangular (lower, column-wise) format.
    fn h2() -> HighsHessian {
        HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 2, 3],
            index_: vec![0, 1, 1],
            value_: vec![2.0, 1.0, 4.0],
        }
    }

    #[test]
    fn objective_and_gradient() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.col_cost_ = vec![1.0, -1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf];
        lp.offset_ = 0.5;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.start_ = vec![0, 0, 0];
        let x = [1.0, 1.0];
        // 0.5 + (1 - 1) + 4
        assert_eq!(model_objective_value(&lp, &h2(), &x, false), 4.5);
        let mut g = Vec::new();
        model_objective_gradient(&lp, &h2(), &x, &mut g);
        assert_eq!(g, vec![4.0, 4.0]);
        let model = HighsModel {
            lp_: lp,
            hessian_: h2(),
        };
        assert!(model.is_qp());
        assert!(!model.is_mip());
        assert!(!model.is_empty());
        assert_eq!(model.objective_value(&x, false), 4.5);
    }

    #[test]
    fn gradient_no_hessian() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.col_cost_ = vec![1.0, 2.0, 3.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf, kHighsInf];
        lp.offset_ = 0.0;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.start_ = vec![0, 0, 0, 0];
        let h = HighsHessian::default();
        let x = [1.0, 2.0, 3.0];
        let mut g = Vec::new();
        model_objective_gradient(&lp, &h, &x, &mut g);
        assert_eq!(g, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn model_clear() {
        let mut model = HighsModel {
            lp_: HighsLp::default(),
            hessian_: h2(),
        };
        model.clear();
        assert_eq!(model.hessian_.dim_, 0);
        assert_eq!(model.lp_.num_col_, 0);
    }
}
