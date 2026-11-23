//! HiGHS source: ipm/hipo/ipm/PreProcess.cpp (lines 289-684)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Scaling (Curtis-Reid style power-of-two scaling), formulation (slack
//! columns for inequality rows), free variables, and the `Preprocessor`
//! stack. `APPLY_ACTION(T)` pushes `Box::new(T::default())` and applies it.

use super::{
    PreprocessAction, PreprocessActionBase, PreprocessEmptyRows, PreprocessFixedVars,
    PreprocessorPoint,
};
use crate::ipm::hipo::Int;
use crate::ipm::hipo::auxiliary::vector_operations::{isInfVector, isNanVector};
use crate::ipm::hipo::ipm::iterate::Iterate;
use crate::ipm::hipo::ipm::model::Model;
use crate::ipm::hipo::ipm::parameters::{
    kFreeVarsInitialBound, kLargeScalingCoeff, kSmallBoundDiff, kSmallScalingCoeff,
};
use crate::lp_data::hconst::kHighsInf;
use crate::util::highs_cdouble::ldexp_f64;
use crate::util::highs_hash::frexp;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// HiGHS: PreProcess.cpp::roundToPowerOf2 (lines 289-300)
/// `std::frexp` / `std::ldexp`; the exponent is a C `int`.
pub fn roundToPowerOf2(d: f64) -> f64 {
    // HiGHS: PreProcess.cpp::roundToPowerOf2 (lines 289-300)
    // int exp;  // int, not Int

    // d = x * 2^exp, d\in[0.5,1)
    let (x, exp) = frexp(d);

    // d is already exact power of 2
    if x == 0.5 {
        return d;
    }

    // return 1*2^exp
    ldexp_f64(1.0, exp)
}

/// HiGHS: PreProcess.h `struct PreprocessScaling`
#[derive(Clone, Debug, Default)]
pub struct PreprocessScaling {
    pub base: PreprocessActionBase,
    pub scaled: bool,
}

impl PreprocessScaling {
    /// HiGHS: PreProcess.cpp::PreprocessScaling::apply (lines 336-376) colScaling lambda
    fn col_scaling(model: &mut Model) {
        let n = model.n_;
        let m = model.m_;
        let qp = model.qp();

        // infinity norm of columns
        let mut norm_cols = vec![0.0f64; n as usize];
        for col in 0..n {
            for el in model.A_.start_[col as usize]..model.A_.start_[(col + 1) as usize] {
                let row = model.A_.index_[el as usize];
                let mut value = model.A_.value_[el as usize];
                value *= model.colscale_[col as usize];
                value *= model.rowscale_[row as usize];
                norm_cols[col as usize] = cpp_max(norm_cols[col as usize], value.abs());
            }
        }

        // apply col scaling
        for i in 0..n {
            let mut coeff = 1.0 / norm_cols[i as usize].sqrt();

            // use geomean of norm scaling and scaling that makes diagonal of Q equal
            // to 1
            if qp {
                let qii =
                    model.Q_.diag(i) * model.colscale_[i as usize] * model.colscale_[i as usize];
                coeff = (coeff * 1.0 / qii.abs().sqrt()).sqrt();
            }

            // if bounds interval is too small, use geomean of norm scaling
            // and maximum scaling allowed by bounds
            if model.lower_[i as usize].is_finite() && model.upper_[i as usize].is_finite() {
                let l = model.lower_[i as usize] / model.colscale_[i as usize];
                let u = model.upper_[i as usize] / model.colscale_[i as usize];
                let diff = (u - l).abs();

                if diff / coeff < kSmallBoundDiff {
                    coeff = (coeff * diff / kSmallBoundDiff).sqrt();
                }
            }

            if !coeff.is_infinite() && !coeff.is_nan() {
                model.colscale_[i as usize] *= coeff;
            }

            model.colscale_[i as usize] = cpp_max(model.colscale_[i as usize], kSmallScalingCoeff);
            model.colscale_[i as usize] = cpp_min(model.colscale_[i as usize], kLargeScalingCoeff);
        }
    }

    /// HiGHS: PreProcess.cpp::PreprocessScaling::apply (lines 377-396) rowScaling lambda
    fn row_scaling(model: &mut Model) {
        let n = model.n_;
        let m = model.m_;

        // infinity norm of rows
        let mut norm_rows = vec![0.0f64; m as usize];
        for col in 0..n {
            for el in model.A_.start_[col as usize]..model.A_.start_[(col + 1) as usize] {
                let row = model.A_.index_[el as usize];
                let mut value = model.A_.value_[el as usize];
                value *= model.colscale_[col as usize];
                value *= model.rowscale_[row as usize];
                norm_rows[row as usize] = cpp_max(norm_rows[row as usize], value.abs());
            }
        }

        // apply row scaling
        for i in 0..m {
            if norm_rows[i as usize] > 0.0 {
                model.rowscale_[i as usize] *= 1.0 / norm_rows[i as usize].sqrt();
            }
            model.rowscale_[i as usize] = cpp_max(model.rowscale_[i as usize], kSmallScalingCoeff);
            model.rowscale_[i as usize] = cpp_min(model.rowscale_[i as usize], kLargeScalingCoeff);
        }
    }
}

impl PreprocessAction for PreprocessScaling {
    fn base(&self) -> &PreprocessActionBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut PreprocessActionBase {
        &mut self.base
    }

    /// HiGHS: PreProcess.cpp::PreprocessScaling::apply (lines 302-463)
    fn apply(&mut self, model: &mut Model) {
        let n = model.n_;
        let m = model.m_;

        self.base.n_pre = n;
        self.base.m_pre = m;
        self.base.n_post = n;
        self.base.m_post = m;

        // *********************************************************************
        // Compute scaling
        // *********************************************************************
        // Transformation:
        // A -> R * A * C
        // b -> R * b
        // c -> C * c
        // x -> C^-1 * x
        // y -> R^-1 * y
        // z -> C * z
        // Q -> C * Q * C
        // where R is row scaling, C is col scaling.

        model.colscale_.clear();
        model.colscale_.resize(n as usize, 1.0);
        model.rowscale_.clear();
        model.rowscale_.resize(m as usize, 1.0);

        let num_passes: Int = 10;
        for _pass in 0..num_passes {
            Self::row_scaling(model);
            Self::col_scaling(model);
        }

        let scaling_failed = isInfVector(&model.colscale_)
            || isNanVector(&model.colscale_)
            || isInfVector(&model.rowscale_)
            || isNanVector(&model.rowscale_);
        if scaling_failed {
            model.colscale_.clear();
            model.rowscale_.clear();
            return;
        }

        self.scaled = true;

        // *********************************************************************
        // Adjust scaling
        // *********************************************************************

        for i in 0..n {
            model.colscale_[i as usize] = cpp_max(model.colscale_[i as usize], kSmallScalingCoeff);
            model.colscale_[i as usize] = cpp_min(model.colscale_[i as usize], kLargeScalingCoeff);
            model.colscale_[i as usize] = roundToPowerOf2(model.colscale_[i as usize]);
        }
        for i in 0..m {
            model.rowscale_[i as usize] = cpp_max(model.rowscale_[i as usize], kSmallScalingCoeff);
            model.rowscale_[i as usize] = cpp_min(model.rowscale_[i as usize], kLargeScalingCoeff);
            model.rowscale_[i as usize] = roundToPowerOf2(model.rowscale_[i as usize]);
        }

        // *********************************************************************
        // Apply scaling
        // *********************************************************************

        // Column has been scaled up by colscale_[col], so cost is scaled up and
        // bounds are scaled down
        for col in 0..n {
            model.c_[col as usize] *= model.colscale_[col as usize];
            model.lower_[col as usize] /= model.colscale_[col as usize];
            model.upper_[col as usize] /= model.colscale_[col as usize];
        }

        // Row has been scaled up by rowscale_[row], so b is scaled up
        for row in 0..m {
            model.b_[row as usize] *= model.rowscale_[row as usize];
        }

        // Each entry of the matrix is scaled by the corresponding row and col
        // factor
        for col in 0..n {
            for el in model.A_.start_[col as usize]..model.A_.start_[(col + 1) as usize] {
                let row = model.A_.index_[el as usize];
                model.A_.value_[el as usize] *= model.rowscale_[row as usize];
                model.A_.value_[el as usize] *= model.colscale_[col as usize];
            }
        }

        if model.qp() {
            for col in 0..model.Q_.dim_ {
                for el in model.Q_.start_[col as usize]..model.Q_.start_[(col + 1) as usize] {
                    let row = model.Q_.index_[el as usize];
                    model.Q_.value_[el as usize] *= model.colscale_[row as usize];
                    model.Q_.value_[el as usize] *= model.colscale_[col as usize];
                }
            }
        }
    }

    /// HiGHS: PreProcess.cpp::PreprocessScaling::undo (lines 465-503)
    fn undo(&self, point: &mut PreprocessorPoint<'_>, model: &Model, it: &Iterate) {
        point.assertConsistency(self.base.n_post, self.base.m_post);
        if model.scaled() {
            let colscale = &model.colscale_;
            let rowscale = &model.rowscale_;

            for i in 0..self.base.n_post {
                point.x[i as usize] *= colscale[i as usize];
                point.xl[i as usize] *= colscale[i as usize];
                point.xu[i as usize] *= colscale[i as usize];
                point.zl[i as usize] /= colscale[i as usize];
                point.zu[i as usize] /= colscale[i as usize];
            }
            for i in 0..self.base.m_post {
                point.y[i as usize] *= rowscale[i as usize];
                point.slack[i as usize] /= rowscale[i as usize];
            }

            // set variables that were ignored
            for i in 0..self.base.n_post {
                if !model.hasLb(i) {
                    point.xl[i as usize] = kHighsInf;
                    point.zl[i as usize] = 0.0;
                }
                if !model.hasUb(i) {
                    point.xu[i as usize] = kHighsInf;
                    point.zu[i as usize] = 0.0;
                }
            }
            for &i in model.freeVars() {
                point.xl[i as usize] = kHighsInf;
                point.zl[i as usize] = 0.0;
                point.xu[i as usize] = kHighsInf;
                point.zu[i as usize] = 0.0;
            }
        }
        point.assertConsistency(self.base.n_pre, self.base.m_pre);
    }

    /// HiGHS: PreProcess.cpp::PreprocessScaling::print (lines 505-505)
    fn print(&self, stream: &mut String) {
        // empty body in C++
    }
}

/// HiGHS: PreProcess.h `struct PreprocessFormulation`
#[derive(Clone, Debug, Default)]
pub struct PreprocessFormulation {
    pub base: PreprocessActionBase,
}

impl PreprocessAction for PreprocessFormulation {
    fn base(&self) -> &PreprocessActionBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut PreprocessActionBase {
        &mut self.base
    }

    /// HiGHS: PreProcess.cpp::PreprocessFormulation::apply (lines 507-555)
    fn apply(&mut self, model: &mut Model) {
        let m = model.m_;

        self.base.n_pre = model.n_;
        self.base.m_pre = m;

        let _annz = model.A_.num_nz();

        for i in 0..m {
            if model.constraints_[i as usize] != b'=' {
                // inequality constraint, add slack variable

                model.n_ += 1;

                // lower/upper bound for new slack
                if model.constraints_[i as usize] == b'>' {
                    model.lower_.push(-kHighsInf);
                    model.upper_.push(0.0);
                } else {
                    model.lower_.push(0.0);
                    model.upper_.push(kHighsInf);
                }

                // cost for new slack
                model.c_.push(0.0);

                // add column of identity to A_
                let temp_ind: Vec<Int> = vec![i];
                let temp_val: Vec<f64> = vec![1.0];
                model.A_.add_vec(1, &temp_ind, &temp_val, 1.0);

                // set scaling to 1
                if model.scaled() {
                    model.colscale_.push(1.0);
                }
            }
        }

        if model.qp() {
            crate::model::highs_hessian_utils::complete_hessian(model.n_, &mut model.Q_);
        }

        self.base.n_post = model.n_;
        self.base.m_post = m;
    }

    /// HiGHS: PreProcess.cpp::PreprocessFormulation::undo (lines 557-626)
    fn undo(&self, point: &mut PreprocessorPoint<'_>, model: &Model, it: &Iterate) {
        it.assertConsistency(self.base.n_post, self.base.m_post);

        // Copy x, xl, xu, zl, zu without slacks
        *point.x = it.x[..self.base.n_pre as usize].to_vec();
        *point.xl = it.xl[..self.base.n_pre as usize].to_vec();
        *point.xu = it.xu[..self.base.n_pre as usize].to_vec();
        *point.zl = it.zl[..self.base.n_pre as usize].to_vec();
        *point.zu = it.zu[..self.base.n_pre as usize].to_vec();

        // force unused entries to have correct value
        for i in 0..self.base.n_pre {
            if !model.hasLb(i) {
                point.xl[i as usize] = kHighsInf;
                point.zl[i as usize] = 0.0;
            }
            if !model.hasUb(i) {
                point.xu[i as usize] = kHighsInf;
                point.zu[i as usize] = 0.0;
            }
        }
        for &i in model.freeVars() {
            point.xl[i as usize] = kHighsInf;
            point.zl[i as usize] = 0.0;
            point.xu[i as usize] = kHighsInf;
            point.zu[i as usize] = 0.0;
        }

        // For the Lagrange multipliers, use slacks from zl and zu, to get correct
        // sign. NB: there is no explicit slack stored for equality constraints.
        point.y.resize(self.base.m_pre as usize, 0.0);
        let mut slack_pos: Int = 0;
        for i in 0..self.base.m_pre {
            match model.constraint(i) {
                b'=' => {
                    point.y[i as usize] = it.y[i as usize];
                }
                b'>' => {
                    point.y[i as usize] = it.zu[(self.base.n_pre + slack_pos) as usize];
                    slack_pos += 1;
                }
                b'<' => {
                    point.y[i as usize] = -it.zl[(self.base.n_pre + slack_pos) as usize];
                    slack_pos += 1;
                }
                _ => {}
            }
        }

        // For x-slacks, use slacks from xl and xu, to get correct sign.
        // NB: there is no explicit slack stored for equality constraints.
        point.slack.resize(self.base.m_pre as usize, 0.0);
        slack_pos = 0;
        for i in 0..self.base.m_pre {
            match model.constraint(i) {
                b'=' => {
                    point.slack[i as usize] = 0.0;
                }
                b'>' => {
                    point.slack[i as usize] = -it.xu[(self.base.n_pre + slack_pos) as usize];
                    slack_pos += 1;
                }
                b'<' => {
                    point.slack[i as usize] = it.xl[(self.base.n_pre + slack_pos) as usize];
                    slack_pos += 1;
                }
                _ => {}
            }
        }
        point.assertConsistency(self.base.n_pre, self.base.m_pre);
    }

    /// HiGHS: PreProcess.cpp::PreprocessFormulation::print (lines 628-630)
    fn print(&self, stream: &mut String) {
        stream.push_str(&format!(
            "Added {} slacks\n",
            self.base.n_post - self.base.n_pre
        ));
    }
}

/// HiGHS: PreProcess.h `struct PreprocessFreeVars`
#[derive(Clone, Debug, Default)]
pub struct PreprocessFreeVars {
    pub base: PreprocessActionBase,
    pub free_vars_count: Int,
}

impl PreprocessAction for PreprocessFreeVars {
    fn base(&self) -> &PreprocessActionBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut PreprocessActionBase {
        &mut self.base
    }

    /// HiGHS: PreProcess.cpp::PreprocessFreeVars::apply (lines 632-649)
    fn apply(&mut self, model: &mut Model) {
        let n = model.n_;

        model.free_variables_.clear();
        for i in 0..n {
            if !model.lower_[i as usize].is_finite()
                && !model.upper_[i as usize].is_finite()
                && model.lower_[i as usize] != model.upper_[i as usize]
            {
                // free variable
                model.free_variables_.push(i);
                self.free_vars_count += 1;
                model.lower_[i as usize] = -kFreeVarsInitialBound;
                model.upper_[i as usize] = kFreeVarsInitialBound;
            }
        }
    }

    /// HiGHS: PreProcess.cpp::PreprocessFreeVars::undo (lines 651-652)
    fn undo(&self, point: &mut PreprocessorPoint<'_>, model: &Model, it: &Iterate) {
        // empty body in C++
    }

    /// HiGHS: PreProcess.cpp::PreprocessFreeVars::print (lines 654-657)
    fn print(&self, stream: &mut String) {
        if self.free_vars_count > 0 {
            stream.push_str(&format!("Found {} free variables\n", self.free_vars_count));
        }
    }
}

/// HiGHS: PreProcess.h `struct Preprocessor`
#[derive(Default)]
pub struct Preprocessor {
    pub stack: Vec<Box<dyn PreprocessAction>>,
}

impl Preprocessor {
    /// HiGHS: PreProcess.cpp::Preprocessor::apply (lines 663-672)
    /// Order: FixedVars, EmptyRows, Scaling, Formulation, FreeVars.
    pub fn apply(&mut self, model: &mut Model) {
        // Remove fixed variables before removing empty rows, because removing columns
        // may create empty rows

        // APPLY_ACTION(PreprocessFixedVars);
        {
            let mut action = Box::new(PreprocessFixedVars::default());
            action.apply(model);
            self.stack.push(action);
        }
        // APPLY_ACTION(PreprocessEmptyRows);
        {
            let mut action = Box::new(PreprocessEmptyRows::default());
            action.apply(model);
            self.stack.push(action);
        }
        // APPLY_ACTION(PreprocessScaling);
        {
            let mut action = Box::new(PreprocessScaling::default());
            action.apply(model);
            self.stack.push(action);
        }
        // APPLY_ACTION(PreprocessFormulation);
        {
            let mut action = Box::new(PreprocessFormulation::default());
            action.apply(model);
            self.stack.push(action);
        }
        // APPLY_ACTION(PreprocessFreeVars);
        {
            let mut action = Box::new(PreprocessFreeVars::default());
            action.apply(model);
            self.stack.push(action);
        }
    }

    /// HiGHS: PreProcess.cpp::Preprocessor::undo (lines 673-678)
    pub fn undo(&self, point: &mut PreprocessorPoint<'_>, model: &Model, it: &Iterate) {
        for iterator in self.stack.iter().rev() {
            iterator.undo(point, model, it);
        }
    }

    /// HiGHS: PreProcess.cpp::Preprocessor::print (lines 680-684)
    pub fn print(&self, log_stream: &mut String) {
        for iterator in self.stack.iter() {
            iterator.print(log_stream);
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    #[ignore = "not yet filled"]
    fn round_to_power_of_2() {
        assert_eq!(roundToPowerOf2(0.5), 0.5);
        assert_eq!(roundToPowerOf2(3.0), 4.0);
        assert_eq!(roundToPowerOf2(4.0), 4.0);
        assert_eq!(roundToPowerOf2(0.3), 0.5);
    }
}
