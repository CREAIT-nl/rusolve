//! HiGHS source: ipm/ipx/model.cc (lines 10-497)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::super::control::{Control, textline};
use super::super::sparse_matrix::{dot_column, scatter_column, transpose};
use super::super::utils::{dot, infnorm};
use super::*;
use crate::util::highs_utils::cpp_max;

impl Model {
    /// HiGHS: ipx/model.cc::Model::Load (lines 10-73)
    ///
    /// Nullable C pointers are `Option` (see the module doc): any `None`
    /// returns `IPX_ERROR_argument_null` from `CopyInput`.
    pub fn load(
        &mut self,
        control: &Control,
        num_constr: Int,
        num_var: Int,
        Ap: Option<&[Int]>,
        Ai: Option<&[Int]>,
        Ax: Option<&[f64]>,
        rhs: Option<&[f64]>,
        constr_type: Option<&[u8]>,
        offset: f64,
        obj: Option<&[f64]>,
        lbuser: Option<&[f64]>,
        ubuser: Option<&[f64]>,
    ) -> Int {
        self.clear();
        let errflag = self.copy_input(
            num_constr,
            num_var,
            Ap,
            Ai,
            Ax,
            rhs,
            constr_type,
            offset,
            obj,
            lbuser,
            ubuser,
        );
        if errflag != 0 {
            return errflag;
        }
        let mut logging = String::new();
        logging.push_str(&format!(
            "Input\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n",
            textline("Number of variables:"),
            self.num_var_,
            textline("Number of free variables:"),
            self.num_free_var_,
            textline("Number of constraints:"),
            self.num_constr_,
            textline("Number of equality constraints:"),
            self.num_eqconstr_,
            textline("Number of matrix entries:"),
            self.num_entries_,
        ));
        control.h_log(&mut logging);
        self.print_coefficient_range(control);

        self.scale_model(control);

        // Make an automatic decision for dualization if not specified by user.
        let mut dualize = control.dualize();
        // dualize = -2 => Possibly dualize - Filippo style
        // dualize = -1 => Possibly dualize - Lukas style
        // dualize = 0 => No dualization
        // dualize = 1 => Perform dualization
        debug_assert!(dualize == 0 || dualize == -1);
        let dualize_lukas = num_constr > 2 * num_var;
        let dualize_filippo = self.filippo_dualization_test();
        if dualize == -1 {
            dualize = Int::from(dualize_lukas);
        } else if dualize == -2 {
            dualize = Int::from(dualize_filippo);
        }
        if dualize != 0 {
            self.load_dual();
        } else {
            self.load_primal();
        }

        self.A_.clear();
        self.AIt_ = transpose(&self.AI_);
        debug_assert_eq!(
            self.AI_.begin(self.num_cols_ + self.num_rows_),
            self.AIt_.begin(self.num_rows_)
        );
        self.find_dense_columns();
        self.norm_c_ = infnorm(&self.c_);
        self.norm_bounds_ = infnorm(&self.b_);
        for &x in self.lb_.iter() {
            if x.is_finite() {
                self.norm_bounds_ = cpp_max(self.norm_bounds_, x.abs());
            }
        }
        for &x in self.ub_.iter() {
            if x.is_finite() {
                self.norm_bounds_ = cpp_max(self.norm_bounds_, x.abs());
            }
        }

        self.print_preprocessing_log(control);
        0
    }

    /// HiGHS: ipx/model.cc::Model::filippoDualizationTest (lines 75-77)
    pub fn filippo_dualization_test(&self) -> bool {
        false
    }

    /// HiGHS: ipx/model.cc::Model::GetInfo (lines 79-88)
    pub fn get_info(&self, info: &mut Info) {
        info.num_var = self.num_var_;
        info.num_constr = self.num_constr_;
        info.num_entries = self.num_entries_;
        info.num_rows_solver = self.num_rows_;
        info.num_cols_solver = self.num_cols_ + self.num_rows_; // including slack columns
        info.num_entries_solver = self.AI_.entries();
        info.dualized = Int::from(self.dualized_);
        info.dense_cols = self.num_dense_cols();
    }

    /// HiGHS: ipx/model.cc::Model::clear (lines 90-125)
    ///
    /// IPX quirk: `offset_` is not reset (C++ `clear()` leaves it alone).
    pub fn clear(&mut self) {
        // clear computational form model
        self.dualized_ = false;
        self.num_rows_ = 0;
        self.num_cols_ = 0;
        self.num_dense_cols_ = 0;
        self.nz_dense_ = 0;
        self.AI_.clear();
        self.AIt_.clear();
        self.b_.resize(0);
        self.c_.resize(0);
        self.lb_.resize(0);
        self.ub_.resize(0);
        self.norm_bounds_ = 0.0;
        self.norm_c_ = 0.0;

        // clear user model
        self.num_constr_ = 0;
        self.num_eqconstr_ = 0;
        self.num_var_ = 0;
        self.num_free_var_ = 0;
        self.num_entries_ = 0;
        self.boxed_vars_.clear();
        self.constr_type_.clear();
        self.norm_obj_ = 0.0;
        self.norm_rhs_ = 0.0;
        self.scaled_obj_.resize(0);
        self.scaled_rhs_.resize(0);
        self.scaled_lbuser_.resize(0);
        self.scaled_ubuser_.resize(0);
        self.A_.clear();

        self.flipped_vars_.clear();
        self.colscale_.resize(0);
        self.rowscale_.resize(0);
    }

    /// HiGHS: ipx/model.cc::Model::PresolveStartingPoint (lines 127-155)
    ///
    /// Each `None` user vector is taken as all zeros.
    pub fn presolve_starting_point(
        &self,
        x_user: Option<&[f64]>,
        slack_user: Option<&[f64]>,
        y_user: Option<&[f64]>,
        z_user: Option<&[f64]>,
        x_solver: &mut Vector,
        y_solver: &mut Vector,
        z_solver: &mut Vector,
    ) {
        let m = self.rows();
        let n = self.cols();
        debug_assert_eq!(x_solver.size(), n + m);
        debug_assert_eq!(y_solver.size(), m);
        debug_assert_eq!(z_solver.size(), n + m);

        let mut x_temp = Vector::new(self.num_var_);
        let mut slack_temp = Vector::new(self.num_constr_);
        let mut y_temp = Vector::new(self.num_constr_);
        let mut z_temp = Vector::new(self.num_var_);
        if let Some(x) = x_user {
            copy_n(x, self.num_var_, &mut x_temp[..]);
        }
        if let Some(slack) = slack_user {
            copy_n(slack, self.num_constr_, &mut slack_temp[..]);
        }
        if let Some(y) = y_user {
            copy_n(y, self.num_constr_, &mut y_temp[..]);
        }
        if let Some(z) = z_user {
            copy_n(z, self.num_var_, &mut z_temp[..]);
        }
        self.scale_point(&mut x_temp, &mut slack_temp, &mut y_temp, &mut z_temp);
        self.dualize_basic_solution(
            &x_temp,
            &slack_temp,
            &y_temp,
            &z_temp,
            x_solver,
            y_solver,
            z_solver,
        );
    }

    /// HiGHS: ipx/model.cc::Model::PresolveIPMStartingPoint (lines 157-232)
    ///
    /// Any `None` user vector returns `IPX_ERROR_argument_null` (line 172).
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    pub fn presolve_ipm_starting_point(
        &self,
        x_user: Option<&[f64]>,
        xl_user: Option<&[f64]>,
        xu_user: Option<&[f64]>,
        slack_user: Option<&[f64]>,
        y_user: Option<&[f64]>,
        zl_user: Option<&[f64]>,
        zu_user: Option<&[f64]>,
        x_solver: &mut Vector,
        xl_solver: &mut Vector,
        xu_solver: &mut Vector,
        y_solver: &mut Vector,
        zl_solver: &mut Vector,
        zu_solver: &mut Vector,
    ) -> Int {
        let (
            Some(x_user),
            Some(xl_user),
            Some(xu_user),
            Some(slack_user),
            Some(y_user),
            Some(zl_user),
            Some(zu_user),
        ) = (
            x_user, xl_user, xu_user, slack_user, y_user, zl_user, zu_user,
        )
        else {
            return IPX_ERROR_argument_null;
        };
        if self.dualized_ {
            return IPX_ERROR_not_implemented;
        }
        let nv = self.num_var_ as usize;
        let nc = self.num_constr_ as usize;

        // Copy user point into workspace and apply model scaling.
        let mut x_temp = Vector::from_slice(&x_user[..nv]);
        let mut xl_temp = Vector::from_slice(&xl_user[..nv]);
        let mut xu_temp = Vector::from_slice(&xu_user[..nv]);
        let mut slack_temp = Vector::from_slice(&slack_user[..nc]);
        let mut y_temp = Vector::from_slice(&y_user[..nc]);
        let mut zl_temp = Vector::from_slice(&zl_user[..nv]);
        let mut zu_temp = Vector::from_slice(&zu_user[..nv]);
        self.scale_point_bounds(
            &mut x_temp,
            &mut xl_temp,
            &mut xu_temp,
            &mut slack_temp,
            &mut y_temp,
            &mut zl_temp,
            &mut zu_temp,
        );

        // Check that point is compatible with bounds.
        for j in 0..nv {
            if !x_temp[j].is_finite() {
                return IPX_ERROR_invalid_vector;
            }
        }
        for j in 0..nv {
            if !(xl_temp[j] >= 0.0)
                || (self.scaled_lbuser_[j] == f64::NEG_INFINITY && xl_temp[j] != f64::INFINITY)
                || (self.scaled_lbuser_[j] != f64::NEG_INFINITY && xl_temp[j] == f64::INFINITY)
            {
                return IPX_ERROR_invalid_vector;
            }
        }
        for j in 0..nv {
            if !(xu_temp[j] >= 0.0)
                || (self.scaled_ubuser_[j] == f64::INFINITY && xu_temp[j] != f64::INFINITY)
                || (self.scaled_ubuser_[j] != f64::INFINITY && xu_temp[j] == f64::INFINITY)
            {
                return IPX_ERROR_invalid_vector;
            }
        }
        for i in 0..nc {
            if !slack_temp[i].is_finite()
                || (self.constr_type_[i] == b'=' && !(slack_temp[i] == 0.0))
                || (self.constr_type_[i] == b'<' && !(slack_temp[i] >= 0.0))
                || (self.constr_type_[i] == b'>' && !(slack_temp[i] <= 0.0))
            {
                return IPX_ERROR_invalid_vector;
            }
        }
        for i in 0..nc {
            if !y_temp[i].is_finite()
                || (self.constr_type_[i] == b'<' && !(y_temp[i] <= 0.0))
                || (self.constr_type_[i] == b'>' && !(y_temp[i] >= 0.0))
            {
                return IPX_ERROR_invalid_vector;
            }
        }
        for j in 0..nv {
            if !(zl_temp[j] >= 0.0 && zl_temp[j] < f64::INFINITY)
                || (self.scaled_lbuser_[j] == f64::NEG_INFINITY && zl_temp[j] != 0.0)
            {
                return IPX_ERROR_invalid_vector;
            }
        }
        for j in 0..nv {
            if !(zu_temp[j] >= 0.0 && zu_temp[j] < f64::INFINITY)
                || (self.scaled_ubuser_[j] == f64::INFINITY && zu_temp[j] != 0.0)
            {
                return IPX_ERROR_invalid_vector;
            }
        }

        self.dualize_ipm_starting_point(
            &x_temp,
            &xl_temp,
            &xu_temp,
            &slack_temp,
            &y_temp,
            &zl_temp,
            &zu_temp,
            x_solver,
            xl_solver,
            xu_solver,
            y_solver,
            zl_solver,
            zu_solver,
        );
        0
    }

    /// HiGHS: ipx/model.cc::Model::PostsolveInteriorSolution (lines 235-281)
    ///
    /// Each `None` output is not returned.
    pub fn postsolve_interior_solution(
        &self,
        x_solver: &Vector,
        xl_solver: &Vector,
        xu_solver: &Vector,
        y_solver: &Vector,
        zl_solver: &Vector,
        zu_solver: &Vector,
        x_user: Option<&mut [f64]>,
        xl_user: Option<&mut [f64]>,
        xu_user: Option<&mut [f64]>,
        slack_user: Option<&mut [f64]>,
        y_user: Option<&mut [f64]>,
        zl_user: Option<&mut [f64]>,
        zu_user: Option<&mut [f64]>,
    ) {
        let m = self.rows();
        let n = self.cols();
        debug_assert_eq!(x_solver.size(), n + m);
        debug_assert_eq!(xl_solver.size(), n + m);
        debug_assert_eq!(xu_solver.size(), n + m);
        debug_assert_eq!(y_solver.size(), m);
        debug_assert_eq!(zl_solver.size(), n + m);
        debug_assert_eq!(zu_solver.size(), n + m);

        let mut x_temp = Vector::new(self.num_var_);
        let mut xl_temp = Vector::new(self.num_var_);
        let mut xu_temp = Vector::new(self.num_var_);
        let mut slack_temp = Vector::new(self.num_constr_);
        let mut y_temp = Vector::new(self.num_constr_);
        let mut zl_temp = Vector::new(self.num_var_);
        let mut zu_temp = Vector::new(self.num_var_);
        self.dualize_back_interior_solution(
            x_solver,
            xl_solver,
            xu_solver,
            y_solver,
            zl_solver,
            zu_solver,
            &mut x_temp,
            &mut xl_temp,
            &mut xu_temp,
            &mut slack_temp,
            &mut y_temp,
            &mut zl_temp,
            &mut zu_temp,
        );
        self.scale_back_interior_solution(
            &mut x_temp,
            &mut xl_temp,
            &mut xu_temp,
            &mut slack_temp,
            &mut y_temp,
            &mut zl_temp,
            &mut zu_temp,
        );
        if let Some(out) = x_user {
            copy_all(&x_temp[..], out);
        }
        if let Some(out) = xl_user {
            copy_all(&xl_temp[..], out);
        }
        if let Some(out) = xu_user {
            copy_all(&xu_temp[..], out);
        }
        if let Some(out) = slack_user {
            copy_all(&slack_temp[..], out);
        }
        if let Some(out) = y_user {
            copy_all(&y_temp[..], out);
        }
        if let Some(out) = zl_user {
            copy_all(&zl_temp[..], out);
        }
        if let Some(out) = zu_user {
            copy_all(&zu_temp[..], out);
        }
    }

    /// HiGHS: ipx/model.cc::Model::EvaluateInteriorSolution (lines 283-383)
    pub fn evaluate_interior_solution(
        &self,
        x_solver: &Vector,
        xl_solver: &Vector,
        xu_solver: &Vector,
        y_solver: &Vector,
        zl_solver: &Vector,
        zu_solver: &Vector,
        info: &mut Info,
    ) {
        let m = self.rows();
        let n = self.cols();
        debug_assert_eq!(x_solver.size(), n + m);
        debug_assert_eq!(xl_solver.size(), n + m);
        debug_assert_eq!(xu_solver.size(), n + m);
        debug_assert_eq!(y_solver.size(), m);
        debug_assert_eq!(zl_solver.size(), n + m);
        debug_assert_eq!(zu_solver.size(), n + m);
        let nv = self.num_var_ as usize;
        let nc = self.num_constr_ as usize;

        // Build solution to scaled user model.
        let mut x = Vector::new(self.num_var_);
        let mut xl = Vector::new(self.num_var_);
        let mut xu = Vector::new(self.num_var_);
        let mut slack = Vector::new(self.num_constr_);
        let mut y = Vector::new(self.num_constr_);
        let mut zl = Vector::new(self.num_var_);
        let mut zu = Vector::new(self.num_var_);
        self.dualize_back_interior_solution(
            x_solver, xl_solver, xu_solver, y_solver, zl_solver, zu_solver, &mut x, &mut xl,
            &mut xu, &mut slack, &mut y, &mut zl, &mut zu,
        );

        // Build residuals for scaled model.
        // rl = lb-x+xl
        let mut rl = Vector::new(self.num_var_);
        for j in 0..nv {
            if self.scaled_lbuser_[j].is_finite() {
                rl[j] = self.scaled_lbuser_[j] - x[j] + xl[j];
            } else {
                debug_assert!(xl[j] == f64::INFINITY);
            }
        }
        // ru = ub-x-xu
        let mut ru = Vector::new(self.num_var_);
        for j in 0..nv {
            if self.scaled_ubuser_[j].is_finite() {
                ru[j] = self.scaled_ubuser_[j] - x[j] - xu[j];
            } else {
                debug_assert!(xu[j] == f64::INFINITY);
            }
        }
        // rb = rhs-slack-A*x
        // Add rhs at the end to avoid losing digits when x is huge.
        let mut rb = Vector::new(self.num_constr_);
        self.multiply_with_scaled_matrix(&x, -1.0, &mut rb, b'N');
        for i in 0..nc {
            rb[i] -= slack[i];
        }
        for i in 0..nc {
            rb[i] += self.scaled_rhs_[i];
        }
        // rc = obj-zl+zu-A'y
        // Add obj at the end to avoid losing digits when y, z are huge.
        let mut rc = Vector::new(self.num_var_);
        self.multiply_with_scaled_matrix(&y, -1.0, &mut rc, b'T');
        for j in 0..nv {
            rc[j] -= zl[j] - zu[j];
        }
        for j in 0..nv {
            rc[j] += self.scaled_obj_[j];
        }

        self.scale_back_residuals(&mut rb, &mut rc, &mut rl, &mut ru);
        let mut presidual = infnorm(&rb);
        presidual = cpp_max(presidual, infnorm(&rl));
        presidual = cpp_max(presidual, infnorm(&ru));
        let dresidual = infnorm(&rc);

        let pobjective = self.offset_ + dot(&self.scaled_obj_, &x);
        let mut dobjective = self.offset_ + dot(&self.scaled_rhs_, &y);
        for j in 0..nv {
            if self.scaled_lbuser_[j].is_finite() {
                dobjective += self.scaled_lbuser_[j] * zl[j];
            }
            if self.scaled_ubuser_[j].is_finite() {
                dobjective -= self.scaled_ubuser_[j] * zu[j];
            }
        }
        debug_assert!(dobjective.is_finite());
        let objective_gap =
            (pobjective - dobjective) / (1.0 + 0.5 * (pobjective + dobjective).abs());

        let mut complementarity = 0.0;
        for j in 0..nv {
            if self.scaled_lbuser_[j].is_finite() {
                complementarity += xl[j] * zl[j];
            }
            if self.scaled_ubuser_[j].is_finite() {
                complementarity += xu[j] * zu[j];
            }
        }
        for i in 0..nc {
            complementarity -= y[i] * slack[i];
        }

        // For computing the norms of the user variables, we have to scale back.
        self.scale_back_interior_solution(
            &mut x, &mut xl, &mut xu, &mut slack, &mut y, &mut zl, &mut zu,
        );

        info.abs_presidual = presidual;
        info.abs_dresidual = dresidual;
        info.rel_presidual = presidual / (1.0 + self.norm_rhs_);
        info.rel_dresidual = dresidual / (1.0 + self.norm_obj_);
        info.pobjval = pobjective;
        info.dobjval = dobjective;
        info.rel_objgap = objective_gap;
        info.complementarity = complementarity;
        info.normx = infnorm(&x);
        info.normy = infnorm(&y);
        info.normz = cpp_max(infnorm(&zl), infnorm(&zu));
    }

    /// HiGHS: ipx/model.cc::Model::PostsolveBasicSolution (lines 385-418)
    ///
    /// Each `None` output is not returned.
    pub fn postsolve_basic_solution(
        &self,
        x_solver: &Vector,
        y_solver: &Vector,
        z_solver: &Vector,
        basic_status_solver: &[Int],
        x_user: Option<&mut [f64]>,
        slack_user: Option<&mut [f64]>,
        y_user: Option<&mut [f64]>,
        z_user: Option<&mut [f64]>,
    ) {
        let m = self.rows();
        let n = self.cols();
        debug_assert_eq!(x_solver.size(), n + m);
        debug_assert_eq!(y_solver.size(), m);
        debug_assert_eq!(z_solver.size(), n + m);
        debug_assert_eq!(basic_status_solver.len() as Int, n + m);

        let mut x_temp = Vector::new(self.num_var_);
        let mut slack_temp = Vector::new(self.num_constr_);
        let mut y_temp = Vector::new(self.num_constr_);
        let mut z_temp = Vector::new(self.num_var_);
        let mut cbasis_temp: Vec<Int> = vec![0; self.num_constr_ as usize];
        let mut vbasis_temp: Vec<Int> = vec![0; self.num_var_ as usize];
        self.dualize_back_basic_solution(
            x_solver,
            y_solver,
            z_solver,
            &mut x_temp,
            &mut slack_temp,
            &mut y_temp,
            &mut z_temp,
        );
        self.dualize_back_basis(basic_status_solver, &mut cbasis_temp, &mut vbasis_temp);
        self.correct_scaled_basic_solution(
            &mut x_temp,
            &mut slack_temp,
            &mut y_temp,
            &mut z_temp,
            &cbasis_temp,
            &vbasis_temp,
        );
        self.scale_back_basic_solution(&mut x_temp, &mut slack_temp, &mut y_temp, &mut z_temp);
        if let Some(out) = x_user {
            copy_all(&x_temp[..], out);
        }
        if let Some(out) = slack_user {
            copy_all(&slack_temp[..], out);
        }
        if let Some(out) = y_user {
            copy_all(&y_temp[..], out);
        }
        if let Some(out) = z_user {
            copy_all(&z_temp[..], out);
        }
    }

    /// HiGHS: ipx/model.cc::Model::EvaluateBasicSolution (lines 420-481)
    pub fn evaluate_basic_solution(
        &self,
        x_solver: &Vector,
        y_solver: &Vector,
        z_solver: &Vector,
        basic_status_solver: &[Int],
        info: &mut Info,
    ) {
        let m = self.rows();
        let n = self.cols();
        debug_assert_eq!(x_solver.size(), n + m);
        debug_assert_eq!(y_solver.size(), m);
        debug_assert_eq!(z_solver.size(), n + m);
        debug_assert_eq!(basic_status_solver.len() as Int, n + m);
        let nv = self.num_var_ as usize;
        let nc = self.num_constr_ as usize;

        // Build basic solution to scaled user model.
        let mut x = Vector::new(self.num_var_);
        let mut slack = Vector::new(self.num_constr_);
        let mut y = Vector::new(self.num_constr_);
        let mut z = Vector::new(self.num_var_);
        let mut cbasis: Vec<Int> = vec![0; nc];
        let mut vbasis: Vec<Int> = vec![0; nv];
        self.dualize_back_basic_solution(
            x_solver, y_solver, z_solver, &mut x, &mut slack, &mut y, &mut z,
        );
        self.dualize_back_basis(basic_status_solver, &mut cbasis, &mut vbasis);
        self.correct_scaled_basic_solution(&mut x, &mut slack, &mut y, &mut z, &cbasis, &vbasis);
        let pobj = dot(&self.scaled_obj_, &x);

        // Build infeasibilities in scaled user model.
        let mut xinfeas = Vector::new(self.num_var_);
        let mut sinfeas = Vector::new(self.num_constr_);
        let mut yinfeas = Vector::new(self.num_constr_);
        let mut zinfeas = Vector::new(self.num_var_);
        for j in 0..nv {
            if x[j] < self.scaled_lbuser_[j] {
                xinfeas[j] = x[j] - self.scaled_lbuser_[j];
            }
            if x[j] > self.scaled_ubuser_[j] {
                xinfeas[j] = x[j] - self.scaled_ubuser_[j];
            }
            if vbasis[j] != IPX_nonbasic_lb && z[j] > 0.0 {
                zinfeas[j] = z[j];
            }
            if vbasis[j] != IPX_nonbasic_ub && z[j] < 0.0 {
                zinfeas[j] = z[j];
            }
        }
        for i in 0..nc {
            if self.constr_type_[i] == b'<' {
                if slack[i] < 0.0 {
                    sinfeas[i] = slack[i];
                }
                if y[i] > 0.0 {
                    yinfeas[i] = y[i];
                }
            }
            if self.constr_type_[i] == b'>' {
                if slack[i] > 0.0 {
                    sinfeas[i] = slack[i];
                }
                if y[i] < 0.0 {
                    yinfeas[i] = y[i];
                }
            }
        }

        // Scale back basic solution and infeasibilities.
        self.scale_back_basic_solution(&mut x, &mut slack, &mut y, &mut z);
        self.scale_back_basic_solution(&mut xinfeas, &mut sinfeas, &mut yinfeas, &mut zinfeas);

        info.primal_infeas = cpp_max(infnorm(&xinfeas), infnorm(&sinfeas));
        info.dual_infeas = cpp_max(infnorm(&zinfeas), infnorm(&yinfeas));
        info.objval = pobj;
    }

    /// HiGHS: ipx/model.cc::Model::PostsolveBasis (lines 483-497)
    ///
    /// Each `None` output is not returned.
    pub fn postsolve_basis(
        &self,
        basic_status_solver: &[Int],
        cbasis_user: Option<&mut [Int]>,
        vbasis_user: Option<&mut [Int]>,
    ) {
        let m = self.rows();
        let n = self.cols();
        debug_assert_eq!(basic_status_solver.len() as Int, n + m);

        let mut cbasis_temp: Vec<Int> = vec![0; self.num_constr_ as usize];
        let mut vbasis_temp: Vec<Int> = vec![0; self.num_var_ as usize];
        self.dualize_back_basis(basic_status_solver, &mut cbasis_temp, &mut vbasis_temp);
        self.scale_back_basis(&mut cbasis_temp, &mut vbasis_temp);
        if let Some(out) = cbasis_user {
            copy_all(&cbasis_temp[..], out);
        }
        if let Some(out) = vbasis_user {
            copy_all(&vbasis_temp[..], out);
        }
    }
}
