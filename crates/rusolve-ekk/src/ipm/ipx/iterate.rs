//! HiGHS source: ipm/ipx/iterate.h, iterate.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! C++ `mutable` members (the residuals and objectives that the `const`
//! method `Evaluate()` computes lazily) are `RefCell<Vector>` / `Cell<..>`;
//! the residual getters return `Ref<'_, Vector>` (drop it before anything
//! that can re-evaluate). Header inline members are implemented.

use super::model::Model;
use super::sparse_matrix::{dot_column, multiply_add};
use super::utils::{all_finite, dot, infnorm};
use super::*;
use crate::util::highs_utils::{cpp_max, cpp_min};
use std::cell::{Cell, Ref, RefCell};
use std::rc::Rc;

/// HiGHS: ipx/iterate.h::Iterate::State
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    fixed,
    free,
    barrier,
}

/// HiGHS: ipx/iterate.h::Iterate::StateDetail
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StateDetail {
    BARRIER_LB,
    BARRIER_UB,
    BARRIER_BOXED,
    FREE,
    FIXED,
    IMPLIED_LB,
    IMPLIED_UB,
    IMPLIED_EQ,
}

/// HiGHS: ipx/iterate.h::Iterate
pub struct Iterate {
    pub bounds_measure_: f64,
    pub costs_measure_: f64,
    model_: Rc<Model>,
    x_: Vector,
    xl_: Vector,
    xu_: Vector,
    y_: Vector,
    zl_: Vector,
    zu_: Vector,
    variable_state_: Vec<StateDetail>,
    rb_: RefCell<Vector>,
    rl_: RefCell<Vector>,
    ru_: RefCell<Vector>,
    rc_: RefCell<Vector>,
    pobjective_: Cell<f64>,
    dobjective_: Cell<f64>,
    presidual_: Cell<f64>,
    dresidual_: Cell<f64>,
    offset_: Cell<f64>,
    complementarity_: Cell<f64>,
    mu_: Cell<f64>,
    mu_min_: Cell<f64>,
    mu_max_: Cell<f64>,
    evaluated_: Cell<bool>,
    postprocessed_: bool,
    feasibility_tol_: f64,
    optimality_tol_: f64,
    start_crossover_tol_: f64,
    /// rusolve (HiGHS 1.16 backport, default off): IPX gap-offset flag, default false
    h16_gap_offset_: bool,
}

impl Iterate {
    pub const kBarrierMin: f64 = 1e-30;

    /// HiGHS: ipx/iterate.cc::Iterate::Iterate (lines 13-61)
    pub fn new(model: Rc<Model>) -> Self {
        let m = model.rows();
        let n = model.cols();
        let lb = model.lb();
        let ub = model.ub();
        let mut x_ = Vector::new(n + m);
        let mut xl_ = Vector::new(n + m);
        let mut xu_ = Vector::new(n + m);
        let y_ = Vector::new(m);
        let mut zl_ = Vector::new(n + m);
        let mut zu_ = Vector::new(n + m);
        let rb_ = Vector::new(m);
        let rl_ = Vector::new(n + m);
        let ru_ = Vector::new(n + m);
        let rc_ = Vector::new(n + m);
        let mut variable_state_ = vec![StateDetail::FREE; (n + m) as usize];

        for j in 0..n + m {
            let ju = j as usize;
            if lb[ju].is_finite() && ub[ju].is_finite() {
                variable_state_[ju] = StateDetail::BARRIER_BOXED;
                xl_[ju] = 1.0;
                xu_[ju] = 1.0;
                zl_[ju] = 1.0;
                zu_[ju] = 1.0;
            } else if lb[ju].is_finite() {
                variable_state_[ju] = StateDetail::BARRIER_LB;
                xl_[ju] = 1.0;
                xu_[ju] = f64::INFINITY;
                zl_[ju] = 1.0;
                zu_[ju] = 0.0;
            } else if ub[ju].is_finite() {
                variable_state_[ju] = StateDetail::BARRIER_UB;
                xl_[ju] = f64::INFINITY;
                xu_[ju] = 1.0;
                zl_[ju] = 0.0;
                zu_[ju] = 1.0;
            } else {
                variable_state_[ju] = StateDetail::FREE;
                xl_[ju] = f64::INFINITY;
                xu_[ju] = f64::INFINITY;
                zl_[ju] = 0.0;
                zu_[ju] = 0.0;
            }
        }
        let bounds_measure_ = 1.0 + model.norm_bounds();
        let costs_measure_ = 1.0 + model.norm_c();
        let mut it = Iterate {
            bounds_measure_,
            costs_measure_,
            model_: model,
            x_,
            xl_,
            xu_,
            y_,
            zl_,
            zu_,
            variable_state_,
            rb_: RefCell::new(rb_),
            rl_: RefCell::new(rl_),
            ru_: RefCell::new(ru_),
            rc_: RefCell::new(rc_),
            pobjective_: Cell::new(0.0),
            dobjective_: Cell::new(0.0),
            presidual_: Cell::new(0.0),
            dresidual_: Cell::new(0.0),
            offset_: Cell::new(0.0),
            complementarity_: Cell::new(0.0),
            mu_: Cell::new(0.0),
            mu_min_: Cell::new(0.0),
            mu_max_: Cell::new(0.0),
            evaluated_: Cell::new(false),
            postprocessed_: false,
            feasibility_tol_: 1e-6,
            optimality_tol_: 1e-8,
            start_crossover_tol_: -1.0,
            h16_gap_offset_: false,
        };
        if cfg!(debug_assertions) {
            it.assert_consistency();
        }
        it
    }

    /// HiGHS: ipx/iterate.cc::Iterate::Initialize (lines 63-94)
    pub fn initialize(
        &mut self,
        x: &Vector,
        xl: &Vector,
        xu: &Vector,
        y: &Vector,
        zl: &Vector,
        zu: &Vector,
    ) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        let lb = model.lb();
        let ub = model.ub();
        debug_assert!(x.size() == n + m);
        debug_assert!(xl.size() == n + m);
        debug_assert!(xu.size() == n + m);
        debug_assert!(y.size() == m);
        debug_assert!(zl.size() == n + m);
        debug_assert!(zu.size() == n + m);
        self.x_.assign(x);
        self.xl_.assign(xl);
        self.xu_.assign(xu);
        self.y_.assign(y);
        self.zl_.assign(zl);
        self.zu_.assign(zu);

        // Set variable statuses.
        for j in 0..n + m {
            let ju = j as usize;
            if lb[ju] == ub[ju] {
                self.variable_state_[ju] = StateDetail::BARRIER_BOXED;
            } else if lb[ju].is_finite() && ub[ju].is_finite() {
                self.variable_state_[ju] = StateDetail::BARRIER_BOXED;
            } else if lb[ju].is_finite() {
                self.variable_state_[ju] = StateDetail::BARRIER_LB;
            } else if ub[ju].is_finite() {
                self.variable_state_[ju] = StateDetail::BARRIER_UB;
            } else {
                self.variable_state_[ju] = StateDetail::FREE;
            }
        }
        if cfg!(debug_assertions) {
            self.assert_consistency();
        }
        self.evaluated_.set(false);
        self.postprocessed_ = false;
    }

    /// Test hook (not in HiGHS): `variable_state_` as the C++ enum
    /// values, for the golden replay (the oracle reads the private member
    /// with `#define private public`).
    #[doc(hidden)]
    pub fn oracle_state(&self) -> Vec<Int> {
        self.variable_state_.iter().map(|&s| s as Int).collect()
    }

    /// Test hook (not in HiGHS): the private `postprocessed_` flag.
    #[doc(hidden)]
    pub fn oracle_postprocessed(&self) -> bool {
        self.postprocessed_
    }

    /// Test hook (not in HiGHS): sets the private state of an iterate the
    /// oracle dumped (vectors, `variable_state_` as C++ enum values,
    /// `postprocessed_`), so that a golden replay can start from a state
    /// produced by the C++ IPM. Cached evaluations are invalidated.
    #[doc(hidden)]
    #[allow(clippy::too_many_arguments)]
    pub fn oracle_load(
        &mut self,
        x: &[f64],
        xl: &[f64],
        xu: &[f64],
        y: &[f64],
        zl: &[f64],
        zu: &[f64],
        state: &[Int],
        postprocessed: bool,
    ) {
        const STATES: [StateDetail; 8] = [
            StateDetail::BARRIER_LB,
            StateDetail::BARRIER_UB,
            StateDetail::BARRIER_BOXED,
            StateDetail::FREE,
            StateDetail::FIXED,
            StateDetail::IMPLIED_LB,
            StateDetail::IMPLIED_UB,
            StateDetail::IMPLIED_EQ,
        ];
        assert_eq!(x.len(), self.x_.len());
        assert_eq!(y.len(), self.y_.len());
        assert_eq!(state.len(), self.variable_state_.len());
        self.x_.0.copy_from_slice(x);
        self.xl_.0.copy_from_slice(xl);
        self.xu_.0.copy_from_slice(xu);
        self.y_.0.copy_from_slice(y);
        self.zl_.0.copy_from_slice(zl);
        self.zu_.0.copy_from_slice(zu);
        for (d, &s) in self.variable_state_.iter_mut().zip(state) {
            *d = STATES[s as usize];
        }
        self.postprocessed_ = postprocessed;
        self.evaluated_.set(false);
    }

    /// HiGHS: ipx/iterate.cc::Iterate::Update (lines 96-141)
    ///
    /// Every direction may be `nullptr` (starting_basis.cc:112).
    #[allow(clippy::too_many_arguments)]
    pub fn update(
        &mut self,
        sp: f64,
        dx: Option<&[f64]>,
        dxl: Option<&[f64]>,
        dxu: Option<&[f64]>,
        sd: f64,
        dy: Option<&[f64]>,
        dzl: Option<&[f64]>,
        dzu: Option<&[f64]>,
    ) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();

        if let Some(dx) = dx {
            for j in 0..n + m {
                if self.state_of(j) != State::fixed {
                    self.x_[j as usize] += sp * dx[j as usize];
                }
            }
        }
        if let Some(dxl) = dxl {
            for j in 0..n + m {
                if self.has_barrier_lb(j) {
                    let ju = j as usize;
                    self.xl_[ju] += sp * dxl[ju];
                    self.xl_[ju] = cpp_max(self.xl_[ju], Self::kBarrierMin);
                }
            }
        }
        if let Some(dxu) = dxu {
            for j in 0..n + m {
                if self.has_barrier_ub(j) {
                    let ju = j as usize;
                    self.xu_[ju] += sp * dxu[ju];
                    self.xu_[ju] = cpp_max(self.xu_[ju], Self::kBarrierMin);
                }
            }
        }
        if let Some(dy) = dy {
            for i in 0..m {
                let iu = i as usize;
                self.y_[iu] += sd * dy[iu];
            }
        }
        if let Some(dzl) = dzl {
            for j in 0..n + m {
                if self.has_barrier_lb(j) {
                    let ju = j as usize;
                    self.zl_[ju] += sd * dzl[ju];
                    self.zl_[ju] = cpp_max(self.zl_[ju], Self::kBarrierMin);
                }
            }
        }
        if let Some(dzu) = dzu {
            for j in 0..n + m {
                if self.has_barrier_ub(j) {
                    let ju = j as usize;
                    self.zu_[ju] += sd * dzu[ju];
                    self.zu_[ju] = cpp_max(self.zu_[ju], Self::kBarrierMin);
                }
            }
        }
        if cfg!(debug_assertions) {
            self.assert_consistency();
        }
        self.evaluated_.set(false);
    }

    pub fn model(&self) -> &Model {
        &self.model_
    }
    /// The shared handle (for objects that keep `const Model&`).
    pub fn model_rc(&self) -> Rc<Model> {
        Rc::clone(&self.model_)
    }
    pub fn x(&self) -> &Vector {
        &self.x_
    }
    pub fn xl(&self) -> &Vector {
        &self.xl_
    }
    pub fn xu(&self) -> &Vector {
        &self.xu_
    }
    pub fn y(&self) -> &Vector {
        &self.y_
    }
    pub fn zl(&self) -> &Vector {
        &self.zl_
    }
    pub fn zu(&self) -> &Vector {
        &self.zu_
    }
    /// HiGHS: ipx/iterate.h::Iterate::x(Int)
    pub fn x_at(&self, j: Int) -> f64 {
        self.x_[j as usize]
    }
    pub fn xl_at(&self, j: Int) -> f64 {
        self.xl_[j as usize]
    }
    pub fn xu_at(&self, j: Int) -> f64 {
        self.xu_[j as usize]
    }
    pub fn y_at(&self, j: Int) -> f64 {
        self.y_[j as usize]
    }
    pub fn zl_at(&self, j: Int) -> f64 {
        self.zl_[j as usize]
    }
    pub fn zu_at(&self, j: Int) -> f64 {
        self.zu_[j as usize]
    }

    /// HiGHS: ipx/iterate.cc::Iterate::rb (line 143)
    pub fn rb(&self) -> Ref<'_, Vector> {
        self.evaluate();
        self.rb_.borrow()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::rl (line 144)
    pub fn rl(&self) -> Ref<'_, Vector> {
        self.evaluate();
        self.rl_.borrow()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::ru (line 145)
    pub fn ru(&self) -> Ref<'_, Vector> {
        self.evaluate();
        self.ru_.borrow()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::rc (line 146)
    pub fn rc(&self) -> Ref<'_, Vector> {
        self.evaluate();
        self.rc_.borrow()
    }

    /// HiGHS: ipx/iterate.h::Iterate::StateOf (inline)
    pub fn state_of(&self, j: Int) -> State {
        match self.variable_state_[j as usize] {
            StateDetail::FIXED => State::fixed,
            StateDetail::FREE
            | StateDetail::IMPLIED_LB
            | StateDetail::IMPLIED_UB
            | StateDetail::IMPLIED_EQ => State::free,
            _ => State::barrier,
        }
    }
    /// HiGHS: ipx/iterate.h::Iterate::has_barrier_lb (inline)
    pub fn has_barrier_lb(&self, j: Int) -> bool {
        let s = self.variable_state_[j as usize];
        s == StateDetail::BARRIER_LB || s == StateDetail::BARRIER_BOXED
    }
    /// HiGHS: ipx/iterate.h::Iterate::has_barrier_ub (inline)
    pub fn has_barrier_ub(&self, j: Int) -> bool {
        let s = self.variable_state_[j as usize];
        s == StateDetail::BARRIER_UB || s == StateDetail::BARRIER_BOXED
    }
    /// HiGHS: ipx/iterate.h::Iterate::is_implied (inline)
    pub fn is_implied(&self, j: Int) -> bool {
        let s = self.variable_state_[j as usize];
        s == StateDetail::IMPLIED_LB || s == StateDetail::IMPLIED_UB || s == StateDetail::IMPLIED_EQ
    }

    /// HiGHS: ipx/iterate.cc::Iterate::make_fixed(Int) (lines 148-155)
    pub fn make_fixed(&mut self, j: Int) {
        let ju = j as usize;
        self.xl_[ju] = 0.0;
        self.xu_[ju] = 0.0;
        self.zl_[ju] = 0.0;
        self.zu_[ju] = 0.0;
        self.variable_state_[ju] = StateDetail::FIXED;
        self.evaluated_.set(false);
    }
    /// HiGHS: ipx/iterate.cc::Iterate::make_fixed(Int, double) (lines 157-160)
    pub fn make_fixed_value(&mut self, j: Int, value: f64) {
        self.x_[j as usize] = value;
        self.make_fixed(j);
    }
    /// HiGHS: ipx/iterate.cc::Iterate::make_implied_lb (lines 162-167)
    pub fn make_implied_lb(&mut self, j: Int) {
        let ju = j as usize;
        self.xl_[ju] = f64::INFINITY;
        self.xu_[ju] = f64::INFINITY;
        self.variable_state_[ju] = StateDetail::IMPLIED_LB;
        self.evaluated_.set(false);
    }
    /// HiGHS: ipx/iterate.cc::Iterate::make_implied_ub (lines 169-174)
    pub fn make_implied_ub(&mut self, j: Int) {
        let ju = j as usize;
        self.xl_[ju] = f64::INFINITY;
        self.xu_[ju] = f64::INFINITY;
        self.variable_state_[ju] = StateDetail::IMPLIED_UB;
        self.evaluated_.set(false);
    }
    /// HiGHS: ipx/iterate.cc::Iterate::make_implied_eq (lines 176-183)
    pub fn make_implied_eq(&mut self, j: Int) {
        let ju = j as usize;
        self.xl_[ju] = f64::INFINITY;
        self.xu_[ju] = f64::INFINITY;
        self.zl_[ju] = 0.0;
        self.zu_[ju] = 0.0;
        self.variable_state_[ju] = StateDetail::IMPLIED_EQ;
        self.evaluated_.set(false);
    }
    /// HiGHS: ipx/iterate.cc::Iterate::ScalingFactor (lines 185-200)
    pub fn scaling_factor(&self, j: Int) -> f64 {
        match self.state_of(j) {
            State::fixed => 0.0,
            State::free => f64::INFINITY,
            _ => {
                let ju = j as usize;
                debug_assert!(self.xl_[ju] > 0.0);
                debug_assert!(self.xu_[ju] > 0.0);
                let g = self.zl_[ju] / self.xl_[ju] + self.zu_[ju] / self.xu_[ju];
                let d = 1.0 / g.sqrt();
                debug_assert!(d.is_finite());
                debug_assert!(d > 0.0);
                d
            }
        }
    }

    /// HiGHS: ipx/iterate.cc::Iterate::pobjective (line 202)
    pub fn pobjective(&self) -> f64 {
        self.evaluate();
        self.pobjective_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::dobjective (line 203)
    pub fn dobjective(&self) -> f64 {
        self.evaluate();
        self.dobjective_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::pobjective_after_postproc (lines 205-208)
    pub fn pobjective_after_postproc(&self) -> f64 {
        self.evaluate();
        self.pobjective_.get() + self.offset_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::dobjective_after_postproc (lines 210-213)
    pub fn dobjective_after_postproc(&self) -> f64 {
        self.evaluate();
        self.dobjective_.get() + self.offset_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::presidual (line 215)
    pub fn presidual(&self) -> f64 {
        self.evaluate();
        self.presidual_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::dresidual (line 216)
    pub fn dresidual(&self) -> f64 {
        self.evaluate();
        self.dresidual_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::complementarity (line 218)
    pub fn complementarity(&self) -> f64 {
        self.evaluate();
        self.complementarity_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::mu (line 219)
    pub fn mu(&self) -> f64 {
        self.evaluate();
        self.mu_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::mu_min (line 220)
    pub fn mu_min(&self) -> f64 {
        self.evaluate();
        self.mu_min_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::mu_max (line 221)
    pub fn mu_max(&self) -> f64 {
        self.evaluate();
        self.mu_max_.get()
    }
    /// HiGHS: ipx/iterate.cc::Iterate::feasible (lines 223-229)
    pub fn feasible(&self) -> bool {
        self.evaluate();
        let primal_feasible = self.presidual_.get() <= self.feasibility_tol_ * self.bounds_measure_;
        let dual_feasible = self.dresidual_.get() <= self.feasibility_tol_ * self.costs_measure_;
        primal_feasible && dual_feasible
    }
    /// HiGHS: ipx/iterate.cc::Iterate::optimal (lines 231-241)
    /// HiGHS 1.16 (#3427): ipm/ipx/iterate.cc::Iterate::optimal (lines 231-241)
    pub fn optimal(&self) -> bool {
        self.evaluate();
        let pobj = self.pobjective_after_postproc();
        let dobj = self.dobjective_after_postproc();
        // HiGHS 1.16 #3427 (h16_ipm_gap_offset): ipm/ipx/iterate.cc line 235 (1.16.0)
        let ave_obj = if self.h16_gap_offset_ {
            0.5 * (pobj + dobj - 2.0 * self.offset_.get())
        } else {
            0.5 * (pobj + dobj)
        };
        let gap = pobj - dobj;
        let abs_gap = gap.abs();
        let obj_measure = 1.0 + ave_obj.abs();
        abs_gap <= self.optimality_tol_ * obj_measure
    }
    /// HiGHS: ipx/iterate.cc::Iterate::term_crit_reached (lines 243-254)
    pub fn term_crit_reached(&self) -> bool {
        if self.feasible() && self.optimal() {
            if self.start_crossover_tol_ <= 0.0 {
                return true;
            }
            let mut pres = 0.0;
            let mut dres = 0.0;
            self.residuals_from_dropping(&mut pres, &mut dres);
            let model = &*self.model_;
            if pres <= self.start_crossover_tol_ * (1.0 + model.norm_bounds())
                && dres <= self.start_crossover_tol_ * (1.0 + model.norm_c())
            {
                return true;
            }
        }
        false
    }

    pub fn feasibility_tol(&self) -> f64 {
        self.feasibility_tol_
    }
    pub fn optimality_tol(&self) -> f64 {
        self.optimality_tol_
    }
    pub fn start_crossover_tol(&self) -> f64 {
        self.start_crossover_tol_
    }
    /// HiGHS: ipx/iterate.h::Iterate::feasibility_tol(double)
    pub fn set_feasibility_tol(&mut self, new_tol: f64) {
        self.feasibility_tol_ = new_tol;
    }
    /// HiGHS: ipx/iterate.h::Iterate::optimality_tol(double)
    pub fn set_optimality_tol(&mut self, new_tol: f64) {
        self.optimality_tol_ = new_tol;
    }
    /// HiGHS: ipx/iterate.h::Iterate::start_crossover_tol(double)
    pub fn set_start_crossover_tol(&mut self, new_tol: f64) {
        self.start_crossover_tol_ = new_tol;
    }
    /// rusolve (HiGHS 1.16 backport, default off): setter for `h16_gap_offset_`
    pub fn set_h16_gap_offset(&mut self, value: bool) {
        self.h16_gap_offset_ = value;
    }
    /// rusolve (HiGHS 1.16 backport, default off): getter for `h16_gap_offset_`
    pub fn h16_gap_offset(&self) -> bool {
        self.h16_gap_offset_
    }

    /// HiGHS: ipx/iterate.cc::Iterate::Postprocess (lines 256-319)
    #[allow(clippy::needless_range_loop)]
    pub fn postprocess(&mut self) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        let c = model.c();
        let lb = model.lb();
        let ub = model.ub();
        let ai = model.AI();

        // For fixed variables compute xl[j] and xu[j] from x[j]. If the lower and
        // upper bound are equal, set zl[j] or zu[j] such that the variable is dual
        // feasible. Otherwise leave them zero.
        for j in 0..n + m {
            let ju = j as usize;
            if self.state_of(j) == State::fixed {
                self.xl_[ju] = self.x_[ju] - lb[ju];
                self.xu_[ju] = ub[ju] - self.x_[ju];
                debug_assert!(self.zl_[ju] == 0.0);
                debug_assert!(self.zu_[ju] == 0.0);
                if lb[ju] == ub[ju] {
                    let z = c[ju] - dot_column(ai, j, &self.y_);
                    if z >= 0.0 {
                        self.zl_[ju] = z;
                    } else {
                        self.zu_[ju] = -z;
                    }
                }
            }
        }
        // For implied variables set x[j] to the bound at which it was implied and
        // compute zl[j] or zu[j]. If the variable was implied at both bounds,
        // choose between zl and zu depending on sign.
        for j in 0..n + m {
            let ju = j as usize;
            if self.is_implied(j) {
                let z = c[ju] - dot_column(ai, j, &self.y_);
                match self.variable_state_[ju] {
                    StateDetail::IMPLIED_EQ => {
                        debug_assert!(lb[ju] == ub[ju]);
                        if z >= 0.0 {
                            self.zl_[ju] = z;
                            self.zu_[ju] = 0.0;
                        } else {
                            self.zl_[ju] = 0.0;
                            self.zu_[ju] = -z;
                        }
                        self.x_[ju] = lb[ju];
                    }
                    StateDetail::IMPLIED_LB => {
                        self.zl_[ju] = z;
                        self.zu_[ju] = 0.0;
                        self.x_[ju] = lb[ju];
                    }
                    StateDetail::IMPLIED_UB => {
                        self.zl_[ju] = 0.0;
                        self.zu_[ju] = -z;
                        self.x_[ju] = ub[ju];
                    }
                    _ => debug_assert!(false),
                }
                self.xl_[ju] = self.x_[ju] - lb[ju];
                self.xu_[ju] = ub[ju] - self.x_[ju];
            }
        }
        self.postprocessed_ = true;
        self.evaluated_.set(false);
    }
    /// HiGHS: ipx/iterate.cc::Iterate::EvaluatePostsolved (lines 321-323)
    pub fn evaluate_postsolved(&self, info: &mut Info) {
        self.model_.evaluate_interior_solution(
            &self.x_, &self.xl_, &self.xu_, &self.y_, &self.zl_, &self.zu_, info,
        );
    }
    /// HiGHS: ipx/iterate.cc::Iterate::DropToComplementarity (lines 325-398)
    #[allow(clippy::needless_range_loop)]
    pub fn drop_to_complementarity(&self, x: &mut Vector, y: &mut Vector, z: &mut Vector) {
        let model: &Model = &self.model_;
        let m = model.rows();
        let n = model.cols();
        let lb = model.lb();
        let ub = model.ub();
        debug_assert!(self.postprocessed_);

        debug_assert!(x.size() == n + m);
        debug_assert!(y.size() == m);
        debug_assert!(z.size() == n + m);
        y.assign(&self.y_);
        for j in 0..n + m {
            let ju = j as usize;
            let xlj = self.xl_[ju];
            let xuj = self.xu_[ju];
            let zlj = self.zl_[ju];
            let zuj = self.zu_[ju];
            debug_assert!(xlj >= 0.0);
            debug_assert!(xuj >= 0.0);
            debug_assert!(zlj >= 0.0);
            debug_assert!(zuj >= 0.0);
            let mut xj = self.x_[ju];
            xj = cpp_max(xj, lb[ju]);
            xj = cpp_min(xj, ub[ju]);

            if lb[ju] == ub[ju] {
                // fixed variable
                x[ju] = lb[ju];
                z[ju] = zlj - zuj;
            } else if lb[ju].is_finite() && ub[ju].is_finite() {
                // boxed variable
                if zlj * xuj >= zuj * xlj {
                    // either active at lower bound or inactive
                    if zlj >= xlj {
                        x[ju] = lb[ju];
                        z[ju] = cpp_max(0.0, zlj - zuj);
                    } else {
                        x[ju] = xj;
                        z[ju] = 0.0;
                    }
                } else {
                    // either active at upper bound or inactive
                    if zuj >= xuj {
                        x[ju] = ub[ju];
                        z[ju] = cpp_min(0.0, zlj - zuj);
                    } else {
                        x[ju] = xj;
                        z[ju] = 0.0;
                    }
                }
            } else if lb[ju].is_finite() {
                // lower bound only
                if zlj >= xlj {
                    x[ju] = lb[ju];
                    z[ju] = cpp_max(0.0, zlj - zuj);
                } else {
                    x[ju] = xj;
                    z[ju] = 0.0;
                }
            } else if ub[ju].is_finite() {
                // upper bound only
                if zuj >= xuj {
                    x[ju] = ub[ju];
                    z[ju] = cpp_min(0.0, zlj - zuj);
                } else {
                    x[ju] = xj;
                    z[ju] = 0.0;
                }
            } else {
                // free variable
                x[ju] = xj;
                z[ju] = 0.0;
            }
        }
    }
    /// HiGHS: ipx/iterate.cc::Iterate::ResidualsFromDropping (lines 400-450)
    #[allow(clippy::needless_range_loop)]
    fn residuals_from_dropping(&self, pres: &mut f64, dres: &mut f64) {
        let model: &Model = &self.model_;
        let m = model.rows();
        let n = model.cols();
        let ai = model.AI();
        let lb = model.lb();
        let ub = model.ub();
        let mut presmax = 0.0;
        let mut dresmax = 0.0;

        for j in 0..n + m {
            let ju = j as usize;
            let mut xdrop = 0.0; // xnew = xold - xdrop
            let mut zdrop = 0.0;
            let xl = self.xl_[ju];
            let xu = self.xu_[ju];
            let zl = self.zl_[ju];
            let zu = self.zu_[ju];
            let x = self.x_[ju];
            match self.variable_state_[ju] {
                StateDetail::BARRIER_LB => {
                    if zl >= xl {
                        xdrop = x - lb[ju]; // active at lower bound
                    } else {
                        zdrop = zl - zu; // inactive
                    }
                }
                StateDetail::BARRIER_UB => {
                    if zu >= xu {
                        xdrop = x - ub[ju]; // active at upper bound
                    } else {
                        zdrop = zl - zu; // inactive
                    }
                }
                StateDetail::BARRIER_BOXED => {
                    if zl / xl >= zu / xu {
                        if zl >= xl {
                            xdrop = x - lb[ju]; // active at lower bound
                        } else {
                            zdrop = zl - zu; // inactive
                        }
                    } else if zu >= xu {
                        xdrop = x - ub[ju]; // active at upper bound
                    } else {
                        zdrop = zl - zu; // inactive
                    }
                }
                _ => {}
            }
            let mut amax = 0.0;
            for p in ai.begin(j)..ai.begin(j + 1) {
                amax = cpp_max(amax, ai.value(p).abs());
            }
            presmax = cpp_max(presmax, xdrop.abs() * amax);
            dresmax = cpp_max(dresmax, zdrop.abs());
        }
        *pres = presmax;
        *dres = dresmax;
    }
    /// HiGHS: ipx/iterate.cc::Iterate::assert_consistency (lines 452-527)
    ///
    /// Body is `assert`s only: run under `if cfg!(debug_assertions)`.
    fn assert_consistency(&mut self) {
        let model: &Model = &self.model_;
        let m = model.rows();
        let n = model.cols();
        let lb = model.lb();
        let ub = model.ub();
        debug_assert!(all_finite(&self.x_));
        debug_assert!(all_finite(&self.y_));
        debug_assert!(all_finite(&self.zl_));
        debug_assert!(all_finite(&self.zu_));

        for j in 0..n + m {
            let ju = j as usize;
            let (l, u) = (lb[ju], ub[ju]);
            let (xl, xu) = (self.xl_[ju], self.xu_[ju]);
            let (zl, zu) = (self.zl_[ju], self.zu_[ju]);
            match self.variable_state_[ju] {
                StateDetail::BARRIER_LB => {
                    debug_assert!(l.is_finite());
                    debug_assert!(u.is_infinite());
                    debug_assert!(xl.is_finite() && xl > 0.0);
                    debug_assert!(xu.is_infinite());
                    debug_assert!(zl > 0.0);
                    debug_assert!(zu == 0.0);
                }
                StateDetail::BARRIER_UB => {
                    debug_assert!(l.is_infinite());
                    debug_assert!(u.is_finite());
                    debug_assert!(xl.is_infinite());
                    debug_assert!(xu.is_finite() && xu > 0.0);
                    debug_assert!(zl == 0.0);
                    debug_assert!(zu > 0.0);
                }
                StateDetail::BARRIER_BOXED => {
                    debug_assert!(l.is_finite());
                    debug_assert!(u.is_finite());
                    debug_assert!(xl.is_finite() && xl > 0.0);
                    debug_assert!(xu.is_finite() && xu > 0.0);
                    debug_assert!(zl > 0.0);
                    debug_assert!(zu > 0.0);
                }
                StateDetail::FREE => {
                    debug_assert!(l.is_infinite());
                    debug_assert!(u.is_infinite());
                    debug_assert!(xl.is_infinite());
                    debug_assert!(xu.is_infinite());
                    debug_assert!(zl == 0.0);
                    debug_assert!(zu == 0.0);
                }
                StateDetail::FIXED => {
                    debug_assert!(xl == 0.0);
                    debug_assert!(xu == 0.0);
                    debug_assert!(zl == 0.0);
                    debug_assert!(zu == 0.0);
                }
                StateDetail::IMPLIED_LB => {
                    debug_assert!(l.is_finite());
                    debug_assert!(xl.is_infinite());
                    debug_assert!(xu.is_infinite());
                    debug_assert!(zl >= 0.0);
                    debug_assert!(zu >= 0.0);
                }
                StateDetail::IMPLIED_UB => {
                    debug_assert!(u.is_finite());
                    debug_assert!(xl.is_infinite());
                    debug_assert!(xu.is_infinite());
                    debug_assert!(zl >= 0.0);
                    debug_assert!(zu >= 0.0);
                }
                StateDetail::IMPLIED_EQ => {
                    debug_assert!(l == u);
                    debug_assert!(xl.is_infinite());
                    debug_assert!(xu.is_infinite());
                    debug_assert!(zl == 0.0);
                    debug_assert!(zu == 0.0);
                }
            }
        }
    }
    /// HiGHS: ipx/iterate.cc::Iterate::Evaluate (lines 529-536)
    fn evaluate(&self) {
        if !self.evaluated_.get() {
            self.compute_residuals();
            self.compute_objectives();
            self.compute_complementarity();
            self.evaluated_.set(true);
        }
    }
    /// HiGHS: ipx/iterate.cc::Iterate::ComputeResiduals (lines 538-599)
    #[allow(clippy::needless_range_loop)]
    fn compute_residuals(&self) {
        let model: &Model = &self.model_;
        let m = model.rows();
        let n = model.cols();
        let lb = model.lb();
        let ub = model.ub();
        let ai = model.AI();

        // Primal residual: rb = b-AI*x.
        let mut rb = self.rb_.borrow_mut();
        rb.assign(model.b());
        multiply_add(ai, &self.x_, -1.0, &mut rb, b'N');

        // Dual residual: rc = c-AI'y-zl+zu. If the iterate has not been
        // postprocessed, then the dual residual for fixed variables is zero
        // because these variables are treated as non-existent by the IPM.
        let mut rc = self.rc_.borrow_mut();
        let c = model.c();
        for j in 0..n + m {
            let ju = j as usize;
            rc[ju] = (c[ju] - self.zl_[ju]) + self.zu_[ju];
        }
        multiply_add(ai, &self.y_, -1.0, &mut rc, b'T');
        if !self.postprocessed_ {
            for j in 0..n + m {
                if self.state_of(j) == State::fixed {
                    rc[j as usize] = 0.0;
                }
            }
        }

        // Bound residuals: rl = lb-x+xl and ru = ub-x-xu. If a variable has no
        // lower or upper bound, then the corresponding bound residual is zero.
        // Also, if the variable is implied (i.e. treated as free by the IPM) or
        // fixed (i.e. treated as non-existent by the IPM), the bound residual is
        // zero.
        //
        // Notice that postprocessing sets xl and xu for fixed/implied variables
        // such that rl and ru are zero; hence the residuals computed here are also
        // correct if the iterate was postprocessed.
        let mut rl = self.rl_.borrow_mut();
        let mut ru = self.ru_.borrow_mut();
        for j in 0..n + m {
            let ju = j as usize;
            rl[ju] = if self.has_barrier_lb(j) {
                lb[ju] - self.x_[ju] + self.xl_[ju]
            } else {
                0.0
            };
            ru[ju] = if self.has_barrier_ub(j) {
                ub[ju] - self.x_[ju] - self.xu_[ju]
            } else {
                0.0
            };
        }

        debug_assert!(all_finite(&rb));
        debug_assert!(all_finite(&rc));
        debug_assert!(all_finite(&rl));
        debug_assert!(all_finite(&ru));

        let mut presidual = infnorm(&rb);
        let dresidual = infnorm(&rc);
        presidual = cpp_max(presidual, infnorm(&rl));
        presidual = cpp_max(presidual, infnorm(&ru));
        self.presidual_.set(presidual);
        self.dresidual_.set(dresidual);
    }
    /// HiGHS: ipx/iterate.cc::Iterate::ComputeObjectives (lines 601-651)
    #[allow(clippy::needless_range_loop)]
    fn compute_objectives(&self) {
        let model: &Model = &self.model_;
        let m = model.rows();
        let n = model.cols();
        let b = model.b();
        let c = model.c();
        let lb = model.lb();
        let ub = model.ub();
        let ai = model.AI();

        if self.postprocessed_ {
            // Compute objective values as defined for the LP model.
            let mut dobjective = model.offset() + dot(b, &self.y_);
            let pobjective = model.offset() + dot(c, &self.x_);
            for j in 0..n + m {
                let ju = j as usize;
                if lb[ju].is_finite() {
                    dobjective += lb[ju] * self.zl_[ju];
                }
                if ub[ju].is_finite() {
                    dobjective -= ub[ju] * self.zu_[ju];
                }
            }
            self.offset_.set(0.0);
            self.pobjective_.set(pobjective);
            self.dobjective_.set(dobjective);
        } else {
            // Compute objective values for the LP that is solved at the very moment
            // (after fixing and implying variables). The offset is such that
            // pobjective_ + offset_ is the primal objective after postprocessing.
            let mut offset = 0.0;
            let mut pobjective = model.offset();
            for j in 0..n + m {
                let ju = j as usize;
                if self.state_of(j) != State::fixed {
                    pobjective += c[ju] * self.x_[ju];
                } else {
                    offset += c[ju] * self.x_[ju];
                }
                if self.is_implied(j) {
                    // At the moment, we are solving an LP with the cost coefficient
                    // for variable j decreased by zl[j]-zu[j].
                    pobjective -= (self.zl_[ju] - self.zu_[ju]) * self.x_[ju];
                    offset += (self.zl_[ju] - self.zu_[ju]) * self.x_[ju];
                }
            }
            let mut dobjective = model.offset() + dot(b, &self.y_);
            for j in 0..n + m {
                let ju = j as usize;
                if self.has_barrier_lb(j) {
                    dobjective += lb[ju] * self.zl_[ju];
                }
                if self.has_barrier_ub(j) {
                    dobjective -= ub[ju] * self.zu_[ju];
                }
                if self.state_of(j) == State::fixed {
                    // At the moment, we are solving the LP without variable j,
                    // but with the RHS decreased by AI[:,j]*x[j].
                    dobjective -= self.x_[ju] * dot_column(ai, j, &self.y_);
                }
            }
            self.offset_.set(offset);
            self.pobjective_.set(pobjective);
            self.dobjective_.set(dobjective);
        }
    }
    /// HiGHS: ipx/iterate.cc::Iterate::ComputeComplementarity (lines 653-681)
    #[allow(clippy::needless_range_loop)]
    fn compute_complementarity(&self) {
        let m = self.model_.rows();
        let n = self.model_.cols();

        let mut complementarity = 0.0;
        let mut mu_min = f64::INFINITY;
        let mut mu_max = 0.0;
        let mut num_finite: Int = 0;
        for j in 0..n + m {
            let ju = j as usize;
            if self.has_barrier_lb(j) {
                complementarity += self.xl_[ju] * self.zl_[ju];
                mu_min = cpp_min(mu_min, self.xl_[ju] * self.zl_[ju]);
                mu_max = cpp_max(mu_max, self.xl_[ju] * self.zl_[ju]);
                num_finite += 1;
            }
        }
        for j in 0..n + m {
            let ju = j as usize;
            if self.has_barrier_ub(j) {
                complementarity += self.xu_[ju] * self.zu_[ju];
                mu_min = cpp_min(mu_min, self.xu_[ju] * self.zu_[ju]);
                mu_max = cpp_max(mu_max, self.xu_[ju] * self.zu_[ju]);
                num_finite += 1;
            }
        }
        self.complementarity_.set(complementarity);
        self.mu_max_.set(mu_max);
        if num_finite > 0 {
            self.mu_min_.set(mu_min);
            self.mu_.set(complementarity / num_finite as f64);
        } else {
            self.mu_min_.set(0.0);
            self.mu_.set(0.0);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iterate_constructor_defaults() {
        // Empty model: the constructor defaults (C++ member initialisers and
        // the measures `1 + norm`) without touching any model data.
        let it = Iterate::new(Rc::new(Model::new()));
        assert_eq!(it.bounds_measure_, 1.0);
        assert_eq!(it.costs_measure_, 1.0);
        assert_eq!(it.feasibility_tol(), 1e-6);
        assert_eq!(it.optimality_tol(), 1e-8);
        assert_eq!(it.start_crossover_tol(), -1.0);
    }

    /// Iterate with a large objective offset (1e6), a true objective near 1,
    /// and a primal-dual objective gap of 1e-3. Returns the iterate with its
    /// objectives set and evaluation marked done, so `optimal()` reads them.
    fn h16_gap_offset_iterate(on: bool) -> Iterate {
        let mut it = Iterate::new(Rc::new(Model::new()));
        it.offset_.set(1e6);
        it.pobjective_.set(1.0);
        it.dobjective_.set(1.0 - 1e-3);
        it.evaluated_.set(true);
        it.set_h16_gap_offset(on);
        it
    }

    #[test]
    fn h16_ipm_gap_offset_ipx_optimal_off() {
        // Flag off: 1.15.1 measure 0.5*(pobj+dobj) ~ 1e6, so the 1e-3 gap is
        // within optimality_tol * 1e6 = 1e-2 and the iterate is optimal.
        let it = h16_gap_offset_iterate(false);
        assert!(it.optimal());
    }

    #[test]
    fn h16_ipm_gap_offset_ipx_optimal_on() {
        // Flag on (HiGHS 1.16 iterate.cc line 235): the average objective is
        // offset-free, ~1, so the same 1e-3 gap is far above 1e-8 * 2 and the
        // iterate is not optimal yet.
        let it = h16_gap_offset_iterate(true);
        assert!(!it.optimal());
        // A gap of 1e-9 still terminates with the flag on.
        let it = h16_gap_offset_iterate(true);
        it.dobjective_.set(1.0 - 1e-9);
        assert!(it.optimal());
    }
}
