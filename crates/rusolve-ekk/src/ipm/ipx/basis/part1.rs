//! HiGHS source: ipm/ipx/basis.cc (lines 17-482)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::super::basiclu_kernel::BasicLuKernel;
use super::super::basiclu_wrapper::BasicLu;
use super::super::control::{
    fixed, format_int, format_str, ostream_double, sci2, textline, time_fmt,
};
use super::super::forrest_tomlin::ForrestTomlin;
use super::super::indexed_vector::IndexedVector;
use super::super::lu_factorization::LuFactorization;
use super::super::power_method::power_method;
use super::super::sparse_matrix::{SparseMatrix, dot_column, scatter_column};
use super::super::symbolic_invert::symbolic_invert;
use super::super::timer::Timer;
use super::*;
use crate::ipm::basiclu::lu_internal::lu_time;

impl Basis {
    /// HiGHS: ipx/basis.cc::Basis::Basis (lines 17-31)
    pub fn new(control: Rc<Control>, model: Rc<Model>) -> Self {
        let m = model.rows();
        let n = model.cols();
        let lu: Box<dyn LuUpdate> = if control.lu_kernel() <= 0 {
            Box::new(BasicLu::new(Rc::clone(&control), m))
        } else {
            let kernel: Box<dyn LuFactorization> = Box::new(BasicLuKernel::default());
            Box::new(ForrestTomlin::new(Rc::clone(&control), m, kernel))
        };
        let mut basis = Basis {
            control_: control,
            model_: model,
            basis_: vec![0; m as usize],
            map2basis_: vec![0; (n + m) as usize],
            lu_: RefCell::new(lu),
            factorization_is_fresh_: false,
            num_factorizations_: 0,
            num_updates_: 0,
            num_ftran_: 0,
            num_btran_: 0,
            num_ftran_sparse_: 0,
            num_btran_sparse_: 0,
            time_ftran_: 0.0,
            time_btran_: 0.0,
            time_update_: 0.0,
            time_factorize_: 0.0,
            fill_factors_: Vec::new(),
            sum_ftran_density_: 0.0,
            sum_btran_density_: 0.0,
        };
        let pivottol = basis.control_.lu_pivottol();
        basis.lu_.borrow_mut().set_pivottol(pivottol);
        basis.set_to_slack_basis();
        basis
    }

    /// HiGHS: ipx/basis.cc::Basis::FixNonbasicVariable (lines 33-39)
    pub fn fix_nonbasic_variable(&mut self, j: Int) {
        if self.status_of(j) == NONBASIC_FIXED {
            return;
        }
        debug_assert!(self.status_of(j) == NONBASIC);
        debug_assert!(self.map2basis_[j as usize] == -1);
        self.map2basis_[j as usize] = -2;
    }

    /// HiGHS: ipx/basis.cc::Basis::FreeBasicVariable (lines 41-48)
    pub fn free_basic_variable(&mut self, j: Int) {
        let m = self.model_.rows();
        if self.status_of(j) == BASIC_FREE {
            return;
        }
        debug_assert!(self.status_of(j) == BASIC);
        debug_assert!(self.map2basis_[j as usize] >= 0 && self.map2basis_[j as usize] < m);
        self.map2basis_[j as usize] += m;
    }

    /// HiGHS: ipx/basis.cc::Basis::UnfixVariables (lines 50-56)
    pub fn unfix_variables(&mut self) {
        let m = self.model_.rows();
        let n = self.model_.cols();
        for j in 0..n + m {
            if self.map2basis_[j as usize] == -2 {
                self.map2basis_[j as usize] = -1;
            }
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::UnfreeVariables (lines 58-64)
    pub fn unfree_variables(&mut self) {
        let m = self.model_.rows();
        let n = self.model_.cols();
        for j in 0..n + m {
            if self.map2basis_[j as usize] >= m {
                self.map2basis_[j as usize] -= m;
            }
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::SetToSlackBasis (lines 66-78)
    pub fn set_to_slack_basis(&mut self) {
        let m = self.model_.rows();
        let n = self.model_.cols();
        for i in 0..m {
            self.basis_[i as usize] = n + i;
        }
        for j in 0..n {
            self.map2basis_[j as usize] = -1;
        }
        for i in 0..m {
            self.map2basis_[(n + i) as usize] = i;
        }
        let allow_timeout = false;
        let err = self.factorize_with(allow_timeout);
        // factorization of slack basis cannot fail other than out of memory
        debug_assert!(err == 0);
    }

    /// HiGHS: ipx/basis.cc::Basis::Load (lines 80-113)
    ///
    /// `const int* basic_status`: C `int`, so `&[i32]`.
    pub fn load(&mut self, basic_status: &[i32]) -> Int {
        let m = self.model_.rows();
        let n = self.model_.cols();

        // Change member variables only when basis is valid.
        let mut basis: Vec<Int> = Vec::new();
        let mut map2basis: Vec<Int> = vec![0; (n + m) as usize];
        let mut p: Int = 0;
        for j in 0..n + m {
            let ju = j as usize;
            match basic_status[ju] {
                NONBASIC_FIXED => map2basis[ju] = -2,
                NONBASIC => map2basis[ju] = -1,
                BASIC => {
                    basis.push(j);
                    map2basis[ju] = p;
                    p += 1;
                }
                BASIC_FREE => {
                    basis.push(j);
                    map2basis[ju] = p + m;
                    p += 1;
                }
                _ => return IPX_ERROR_invalid_basis,
            }
        }
        if p != m {
            return IPX_ERROR_invalid_basis;
        }

        self.basis_.copy_from_slice(&basis);
        self.map2basis_.copy_from_slice(&map2basis);
        self.factorize()
    }

    /// HiGHS: ipx/basis.cc::Basis::Factorize (lines 115-171)
    pub fn factorize(&mut self) -> Int {
        self.factorize_with(true)
    }

    /// rusolve BASICLU time limit (not in the HiGHS reference): ipx/basis.cc::Basis::Factorize(allow_timeout); the flag is
    /// false for a slack basis.
    pub fn factorize_with(&mut self, allow_timeout: bool) -> Int {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let ai = model.AI();
        let timer = Timer::new(0.0);

        // Build column pointers for passing to LU factorization.
        let mut begin = vec![0; m as usize];
        let mut end = vec![0; m as usize];
        let mut basis_num_nz: Int = 0;
        for i in 0..m {
            let iu = i as usize;
            debug_assert!(self.basis_[iu] >= 0);
            begin[iu] = ai.begin(self.basis_[iu]);
            end[iu] = ai.end(self.basis_[iu]);
            basis_num_nz += end[iu] - begin[iu];
        }

        let mut h = format!(
            "    Start  factorization {}: nonzeros in basis = {}{}{}\n",
            format_int(self.num_factorizations_ + 1, 3),
            format_int(basis_num_nz, 9),
            format_str("", 14),
            time_fmt(self.control_.elapsed())
        );
        self.control_.h_interval_log(&mut h);

        let mut err: Int = 0; // return code
        loop {
            let highs_time_limit = self.control_.time_limit();
            debug_assert!(highs_time_limit >= 0.0);
            self.lu_.borrow_mut().time_start(lu_time());
            let mut basiclu_time_limit = f64::INFINITY;
            if allow_timeout && highs_time_limit < f64::INFINITY {
                basiclu_time_limit = highs_time_limit - self.control_.elapsed();
                if basiclu_time_limit <= 0.0 {
                    return IPX_ERROR_time_interrupt;
                }
            }
            self.lu_.borrow_mut().time_limit(basiclu_time_limit);
            let flag =
                self.lu_
                    .borrow_mut()
                    .factorize(&begin, &end, ai.rowidx(), ai.values(), false);
            if flag == IPX_ERROR_time_interrupt {
                return flag;
            }
            self.num_factorizations_ += 1;
            let fill = self.lu_.borrow().fill_factor();
            self.fill_factors_.push(fill);
            if (flag & 2) != 0 {
                self.adapt_to_singular_factorization();
                err = IPX_ERROR_basis_singular;
                break;
            }
            if (flag & 1) != 0 && self.tighten_lu_pivot_tol() {
                // The factorization was numerically unstable and the pivot
                // tolerance could be tightened. Repeat.
                continue;
            }
            if (flag & 1) != 0 && self.control_.debug_enabled(3) {
                let pivottol = self.lu_.borrow().pivottol();
                self.control_.debug(
                    3,
                    &format!(
                        " LU factorization unstable with pivot tolerance {}\n",
                        ostream_double(pivottol)
                    ),
                );
                // ignore instability
            }
            break;
        }
        self.time_factorize_ += timer.elapsed();
        self.factorization_is_fresh_ = true;
        let mut h = format!(
            "    Finish factorization {}: fill factor = {}{}{}\n",
            format_int(self.num_factorizations_, 3),
            fixed(self.lu_.borrow().fill_factor(), 6, 2),
            format_str("", 23),
            time_fmt(self.control_.elapsed())
        );
        self.control_.h_interval_log(&mut h);
        err
    }

    /// HiGHS: ipx/basis.cc::Basis::FactorizationIsFresh (lines 173-175)
    pub fn factorization_is_fresh(&self) -> bool {
        self.factorization_is_fresh_
    }

    /// HiGHS: ipx/basis.cc::Basis::GetLuFactors (lines 177-181)
    #[allow(non_snake_case)]
    pub fn get_lu_factors(
        &self,
        L: &mut SparseMatrix,
        U: &mut SparseMatrix,
        rowperm: &mut [Int],
        colperm: &mut [Int],
    ) {
        debug_assert!(self.factorization_is_fresh());
        self.lu_
            .borrow_mut()
            .get_factors(Some(L), Some(U), Some(rowperm), Some(colperm), None);
    }

    /// HiGHS: ipx/basis.cc::Basis::SolveDense (lines 183-185)
    pub fn solve_dense(&self, rhs: &Vector, lhs: &mut Vector, trans: u8) {
        self.lu_.borrow_mut().solve_dense(rhs, lhs, trans);
    }

    /// HiGHS: ipx/basis.cc::Basis::SolveForUpdate(Int, IndexedVector&) (lines 187-217)
    pub fn solve_for_update_lhs(&mut self, j: Int, lhs: &mut IndexedVector) {
        let p = self.position_of(j);
        let timer = Timer::new(0.0);
        let model = Rc::clone(&self.model_);
        let dim = model.rows();
        debug_assert!(dim > 0);
        if p < 0 {
            // ftran
            let ai = model.AI();
            let begin = ai.begin(j);
            let end = ai.end(j);
            let nz = end - begin;
            let bi = &ai.rowidx()[begin as usize..end as usize];
            let bx = &ai.values()[begin as usize..end as usize];
            self.lu_.borrow_mut().ftran_for_update_lhs(nz, bi, bx, lhs);
            self.num_ftran_ += 1;
            let density = (1.0 * lhs.nnz() as f64) / dim as f64;
            self.sum_ftran_density_ += density;
            if lhs.sparse() {
                self.num_ftran_sparse_ += 1;
            }
            self.time_ftran_ += timer.elapsed();
        } else {
            // btran
            self.lu_.borrow_mut().btran_for_update_lhs(p, lhs);
            self.num_btran_ += 1;
            let density = (1.0 * lhs.nnz() as f64) / dim as f64;
            self.sum_btran_density_ += density;
            if lhs.sparse() {
                self.num_btran_sparse_ += 1;
            }
            self.time_btran_ += timer.elapsed();
        }
    }

    /// rusolve `lp_crossover_batch_push`: FTRAN of a caller-built sparse right-hand side
    /// (the sum of several nonbasic columns); `idx` must hold distinct row indices.
    pub fn solve_rhs_lhs(&mut self, idx: &[Int], val: &[f64], lhs: &mut IndexedVector) {
        let timer = Timer::new(0.0);
        let dim = self.model_.rows();
        self.lu_
            .borrow_mut()
            .ftran_for_update_lhs(idx.len() as Int, idx, val, lhs);
        self.num_ftran_ += 1;
        self.sum_ftran_density_ += (1.0 * lhs.nnz() as f64) / dim as f64;
        if lhs.sparse() {
            self.num_ftran_sparse_ += 1;
        }
        self.time_ftran_ += timer.elapsed();
    }

    /// HiGHS: ipx/basis.cc::Basis::SolveForUpdate(Int) (lines 219-236)
    pub fn solve_for_update(&mut self, j: Int) {
        let p = self.position_of(j);
        let timer = Timer::new(0.0);
        let model = Rc::clone(&self.model_);
        if p < 0 {
            let ai = model.AI();
            let begin = ai.begin(j);
            let end = ai.end(j);
            let nz = end - begin;
            let bi = &ai.rowidx()[begin as usize..end as usize];
            let bx = &ai.values()[begin as usize..end as usize];
            self.lu_.borrow_mut().ftran_for_update(nz, bi, bx);
            self.time_ftran_ += timer.elapsed();
        } else {
            self.lu_.borrow_mut().btran_for_update(p);
            self.time_btran_ += timer.elapsed();
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::TableauRow (lines 238-305)
    ///
    /// C++ default arguments: `ignore_fixed = false`.
    pub fn tableau_row(
        &mut self,
        jb: Int,
        btran: &mut IndexedVector,
        row: &mut IndexedVector,
        ignore_fixed: bool,
    ) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        debug_assert!(self.is_basic(jb));
        self.solve_for_update_lhs(jb, btran);

        // Estimate if tableau row is sparse.
        let mut is_sparse = btran.sparse();
        if is_sparse {
            let ait = model.AIt();
            let bi = btran.pattern();
            let mut nz: Int = 0;
            for k in 0..btran.nnz() {
                let i = bi[k as usize];
                nz += ait.end(i) - ait.begin(i);
            }
            nz /= 2; // guess for overlap
            if (nz as f64) > kHypersparseThreshold * (n as f64) {
                is_sparse = false;
            }
        }

        if is_sparse {
            // sparse vector * sparse matrix: accesses A rowwise
            let ait = model.AIt();
            let ati = ait.rowidx();
            let atx = ait.values();
            row.set_to_zero();
            let mut nz: Int = 0;
            {
                let (row_elems, row_pattern) = row.elements_and_pattern_mut();
                for k in 0..btran.nnz() {
                    let i = btran.pattern()[k as usize];
                    let temp = btran.get(i);
                    let begin = ait.begin(i);
                    let end = ait.end(i);
                    for p in begin..end {
                        let j = ati[p as usize];
                        let ju = j as usize;
                        if self.map2basis_[ju] == -1 || (self.map2basis_[ju] == -2 && !ignore_fixed)
                        {
                            self.map2basis_[ju] -= 2; // mark column
                            row_pattern[nz as usize] = j;
                            nz += 1;
                        }
                        if self.map2basis_[ju] < -2 {
                            // marked column
                            row_elems[ju] += temp * atx[p as usize];
                        }
                    }
                }
                for k in 0..nz {
                    // reset marked
                    self.map2basis_[row_pattern[k as usize] as usize] += 2;
                }
            }
            row.set_nnz(nz);
        } else {
            // dense vector * sparse matrix: accesses A columnwise
            let ai = model.AI();
            let ai_idx = ai.rowidx();
            let ax = ai.values();
            for j in 0..n + m {
                let ju = j as usize;
                let mut result: f64 = 0.0;
                if self.map2basis_[ju] == -1 || (self.map2basis_[ju] == -2 && !ignore_fixed) {
                    let begin = ai.begin(j);
                    let end = ai.end(j);
                    for p in begin..end {
                        result += ax[p as usize] * btran.get(ai_idx[p as usize]);
                    }
                }
                *row.get_mut(j) = result;
            }
            row.invalidate_pattern();
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::ExchangeIfStable (lines 307-343)
    pub fn exchange_if_stable(
        &mut self,
        jb: Int,
        jn: Int,
        tableau_entry: f64,
        sys: i32,
        exchanged: &mut bool,
    ) -> Int {
        debug_assert!(self.is_basic(jb));
        debug_assert!(self.is_nonbasic(jn));
        if sys > 0 {
            // forward system needs to be solved
            self.solve_for_update(jn);
        }
        if sys < 0 {
            // transposed system needs to be solved
            self.solve_for_update(jb);
        }

        // Update factorization.
        *exchanged = false;
        let timer = Timer::new(0.0);
        let err = self.lu_.borrow_mut().update(tableau_entry);
        self.time_update_ += timer.elapsed();
        if err != 0 {
            if self.factorization_is_fresh() && !self.tighten_lu_pivot_tol() {
                return IPX_ERROR_basis_too_ill_conditioned;
            }
            if self.control_.debug_enabled(3) {
                let updates = self.lu_.borrow().updates();
                self.control_.debug(
                    3,
                    &format!(
                        " stability check forced refactorization after {} updates\n",
                        updates - 1
                    ),
                );
            }
            return self.factorize(); // refactorizes the old basis
        }

        // Update basis.
        let ib = self.position_of(jb);
        debug_assert!(self.basis_[ib as usize] == jb);
        self.basis_[ib as usize] = jn;
        self.map2basis_[jn as usize] = ib; // status now BASIC
        self.map2basis_[jb as usize] = -1; // status now NONBASIC
        self.num_updates_ += 1;
        self.factorization_is_fresh_ = false;
        *exchanged = true;

        let need_fresh = self.lu_.borrow_mut().need_fresh_factorization();
        if need_fresh {
            return self.factorize();
        }
        0
    }

    /// HiGHS: ipx/basis.cc::Basis::ComputeBasicSolution (lines 345-369)
    pub fn compute_basic_solution(&self, x: &mut Vector, y: &mut Vector, z: &mut Vector) {
        let model = self.model();
        let m = model.rows();
        let n = model.cols();
        let b = model.b();
        let cost = model.c();
        let ai = model.AI();

        // Compute x[basic] so that Ax=b. Use y as workspace.
        y.assign(b);
        for j in 0..n + m {
            if self.is_nonbasic(j) {
                scatter_column(ai, j, -x[j as usize], &mut y[..]);
            }
        }
        let rhs = y.clone();
        self.solve_dense(&rhs, y, b'N');
        for p in 0..m {
            x[self.basis_[p as usize] as usize] = y[p as usize];
        }

        // Compute y and z[nonbasic] so that AI'y+z=c.
        for p in 0..m {
            let jb = self.basis_[p as usize] as usize;
            y[p as usize] = cost[jb] - z[jb];
        }
        let rhs = y.clone();
        self.solve_dense(&rhs, y, b'T');
        for j in 0..n + m {
            if self.is_nonbasic(j) {
                z[j as usize] = cost[j as usize] - dot_column(ai, j, &y[..]);
            }
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::ConstructBasisFromWeights (lines 371-403)
    // IPX quirk: the trailing `if info.errflag != 0 { return; }` mirrors the C++ source.
    #[allow(clippy::needless_return)]
    pub fn construct_basis_from_weights(&mut self, colscale: &[f64], info: &mut Info) {
        info.errflag = 0;
        info.dependent_rows = 0;
        info.dependent_cols = 0;

        if self.control_.crash_basis() != 0 {
            let mut interrupt = false;
            self.crash_basis(colscale, &mut interrupt);
            if interrupt {
                info.errflag = IPX_ERROR_time_interrupt;
                return;
            }
            let mut sigma = self.min_singular_value();
            self.control_.debug(
                1,
                &format!(
                    "{}{}\n",
                    textline("Minimum singular value of crash basis:"),
                    sci2(sigma)
                ),
            );
            self.repair(info, &mut interrupt);
            if interrupt {
                info.errflag = IPX_ERROR_time_interrupt;
                return;
            }
            if info.basis_repairs < 0 {
                self.control_.h_log_str(" discarding crash basis\n");
                self.set_to_slack_basis();
            } else if info.basis_repairs > 0 {
                sigma = self.min_singular_value();
                self.control_.debug(
                    1,
                    &format!(
                        "{}{}\n",
                        textline("Minimum singular value of repaired crash basis:"),
                        sci2(sigma)
                    ),
                );
            }
        } else {
            self.set_to_slack_basis();
        }
        self.pivot_free_variables_into_basis(colscale, info);
        if info.errflag != 0 {
            return;
        }
        self.pivot_fixed_variables_out_of_basis(colscale, info);
        if info.errflag != 0 {
            return;
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::MinSingularValue (lines 405-415)
    pub fn min_singular_value(&self) -> f64 {
        let m = self.model_.rows();
        let mut v = Vector::new(m);

        // Computes maximum eigenvalue of inverse(B*B').
        let lambda = power_method(
            |x: &Vector, fx: &mut Vector| {
                self.solve_dense(x, fx, b'N');
                let tmp = fx.clone();
                self.solve_dense(&tmp, fx, b'T');
            },
            &mut v,
        );
        (1.0 / lambda).sqrt()
    }

    /// HiGHS: ipx/basis.cc::Basis::SymbolicInvert (lines 417-419)
    ///
    /// `colcounts` is `nullptr` from `DensityInverse` (line 424).
    pub fn symbolic_invert(&self, rowcounts: Option<&mut [Int]>, colcounts: Option<&mut [Int]>) {
        symbolic_invert(&self.model_, &self.basis_, rowcounts, colcounts);
    }

    /// HiGHS: ipx/basis.cc::Basis::DensityInverse (lines 421-430)
    pub fn density_inverse(&self) -> f64 {
        let m = self.model_.rows();
        let mut rowcounts: Vec<Int> = vec![0; m as usize];
        self.symbolic_invert(Some(&mut rowcounts), None);
        // Accumulating rowcounts would result in overflow for large LPs.
        let mut density = 0.0;
        for &rc in &rowcounts {
            density += 1.0 * rc as f64 / m as f64;
        }
        density / m as f64
    }

    /// HiGHS: ipx/basis.cc::Basis::model (lines 432-434)
    pub fn model(&self) -> &Model {
        &self.model_
    }

    /// The shared handle (for objects that keep `const Model&`).
    pub fn model_rc(&self) -> Rc<Model> {
        Rc::clone(&self.model_)
    }

    /// HiGHS: ipx/basis.cc::Basis::factorizations (lines 436-438)
    pub fn factorizations(&self) -> Int {
        self.num_factorizations_
    }

    /// HiGHS: ipx/basis.cc::Basis::updates_total (lines 440-442)
    pub fn updates_total(&self) -> Int {
        self.num_updates_
    }

    /// HiGHS: ipx/basis.cc::Basis::frac_ftran_sparse (lines 444-446)
    pub fn frac_ftran_sparse(&self) -> f64 {
        1.0 * self.num_ftran_sparse_ as f64 / self.num_ftran_ as f64
    }

    /// HiGHS: ipx/basis.cc::Basis::frac_btran_sparse (lines 448-450)
    pub fn frac_btran_sparse(&self) -> f64 {
        1.0 * self.num_btran_sparse_ as f64 / self.num_btran_ as f64
    }

    /// HiGHS: ipx/basis.cc::Basis::time_factorize (lines 452-454)
    pub fn time_factorize(&self) -> f64 {
        self.time_factorize_
    }

    /// HiGHS: ipx/basis.cc::Basis::time_ftran (lines 456-458)
    pub fn time_ftran(&self) -> f64 {
        self.time_ftran_
    }

    /// HiGHS: ipx/basis.cc::Basis::time_btran (lines 460-462)
    pub fn time_btran(&self) -> f64 {
        self.time_btran_
    }

    /// HiGHS: ipx/basis.cc::Basis::time_update (lines 464-466)
    pub fn time_update(&self) -> f64 {
        self.time_update_
    }

    /// HiGHS: ipx/basis.cc::Basis::mean_fill (lines 468-476)
    pub fn mean_fill(&self) -> f64 {
        if self.fill_factors_.is_empty() {
            return 0.0;
        }
        let mut mean = 1.0;
        let num_factors = self.fill_factors_.len() as Int;
        for &f in &self.fill_factors_ {
            mean *= f.powf(1.0 / num_factors as f64);
        }
        mean
    }

    /// HiGHS: ipx/basis.cc::Basis::max_fill (lines 478-482)
    ///
    /// `*std::max_element` (first maximal element; empty is UB in C++).
    pub fn max_fill(&self) -> f64 {
        if self.fill_factors_.is_empty() {
            return 0.0;
        }
        let v = &self.fill_factors_;
        let mut largest = v[0];
        for &f in &v[1..] {
            if largest < f {
                largest = f;
            }
        }
        largest
    }
}
