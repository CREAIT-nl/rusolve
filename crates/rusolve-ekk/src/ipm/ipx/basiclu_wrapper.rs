//! HiGHS source: ipm/ipx/basiclu_wrapper.h, basiclu_wrapper.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `BasicLu` is the default `LuUpdate` of `Basis` (`lu_kernel <= 0`):
//! BASICLU with its own Forrest-Tomlin update. It owns the BASICLU store
//! (`istore_`, `xstore_`, `Li_` ..) and calls `crate::ipm::basiclu`.

use super::control::{Control, sci2, textline};
use super::indexed_vector::IndexedVector;
use super::lu_update::LuUpdate;
use super::sparse_matrix::{SparseMatrix, remove_diagonal};
use super::*;
use crate::ipm::basiclu::basiclu_factorize::basiclu_factorize;
use crate::ipm::basiclu::basiclu_get_factors::basiclu_get_factors;
use crate::ipm::basiclu::basiclu_initialize::basiclu_initialize;
use crate::ipm::basiclu::basiclu_solve_dense::basiclu_solve_dense;
use crate::ipm::basiclu::basiclu_solve_for_update::basiclu_solve_for_update;
use crate::ipm::basiclu::basiclu_update::basiclu_update;
use crate::ipm::basiclu::*;
use std::rc::Rc;

/// HiGHS: ipx/basiclu_wrapper.h::BasicLu
pub struct BasicLu {
    /// base member `LuUpdate::updates_`
    updates_: Int,
    control_: Rc<Control>,
    istore_: Vec<Int>,
    xstore_: Vec<f64>,
    Li_: Vec<Int>,
    Ui_: Vec<Int>,
    Wi_: Vec<Int>,
    Lx_: Vec<f64>,
    Ux_: Vec<f64>,
    Wx_: Vec<f64>,
    fill_factor_: f64,
}

impl BasicLu {
    pub const kReallocFactor: f64 = 1.5;

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::BasicLu (lines 10-32)
    pub fn new(control: Rc<Control>, dim: Int) -> Self {
        let size = (BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * dim) as usize;
        let mut istore_: Vec<Int> = vec![0; size];
        let mut xstore_: Vec<f64> = vec![0.0; size];

        let status = basiclu_initialize(dim, &mut istore_, &mut xstore_);
        if status != BASICLU_OK {
            panic!("basiclu_initialize failed");
        }

        // Set initial size of BASICLU work arrays to 1 element, so that data()
        // does not return NULL (not sure if it would if size is 0).
        let mut this = BasicLu {
            updates_: 0,
            control_: control,
            istore_,
            xstore_,
            Li_: vec![0; 1],
            Ui_: vec![0; 1],
            Wi_: vec![0; 1],
            Lx_: vec![0.0; 1],
            Ux_: vec![0.0; 1],
            Wx_: vec![0.0; 1],
            fill_factor_: 0.0,
        };
        this.xstore_[BASICLU_MEMORYL] = 1.0;
        this.xstore_[BASICLU_MEMORYU] = 1.0;
        this.xstore_[BASICLU_MEMORYW] = 1.0;
        this
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::Reallocate (lines 268-297)
    fn reallocate(&mut self) {
        debug_assert_eq!(self.Li_.len(), self.xstore_[BASICLU_MEMORYL] as usize);
        debug_assert_eq!(self.Lx_.len(), self.xstore_[BASICLU_MEMORYL] as usize);
        debug_assert_eq!(self.Ui_.len(), self.xstore_[BASICLU_MEMORYU] as usize);
        debug_assert_eq!(self.Ux_.len(), self.xstore_[BASICLU_MEMORYU] as usize);
        debug_assert_eq!(self.Wi_.len(), self.xstore_[BASICLU_MEMORYW] as usize);
        debug_assert_eq!(self.Wx_.len(), self.xstore_[BASICLU_MEMORYW] as usize);

        // C: `Int new_size = double sum;` then `new_size *= kReallocFactor;`
        // (double product converted back to Int, truncating).
        if self.xstore_[BASICLU_ADD_MEMORYL] > 0.0 {
            let mut new_size =
                (self.xstore_[BASICLU_MEMORYL] + self.xstore_[BASICLU_ADD_MEMORYL]) as Int;
            new_size = (new_size as f64 * Self::kReallocFactor) as Int;
            self.Li_.resize(new_size as usize, 0);
            self.Lx_.resize(new_size as usize, 0.0);
            self.xstore_[BASICLU_MEMORYL] = new_size as f64;
        }
        if self.xstore_[BASICLU_ADD_MEMORYU] > 0.0 {
            let mut new_size =
                (self.xstore_[BASICLU_MEMORYU] + self.xstore_[BASICLU_ADD_MEMORYU]) as Int;
            new_size = (new_size as f64 * Self::kReallocFactor) as Int;
            self.Ui_.resize(new_size as usize, 0);
            self.Ux_.resize(new_size as usize, 0.0);
            self.xstore_[BASICLU_MEMORYU] = new_size as f64;
        }
        if self.xstore_[BASICLU_ADD_MEMORYW] > 0.0 {
            let mut new_size =
                (self.xstore_[BASICLU_MEMORYW] + self.xstore_[BASICLU_ADD_MEMORYW]) as Int;
            new_size = (new_size as f64 * Self::kReallocFactor) as Int;
            self.Wi_.resize(new_size as usize, 0);
            self.Wx_.resize(new_size as usize, 0.0);
            self.xstore_[BASICLU_MEMORYW] = new_size as f64;
        }
    }
}

impl LuUpdate for BasicLu {
    fn oracle_stats(&self) -> super::lu_update::LuOracleStats {
        super::lu_update::LuOracleStats::BasicLu {
            lmem: self.xstore_[BASICLU_MEMORYL] as Int,
            umem: self.xstore_[BASICLU_MEMORYU] as Int,
            wmem: self.xstore_[BASICLU_MEMORYW] as Int,
            istore_size: self.istore_.len(),
            xstore_size: self.xstore_.len(),
            li_size: self.Li_.len(),
            ui_size: self.Ui_.len(),
            wi_size: self.Wi_.len(),
            fill_factor: self.fill_factor_,
        }
    }

    /// HiGHS: ipx/lu_update.h::LuUpdate::updates_ (base member)
    fn updates_ref(&self) -> Int {
        self.updates_
    }
    fn updates_mut(&mut self) -> &mut Int {
        &mut self.updates_
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_Factorize (lines 34-77)
    fn _factorize(
        &mut self,
        Bbegin: &[Int],
        Bend: &[Int],
        Bi: &[Int],
        Bx: &[f64],
        strict_abs_pivottol: bool,
    ) -> Int {
        if strict_abs_pivottol {
            self.xstore_[BASICLU_REMOVE_COLUMNS] = 1.0;
            self.xstore_[BASICLU_ABS_PIVOT_TOLERANCE] = kLuDependencyTol;
        } else {
            self.xstore_[BASICLU_REMOVE_COLUMNS] = 0.0;
            self.xstore_[BASICLU_ABS_PIVOT_TOLERANCE] = 1e-14; // BASICLU default
        }
        let mut ncall: Int = 0;
        let status = loop {
            let s = basiclu_factorize(
                &mut self.istore_,
                &mut self.xstore_,
                &mut self.Li_,
                &mut self.Lx_,
                &mut self.Ui_,
                &mut self.Ux_,
                &mut self.Wi_,
                &mut self.Wx_,
                Bbegin,
                Bend,
                Bi,
                Bx,
                ncall,
            );
            if s != BASICLU_REALLOCATE {
                break s;
            }
            self.reallocate();
            ncall += 1;
        };
        // rusolve BASICLU time limit (not in the HiGHS reference): a BASICLU time-out is reported to the caller
        if status == BASICLU_WARNING_timeout {
            return IPX_ERROR_time_interrupt;
        }
        if status != BASICLU_OK && status != BASICLU_WARNING_singular_matrix {
            panic!("basiclu_factorize failed");
        }

        let matrix_nz = self.xstore_[BASICLU_MATRIX_NZ] as Int;
        let lnz = self.xstore_[BASICLU_LNZ] as Int;
        let unz = self.xstore_[BASICLU_UNZ] as Int;
        let dim = self.xstore_[BASICLU_DIM] as Int;
        self.fill_factor_ = 1.0 * ((lnz + unz + dim) as f64) / (matrix_nz as f64);

        let normLinv = self.xstore_[BASICLU_NORMEST_LINV];
        let normUinv = self.xstore_[BASICLU_NORMEST_UINV];
        let stability = self.xstore_[BASICLU_RESIDUAL_TEST];
        self.control_.debug(
            3,
            &format!(
                " normLinv = {}, normUinv = {}, stability = {}\n",
                sci2(normLinv),
                sci2(normUinv),
                sci2(stability)
            ),
        );

        let mut ret: Int = 0;
        if stability > kLuStabilityThreshold {
            ret |= 1;
        }
        if status == BASICLU_WARNING_singular_matrix {
            ret |= 2;
        }
        ret
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_GetFactors (lines 79-122)
    fn _get_factors(
        &mut self,
        L: Option<&mut SparseMatrix>,
        U: Option<&mut SparseMatrix>,
        rowperm: Option<&mut [Int]>,
        colperm: Option<&mut [Int]>,
        dependent_cols: Option<&mut Vec<Int>>,
    ) {
        let dim = self.xstore_[BASICLU_DIM] as Int;
        let mut L = L;
        let mut U = U;

        if let Some(l) = L.as_deref_mut() {
            let lnz = self.xstore_[BASICLU_LNZ] as Int;
            l.resize(dim, dim, dim + lnz);
        }
        if let Some(u) = U.as_deref_mut() {
            let unz = self.xstore_[BASICLU_UNZ] as Int;
            u.resize(dim, dim, dim + unz);
        }
        let (lcolptr, lrowidx, lvalue) = match L.as_deref_mut() {
            Some(l) => {
                let (a, b, c) = l.arrays_mut();
                (Some(a), Some(b), Some(c))
            }
            None => (None, None, None),
        };
        let (ucolptr, urowidx, uvalue) = match U.as_mut() {
            Some(u) => {
                let (a, b, c) = u.arrays_mut();
                (Some(a), Some(b), Some(c))
            }
            None => (None, None, None),
        };

        let status = basiclu_get_factors(
            &mut self.istore_,
            &mut self.xstore_,
            &mut self.Li_,
            &mut self.Lx_,
            &mut self.Ui_,
            &mut self.Ux_,
            &mut self.Wi_,
            &mut self.Wx_,
            rowperm,
            colperm,
            lcolptr,
            lrowidx,
            lvalue,
            ucolptr,
            urowidx,
            uvalue,
        );
        if status != BASICLU_OK {
            panic!("basiclu_get_factors failed");
        }

        if let Some(l) = L.as_mut() {
            // Remove unit diagonal from L.
            let num_dropped = remove_diagonal(l, None);
            debug_assert_eq!(num_dropped, dim);
        }
        if let Some(dc) = dependent_cols {
            // Dependent columns are at the end of the BASICLU pivot sequence.
            let rank = self.xstore_[BASICLU_RANK] as Int;
            dc.clear();
            for k in rank..dim {
                dc.push(k);
            }
        }
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_SolveDense (lines 124-132)
    fn _solve_dense(&mut self, rhs: &Vector, lhs: &mut Vector, trans: u8) {
        let status = basiclu_solve_dense(
            &mut self.istore_,
            &mut self.xstore_,
            &mut self.Li_,
            &mut self.Lx_,
            &mut self.Ui_,
            &mut self.Ux_,
            &mut self.Wi_,
            &mut self.Wx_,
            &rhs[..],
            &mut lhs[..],
            trans,
        );
        if status != BASICLU_OK {
            panic!("basiclu_solve_dense failed");
        }
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_FtranForUpdate(nz, bi, bx) (lines 134-150)
    fn _ftran_for_update(&mut self, nzrhs: Int, bi: &[Int], bx: &[f64]) {
        let status = loop {
            let s = basiclu_solve_for_update(
                &mut self.istore_,
                &mut self.xstore_,
                &mut self.Li_,
                &mut self.Lx_,
                &mut self.Ui_,
                &mut self.Ux_,
                &mut self.Wi_,
                &mut self.Wx_,
                nzrhs,
                bi,
                Some(bx),
                None,
                None,
                None,
                b'N',
            );
            if s != BASICLU_REALLOCATE {
                break s;
            }
            self.reallocate();
        };
        if status != BASICLU_OK {
            panic!("basiclu_solve_for_update (ftran without lhs) failed");
        }
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_FtranForUpdate(nz, bi, bx, lhs) (lines 152-173)
    fn _ftran_for_update_lhs(
        &mut self,
        nzrhs: Int,
        bi: &[Int],
        bx: &[f64],
        lhs: &mut IndexedVector,
    ) {
        let mut nzlhs: Int = 0;
        lhs.set_to_zero();
        let status = loop {
            let (elements, pattern) = lhs.elements_and_pattern_mut();
            let s = basiclu_solve_for_update(
                &mut self.istore_,
                &mut self.xstore_,
                &mut self.Li_,
                &mut self.Lx_,
                &mut self.Ui_,
                &mut self.Ux_,
                &mut self.Wi_,
                &mut self.Wx_,
                nzrhs,
                bi,
                Some(bx),
                Some(&mut nzlhs),
                Some(pattern),
                Some(elements),
                b'N',
            );
            if s != BASICLU_REALLOCATE {
                break s;
            }
            self.reallocate();
        };
        if status != BASICLU_OK {
            panic!("basiclu_solve_for_update (ftran with lhs) failed");
        }
        lhs.set_nnz(nzlhs);
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_BtranForUpdate(j) (lines 175-191)
    fn _btran_for_update(&mut self, j: Int) {
        let status = loop {
            let s = basiclu_solve_for_update(
                &mut self.istore_,
                &mut self.xstore_,
                &mut self.Li_,
                &mut self.Lx_,
                &mut self.Ui_,
                &mut self.Ux_,
                &mut self.Wi_,
                &mut self.Wx_,
                0,
                &[j],
                None,
                None,
                None,
                None,
                b'T',
            );
            if s != BASICLU_REALLOCATE {
                break s;
            }
            self.reallocate();
        };
        if status != BASICLU_OK {
            panic!("basiclu_solve_for_update (btran without lhs) failed");
        }
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_BtranForUpdate(j, lhs) (lines 193-213)
    fn _btran_for_update_lhs(&mut self, j: Int, lhs: &mut IndexedVector) {
        let mut nzlhs: Int = 0;
        lhs.set_to_zero();
        let status = loop {
            let (elements, pattern) = lhs.elements_and_pattern_mut();
            let s = basiclu_solve_for_update(
                &mut self.istore_,
                &mut self.xstore_,
                &mut self.Li_,
                &mut self.Lx_,
                &mut self.Ui_,
                &mut self.Ux_,
                &mut self.Wi_,
                &mut self.Wx_,
                0,
                &[j],
                None,
                Some(&mut nzlhs),
                Some(pattern),
                Some(elements),
                b'T',
            );
            if s != BASICLU_REALLOCATE {
                break s;
            }
            self.reallocate();
        };
        if status != BASICLU_OK {
            panic!("basiclu_solve_for_update (btran with lhs) failed");
        }
        lhs.set_nnz(nzlhs);
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_Update (lines 215-246)
    fn _update(&mut self, pivot: f64) -> Int {
        let max_eta_old = self.xstore_[BASICLU_MAX_ETA];
        let status = loop {
            let s = basiclu_update(
                &mut self.istore_,
                &mut self.xstore_,
                &mut self.Li_,
                &mut self.Lx_,
                &mut self.Ui_,
                &mut self.Ux_,
                &mut self.Wi_,
                &mut self.Wx_,
                pivot,
            );
            if s != BASICLU_REALLOCATE {
                break s;
            }
            self.reallocate();
        };
        if status != BASICLU_OK && status != BASICLU_ERROR_singular_update {
            panic!("basiclu_update failed");
        }
        if status == BASICLU_ERROR_singular_update {
            return -1;
        }

        // Print a debugging message if a new eta entry is large.
        let max_eta = self.xstore_[BASICLU_MAX_ETA];
        if max_eta > 1e10 && max_eta > max_eta_old {
            self.control_
                .debug(3, &format!(" max eta = {}\n", sci2(max_eta)));
        }

        // stability check
        let pivot_error = self.xstore_[BASICLU_PIVOT_ERROR];
        if pivot_error > kFtDiagErrorTol {
            self.control_.debug(
                3,
                &format!(
                    " relative error in new diagonal entry of U = {}\n",
                    sci2(pivot_error)
                ),
            );
            return 1;
        }
        0
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_NeedFreshFactorization (lines 248-254)
    fn _need_fresh_factorization(&mut self) -> bool {
        let dim = self.xstore_[BASICLU_DIM] as Int;
        let nforrest = self.xstore_[BASICLU_NFORREST] as Int;
        let update_cost = self.xstore_[BASICLU_UPDATE_COST];

        nforrest == dim || update_cost > 1.0
    }

    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_fill_factor (lines 256-258)
    fn _fill_factor(&self) -> f64 {
        self.fill_factor_
    }
    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_pivottol() (lines 260-262)
    fn _pivottol(&self) -> f64 {
        self.xstore_[BASICLU_REL_PIVOT_TOLERANCE]
    }
    /// HiGHS: ipx/basiclu_wrapper.cc::BasicLu::_pivottol(double) (lines 264-266)
    fn _set_pivottol(&mut self, new_pivottol: f64) {
        self.xstore_[BASICLU_REL_PIVOT_TOLERANCE] = new_pivottol;
    }
    /// rusolve BASICLU time limit (not in the HiGHS reference): ipx/basiclu_wrapper.cc::BasicLu::_timeStart
    fn _time_start(&mut self, new_time_start: f64) {
        self.xstore_[BASICLU_TIME_START] = new_time_start;
    }
    /// rusolve BASICLU time limit (not in the HiGHS reference): ipx/basiclu_wrapper.cc::BasicLu::_timeLimit
    fn _time_limit(&mut self, new_time_limit: f64) {
        self.xstore_[BASICLU_TIME_LIMIT] = new_time_limit;
    }
}

#[cfg(test)]
mod time_limit_tests {
    use super::*;
    use crate::ipm::basiclu::lu_internal::lu_time;

    /// rusolve BASICLU time limit (not in the HiGHS reference): a factorization past its limit reports
    /// `IPX_ERROR_time_interrupt` instead of running to the end
    /// (bdry2 starting basis).
    #[test]
    fn factorize_times_out_and_unlimited_succeeds() {
        let mut lu = BasicLu::new(Rc::new(Control::new()), 3);
        let begin = [0, 1, 2];
        let end = [1, 2, 3];
        let bi = [0, 1, 2];
        let bx = [1.0, 1.0, 1.0];
        lu.time_start(lu_time());
        lu.time_limit(-1.0);
        assert_eq!(
            lu.factorize(&begin, &end, &bi, &bx, false),
            IPX_ERROR_time_interrupt
        );
        lu.time_start(lu_time());
        lu.time_limit(f64::INFINITY);
        assert_eq!(lu.factorize(&begin, &end, &bi, &bx, false), 0);
    }
}
