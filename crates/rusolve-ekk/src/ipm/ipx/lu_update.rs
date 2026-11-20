//! HiGHS source: ipm/ipx/lu_update.h, lu_update.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Abstract base class -> trait. The base member `updates_` lives in each
//! implementor and is reached through `updates_ref`/`updates_mut`. The
//! public non-virtual wrappers are provided methods; each implementor
//! supplies the private virtuals (`_factorize`, ...). Overloads:
//! `FtranForUpdate(nz, bi, bx)` -> `ftran_for_update`,
//! `FtranForUpdate(nz, bi, bx, lhs)` -> `ftran_for_update_lhs`; same for
//! `BtranForUpdate`; the setter `pivottol(double)` -> `set_pivottol`.

use super::indexed_vector::IndexedVector;
use super::sparse_matrix::SparseMatrix;
use super::*;

/// HiGHS: ipx/lu_update.h::LuUpdate (abstract)
pub trait LuUpdate {
    /// base member `updates_` (read)
    fn updates_ref(&self) -> Int;
    /// base member `updates_` (write)
    fn updates_mut(&mut self) -> &mut Int;

    /// HiGHS: ipx/lu_update.cc::LuUpdate::Factorize (lines 5-9)
    fn factorize(
        &mut self,
        Bbegin: &[Int],
        Bend: &[Int],
        Bi: &[Int],
        Bx: &[f64],
        strict_abs_pivottol: bool,
    ) -> Int {
        *self.updates_mut() = 0;
        self._factorize(Bbegin, Bend, Bi, Bx, strict_abs_pivottol)
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::GetFactors (lines 11-14)
    ///
    /// Every argument may be `nullptr`.
    fn get_factors(
        &mut self,
        L: Option<&mut SparseMatrix>,
        U: Option<&mut SparseMatrix>,
        rowperm: Option<&mut [Int]>,
        colperm: Option<&mut [Int]>,
        dependent_cols: Option<&mut Vec<Int>>,
    ) {
        self._get_factors(L, U, rowperm, colperm, dependent_cols);
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::SolveDense (lines 16-18)
    fn solve_dense(&mut self, rhs: &Vector, lhs: &mut Vector, trans: u8) {
        self._solve_dense(rhs, lhs, trans);
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::FtranForUpdate(nz, bi, bx) (lines 20-22)
    fn ftran_for_update(&mut self, nz: Int, bi: &[Int], bx: &[f64]) {
        self._ftran_for_update(nz, bi, bx);
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::FtranForUpdate(nz, bi, bx, lhs) (lines 24-27)
    fn ftran_for_update_lhs(&mut self, nz: Int, bi: &[Int], bx: &[f64], lhs: &mut IndexedVector) {
        self._ftran_for_update_lhs(nz, bi, bx, lhs);
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::BtranForUpdate(j) (lines 29-31)
    fn btran_for_update(&mut self, j: Int) {
        self._btran_for_update(j);
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::BtranForUpdate(j, lhs) (lines 33-35)
    fn btran_for_update_lhs(&mut self, j: Int, lhs: &mut IndexedVector) {
        self._btran_for_update_lhs(j, lhs);
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::Update (lines 37-40)
    fn update(&mut self, pivot: f64) -> Int {
        *self.updates_mut() += 1;
        self._update(pivot)
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::NeedFreshFactorization (lines 42-44)
    fn need_fresh_factorization(&mut self) -> bool {
        self._need_fresh_factorization()
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::fill_factor (lines 46-48)
    fn fill_factor(&self) -> f64 {
        self._fill_factor()
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::pivottol() (lines 50-52)
    fn pivottol(&self) -> f64 {
        self._pivottol()
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::pivottol(double) (lines 54-56)
    fn set_pivottol(&mut self, new_pivottol: f64) {
        self._set_pivottol(new_pivottol);
    }
    /// rusolve BASICLU time limit (not in the HiGHS reference): ipx/lu_update.cc::LuUpdate::timeStart
    fn time_start(&mut self, new_time_start: f64) {
        self._time_start(new_time_start);
    }
    /// rusolve BASICLU time limit (not in the HiGHS reference): ipx/lu_update.cc::LuUpdate::timeLimit
    fn time_limit(&mut self, new_time_limit: f64) {
        self._time_limit(new_time_limit);
    }
    /// HiGHS: ipx/lu_update.cc::LuUpdate::updates (lines 58-60)
    fn updates(&self) -> Int {
        self.updates_ref()
    }

    /// Test hook (not in HiGHS): the implementor's internals that the
    /// golden replay compares against the oracle, which reads them with
    /// `#define private public` (tools/highs-oracle/src/ipx.cpp basisObj).
    #[doc(hidden)]
    fn oracle_stats(&self) -> LuOracleStats {
        LuOracleStats::Other
    }

    // private virtuals (lu_update.h)
    fn _factorize(
        &mut self,
        Bbegin: &[Int],
        Bend: &[Int],
        Bi: &[Int],
        Bx: &[f64],
        strict_abs_pivottol: bool,
    ) -> Int;
    fn _get_factors(
        &mut self,
        L: Option<&mut SparseMatrix>,
        U: Option<&mut SparseMatrix>,
        rowperm: Option<&mut [Int]>,
        colperm: Option<&mut [Int]>,
        dependent_cols: Option<&mut Vec<Int>>,
    );
    fn _solve_dense(&mut self, rhs: &Vector, lhs: &mut Vector, trans: u8);
    fn _ftran_for_update(&mut self, nz: Int, bi: &[Int], bx: &[f64]);
    fn _ftran_for_update_lhs(&mut self, nz: Int, bi: &[Int], bx: &[f64], lhs: &mut IndexedVector);
    fn _btran_for_update(&mut self, p: Int);
    fn _btran_for_update_lhs(&mut self, p: Int, lhs: &mut IndexedVector);
    fn _update(&mut self, pivot: f64) -> Int;
    fn _need_fresh_factorization(&mut self) -> bool;
    fn _fill_factor(&self) -> f64;
    fn _pivottol(&self) -> f64;
    fn _set_pivottol(&mut self, new_pivottol: f64);
    fn _time_start(&mut self, new_time_start: f64);
    fn _time_limit(&mut self, new_time_limit: f64);
}

/// See [`LuUpdate::oracle_stats`].
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq)]
pub enum LuOracleStats {
    /// `BasicLu`: BASICLU memory sizes after the Reallocate loop and the
    /// vector sizes.
    BasicLu {
        lmem: Int,
        umem: Int,
        wmem: Int,
        istore_size: usize,
        xstore_size: usize,
        li_size: usize,
        ui_size: usize,
        wi_size: usize,
        fill_factor: f64,
    },
    /// `ForrestTomlin`: number of updates in the eta file.
    ForrestTomlin {
        ft_updates: usize,
        fill_factor: f64,
    },
    Other,
}
