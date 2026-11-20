//! HiGHS source: ipm/ipx/basis.h, basis.cc (class Basis)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `basis.cc` (965 lines) is split: [`part1`] lines 17-482 (public API),
//! [`part2`] lines 484-965 (repair, crash, pivoting free/fixed variables,
//! `CopyBasic`, `reportBasisData`). Header inline members are implemented
//! here.
//!
//! C++ `mutable std::unique_ptr<LuUpdate> lu_` (the `const` solve methods
//! update LU workspace) is `RefCell<Box<dyn LuUpdate>>`.

use super::control::Control;
use super::lu_update::LuUpdate;
use super::model::Model;
use super::*;
use std::cell::RefCell;
use std::rc::Rc;

mod part1;
mod part2;

pub use part2::copy_basic;

/// HiGHS: ipx/basis.h::Basis::BasicStatus (plain enum, used as integers)
pub type BasicStatus = Int;
pub const NONBASIC_FIXED: BasicStatus = -2;
pub const NONBASIC: BasicStatus = -1;
pub const BASIC: BasicStatus = 0;
pub const BASIC_FREE: BasicStatus = 1;

/// HiGHS: ipx/basis.h::Basis
pub struct Basis {
    control_: Rc<Control>,
    model_: Rc<Model>,
    basis_: Vec<Int>,
    map2basis_: Vec<Int>,
    lu_: RefCell<Box<dyn LuUpdate>>,
    factorization_is_fresh_: bool,
    num_factorizations_: Int,
    num_updates_: Int,
    num_ftran_: Int,
    num_btran_: Int,
    num_ftran_sparse_: Int,
    num_btran_sparse_: Int,
    time_ftran_: f64,
    time_btran_: f64,
    time_update_: f64,
    time_factorize_: f64,
    fill_factors_: Vec<f64>,
    sum_ftran_density_: f64,
    sum_btran_density_: f64,
}

impl Basis {
    pub const kBasisRepairThreshold: f64 = 1e5;
    pub const kMaxBasisRepair: Int = 200;

    /// HiGHS: ipx/basis.h::Basis::operator[] (inline)
    pub fn at(&self, p: Int) -> Int {
        self.basis_[p as usize]
    }

    /// HiGHS: ipx/basis.h::Basis::StatusOf (inline)
    pub fn status_of(&self, j: Int) -> BasicStatus {
        let m = self.model_.rows();
        let p = self.map2basis_[j as usize];
        debug_assert!(p >= -2 && p < 2 * m);
        if p < 0 {
            if p == -1 { NONBASIC } else { NONBASIC_FIXED }
        } else if p < m {
            BASIC
        } else {
            BASIC_FREE
        }
    }

    /// HiGHS: ipx/basis.h::Basis::PositionOf (inline)
    pub fn position_of(&self, j: Int) -> Int {
        let m = self.model_.rows();
        let p = self.map2basis_[j as usize];
        debug_assert!(p >= -2 && p < 2 * m);
        if p < 0 {
            -1
        } else if p < m {
            p
        } else {
            p - m
        }
    }

    /// HiGHS: ipx/basis.h::Basis::IsBasic (inline)
    pub fn is_basic(&self, j: Int) -> bool {
        self.status_of(j) == BASIC || self.status_of(j) == BASIC_FREE
    }

    /// HiGHS: ipx/basis.h::Basis::IsNonbasic (inline)
    pub fn is_nonbasic(&self, j: Int) -> bool {
        self.status_of(j) == NONBASIC || self.status_of(j) == NONBASIC_FIXED
    }

    /// Test hook (not in HiGHS): the private state the golden replay
    /// compares against the oracle, which reads it with `#define private
    /// public` (tools/highs-oracle/src/ipx.cpp basisObj).
    #[doc(hidden)]
    pub fn oracle_view(&self) -> BasisOracleView<'_> {
        let lu = self.lu_.borrow();
        BasisOracleView {
            basis: &self.basis_,
            map2basis: &self.map2basis_,
            fresh: self.factorization_is_fresh_,
            factorizations: self.num_factorizations_,
            updates: self.num_updates_,
            ftran: self.num_ftran_,
            btran: self.num_btran_,
            ftran_sparse: self.num_ftran_sparse_,
            btran_sparse: self.num_btran_sparse_,
            fill_factors: &self.fill_factors_,
            sum_ftran_density: self.sum_ftran_density_,
            sum_btran_density: self.sum_btran_density_,
            pivottol: lu.pivottol(),
            lu: lu.oracle_stats(),
        }
    }
}

/// See [`Basis::oracle_view`].
#[doc(hidden)]
pub struct BasisOracleView<'a> {
    pub basis: &'a [Int],
    pub map2basis: &'a [Int],
    pub fresh: bool,
    pub factorizations: Int,
    pub updates: Int,
    pub ftran: Int,
    pub btran: Int,
    pub ftran_sparse: Int,
    pub btran_sparse: Int,
    pub fill_factors: &'a [f64],
    pub sum_ftran_density: f64,
    pub sum_btran_density: f64,
    pub pivottol: f64,
    pub lu: super::lu_update::LuOracleStats,
}
