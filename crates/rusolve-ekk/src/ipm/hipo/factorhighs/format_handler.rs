//! HiGHS source: ipm/hipo/factorhighs/FormatHandler.h, ipm/hipo/factorhighs/FormatHandler.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Base of the frontal-matrix format handlers (virtual-dispatch mapping in
//! `hipo/mod.rs`). `clique_ptr_` is an offset into the serial `CliqueStack`
//! buffer (`None` = nullptr: the clique lives in the local `clique_`); the
//! methods that touch the clique take that buffer as `stack`.

use super::super::{Int, Int64};
use super::Regul;
use super::symbolic::Symbolic;

/// HiGHS: FormatHandler.h `const Int extra_space_frontal`
pub const extra_space_frontal: Int = 10;

/// HiGHS: FormatHandler.h `class FormatHandler` (data members)
pub struct FormatHandler<'a> {
    pub S_: &'a Symbolic,
    /// C++ `const Regul&`; `Regul` is two doubles, copied.
    pub regul_: Regul,
    pub sn_: Int,
    pub nb_: Int,
    pub sn_size_: Int,
    pub ldf_: Int,
    pub ldc_: Int,
    pub frontal_: &'a mut Vec<f64>,
    pub clique_: Vec<f64>,
    pub clique_ptr_: Option<Int64>,
    pub local_reg_: Vec<f64>,
    pub swaps_: Vec<Int>,
    pub pivot_2x2_: Vec<f64>,
}

impl<'a> FormatHandler<'a> {
    /// HiGHS: FormatHandler.cpp::FormatHandler::FormatHandler (lines 9-24)
    pub fn new(
        S: &'a Symbolic,
        sn: Int,
        regul: &Regul,
        frontal: &'a mut Vec<f64>,
        clique_ptr: Option<Int64>,
    ) -> Self {
        let sn_size_ = S.snStart(sn + 1) - S.snStart(sn);
        let ldf_ = (S.ptr(sn + 1) - S.ptr(sn)) as Int;
        let ldc_ = ldf_ - sn_size_;
        let nb_ = S.blockSize();

        let local_reg_ = vec![0.0f64; sn_size_ as usize];
        let swaps_ = vec![0 as Int; sn_size_ as usize];
        let pivot_2x2_ = vec![0.0f64; sn_size_ as usize];

        FormatHandler {
            S_: S,
            regul_: *regul,
            sn_: sn,
            nb_,
            sn_size_,
            ldf_,
            ldc_,
            frontal_: frontal,
            clique_: Vec::new(),
            clique_ptr_: clique_ptr,
            local_reg_,
            swaps_,
            pivot_2x2_,
        }
    }

    /// HiGHS: FormatHandler.cpp::FormatHandler::terminate (lines 25-44)
    pub fn terminate(
        &mut self,
        clique: &mut Vec<f64>,
        total_reg: &mut Vec<f64>,
        swaps: &mut Vec<Int>,
        pivot_2x2: &mut Vec<f64>,
    ) {
        // Move local copies of data into their final position.
        // In this way, the shared objects sn_columns_ and schur_contribution_ are
        // accessed only here, while a local copy is used for the assembly and dense
        // factorisation. This should avoid the problem of false sharing.

        *clique = std::mem::take(&mut self.clique_);
        *swaps = std::mem::take(&mut self.swaps_);
        *pivot_2x2 = std::mem::take(&mut self.pivot_2x2_);

        // Move local regularisation into total regularisation.
        for i in 0..self.sn_size_ {
            total_reg[(self.S_.snStart(self.sn_) + i) as usize] = self.local_reg_[i as usize];
        }

        // This function should not require a lock, since all threads access different
        // locations within total_reg, sn_columns and schur_contribution.
    }
}

/// HiGHS: FormatHandler.h `class FormatHandler` (virtual methods)
pub trait FormatHandlerOps<'a> {
    fn base(&self) -> &FormatHandler<'a>;
    fn base_mut(&mut self) -> &mut FormatHandler<'a>;
    fn initFrontal(&mut self);
    fn initClique(&mut self);
    fn assembleFrontal(&mut self, i: Int, j: Int, val: f64);
    /// `child` starts at the child clique.
    #[allow(clippy::too_many_arguments)]
    fn assembleFrontalMultiple(
        &mut self,
        num: &mut Int,
        child: &[f64],
        nc: Int,
        child_sn: Int,
        row: Int,
        col: Int,
        i: Int,
        j: Int,
    );
    /// The child clique is `stack[child_off..]`.
    fn assembleClique(&mut self, stack: &mut [f64], child_off: Int64, nc: Int, child_sn: Int);
    fn denseFactorise(&mut self, stack: &mut [f64], reg_thresh: f64) -> Int;
    /// HiGHS: FormatHandler.h `virtual void extremeEntries() {}`
    fn extremeEntries(&mut self) {}
}
