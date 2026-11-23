//! HiGHS source: ipm/ipx/forrest_tomlin.h, forrest_tomlin.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Used only when `lu_kernel > 0`; not on the HiGHS default path.

use super::control::{Control, sci2};
use super::indexed_vector::IndexedVector;
use super::lu_factorization::LuFactorization;
use super::lu_update::LuUpdate;
use super::sparse_matrix::{
    SparseMatrix, dot_column, normest_inverse, scatter_column, triangular_solve,
};
use super::utils::{inverse_perm, permute, permute_back};
use super::*;
use crate::util::highs_utils::cpp_max;
use std::rc::Rc;

/// HiGHS: ipx/forrest_tomlin.h::ForrestTomlin
pub struct ForrestTomlin {
    /// base member `LuUpdate::updates_`
    updates_: Int,
    control_: Rc<Control>,
    dim_: Int,
    lu_: Box<dyn LuFactorization>,
    rowperm_: Vec<Int>,
    colperm_: Vec<Int>,
    rowperm_inv_: Vec<Int>,
    colperm_inv_: Vec<Int>,
    dependent_cols_: Vec<Int>,
    L_: SparseMatrix,
    U_: SparseMatrix,
    R_: SparseMatrix,
    replaced_: Vec<Int>,
    replace_next_: Int,
    have_btran_: bool,
    have_ftran_: bool,
    fill_factor_: f64,
    pivottol_: f64,
    /// rusolve BASICLU time limit (not in the HiGHS reference): basiclu_time_start_ / basiclu_time_limit_
    basiclu_time_start_: f64,
    basiclu_time_limit_: f64,
    work_: Vector,
}

impl ForrestTomlin {
    pub const kMaxUpdates: Int = 5000;

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::ForrestTomlin (lines 9-14)
    ///
    /// C++ moves the `unique_ptr<LuFactorization>&` argument in.
    pub fn new(control: Rc<Control>, dim: Int, lu: Box<dyn LuFactorization>) -> Self {
        ForrestTomlin {
            updates_: 0,
            control_: control,
            dim_: dim,
            lu_: lu,
            rowperm_: Vec::new(),
            colperm_: Vec::new(),
            rowperm_inv_: Vec::new(),
            colperm_inv_: Vec::new(),
            dependent_cols_: Vec::new(),
            L_: SparseMatrix::new(),
            U_: SparseMatrix::new(),
            R_: SparseMatrix::new(),
            replaced_: Vec::new(),
            // PORT NOTE: C++ leaves `replace_next_` uninitialized; it is set by
            // ComputeEta before any read.
            replace_next_: -1,
            have_btran_: false,
            have_ftran_: false,
            fill_factor_: 0.0,
            pivottol_: 0.1,
            basiclu_time_start_: 0.0,
            basiclu_time_limit_: f64::INFINITY,
            work_: Vector::new(dim + Self::kMaxUpdates),
        }
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::SolvePermuted (lines 250-300)
    fn solve_permuted(&mut self, x: &mut [f64], trans: u8) {
        let num_updates = self.replaced_.len() as Int;
        let dim = self.dim_;
        debug_assert_eq!(self.U_.cols(), dim + num_updates);

        // Require num_updates elements workspace at end of lhs.
        debug_assert!(x.len() as Int >= dim + num_updates);

        if trans == b't' || trans == b'T' {
            // Move replaced entries to the end of the pivot sequence and zero out
            // their old position. Because the corresponding columns of U are unit
            // columns now, the replaced positions remain zero when solving with U'.
            // This is crucial because we neglected to zero out the corresponding
            // rows of U in the update, so there are entries in U which actually
            // should not be there.
            for k in 0..num_updates {
                let r = self.replaced_[k as usize] as usize;
                x[(dim + k) as usize] = x[r];
                x[r] = 0.0;
            }
            triangular_solve(&self.U_, x, b't', "upper", 0);
            // Solve backwards with row eta matrices (leading scatter operations)
            // and put the entry from the end of the pivot sequence back into the
            // position that its update replaced.
            for k in (0..num_updates).rev() {
                let alpha = -x[(dim + k) as usize];
                scatter_column(&self.R_, k, alpha, x);
                let r = self.replaced_[k as usize] as usize;
                x[r] = x[(dim + k) as usize];
                x[(dim + k) as usize] = 0.0;
            }
            triangular_solve(&self.L_, x, b't', "lower", 1);
        } else {
            triangular_solve(&self.L_, x, b'n', "lower", 1);
            // Solve forward with row eta matrices (leading gather operations) and
            // put the newly computed entry at the end of the pivot sequence. We
            // actually would not need to zero out the replaced positions here
            // because they will not be accessed by any subsequent eta matrix.
            for k in 0..num_updates {
                let r = self.replaced_[k as usize] as usize;
                let v = x[r] - dot_column(&self.R_, k, x);
                x[(dim + k) as usize] = v;
                x[r] = 0.0;
            }
            // The triangular solve with U fills the replaced position with garbage.
            // This does not hurt because by the fact that the corresponding columns
            // of U have no nonzero off-diagonal entries, the computed values are
            // not propagated further. We finally overwrite them with their actual
            // values from the end of the pivot sequence.
            triangular_solve(&self.U_, x, b'n', "upper", 0);
            for k in (0..num_updates).rev() {
                let r = self.replaced_[k as usize] as usize;
                x[r] = x[(dim + k) as usize];
                x[(dim + k) as usize] = 0.0;
            }
        }
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::ComputeSpike (lines 302-325)
    fn compute_spike(&mut self, nb: Int, bi: &[Int], bx: &[f64]) {
        let num_updates = self.replaced_.len() as Int;
        let dim = self.dim_;

        // Solve L*lhs=b.
        self.work_.fill(0.0);
        for p in 0..nb {
            let r = self.rowperm_inv_[bi[p as usize] as usize] as usize;
            self.work_[r] = bx[p as usize];
        }
        triangular_solve(&self.L_, &mut self.work_, b'n', "lower", 1);

        // Apply update etas.
        for k in 0..num_updates {
            let r = self.replaced_[k as usize] as usize;
            let v = self.work_[r] - dot_column(&self.R_, k, &self.work_);
            self.work_[(dim + k) as usize] = v;
            self.work_[r] = 0.0;
        }

        // Store spike in U. Indices are sorted, which is required for the sparse
        // dot product in Update().
        self.U_.clear_queue();
        for p in 0..(dim + num_updates) {
            if self.work_[p as usize] != 0.0 {
                self.U_.push_back(p, self.work_[p as usize]);
            }
        }
        self.have_ftran_ = true;
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::ComputeEta (lines 327-358)
    fn compute_eta(&mut self, j: Int) {
        let num_updates = self.replaced_.len() as Int;
        let dim = self.dim_;
        debug_assert_eq!(self.U_.cols(), dim + num_updates);

        // Find permuted position of j.
        let mut pos = self.colperm_inv_[j as usize];
        for k in 0..num_updates {
            if self.replaced_[k as usize] == pos {
                pos = dim + k;
            }
        }

        // Solve lhs'U=e_pos'. Replaced position must remain zero (by the fact that
        // the corresponding columns of U are unit), so replaced indices cannot
        // occur in the new eta matrix.
        self.work_.fill(0.0);
        self.work_[pos as usize] = 1.0;
        triangular_solve(&self.U_, &mut self.work_, b't', "upper", 0);
        for k in 0..num_updates {
            debug_assert!(self.work_[self.replaced_[k as usize] as usize] == 0.0);
        }

        // Queue eta at end of R. Indices are sorted, which is required for the
        // sparse dot product in Update().
        self.R_.clear_queue();
        let pivot = self.work_[pos as usize];
        for i in (pos + 1)..(dim + num_updates) {
            if self.work_[i as usize] != 0.0 {
                self.R_.push_back(i, -self.work_[i as usize] / pivot);
            }
        }
        self.have_btran_ = true;
        self.replace_next_ = pos;
    }
}

/// HiGHS: ipx/forrest_tomlin.cc::SparseDot (static, lines 122-141)
///
/// Dot product of entries in queues of A1 and A2. Indices must be sorted in
/// increasing order.
fn sparse_dot(A1: &SparseMatrix, A2: &SparseMatrix) -> f64 {
    let q1 = A1.queue_size();
    let q2 = A2.queue_size();
    let mut p1: Int = 0;
    let mut p2: Int = 0;
    let mut d = 0.0;
    while p1 < q1 && p2 < q2 {
        let i1 = A1.qindex(p1);
        let i2 = A2.qindex(p2);
        if i1 == i2 {
            d += A1.qvalue(p1) * A2.qvalue(p2);
            p1 += 1;
            p2 += 1;
        } else if i1 < i2 {
            p1 += 1;
        } else {
            p2 += 1;
        }
    }
    d
}

impl LuUpdate for ForrestTomlin {
    fn oracle_stats(&self) -> super::lu_update::LuOracleStats {
        super::lu_update::LuOracleStats::ForrestTomlin {
            ft_updates: self.replaced_.len(),
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

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_Factorize (lines 16-51)
    fn _factorize(
        &mut self,
        Bbegin: &[Int],
        Bend: &[Int],
        Bi: &[Int],
        Bx: &[f64],
        strict_abs_pivottol: bool,
    ) -> Int {
        // Reset updates.
        let dim = self.dim_;
        self.R_.resize(dim, 0, 0);
        self.replaced_.clear();
        self.replace_next_ = -1;
        self.have_btran_ = false;
        self.have_ftran_ = false;

        self.lu_.factorize(
            dim,
            Bbegin,
            Bend,
            Bi,
            Bx,
            self.pivottol_,
            strict_abs_pivottol,
            &mut self.L_,
            &mut self.U_,
            &mut self.rowperm_,
            &mut self.colperm_,
            &mut self.dependent_cols_,
        );
        self.rowperm_inv_ = inverse_perm(&self.rowperm_);
        self.colperm_inv_ = inverse_perm(&self.colperm_);

        // Compute fill factor.
        let mut bnz: Int = 0;
        for j in 0..dim {
            bnz += Bend[j as usize] - Bbegin[j as usize];
        }
        self.fill_factor_ = 1.0 * ((self.L_.entries() + self.U_.entries()) as f64) / (bnz as f64);

        if self.control_.debug_enabled(3) {
            let normLinv = normest_inverse(&self.L_, "lower", 1);
            let normUinv = normest_inverse(&self.U_, "upper", 0);
            self.control_.debug(
                3,
                &format!(
                    " normLinv = {}, normUinv = {}, stability = {}\n",
                    sci2(normLinv),
                    sci2(normUinv),
                    sci2(self.lu_.stability())
                ),
            );
        }

        let mut ret: Int = 0;
        if self.lu_.stability() > kLuStabilityThreshold {
            ret |= 1;
        }
        if !self.dependent_cols_.is_empty() {
            ret |= 2;
        }
        ret
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_GetFactors (lines 53-63)
    fn _get_factors(
        &mut self,
        L: Option<&mut SparseMatrix>,
        U: Option<&mut SparseMatrix>,
        rowperm: Option<&mut [Int]>,
        colperm: Option<&mut [Int]>,
        dependent_cols: Option<&mut Vec<Int>>,
    ) {
        if let Some(l) = L {
            *l = self.L_.clone();
        }
        if let Some(u) = U {
            *u = self.U_.clone();
        }
        // C++ `std::copy` into caller buffers of at least dim entries: only
        // the first dim are written (a longer slice is allowed).
        if let Some(rp) = rowperm {
            rp[..self.rowperm_.len()].copy_from_slice(&self.rowperm_);
        }
        if let Some(cp) = colperm {
            cp[..self.colperm_.len()].copy_from_slice(&self.colperm_);
        }
        if let Some(dc) = dependent_cols {
            *dc = self.dependent_cols_.clone();
        }
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_SolveDense (lines 65-76)
    fn _solve_dense(&mut self, rhs: &Vector, lhs: &mut Vector, trans: u8) {
        if trans == b't' || trans == b'T' {
            permute_back(&self.colperm_, rhs, &mut self.work_);
            let mut w = std::mem::take(&mut self.work_);
            self.solve_permuted(&mut w, b'T');
            self.work_ = w;
            permute(&self.rowperm_, &self.work_, lhs);
        } else {
            permute_back(&self.rowperm_, rhs, &mut self.work_);
            let mut w = std::mem::take(&mut self.work_);
            self.solve_permuted(&mut w, b'N');
            self.work_ = w;
            permute(&self.colperm_, &self.work_, lhs);
        }
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_FtranForUpdate(nz, bi, bx) (lines 78-80)
    fn _ftran_for_update(&mut self, nz: Int, bi: &[Int], bx: &[f64]) {
        self.compute_spike(nz, bi, bx);
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_FtranForUpdate(nz, bi, bx, lhs) (lines 82-96)
    fn _ftran_for_update_lhs(&mut self, nb: Int, bi: &[Int], bx: &[f64], lhs: &mut IndexedVector) {
        let dim = self.dim_;
        self.compute_spike(nb, bi, bx);
        triangular_solve(&self.U_, &mut self.work_, b'n', "upper", 0);

        // Move extra variables from updates to replaced positions.
        let num_updates = self.replaced_.len() as Int;
        for k in (0..num_updates).rev() {
            let r = self.replaced_[k as usize] as usize;
            self.work_[r] = self.work_[(dim + k) as usize];
        }

        // Return lhs without pattern.
        for p in 0..dim {
            let cp = self.colperm_[p as usize];
            *lhs.get_mut(cp) = self.work_[p as usize];
        }
        lhs.invalidate_pattern();
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_BtranForUpdate(p) (lines 98-100)
    fn _btran_for_update(&mut self, p: Int) {
        self.compute_eta(p);
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_BtranForUpdate(p, lhs) (lines 102-118)
    fn _btran_for_update_lhs(&mut self, p: Int, lhs: &mut IndexedVector) {
        let dim = self.dim_;
        self.compute_eta(p);

        // Apply update etas and solve with L'.
        let num_updates = self.replaced_.len() as Int;
        for k in (0..num_updates).rev() {
            let alpha = -self.work_[(dim + k) as usize];
            scatter_column(&self.R_, k, alpha, &mut self.work_);
            let r = self.replaced_[k as usize] as usize;
            self.work_[r] = self.work_[(dim + k) as usize];
            self.work_[(dim + k) as usize] = 0.0;
        }
        triangular_solve(&self.L_, &mut self.work_, b't', "lower", 1);

        // Return lhs without pattern.
        for pp in 0..dim {
            let rp = self.rowperm_[pp as usize];
            *lhs.get_mut(rp) = self.work_[pp as usize];
        }
        lhs.invalidate_pattern();
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_Update (lines 143-210)
    fn _update(&mut self, pivot: f64) -> Int {
        let dim = self.dim_;
        let num_updates = self.replaced_.len() as Int;
        debug_assert!(self.have_ftran_);
        debug_assert!(self.have_btran_);

        // Find the position of the entry with row index replace_next_ in spike.
        // If not present, we will have where == qend.
        debug_assert!(self.replace_next_ >= 0);
        let qend = self.U_.queue_size();
        let mut where_: Int = 0;
        while where_ < qend && self.U_.qindex(where_) != self.replace_next_ {
            where_ += 1;
        }

        // Compute new diagonal entry of U. newdiag1 will be inserted into U;
        // newdiag2 would be the same in exact arithmetic and is for monitoring
        // numerical stability.
        let olddiag = self.U_.value(self.U_.end(self.replace_next_) - 1);
        let newdiag1 = pivot * olddiag;
        let qv = if where_ == qend {
            0.0
        } else {
            self.U_.qvalue(where_)
        };
        let newdiag2 = qv - sparse_dot(&self.U_, &self.R_);
        let newdiag_err = (newdiag1 - newdiag2).abs();
        let rel_newdiag_err = newdiag_err / newdiag1.abs();

        // Put new diagonal entry at end of spike.
        if where_ < qend {
            for l in where_..(qend - 1) {
                let qi = self.U_.qindex(l + 1);
                *self.U_.qindex_mut(l) = qi;
                let qval = self.U_.qvalue(l + 1);
                *self.U_.qvalue_mut(l) = qval;
            }
            *self.U_.qindex_mut(qend - 1) = dim + num_updates;
            *self.U_.qvalue_mut(qend - 1) = newdiag1;
        } else {
            self.U_.push_back(dim + num_updates, newdiag1);
        }

        // Overwrite replaced column by unit column in U.
        let end = self.U_.end(self.replace_next_);
        let begin = self.U_.begin(self.replace_next_);
        for l in begin..(end - 1) {
            *self.U_.value_mut(l) = 0.0;
        }
        *self.U_.value_mut(end - 1) = 1.0;

        // Finish update.
        self.U_.add_column();
        self.R_.add_column();
        self.replaced_.push(self.replace_next_);
        self.replace_next_ = -1;
        self.have_btran_ = false;
        self.have_ftran_ = false;

        if newdiag1 == 0.0 {
            return -1;
        }

        // Print a debugging message if a new eta entry is large.
        let mut max_eta = 0.0;
        for l in self.R_.begin(num_updates)..self.R_.end(num_updates) {
            max_eta = cpp_max(max_eta, self.R_.value(l).abs());
        }
        if max_eta > 1e10 {
            self.control_
                .debug(3, &format!(" max eta = {}\n", sci2(max_eta)));
        }

        // stability check
        if rel_newdiag_err > kFtDiagErrorTol {
            self.control_.debug(
                3,
                &format!(
                    " relative error in new diagonal entry of U = {}\n",
                    sci2(rel_newdiag_err)
                ),
            );
            return 1;
        }
        0
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_NeedFreshFactorization (lines 212-236)
    fn _need_fresh_factorization(&mut self) -> bool {
        let num_updates = self.replaced_.len() as Int;
        let Rnz = self.R_.entries(); // nnz in accumulated row etas
        let Lnz = self.L_.entries() + self.dim_; // nnz(L) incl. diagonal
        let Unz = self.U_.entries(); // nnz(U) incl. zeroed out columns
        let U0nz = self.U_.begin(self.dim_); // nnz(U) after factorization

        if num_updates == Self::kMaxUpdates {
            return true;
        }
        if num_updates < 100 {
            return false;
        }
        if (Rnz as f64) > 1.0 * (Lnz as f64) {
            return true;
        }
        if (Unz as f64) > 1.7 * (U0nz as f64) {
            return true;
        }
        false
    }

    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_fill_factor (lines 238-240)
    fn _fill_factor(&self) -> f64 {
        self.fill_factor_
    }
    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_pivottol() (lines 242-244)
    fn _pivottol(&self) -> f64 {
        self.pivottol_
    }
    /// HiGHS: ipx/forrest_tomlin.cc::ForrestTomlin::_pivottol(double) (lines 246-248)
    fn _set_pivottol(&mut self, new_pivottol: f64) {
        self.pivottol_ = new_pivottol;
    }
    /// rusolve BASICLU time limit (not in the HiGHS reference): ipx/forrest_tomlin.cc::ForrestTomlin::_timeStart
    fn _time_start(&mut self, new_time_start: f64) {
        self.basiclu_time_start_ = new_time_start;
    }
    /// rusolve BASICLU time limit (not in the HiGHS reference): ipx/forrest_tomlin.cc::ForrestTomlin::_timeLimit
    fn _time_limit(&mut self, new_time_limit: f64) {
        self.basiclu_time_limit_ = new_time_limit;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Queue-only matrix: `sparse_dot` reads the queues, not the CSC arrays.
    fn queued(entries: &[(Int, f64)]) -> SparseMatrix {
        let mut A = SparseMatrix::default();
        for &(i, x) in entries {
            A.push_back(i, x);
        }
        A
    }

    #[test]
    fn sparse_dot_merges_sorted_queues() {
        // Overlap at indices 2 and 5 only: 3*4 + 1*(-2) = 10.
        let A1 = queued(&[(0, 2.0), (2, 3.0), (5, 1.0)]);
        let A2 = queued(&[(2, 4.0), (3, 1.0), (5, -2.0)]);
        assert_eq!(sparse_dot(&A1, &A2), 10.0);
    }

    #[test]
    fn sparse_dot_disjoint_and_empty() {
        let A1 = queued(&[(0, 2.0), (4, 3.0)]);
        let A2 = queued(&[(1, 5.0), (3, 7.0)]);
        assert_eq!(sparse_dot(&A1, &A2), 0.0);
        assert_eq!(sparse_dot(&A1, &SparseMatrix::default()), 0.0);
    }
}
