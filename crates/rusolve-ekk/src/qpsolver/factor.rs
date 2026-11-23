//! HiGHS source: qpsolver/factor.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//!
//!
//! PORT NOTE (ownership, lib.rs rule (a)): C++ `CholeskyFactor` stores
//! `Runtime& runtime` and `Basis& basis`. Here they are parameters of the
//! methods that read them: the constructor (`rt.instance.num_var`,
//! `basis.getnuminactive()`), `recompute` (`Q.mat.extractcol`,
//! `basis.Ztprod`, `basis.getinactive`), `solveL` (line 190 may call
//! `recompute`) and `solve` (line 222 reads `runtime.instance.num_var` and
//! may call `recompute`).
//! NOTE: `eliminate` loops `k < current_k` on the MEMBER, not on its
//! `currentk` parameter; keep that.

use crate::prelude::*;
use crate::qpsolver::*;
use crate::util::c_printf::c_fmt_f;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// HiGHS: factor.hpp::CholeskyFactor (lines 20-418)
#[derive(Clone, Debug, Default)]
pub struct CholeskyFactor {
    uptodate: bool,
    numberofreduces: HighsInt,

    current_k: HighsInt,
    current_k_max: HighsInt,
    L: Vec<f64>,

    has_negative_eigenvalue: bool,
    a: Vec<f64>,
    /// HiGHS 1.16 #3325 (h16_qp_convexity): set from `rt.settings` in `new`; false
    /// reproduces 1.15.1. Read by `resize`, `expand`.
    h16_qp_convexity: bool,
}

impl CholeskyFactor {
    /// HiGHS 1.16 (#3458): factor.hpp::CholeskyFactor::resize (lines 36-56)
    /// HiGHS: factor.hpp::CholeskyFactor::resize (lines 36-51)
    fn resize(&mut self, new_k_max: HighsInt) {
        // HiGHS 1.16 #3458 (h16_qp_convexity): factor.hpp lines 37-42 (1.16.0)
        let new_k_max = if self.h16_qp_convexity {
            cpp_max(new_k_max, 1)
        } else {
            new_k_max
        };
        let l_old = self.L.clone();
        // C++ `L.clear(); L.resize(n)` keeps the vector's capacity: slots past
        // n keep their stale values, and expand() writes L[0] when
        // current_k_max is 0 (UB in C++ that lands in that spare capacity).
        // So L never shrinks here and only its first n entries are zeroed.
        let n = (new_k_max * new_k_max) as usize;
        if self.L.len() < n {
            self.L.resize(n, 0.0);
        }
        self.L[..n].fill(0.0);
        let l_size = n as HighsInt;
        // Driven by #958, changes made in following lines to avoid array
        // bound error when new_k_max < current_k_max
        let min_k_max = cpp_min(new_k_max, self.current_k_max);
        for i in 0..min_k_max {
            for j in 0..min_k_max {
                debug_assert!(i * new_k_max + j < l_size);
                self.L[(i * new_k_max + j) as usize] = l_old[(i * self.current_k_max + j) as usize];
            }
        }
        self.current_k_max = new_k_max;
    }

    /// HiGHS: factor.hpp::CholeskyFactor::CholeskyFactor (lines 54-60)
    pub fn new(rt: &Runtime, bas: &Basis) -> Self {
        let mut s = CholeskyFactor {
            uptodate: false,
            numberofreduces: 0,
            current_k: 0,
            current_k_max: 0,
            L: Vec::new(),
            has_negative_eigenvalue: false,
            a: Vec::new(),
            h16_qp_convexity: rt.settings.h16_qp_convexity,
        };
        s.uptodate = false;
        s.current_k_max = cpp_max(
            cpp_min(
                (rt.instance.num_var as f64 / 16.0).ceil() as HighsInt,
                1000 as HighsInt,
            ),
            bas.getnuminactive(),
        );
        s.L.resize((s.current_k_max * s.current_k_max) as usize, 0.0);
        s
    }

    /// HiGHS: factor.hpp::CholeskyFactor::recompute (lines 62-108)
    pub fn recompute(&mut self, rt: &Runtime, basis: &mut Basis) -> QpSolverStatus {
        let mut orig: Vec<Vec<f64>> = Vec::new();
        let dim_ns = basis.getinactive().len() as HighsInt;
        self.numberofreduces = 0;

        orig.resize(dim_ns as usize, vec![0.0; dim_ns as usize]);
        self.resize(dim_ns);

        let mut temp = Matrix::new(dim_ns, 0);

        let mut buffer_qcol = QpVector::new(rt.instance.num_var);
        let mut buffer_ztqi = QpVector::new(dim_ns);
        for i in 0..rt.instance.num_var {
            rt.instance.Q.mat.extractcol_into(i, &mut buffer_qcol);
            basis.Ztprod(&buffer_qcol, &mut buffer_ztqi, false, -1);
            temp.append(&buffer_ztqi);
        }
        let temp_t = temp.t();
        for i in 0..dim_ns {
            temp_t.extractcol_into(i, &mut buffer_qcol);
            basis.Ztprod(&buffer_qcol, &mut buffer_ztqi, false, -1);
            for j in 0..buffer_ztqi.num_nz {
                orig[i as usize][buffer_ztqi.index[j as usize] as usize] =
                    buffer_ztqi.value[buffer_ztqi.index[j as usize] as usize];
            }
        }

        let ckm = self.current_k_max as usize;
        for col in 0..orig.len() {
            for row in 0..=col {
                let mut sum = 0.0;
                if row == col {
                    for k in 0..row {
                        sum += self.L[k * ckm + row] * self.L[k * ckm + row];
                    }
                    let d_value = orig[row][row] - sum;
                    if d_value <= 0.0 {
                        return QpSolverStatus::NOTPOSITIVDEFINITE;
                    }
                    self.L[row * ckm + row] = d_value.sqrt();
                } else {
                    for k in 0..row {
                        sum += self.L[k * ckm + col] * self.L[k * ckm + row];
                    }
                    self.L[row * ckm + col] = (orig[col][row] - sum) / self.L[row * ckm + row];
                }
            }
        }
        self.current_k = dim_ns;
        self.uptodate = true;

        QpSolverStatus::OK
    }

    /// HiGHS 1.16 (#3325): factor.hpp::CholeskyFactor::invalidate (line 68).
    /// Called only on the on path (quass/solve.rs).
    pub fn invalidate(&mut self) {
        self.uptodate = false;
    }

    /// HiGHS 1.16 (#3325): factor.hpp::CholeskyFactor::expand (lines 127-156)
    /// HiGHS: factor.hpp::CholeskyFactor::expand (lines 110-185)
    ///
    /// `curvature` is p'Qp from the caller; it is read only when
    /// `h16_qp_convexity` is set, and the off path passes 0.0.
    pub fn expand(
        &mut self,
        yp: &QpVector,
        gyp: &mut QpVector,
        l: &mut QpVector,
        m: &mut QpVector,
        curvature: f64,
    ) -> QpSolverStatus {
        if !self.uptodate {
            return QpSolverStatus::OK;
        }
        let mu = gyp.dot(yp);
        l.resparsify();
        // HiGHS 1.16 #3325 (h16_qp_convexity): factor.hpp lines 134 (1.16.0)
        let lambda = if self.h16_qp_convexity {
            if l.dim > 0 { curvature } else { mu }
        } else {
            mu - l.norm2()
        };
        if lambda > 0.0 {
            if self.current_k_max <= self.current_k + 1 {
                self.resize(self.current_k_max * 2);
            }

            for i in 0..self.current_k {
                self.L[(i * self.current_k_max + self.current_k) as usize] = l.value[i as usize];
            }
            self.L[(self.current_k * self.current_k_max + self.current_k) as usize] = lambda.sqrt();

            self.current_k += 1;
        } else {
            return QpSolverStatus::NOTPOSITIVDEFINITE;
            // Commented out unreachable code (C++ lines 132-183)
        }
        QpSolverStatus::OK
    }

    /// HiGHS: factor.hpp::CholeskyFactor::solveL (lines 187-207)
    ///
    /// PORT NOTE (rule (a)): takes `rt` and `basis` because the `!uptodate`
    /// branch calls `recompute()` (line 190); the skeleton signature omitted
    /// them. The caller (quass.cpp:117) holds both.
    pub fn solve_l(
        &mut self,
        rt: &Runtime,
        basis: &mut Basis,
        rhs: &mut QpVector,
    ) -> QpSolverStatus {
        if !self.uptodate {
            // Recompute Cholesky
            let status = self.recompute(rt, basis);
            if status != QpSolverStatus::OK {
                return status;
            }
        }

        if self.current_k != rhs.dim {
            debug_assert!(self.current_k == rhs.dim);
            return QpSolverStatus::kError;
        }

        for r in 0..rhs.dim {
            for j in 0..r {
                rhs.value[r as usize] -=
                    rhs.value[j as usize] * self.L[(j * self.current_k_max + r) as usize];
            }

            rhs.value[r as usize] /= self.L[(r * self.current_k_max + r) as usize];
        }
        QpSolverStatus::OK
    }

    /// HiGHS: factor.hpp::CholeskyFactor::solveLT (lines 210-218)
    pub fn solve_lt(&mut self, rhs: &mut QpVector) {
        let mut i = rhs.dim - 1;
        while i >= 0 {
            let mut sum = 0.0;
            let mut j = rhs.dim - 1;
            while j > i {
                sum += rhs.value[j as usize] * self.L[(i * self.current_k_max + j) as usize];
                j -= 1;
            }
            rhs.value[i as usize] =
                (rhs.value[i as usize] - sum) / self.L[(i * self.current_k_max + i) as usize];
            i -= 1;
        }
    }

    /// HiGHS: factor.hpp::CholeskyFactor::solve (lines 220-234)
    pub fn solve(&mut self, rt: &Runtime, basis: &mut Basis, rhs: &mut QpVector) -> QpSolverStatus {
        if !self.uptodate
            || (self.numberofreduces >= rt.instance.num_var / 2 && !self.has_negative_eigenvalue)
        {
            // Recompute Cholesky
            let status = self.recompute(rt, basis);
            if status != QpSolverStatus::OK {
                return status;
            }
        }
        let status = self.solve_l(rt, basis, rhs);
        if status != QpSolverStatus::OK {
            return status;
        }
        self.solve_lt(rhs);

        rhs.resparsify();
        QpSolverStatus::OK
    }

    /// HiGHS: factor.hpp::CholeskyFactor::eliminate (lines 236-294)
    ///
    /// C++ passes the member `L` (and in `reduce`, `L` too): take it with
    /// `std::mem::take(&mut self.L)` around the call. The loop bound reads
    /// `self.current_k`, not `currentk`.
    pub fn eliminate(
        &self,
        m: &mut Vec<f64>,
        i: HighsInt,
        j: HighsInt,
        kmax: HighsInt,
        currentk: HighsInt,
    ) {
        // i = col, j = row
        if m[(j * kmax + i) as usize] == 0.0 {
            return;
        }
        let z = (m[(i * kmax + i) as usize] * m[(i * kmax + i) as usize]
            + m[(j * kmax + i) as usize] * m[(j * kmax + i) as usize])
            .sqrt();
        let (cos_, sin_) = if z == 0.0 {
            (1.0, 0.0)
        } else {
            (
                m[(i * kmax + i) as usize] / z,
                -m[(j * kmax + i) as usize] / z,
            )
        };

        if sin_ == 0.0 {
            if cos_ > 0.0 {
                // nothing
            } else {
                for k in 0..self.current_k {
                    // update entry i and j of column k
                    let a_ik = m[(i * kmax + k) as usize];
                    // entry i
                    m[(i * kmax + k) as usize] = -a_ik;
                    m[(j * kmax + k) as usize] = -m[(j * kmax + k) as usize];
                }
            }
        } else if cos_ == 0.0 {
            if sin_ > 0.0 {
                for k in 0..self.current_k {
                    // update entry i and j of column k
                    let a_ik = m[(i * kmax + k) as usize];
                    // entry i
                    m[(i * kmax + k) as usize] = -m[(j * kmax + k) as usize];
                    m[(j * kmax + k) as usize] = a_ik;
                }
            } else {
                for k in 0..self.current_k {
                    // update entry i and j of column k
                    let a_ik = m[(i * kmax + k) as usize];
                    // entry i
                    m[(i * kmax + k) as usize] = m[(j * kmax + k) as usize];
                    m[(j * kmax + k) as usize] = -a_ik;
                }
            }
        } else {
            // #pragma omp parallel for
            for k in 0..self.current_k {
                // update entry i and j of column k
                let a_ik = m[(i * kmax + k) as usize];
                // entry i
                m[(i * kmax + k) as usize] = cos_ * a_ik - sin_ * m[(j * kmax + k) as usize];
                m[(j * kmax + k) as usize] = sin_ * a_ik + cos_ * m[(j * kmax + k) as usize];
            }
        }
        m[(j * kmax + i) as usize] = 0.0;
    }

    /// HiGHS: factor.hpp::CholeskyFactor::reduce (lines 296-391)
    ///
    /// NOTE: `unsigned p = maxabsd;` (line 305): keep the unsigned arithmetic.
    pub fn reduce(&mut self, buffer_d: &QpVector, maxabsd: HighsInt, p_in_v: bool) {
        if self.current_k == 0 {
            return;
        }
        if !self.uptodate {
            return;
        }
        self.numberofreduces += 1;

        let p: u32 = maxabsd as u32; // col we push to the right and remove

        // start situation: p=3, current_k = 5
        // |1 x  | |x    |       |1   | |xxxxx|
        // | 1x  | |xx   |  ===  | 1  | | xxxx|
        // |  x1 | |xxx  |       |xxxx| |  xxx|
        // |  x 1| |xxxx |       |  1 | |   xx|
        //         |xxxxx|       |   1| |    x|
        // next step: move row/col p to the bottom/right

        //> save row p
        let mut row_p = vec![0.0; self.current_k as usize];
        for i in 0..self.current_k {
            row_p[i as usize] = self.L[(p as HighsInt * self.current_k_max + i) as usize];
        }

        //> move all rows > p up by one row
        let mut row = p as HighsInt;
        while row < self.current_k - 1 {
            for i in 0..self.current_k {
                self.L[(row * self.current_k_max + i) as usize] =
                    self.L[((row + 1) * self.current_k_max + i) as usize];
            }
            row += 1;
        }

        //> load row p
        for i in 0..self.current_k {
            self.L[((self.current_k - 1) * self.current_k_max + i) as usize] = row_p[i as usize];
        }

        //> now move col p to the right in each row
        for row in 0..self.current_k {
            let p_entry = self.L[(row * self.current_k_max + p as HighsInt) as usize];
            let mut col = p as HighsInt;
            while col < self.current_k - 1 {
                self.L[(row * self.current_k_max + col) as usize] =
                    self.L[(row * self.current_k_max + col + 1) as usize];
                col += 1;
            }
            self.L[(row * self.current_k_max + self.current_k - 1) as usize] = p_entry;
        }

        if self.current_k == 1 {
            self.current_k -= 1;
            return;
        }

        if !p_in_v {
            // situation now:
            // |1   x| |x    |       |1   | |xxxxx|
            // | 1  x| |xx   |  ===  | 1  | | xxxx|
            // |  1 x| |xxx x|       |  1 | |  xx |
            // |   1x| |xxxxx|       |   1| |   x |
            //         |xx  x|       |xxxx| |  xxx|
            // next: remove nonzero entries in last column except for diagonal element
            let mut r = p as HighsInt - 1;
            while r >= 0 {
                // to current_k-1
                let mut l = std::mem::take(&mut self.L);
                self.eliminate(
                    &mut l,
                    self.current_k - 1,
                    r,
                    self.current_k_max,
                    self.current_k,
                );
                self.L = l;
                r -= 1;
            }

            // situation now:
            // |1   x| |x   x|        |xxxx | |1   |
            // | 1  x| |xx  x|  ===   | xxx | | 1  |
            // |  1 x| |xxx x|        |  xx | |  1 |
            // |   1x| |xxxxx|        |   x | |   1|
            //         |    x|        |xxxxx| |xxxx|
            // next: multiply product
            // new last row: old last row (first current_k-1 elements) + r *
            // R_current_k_current_k

            for i in 0..buffer_d.num_nz {
                let idx = buffer_d.index[i as usize];
                if idx == maxabsd {
                    continue;
                }
                if idx < maxabsd {
                    self.L[((self.current_k - 1) * self.current_k_max + idx) as usize] +=
                        -buffer_d.value[idx as usize] / buffer_d.value[maxabsd as usize]
                            * self.L[((self.current_k - 1) * self.current_k_max + self.current_k
                                - 1) as usize];
                } else {
                    self.L[((self.current_k - 1) * self.current_k_max + idx - 1) as usize] +=
                        -buffer_d.value[idx as usize] / buffer_d.value[maxabsd as usize]
                            * self.L[((self.current_k - 1) * self.current_k_max + self.current_k
                                - 1) as usize];
                }
            }
            // situation now: as above, but no more product
        }
        // next: eliminate last row
        for i in 0..self.current_k - 1 {
            let mut l = std::mem::take(&mut self.L);
            self.eliminate(
                &mut l,
                i,
                self.current_k - 1,
                self.current_k_max,
                self.current_k,
            );
            self.L = l;
        }
        self.current_k -= 1;
    }

    /// HiGHS: factor.hpp::CholeskyFactor::report (lines 393-401)
    ///
    /// C++ default arguments: name = ""
    pub fn report(&self, name: &str) {
        println!("{}", name);
        for i in 0..self.current_k {
            for j in 0..self.current_k {
                print!(
                    "{} ",
                    c_fmt_f(self.L[(i * self.current_k_max + j) as usize], 6)
                );
            }
            println!();
        }
    }

    /// HiGHS: factor.hpp::CholeskyFactor::density (lines 403-417)
    pub fn density(&self) -> f64 {
        if self.current_k == 0 {
            return 0.0;
        }

        let mut num_nz: HighsInt = 0;
        for i in 0..self.current_k {
            for j in 0..self.current_k {
                if self.L[(i * self.current_k_max + j) as usize].abs() > 1e-7 {
                    num_nz += 1;
                }
            }
        }
        num_nz as f64 / ((self.current_k * (self.current_k + 1)) as f64 / 2.0)
    }
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
    use crate::qpsolver::*;

    /// min 0.5 (2 x0^2 + 4 x1^2) - 2 x0 - 4 x1, 0 <= x <= 10 and, when
    /// `with_row`, x0 + x1 <= 10. Optimum x = (1, 1), objective -3.
    fn box_qp(with_row: bool) -> Instance {
        let nc = if with_row { 1 } else { 0 };
        let mut inst = Instance::new(2, nc);
        inst.c = QpVector {
            num_nz: 2,
            dim: 2,
            index: vec![0, 1],
            value: vec![-2.0, -4.0],
        };
        inst.Q.mat = MatrixBase {
            num_row: 2,
            num_col: 2,
            start: vec![0, 1, 2],
            index: vec![0, 1],
            value: vec![2.0, 4.0],
        };
        inst.A.mat = if with_row {
            MatrixBase {
                num_row: 1,
                num_col: 2,
                start: vec![0, 1, 2],
                index: vec![0, 0],
                value: vec![1.0, 1.0],
            }
        } else {
            MatrixBase {
                num_row: 0,
                num_col: 2,
                start: vec![0, 0, 0],
                index: vec![],
                value: vec![],
            }
        };
        inst.con_lo = vec![-kHighsInf; nc as usize];
        inst.con_up = vec![10.0; nc as usize];
        inst.var_lo = vec![0.0, 0.0];
        inst.var_up = vec![10.0, 10.0];
        inst
    }

    /// Both variables active at their lower bound (constraint ids
    /// `num_con + i`); the row (id 0) is not in the basis.
    fn at_lower(rt: &mut Runtime) -> Basis {
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kActiveAtLower,
            BasisStatus::kActiveAtLower,
        ];
        Basis::new(rt, vec![1, 2], status, vec![])
    }

    /// Both variables in the nullspace (inactive in basis): Z = I.
    fn all_free(rt: &mut Runtime) -> Basis {
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kInactiveInBasis,
            BasisStatus::kInactiveInBasis,
        ];
        Basis::new(rt, vec![], status, vec![1, 2])
    }

    fn close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-9)
    }

    #[test]
    fn reduced_hessian_solve() {
        // Z = I, so the factor is the Cholesky factor of Q = diag(2, 4).
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut basis = all_free(&mut rt);
        let mut f = CholeskyFactor::new(&rt, &basis);
        assert_eq!(f.recompute(&rt, &mut basis), QpSolverStatus::OK);
        let mut rhs = QpVector::new(2);
        rhs.value = vec![2.0, 4.0];
        rhs.resparsify();
        assert_eq!(f.solve(&rt, &mut basis, &mut rhs), QpSolverStatus::OK);
        assert!(close(&rhs.value, &[1.0, 1.0]));
    }

    #[test]
    fn solve_l_recomputes_when_not_uptodate() {
        // factor.hpp:188-192: a fresh factor recomputes before solving L y = b.
        // L = diag(sqrt 2, 2), so y = (2 / sqrt 2, 4 / 2).
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut basis = all_free(&mut rt);
        let mut f = CholeskyFactor::new(&rt, &basis);
        let mut rhs = QpVector::new(2);
        rhs.value = vec![2.0, 4.0];
        rhs.resparsify();
        assert_eq!(f.solve_l(&rt, &mut basis, &mut rhs), QpSolverStatus::OK);
        assert_eq!(rhs.value, vec![2.0 / 2.0f64.sqrt(), 2.0]);
        assert_eq!(f.current_k, 2);
        assert!(f.uptodate);
    }

    #[test]
    fn recompute_not_positive_definite() {
        // factor.hpp:93-94: d_value = 0 on the first pivot.
        let mut inst = box_qp(true);
        inst.Q.mat.value[0] = 0.0;
        let mut rt = Runtime::new(&inst, Statistics::default());
        let mut basis = all_free(&mut rt);
        let mut f = CholeskyFactor::new(&rt, &basis);
        assert_eq!(
            f.recompute(&rt, &mut basis),
            QpSolverStatus::NOTPOSITIVDEFINITE
        );
        assert!(!f.uptodate);
    }

    #[test]
    fn expand_and_density() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut basis = at_lower(&mut rt);
        let mut f = CholeskyFactor::new(&rt, &basis);
        let mut yp = QpVector::new(2);
        yp.value = vec![1.0, 0.0];
        yp.resparsify();
        let mut gyp = QpVector::new(2);
        gyp.value = vec![2.0, 0.0];
        gyp.resparsify();
        let mut l = QpVector::new(0);
        let mut m = QpVector::new(0);
        // !uptodate: no-op returning OK (lines 112-114).
        assert_eq!(
            f.expand(&yp, &mut gyp, &mut l, &mut m, 0.0),
            QpSolverStatus::OK
        );
        assert_eq!(f.current_k, 0);
        assert_eq!(f.density(), 0.0);
        // Empty nullspace: recompute gives current_k = 0, current_k_max = 0.
        assert_eq!(f.recompute(&rt, &mut basis), QpSolverStatus::OK);
        assert_eq!(f.current_k_max, 0);
        // With current_k_max = 0, line 120 resizes to 0 * 2 = 0 and the write
        // at line 126 is out of bounds (UB in C++ too), so start from 1.
        f.resize(1);
        // lambda = 2 > 0: resize to 2 (line 120), L[0] = sqrt 2 (line 126).
        assert_eq!(
            f.expand(&yp, &mut gyp, &mut l, &mut m, 0.0),
            QpSolverStatus::OK
        );
        assert_eq!(f.current_k, 1);
        assert_eq!(f.current_k_max, 2);
        assert_eq!(f.L[0], 2.0f64.sqrt());
        assert_eq!(f.density(), 1.0);
        // lambda = mu - |l|^2 <= 0 -> NOTPOSITIVDEFINITE (line 130).
        let mut zero = QpVector::new(2);
        assert_eq!(
            f.expand(&zero.clone(), &mut zero, &mut l, &mut m, 0.0),
            QpSolverStatus::NOTPOSITIVDEFINITE
        );
        assert_eq!(f.current_k, 1);
    }

    #[test]
    fn reduce_moves_column_to_end() {
        // Z = I, L = diag(sqrt 2, 2). Removing p = 0 with p_in_v moves row
        // and column 0 to the end (lines 316-340); the last row has no
        // off-diagonal entry, so eliminate returns early (line 239).
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut basis = all_free(&mut rt);
        let mut f = CholeskyFactor::new(&rt, &basis);
        assert_eq!(f.recompute(&rt, &mut basis), QpSolverStatus::OK);
        let mut d = QpVector::new(2);
        d.value = vec![1.0, 0.0];
        d.resparsify();
        f.reduce(&d, 0, true);
        assert_eq!(f.current_k, 1);
        assert_eq!(f.numberofreduces, 1);
        assert_eq!(f.L[0], 2.0);
        assert_eq!(f.L[3], 2.0f64.sqrt());
    }

    #[test]
    fn eliminate_rotates_rows() {
        // Givens rotation of rows i = 0 and j = 1 on column 0 (lines 242-293):
        // z = 5, cos = 3 / 5, sin = -4 / 5.
        let mut f = CholeskyFactor {
            current_k: 2,
            current_k_max: 2,
            ..Default::default()
        };
        let mut m = vec![3.0, 1.0, 4.0, 2.0];
        f.eliminate(&mut m, 0, 1, 2, 2);
        let (cos_, sin_) = (3.0 / 5.0, -4.0 / 5.0);
        assert_eq!(m[0], cos_ * 3.0 - sin_ * 4.0);
        assert_eq!(m[1], cos_ * 1.0 - sin_ * 2.0);
        assert_eq!(m[2], 0.0);
        assert_eq!(m[3], sin_ * 1.0 + cos_ * 2.0);
        // m[j][i] == 0: early return (line 239).
        let mut m = vec![-1.0, 2.0, 0.0, 3.0];
        f.eliminate(&mut m, 0, 1, 2, 2);
        assert_eq!(m, vec![-1.0, 2.0, 0.0, 3.0]);
        // cos == 0, sin = -1 < 0: swap with sign (lines 275-281). The loop
        // bound is the member current_k = 1, not the currentk argument 2.
        f.current_k = 1;
        let mut m = vec![0.0, 5.0, 2.0, 7.0];
        f.eliminate(&mut m, 0, 1, 2, 2);
        assert_eq!(m, vec![2.0, 5.0, 0.0, 7.0]);
    }

    /// HiGHS 1.16 #3458: the 1.16 `resize` raises a zero `new_k_max` to 1; the off path
    /// keeps the 1.15.1 value 0 (the UB emulation documented above).
    #[test]
    fn h16_qpconv_resize_zero_gives_k_max_one() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        rt.settings.h16_qp_convexity = true;
        let basis = all_free(&mut rt);
        let mut f = CholeskyFactor::new(&rt, &basis);
        f.resize(0);
        assert_eq!(f.current_k_max, 1);
    }

    #[test]
    fn h16_qpconv_off_resize_zero_keeps_k_max_zero() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let basis = all_free(&mut rt);
        let mut f = CholeskyFactor::new(&rt, &basis);
        f.resize(0);
        assert_eq!(f.current_k_max, 0);
    }

    /// HiGHS 1.16 #3325: with the flag on, the new diagonal entry is the curvature
    /// p'Qp passed in (here 1.0), not mu - l'l (here 2 - 0 = 2). With the flag
    /// off the same call ignores the curvature argument.
    fn expand_with_curvature(flag: bool, curvature: f64) -> (HighsInt, f64) {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        rt.settings.h16_qp_convexity = flag;
        let mut basis = all_free(&mut rt);
        let mut f = CholeskyFactor::new(&rt, &basis);
        assert_eq!(f.recompute(&rt, &mut basis), QpSolverStatus::OK);
        let mut yp = QpVector::new(2);
        yp.value = vec![1.0, 0.0];
        yp.resparsify();
        let mut gyp = QpVector::new(2);
        gyp.value = vec![2.0, 0.0];
        gyp.resparsify();
        // l = 0 with dim 2: l'l = 0, and l.dim > 0 selects the curvature branch.
        let mut l = QpVector::new(2);
        let mut m = QpVector::new(2);
        assert_eq!(
            f.expand(&yp, &mut gyp, &mut l, &mut m, curvature),
            QpSolverStatus::OK
        );
        assert_eq!(f.current_k, 3);
        // New diagonal entry sits at row 2, column 2 of the k_max = 4 layout.
        (f.current_k_max, f.L[2 * 4 + 2])
    }

    #[test]
    fn h16_qpconv_expand_uses_curvature_when_on() {
        let (k_max, diag) = expand_with_curvature(true, 1.0);
        assert_eq!(k_max, 4);
        assert_eq!(diag, 1.0f64.sqrt());
    }

    #[test]
    fn h16_qpconv_off_expand_ignores_curvature() {
        let (k_max, diag) = expand_with_curvature(false, 1.0);
        assert_eq!(k_max, 4);
        assert_eq!(diag, 2.0f64.sqrt());
    }
}
