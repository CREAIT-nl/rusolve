//! HiGHS source: ipm/hipo/factorhighs/FactorHiGHS.h, ipm/hipo/factorhighs/FactorHiGHS.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `FHsolver`, the public interface of FactorHiGHS. `logger_` is a clone of
//! the IPM logger. The destructor (lines 23-28) only calls the empty
//! `printTimes` / `printIter` and is not ported. `solve` and
//! `getRegularisation` take the `Symbolic` the C++ `Numeric` points at.

use super::super::auxiliary::logger::Logger;
use super::super::{Int, Int64};
use super::analyse::Analyse;
use super::clique_stack::CliqueStack;
use super::data_collector::DataCollector;
use super::factorise::Factorise;
use super::numeric::Numeric;
use super::symbolic::Symbolic;
use super::*;

/// HiGHS: FactorHiGHS.h `class FHsolver`
#[derive(Clone)]
pub struct FHsolver {
    pub logger_: Option<Logger>,
    pub data_: DataCollector,
    pub regul_: Regul,
    pub N_: Numeric,
    pub serial_stack_: CliqueStack,
    pub nb_: Int,
    pub sn_columns_: Vec<Vec<f64>>,
    /// rusolve lever `ipm_hipo_cache_perm`: replayed permutation of the matrix structure.
    pub perm_cache_: crate::ipm::hipo::factorhighs::factorise::PermCache,
    /// rusolve: deadline polled between supernodes; `timed_out_` is set when it stopped a factorisation
    pub deadline_: Option<std::time::Instant>,
    pub timed_out_: bool,
}

impl FHsolver {
    /// HiGHS: FactorHiGHS.h `static const Int default_nb_ = 128`
    pub const default_nb_: Int = 128;

    /// HiGHS: FactorHiGHS.cpp::FHsolver::FHsolver (lines 10-21)
    /// C++ default arguments: `logger = nullptr`, `block_size = default_nb_`.
    pub fn new(logger: Option<&Logger>, block_size: Int) -> Self {
        // C++ lines 12-20: #ifdef HIPO_COLLECT_EXPENSIVE_DATA and
        // #if HIPO_TIMING_LEVEL > 0 are not defined at level 0: dropped.
        Self {
            logger_: logger.cloned(),
            data_: DataCollector::new(),
            regul_: Regul::default(),
            N_: Numeric::default(),
            serial_stack_: CliqueStack::default(),
            nb_: if block_size > 0 {
                block_size
            } else {
                Self::default_nb_
            },
            sn_columns_: Vec::new(),
            perm_cache_: Default::default(),
            deadline_: None,
            timed_out_: false,
        }
    }

    /// HiGHS: FactorHiGHS.cpp::FHsolver::newIter (lines 30-30)
    pub fn newIter(&mut self) {
        self.data_.append();
    }

    /// HiGHS: FactorHiGHS.cpp::FHsolver::setRegularisation (lines 32-35)
    pub fn setRegularisation(&mut self, reg_p: f64, reg_d: f64) {
        self.regul_.primal = reg_p;
        self.regul_.dual = reg_d;
    }

    /// HiGHS: FactorHiGHS.cpp::FHsolver::analyse (lines 37-43)
    pub fn analyse(
        &mut self,
        S: &mut Symbolic,
        rows: &[Int],
        ptr: &[Int],
        signs: &[Int],
        perm: &[Int],
    ) -> Int {
        let mut an_obj = Analyse::new(
            rows,
            ptr,
            signs,
            self.nb_,
            self.logger_.as_ref(),
            &self.data_,
            perm,
        );
        an_obj.run(S)
    }

    /// HiGHS: FactorHiGHS.cpp::FHsolver::factorise (lines 45-51)
    pub fn factorise(&mut self, S: &Symbolic, rows: &[Int], ptr: &[Int], vals: &[f64]) -> Int {
        let use_cache =
            crate::ipm::hipo::blas::CACHE_PERM.load(std::sync::atomic::Ordering::Relaxed);
        let mut fact_obj = Factorise::new_cached(
            S,
            rows,
            ptr,
            vals,
            &self.regul_,
            self.logger_.as_ref(),
            &self.data_,
            &mut self.sn_columns_,
            Some(&mut self.serial_stack_),
            if use_cache {
                Some(&mut self.perm_cache_)
            } else {
                None
            },
        );
        fact_obj.deadline_ = self.deadline_;
        let failed = fact_obj.run(&mut self.N_);
        self.timed_out_ = fact_obj
            .timed_out_
            .load(std::sync::atomic::Ordering::Relaxed);
        // C++: return fact_obj.run(N_); -- bool implicitly converts to Int
        // (true = 1 = failure, false = 0 = kRetOk).
        if failed { 1 } else { 0 }
    }

    /// HiGHS: FactorHiGHS.cpp::FHsolver::solve (lines 53-53)
    pub fn solve(&mut self, S: &Symbolic, x: &mut Vec<f64>) -> Int {
        self.N_.solve(S, &self.sn_columns_, &self.data_, x)
    }

    /// rusolve (ipm_hipo_dense_schur): `solve` through a shared borrow, so
    /// several right-hand sides can be solved at once against one factor (Numeric::solve is `&self`).
    pub fn solve_shared(&self, S: &Symbolic, x: &mut Vec<f64>) -> Int {
        self.N_.solve(S, &self.sn_columns_, &self.data_, x)
    }

    /// HiGHS: FactorHiGHS.cpp::FHsolver::getRegularisation (lines 55-55)
    pub fn getRegularisation(&mut self, S: &Symbolic, reg: &mut Vec<f64>) {
        self.N_.getReg(S, reg);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn solve_small_quasidefinite() {
        // [2 1; 1 -3] lower CSC, signs (+, -), identity ordering
        let (rows, ptr, vals) = (vec![0, 1, 1], vec![0, 2, 3], vec![2.0, 1.0, -3.0]);
        let mut fh = FHsolver::new(None, FHsolver::default_nb_);
        let mut s = Symbolic::new();
        assert_eq!(fh.analyse(&mut s, &rows, &ptr, &[1, -1], &[0, 1]), kRetOk);
        assert_eq!(fh.factorise(&s, &rows, &ptr, &vals), kRetOk);
        let mut x = vec![3.0, -2.0];
        fh.solve(&s, &mut x);
        assert!((x[0] - 1.0).abs() < 1e-12 && (x[1] - 1.0).abs() < 1e-12);
    }

    /// Nested-dissection order of a `nr` x `nc` grid block (separator last).
    fn nd_order(r0: Int, r1: Int, c0: Int, c1: Int, ncols: Int, out: &mut Vec<Int>) {
        if (r1 - r0) * (c1 - c0) <= 4 {
            for r in r0..r1 {
                for c in c0..c1 {
                    out.push(r * ncols + c);
                }
            }
            return;
        }
        if r1 - r0 >= c1 - c0 {
            let m = (r0 + r1) / 2;
            nd_order(r0, m, c0, c1, ncols, out);
            nd_order(m + 1, r1, c0, c1, ncols, out);
            nd_order(m, m + 1, c0, c1, ncols, out);
        } else {
            let m = (c0 + c1) / 2;
            nd_order(r0, r1, c0, m, ncols, out);
            nd_order(r0, r1, m + 1, c1, ncols, out);
            nd_order(r0, r1, m, m + 1, ncols, out);
        }
    }

    #[test]
    fn solve_grid_laplacian_nested_dissection() {
        // 5-point Laplacian (+ shift) on a 40 x 40 grid, lower CSC, nested
        // dissection ordering: many supernodes whose children sit at nonzero
        // offsets of the clique stack (assembleClique, Factorise.cpp 239-241).
        let g: Int = 40;
        let n = g * g;
        let (mut rows, mut ptr, mut vals) = (Vec::new(), vec![0 as Int], Vec::new());
        for j in 0..n {
            let (r, c) = (j / g, j % g);
            rows.push(j);
            vals.push(4.5);
            if c + 1 < g {
                rows.push(j + 1);
                vals.push(-1.0);
            }
            if r + 1 < g {
                rows.push(j + g);
                vals.push(-1.0);
            }
            ptr.push(rows.len() as Int);
        }
        let mut perm = Vec::new();
        nd_order(0, g, 0, g, g, &mut perm);
        let signs = vec![1 as Int; n as usize];
        let mut fh = FHsolver::new(None, 8);
        let mut s = Symbolic::new();
        assert_eq!(fh.analyse(&mut s, &rows, &ptr, &signs, &perm), kRetOk);
        assert!(s.sn() > 4);
        assert_eq!(fh.factorise(&s, &rows, &ptr, &vals), kRetOk);
        let b: Vec<f64> = (0..n).map(|i| 1.0 + (i % 7) as f64).collect();
        let mut x = b.clone();
        fh.solve(&s, &mut x);
        // residual of the full symmetric matrix
        let mut res = b.clone();
        for j in 0..n as usize {
            for el in ptr[j] as usize..ptr[j + 1] as usize {
                let i = rows[el] as usize;
                res[i] -= vals[el] * x[j];
                if i != j {
                    res[j] -= vals[el] * x[i];
                }
            }
        }
        let err = res.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(err < 1e-10, "residual {err}");
    }
}
