//! HiGHS source: ipm/ipx/basiclu_kernel.h, basiclu_kernel.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Used only when `lu_kernel > 0` (`ForrestTomlin` over BASICLU); not on
//! the HiGHS default path. `BasicLuHelper` (RAII wrapper) becomes a value
//! holding a `BasicluObject` whose `Drop` is `basiclu_obj_free`.

use super::lu_factorization::LuFactorization;
use super::sparse_matrix::{SparseMatrix, remove_diagonal};
use super::*;
use crate::ipm::basiclu::basiclu_object::{
    BasicluObject, basiclu_obj_factorize, basiclu_obj_free, basiclu_obj_get_factors,
    basiclu_obj_initialize,
};
use crate::ipm::basiclu::*;

/// HiGHS: ipx/basiclu_kernel.cc::BasicLuHelper (lines 11-28)
struct BasicLuHelper {
    obj: BasicluObject,
}

impl BasicLuHelper {
    /// HiGHS: ipx/basiclu_kernel.cc::BasicLuHelper::BasicLuHelper (lines 15-21)
    ///
    /// `throw bad_alloc` / `throw logic_error` -> `panic!`.
    fn new(dim: Int) -> Self {
        let mut obj = BasicluObject::default();
        let err = basiclu_obj_initialize(&mut obj, dim);
        if err == BASICLU_ERROR_out_of_memory {
            panic!("std::bad_alloc");
        }
        if err != BASICLU_OK {
            panic!("basiclu_obj_initialize failed");
        }
        BasicLuHelper { obj }
    }
}

impl Drop for BasicLuHelper {
    /// HiGHS: ipx/basiclu_kernel.cc::BasicLuHelper::~BasicLuHelper (lines 23-25)
    fn drop(&mut self) {
        basiclu_obj_free(&mut self.obj);
    }
}

/// HiGHS: ipx/basiclu_kernel.h::BasicLuKernel
#[derive(Debug, Default)]
pub struct BasicLuKernel {
    /// base member `LuFactorization::stability_`
    stability_: f64,
}

impl LuFactorization for BasicLuKernel {
    fn stability_ref(&self) -> f64 {
        self.stability_
    }
    fn stability_mut(&mut self) -> &mut f64 {
        &mut self.stability_
    }

    /// HiGHS: ipx/basiclu_kernel.cc::BasicLuKernel::_Factorize (lines 30-69)
    fn _factorize(
        &mut self,
        dim: Int,
        Bbegin: &[Int],
        Bend: &[Int],
        Bi: &[Int],
        Bx: &[f64],
        pivottol: f64,
        strict_abs_pivottol: bool,
        L: &mut SparseMatrix,
        U: &mut SparseMatrix,
        rowperm: &mut Vec<Int>,
        colperm: &mut Vec<Int>,
        dependent_cols: &mut Vec<Int>,
    ) {
        let mut lu = BasicLuHelper::new(dim);
        lu.obj.xstore[BASICLU_REL_PIVOT_TOLERANCE] = pivottol;
        if strict_abs_pivottol {
            lu.obj.xstore[BASICLU_ABS_PIVOT_TOLERANCE] = kLuDependencyTol;
            lu.obj.xstore[BASICLU_REMOVE_COLUMNS] = 1.0;
        }
        let err = basiclu_obj_factorize(&mut lu.obj, Bbegin, Bend, Bi, Bx);
        if err == BASICLU_ERROR_out_of_memory {
            panic!("std::bad_alloc");
        }
        if err != BASICLU_OK && err != BASICLU_WARNING_singular_matrix {
            panic!("basiclu_obj_factorize failed");
        }

        // Extract factors. Dependent columns are at the end of the BASICLU pivot
        // sequence.
        let rank = lu.obj.xstore[BASICLU_RANK] as Int;
        dependent_cols.clear();
        for k in rank..dim {
            dependent_cols.push(k);
        }
        // C: `dim + lu.obj.xstore[...]` is double, converted to Int by resize().
        L.resize(dim, dim, (dim as f64 + lu.obj.xstore[BASICLU_LNZ]) as Int);
        U.resize(dim, dim, (dim as f64 + lu.obj.xstore[BASICLU_UNZ]) as Int);
        rowperm.resize(dim as usize, 0);
        colperm.resize(dim as usize, 0);
        let err = {
            let (lcolptr, lrowidx, lvalue) = L.arrays_mut();
            let (ucolptr, urowidx, uvalue) = U.arrays_mut();
            basiclu_obj_get_factors(
                &mut lu.obj,
                Some(rowperm.as_mut_slice()),
                Some(colperm.as_mut_slice()),
                Some(lcolptr),
                Some(lrowidx),
                Some(lvalue),
                Some(ucolptr),
                Some(urowidx),
                Some(uvalue),
            )
        };
        if err != BASICLU_OK {
            panic!("basiclu_obj_get_factors failed");
        }

        // Remove unit diagonal from L.
        let num_dropped = remove_diagonal(L, None);
        debug_assert_eq!(num_dropped, dim);
        debug_assert!(L.entries() as f64 == lu.obj.xstore[BASICLU_LNZ]);
    }
}
