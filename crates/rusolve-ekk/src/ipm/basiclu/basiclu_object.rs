//! HiGHS source: ipm/basiclu/basiclu_object.h, basiclu_object.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Used only by `ipx::basiclu_kernel::BasicLuKernel` (`lu_kernel > 0`, not
//! the IPX default). `malloc`/`free`/`realloc` storage becomes `Vec`s; an
//! empty `istore` or `xstore` is the C NULL "freed object" state (`isvalid`).

use super::basiclu_factorize::basiclu_factorize;
use super::basiclu_get_factors::basiclu_get_factors;
use super::basiclu_solve_dense::basiclu_solve_dense;
use super::basiclu_solve_for_update::basiclu_solve_for_update;
use super::basiclu_solve_sparse::basiclu_solve_sparse;
use super::basiclu_update::basiclu_update;
use super::lu_initialize::lu_initialize;
use super::*;

/// HiGHS: basiclu/basiclu_object.h::struct basiclu_object
#[derive(Debug, Clone, Default)]
pub struct BasicluObject {
    pub istore: Vec<LuInt>,
    pub xstore: Vec<f64>,
    pub Li: Vec<LuInt>,
    pub Ui: Vec<LuInt>,
    pub Wi: Vec<LuInt>,
    pub Lx: Vec<f64>,
    pub Ux: Vec<f64>,
    pub Wx: Vec<f64>,
    pub lhs: Vec<f64>,
    pub ilhs: Vec<LuInt>,
    pub nzlhs: LuInt,
    pub realloc_factor: f64,
}

/// HiGHS: basiclu/basiclu_object.c::lu_reallocix (lines 32-45)
///
/// PORT NOTE: `realloc` failure (`BASICLU_ERROR_out_of_memory`) cannot be
/// observed with `Vec::resize`. A negative `nz` still reports out of memory,
/// as in C, where `nz * sizeof(lu_int)` converts to a huge `size_t`.
fn lu_reallocix(nz: LuInt, p_ai: &mut Vec<LuInt>, p_ax: &mut Vec<f64>) -> LuInt {
    if nz < 0 {
        return BASICLU_ERROR_out_of_memory;
    }
    p_ai.resize(nz as usize, 0);
    p_ax.resize(nz as usize, 0.0);
    BASICLU_OK
}

/// HiGHS: basiclu/basiclu_object.c::lu_realloc_obj (lines 54-89)
fn lu_realloc_obj(obj: &mut BasicluObject) -> LuInt {
    let addmem_l = obj.xstore[BASICLU_ADD_MEMORYL] as LuInt;
    let addmem_u = obj.xstore[BASICLU_ADD_MEMORYU] as LuInt;
    let addmem_w = obj.xstore[BASICLU_ADD_MEMORYW] as LuInt;
    let realloc_factor = f64::max(1.0, obj.realloc_factor);
    let mut status = BASICLU_OK;

    if status == BASICLU_OK && addmem_l > 0 {
        // C: nelem = xstore + addmem (double sum, truncated); nelem *= factor
        let mut nelem = (obj.xstore[BASICLU_MEMORYL] + addmem_l as f64) as LuInt;
        nelem = (nelem as f64 * realloc_factor) as LuInt;
        status = lu_reallocix(nelem, &mut obj.Li, &mut obj.Lx);
        if status == BASICLU_OK {
            obj.xstore[BASICLU_MEMORYL] = nelem as f64;
        }
    }
    if status == BASICLU_OK && addmem_u > 0 {
        let mut nelem = (obj.xstore[BASICLU_MEMORYU] + addmem_u as f64) as LuInt;
        nelem = (nelem as f64 * realloc_factor) as LuInt;
        status = lu_reallocix(nelem, &mut obj.Ui, &mut obj.Ux);
        if status == BASICLU_OK {
            obj.xstore[BASICLU_MEMORYU] = nelem as f64;
        }
    }
    if status == BASICLU_OK && addmem_w > 0 {
        let mut nelem = (obj.xstore[BASICLU_MEMORYW] + addmem_w as f64) as LuInt;
        nelem = (nelem as f64 * realloc_factor) as LuInt;
        status = lu_reallocix(nelem, &mut obj.Wi, &mut obj.Wx);
        if status == BASICLU_OK {
            obj.xstore[BASICLU_MEMORYW] = nelem as f64;
        }
    }
    status
}

/// HiGHS: basiclu/basiclu_object.c::isvalid (lines 94-97)
///
/// C tests the pointers `obj->istore && obj->xstore`; the empty `Vec` is NULL.
fn isvalid(obj: &BasicluObject) -> bool {
    !obj.istore.is_empty() && !obj.xstore.is_empty()
}

/// HiGHS: basiclu/basiclu_object.c::lu_clear_lhs (lines 102-118)
fn lu_clear_lhs(obj: &mut BasicluObject) {
    let m = obj.xstore[BASICLU_DIM] as LuInt;
    let nzsparse = (obj.xstore[BASICLU_SPARSE_THRESHOLD] * m as f64) as LuInt;
    let nz = obj.nzlhs;

    if nz != 0 {
        if nz <= nzsparse {
            for p in 0..nz {
                let k = obj.ilhs[p as usize] as usize;
                obj.lhs[k] = 0.0;
            }
        } else {
            // C: memset(obj->lhs, 0, m*sizeof(double))
            obj.lhs[..m as usize].fill(0.0);
        }
        obj.nzlhs = 0;
    }
}

/// HiGHS: basiclu/basiclu_object.c::basiclu_obj_initialize (lines 123-177)
///
/// PORT NOTE: C's `!obj` check (`BASICLU_ERROR_argument_missing`) has no
/// counterpart for `&mut BasicluObject`. Allocation failure is not
/// observable with `vec!`, so the `out_of_memory` branch is not reproduced.
pub fn basiclu_obj_initialize(obj: &mut BasicluObject, m: LuInt) -> LuInt {
    if m < 0 {
        return BASICLU_ERROR_invalid_argument;
    }

    if m == 0 {
        obj.istore = Vec::new();
        obj.xstore = Vec::new();
        obj.Li = Vec::new();
        obj.Lx = Vec::new();
        obj.Ui = Vec::new();
        obj.Ux = Vec::new();
        obj.Wi = Vec::new();
        obj.Wx = Vec::new();
        obj.lhs = Vec::new();
        obj.ilhs = Vec::new();
        obj.nzlhs = 0;
        return BASICLU_OK;
    }

    let imemsize = (BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * m) as usize;
    let xmemsize = (BASICLU_SIZE_XSTORE_1 + BASICLU_SIZE_XSTORE_M * m) as usize;
    let fmemsize = m; // initial length of Li, Lx, Ui, Ux, Wi, Wx
    let mu = m as usize;

    obj.istore = vec![0; imemsize];
    obj.xstore = vec![0.0; xmemsize];
    obj.Li = vec![0; mu];
    obj.Lx = vec![0.0; mu];
    obj.Ui = vec![0; mu];
    obj.Ux = vec![0.0; mu];
    obj.Wi = vec![0; mu];
    obj.Wx = vec![0.0; mu];
    obj.lhs = vec![0.0; mu]; // calloc
    obj.ilhs = vec![0; mu];
    obj.nzlhs = 0;
    obj.realloc_factor = 1.5;

    lu_initialize(m, &mut obj.istore, &mut obj.xstore);
    obj.xstore[BASICLU_MEMORYL] = fmemsize as f64;
    obj.xstore[BASICLU_MEMORYU] = fmemsize as f64;
    obj.xstore[BASICLU_MEMORYW] = fmemsize as f64;
    BASICLU_OK
}

/// HiGHS: basiclu/basiclu_object.c::basiclu_obj_free (lines 182-198)
pub fn basiclu_obj_free(obj: &mut BasicluObject) {
    obj.istore = Vec::new();
    obj.xstore = Vec::new();
    obj.Li = Vec::new();
    obj.Lx = Vec::new();
    obj.Ui = Vec::new();
    obj.Ux = Vec::new();
    obj.Wi = Vec::new();
    obj.Wx = Vec::new();
    obj.lhs = Vec::new();
    obj.ilhs = Vec::new();
    obj.nzlhs = -1;
}

/// HiGHS: basiclu/basiclu_object.c::basiclu_obj_factorize (lines 203-227)
pub fn basiclu_obj_factorize(
    obj: &mut BasicluObject,
    bbegin: &[LuInt],
    bend: &[LuInt],
    bi: &[LuInt],
    bx: &[f64],
) -> LuInt {
    if !isvalid(obj) {
        return BASICLU_ERROR_invalid_object;
    }

    let mut status = basiclu_factorize(
        &mut obj.istore,
        &mut obj.xstore,
        &mut obj.Li,
        &mut obj.Lx,
        &mut obj.Ui,
        &mut obj.Ux,
        &mut obj.Wi,
        &mut obj.Wx,
        bbegin,
        bend,
        bi,
        bx,
        0,
    );

    while status == BASICLU_REALLOCATE {
        status = lu_realloc_obj(obj);
        if status != BASICLU_OK {
            break;
        }
        status = basiclu_factorize(
            &mut obj.istore,
            &mut obj.xstore,
            &mut obj.Li,
            &mut obj.Lx,
            &mut obj.Ui,
            &mut obj.Ux,
            &mut obj.Wi,
            &mut obj.Wx,
            bbegin,
            bend,
            bi,
            bx,
            1,
        );
    }

    status
}

/// HiGHS: basiclu/basiclu_object.c::basiclu_obj_get_factors (lines 232-245)
#[allow(clippy::too_many_arguments)]
pub fn basiclu_obj_get_factors(
    obj: &mut BasicluObject,
    rowperm: Option<&mut [LuInt]>,
    colperm: Option<&mut [LuInt]>,
    lcolptr: Option<&mut [LuInt]>,
    lrowidx: Option<&mut [LuInt]>,
    lvalue_: Option<&mut [f64]>,
    ucolptr: Option<&mut [LuInt]>,
    urowidx: Option<&mut [LuInt]>,
    uvalue_: Option<&mut [f64]>,
) -> LuInt {
    if !isvalid(obj) {
        return BASICLU_ERROR_invalid_object;
    }

    basiclu_get_factors(
        &mut obj.istore,
        &mut obj.xstore,
        &mut obj.Li,
        &mut obj.Lx,
        &mut obj.Ui,
        &mut obj.Ux,
        &mut obj.Wi,
        &mut obj.Wx,
        rowperm,
        colperm,
        lcolptr,
        lrowidx,
        lvalue_,
        ucolptr,
        urowidx,
        uvalue_,
    )
}

/// HiGHS: basiclu/basiclu_object.c::basiclu_obj_solve_dense (lines 250-259)
pub fn basiclu_obj_solve_dense(
    obj: &mut BasicluObject,
    rhs: &[f64],
    lhs: &mut [f64],
    trans: u8,
) -> LuInt {
    if !isvalid(obj) {
        return BASICLU_ERROR_invalid_object;
    }

    basiclu_solve_dense(
        &mut obj.istore,
        &mut obj.xstore,
        &mut obj.Li,
        &mut obj.Lx,
        &mut obj.Ui,
        &mut obj.Ux,
        &mut obj.Wi,
        &mut obj.Wx,
        rhs,
        lhs,
        trans,
    )
}

/// HiGHS: basiclu/basiclu_object.c::basiclu_obj_solve_sparse (lines 264-275)
pub fn basiclu_obj_solve_sparse(
    obj: &mut BasicluObject,
    nzrhs: LuInt,
    irhs: &[LuInt],
    xrhs: &[f64],
    trans: u8,
) -> LuInt {
    if !isvalid(obj) {
        return BASICLU_ERROR_invalid_object;
    }

    lu_clear_lhs(obj);
    basiclu_solve_sparse(
        &mut obj.istore,
        &mut obj.xstore,
        &mut obj.Li,
        &mut obj.Lx,
        &mut obj.Ui,
        &mut obj.Ux,
        &mut obj.Wi,
        &mut obj.Wx,
        nzrhs,
        irhs,
        xrhs,
        &mut obj.nzlhs,
        &mut obj.ilhs,
        &mut obj.lhs,
        trans,
    )
}

/// HiGHS: basiclu/basiclu_object.c::basiclu_obj_solve_for_update (lines 280-303)
pub fn basiclu_obj_solve_for_update(
    obj: &mut BasicluObject,
    nzrhs: LuInt,
    irhs: &[LuInt],
    xrhs: Option<&[f64]>,
    trans: u8,
    want_solution: LuInt,
) -> LuInt {
    let mut status = BASICLU_OK;

    if !isvalid(obj) {
        return BASICLU_ERROR_invalid_object;
    }

    lu_clear_lhs(obj);
    while status == BASICLU_OK {
        // C: want_solution ? &obj->nzlhs : NULL
        let p_nzlhs = if want_solution != 0 {
            Some(&mut obj.nzlhs)
        } else {
            None
        };
        status = basiclu_solve_for_update(
            &mut obj.istore,
            &mut obj.xstore,
            &mut obj.Li,
            &mut obj.Lx,
            &mut obj.Ui,
            &mut obj.Ux,
            &mut obj.Wi,
            &mut obj.Wx,
            nzrhs,
            irhs,
            xrhs,
            p_nzlhs,
            Some(&mut obj.ilhs[..]),
            Some(&mut obj.lhs[..]),
            trans,
        );
        if status != BASICLU_REALLOCATE {
            break;
        }
        status = lu_realloc_obj(obj);
    }

    status
}

/// HiGHS: basiclu/basiclu_object.c::basiclu_obj_update (lines 308-325)
pub fn basiclu_obj_update(obj: &mut BasicluObject, xtbl: f64) -> LuInt {
    let mut status = BASICLU_OK;

    if !isvalid(obj) {
        return BASICLU_ERROR_invalid_object;
    }

    while status == BASICLU_OK {
        status = basiclu_update(
            &mut obj.istore,
            &mut obj.xstore,
            &mut obj.Li,
            &mut obj.Lx,
            &mut obj.Ui,
            &mut obj.Ux,
            &mut obj.Wi,
            &mut obj.Wx,
            xtbl,
        );
        if status != BASICLU_REALLOCATE {
            break;
        }
        status = lu_realloc_obj(obj);
    }

    status
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn realloc_obj_grows_li_by_factor_and_records_memory() {
        // memoryL = 10, addmemL = 3, factor 1.5: nelem = (10+3)*1.5 = 19.5 -> 19
        let mut obj = BasicluObject {
            xstore: vec![0.0; 80],
            Li: vec![0; 5],
            Lx: vec![0.0; 5],
            Ui: vec![0; 5],
            Ux: vec![0.0; 5],
            Wi: vec![0; 5],
            Wx: vec![0.0; 5],
            realloc_factor: 1.5,
            ..Default::default()
        };
        obj.xstore[BASICLU_MEMORYL] = 10.0;
        obj.xstore[BASICLU_ADD_MEMORYL] = 3.0;
        assert_eq!(lu_realloc_obj(&mut obj), BASICLU_OK);
        assert_eq!(obj.Li.len(), 19);
        assert_eq!(obj.Lx.len(), 19);
        assert_eq!(obj.xstore[BASICLU_MEMORYL], 19.0);
        // no request for U or W: untouched
        assert_eq!(obj.Ui.len(), 5);
        assert_eq!(obj.xstore[BASICLU_MEMORYU], 0.0);
    }

    #[test]
    fn clear_lhs_sparse_branch_zeroes_listed_entries_only() {
        // m = 4, sparse threshold 0.5 -> nzsparse = 2; nz = 2 <= 2: sparse clear
        let mut xstore = vec![0.0; 80];
        xstore[BASICLU_DIM] = 4.0;
        xstore[BASICLU_SPARSE_THRESHOLD] = 0.5;
        let mut obj = BasicluObject {
            xstore,
            lhs: vec![1.0, 2.0, 3.0, 4.0],
            ilhs: vec![0, 2, 1, 3],
            nzlhs: 2,
            ..Default::default()
        };
        lu_clear_lhs(&mut obj);
        assert_eq!(obj.lhs, vec![0.0, 2.0, 0.0, 4.0]);
        assert_eq!(obj.nzlhs, 0);
    }

    #[test]
    fn clear_lhs_dense_branch_memsets_first_m_entries() {
        // nz = 3 > nzsparse = 2: whole first m = 4 entries are cleared
        let mut xstore = vec![0.0; 80];
        xstore[BASICLU_DIM] = 4.0;
        xstore[BASICLU_SPARSE_THRESHOLD] = 0.5;
        let mut obj = BasicluObject {
            xstore,
            lhs: vec![1.0, 2.0, 3.0, 4.0],
            ilhs: vec![1, 0, 3, 2],
            nzlhs: 3,
            ..Default::default()
        };
        lu_clear_lhs(&mut obj);
        assert_eq!(obj.lhs, vec![0.0, 0.0, 0.0, 0.0]);
        assert_eq!(obj.nzlhs, 0);
    }

    #[test]
    fn initialize_rejects_negative_and_empties_for_zero() {
        let mut obj = BasicluObject::default();
        assert_eq!(
            basiclu_obj_initialize(&mut obj, -1),
            BASICLU_ERROR_invalid_argument
        );
        assert_eq!(basiclu_obj_initialize(&mut obj, 0), BASICLU_OK);
        assert!(obj.istore.is_empty());
        assert!(obj.xstore.is_empty());
        assert_eq!(obj.nzlhs, 0);
        // a freed or empty object is invalid
        assert_eq!(
            basiclu_obj_update(&mut obj, 0.0),
            BASICLU_ERROR_invalid_object
        );
        basiclu_obj_free(&mut obj);
        assert_eq!(obj.nzlhs, -1);
    }
}
