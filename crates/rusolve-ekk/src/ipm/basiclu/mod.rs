//! HiGHS source: ipm/basiclu/basiclu.h, lu_def.h (constants), plus one Rust
//! file per BASICLU `.c` file
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! # Storage model
//!
//! BASICLU keeps its whole state in two caller-owned arrays, `istore`
//! (`lu_int`) and `xstore` (`double`), plus six caller-owned factor arrays
//! `Li/Lx/Ui/Ux/Wi/Wx`. Every API routine starts with `lu_load`, which
//! copies the scalars out of `xstore` into a `struct lu` and points the
//! array members of `struct lu` into `istore`/`xstore`, and ends with
//! `lu_save`, which writes the scalars back.
//!
//! In Rust, [`lu_internal::lu_load`] returns an [`lu_internal::Lu`] that
//! holds the scalars by value and one `&mut [T]` slice per disjoint region
//! of `istore`/`xstore` (carved with `split_at_mut` in the C partition
//! order), plus the six factor arrays. The C "second list" aliases (arrays
//! that share memory with a first-list array, used in another phase) are
//! not fields: they are written as the region they live in, with the C
//! offset (`Lbegin` = `wbegin[m+1..]`, `pivotcol` = `colcount_flink`, ...),
//! see the table on [`lu_internal::Lu`].
//!
//! C pointer arguments become slices that start at the C pointer
//! (`Wbegin+m` -> `&wbegin[m..]`). A nullable C pointer becomes `Option`.
//! The `(void *)` type puns of BASICLU (`pstack` over `work1`, `cancelled`
//! over `row_pivot`) keep their C storage: the integer is stored in an
//! `f64` slot through [`pun_set_int`]/[`pun_get_int`] (bit casts, no
//! rounding), one integer per slot.

pub mod basiclu_factorize;
pub mod basiclu_get_factors;
pub mod basiclu_initialize;
pub mod basiclu_object;
pub mod basiclu_solve_dense;
pub mod basiclu_solve_for_update;
pub mod basiclu_solve_sparse;
pub mod basiclu_update;
pub mod lu_build_factors;
pub mod lu_condest;
pub mod lu_dfs;
pub mod lu_factorize_bump;
pub mod lu_file;
pub mod lu_garbage_perm;
pub mod lu_initialize;
pub mod lu_internal;
pub mod lu_list;
pub mod lu_markowitz;
pub mod lu_matrix_norm;
pub mod lu_pivot;
pub mod lu_residual_test;
pub mod lu_setup_bump;
pub mod lu_singletons;
pub mod lu_solve_dense;
pub mod lu_solve_for_update;
pub mod lu_solve_sparse;
pub mod lu_solve_symbolic;
pub mod lu_solve_triangular;
pub mod lu_update;

/// HiGHS: basiclu/lu_def.h `lu_int` (= `HighsInt` = `int` in the default build)
pub type LuInt = i32;

/// HiGHS: basiclu/lu_def.h::LU_INT_MAX
pub const LU_INT_MAX: LuInt = LuInt::MAX;

// basiclu.h: size of istore and xstore
pub const BASICLU_SIZE_ISTORE_1: LuInt = 1024;
pub const BASICLU_SIZE_ISTORE_M: LuInt = 21;
pub const BASICLU_SIZE_XSTORE_1: LuInt = 1024;
pub const BASICLU_SIZE_XSTORE_M: LuInt = 4;

// basiclu.h: status codes
pub const BASICLU_OK: LuInt = 0;
pub const BASICLU_REALLOCATE: LuInt = 1;
pub const BASICLU_WARNING_singular_matrix: LuInt = 2;
/// rusolve BASICLU time limit (not in the HiGHS reference): basiclu.h::BASICLU_WARNING_timeout
pub const BASICLU_WARNING_timeout: LuInt = 3;
pub const BASICLU_ERROR_invalid_store: LuInt = -1;
pub const BASICLU_ERROR_invalid_call: LuInt = -2;
pub const BASICLU_ERROR_argument_missing: LuInt = -3;
pub const BASICLU_ERROR_invalid_argument: LuInt = -4;
pub const BASICLU_ERROR_maximum_updates: LuInt = -5;
pub const BASICLU_ERROR_singular_update: LuInt = -6;
pub const BASICLU_ERROR_invalid_object: LuInt = -8;
pub const BASICLU_ERROR_out_of_memory: LuInt = -9;

// basiclu.h: user parameters in xstore (index)
pub const BASICLU_MEMORYL: usize = 1;
pub const BASICLU_MEMORYU: usize = 2;
pub const BASICLU_MEMORYW: usize = 3;
pub const BASICLU_DROP_TOLERANCE: usize = 4;
pub const BASICLU_ABS_PIVOT_TOLERANCE: usize = 5;
pub const BASICLU_REL_PIVOT_TOLERANCE: usize = 6;
pub const BASICLU_BIAS_NONZEROS: usize = 7;
pub const BASICLU_MAXN_SEARCH_PIVOT: usize = 8;
pub const BASICLU_PAD: usize = 9;
pub const BASICLU_STRETCH: usize = 10;
pub const BASICLU_COMPRESSION_THRESHOLD: usize = 11;
pub const BASICLU_SPARSE_THRESHOLD: usize = 12;
pub const BASICLU_REMOVE_COLUMNS: usize = 13;
// rusolve BASICLU time limit (not in the HiGHS reference): basiclu.h::BASICLU_TIME_START / BASICLU_TIME_LIMIT
pub const BASICLU_TIME_START: usize = 15;
pub const BASICLU_TIME_LIMIT: usize = 16;
pub const BASICLU_SEARCH_ROWS: usize = 14;

// basiclu.h: user readable info in xstore (index)
pub const BASICLU_DIM: usize = 64;
pub const BASICLU_STATUS: usize = 65;
pub const BASICLU_ADD_MEMORYL: usize = 66;
pub const BASICLU_ADD_MEMORYU: usize = 67;
pub const BASICLU_ADD_MEMORYW: usize = 68;
pub const BASICLU_NUPDATE: usize = 70;
pub const BASICLU_NFORREST: usize = 71;
pub const BASICLU_NFACTORIZE: usize = 72;
pub const BASICLU_NUPDATE_TOTAL: usize = 73;
pub const BASICLU_NFORREST_TOTAL: usize = 74;
pub const BASICLU_NSYMPERM_TOTAL: usize = 75;
pub const BASICLU_LNZ: usize = 76;
pub const BASICLU_UNZ: usize = 77;
pub const BASICLU_RNZ: usize = 78;
pub const BASICLU_MIN_PIVOT: usize = 79;
pub const BASICLU_MAX_PIVOT: usize = 80;
pub const BASICLU_UPDATE_COST: usize = 81;
pub const BASICLU_TIME_FACTORIZE: usize = 82;
pub const BASICLU_TIME_SOLVE: usize = 83;
pub const BASICLU_TIME_UPDATE: usize = 84;
pub const BASICLU_TIME_FACTORIZE_TOTAL: usize = 85;
pub const BASICLU_TIME_SOLVE_TOTAL: usize = 86;
pub const BASICLU_TIME_UPDATE_TOTAL: usize = 87;
pub const BASICLU_LFLOPS: usize = 88;
pub const BASICLU_UFLOPS: usize = 89;
pub const BASICLU_RFLOPS: usize = 90;
pub const BASICLU_CONDEST_L: usize = 91;
pub const BASICLU_CONDEST_U: usize = 92;
pub const BASICLU_MAX_ETA: usize = 93;
pub const BASICLU_NORM_L: usize = 94;
pub const BASICLU_NORM_U: usize = 95;
pub const BASICLU_NORMEST_LINV: usize = 96;
pub const BASICLU_NORMEST_UINV: usize = 97;
pub const BASICLU_MATRIX_ONENORM: usize = 98;
pub const BASICLU_MATRIX_INFNORM: usize = 99;
pub const BASICLU_RESIDUAL_TEST: usize = 111;
pub const BASICLU_MATRIX_NZ: usize = 100;
pub const BASICLU_RANK: usize = 101;
pub const BASICLU_BUMP_SIZE: usize = 102;
pub const BASICLU_BUMP_NZ: usize = 103;
pub const BASICLU_NSEARCH_PIVOT: usize = 104;
pub const BASICLU_NEXPAND: usize = 105;
pub const BASICLU_NGARBAGE: usize = 106;
pub const BASICLU_FACTOR_FLOPS: usize = 107;
pub const BASICLU_TIME_SINGLETONS: usize = 108;
pub const BASICLU_TIME_SEARCH_PIVOT: usize = 109;
pub const BASICLU_TIME_ELIM_PIVOT: usize = 110;
pub const BASICLU_PIVOT_ERROR: usize = 120;

/// HiGHS: basiclu/lu_def.h::BASICLU_HASH (hash in istore[0], xstore[0])
pub const BASICLU_HASH: LuInt = 7743090;

// lu_def.h: `enum { NO_TASK, SINGLETONS, SETUP_BUMP, FACTORIZE_BUMP, BUILD_FACTORS }`
pub const NO_TASK: LuInt = 0;
pub const SINGLETONS: LuInt = 1;
pub const SETUP_BUMP: LuInt = 2;
pub const FACTORIZE_BUMP: LuInt = 3;
pub const BUILD_FACTORS: LuInt = 4;

/// HiGHS: basiclu/lu_def.h::lu_iswap (static inline)
#[inline]
pub(crate) fn lu_iswap(x: &mut [LuInt], i: usize, j: usize) {
    x.swap(i, j);
}

/// HiGHS: basiclu/lu_def.h::lu_fswap (static inline)
#[inline]
pub(crate) fn lu_fswap(x: &mut [f64], i: usize, j: usize) {
    x.swap(i, j);
}

/// Store an integer in an `f64` workspace slot (C `(void *)` pun, e.g.
/// `lu_int *pstack = (void *) this->work1`). Exact bit cast, so any value
/// round-trips through [`pun_get_int`]. One integer per slot.
#[inline]
pub fn pun_set_int(work: &mut [f64], k: usize, value: i64) {
    work[k] = f64::from_bits(value as u64);
}

/// Read an integer stored by [`pun_set_int`].
#[inline]
pub fn pun_get_int(work: &[f64], k: usize) -> i64 {
    work[k].to_bits() as i64
}
