//! HiGHS source: extern/amd/amd.h, extern/amd/amd_internal.h
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! SuiteSparse AMD 3.3.4 (approximate minimum degree), trimmed to the call
//! graph of `Highs_amd_defaults` + `Highs_amd_order`, the two entry points
//! HiPO uses through `HighsExtras::amd`. `amd_int` is the HiGHS integer
//! (`int`, or `int64_t` under `HIGHSINT64`), `amd_uint` its unsigned twin.
//! `Control` / `Info` are `double` arrays of `AMD_CONTROL` / `AMD_INFO`
//! entries; the C `NULL` is `None`.
//!
//! Not ported: `Highs_amd_control` / `Highs_amd_info` (printers, never called
//! from HiGHS), `Highs_SuiteSparse_malloc` / `_free` (replaced by `Vec`; the
//! size-overflow check of `SuiteSparse_malloc` maps to `AMD_OUT_OF_MEMORY` only
//! for sizes Rust cannot allocate either), the `NDEBUG`-off debug checks
//! (`AMD_debug*`).

pub mod aat;
pub mod amd_2;
pub mod order;
pub mod postorder;

pub use aat::{Highs_amd_1, Highs_amd_aat};
pub use amd_2::Highs_amd_2;
pub use order::{Highs_amd_defaults, Highs_amd_order, Highs_amd_preprocess, Highs_amd_valid};
pub use postorder::{Highs_amd_post_tree, Highs_amd_postorder};

/// HiGHS: amd.h `amd_int`
pub type amd_int = crate::HighsInt;
/// HiGHS: amd.h `amd_uint`
pub type amd_uint = crate::HighsUInt;
/// HiGHS: amd.h `amd_int_max`
pub const amd_int_max: amd_int = amd_int::MAX;

// HiGHS: amd.h Control / Info layout
pub const AMD_CONTROL: usize = 5;
pub const AMD_INFO: usize = 20;
pub const AMD_DENSE: usize = 0;
pub const AMD_AGGRESSIVE: usize = 1;
pub const AMD_DEFAULT_DENSE: f64 = 10.0;
pub const AMD_DEFAULT_AGGRESSIVE: f64 = 1.0;
pub const AMD_STATUS: usize = 0;
pub const AMD_N: usize = 1;
pub const AMD_NZ: usize = 2;
pub const AMD_SYMMETRY: usize = 3;
pub const AMD_NZDIAG: usize = 4;
pub const AMD_NZ_A_PLUS_AT: usize = 5;
pub const AMD_NDENSE: usize = 6;
pub const AMD_MEMORY: usize = 7;
pub const AMD_NCMPA: usize = 8;
pub const AMD_LNZ: usize = 9;
pub const AMD_NDIV: usize = 10;
pub const AMD_NMULTSUBS_LDL: usize = 11;
pub const AMD_NMULTSUBS_LU: usize = 12;
pub const AMD_DMAX: usize = 13;

// HiGHS: amd.h return codes
pub const AMD_OK: i32 = 0;
pub const AMD_OUT_OF_MEMORY: i32 = -1;
pub const AMD_INVALID: i32 = -2;
pub const AMD_OK_BUT_JUMBLED: i32 = 1;

/// HiGHS: amd_internal.h `EMPTY`
pub const EMPTY: amd_int = -1;

/// HiGHS: amd_internal.h `FLIP(i)`
#[inline]
pub fn FLIP(i: amd_int) -> amd_int {
    -i - 2
}

/// HiGHS: amd_internal.h `UNFLIP(i)`
#[inline]
pub fn UNFLIP(i: amd_int) -> amd_int {
    if i < EMPTY { FLIP(i) } else { i }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn flip_is_involution() {
        for i in [-5, -2, -1, 0, 3] {
            assert_eq!(FLIP(FLIP(i)), i);
        }
        assert_eq!(UNFLIP(FLIP(4)), 4);
        assert_eq!(UNFLIP(EMPTY), EMPTY);
    }
}
