//! HiGHS source: simplex/ (subset)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiGHS `simplex/` subset: the EKK simplex solver (HEkk, dual and primal simplex, NLA).

pub mod analysis;
pub mod dual;
pub mod dual_rhs;
pub mod dual_row;
pub mod ekk;
pub mod happ;
pub mod hsimplex;
pub mod hsimplex_debug;
pub mod nla;
pub mod primal;
pub mod simplex_const;
pub mod simplex_struct;
pub mod simplex_timer;
