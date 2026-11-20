//! HiGHS source: extern/HighsExtrasExternalDeps.h (amd / metis / rcm feature wrappers), extern/{metis,amd,rcm}
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Fill-reducing orderings for FactorHiGHS. HiGHS 1.15.1 loads them through
//! the `HighsExtras` feature API (`HighsExtras::amd::order`,
//! `HighsExtras::metis::nodeND`, `HighsExtras::rcm::genrcm`); when the extras
//! library is absent HiPO is unavailable, which is what
//! `HighsOptions::hipo_available` models. [`extras`] is that facade (complete:
//! plain forwarding), [`metis`], [`amd`] and [`rcm`] are the vendored
//! libraries ported with their own licenses (THIRD_PARTY_NOTICES.md).

pub mod amd;
pub mod metis;
pub mod rcm;

/// HiGHS: HighsExtrasExternalDeps.h `namespace HighsExtras` (the ordering
/// features; BLAS is [`crate::ipm::hipo::blas`]).
pub mod extras {
    /// HiGHS: HighsExtrasExternalDeps.h `struct amd` (lines 100-111)
    pub mod amd {
        use super::super::amd::{Highs_amd_defaults, Highs_amd_order, amd_int};

        pub fn set_defaults(Control: &mut [f64]) {
            Highs_amd_defaults(Control)
        }

        pub fn order(
            n: amd_int,
            Ap: &[amd_int],
            Ai: &[amd_int],
            P: &mut [amd_int],
            Control: Option<&[f64]>,
            Info: Option<&mut [f64]>,
        ) -> i32 {
            Highs_amd_order(n, Ap, Ai, P, Control, Info)
        }
    }

    /// HiGHS: HighsExtrasExternalDeps.h `struct metis` (lines 192-204)
    pub mod metis {
        use super::super::metis::idx_t;
        use super::super::metis::ometis::Highs_METIS_NodeND;
        use super::super::metis::options::Highs_METIS_SetDefaultOptions;

        pub fn set_default_options(options: &mut [idx_t]) -> i32 {
            Highs_METIS_SetDefaultOptions(options)
        }

        pub fn nodeND(
            nvtxs: &mut idx_t,
            xadj: &[idx_t],
            adjncy: &[idx_t],
            vwgt: Option<&[idx_t]>,
            options: Option<&[idx_t]>,
            perm: &mut [idx_t],
            iperm: &mut [idx_t],
        ) -> i32 {
            Highs_METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm)
        }
    }

    /// HiGHS: HighsExtrasExternalDeps.h `struct rcm` (lines 206-214)
    pub mod rcm {
        use crate::HighsInt;

        pub fn genrcm(
            node_num: HighsInt,
            adj_num: HighsInt,
            adj_row: &[HighsInt],
            adj: &[HighsInt],
            perm: &mut [HighsInt],
        ) -> HighsInt {
            super::super::rcm::Highs_genrcm(node_num, adj_num, adj_row, adj, perm)
        }
    }
}
