//! HiGHS source: ipm/hipo/ipm/Info.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiPO solve statistics (complete). Fields the C++ leaves uninitialised
//! (`m_solver`, `p_res_rel`, `correctors`, ...) start at zero; every one of
//! them is written before it is read in 1.15.1.

use super::status::{Status, kStatusNotRun};
use crate::ipm::hipo::Int;
use crate::ipm::ipx;

/// HiGHS: Info.h `struct Info`
#[derive(Clone, Copy, Debug)]
pub struct Info {
    pub m_solver: Int,
    pub n_solver: Int,
    pub m_original: Int,
    pub n_original: Int,
    pub status: Status,
    pub p_res_rel: f64,
    pub p_res_abs: f64,
    pub d_res_rel: f64,
    pub d_res_abs: f64,
    pub p_obj: f64,
    pub d_obj: f64,
    pub pd_gap: f64,
    pub ipm_iter: Int,
    pub ipx_used: bool,
    pub ipx_info: ipx::Info,
    pub correctors: Int,
    pub analyse_NE_time: f64,
    pub analyse_AS_time: f64,
    pub matrix_time: f64,
    pub AS_structure_time: f64,
    pub NE_structure_time: f64,
    pub factor_time: f64,
    pub solve_time: f64,
    pub residual_time: f64,
    pub omega_time: f64,
    pub factor_number: Int,
    pub solve_number: Int,
}

impl Default for Info {
    fn default() -> Self {
        Info {
            m_solver: 0,
            n_solver: 0,
            m_original: 0,
            n_original: 0,
            status: kStatusNotRun,
            p_res_rel: 0.0,
            p_res_abs: 0.0,
            d_res_rel: 0.0,
            d_res_abs: 0.0,
            p_obj: 0.0,
            d_obj: 0.0,
            pd_gap: 0.0,
            ipm_iter: 0,
            ipx_used: false,
            ipx_info: ipx::Info::default(),
            correctors: 0,
            analyse_NE_time: 0.0,
            analyse_AS_time: 0.0,
            matrix_time: 0.0,
            AS_structure_time: 0.0,
            NE_structure_time: 0.0,
            factor_time: 0.0,
            solve_time: 0.0,
            residual_time: 0.0,
            omega_time: 0.0,
            factor_number: 0,
            solve_number: 0,
        }
    }
}
