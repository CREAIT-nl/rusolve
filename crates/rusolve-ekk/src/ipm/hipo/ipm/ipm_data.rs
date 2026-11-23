//! HiGHS source: ipm/hipo/ipm/IpmData.h, ipm/hipo/ipm/IpmData.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Per-iteration IPM statistics (complete).

use crate::ipm::hipo::Int;

/// HiGHS: IpmData.h `struct IpmIterData`
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IpmIterData {
    pub sigma_aff: f64,
    pub sigma: f64,
    pub correctors: Int,
    pub num_solves: Int,
    pub min_theta: f64,
    pub max_theta: f64,
    pub min_prod: f64,
    pub max_prod: f64,
    pub num_small_prod: Int,
    pub num_large_prod: Int,
    pub omega: f64,
    pub nw_back_err: f64,
    pub cw_back_err: f64,
    pub large_components_cw: Int,
}

/// HiGHS: IpmData.h `struct IpmData`
#[derive(Clone, Debug, Default)]
pub struct IpmData {
    pub record: Vec<IpmIterData>,
}

impl IpmData {
    /// HiGHS: IpmData.cpp::IpmData::append (lines 5-5)
    pub fn append(&mut self) {
        self.record.push(IpmIterData::default());
    }

    /// HiGHS: IpmData.cpp::IpmData::back (lines 6-6)
    pub fn back(&mut self) -> &mut IpmIterData {
        self.record
            .last_mut()
            .expect("IpmData::back on empty record")
    }
}
