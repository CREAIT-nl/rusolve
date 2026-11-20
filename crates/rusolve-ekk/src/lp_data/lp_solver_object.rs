//! HiGHS source: lp_data/HighsLpSolverObject.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! As in C++, a bundle of mutable references to the caller's objects. The
//! `HEkk& ekk_instance_` member is removed (crate ownership rule (a)): the
//! `HEkk` is the receiver of `HEkk::move_lp`, so it is never inside the bundle.
//! `profiling_` is out of scope and dropped.

use crate::lp_data::hconst::HighsModelStatus;
use crate::lp_data::highs_callback::HighsCallback;
use crate::lp_data::highs_info::HighsInfo;
use crate::lp_data::highs_lp::HighsLp;
use crate::lp_data::highs_options::HighsOptions;
use crate::lp_data::hstruct::{HighsBasis, HighsSolution};
use crate::util::highs_timer::HighsTimer;

/// HiGHS: HighsLpSolverObject.h::HighsLpSolverObject
pub struct HighsLpSolverObject<'a> {
    pub lp_: &'a mut HighsLp,
    pub basis_: &'a mut HighsBasis,
    pub solution_: &'a mut HighsSolution,
    pub highs_info_: &'a mut HighsInfo,
    pub callback_: &'a mut HighsCallback,
    pub options_: &'a mut HighsOptions,
    pub timer_: &'a mut HighsTimer,
    pub model_status_: HighsModelStatus,
}

impl<'a> HighsLpSolverObject<'a> {
    /// HiGHS: HighsLpSolverObject.h::HighsLpSolverObject::HighsLpSolverObject
    /// (`ekk_instance` parameter removed).
    pub fn new(
        lp: &'a mut HighsLp,
        basis: &'a mut HighsBasis,
        solution: &'a mut HighsSolution,
        highs_info: &'a mut HighsInfo,
        callback: &'a mut HighsCallback,
        options: &'a mut HighsOptions,
        timer: &'a mut HighsTimer,
    ) -> Self {
        Self {
            lp_: lp,
            basis_: basis,
            solution_: solution,
            highs_info_: highs_info,
            callback_: callback,
            options_: options,
            timer_: timer,
            model_status_: HighsModelStatus::kNotset,
        }
    }
}
