//! HiGHS source: model/ (HighsHessian, HighsHessianUtils, HighsModel)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `Highs` keeps the LP in
//! `model_lp_` and the Hessian in `model_hessian_` (C++ `model_.lp_` /
//! `model_.hessian_`); a `HighsModel` is only built where C++ passes one by
//! value (`Highs::passModel(HighsModel)`, `HighsIis::model_`).

pub mod highs_hessian;
pub mod highs_hessian_utils;
pub mod highs_model;
