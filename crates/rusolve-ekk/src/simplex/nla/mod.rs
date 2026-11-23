//! HiGHS source: simplex/HSimplexNla.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Structs for the simplex NLA layer. Bodies: `nla.rs` (HSimplexNla.cpp),
//! `product_form.rs` (HSimplexNlaProductForm.cpp + HSimplexNlaFreeze.cpp), `debug.rs`
//! (HSimplexNlaDebug.cpp, stubs).
//!
//! Ownership (crate docs, pattern (b)): the C++ pointer members `lp_`, `scale_`, `basic_index_`,
//! `options_` are replaced by an [`NlaRefs`] argument built from disjoint `HEkk` fields;
//! `timer_` and `analysis_` are dropped (stubs; `analysis_->analyse_factor_time` is `false`).
//! Translate `lp_->x` as `refs.lp_.x`, `scale_` as `refs.scale_` (an `Option`; C++ `scale_ == NULL`
//! is `refs.scale_.is_none()`), `basic_index_[i]` as `refs.basic_index_[i as usize]`,
//! `options_->x` as `refs.options_.x`. Methods that permute the basic index take
//! `basic_index: &mut [HighsInt]` instead of `refs`.

use crate::prelude::*;

mod debug;
mod nla;
mod product_form;

/// HiGHS: HSimplexNla.h kReportItemLimit
pub const kReportItemLimit: HighsInt = 25;
/// HiGHS: HSimplexNla.h kDensityForIndexing
pub const kDensityForIndexing: f64 = 0.4;

/// HiGHS: HSimplexNla.h::ProductFormUpdate
#[derive(Clone, Debug)]
pub struct ProductFormUpdate {
    pub valid_: bool,
    pub num_row_: HighsInt,
    pub update_count_: HighsInt,
    pub pivot_index_: Vec<HighsInt>,
    pub pivot_value_: Vec<f64>,
    pub start_: Vec<HighsInt>,
    pub index_: Vec<HighsInt>,
    pub value_: Vec<f64>,
}

impl Default for ProductFormUpdate {
    fn default() -> Self {
        Self {
            valid_: false,
            num_row_: Default::default(),
            update_count_: Default::default(),
            pivot_index_: Default::default(),
            pivot_value_: Default::default(),
            start_: Default::default(),
            index_: Default::default(),
            value_: Default::default(),
        }
    }
}

/// HiGHS: HSimplexNla.h::SimplexIterate
#[derive(Clone, Debug)]
pub struct SimplexIterate {
    pub valid_: bool,
    pub basis_: SimplexBasis,
    pub invert_: InvertibleRepresentation,
    pub dual_edge_weight_: Vec<f64>,
}

impl Default for SimplexIterate {
    fn default() -> Self {
        Self {
            valid_: false,
            basis_: Default::default(),
            invert_: Default::default(),
            dual_edge_weight_: Default::default(),
        }
    }
}

/// HiGHS: HSimplexNla.h::HSimplexNla
/// C++ pointer members `lp_`, `scale_`, `basic_index_`, `options_`, `timer_`, `analysis_` are
/// not stored; see [`NlaRefs`].
#[derive(Clone, Debug)]
pub struct HSimplexNla {
    pub factor_: HFactor,
    pub report_: bool,
    pub build_synthetic_tick_: f64,
    pub update_: ProductFormUpdate,
    pub simplex_iterate_: SimplexIterate,
}

impl Default for HSimplexNla {
    fn default() -> Self {
        Self {
            factor_: Default::default(),
            report_: Default::default(),
            build_synthetic_tick_: Default::default(),
            update_: Default::default(),
            simplex_iterate_: Default::default(),
        }
    }
}

/// The borrowed context replacing `HSimplexNla`'s C++ pointer members
/// `lp_`, `scale_`, `basic_index_` and `options_` (crate docs, pattern (b)).
#[derive(Clone, Copy, Debug)]
pub struct NlaRefs<'a> {
    pub lp_: &'a HighsLp,
    /// `Some(&lp.scale_)` iff `lp.scale_.has_scaling && !lp.is_scaled_` (C++ `scale_ != NULL`).
    pub scale_: Option<&'a HighsScale>,
    pub basic_index_: &'a [HighsInt],
    pub options_: &'a HighsOptions,
}

impl<'a> NlaRefs<'a> {
    /// HiGHS: HSimplexNla.cpp::HSimplexNla::setLpAndScalePointers (lines 43-48) plus the
    /// `basic_index_`/`options_` assignments of `setup`/`setPointers`.
    pub fn new(lp: &'a HighsLp, basic_index: &'a [HighsInt], options: &'a HighsOptions) -> Self {
        let scale_ = if lp.scale_.has_scaling && !lp.is_scaled_ {
            Some(&lp.scale_)
        } else {
            None
        };
        NlaRefs {
            lp_: lp,
            scale_,
            basic_index_: basic_index,
            options_: options,
        }
    }
}
