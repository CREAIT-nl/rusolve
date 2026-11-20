//! HiGHS source: util/HighsLinearSumBounds.h, util/HighsLinearSumBounds.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Activity bounds of linear sums (HPresolve: implied row bounds, implied
//! dual row bounds).
//!
//! # Design: the bound arrays
//!
//! C++ stores six raw pointers (`varLower`, `varUpper`, `implVarLower`,
//! `implVarUpper`, `implVarLowerSource`, `implVarUpperSource`) set once by
//! `setBoundArrays`, pointing into HPresolve vectors. Rust cannot keep those
//! borrows alive, so every method that dereferences them takes a
//! [`HighsLinearSumBoundArrays`] argument `b` instead, and
//! [`HighsLinearSumBounds::set_bound_arrays`] is a documented no-op.
//! Translation: `varLower[var]` -> `b.var_lower[var as usize]`, and so on.
//! HPresolve builds `b` with the `col_bound_arrays!(self)` /
//! `row_dual_bound_arrays!(self)` macros (see `presolve/hpresolve/mod.rs`).
//!
//! The private C++ overloads get distinct Rust names:
//! * `getImplVarUpper(sum, myVarUpper, myImplVarUpper, myImplVarUpperSource)`
//!   -> [`HighsLinearSumBounds::get_impl_var_upper_from`] (same for lower);
//! * private `updatedImplVarUpper(sum, var, coef, oldVarUpper,
//!   oldImplVarUpper, oldImplVarUpperSource)` ->
//!   [`HighsLinearSumBounds::updated_impl_var_upper_full`] (same for lower);
//! * `update(numInfs, activity, direction, bound, coef)` ->
//!   [`HighsLinearSumBounds::update_dir`], `update(numInfs, activity,
//!   oldBound, newBound, coef)` -> [`HighsLinearSumBounds::update_old_new`].
//!
//! `update`/`residual` are `const` members that touch no member, so they are
//! associated functions. `handleVar*` take `&mut self` and pass
//! `&mut self.num_inf_sum_upper_orig[sum]` / `&mut self.sum_upper_orig[sum]`
//! (disjoint fields) to `update_dir`; keep the `coefficient > 0 ? a : b`
//! selection as an `if` choosing which pair of fields to borrow.
//!
//! `activity += direction * static_cast<HighsCDouble>(bound) * coefficient`
//! is `(double(direction) * HighsCDouble(bound)) * coefficient`:
//! `*activity += (direction as f64) * HighsCDouble::from(bound) * coefficient`
//! (`f64 * HighsCDouble`, then `HighsCDouble * f64`, left to right).

use crate::prelude::*;

/// The six C++ bound-array pointers of [`HighsLinearSumBounds`], borrowed
/// for the duration of one call.
#[derive(Clone, Copy)]
pub struct HighsLinearSumBoundArrays<'a> {
    pub var_lower: &'a [f64],
    pub var_upper: &'a [f64],
    pub impl_var_lower: &'a [f64],
    pub impl_var_upper: &'a [f64],
    pub impl_var_lower_source: &'a [HighsInt],
    pub impl_var_upper_source: &'a [HighsInt],
}

/// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds
///
/// The C++ pointer members are not stored (see the module doc).
#[derive(Clone, Debug, Default)]
pub struct HighsLinearSumBounds {
    pub sum_lower_orig: Vec<HighsCDouble>,
    pub sum_upper_orig: Vec<HighsCDouble>,
    pub num_inf_sum_lower_orig: Vec<HighsInt>,
    pub num_inf_sum_upper_orig: Vec<HighsInt>,
    pub sum_lower: Vec<HighsCDouble>,
    pub sum_upper: Vec<HighsCDouble>,
    pub num_inf_sum_lower: Vec<HighsInt>,
    pub num_inf_sum_upper: Vec<HighsInt>,
}

impl HighsLinearSumBounds {
    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::setNumSums
    pub fn set_num_sums(&mut self, num_sums: HighsInt) {
        let n = num_sums as usize;
        self.num_inf_sum_lower.resize(n, 0);
        self.num_inf_sum_upper.resize(n, 0);
        self.sum_lower.resize(n, HighsCDouble::default());
        self.sum_upper.resize(n, HighsCDouble::default());
        self.num_inf_sum_lower_orig.resize(n, 0);
        self.num_inf_sum_upper_orig.resize(n, 0);
        self.sum_lower_orig.resize(n, HighsCDouble::default());
        self.sum_upper_orig.resize(n, HighsCDouble::default());
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::setBoundArrays
    ///
    /// No-op: the arrays are passed to each call as
    /// [`HighsLinearSumBoundArrays`]. Keep the call so the structure stays
    /// parallel with the C++.
    pub fn set_bound_arrays(&mut self) {}

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::sumScaled
    pub fn sum_scaled(&mut self, sum: HighsInt, scale: f64) {
        let s = sum as usize;
        self.sum_lower_orig[s] *= scale;
        self.sum_upper_orig[s] *= scale;
        self.sum_lower[s] *= scale;
        self.sum_upper[s] *= scale;

        if scale < 0.0 {
            std::mem::swap(&mut self.sum_lower[s], &mut self.sum_upper[s]);
            std::mem::swap(&mut self.sum_lower_orig[s], &mut self.sum_upper_orig[s]);
            std::mem::swap(
                &mut self.num_inf_sum_lower[s],
                &mut self.num_inf_sum_upper[s],
            );
            std::mem::swap(
                &mut self.num_inf_sum_lower_orig[s],
                &mut self.num_inf_sum_upper_orig[s],
            );
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::add (lines 12-17, 6 lines)
    pub fn add(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
    ) {
        let v = var as usize;
        self.handle_var_upper(sum, coefficient, b.var_upper[v], 1);
        self.handle_var_lower(sum, coefficient, b.var_lower[v], 1);
        let impl_var_upper = self.get_impl_var_upper(b, sum, var);
        self.handle_impl_var_upper(sum, coefficient, impl_var_upper, 1);
        let impl_var_lower = self.get_impl_var_lower(b, sum, var);
        self.handle_impl_var_lower(sum, coefficient, impl_var_lower, 1);
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::remove (lines 19-25, 7 lines)
    pub fn remove(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
    ) {
        let v = var as usize;
        self.handle_var_upper(sum, coefficient, b.var_upper[v], -1);
        self.handle_var_lower(sum, coefficient, b.var_lower[v], -1);
        let impl_var_upper = self.get_impl_var_upper(b, sum, var);
        self.handle_impl_var_upper(sum, coefficient, impl_var_upper, -1);
        let impl_var_lower = self.get_impl_var_lower(b, sum, var);
        self.handle_impl_var_lower(sum, coefficient, impl_var_lower, -1);
    }

    /// HiGHS 1.16.0 (#3318): HighsLinearSumBounds.cpp::HighsLinearSumBounds::implBndsValidAfterCoefChange (lines 176-191)
    ///
    /// Replaces the contribution of `var` to `sum` (`old_coef` by `new_coef`, 0.0 meaning
    /// absent) and reports whether neither activity bound was relaxed. Only called under
    /// `mip_coef_change_impl_bounds` (HiGHS #3318).
    pub fn impl_bnds_valid_after_coef_change(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        old_coef: f64,
        new_coef: f64,
    ) -> bool {
        let s = sum as usize;
        let old_num_inf_lower = self.num_inf_sum_lower[s];
        let old_sum_lower = self.sum_lower[s];
        let old_num_inf_upper = self.num_inf_sum_upper[s];
        let old_sum_upper = self.sum_upper[s];

        if old_coef != 0.0 {
            self.remove(b, sum, var, old_coef);
        }
        if new_coef != 0.0 {
            self.add(b, sum, var, new_coef);
        }

        self.num_inf_sum_lower[s] <= old_num_inf_lower
            && self.sum_lower[s] >= old_sum_lower
            && self.num_inf_sum_upper[s] <= old_num_inf_upper
            && self.sum_upper[s] <= old_sum_upper
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::updatedVarUpper (lines 27-34, 8 lines)
    pub fn updated_var_upper(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        old_var_upper: f64,
    ) {
        let v = var as usize;
        self.handle_var_upper(sum, coefficient, old_var_upper, -1);
        self.handle_var_upper(sum, coefficient, b.var_upper[v], 1);
        self.updated_impl_var_upper_full(
            b,
            sum,
            var,
            coefficient,
            old_var_upper,
            b.impl_var_upper[v],
            b.impl_var_upper_source[v],
        );
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::updatedVarLower (lines 36-43, 8 lines)
    pub fn updated_var_lower(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        old_var_lower: f64,
    ) {
        let v = var as usize;
        self.handle_var_lower(sum, coefficient, old_var_lower, -1);
        self.handle_var_lower(sum, coefficient, b.var_lower[v], 1);
        self.updated_impl_var_lower_full(
            b,
            sum,
            var,
            coefficient,
            old_var_lower,
            b.impl_var_lower[v],
            b.impl_var_lower_source[v],
        );
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::updatedImplVarUpper (lines 45-51, 7 lines)
    ///
    /// Public 5-argument overload; forwards to
    /// [`Self::updated_impl_var_upper_full`] with `varUpper[var]`.
    pub fn updated_impl_var_upper(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        old_impl_var_upper: f64,
        old_impl_var_upper_source: HighsInt,
    ) {
        let v = var as usize;
        self.updated_impl_var_upper_full(
            b,
            sum,
            var,
            coefficient,
            b.var_upper[v],
            old_impl_var_upper,
            old_impl_var_upper_source,
        );
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::updatedImplVarLower (lines 53-59, 7 lines)
    ///
    /// Public 5-argument overload; forwards to
    /// [`Self::updated_impl_var_lower_full`] with `varLower[var]`.
    pub fn updated_impl_var_lower(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        old_impl_var_lower: f64,
        old_impl_var_lower_source: HighsInt,
    ) {
        let v = var as usize;
        self.updated_impl_var_lower_full(
            b,
            sum,
            var,
            coefficient,
            b.var_lower[v],
            old_impl_var_lower,
            old_impl_var_lower_source,
        );
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::getResidualSumLower (lines 61-78, 18 lines)
    ///
    /// C++ default arguments: `boundVar = -1, boundVarCoefficient = -kHighsInf,
    /// boundVarValue = -kHighsInf`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_residual_sum_lower(
        &self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        bound_var: HighsInt,
        bound_var_coefficient: f64,
        bound_var_value: f64,
    ) -> f64 {
        let s = sum as usize;
        let mut activity = self.sum_lower[s];
        let mut num_infs = self.num_inf_sum_lower[s];
        let var_bound = if coefficient > 0.0 {
            self.get_impl_var_lower(b, sum, var)
        } else {
            self.get_impl_var_upper(b, sum, var)
        };
        let old_bound_var_bound = if bound_var != -1 {
            if bound_var_coefficient > 0.0 {
                self.get_impl_var_lower(b, sum, bound_var)
            } else {
                self.get_impl_var_upper(b, sum, bound_var)
            }
        } else {
            kHighsInf
        };
        Self::residual(
            &mut num_infs,
            &mut activity,
            var_bound,
            coefficient,
            bound_var,
            bound_var_coefficient,
            old_bound_var_bound,
            bound_var_value,
        );
        if num_infs == 0 {
            f64::from(activity)
        } else {
            -kHighsInf
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::getResidualSumUpper (lines 80-97, 18 lines)
    ///
    /// C++ default arguments: `boundVar = -1, boundVarCoefficient = -kHighsInf,
    /// boundVarValue = -kHighsInf`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_residual_sum_upper(
        &self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        bound_var: HighsInt,
        bound_var_coefficient: f64,
        bound_var_value: f64,
    ) -> f64 {
        let s = sum as usize;
        let mut activity = self.sum_upper[s];
        let mut num_infs = self.num_inf_sum_upper[s];
        let var_bound = if coefficient < 0.0 {
            self.get_impl_var_lower(b, sum, var)
        } else {
            self.get_impl_var_upper(b, sum, var)
        };
        let old_bound_var_bound = if bound_var != -1 {
            if bound_var_coefficient < 0.0 {
                self.get_impl_var_lower(b, sum, bound_var)
            } else {
                self.get_impl_var_upper(b, sum, bound_var)
            }
        } else {
            kHighsInf
        };
        Self::residual(
            &mut num_infs,
            &mut activity,
            var_bound,
            coefficient,
            bound_var,
            bound_var_coefficient,
            old_bound_var_bound,
            bound_var_value,
        );
        if num_infs == 0 {
            f64::from(activity)
        } else {
            kHighsInf
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::getResidualSumLowerOrig (lines 99-111, 13 lines)
    ///
    /// C++ default arguments: `boundVar = -1, boundVarCoefficient = -kHighsInf,
    /// boundVarValue = -kHighsInf`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_residual_sum_lower_orig(
        &self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        bound_var: HighsInt,
        bound_var_coefficient: f64,
        bound_var_value: f64,
    ) -> f64 {
        let s = sum as usize;
        let mut activity = self.sum_lower_orig[s];
        let mut num_infs = self.num_inf_sum_lower_orig[s];
        let v = var as usize;
        let var_bound = if coefficient > 0.0 {
            b.var_lower[v]
        } else {
            b.var_upper[v]
        };
        let old_bound_var_bound = if bound_var != -1 {
            let bv = bound_var as usize;
            if bound_var_coefficient > 0.0 {
                b.var_lower[bv]
            } else {
                b.var_upper[bv]
            }
        } else {
            kHighsInf
        };
        Self::residual(
            &mut num_infs,
            &mut activity,
            var_bound,
            coefficient,
            bound_var,
            bound_var_coefficient,
            old_bound_var_bound,
            bound_var_value,
        );
        if num_infs == 0 {
            f64::from(activity)
        } else {
            -kHighsInf
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::getResidualSumUpperOrig (lines 113-125, 13 lines)
    ///
    /// C++ default arguments: `boundVar = -1, boundVarCoefficient = -kHighsInf,
    /// boundVarValue = -kHighsInf`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_residual_sum_upper_orig(
        &self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        bound_var: HighsInt,
        bound_var_coefficient: f64,
        bound_var_value: f64,
    ) -> f64 {
        let s = sum as usize;
        let mut activity = self.sum_upper_orig[s];
        let mut num_infs = self.num_inf_sum_upper_orig[s];
        let v = var as usize;
        let var_bound = if coefficient < 0.0 {
            b.var_lower[v]
        } else {
            b.var_upper[v]
        };
        let old_bound_var_bound = if bound_var != -1 {
            let bv = bound_var as usize;
            if bound_var_coefficient < 0.0 {
                b.var_lower[bv]
            } else {
                b.var_upper[bv]
            }
        } else {
            kHighsInf
        };
        Self::residual(
            &mut num_infs,
            &mut activity,
            var_bound,
            coefficient,
            bound_var,
            bound_var_coefficient,
            old_bound_var_bound,
            bound_var_value,
        );
        if num_infs == 0 {
            f64::from(activity)
        } else {
            kHighsInf
        }
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::getSumLowerOrig
    ///
    /// C++ `template <typename T = double> (sum, T offset = T())`: pass
    /// `0.0` for the default; a `double` or `HighsCDouble` offset is converted
    /// with `HighsCDouble::from` exactly like `static_cast<HighsCDouble>`.
    pub fn get_sum_lower_orig<T: Into<HighsCDouble>>(&self, sum: HighsInt, offset: T) -> f64 {
        if self.num_inf_sum_lower_orig[sum as usize] == 0 {
            f64::from(self.sum_lower_orig[sum as usize] + offset.into())
        } else {
            -kHighsInf
        }
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::getSumUpperOrig
    ///
    /// C++ default arguments: `offset = T()` (pass `0.0`).
    pub fn get_sum_upper_orig<T: Into<HighsCDouble>>(&self, sum: HighsInt, offset: T) -> f64 {
        if self.num_inf_sum_upper_orig[sum as usize] == 0 {
            f64::from(self.sum_upper_orig[sum as usize] + offset.into())
        } else {
            kHighsInf
        }
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::getSumLower
    ///
    /// C++ default arguments: `offset = T()` (pass `0.0`).
    pub fn get_sum_lower<T: Into<HighsCDouble>>(&self, sum: HighsInt, offset: T) -> f64 {
        if self.num_inf_sum_lower[sum as usize] == 0 {
            f64::from(self.sum_lower[sum as usize] + offset.into())
        } else {
            -kHighsInf
        }
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::getSumUpper
    ///
    /// C++ default arguments: `offset = T()` (pass `0.0`).
    pub fn get_sum_upper<T: Into<HighsCDouble>>(&self, sum: HighsInt, offset: T) -> f64 {
        if self.num_inf_sum_upper[sum as usize] == 0 {
            f64::from(self.sum_upper[sum as usize] + offset.into())
        } else {
            kHighsInf
        }
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::getNumInfSumLower
    pub fn get_num_inf_sum_lower(&self, sum: HighsInt) -> HighsInt {
        self.num_inf_sum_lower[sum as usize]
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::getNumInfSumUpper
    pub fn get_num_inf_sum_upper(&self, sum: HighsInt) -> HighsInt {
        self.num_inf_sum_upper[sum as usize]
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::getNumInfSumLowerOrig
    pub fn get_num_inf_sum_lower_orig(&self, sum: HighsInt) -> HighsInt {
        self.num_inf_sum_lower_orig[sum as usize]
    }

    /// HiGHS: HighsLinearSumBounds.h::HighsLinearSumBounds::getNumInfSumUpperOrig
    pub fn get_num_inf_sum_upper_orig(&self, sum: HighsInt) -> HighsInt {
        self.num_inf_sum_upper_orig[sum as usize]
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::shrink (lines 127-150, 24 lines)
    ///
    /// `resize(newSize)` on the vectors: `resize(new_size as usize, Default)`.
    pub fn shrink(&mut self, new_indices: &[HighsInt], new_size: HighsInt) {
        for (i, &new_index) in new_indices.iter().enumerate() {
            if new_index != -1 {
                let t = new_index as usize;
                self.sum_lower[t] = self.sum_lower[i];
                self.sum_upper[t] = self.sum_upper[i];
                self.num_inf_sum_lower[t] = self.num_inf_sum_lower[i];
                self.num_inf_sum_upper[t] = self.num_inf_sum_upper[i];
                self.sum_lower_orig[t] = self.sum_lower_orig[i];
                self.sum_upper_orig[t] = self.sum_upper_orig[i];
                self.num_inf_sum_lower_orig[t] = self.num_inf_sum_lower_orig[i];
                self.num_inf_sum_upper_orig[t] = self.num_inf_sum_upper_orig[i];
            }
        }

        let n = new_size as usize;
        self.sum_lower.resize(n, HighsCDouble::default());
        self.sum_upper.resize(n, HighsCDouble::default());
        self.num_inf_sum_lower.resize(n, 0);
        self.num_inf_sum_upper.resize(n, 0);
        self.sum_lower_orig.resize(n, HighsCDouble::default());
        self.sum_upper_orig.resize(n, HighsCDouble::default());
        self.num_inf_sum_lower_orig.resize(n, 0);
        self.num_inf_sum_upper_orig.resize(n, 0);
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::getImplVarUpper (lines 152-155, 4 lines)
    pub fn get_impl_var_upper(
        &self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
    ) -> f64 {
        let v = var as usize;
        Self::get_impl_var_upper_from(
            sum,
            b.var_upper[v],
            b.impl_var_upper[v],
            b.impl_var_upper_source[v],
        )
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::getImplVarLower (lines 157-160, 4 lines)
    pub fn get_impl_var_lower(
        &self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
    ) -> f64 {
        let v = var as usize;
        Self::get_impl_var_lower_from(
            sum,
            b.var_lower[v],
            b.impl_var_lower[v],
            b.impl_var_lower_source[v],
        )
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::getImplVarUpper (lines 162-167, 6 lines)
    ///
    /// Private 4-argument overload. `std::min` -> `cpp_min`.
    fn get_impl_var_upper_from(
        sum: HighsInt,
        my_var_upper: f64,
        my_impl_var_upper: f64,
        my_impl_var_upper_source: HighsInt,
    ) -> f64 {
        if my_impl_var_upper_source == sum {
            my_var_upper
        } else {
            cpp_min(my_impl_var_upper, my_var_upper)
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::getImplVarLower (lines 169-174, 6 lines)
    ///
    /// Private 4-argument overload. `std::max` -> `cpp_max`.
    fn get_impl_var_lower_from(
        sum: HighsInt,
        my_var_lower: f64,
        my_impl_var_lower: f64,
        my_impl_var_lower_source: HighsInt,
    ) -> f64 {
        if my_impl_var_lower_source == sum {
            my_var_lower
        } else {
            cpp_max(my_impl_var_lower, my_var_lower)
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::update (lines 176-183, 8 lines)
    ///
    /// The `(numInfs, activity, direction, bound, coefficient)` overload.
    #[inline]
    fn update_dir(
        num_infs: &mut HighsInt,
        activity: &mut HighsCDouble,
        direction: HighsInt,
        bound: f64,
        coefficient: f64,
    ) {
        if bound.abs() != kHighsInf {
            *activity += (direction as f64) * HighsCDouble::from(bound) * coefficient;
        } else {
            *num_infs += direction;
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::update (lines 185-191, 7 lines)
    ///
    /// The `(numInfs, activity, oldBound, newBound, coefficient)` overload.
    #[inline]
    fn update_old_new(
        num_infs: &mut HighsInt,
        activity: &mut HighsCDouble,
        old_bound: f64,
        new_bound: f64,
        coefficient: f64,
    ) {
        if old_bound == new_bound {
            return;
        }
        Self::update_dir(num_infs, activity, -1, old_bound, coefficient);
        Self::update_dir(num_infs, activity, 1, new_bound, coefficient);
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::residual (lines 193-205, 13 lines)
    ///
    /// C++ default arguments: `boundVar = -1, boundVarCoefficient = kHighsInf,
    /// oldBoundVarBound = kHighsInf, newBoundVarBound = kHighsInf` (all call
    /// sites pass every argument).
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn residual(
        num_infs: &mut HighsInt,
        activity: &mut HighsCDouble,
        var_bound: f64,
        coefficient: f64,
        bound_var: HighsInt,
        bound_var_coefficient: f64,
        old_bound_var_bound: f64,
        new_bound_var_bound: f64,
    ) {
        Self::update_dir(num_infs, activity, -1, var_bound, coefficient);
        if bound_var != -1 {
            Self::update_old_new(
                num_infs,
                activity,
                old_bound_var_bound,
                new_bound_var_bound,
                bound_var_coefficient,
            );
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::handleVarUpper (lines 207-213, 7 lines)
    fn handle_var_upper(
        &mut self,
        sum: HighsInt,
        coefficient: f64,
        my_var_upper: f64,
        direction: HighsInt,
    ) {
        let s = sum as usize;
        if coefficient > 0.0 {
            Self::update_dir(
                &mut self.num_inf_sum_upper_orig[s],
                &mut self.sum_upper_orig[s],
                direction,
                my_var_upper,
                coefficient,
            );
        } else {
            Self::update_dir(
                &mut self.num_inf_sum_lower_orig[s],
                &mut self.sum_lower_orig[s],
                direction,
                my_var_upper,
                coefficient,
            );
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::handleVarLower (lines 215-221, 7 lines)
    fn handle_var_lower(
        &mut self,
        sum: HighsInt,
        coefficient: f64,
        my_var_lower: f64,
        direction: HighsInt,
    ) {
        let s = sum as usize;
        if coefficient > 0.0 {
            Self::update_dir(
                &mut self.num_inf_sum_lower_orig[s],
                &mut self.sum_lower_orig[s],
                direction,
                my_var_lower,
                coefficient,
            );
        } else {
            Self::update_dir(
                &mut self.num_inf_sum_upper_orig[s],
                &mut self.sum_upper_orig[s],
                direction,
                my_var_lower,
                coefficient,
            );
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::handleImplVarUpper (lines 223-229, 7 lines)
    fn handle_impl_var_upper(
        &mut self,
        sum: HighsInt,
        coefficient: f64,
        my_impl_var_upper: f64,
        direction: HighsInt,
    ) {
        let s = sum as usize;
        if coefficient > 0.0 {
            Self::update_dir(
                &mut self.num_inf_sum_upper[s],
                &mut self.sum_upper[s],
                direction,
                my_impl_var_upper,
                coefficient,
            );
        } else {
            Self::update_dir(
                &mut self.num_inf_sum_lower[s],
                &mut self.sum_lower[s],
                direction,
                my_impl_var_upper,
                coefficient,
            );
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::handleImplVarLower (lines 231-237, 7 lines)
    fn handle_impl_var_lower(
        &mut self,
        sum: HighsInt,
        coefficient: f64,
        my_impl_var_lower: f64,
        direction: HighsInt,
    ) {
        let s = sum as usize;
        if coefficient > 0.0 {
            Self::update_dir(
                &mut self.num_inf_sum_lower[s],
                &mut self.sum_lower[s],
                direction,
                my_impl_var_lower,
                coefficient,
            );
        } else {
            Self::update_dir(
                &mut self.num_inf_sum_upper[s],
                &mut self.sum_upper[s],
                direction,
                my_impl_var_lower,
                coefficient,
            );
        }
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::updatedImplVarUpper (lines 239-252, 14 lines)
    ///
    /// Private 6-argument overload.
    #[allow(clippy::too_many_arguments)]
    fn updated_impl_var_upper_full(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        old_var_upper: f64,
        old_impl_var_upper: f64,
        old_impl_var_upper_source: HighsInt,
    ) {
        let old_v_upper = Self::get_impl_var_upper_from(
            sum,
            old_var_upper,
            old_impl_var_upper,
            old_impl_var_upper_source,
        );
        let v_upper = self.get_impl_var_upper(b, sum, var);

        if v_upper == old_v_upper {
            return;
        }

        self.handle_impl_var_upper(sum, coefficient, old_v_upper, -1);
        self.handle_impl_var_upper(sum, coefficient, v_upper, 1);
    }

    /// HiGHS: HighsLinearSumBounds.cpp::HighsLinearSumBounds::updatedImplVarLower (lines 254-267, 14 lines)
    ///
    /// Private 6-argument overload.
    #[allow(clippy::too_many_arguments)]
    fn updated_impl_var_lower_full(
        &mut self,
        b: &HighsLinearSumBoundArrays,
        sum: HighsInt,
        var: HighsInt,
        coefficient: f64,
        old_var_lower: f64,
        old_impl_var_lower: f64,
        old_impl_var_lower_source: HighsInt,
    ) {
        let old_v_lower = Self::get_impl_var_lower_from(
            sum,
            old_var_lower,
            old_impl_var_lower,
            old_impl_var_lower_source,
        );
        let v_lower = self.get_impl_var_lower(b, sum, var);

        if v_lower == old_v_lower {
            return;
        }

        self.handle_impl_var_lower(sum, coefficient, old_v_lower, -1);
        self.handle_impl_var_lower(sum, coefficient, v_lower, 1);
    }
}
