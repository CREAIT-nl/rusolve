//! HiGHS source: mip/HighsTransformedLp.h, mip/HighsTransformedLp.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsTransformedLp.h::HighsTransformedLp (lines 29-65)
///
/// Removed C++ members:
/// * C++: `const HighsLpRelaxation & lprelaxation`: removed; methods that use it take
///   `lprelaxation: &HighsLpRelaxation` (the separation round's detached LP).
/// * C++: `const HighsDomain & globaldom_`: kept as the index [`HighsTransformedLp::globaldom_`]
///   into `ms.mipdata_.domains` (the worker's global domain). It must stay attached: the body
///   reads it and `lprelaxation.slack_lower(ms, row, globaldom_)` detaches it for the cut-pool
///   activity replay.
#[derive(Clone, Default, Debug)]
pub struct HighsTransformedLp {
    /// C++: `std::vector<std::pair<HighsInt, HighsImplications::VarBound>> bestVub`
    pub best_vub: Vec<(HighsInt, HighsImplications_VarBound)>,
    /// C++: `std::vector<std::pair<HighsInt, HighsImplications::VarBound>> bestVlb`
    pub best_vlb: Vec<(HighsInt, HighsImplications_VarBound)>,
    /// C++: `std::vector<double> simpleLbDist`
    pub simple_lb_dist: Vec<f64>,
    /// C++: `std::vector<double> simpleUbDist`
    pub simple_ub_dist: Vec<f64>,
    /// C++: `std::vector<double> lbDist`
    pub lb_dist: Vec<f64>,
    /// C++: `std::vector<double> ubDist`
    pub ub_dist: Vec<f64>,
    /// C++: `std::vector<double> boundDist`
    pub bound_dist: Vec<f64>,
    /// C++: `std::vector<BoundType> boundTypes`
    pub bound_types: Vec<HighsTransformedLp_BoundType>,
    /// C++: `HighsSparseVectorSum vectorsum`
    pub vectorsum: HighsSparseVectorSum,
    /// C++: `const HighsDomain & globaldom_` (index into `ms.mipdata_.domains`).
    pub globaldom_: usize,
}

/// HiGHS: HighsTransformedLp.h::HighsTransformedLp::BoundType (enum class)
/// HiGHS latest (#3179, 28339ce4b): HighsTransformedLp.h::HighsTransformedLp::BoundType (lines 41-46)
#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsTransformedLp_BoundType {
    #[default]
    kSimpleUb = 0,
    kSimpleLb = 1,
    kVariableUb = 2,
    kVariableLb = 3,
    kUnused = 4,
}

/// C++ `std::vector::erase` of the `remove` lambda in `transform`: moves the last entry into
/// `position` and zeroes the tail (HighsTransformedLp.cpp:145-150).
///
/// PORT NOTE (perf): the buffers are the `inds`/`vals` slices bound once
/// by [`HighsTransformedLp::transform`], whose length is the initial `num_nz`; the entries the C++
/// `std::vector::erase` names are the same ones.
fn remove_position(
    inds: &mut [HighsInt],
    vals: &mut [f64],
    num_nz: &mut HighsInt,
    position: usize,
) {
    *num_nz -= 1;
    let last = *num_nz as usize;
    // SAFETY: `transform` calls this only with `position < num_nz` as it was before the
    // decrement above, and `num_nz` never exceeds the initial `inds.len()`, which is the length of
    // both slices (`transform` slices them from the same cut row, and C++ `transform` requires
    // `vals.size() == inds.size()`). So `position` and `last = num_nz - 1` are both in range.
    // With `position == last` the sequence writes the entry and then zeroes it, as the C++ does.
    unsafe {
        let moved_ind = *inds.get_unchecked(last);
        let moved_val = *vals.get_unchecked(last);
        *inds.get_unchecked_mut(position) = moved_ind;
        *vals.get_unchecked_mut(position) = moved_val;
        *inds.get_unchecked_mut(last) = 0;
        *vals.get_unchecked_mut(last) = 0.0;
    }
}

/// C++ `getLb` lambda in `transform` (HighsTransformedLp.cpp:150-154). A slack's bound reads
/// `lprelaxation.slackLower`, which detaches the global domain for the cut-pool replay.
///
/// PORT NOTE (perf): `clean` is true when no pool event is pending for `globaldom`. The pending
/// test is then made once per `transform` call instead of once per slack bound. It is recomputed
/// after every replaying `slack_lower` and after every implication cleanup, the only calls in
/// `transform` that can change the pool logs or the cursors. The value is the same bit for bit.
///
/// PORT NOTE (perf): inlined, as the C++ lambda is.
#[inline(always)]
fn transformed_lb(
    ms: &mut HighsMipSolver,
    lprelaxation: &HighsLpRelaxation,
    globaldom: usize,
    slack_offset: HighsInt,
    col: HighsInt,
    clean: &mut bool,
) -> f64 {
    if col < slack_offset {
        return ms.mipdata_.domains[globaldom].col_lower_[col as usize];
    }
    let row = col - slack_offset;
    if *clean {
        return lprelaxation.slack_lower_no_pending(ms, row, globaldom);
    }
    let lb = lprelaxation.slack_lower(ms, row, globaldom);
    *clean = !ms.mipdata_.domains[globaldom].has_pending_pool_events(ms);
    lb
}

/// C++ `getUb` lambda in `transform` (HighsTransformedLp.cpp:155-159).
#[inline(always)]
fn transformed_ub(
    ms: &mut HighsMipSolver,
    lprelaxation: &HighsLpRelaxation,
    globaldom: usize,
    slack_offset: HighsInt,
    col: HighsInt,
) -> f64 {
    if col < slack_offset {
        ms.mipdata_.domains[globaldom].col_upper_[col as usize]
    } else {
        lprelaxation.slack_upper(ms, col - slack_offset, globaldom)
    }
}

/// C++ `implications.getBestVub(col, lpSolution, bestub, globaldom_)`. Both the implications
/// and the global domain (a `const HighsDomain &` argument) are detached for the call.
fn best_vub_of(
    ms: &mut HighsMipSolver,
    globaldom: usize,
    col: HighsInt,
    sol: &HighsSolution,
    best_ub: &mut f64,
) -> (HighsInt, HighsImplications_VarBound) {
    detach!(ms.mipdata_.domains[globaldom], |d| detach!(
        ms.mipdata_.implications,
        |imp| imp.get_best_vub(ms, col, sol, best_ub, d)
    ))
}

/// C++ `implications.getBestVlb(col, lpSolution, bestlb, globaldom_)`; see [`best_vub_of`].
fn best_vlb_of(
    ms: &mut HighsMipSolver,
    globaldom: usize,
    col: HighsInt,
    sol: &HighsSolution,
    best_lb: &mut f64,
) -> (HighsInt, HighsImplications_VarBound) {
    detach!(ms.mipdata_.domains[globaldom], |d| detach!(
        ms.mipdata_.implications,
        |imp| imp.get_best_vlb(ms, col, sol, best_lb, d)
    ))
}

impl HighsTransformedLp {
    /// HiGHS: HighsTransformedLp.cpp::HighsTransformedLp::HighsTransformedLp (lines 15-130)
    /// C++: `HighsTransformedLp(const HighsLpRelaxation & lprelaxation, HighsImplications & implications, const HighsDomain & globaldom)`
    // PORT NOTE: the `implications` parameter is dropped: the only C++ call site
    // (HighsSeparation.cpp:132) passes `mipdata.implications`, and the body mixes it with
    // `mipsolver.mipdata_->implications` freely. Each call is
    // `detach!(ms.mipdata_.implications, |imp| ..)`; `getBestVub/Vlb` also need the global
    // domain as `&HighsDomain`, so the body detaches `domains[globaldom]` around them, while
    // `cleanupVarbounds` reads `domains[0]` attached.
    pub fn new(
        ms: &mut HighsMipSolver,
        lprelaxation: &HighsLpRelaxation,
        globaldom: usize,
    ) -> Self {
        let mut this = HighsTransformedLp {
            globaldom_: globaldom,
            ..Default::default()
        };
        debug_assert!(HighsLpRelaxation::scaled_optimal(lprelaxation.get_status()));
        let sol = lprelaxation.get_lp_solver_2().get_solution();

        let num_transformed_col = lprelaxation.num_cols() + lprelaxation.num_rows();
        let n = num_transformed_col as usize;

        this.bound_dist = vec![0.0; n];
        this.simple_lb_dist = vec![0.0; n];
        this.simple_ub_dist = vec![0.0; n];
        this.lb_dist = vec![0.0; n];
        this.ub_dist = vec![0.0; n];
        this.best_vlb = vec![(-1, HighsImplications_VarBound::default()); n];
        this.best_vub = vec![(-1, HighsImplications_VarBound::default()); n];
        // HiGHS #3179 (mip_path_same_bds): HighsTransformedLp.h lines 41-46 at 28339ce4b put
        // kUnused first in BoundType, so `boundTypes.resize(numTransformedCol)` (cpp line 34)
        // value-initialises every column to kUnused; in 1.15.1 the first enumerator is
        // kSimpleUb. The initial type is visible as the restored `oldBoundType` and in the
        // second loop of `transform` (continuous columns with kUnused get a simple bound).
        let init_type = if ms.options_mip_.mip_path_same_bds {
            HighsTransformedLp_BoundType::kUnused
        } else {
            HighsTransformedLp_BoundType::kSimpleUb
        };
        this.bound_types = vec![init_type; n];
        this.vectorsum.set_dimension(num_transformed_col);

        for col in ms.mipdata_.continuous_cols.clone() {
            let ci = col as usize;
            if !ms.parallel_lock_active() {
                detach!(ms.mipdata_.implications, |imp| imp
                    .cleanup_varbounds(ms, col));
            }

            if ms.mipdata_.domains[globaldom].infeasible() {
                return this;
            }

            if ms.mipdata_.domains[globaldom].is_fixed(col) {
                continue;
            }

            let mut bestub = ms.mipdata_.domains[globaldom].col_upper_[ci];
            this.simple_ub_dist[ci] = bestub - sol.col_value[ci];
            if this.simple_ub_dist[ci] <= ms.mipdata_.feastol {
                this.simple_ub_dist[ci] = 0.0;
            }
            this.best_vub[ci] = best_vub_of(ms, globaldom, col, sol, &mut bestub);

            let mut bestlb = ms.mipdata_.domains[globaldom].col_lower_[ci];
            this.simple_lb_dist[ci] = sol.col_value[ci] - bestlb;
            if this.simple_lb_dist[ci] <= ms.mipdata_.feastol {
                this.simple_lb_dist[ci] = 0.0;
            }
            this.best_vlb[ci] = best_vlb_of(ms, globaldom, col, sol, &mut bestlb);

            this.lb_dist[ci] = sol.col_value[ci] - bestlb;
            if this.lb_dist[ci] <= ms.mipdata_.feastol {
                this.lb_dist[ci] = 0.0;
            }
            this.ub_dist[ci] = bestub - sol.col_value[ci];
            if this.ub_dist[ci] <= ms.mipdata_.feastol {
                this.ub_dist[ci] = 0.0;
            }

            this.bound_dist[ci] = cpp_min(this.lb_dist[ci], this.ub_dist[ci]);
        }

        for col in ms.mipdata_.integral_cols.clone() {
            let ci = col as usize;
            let mut bestub = ms.mipdata_.domains[globaldom].col_upper_[ci];
            let mut bestlb = ms.mipdata_.domains[globaldom].col_lower_[ci];

            if !ms.parallel_lock_active() {
                detach!(ms.mipdata_.implications, |imp| imp
                    .cleanup_varbounds(ms, col));
            }

            if ms.mipdata_.domains[globaldom].infeasible() {
                return this;
            }
            this.simple_ub_dist[ci] = bestub - sol.col_value[ci];
            if this.simple_ub_dist[ci] <= ms.mipdata_.feastol {
                this.simple_ub_dist[ci] = 0.0;
            }

            this.simple_lb_dist[ci] = sol.col_value[ci] - bestlb;
            if this.simple_lb_dist[ci] <= ms.mipdata_.feastol {
                this.simple_lb_dist[ci] = 0.0;
            }
            let simple_bnd_dist = cpp_min(this.simple_lb_dist[ci], this.simple_ub_dist[ci]);
            if simple_bnd_dist > 0.0
                && ((HighsIntegers::nearest_integer(sol.col_value[ci]) as f64) - sol.col_value[ci])
                    .abs()
                    < ms.mipdata_.feastol
            {
                this.best_vub[ci] = best_vub_of(ms, globaldom, col, sol, &mut bestub);
                this.best_vlb[ci] = best_vlb_of(ms, globaldom, col, sol, &mut bestlb);

                this.lb_dist[ci] = sol.col_value[ci] - bestlb;
                if this.lb_dist[ci] <= ms.mipdata_.feastol {
                    this.lb_dist[ci] = 0.0;
                }
                this.ub_dist[ci] = bestub - sol.col_value[ci];
                if this.ub_dist[ci] <= ms.mipdata_.feastol {
                    this.ub_dist[ci] = 0.0;
                }
                this.bound_dist[ci] = cpp_min(this.lb_dist[ci], this.ub_dist[ci]);
                if this.bound_dist[ci] > simple_bnd_dist + ms.mipdata_.feastol {
                    this.lb_dist[ci] = this.simple_lb_dist[ci];
                    this.ub_dist[ci] = this.simple_ub_dist[ci];
                    this.bound_dist[ci] = simple_bnd_dist;
                    this.best_vub[ci].0 = -1;
                    this.best_vlb[ci].0 = -1;
                }
            } else {
                this.lb_dist[ci] = this.simple_lb_dist[ci];
                this.ub_dist[ci] = this.simple_ub_dist[ci];
                this.bound_dist[ci] = simple_bnd_dist;
            }
        }

        // setup information of slackVariables
        let num_lp_row = lprelaxation.num_rows();
        let index_offset = ms.num_col();
        for row in 0..num_lp_row {
            let slack_index = (index_offset + row) as usize;
            let bestub = lprelaxation.slack_upper(ms, row, globaldom);
            let bestlb = lprelaxation.slack_lower(ms, row, globaldom);

            if bestlb == bestub {
                continue;
            }

            this.lb_dist[slack_index] = sol.row_value[row as usize] - bestlb;
            if this.lb_dist[slack_index] <= ms.mipdata_.feastol {
                this.lb_dist[slack_index] = 0.0;
            }
            this.simple_lb_dist[slack_index] = this.lb_dist[slack_index];

            this.ub_dist[slack_index] = bestub - sol.row_value[row as usize];
            if this.ub_dist[slack_index] <= ms.mipdata_.feastol {
                this.ub_dist[slack_index] = 0.0;
            }
            this.simple_ub_dist[slack_index] = this.ub_dist[slack_index];

            this.bound_dist[slack_index] =
                cpp_min(this.lb_dist[slack_index], this.ub_dist[slack_index]);
        }

        this
    }

    /// HiGHS: HighsTransformedLp.cpp::HighsTransformedLp::transform (lines 132-447)
    /// HiGHS latest (#3177, 26e3f354b): HighsTransformedLp.cpp::HighsTransformedLp::transform (lines 195-203, 215-223)
    /// HiGHS latest (#3179, 28339ce4b): HighsTransformedLp.cpp::HighsTransformedLp::transform (lines 132-489)
    /// C++: `bool transform(std::vector<double> & vals, std::vector<double> & upper, std::vector<double> & solval, std::vector<HighsInt> & inds, double & rhs, bool & integersPositive, bool preferVbds = false, bool enforceSameBds = false)`
    /// C++ default arguments: `prefer_vbds = false`, `enforce_same_bds = false`
    // PORT NOTE: `cleanupVub` / `cleanupVlb` detach `ms.mipdata_.implications`.
    #[allow(clippy::too_many_arguments)]
    pub fn transform(
        &mut self,
        ms: &mut HighsMipSolver,
        lprelaxation: &HighsLpRelaxation,
        vals: &mut Vec<f64>,
        upper: &mut Vec<f64>,
        solval: &mut Vec<f64>,
        inds: &mut Vec<HighsInt>,
        rhs: &mut f64,
        integers_positive: &mut bool,
        prefer_vbds: bool,
        enforce_same_bds: bool,
    ) -> bool {
        // vector sum should be empty
        debug_assert!(self.vectorsum.get_nonzeros().is_empty());

        let mut tmp_rhs = HighsCDouble::from(*rhs);

        let slack_offset = lprelaxation.num_cols();
        let globaldom = self.globaldom_;
        let mut clean = !ms.mipdata_.domains[globaldom].has_pending_pool_events(ms);

        let mut num_nz = inds.len() as HighsInt;

        // HiGHS #3177 (mip_redundant_vbd): HighsTransformedLp.cpp lines 199-203, 219-223 at 26e3f354b
        let redundant_vbd = ms.options_mip_.mip_redundant_vbd;
        // HiGHS #3264 (mip_aggregate_implications): inf-safe VarBound::maxValue/minValue
        // (HighsImplications.h lines 56-67 at f2ae2b320) in the cleanup tests below.
        let inf_safe_vbd = ms.options_mip_.mip_aggregate_implications;

        // PORT NOTE (perf): loop invariants. C++ keeps
        // `mip.options_mip_->small_matrix_value` and `mip.mipdata_->feastol` in registers over the
        // whole of `transform`; the port re-read them through `ms` at every nonzero (`feastol`
        // twice in the variable-bound cleanup tests and once per iteration below the
        // bound-type comment). Nothing in `transform` writes an option or `feastol`, so the value
        // is the same bit for bit.
        let small_matrix_value = ms.options_mip_.small_matrix_value;
        let feastol = ms.mipdata_.feastol;

        // PORT NOTE (perf): the first loop's working buffers are bound
        // once as the slices of the initial `inds`/`vals` range, so the per-nonzero accesses are
        // unchecked (see [`remove_position`]). `num_nz` starts at `inds.len()` and only decreases,
        // and `remove` only moves the entry at `num_nz - 1` - inside the range - into position `i`,
        // so `i` and `num_nz - 1` name the same entries as the C++ `inds[i]` / `vals[i]` for the
        // whole loop. `inds` and `vals` are the same cut row at all three call sites, as C++
        // requires, so the two slices have the same length.
        let num_nz0 = num_nz as usize;
        let (inds_buf, vals_buf) = (&mut inds[..num_nz0], &mut vals[..num_nz0]);

        let mut i: HighsInt = 0;
        while i < num_nz {
            let iu = i as usize;
            debug_assert!(iu < num_nz0);
            // SAFETY: `iu < num_nz <= num_nz0 = inds_buf.len()`; see the note above.
            let col = unsafe { *inds_buf.get_unchecked(iu) };
            let cu = col as usize;
            // SAFETY: as for `inds_buf`; `vals` has the same length.
            let mut val = unsafe { *vals_buf.get_unchecked(iu) };

            let lb = transformed_lb(ms, lprelaxation, globaldom, slack_offset, col, &mut clean);
            let ub = transformed_ub(ms, lprelaxation, globaldom, slack_offset, col);

            if ub - lb < small_matrix_value {
                tmp_rhs -= cpp_min(lb, ub) * val;
                remove_position(inds_buf, vals_buf, &mut num_nz, iu);
                continue;
            }

            if lb == -kHighsInf && ub == kHighsInf {
                self.vectorsum.clear();
                return false;
            }

            // the code below uses the difference between the column upper and lower
            // bounds as the upper bound for the slack from the variable upper bound
            // constraint (upper[j] = ub - lb) and thus assumes that the variable upper
            // bound constraints are tight. this assumption may not be satisfied when
            // new bound changes were derived during cut generation and, therefore, we
            // tighten the best variable upper bound.
            // PORT NOTE (perf): `bestVub[col].first` /
            // `bestVlb[col].first` are read at the two cleanup tests, possibly written by the
            // redundant branches below, and read again in the bound-type selection and the match;
            // bound once here and kept in step with the two writes. `cleanupVub`/`cleanupVlb` get
            // the second tuple entry by reference and the first by value, as in C++, so they
            // cannot change it.
            let mut vub0 = self.best_vub[cu].0;
            if vub0 != -1 && self.best_vub[cu].1.max_value_x(inf_safe_vbd) > ub + feastol {
                let mut redundant = false;
                let mut infeasible = false;
                detach!(ms.mipdata_.implications, |imp| imp.cleanup_vub(
                    ms,
                    col,
                    vub0,
                    &mut self.best_vub[cu].1,
                    ub,
                    &mut redundant,
                    &mut infeasible,
                    false
                ));
                clean = !ms.mipdata_.domains[globaldom].has_pending_pool_events(ms);
                if redundant_vbd && redundant {
                    self.best_vub[cu].0 = -1;
                    vub0 = -1;
                    self.ub_dist[cu] = self.simple_ub_dist[cu];
                    self.bound_dist[cu] = cpp_min(self.lb_dist[cu], self.ub_dist[cu]);
                }
            }

            // the code below uses the difference between the column upper and lower
            // bounds as the upper bound for the slack from the variable lower bound
            // constraint (upper[j] = ub - lb) and thus assumes that the variable lower
            // bound constraints are tight. this assumption may not be satisfied when
            // new bound changes were derived during cut generation and, therefore, we
            // tighten the best variable lower bound.
            let mut vlb0 = self.best_vlb[cu].0;
            if vlb0 != -1 && self.best_vlb[cu].1.min_value_x(inf_safe_vbd) < lb - feastol {
                let mut redundant = false;
                let mut infeasible = false;
                detach!(ms.mipdata_.implications, |imp| imp.cleanup_vlb(
                    ms,
                    col,
                    vlb0,
                    &mut self.best_vlb[cu].1,
                    lb,
                    &mut redundant,
                    &mut infeasible,
                    false
                ));
                clean = !ms.mipdata_.domains[globaldom].has_pending_pool_events(ms);
                if redundant_vbd && redundant {
                    self.best_vlb[cu].0 = -1;
                    vlb0 = -1;
                    self.lb_dist[cu] = self.simple_lb_dist[cu];
                    self.bound_dist[cu] = cpp_min(self.lb_dist[cu], self.ub_dist[cu]);
                }
            }

            // store the old bound type so that we can restore it if the continuous
            // column is relaxed out anyways. This allows to correctly transform and
            // then untransform multiple base rows which is useful to compute cuts based
            // on several transformed base rows. It could otherwise lead to bugs if a
            // column is first transformed with a simple bound and not relaxed but for
            // another base row is transformed and relaxed with a variable bound. Should
            // the non-relaxed column now be untransformed we would wrongly use the
            // variable bound even though this is not the correct way to untransform the
            // column.
            let old_bound_type = self.bound_types[cu];

            // PORT NOTE (perf): the bound-distance arrays are read three
            // to six times per nonzero through `self`, each read with its own length check. They
            // are written only in the redundant-variable-bound branches, which ran before this
            // point, and no call in the bound-type selection or in the match touches them, so
            // binding them once here reads the same values as the re-reads did. `bound_types` is
            // written only here and read only by the match below, which now uses the local; it is
            // stored back before the match, exactly where the C++ has already assigned
            // `boundTypes[col]`.
            let lb_dist = self.lb_dist[cu];
            let ub_dist = self.ub_dist[cu];
            let simple_lb_dist = self.simple_lb_dist[cu];
            let simple_ub_dist = self.simple_ub_dist[cu];
            let bound_dist = self.bound_dist[cu];
            let mut bound_type = old_bound_type;

            if enforce_same_bds && old_bound_type != HighsTransformedLp_BoundType::kUnused {
                // Don't allow different bound types to be used in stacked
                // transformations unless they're going to be relaxed out
                // This ensures that a single call to untransform remains valid
                if lprelaxation.is_col_integral(ms, col)
                    && (ub - lb <= 1.5
                        || bound_dist != 0.0
                        || simple_lb_dist == 0.0
                        || simple_ub_dist == 0.0)
                {
                    i += 1;
                    continue;
                }
                if old_bound_type == HighsTransformedLp_BoundType::kSimpleLb
                    || old_bound_type == HighsTransformedLp_BoundType::kSimpleUb
                {
                    if lb_dist < ub_dist && val > 0.0 {
                        bound_type = HighsTransformedLp_BoundType::kSimpleLb;
                    } else if ub_dist < lb_dist && val < 0.0 {
                        bound_type = HighsTransformedLp_BoundType::kSimpleUb;
                    }
                }
            } else if lprelaxation.is_col_integral(ms, col) {
                if ub - lb <= 1.5
                    || bound_dist != 0.0
                    || simple_lb_dist == 0.0
                    || simple_ub_dist == 0.0
                {
                    // since we skip the handling of variable bound constraints for all
                    // binary and some general-integer variables here, the bound type used
                    // should be a simple lower or upper bound
                    bound_type = if simple_lb_dist < simple_ub_dist - feastol {
                        HighsTransformedLp_BoundType::kSimpleLb
                    } else if simple_ub_dist < simple_lb_dist - feastol {
                        HighsTransformedLp_BoundType::kSimpleUb
                    } else if val > 0.0 {
                        HighsTransformedLp_BoundType::kSimpleLb
                    } else {
                        HighsTransformedLp_BoundType::kSimpleUb
                    };
                    self.bound_types[cu] = bound_type;
                    i += 1;
                    continue;
                }
                if vlb0 == -1 || ub_dist < lb_dist - feastol {
                    debug_assert!(vub0 != -1);
                    bound_type = HighsTransformedLp_BoundType::kVariableUb;
                } else if vub0 == -1 || lb_dist < ub_dist - feastol {
                    debug_assert!(vlb0 != -1);
                    bound_type = HighsTransformedLp_BoundType::kVariableLb;
                } else if val > 0.0 {
                    debug_assert!(vub0 != -1);
                    bound_type = HighsTransformedLp_BoundType::kVariableUb;
                } else {
                    debug_assert!(vlb0 != -1);
                    bound_type = HighsTransformedLp_BoundType::kVariableLb;
                }
            } else if lb_dist < ub_dist - feastol {
                bound_type = if vlb0 == -1 {
                    HighsTransformedLp_BoundType::kSimpleLb
                } else if prefer_vbds || val > 0.0 || simple_lb_dist > lb_dist + feastol {
                    HighsTransformedLp_BoundType::kVariableLb
                } else {
                    HighsTransformedLp_BoundType::kSimpleLb
                };
            } else if ub_dist < lb_dist - feastol {
                bound_type = if vub0 == -1 {
                    HighsTransformedLp_BoundType::kSimpleUb
                } else if prefer_vbds || val < 0.0 || simple_ub_dist > ub_dist + feastol {
                    HighsTransformedLp_BoundType::kVariableUb
                } else {
                    HighsTransformedLp_BoundType::kSimpleUb
                };
            } else if val > 0.0 {
                bound_type = if vlb0 != -1 {
                    HighsTransformedLp_BoundType::kVariableLb
                } else if prefer_vbds && vub0 != -1 {
                    HighsTransformedLp_BoundType::kVariableUb
                } else {
                    HighsTransformedLp_BoundType::kSimpleLb
                };
            } else {
                bound_type = if vub0 != -1 {
                    HighsTransformedLp_BoundType::kVariableUb
                } else if prefer_vbds && vlb0 != -1 {
                    HighsTransformedLp_BoundType::kVariableLb
                } else {
                    HighsTransformedLp_BoundType::kSimpleUb
                };
            }

            self.bound_types[cu] = bound_type;

            match bound_type {
                HighsTransformedLp_BoundType::kSimpleLb => {
                    if val > 0.0 {
                        // relax away using lower bound
                        tmp_rhs -= lb * val;
                        self.bound_types[cu] = old_bound_type;
                        remove_position(inds_buf, vals_buf, &mut num_nz, iu);
                        continue;
                    }
                }
                HighsTransformedLp_BoundType::kSimpleUb => {
                    if val < 0.0 {
                        // relax away using upper bound
                        tmp_rhs -= ub * val;
                        self.bound_types[cu] = old_bound_type;
                        remove_position(inds_buf, vals_buf, &mut num_nz, iu);
                        continue;
                    }
                }
                HighsTransformedLp_BoundType::kVariableLb => {
                    let constant = self.best_vlb[cu].1.constant;
                    let coef = self.best_vlb[cu].1.coef;
                    tmp_rhs -= constant * val;
                    self.vectorsum.add(vlb0, val * coef);
                    if val > 0.0 {
                        self.bound_types[cu] = old_bound_type;
                        remove_position(inds_buf, vals_buf, &mut num_nz, iu);
                        continue;
                    }
                }
                HighsTransformedLp_BoundType::kVariableUb => {
                    let constant = self.best_vub[cu].1.constant;
                    let coef = self.best_vub[cu].1.coef;
                    tmp_rhs -= constant * val;
                    self.vectorsum.add(vub0, val * coef);
                    val = -val;
                    // SAFETY: `iu < num_nz0 = vals_buf.len()`; see the note at the loop head.
                    unsafe { *vals_buf.get_unchecked_mut(iu) = val };
                    if val > 0.0 {
                        self.bound_types[cu] = old_bound_type;
                        remove_position(inds_buf, vals_buf, &mut num_nz, iu);
                        continue;
                    }
                }
                HighsTransformedLp_BoundType::kUnused => {
                    debug_assert!(false);
                    self.vectorsum.clear();
                    return false;
                }
            }
            // move to next element
            i += 1;
        }

        if !self.vectorsum.get_nonzeros().is_empty() {
            // PORT NOTE (perf): the live prefix of the first loop's
            // buffers, i.e. the C++ `for (HighsInt j = 0; j < numNz; j++)`.
            let live = num_nz as usize;
            for (&ind, &val) in inds_buf[..live].iter().zip(vals_buf[..live].iter()) {
                if val != 0.0 {
                    self.vectorsum.add(ind, val);
                }
            }

            let max_error = 0.0_f64;
            // HiGHS #3223 (mip_cdouble_proofs): HighsTransformedLp.cpp lines 377-379 at 7db1c8180
            let x1b_cdouble = ms.options_mip_.mip_cdouble_proofs;
            let epsilon = ms.mipdata_.epsilon;
            let is_zero = |col: HighsInt, val: f64| {
                if x1b_cdouble {
                    // C++ getLb / getUb lambdas, as in the bound reads of the first loop
                    let lb =
                        transformed_lb(ms, lprelaxation, globaldom, slack_offset, col, &mut clean);
                    let ub = transformed_ub(ms, lprelaxation, globaldom, slack_offset, col);
                    return HighsDomain::term_is_zero(val, lb, ub, epsilon);
                }
                val.abs() <= small_matrix_value
            };

            self.vectorsum.cleanup(is_zero);
            if max_error > feastol {
                self.vectorsum.clear();
                return false;
            }

            // C++ `inds = vectorsum.getNonzeros()` copy-assigns into the existing capacity.
            inds.clear();
            inds.extend_from_slice(self.vectorsum.get_nonzeros());
            num_nz = inds.len() as HighsInt;

            vals.resize(num_nz as usize, 0.0);
            for (ind, val) in inds.iter().zip(vals.iter_mut()) {
                *val = self.vectorsum.get_value(*ind);
            }

            self.vectorsum.clear();
        } else {
            vals.resize(num_nz as usize, 0.0);
            inds.resize(num_nz as usize, 0);
        }

        // PORT NOTE (perf): both remaining loops walk the slices
        // instead of indexing the `&mut Vec` parameters: after the resizes above, `inds`, `vals`,
        // `upper` and `solval` all have exactly `num_nz` entries, so zipping their iterators
        // visits the same `j` as the C++ `for (HighsInt j = 0; j < numNz; j++)` and drops the
        // per-access length checks. `*integers_positive` is only read in this loop.
        let integers_positive_in = *integers_positive;
        for (&col, &val) in inds.iter().zip(vals.iter()) {
            let cu = col as usize;
            // PORT NOTE (perf): `boundTypes[col]` is read three times
            // here and written once; nothing between the reads writes it, so they use one local.
            let bound_type = self.bound_types[cu];

            if enforce_same_bds && bound_type != HighsTransformedLp_BoundType::kUnused {
                continue;
            }

            // get bounds
            let lb = transformed_lb(ms, lprelaxation, globaldom, slack_offset, col, &mut clean);
            let ub = transformed_ub(ms, lprelaxation, globaldom, slack_offset, col);

            // variable should not be free
            debug_assert!(lb != -kHighsInf || ub != kHighsInf);

            // set bound type for previously unprocessed integer-constrained variables
            if !lprelaxation.is_col_integral(ms, col)
                && bound_type != HighsTransformedLp_BoundType::kUnused
            {
                continue;
            }

            // do not overwrite bound type for integral slacks from vlb / vub
            // constraints
            if matches!(
                bound_type,
                HighsTransformedLp_BoundType::kVariableLb
                    | HighsTransformedLp_BoundType::kVariableUb
            ) {
                continue;
            }

            // complement integers to make coefficients positive if both bounds are
            // finite; otherwise, complement integers with closest bound.
            // take into account 'integersPositive' provided by caller.
            if integers_positive_in {
                self.bound_types[cu] = if (lb != -kHighsInf && val > 0.0) || ub == kHighsInf {
                    HighsTransformedLp_BoundType::kSimpleLb
                } else {
                    HighsTransformedLp_BoundType::kSimpleUb
                };
            } else if self.lb_dist[cu] < self.ub_dist[cu] {
                self.bound_types[cu] = HighsTransformedLp_BoundType::kSimpleLb;
            } else {
                self.bound_types[cu] = HighsTransformedLp_BoundType::kSimpleUb;
            }
        }

        upper.resize(num_nz as usize, 0.0);
        solval.resize(num_nz as usize, 0.0);

        // PORT NOTE (perf): `*integers_positive` is folded per nonzero
        // in C++; held in a local and written back at the end and before the one early return, so
        // the caller sees the same value on both paths.
        let mut ints_positive = *integers_positive;
        for ((col, val), (up, sv)) in inds
            .iter()
            .zip(vals.iter_mut())
            .zip(upper.iter_mut().zip(solval.iter_mut()))
        {
            let col = *col;
            let cu = col as usize;

            let lb = transformed_lb(ms, lprelaxation, globaldom, slack_offset, col, &mut clean);
            let ub = transformed_ub(ms, lprelaxation, globaldom, slack_offset, col);

            *up = ub - lb;

            match self.bound_types[cu] {
                HighsTransformedLp_BoundType::kSimpleLb => {
                    // shift (lower bound)
                    debug_assert!(lb != -kHighsInf);
                    tmp_rhs -= lb * *val;
                    *sv = self.lb_dist[cu];
                }
                HighsTransformedLp_BoundType::kSimpleUb => {
                    // complement (upper bound)
                    debug_assert!(ub != kHighsInf);
                    tmp_rhs -= ub * *val;
                    *val = -*val;
                    *sv = self.ub_dist[cu];
                }
                HighsTransformedLp_BoundType::kVariableLb => {
                    *sv = self.lb_dist[cu];
                }
                HighsTransformedLp_BoundType::kVariableUb => {
                    *sv = self.ub_dist[cu];
                }
                HighsTransformedLp_BoundType::kUnused => {
                    debug_assert!(false);
                    *integers_positive = ints_positive;
                    return false;
                }
            }

            // check if all integer-constrained variables have positive coefficients
            if lprelaxation.is_col_integral(ms, col) {
                ints_positive = ints_positive && *val > 0.0;
            }
        }

        *integers_positive = ints_positive;
        *rhs = f64::from(tmp_rhs);

        if num_nz == 0 && *rhs >= -feastol {
            return false;
        }

        true
    }

    /// HiGHS: HighsTransformedLp.cpp::HighsTransformedLp::untransform (lines 449-570)
    /// HiGHS latest (#3179, 28339ce4b): HighsTransformedLp.cpp::HighsTransformedLp::untransform (lines 491-570)
    /// C++: `bool untransform(std::vector<double> & vals, std::vector<HighsInt> & inds, double & rhs, bool integral = false)`
    /// C++ default arguments: `integral = false`
    pub fn untransform(
        &mut self,
        ms: &mut HighsMipSolver,
        lprelaxation: &HighsLpRelaxation,
        vals: &mut Vec<f64>,
        inds: &mut Vec<HighsInt>,
        rhs: &mut f64,
        integral: bool,
    ) -> bool {
        let mut tmp_rhs = HighsCDouble::from(*rhs);
        let slack_offset = ms.num_col();
        let globaldom = self.globaldom_;

        let num_nz = inds.len() as HighsInt;

        for i in 0..num_nz as usize {
            if vals[i] == 0.0 {
                continue;
            }
            let col = inds[i];
            let cu = col as usize;

            match self.bound_types[cu] {
                HighsTransformedLp_BoundType::kVariableLb => {
                    tmp_rhs += self.best_vlb[cu].1.constant * vals[i];
                    self.vectorsum
                        .add(self.best_vlb[cu].0, -vals[i] * self.best_vlb[cu].1.coef);
                    self.vectorsum.add(col, vals[i]);
                }
                HighsTransformedLp_BoundType::kVariableUb => {
                    tmp_rhs -= self.best_vub[cu].1.constant * vals[i];
                    self.vectorsum
                        .add(self.best_vub[cu].0, vals[i] * self.best_vub[cu].1.coef);
                    self.vectorsum.add(col, -vals[i]);
                }
                HighsTransformedLp_BoundType::kSimpleLb => {
                    if col < slack_offset {
                        tmp_rhs += vals[i] * ms.mipdata_.domains[globaldom].col_lower_[cu];
                        self.vectorsum.add(col, vals[i]);
                    } else {
                        let row = col - slack_offset;
                        tmp_rhs += vals[i] * lprelaxation.slack_lower(ms, row, globaldom);

                        let mut rowlen: HighsInt = 0;
                        let mut rowinds: &[HighsInt] = &[];
                        let mut rowvals: &[f64] = &[];
                        lprelaxation.get_row(ms, row, &mut rowlen, &mut rowinds, &mut rowvals);

                        for j in 0..rowlen as usize {
                            self.vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
                        }
                    }
                }
                HighsTransformedLp_BoundType::kSimpleUb => {
                    if col < slack_offset {
                        tmp_rhs -= vals[i] * ms.mipdata_.domains[globaldom].col_upper_[cu];
                        self.vectorsum.add(col, -vals[i]);
                    } else {
                        let row = col - slack_offset;
                        tmp_rhs -= vals[i] * lprelaxation.slack_upper(ms, row, globaldom);
                        vals[i] = -vals[i];

                        let mut rowlen: HighsInt = 0;
                        let mut rowinds: &[HighsInt] = &[];
                        let mut rowvals: &[f64] = &[];
                        lprelaxation.get_row(ms, row, &mut rowlen, &mut rowinds, &mut rowvals);

                        for j in 0..rowlen as usize {
                            self.vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
                        }
                    }
                }
                HighsTransformedLp_BoundType::kUnused => {
                    debug_assert!(false);
                    self.vectorsum.clear();
                    return false;
                }
            }
        }

        if integral {
            // if the cut is integral, we just round all coefficient values and the
            // right hand side to the nearest integral value, as small deviation
            // only come from numerical errors during resubstitution of slack variables
            let is_zero = |_col: HighsInt, val: f64| val.abs() < 0.5;

            self.vectorsum.cleanup(is_zero);
            *rhs = f64::from(tmp_rhs).round();
        } else {
            let mut abort = false;
            // HiGHS #3223 (mip_cdouble_proofs): HighsTransformedLp.cpp lines 580-583 at 7db1c8180
            let x1b_cdouble = ms.options_mip_.mip_cdouble_proofs;
            let small_matrix_value = ms.options_mip_.small_matrix_value;
            let epsilon = ms.mipdata_.epsilon;
            let feastol = ms.mipdata_.feastol;
            let globaldom_domain = &ms.mipdata_.domains[globaldom];
            let is_zero = |_col: HighsInt, val: f64| {
                if x1b_cdouble {
                    if HighsDomain::term_is_zero(
                        val,
                        globaldom_domain.col_lower_[_col as usize],
                        globaldom_domain.col_upper_[_col as usize],
                        epsilon,
                    ) {
                        return true;
                    }
                }
                let absval = val.abs();
                if !x1b_cdouble && absval <= small_matrix_value {
                    return true;
                }

                if absval <= feastol {
                    if val > 0.0 {
                        if globaldom_domain.col_lower_[_col as usize] == -kHighsInf {
                            abort = true;
                        } else {
                            tmp_rhs -= val * globaldom_domain.col_lower_[_col as usize];
                        }
                    } else if globaldom_domain.col_upper_[_col as usize] == kHighsInf {
                        abort = true;
                    } else {
                        tmp_rhs -= val * globaldom_domain.col_upper_[_col as usize];
                    }
                    return true;
                }
                false
            };

            self.vectorsum.cleanup(is_zero);
            if abort {
                self.vectorsum.clear();
                return false;
            }
            *rhs = f64::from(tmp_rhs);
        }

        // C++ `inds = vectorsum.getNonzeros()` copy-assigns into the existing capacity.
        inds.clear();
        inds.extend_from_slice(self.vectorsum.get_nonzeros());
        let num_nz = inds.len();
        vals.resize(num_nz, 0.0);

        if integral {
            for i in 0..num_nz {
                vals[i] = self.vectorsum.get_value(inds[i]).round();
            }
        } else {
            for i in 0..num_nz {
                vals[i] = self.vectorsum.get_value(inds[i]);
            }
        }
        self.vectorsum.clear();

        true
    }
    /// HiGHS latest (#3179, 28339ce4b): HighsTransformedLp.h::HighsTransformedLp::initMultiRowTransform (lines 67-69)
    /// C++: `void initMultiRowTransform()`
    pub fn init_multi_row_transform(&mut self) {
        self.bound_types.fill(HighsTransformedLp_BoundType::kUnused);
    }

    /// HiGHS: HighsTransformedLp.h::HighsTransformedLp::boundDistance (lines 55-55)
    /// C++: `double boundDistance(HighsInt col) const`
    pub fn bound_distance(&self, col: HighsInt) -> f64 {
        self.bound_dist[col as usize]
    }

    /// HiGHS: HighsTransformedLp.h::HighsTransformedLp::getGlobaldom (lines 64-64)
    /// C++: `const HighsDomain & getGlobaldom() const`
    pub fn get_globaldom<'s>(&self, ms: &'s HighsMipSolver) -> &'s HighsDomain {
        &ms.mipdata_.domains[self.globaldom_]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Columns `x0 in [0,1]` (int), `x1 in [0,3]` (int), `x2 in [1,5]` (continuous) and the
    /// base row `2 x0 + 3 x1 - x2 <= 4`, with LP distances that make `transform` shift x0 and
    /// x2 to their lower bounds and complement x1 at its upper bound.
    #[test]
    fn transform_untransform_round_trip() {
        let mut ms = HighsMipSolver::default();
        {
            let lp = &mut ms.orig_model_;
            lp.num_col_ = 3;
            lp.num_row_ = 0;
            lp.col_cost_ = vec![0.0; 3];
            lp.col_lower_ = vec![0.0, 0.0, 1.0];
            lp.col_upper_ = vec![1.0, 3.0, 5.0];
            lp.integrality_ = vec![
                HighsVarType::kInteger,
                HighsVarType::kInteger,
                HighsVarType::kContinuous,
            ];
            lp.a_matrix_.num_col_ = 3;
            lp.a_matrix_.start_ = vec![0, 0, 0, 0];
        }
        ms.mipdata_.a_rstart_ = vec![0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        let dom = HighsDomain::new(&ms);
        ms.mipdata_.domains = vec![Slot::new(dom)];

        let mut lprel = HighsLpRelaxation::default();
        lprel.lpsolver.model_lp_.num_col_ = 3;

        let mut tlp = HighsTransformedLp {
            globaldom_: 0,
            lb_dist: vec![0.2, 2.5, 0.5],
            ub_dist: vec![0.8, 0.5, 3.5],
            simple_lb_dist: vec![0.2, 2.5, 0.5],
            simple_ub_dist: vec![0.8, 0.5, 3.5],
            bound_dist: vec![0.2, 0.5, 0.5],
            best_vlb: vec![(-1, HighsImplications_VarBound::default()); 3],
            best_vub: vec![(-1, HighsImplications_VarBound::default()); 3],
            bound_types: vec![HighsTransformedLp_BoundType::default(); 3],
            ..Default::default()
        };
        tlp.vectorsum.set_dimension(3);

        let mut inds = vec![0, 1, 2];
        let mut vals = vec![2.0, 3.0, -1.0];
        let (mut upper, mut solval) = (Vec::new(), Vec::new());
        let mut rhs = 4.0;
        let mut integers_positive = false;
        assert!(tlp.transform(
            &mut ms,
            &lprel,
            &mut vals,
            &mut upper,
            &mut solval,
            &mut inds,
            &mut rhs,
            &mut integers_positive,
            false,
            false,
        ));
        use HighsTransformedLp_BoundType::*;
        assert_eq!(tlp.bound_types, vec![kSimpleLb, kSimpleUb, kSimpleLb]);
        assert_eq!(inds, vec![0, 1, 2]);
        assert_eq!(vals, vec![2.0, -3.0, -1.0]);
        assert_eq!(upper, vec![1.0, 3.0, 4.0]);
        assert_eq!(solval, vec![0.2, 0.5, 0.5]);
        assert_eq!(rhs, -4.0);
        assert!(!integers_positive);

        assert!(tlp.untransform(&mut ms, &lprel, &mut vals, &mut inds, &mut rhs, false));
        assert_eq!(inds, vec![0, 1, 2]);
        assert_eq!(vals, vec![2.0, 3.0, -1.0]);
        assert_eq!(rhs, 4.0);
        assert!(tlp.vectorsum.get_nonzeros().is_empty());
    }

    /// #3177: with `mip_redundant_vbd` on, a redundant VUB resets `best_vub`,
    /// `ub_dist` and `bound_dist`.
    #[test]
    fn x1_tlp_redundant_vub_resets_dist() {
        let mut ms = HighsMipSolver::default();
        {
            let lp = &mut ms.orig_model_;
            lp.num_col_ = 1;
            lp.num_row_ = 0;
            lp.col_cost_ = vec![0.0; 1];
            lp.col_lower_ = vec![0.0];
            lp.col_upper_ = vec![10.0];
            lp.integrality_ = vec![HighsVarType::kContinuous];
            lp.a_matrix_.num_col_ = 1;
            lp.a_matrix_.start_ = vec![0, 0];
        }
        ms.mipdata_.a_rstart_ = vec![0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.mip_redundant_vbd = true;
        let dom = HighsDomain::new(&ms);
        ms.mipdata_.domains = vec![Slot::new(dom)];

        let mut lprel = HighsLpRelaxation::default();
        lprel.lpsolver.model_lp_.num_col_ = 1;

        let mut tlp = HighsTransformedLp {
            globaldom_: 0,
            lb_dist: vec![0.5],
            ub_dist: vec![5.0],
            simple_lb_dist: vec![0.5],
            simple_ub_dist: vec![2.0],
            bound_dist: vec![0.5],
            best_vlb: vec![(-1, HighsImplications_VarBound::default()); 1],
            best_vub: vec![(-1, HighsImplications_VarBound::default()); 1],
            bound_types: vec![HighsTransformedLp_BoundType::default(); 1],
            ..Default::default()
        };
        tlp.vectorsum.set_dimension(1);

        // best_vub[0].1.max_value() > ub + feastol triggers cleanup_vub.
        // With a default (empty) implications set, cleanup_vub marks it redundant.
        tlp.best_vub[0] = (
            0,
            HighsImplications_VarBound {
                coef: 1.0,
                constant: 20.0,
                ..Default::default()
            },
        );

        let mut inds = vec![0];
        let mut vals = vec![1.0];
        let (mut upper, mut solval) = (Vec::new(), Vec::new());
        let mut rhs = 0.0;
        let mut integers_positive = false;
        let _ = tlp.transform(
            &mut ms,
            &lprel,
            &mut vals,
            &mut upper,
            &mut solval,
            &mut inds,
            &mut rhs,
            &mut integers_positive,
            false,
            false,
        );

        assert_eq!(tlp.best_vub[0].0, -1);
        assert_eq!(tlp.ub_dist[0], tlp.simple_ub_dist[0]);
        assert_eq!(tlp.bound_dist[0], cpp_min(tlp.lb_dist[0], tlp.ub_dist[0]));
    }

    /// #3179: with `enforce_same_bds` true and a prior bound type that is not
    /// `kUnused`, the first loop keeps the old type for an integral column whose
    /// simple-bound condition holds. After `init_multi_row_transform` every entry
    /// is `kUnused`.
    #[test]
    fn x1_tlp_enforce_same_bds_keeps_type() {
        let mut ms = HighsMipSolver::default();
        {
            let lp = &mut ms.orig_model_;
            lp.num_col_ = 1;
            lp.num_row_ = 0;
            lp.col_cost_ = vec![0.0; 1];
            lp.col_lower_ = vec![0.0];
            lp.col_upper_ = vec![1.0];
            lp.integrality_ = vec![HighsVarType::kInteger];
            lp.a_matrix_.num_col_ = 1;
            lp.a_matrix_.start_ = vec![0, 0];
        }
        ms.mipdata_.a_rstart_ = vec![0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        let dom = HighsDomain::new(&ms);
        ms.mipdata_.domains = vec![Slot::new(dom)];

        let mut lprel = HighsLpRelaxation::default();
        lprel.lpsolver.model_lp_.num_col_ = 1;

        let mut tlp = HighsTransformedLp {
            globaldom_: 0,
            lb_dist: vec![0.2],
            ub_dist: vec![0.8],
            simple_lb_dist: vec![0.2],
            simple_ub_dist: vec![0.8],
            bound_dist: vec![0.2],
            best_vlb: vec![(-1, HighsImplications_VarBound::default()); 1],
            best_vub: vec![(-1, HighsImplications_VarBound::default()); 1],
            bound_types: vec![HighsTransformedLp_BoundType::kSimpleLb; 1],
            ..Default::default()
        };
        tlp.vectorsum.set_dimension(1);

        // ub - lb == 1.0 <= 1.5, so the integral column hits the early-continue
        // branch in the enforce_same_bds path and keeps kSimpleLb.
        let mut inds = vec![0];
        let mut vals = vec![1.0];
        let (mut upper, mut solval) = (Vec::new(), Vec::new());
        let mut rhs = 0.0;
        let mut integers_positive = false;
        let _ = tlp.transform(
            &mut ms,
            &lprel,
            &mut vals,
            &mut upper,
            &mut solval,
            &mut inds,
            &mut rhs,
            &mut integers_positive,
            false,
            true,
        );

        assert_eq!(tlp.bound_types[0], HighsTransformedLp_BoundType::kSimpleLb);

        tlp.init_multi_row_transform();
        assert_eq!(tlp.bound_types[0], HighsTransformedLp_BoundType::kUnused);
    }

    /// #3223 (`mip_cdouble_proofs`): one variable-lower-bound term with coefficient 1e-10
    /// on a column in [0, 1e6]. Runs `transform` (the cleanup `IsZero` lambda) and returns the
    /// surviving indices and values.
    fn x1b_train_b_transform_vlb(on: bool) -> (Vec<HighsInt>, Vec<f64>) {
        let mut ms = HighsMipSolver::default();
        {
            let lp = &mut ms.orig_model_;
            lp.num_col_ = 2;
            lp.num_row_ = 0;
            lp.col_cost_ = vec![0.0; 2];
            lp.col_lower_ = vec![0.0, 0.0];
            lp.col_upper_ = vec![1e6, 1e6];
            lp.integrality_ = vec![HighsVarType::kContinuous; 2];
            lp.a_matrix_.num_col_ = 2;
            lp.a_matrix_.start_ = vec![0, 0, 0];
        }
        ms.mipdata_.a_rstart_ = vec![0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.mip_cdouble_proofs = on;
        let dom = HighsDomain::new(&ms);
        ms.mipdata_.domains = vec![Slot::new(dom)];

        let mut lprel = HighsLpRelaxation::default();
        lprel.lpsolver.model_lp_.num_col_ = 2;

        let mut tlp = HighsTransformedLp {
            globaldom_: 0,
            lb_dist: vec![0.5, 0.5],
            ub_dist: vec![5.0, 5.0],
            simple_lb_dist: vec![0.5, 0.5],
            simple_ub_dist: vec![5.0, 5.0],
            bound_dist: vec![0.5, 0.5],
            best_vlb: vec![
                (
                    1,
                    HighsImplications_VarBound {
                        coef: 1.0,
                        constant: 0.0,
                        ..Default::default()
                    },
                ),
                (-1, HighsImplications_VarBound::default()),
            ],
            best_vub: vec![(-1, HighsImplications_VarBound::default()); 2],
            bound_types: vec![HighsTransformedLp_BoundType::default(); 2],
            ..Default::default()
        };
        tlp.vectorsum.set_dimension(2);

        let mut inds = vec![0];
        let mut vals = vec![1e-10];
        let (mut upper, mut solval) = (Vec::new(), Vec::new());
        let mut rhs = 0.0;
        let mut integers_positive = false;
        let _ = tlp.transform(
            &mut ms,
            &lprel,
            &mut vals,
            &mut upper,
            &mut solval,
            &mut inds,
            &mut rhs,
            &mut integers_positive,
            false,
            false,
        );
        (inds, vals)
    }

    /// #3223: with `mip_cdouble_proofs` on, `HighsDomain::termIsZero` keeps the 1e-10
    /// coefficient (1e-10 * 1e6 > epsilon), so it survives the cleanup.
    #[test]
    fn x1b_train_b_dual_proof_term_is_zero() {
        let (inds, vals) = x1b_train_b_transform_vlb(true);
        assert_eq!(inds, vec![1]);
        assert_eq!(vals.len(), 1);
        // the transform flips the sign of a variable-lower-bound term; the magnitude is kept
        assert!((vals[0].abs() - 1e-10).abs() < 1e-20, "vals {:?}", vals);
    }

    /// #3223 off path: the 1.15.1 `small_matrix_value` (1e-9) test drops the same 1e-10
    /// coefficient, so the transformed row is empty.
    #[test]
    fn x1b_train_b_off_identity() {
        let (inds, vals) = x1b_train_b_transform_vlb(false);
        assert!(inds.is_empty());
        assert!(vals.is_empty());
    }
}
