//! HiGHS source: mip/HighsDomain.cpp (lines 1259-2118)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

impl HighsDomain {
    /// HiGHS: HighsDomain.cpp::HighsDomain::computeMinActivity (lines 1260-1303)
    /// HiGHS latest (#3184, c5ebac196): HighsDomain.cpp::HighsDomain::computeMinActivity (lines 1276-1303)
    /// C++: `void computeMinActivity(HighsInt start, HighsInt end, const HighsInt * ARindex, const double * ARvalue, HighsInt & ninfmin, HighsCDouble & activitymin) const`
    pub fn compute_min_activity(
        &self,
        start: HighsInt,
        end: HighsInt,
        a_rindex: &[HighsInt],
        a_rvalue: &[f64],
        ninfmin: &mut HighsInt,
        activitymin: &mut HighsCDouble,
    ) {
        // HiGHS #3184 (mip_cdouble_act_contrib): HighsDomain.cpp lines 1276-1303 at c5ebac196
        if self.infeasible_ {
            *activitymin = HighsCDouble::from(0.0);
            *ninfmin = 0;
            for j in start..end {
                let col = a_rindex[j as usize];
                let val = a_rvalue[j as usize];

                let mut tmp: HighsInt = 0;
                let lb = self.get_col_lower_pos(col, self.infeasible_pos - 1, &mut tmp);
                let ub = self.get_col_upper_pos(col, self.infeasible_pos - 1, &mut tmp);

                if self.x1_cdouble_act_contrib {
                    let contributionmin = activity_contribution_min_cd(val, lb, ub);

                    if contributionmin == -kHighsInf {
                        *ninfmin += 1;
                    } else {
                        *activitymin += contributionmin;
                    }
                } else {
                    let contributionmin = activity_contribution_min(val, lb, ub);

                    if contributionmin == -kHighsInf {
                        *ninfmin += 1;
                    } else {
                        *activitymin += contributionmin;
                    }
                }
            }
        } else {
            *activitymin = HighsCDouble::from(0.0);
            *ninfmin = 0;
            for j in start..end {
                let col = a_rindex[j as usize];
                let val = a_rvalue[j as usize];

                if self.x1_cdouble_act_contrib {
                    let contributionmin = activity_contribution_min_cd(
                        val,
                        self.col_lower_[col as usize],
                        self.col_upper_[col as usize],
                    );

                    if contributionmin == -kHighsInf {
                        *ninfmin += 1;
                    } else {
                        *activitymin += contributionmin;
                    }
                } else {
                    let contributionmin = activity_contribution_min(
                        val,
                        self.col_lower_[col as usize],
                        self.col_upper_[col as usize],
                    );

                    if contributionmin == -kHighsInf {
                        *ninfmin += 1;
                    } else {
                        *activitymin += contributionmin;
                    }
                }
            }
        }
        activitymin.renormalize();
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::computeMaxActivity (lines 1305-1347)
    /// HiGHS latest (#3184, c5ebac196): HighsDomain.cpp::HighsDomain::computeMaxActivity (lines 1321-1347)
    /// C++: `void computeMaxActivity(HighsInt start, HighsInt end, const HighsInt * ARindex, const double * ARvalue, HighsInt & ninfmax, HighsCDouble & activitymax) const`
    pub fn compute_max_activity(
        &self,
        start: HighsInt,
        end: HighsInt,
        a_rindex: &[HighsInt],
        a_rvalue: &[f64],
        ninfmax: &mut HighsInt,
        activitymax: &mut HighsCDouble,
    ) {
        // HiGHS #3184 (mip_cdouble_act_contrib): HighsDomain.cpp lines 1321-1347 at c5ebac196
        if self.infeasible_ {
            *activitymax = HighsCDouble::from(0.0);
            *ninfmax = 0;
            for j in start..end {
                let col = a_rindex[j as usize];
                let val = a_rvalue[j as usize];

                let mut tmp: HighsInt = 0;
                let lb = self.get_col_lower_pos(col, self.infeasible_pos - 1, &mut tmp);
                let ub = self.get_col_upper_pos(col, self.infeasible_pos - 1, &mut tmp);

                if self.x1_cdouble_act_contrib {
                    let contributionmin = activity_contribution_max_cd(val, lb, ub);

                    if contributionmin == kHighsInf {
                        *ninfmax += 1;
                    } else {
                        *activitymax += contributionmin;
                    }
                } else {
                    let contributionmin = activity_contribution_max(val, lb, ub);

                    if contributionmin == kHighsInf {
                        *ninfmax += 1;
                    } else {
                        *activitymax += contributionmin;
                    }
                }
            }
        } else {
            *activitymax = HighsCDouble::from(0.0);
            *ninfmax = 0;
            for j in start..end {
                let col = a_rindex[j as usize];
                let val = a_rvalue[j as usize];

                if self.x1_cdouble_act_contrib {
                    let contributionmin = activity_contribution_max_cd(
                        val,
                        self.col_lower_[col as usize],
                        self.col_upper_[col as usize],
                    );

                    if contributionmin == kHighsInf {
                        *ninfmax += 1;
                    } else {
                        *activitymax += contributionmin;
                    }
                } else {
                    let contributionmin = activity_contribution_max(
                        val,
                        self.col_lower_[col as usize],
                        self.col_upper_[col as usize],
                    );

                    if contributionmin == kHighsInf {
                        *ninfmax += 1;
                    } else {
                        *activitymax += contributionmin;
                    }
                }
            }
        }
        activitymax.renormalize();
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::adjustedUb (lines 1349-1379)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::adjustedUb (lines 1608-1636)
    /// C++: `double adjustedUb(HighsInt col, HighsCDouble boundVal, bool & accept) const`
    pub fn adjusted_ub(
        &self,
        ms: &HighsMipSolver,
        col: HighsInt,
        bound_val: HighsCDouble,
        accept: &mut bool,
    ) -> f64 {
        let cu = col as usize;
        let bound: f64;

        if ms.is_col_integral(col) {
            bound = f64::from((bound_val + ms.mipdata_.feastol).floor());
            *accept = bound < self.col_upper_[cu]
                && self.col_upper_[cu] - bound > 1000.0 * ms.mipdata_.feastol * bound.abs();
        } else {
            // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1617 (1.16.0)
            let snap = if ms.options_mip_.h16_cdouble {
                (bound_val - self.col_lower_[cu]).abs() <= ms.mipdata_.epsilon
            } else {
                (f64::from(bound_val) - self.col_lower_[cu]).abs() <= ms.mipdata_.epsilon
            };
            if snap {
                bound = self.col_lower_[cu];
            } else {
                bound = f64::from(bound_val);
            }
            if self.col_upper_[cu] == kHighsInf {
                *accept = true;
            } else if bound + 1000.0 * ms.mipdata_.feastol < self.col_upper_[cu] {
                let mut relative_improve = self.col_upper_[cu] - bound;
                if self.col_lower_[cu] != -kHighsInf {
                    relative_improve /= self.col_upper_[cu] - self.col_lower_[cu];
                } else {
                    relative_improve /= cpp_max(self.col_upper_[cu].abs(), bound.abs());
                }
                *accept = relative_improve
                    >= self
                        .cont_min_improve
                        .min(ms.options_mip_.mip_cont_prop_improve);
            } else {
                *accept = false;
            }
        }

        bound
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::adjustedLb (lines 1381-1411)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::adjustedLb (lines 1638-1666)
    /// C++: `double adjustedLb(HighsInt col, HighsCDouble boundVal, bool & accept) const`
    pub fn adjusted_lb(
        &self,
        ms: &HighsMipSolver,
        col: HighsInt,
        bound_val: HighsCDouble,
        accept: &mut bool,
    ) -> f64 {
        let cu = col as usize;
        let bound: f64;

        if ms.is_col_integral(col) {
            bound = f64::from((bound_val - ms.mipdata_.feastol).ceil());
            *accept = bound > self.col_lower_[cu]
                && bound - self.col_lower_[cu] > 1000.0 * ms.mipdata_.feastol * bound.abs();
        } else {
            // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1647 (1.16.0)
            let snap = if ms.options_mip_.h16_cdouble {
                (self.col_upper_[cu] - bound_val).abs() <= ms.mipdata_.epsilon
            } else {
                (self.col_upper_[cu] - f64::from(bound_val)).abs() <= ms.mipdata_.epsilon
            };
            if snap {
                bound = self.col_upper_[cu];
            } else {
                bound = f64::from(bound_val);
            }
            if self.col_lower_[cu] == -kHighsInf {
                *accept = true;
            } else if bound - 1000.0 * ms.mipdata_.feastol > self.col_lower_[cu] {
                let mut relative_improve = bound - self.col_lower_[cu];
                if self.col_upper_[cu] != kHighsInf {
                    relative_improve /= self.col_upper_[cu] - self.col_lower_[cu];
                } else {
                    relative_improve /= cpp_max(self.col_lower_[cu].abs(), bound.abs());
                }
                *accept = relative_improve
                    >= self
                        .cont_min_improve
                        .min(ms.options_mip_.mip_cont_prop_improve);
            } else {
                *accept = false;
            }
        }

        bound
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::propagateRowUpper (lines 1413-1455)
    /// HiGHS latest (#3184, c5ebac196): HighsDomain.cpp::HighsDomain::propagateRowUpper (lines 1420-1455)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::propagateRowUpper (lines 1668-1709)
    /// C++: `HighsInt propagateRowUpper(const HighsInt * Rindex, const double * Rvalue, HighsInt Rlen, double Rupper, const HighsCDouble & minactivity, HighsInt ninfmin, HighsDomainChange * boundchgs) const`
    // PORT NOTE: signature changed: `HighsDomainChange* boundchgs` (room for `Rlen` entries
    // starting at the row's offset in `changedbounds`) is `&mut changedbounds[start..]`;
    // `ms` is only read (`mipdata_->feastol`, `isColIntegral`).
    pub fn propagate_row_upper(
        &self,
        ms: &HighsMipSolver,
        rindex: &[HighsInt],
        rvalue: &[f64],
        rlen: HighsInt,
        rupper: f64,
        minactivity: &HighsCDouble,
        ninfmin: HighsInt,
        boundchgs: &mut [HighsDomainChange],
    ) -> HighsInt {
        if ninfmin > 1 {
            return 0;
        }
        let mut numchgs: HighsInt = 0;
        for i in 0..rlen as usize {
            // HiGHS #3184 (mip_cdouble_act_contrib): HighsDomain.cpp line 1423 at c5ebac196
            let minresact: HighsCDouble = if self.x1_cdouble_act_contrib {
                let actcontribution = activity_contribution_min_cd(
                    rvalue[i],
                    self.col_lower_[rindex[i] as usize],
                    self.col_upper_[rindex[i] as usize],
                );
                if ninfmin == 1 {
                    if actcontribution != -kHighsInf {
                        continue;
                    }

                    *minactivity
                } else {
                    *minactivity - actcontribution
                }
            } else {
                let actcontribution = activity_contribution_min(
                    rvalue[i],
                    self.col_lower_[rindex[i] as usize],
                    self.col_upper_[rindex[i] as usize],
                );
                if ninfmin == 1 {
                    if actcontribution != -kHighsInf {
                        continue;
                    }

                    *minactivity
                } else {
                    *minactivity - actcontribution
                }
            };

            let bound_val: HighsCDouble = (rupper - minresact) / rvalue[i];
            // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1690 (1.16.0)
            let skip = if ms.options_mip_.h16_cdouble {
                bound_val.abs() * kHighsTiny > ms.mipdata_.feastol
            } else {
                (f64::from(bound_val) * kHighsTiny).abs() > ms.mipdata_.feastol
            };
            if skip {
                continue;
            }

            if rvalue[i] > 0.0 {
                let mut accept = false;

                let bound = self.adjusted_ub(ms, rindex[i], bound_val, &mut accept);
                if accept {
                    boundchgs[numchgs as usize] = HighsDomainChange {
                        boundval: bound,
                        column: rindex[i],
                        boundtype: HighsBoundType::kUpper,
                    };
                    numchgs += 1;
                }
            } else {
                let mut accept = false;

                let bound = self.adjusted_lb(ms, rindex[i], bound_val, &mut accept);
                if accept {
                    boundchgs[numchgs as usize] = HighsDomainChange {
                        boundval: bound,
                        column: rindex[i],
                        boundtype: HighsBoundType::kLower,
                    };
                    numchgs += 1;
                }
            }
        }

        numchgs
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::propagateRowLower (lines 1457-1498)
    /// HiGHS latest (#3184, c5ebac196): HighsDomain.cpp::HighsDomain::propagateRowLower (lines 1463-1498)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::propagateRowLower (lines 1711-1751)
    /// C++: `HighsInt propagateRowLower(const HighsInt * Rindex, const double * Rvalue, HighsInt Rlen, double Rlower, const HighsCDouble & maxactivity, HighsInt ninfmax, HighsDomainChange * boundchgs) const`
    // PORT NOTE: as `propagate_row_upper`.
    pub fn propagate_row_lower(
        &self,
        ms: &HighsMipSolver,
        rindex: &[HighsInt],
        rvalue: &[f64],
        rlen: HighsInt,
        rlower: f64,
        maxactivity: &HighsCDouble,
        ninfmax: HighsInt,
        boundchgs: &mut [HighsDomainChange],
    ) -> HighsInt {
        if ninfmax > 1 {
            return 0;
        }
        let mut numchgs: HighsInt = 0;
        for i in 0..rlen as usize {
            // HiGHS #3184 (mip_cdouble_act_contrib): HighsDomain.cpp line 1466 at c5ebac196
            let maxresact: HighsCDouble = if self.x1_cdouble_act_contrib {
                let actcontribution = activity_contribution_max_cd(
                    rvalue[i],
                    self.col_lower_[rindex[i] as usize],
                    self.col_upper_[rindex[i] as usize],
                );
                if ninfmax == 1 {
                    if actcontribution != kHighsInf {
                        continue;
                    }

                    *maxactivity
                } else {
                    *maxactivity - actcontribution
                }
            } else {
                let actcontribution = activity_contribution_max(
                    rvalue[i],
                    self.col_lower_[rindex[i] as usize],
                    self.col_upper_[rindex[i] as usize],
                );
                if ninfmax == 1 {
                    if actcontribution != kHighsInf {
                        continue;
                    }

                    *maxactivity
                } else {
                    *maxactivity - actcontribution
                }
            };

            let bound_val: HighsCDouble = (rlower - maxresact) / rvalue[i];
            // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1733 (1.16.0)
            let skip = if ms.options_mip_.h16_cdouble {
                bound_val.abs() * kHighsTiny > ms.mipdata_.feastol
            } else {
                (f64::from(bound_val) * kHighsTiny).abs() > ms.mipdata_.feastol
            };
            if skip {
                continue;
            }

            if rvalue[i] < 0.0 {
                let mut accept = false;

                let bound = self.adjusted_ub(ms, rindex[i], bound_val, &mut accept);
                if accept {
                    boundchgs[numchgs as usize] = HighsDomainChange {
                        boundval: bound,
                        column: rindex[i],
                        boundtype: HighsBoundType::kUpper,
                    };
                    numchgs += 1;
                }
            } else {
                let mut accept = false;

                let bound = self.adjusted_lb(ms, rindex[i], bound_val, &mut accept);
                if accept {
                    boundchgs[numchgs as usize] = HighsDomainChange {
                        boundval: bound,
                        column: rindex[i],
                        boundtype: HighsBoundType::kLower,
                    };
                    numchgs += 1;
                }
            }
        }

        numchgs
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::updateThresholdLbChange (lines 1500-1513)
    /// HiGHS 1.16 (#3525): HighsDomain.cpp::updateThresholdLbChange (lines 1753-1764)
    /// C++: `void updateThresholdLbChange(HighsInt col, double newbound, double val, double & threshold) const`
    pub fn update_threshold_lb_change(
        &self,
        ms: &HighsMipSolver,
        col: HighsInt,
        newbound: f64,
        val: f64,
        threshold: &mut f64,
    ) {
        let cu = col as usize;
        if newbound != self.col_upper_[cu] {
            // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp line 1758 (1.16.0)
            let threshold_new = val.abs()
                * bound_range(
                    self.col_upper_[cu],
                    newbound,
                    ms.mipdata_.feastol,
                    self.variable_type(ms, col),
                    ms.options_mip_.h16_bound_range_inf,
                    ms.options_mip_.mip_cont_prop_improve,
                );

            // the new threshold is now the maximum of the new threshold and the current
            // one
            *threshold = cpp_max(cpp_max(*threshold, threshold_new), ms.mipdata_.feastol);
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::updateThresholdUbChange (lines 1515-1528)
    /// HiGHS 1.16 (#3525): HighsDomain.cpp::updateThresholdUbChange (lines 1766-1777)
    /// C++: `void updateThresholdUbChange(HighsInt col, double newbound, double val, double & threshold) const`
    pub fn update_threshold_ub_change(
        &self,
        ms: &HighsMipSolver,
        col: HighsInt,
        newbound: f64,
        val: f64,
        threshold: &mut f64,
    ) {
        let cu = col as usize;
        if newbound != self.col_lower_[cu] {
            // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp line 1771 (1.16.0)
            let threshold_new = val.abs()
                * bound_range(
                    newbound,
                    self.col_lower_[cu],
                    ms.mipdata_.feastol,
                    self.variable_type(ms, col),
                    ms.options_mip_.h16_bound_range_inf,
                    ms.options_mip_.mip_cont_prop_improve,
                );

            // the new threshold is now the maximum of the new threshold and the current
            // one
            *threshold = cpp_max(cpp_max(*threshold, threshold_new), ms.mipdata_.feastol);
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::updateActivityLbChange (lines 1530-1695)
    /// C++: `void updateActivityLbChange(HighsInt col, double oldbound, double newbound)`
    // PORT NOTE: signature changed (`ms`, read only). The cut-pool loop (C++ 1650-1684, with
    // its forward/backward infeasibility unwinding) runs on
    // `let mut props = std::mem::take(&mut self.cutpoolpropagation)`, calling
    // `props[i].update_activity_lb_change(ms, self, ...)`, then puts `props` back; the conflict
    // loop does the same with `conflict_pool_propagation` (`update_activity_lb_change(self, ..)`).
    // Put the vector back before any early return. `objProp_.updateActivityLbChange` is
    // `self.obj_prop_update_activity_lb_change(ms, ..)` (guarded by `self.obj_prop_.is_active()`).
    pub fn update_activity_lb_change(
        &mut self,
        ms: &HighsMipSolver,
        col: HighsInt,
        oldbound: f64,
        newbound: f64,
    ) {
        let (mut oldbound, mut newbound) = (oldbound, newbound);
        let mip = crate::mip_model!(ms);
        let colu = col as usize;
        let start = mip.a_matrix_.start_[colu];
        let mut end = mip.a_matrix_.start_[colu + 1];
        crate::work::add_propagation_units((end - start) as u64); // reporting-only

        if self.obj_prop_.is_active() {
            self.obj_prop_update_activity_lb_change(ms, col, oldbound, newbound);
            if self.infeasible_ {
                return;
            }
        }

        // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp line 1783 at 8a5f9c1d8
        let track_redundancy = newbound > oldbound && self.dual_fix_probing_active;

        // PORT NOTE (perf): the column slice and the row-bound vectors
        // are bound once, as the C++ loop keeps their pointers in registers; indexing through
        // `mip` reloads every vector header after each store to `self`.
        let a_index = &mip.a_matrix_.index_[start as usize..end as usize];
        let a_value = &mip.a_matrix_.value_[start as usize..end as usize];
        let row_lower = &mip.row_lower_[..];
        let row_upper = &mip.row_upper_[..];
        let feastol = ms.mipdata_.feastol;
        for (k, (&row, &val)) in a_index.iter().zip(a_value.iter()).enumerate() {
            let i = start + k as HighsInt;
            let rowu = row as usize;
            if val > 0.0 {
                let deltamin = compute_delta(
                    val,
                    oldbound,
                    newbound,
                    -kHighsInf,
                    &mut self.activitymininf_[rowu],
                );
                self.activitymin_[rowu] += deltamin;

                // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 1806-1810 at 8a5f9c1d8
                let record_redundant = if self.dual_fix_probing_active {
                    self.record_redundant_rows_
                        && !self.dual_fix_probing_.is_zero_cost_fixing_active()
                        && row_lower[rowu] != -kHighsInf
                        && row_upper[rowu] == kHighsInf
                } else {
                    self.record_redundant_rows_
                        && row_lower[rowu] != -kHighsInf
                        && row_upper[rowu] == kHighsInf
                };
                if record_redundant {
                    self.update_redundant_rows(ms, row);
                }

                // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 1812-1813 at 8a5f9c1d8
                if track_redundancy {
                    self.dfp_update_lhs_redundant(ms, row);
                }

                if deltamin <= 0.0 {
                    let mut threshold = self.capacity_threshold_[rowu];
                    self.update_threshold_lb_change(ms, col, newbound, val, &mut threshold);
                    self.capacity_threshold_[rowu] = threshold;
                    continue;
                }

                if row_upper[rowu] != kHighsInf
                    && self.activitymininf_[rowu] == 0
                    && self.activitymin_[rowu] - row_upper[rowu] > feastol
                {
                    ms.mipdata_.debug_solution.node_pruned(&*self);
                    self.infeasible_ = true;
                    self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                    self.infeasible_reason = HighsDomain_Reason::model_row_upper(row);
                    end = i + 1;
                    break;
                }

                if self.activitymininf_[rowu] <= 1
                    && self.propagateflags_[rowu] == 0
                    && row_upper[rowu] != kHighsInf
                {
                    self.mark_propagate(ms, row);
                }
            } else {
                let deltamax = compute_delta(
                    val,
                    oldbound,
                    newbound,
                    -kHighsInf,
                    &mut self.activitymaxinf_[rowu],
                );
                self.activitymax_[rowu] += deltamax;

                // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 1858-1862 at 8a5f9c1d8
                let record_redundant = if self.dual_fix_probing_active {
                    self.record_redundant_rows_
                        && !self.dual_fix_probing_.is_zero_cost_fixing_active()
                        && row_lower[rowu] == -kHighsInf
                        && row_upper[rowu] != kHighsInf
                } else {
                    self.record_redundant_rows_
                        && row_lower[rowu] == -kHighsInf
                        && row_upper[rowu] != kHighsInf
                };
                if record_redundant {
                    self.update_redundant_rows(ms, row);
                }

                // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 1864-1865 at 8a5f9c1d8
                if track_redundancy {
                    self.dfp_update_rhs_redundant(ms, row);
                }

                if deltamax >= 0.0 {
                    let mut threshold = self.capacity_threshold_[rowu];
                    self.update_threshold_lb_change(ms, col, newbound, val, &mut threshold);
                    self.capacity_threshold_[rowu] = threshold;
                    continue;
                }

                if row_lower[rowu] != -kHighsInf
                    && self.activitymaxinf_[rowu] == 0
                    && row_lower[rowu] - self.activitymax_[rowu] > feastol
                {
                    ms.mipdata_.debug_solution.node_pruned(&*self);
                    self.infeasible_ = true;
                    self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                    self.infeasible_reason = HighsDomain_Reason::model_row_lower(row);
                    end = i + 1;
                    break;
                }

                if self.activitymaxinf_[rowu] <= 1
                    && self.propagateflags_[rowu] == 0
                    && row_lower[rowu] != -kHighsInf
                {
                    self.mark_propagate(ms, row);
                }
            }
        }

        if self.infeasible_ {
            // C++ debug asserts on `infeasible_reason` (type kModelRowLower / kModelRowUpper and
            // index `a_matrix_.index_[end - 1]`) are not translated (see PORT NOTE above).
        } else {
            // Explanation: If cutpoolpropagation[i] returns infeasible,
            // we still need to update the domain threshold values using
            // cutpoolpropagation[j], j > i, and reverse the activity
            // changes made by cutpoolpropagation[j] j < i.
            let mut infeascutpool: HighsInt = -1;
            let ncutpoolprop = self.cutpoolpropagation.len() as HighsInt;
            let mut props = std::mem::take(&mut self.cutpoolpropagation);
            for i in 0..ncutpoolprop {
                let iu = i as usize;
                if !self.infeasible_ {
                    let infeas = self.infeasible_;
                    props[iu].update_activity_lb_change(
                        ms, self, col, oldbound, newbound, true, true, infeas,
                    );
                    if self.infeasible_ {
                        infeascutpool = i;
                    }
                } else {
                    let infeas = self.infeasible_;
                    props[iu].update_activity_lb_change(
                        ms, self, col, oldbound, newbound, true, false, infeas,
                    );
                }
            }
            for i in 0..infeascutpool {
                let infeas = self.infeasible_;
                props[i as usize].update_activity_lb_change(
                    ms, self, col, oldbound, newbound, false, true, infeas,
                );
            }
            self.cutpoolpropagation = props;
        }

        if self.infeasible_ {
            std::mem::swap(&mut oldbound, &mut newbound);
            for i in start..end {
                let ii = i as usize;
                let rowu = mip.a_matrix_.index_[ii] as usize;
                let val = mip.a_matrix_.value_[ii];
                if val > 0.0 {
                    let d = compute_delta(
                        val,
                        oldbound,
                        newbound,
                        -kHighsInf,
                        &mut self.activitymininf_[rowu],
                    );
                    self.activitymin_[rowu] += d;
                } else {
                    let d = compute_delta(
                        val,
                        oldbound,
                        newbound,
                        -kHighsInf,
                        &mut self.activitymaxinf_[rowu],
                    );
                    self.activitymax_[rowu] += d;
                }
            }

            if self.obj_prop_.is_active() {
                self.obj_prop_update_activity_lb_change(ms, col, oldbound, newbound);
            }
        } else {
            let mut conflict_props = std::mem::take(&mut self.conflict_pool_propagation);
            for conflictprop in conflict_props.iter_mut() {
                conflictprop.update_activity_lb_change(&*self, col, oldbound, newbound);
            }
            self.conflict_pool_propagation = conflict_props;
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::updateActivityUbChange (lines 1697-1868)
    /// C++: `void updateActivityUbChange(HighsInt col, double oldbound, double newbound)`
    // PORT NOTE: as `update_activity_lb_change`.
    pub fn update_activity_ub_change(
        &mut self,
        ms: &HighsMipSolver,
        col: HighsInt,
        oldbound: f64,
        newbound: f64,
    ) {
        let (mut oldbound, mut newbound) = (oldbound, newbound);
        let mip = crate::mip_model!(ms);
        let colu = col as usize;
        let start = mip.a_matrix_.start_[colu];
        let mut end = mip.a_matrix_.start_[colu + 1];
        crate::work::add_propagation_units((end - start) as u64); // reporting-only

        if self.obj_prop_.is_active() {
            self.obj_prop_update_activity_ub_change(ms, col, oldbound, newbound);
            if self.infeasible_ {
                return;
            }
        }

        // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp line 1958 at 8a5f9c1d8
        let track_redundancy = newbound < oldbound && self.dual_fix_probing_active;

        // PORT NOTE (perf): the column slice and the row-bound vectors
        // are bound once, as the C++ loop keeps their pointers in registers; indexing through
        // `mip` reloads every vector header after each store to `self`.
        let a_index = &mip.a_matrix_.index_[start as usize..end as usize];
        let a_value = &mip.a_matrix_.value_[start as usize..end as usize];
        let row_lower = &mip.row_lower_[..];
        let row_upper = &mip.row_upper_[..];
        let feastol = ms.mipdata_.feastol;
        for (k, (&row, &val)) in a_index.iter().zip(a_value.iter()).enumerate() {
            let i = start + k as HighsInt;
            let rowu = row as usize;
            if val > 0.0 {
                let deltamax = compute_delta(
                    val,
                    oldbound,
                    newbound,
                    kHighsInf,
                    &mut self.activitymaxinf_[rowu],
                );
                self.activitymax_[rowu] += deltamax;

                // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 1981-1985 at 8a5f9c1d8
                let record_redundant = if self.dual_fix_probing_active {
                    self.record_redundant_rows_
                        && !self.dual_fix_probing_.is_zero_cost_fixing_active()
                        && row_lower[rowu] == -kHighsInf
                        && row_upper[rowu] != kHighsInf
                } else {
                    self.record_redundant_rows_
                        && row_lower[rowu] == -kHighsInf
                        && row_upper[rowu] != kHighsInf
                };
                if record_redundant {
                    self.update_redundant_rows(ms, row);
                }

                // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 1987-1988 at 8a5f9c1d8
                if track_redundancy {
                    self.dfp_update_rhs_redundant(ms, row);
                }

                if deltamax >= 0.0 {
                    let mut threshold = self.capacity_threshold_[rowu];
                    self.update_threshold_ub_change(ms, col, newbound, val, &mut threshold);
                    self.capacity_threshold_[rowu] = threshold;
                    continue;
                }

                if row_lower[rowu] != -kHighsInf
                    && self.activitymaxinf_[rowu] == 0
                    && row_lower[rowu] - self.activitymax_[rowu] > feastol
                {
                    ms.mipdata_.debug_solution.node_pruned(&*self);
                    self.infeasible_ = true;
                    self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                    self.infeasible_reason = HighsDomain_Reason::model_row_lower(row);
                    end = i + 1;
                    break;
                }

                if self.activitymaxinf_[rowu] <= 1
                    && self.propagateflags_[rowu] == 0
                    && row_lower[rowu] != -kHighsInf
                {
                    self.mark_propagate(ms, row);
                }
            } else {
                let deltamin = compute_delta(
                    val,
                    oldbound,
                    newbound,
                    kHighsInf,
                    &mut self.activitymininf_[rowu],
                );
                self.activitymin_[rowu] += deltamin;

                // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 2036-2040 at 8a5f9c1d8
                let record_redundant = if self.dual_fix_probing_active {
                    self.record_redundant_rows_
                        && !self.dual_fix_probing_.is_zero_cost_fixing_active()
                        && row_lower[rowu] != -kHighsInf
                        && row_upper[rowu] == kHighsInf
                } else {
                    self.record_redundant_rows_
                        && row_lower[rowu] != -kHighsInf
                        && row_upper[rowu] == kHighsInf
                };
                if record_redundant {
                    self.update_redundant_rows(ms, row);
                }

                // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 2042-2043 at 8a5f9c1d8
                if track_redundancy {
                    self.dfp_update_lhs_redundant(ms, row);
                }

                if deltamin <= 0.0 {
                    let mut threshold = self.capacity_threshold_[rowu];
                    self.update_threshold_ub_change(ms, col, newbound, val, &mut threshold);
                    self.capacity_threshold_[rowu] = threshold;
                    continue;
                }

                if row_upper[rowu] != kHighsInf
                    && self.activitymininf_[rowu] == 0
                    && self.activitymin_[rowu] - row_upper[rowu] > feastol
                {
                    ms.mipdata_.debug_solution.node_pruned(&*self);
                    self.infeasible_ = true;
                    self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                    self.infeasible_reason = HighsDomain_Reason::model_row_upper(row);
                    end = i + 1;
                    break;
                }

                if self.activitymininf_[rowu] <= 1
                    && self.propagateflags_[rowu] == 0
                    && row_upper[rowu] != kHighsInf
                {
                    self.mark_propagate(ms, row);
                }
            }
        }

        if self.infeasible_ {
            // C++ debug asserts on `infeasible_reason` are not translated (see
            // `update_activity_lb_change`).
        } else {
            // Explanation: as in `update_activity_lb_change`.
            let mut infeascutpool: HighsInt = -1;
            let ncutpoolprop = self.cutpoolpropagation.len() as HighsInt;
            let mut props = std::mem::take(&mut self.cutpoolpropagation);
            for i in 0..ncutpoolprop {
                let iu = i as usize;
                if !self.infeasible_ {
                    let infeas = self.infeasible_;
                    props[iu].update_activity_ub_change(
                        ms, self, col, oldbound, newbound, true, true, infeas,
                    );
                    if self.infeasible_ {
                        infeascutpool = i;
                    }
                } else {
                    let infeas = self.infeasible_;
                    props[iu].update_activity_ub_change(
                        ms, self, col, oldbound, newbound, true, false, infeas,
                    );
                }
            }
            for i in 0..infeascutpool {
                let infeas = self.infeasible_;
                props[i as usize].update_activity_ub_change(
                    ms, self, col, oldbound, newbound, false, true, infeas,
                );
            }
            self.cutpoolpropagation = props;
        }

        if self.infeasible_ {
            std::mem::swap(&mut oldbound, &mut newbound);
            for i in start..end {
                let ii = i as usize;
                let rowu = mip.a_matrix_.index_[ii] as usize;
                let val = mip.a_matrix_.value_[ii];
                if val > 0.0 {
                    let d = compute_delta(
                        val,
                        oldbound,
                        newbound,
                        kHighsInf,
                        &mut self.activitymaxinf_[rowu],
                    );
                    self.activitymax_[rowu] += d;
                } else {
                    let d = compute_delta(
                        val,
                        oldbound,
                        newbound,
                        kHighsInf,
                        &mut self.activitymininf_[rowu],
                    );
                    self.activitymin_[rowu] += d;
                }
            }

            if self.obj_prop_.is_active() {
                self.obj_prop_update_activity_ub_change(ms, col, oldbound, newbound);
            }
        } else {
            let mut conflict_props = std::mem::take(&mut self.conflict_pool_propagation);
            for conflictprop in conflict_props.iter_mut() {
                conflictprop.update_activity_ub_change(&*self, col, oldbound, newbound);
            }
            self.conflict_pool_propagation = conflict_props;
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::recomputeCapacityThreshold (lines 1870-1887)
    /// HiGHS 1.16 (#3525): HighsDomain.cpp::recomputeCapacityThreshold (lines 2119-2136)
    /// C++: `void recomputeCapacityThreshold(HighsInt row)`
    pub fn recompute_capacity_threshold(&mut self, ms: &HighsMipSolver, row: HighsInt) {
        let rowu = row as usize;
        let start = ms.mipdata_.a_rstart_[rowu];
        let end = ms.mipdata_.a_rstart_[rowu + 1];

        self.capacity_threshold_[rowu] = -self.feastol(ms);
        for i in start..end {
            let iu = i as usize;
            let col = ms.mipdata_.a_rindex_[iu];
            let cu = col as usize;

            if self.col_upper_[cu] == self.col_lower_[cu] {
                continue;
            }

            // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp line 2130 (1.16.0)
            let threshold = ms.mipdata_.a_rvalue_[iu].abs()
                * bound_range(
                    self.col_upper_[cu],
                    self.col_lower_[cu],
                    self.feastol(ms),
                    self.variable_type(ms, col),
                    ms.options_mip_.h16_bound_range_inf,
                    ms.options_mip_.mip_cont_prop_improve,
                );

            self.capacity_threshold_[rowu] = cpp_max(
                cpp_max(self.capacity_threshold_[rowu], threshold),
                self.feastol(ms),
            );
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::updateRedundantRows (lines 1889-1897)
    /// C++: `void updateRedundantRows(HighsInt row)`
    pub fn update_redundant_rows(&mut self, ms: &HighsMipSolver, row: HighsInt) {
        if !self.is_redundant_row(ms, row) {
            // row that was found to be redundant should not be non-redundant
            // (C++ debug assert: `redundantRows_.find(row) == nullptr`)
            return;
        }
        // row is redundant
        self.redundant_rows_
            .insert(HighsHashTableEntry::new(row, ()));
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::getRedundantRowValue (lines 1899-1907)
    /// C++: `double getRedundantRowValue(HighsInt row) const`
    pub fn get_redundant_row_value(&self, ms: &HighsMipSolver, row: HighsInt) -> f64 {
        let rowu = row as usize;
        if ms.row_lower_2(row) != -kHighsInf {
            f64::from(self.activitymin_[rowu] - ms.row_lower_2(row))
        } else {
            f64::from(self.activitymax_[rowu] - ms.row_upper_2(row))
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::isRedundantRow (lines 1909-1912)
    /// HiGHS 1.16 (#3575): HighsDomain.cpp::isRedundantRow (line 2159)
    /// C++: `bool isRedundantRow(HighsInt row) const`
    pub fn is_redundant_row(&self, ms: &HighsMipSolver, row: HighsInt) -> bool {
        // HiGHS 1.16 #3575 (h16_row_redundant): HighsDomain.cpp line 2159 (1.16.0)
        let h16_row_redundant = ms.options_mip_.h16_row_redundant;
        if h16_row_redundant {
            self.is_row_lower_redundant(ms, row) && self.is_row_upper_redundant(ms, row)
        } else {
            self.get_min_activity(row) >= ms.row_lower_2(row) - self.feastol(ms)
                && self.get_max_activity(row) <= ms.row_upper_2(row) + self.feastol(ms)
        }
    }

    /// HiGHS 1.16 (#3575): HighsDomain.h::isRowLowerRedundant (lines 817-821)
    /// C++: `bool isRowLowerRedundant(HighsInt row) const`
    #[inline]
    pub fn is_row_lower_redundant(&self, ms: &HighsMipSolver, row: HighsInt) -> bool {
        let rowu = row as usize;
        ms.row_lower_2(row) == -kHighsInf
            || (self.activitymininf_[rowu] == 0
                && ms.row_lower_2(row) - self.activitymin_[rowu] <= self.feastol(ms))
    }

    /// HiGHS 1.16 (#3575): HighsDomain.h::isRowUpperRedundant (lines 823-827)
    /// C++: `bool isRowUpperRedundant(HighsInt row) const`
    #[inline]
    pub fn is_row_upper_redundant(&self, ms: &HighsMipSolver, row: HighsInt) -> bool {
        let rowu = row as usize;
        ms.row_upper_2(row) == kHighsInf
            || (self.activitymaxinf_[rowu] == 0
                && self.activitymax_[rowu] - ms.row_upper_2(row) <= self.feastol(ms))
    }

    /// HiGHS 1.16 (#3575): HighsDomain.h::isRowLowerInfeasible (lines 829-833)
    /// C++: `bool isRowLowerInfeasible(HighsInt row) const`
    pub fn is_row_lower_infeasible(&self, ms: &HighsMipSolver, row: HighsInt) -> bool {
        let rowu = row as usize;
        ms.row_lower_2(row) != -kHighsInf
            && self.activitymaxinf_[rowu] == 0
            && ms.row_lower_2(row) - self.activitymax_[rowu] > self.feastol(ms)
    }

    /// HiGHS 1.16 (#3575): HighsDomain.h::isRowUpperInfeasible (lines 835-838)
    /// C++: `bool isRowUpperInfeasible(HighsInt row) const`
    pub fn is_row_upper_infeasible(&self, ms: &HighsMipSolver, row: HighsInt) -> bool {
        let rowu = row as usize;
        ms.row_upper_2(row) != kHighsInf
            && self.activitymininf_[rowu] == 0
            && self.activitymin_[rowu] - ms.row_upper_2(row) > self.feastol(ms)
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::markPropagateCut (lines 1914-1934)
    /// C++: `void markPropagateCut(Reason reason)`
    // PORT NOTE: `cutpoolpropagation[t].markPropagateCut(i)` is
    // `self.cutpoolpropagation[t].mark_propagate_cut(rhs, i)` with
    // `rhs = ms.mipdata_.cutpools[self.cutpoolpropagation[t].cutpool].get_rhs()[i]`.
    pub fn mark_propagate_cut(&mut self, ms: &HighsMipSolver, reason: HighsDomain_Reason) {
        match reason.r#type {
            HighsDomain_Reason::kUnknown
            | HighsDomain_Reason::kCliqueTable
            | HighsDomain_Reason::kBranching
            | HighsDomain_Reason::kModelRowLower
            | HighsDomain_Reason::kModelRowUpper
            | HighsDomain_Reason::kConflictingBounds
            | HighsDomain_Reason::kObjective => {}
            t => {
                let ncut = self.cutpoolpropagation.len() as HighsInt;
                if t < ncut {
                    let k = t as usize;
                    let pool = self.cutpoolpropagation[k].cutpool;
                    let rhs = ms.mipdata_.cutpools[pool].get_rhs()[reason.index as usize];
                    self.cutpoolpropagation[k].mark_propagate_cut(
                        rhs,
                        reason.index,
                        ms.options_mip_.h16_cdouble,
                    );
                } else {
                    self.conflict_pool_propagation[(t - ncut) as usize]
                        .mark_propagate_conflict(reason.index);
                }
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::markPropagate (lines 1936-1958)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::markPropagate (lines 2185-2201)
    /// HiGHS 1.16 (#3575): HighsDomain.cpp::markPropagate (lines 2185-2201)
    /// C++: `void markPropagate(HighsInt row)`
    pub fn mark_propagate(&mut self, ms: &HighsMipSolver, row: HighsInt) {
        let rowu = row as usize;
        if self.propagateflags_[rowu] == 0 {
            let row_lower = ms.row_lower_2(row);
            let row_upper = ms.row_upper_2(row);
            // HiGHS 1.16 #3575 (h16_row_redundant): HighsDomain.cpp lines 2187 and 2191 (1.16.0)
            // The 1.16 text includes its CDouble compares, so it subsumes h16_cdouble here.
            let h16_row_redundant = ms.options_mip_.h16_row_redundant;
            let proplower = if h16_row_redundant {
                !self.is_row_lower_redundant(ms, row)
                    && (self.activitymaxinf_[rowu] == 1
                        || (self.activitymax_[rowu] - row_lower) <= self.capacity_threshold_[rowu])
            } else {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 2189 and 2193 (1.16.0)
                row_lower != -kHighsInf
                    && (self.activitymininf_[rowu] != 0
                        || self.activitymin_[rowu] < row_lower - ms.mipdata_.feastol)
                    && (self.activitymaxinf_[rowu] == 1
                        || if ms.options_mip_.h16_cdouble {
                            (self.activitymax_[rowu] - row_lower) <= self.capacity_threshold_[rowu]
                        } else {
                            (f64::from(self.activitymax_[rowu]) - row_lower)
                                <= self.capacity_threshold_[rowu]
                        })
            };
            let propupper = if h16_row_redundant {
                !self.is_row_upper_redundant(ms, row)
                    && (self.activitymininf_[rowu] == 1
                        || (row_upper - self.activitymin_[rowu]) <= self.capacity_threshold_[rowu])
            } else {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 2189 and 2193 (1.16.0)
                row_upper != kHighsInf
                    && (self.activitymaxinf_[rowu] != 0
                        || self.activitymax_[rowu] > row_upper + ms.mipdata_.feastol)
                    && (self.activitymininf_[rowu] == 1
                        || if ms.options_mip_.h16_cdouble {
                            (row_upper - self.activitymin_[rowu]) <= self.capacity_threshold_[rowu]
                        } else {
                            (row_upper - f64::from(self.activitymin_[rowu]))
                                <= self.capacity_threshold_[rowu]
                        })
            };

            if proplower || propupper {
                self.propagateinds_.push(row);
                self.propagateflags_[rowu] = 1;
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::computeRowActivities (lines 1960-1989)
    /// C++: `void computeRowActivities()`
    pub fn compute_row_activities(&mut self, ms: &HighsMipSolver) {
        let numrow = ms.num_row() as usize;
        self.activitymin_.resize(numrow, HighsCDouble::default());
        self.activitymininf_.resize(numrow, 0);
        self.activitymax_.resize(numrow, HighsCDouble::default());
        self.activitymaxinf_.resize(numrow, 0);
        self.capacity_threshold_.resize(numrow, 0.0);
        self.propagateflags_.resize(numrow, 0);
        self.propagateinds_.reserve(numrow);

        for i in 0..numrow {
            let ii = i as HighsInt;
            let start = ms.mipdata_.a_rstart_[i];
            let end = ms.mipdata_.a_rstart_[i + 1];

            // C++ passes `activitymininf_[i]` / `activitymin_[i]` by reference: taken out and
            // written back, because the body only reads the row data and `self` fields it
            // overwrites first.
            let mut ninf = self.activitymininf_[i];
            let mut act = self.activitymin_[i];
            self.compute_min_activity(
                start,
                end,
                &ms.mipdata_.a_rindex_,
                &ms.mipdata_.a_rvalue_,
                &mut ninf,
                &mut act,
            );
            self.activitymininf_[i] = ninf;
            self.activitymin_[i] = act;

            let mut ninf = self.activitymaxinf_[i];
            let mut act = self.activitymax_[i];
            self.compute_max_activity(
                start,
                end,
                &ms.mipdata_.a_rindex_,
                &ms.mipdata_.a_rvalue_,
                &mut ninf,
                &mut act,
            );
            self.activitymaxinf_[i] = ninf;
            self.activitymax_[i] = act;

            self.recompute_capacity_threshold(ms, ii);

            if (self.activitymininf_[i] <= 1 && ms.row_upper_2(ii) != kHighsInf)
                || (self.activitymaxinf_[i] <= 1 && ms.row_lower_2(ii) != -kHighsInf)
            {
                self.mark_propagate(ms, ii);
                // propagateflags_[i] = 1;
                // propagateinds_.push_back(i);
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::doChangeBound (lines 1991-2021)
    /// C++: `double doChangeBound(const HighsDomainChange & boundchg)`
    // PORT NOTE: signature changed (`ms`, read only, for the activity updates). Private in
    // C++; also called from `backtrack*` / `set_domain_change_stack*` (stack.rs), which replay
    // pool events at their entry.
    pub fn do_change_bound(&mut self, ms: &HighsMipSolver, boundchg: &HighsDomainChange) -> f64 {
        let oldbound: f64;
        let colu = boundchg.column as usize;

        if boundchg.boundtype == HighsBoundType::kLower {
            oldbound = self.col_lower_[colu];
            self.col_lower_[colu] = boundchg.boundval;
            if oldbound != boundchg.boundval {
                if !self.infeasible_ {
                    self.update_activity_lb_change(
                        ms,
                        boundchg.column,
                        oldbound,
                        boundchg.boundval,
                    );
                }

                if !self.is_changed_col(boundchg.column) {
                    self.changedcolsflags_[colu] = 1;
                    self.changedcols_.push(boundchg.column);
                }
            }
        } else {
            oldbound = self.col_upper_[colu];
            self.col_upper_[colu] = boundchg.boundval;
            if oldbound != boundchg.boundval {
                if !self.infeasible_ {
                    self.update_activity_ub_change(
                        ms,
                        boundchg.column,
                        oldbound,
                        boundchg.boundval,
                    );
                }

                if !self.is_changed_col(boundchg.column) {
                    self.changedcolsflags_[colu] = 1;
                    self.changedcols_.push(boundchg.column);
                }
            }
        }

        oldbound
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::changeBound (lines 2023-2105)
    /// HiGHS 1.16 (#3320): HighsDomain.cpp::changeBound (lines 2266-2352)
    /// C++: `void changeBound(HighsDomainChange boundchg, Reason reason = Reason :: branching ( ))`
    /// C++ default arguments: `reason = Reason :: branching ( )`
    ///
    /// PORT NOTE: `self.replay_pool_events(ms)` (M3) runs before the first write to `self`,
    /// right after the `boundval <= col_lower_` / `>= col_upper_` early return. That return
    /// leaves the domain untouched, so a pending event is replayed by the next domain method
    /// against the same bounds (replay is exact between two domain methods); skipping it
    /// there saves the pending-event check on every rejected bound change (perf-gap3-116.md).
    /// The re-entrant hooks at the end (C++ 2096-2102) are associated fns that read
    /// `ms.mipdata_` afresh:
    /// `HighsCliqueTable::add_implications(ms, self, col, val)` and
    /// `HighsImplications::apply_implications(ms, self, col, val)` with
    /// `val = (self.col_lower_[col] > 0.5) as HighsInt`. They require that the clique table,
    /// the implications and every pool this domain propagates are attached (a caller running
    /// inside one of their methods lends it back with `lend!`).
    pub fn change_bound(
        &mut self,
        ms: &mut HighsMipSolver,
        boundchg: HighsDomainChange,
        reason: HighsDomain_Reason,
    ) {
        let mut boundchg = boundchg;
        let colu = boundchg.column as usize;
        // HiGHS 1.16 #3320 (h16_int_bound_round): HighsDomain.cpp lines 2279-2281, 2303-2305 (1.16.0)
        let h16_int_bound_round = ms.options_mip_.h16_int_bound_round;

        let prev_pos: HighsInt;
        if boundchg.boundtype == HighsBoundType::kLower {
            if h16_int_bound_round && ms.is_col_integral(boundchg.column) {
                boundchg.boundval = (boundchg.boundval - ms.mipdata_.feastol).ceil();
            }
            if boundchg.boundval <= self.col_lower_[colu] {
                if reason.r#type != HighsDomain_Reason::kBranching {
                    return;
                }
                boundchg.boundval = self.col_lower_[colu];
            }
            // PORT NOTE (perf): M3 replay deferred past the no-op return above (see doc).
            self.replay_pool_events(ms);
            if boundchg.boundval > self.col_upper_[colu] {
                if boundchg.boundval - self.col_upper_[colu] > ms.mipdata_.feastol {
                    ms.mipdata_.debug_solution.node_pruned(&*self);
                    if !self.infeasible_ {
                        self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                        self.infeasible_ = true;
                        self.infeasible_reason = HighsDomain_Reason::conflicting_bounds(
                            self.domchgstack_.len() as HighsInt,
                        );
                    }
                } else {
                    boundchg.boundval = self.col_upper_[colu];
                    if boundchg.boundval == self.col_lower_[colu] {
                        return;
                    }
                }
            }

            prev_pos = self.col_lower_pos_[colu];
            self.col_lower_pos_[colu] = self.domchgstack_.len() as HighsInt;
        } else {
            if h16_int_bound_round && ms.is_col_integral(boundchg.column) {
                boundchg.boundval = (boundchg.boundval + ms.mipdata_.feastol).floor();
            }
            if boundchg.boundval >= self.col_upper_[colu] {
                if reason.r#type != HighsDomain_Reason::kBranching {
                    return;
                }
                boundchg.boundval = self.col_upper_[colu];
            }
            // PORT NOTE (perf): M3 replay deferred past the no-op return above (see doc).
            self.replay_pool_events(ms);
            if boundchg.boundval < self.col_lower_[colu] {
                if self.col_lower_[colu] - boundchg.boundval > ms.mipdata_.feastol {
                    ms.mipdata_.debug_solution.node_pruned(&*self);
                    if !self.infeasible_ {
                        self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                        self.infeasible_ = true;
                        self.infeasible_reason = HighsDomain_Reason::conflicting_bounds(
                            self.domchgstack_.len() as HighsInt,
                        );
                    }
                } else {
                    boundchg.boundval = self.col_lower_[colu];
                    if boundchg.boundval == self.col_upper_[colu] {
                        return;
                    }
                }
            }

            prev_pos = self.col_upper_pos_[colu];
            self.col_upper_pos_[colu] = self.domchgstack_.len() as HighsInt;
        }

        ms.mipdata_.debug_solution.bound_change_added(
            &*self,
            &boundchg,
            reason.r#type == HighsDomain_Reason::kBranching,
            Some(&reason),
        );

        if reason.r#type == HighsDomain_Reason::kBranching {
            self.branch_pos_.push(self.domchgstack_.len() as HighsInt);
        }

        let binary = self.is_binary(ms, boundchg.column);

        let oldbound = self.do_change_bound(ms, &boundchg);

        self.prevboundval_.push((oldbound, prev_pos));
        self.domchgstack_.push(boundchg.clone());
        self.domchgreason_.push(reason);

        if binary && !self.infeasible_ && self.is_fixed(boundchg.column) {
            let val = (self.col_lower_[colu] > 0.5) as HighsInt;
            HighsCliqueTable::add_implications(ms, self, boundchg.column, val);
            if !self.infeasible_ {
                let val = (self.col_lower_[colu] > 0.5) as HighsInt;
                HighsImplications::apply_implications(ms, self, boundchg.column, val);
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::checkChangeBound (lines 2107-2118)
    /// C++: `bool checkChangeBound(HighsBoundType boundtype, HighsInt col, HighsCDouble boundval, Reason reason)`
    // PORT NOTE: `adjustedLb/Ub` and `feastol()` read `ms`; `changeBound` is `change_bound_2`.
    pub fn check_change_bound(
        &mut self,
        ms: &mut HighsMipSolver,
        boundtype: HighsBoundType,
        col: HighsInt,
        boundval: HighsCDouble,
        reason: HighsDomain_Reason,
    ) -> bool {
        if f64::from(boundval * kHighsTiny).abs() > self.feastol(ms) {
            return false;
        }
        let mut accept = false;
        let bound = if boundtype == HighsBoundType::kLower {
            self.adjusted_lb(ms, col, boundval, &mut accept)
        } else {
            self.adjusted_ub(ms, col, boundval, &mut accept)
        };
        if !accept {
            return false;
        }
        self.change_bound_2(ms, boundtype, col, bound, reason);
        true
    }
}

#[cfg(test)]
mod h16_cdouble_prop_tests {
    use super::*;

    // One continuous column in [0, inf), one row x0 <= 4; feastol 1e-6, epsilon 1e-9.
    fn cont_mip(h16_cdouble: bool) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![kHighsInf];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![4.0];
        lp.integrality_ = vec![HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        ms.mipdata_.a_rstart_ = vec![0, 1];
        ms.mipdata_.a_rindex_ = vec![0];
        ms.mipdata_.a_rvalue_ = vec![1.0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.h16_cdouble = h16_cdouble;
        ms
    }

    // Bound value bv = hi + lo with hi = 1 + 4503600 * 2^-52 (a double) and lo = -0.473 * 2^-52.
    // f64(bv) = hi, so |f64(bv) - 1| = 4503600 * 2^-52 > epsilon (1e-9 = 4503599.627 * 2^-52).
    // The exact |bv - 1| = 4503599.527 * 2^-52 < epsilon, so the 1.16 CDouble compare snaps.
    fn snap_case() -> (f64, HighsCDouble) {
        let eps = f64::EPSILON; // 2^-52
        let hi = 1.0 + 4503600.0 * eps;
        (hi, HighsCDouble::from(hi) + (-0.473 * eps))
    }

    #[test]
    fn h16_cdouble_prop_on_adjusted_bounds_snap() {
        let (_hi, bv) = snap_case();
        let ms = cont_mip(true);
        let mut dom = HighsDomain::new(&ms);
        // adjustedUb: lower bound 1.0 and no upper bound -> accepted, snapped to the lower bound.
        dom.col_lower_[0] = 1.0;
        dom.col_upper_[0] = kHighsInf;
        let mut accept = false;
        let b = dom.adjusted_ub(&ms, 0, bv, &mut accept);
        assert!(accept);
        assert_eq!(b, 1.0);
        // adjustedLb: upper bound 1.0 and no lower bound -> accepted, snapped to the upper bound.
        dom.col_lower_[0] = -kHighsInf;
        dom.col_upper_[0] = 1.0;
        let mut accept = false;
        let b = dom.adjusted_lb(&ms, 0, bv, &mut accept);
        assert!(accept);
        assert_eq!(b, 1.0);
    }

    #[test]
    fn h16_cdouble_prop_off_adjusted_bounds_keep_hi() {
        let (hi, bv) = snap_case();
        let ms = cont_mip(false);
        let mut dom = HighsDomain::new(&ms);
        // 1.15.1: f64(bv) = hi is not within epsilon of the bound, so it is kept.
        dom.col_lower_[0] = 1.0;
        dom.col_upper_[0] = kHighsInf;
        let mut accept = false;
        let b = dom.adjusted_ub(&ms, 0, bv, &mut accept);
        assert!(accept);
        assert_eq!(b, hi);
        dom.col_lower_[0] = -kHighsInf;
        dom.col_upper_[0] = 1.0;
        let mut accept = false;
        let b = dom.adjusted_lb(&ms, 0, bv, &mut accept);
        assert!(accept);
        assert_eq!(b, hi);
    }

    // U = 1 + 4503600 * 2^-52; obj = 1 + 0.45 * 2^-52 (hi 1.0, lo 0.45 ulp); capacity 1e-9.
    // Off: U - f64(obj) = 4503600 * 2^-52 > capacity, so the objective is not propagated.
    // On: U - obj = 4503599.55 * 2^-52 <= capacity, so it is.
    fn should_be_propagated_case(h16_cdouble: bool) -> bool {
        let eps = f64::EPSILON;
        let upper = 1.0 + 4503600.0 * eps;
        let mut ms = cont_mip(h16_cdouble);
        ms.mipdata_.upper_limit = upper;
        let mut dom = HighsDomain::new(&ms);
        dom.obj_prop_.is_propagated = false;
        dom.obj_prop_.num_inf_obj_lower = 0;
        dom.obj_prop_.objective_lower = HighsCDouble::from(1.0) + (0.45 * eps);
        dom.obj_prop_.capacity_threshold = 1e-9;
        dom.obj_prop_should_be_propagated(&ms)
    }

    #[test]
    fn h16_cdouble_prop_on_should_be_propagated() {
        assert!(should_be_propagated_case(true));
    }

    #[test]
    fn h16_cdouble_prop_off_should_be_propagated() {
        assert!(!should_be_propagated_case(false));
    }

    // Row 1 <= row activity; activitymax = hi + lo with hi = 1 + 4503600 * 2^-52 and
    // lo = -0.45 * 2^-52. Off: f64(activitymax) - 1 = 4503600 * 2^-52 > capacity 1e-9, so the
    // row is not marked. On: the CDouble difference 4503599.55 * 2^-52 <= capacity, so it is.
    fn mark_propagate_case(h16_cdouble: bool) -> Vec<HighsInt> {
        let eps = f64::EPSILON;
        let hi = 1.0 + 4503600.0 * eps;
        let mut ms = cont_mip(h16_cdouble);
        ms.orig_model_.row_lower_ = vec![1.0];
        ms.orig_model_.row_upper_ = vec![kHighsInf];
        let mut dom = HighsDomain::new(&ms);
        dom.activitymin_ = vec![HighsCDouble::from(-1.0)];
        dom.activitymax_ = vec![HighsCDouble::from(hi) + (-0.45 * eps)];
        dom.activitymininf_ = vec![0];
        dom.activitymaxinf_ = vec![0];
        dom.capacity_threshold_ = vec![1e-9];
        dom.propagateflags_ = vec![0];
        dom.propagateinds_.clear();
        dom.mark_propagate(&ms, 0);
        dom.propagateinds_.clone()
    }

    #[test]
    fn h16_cdouble_prop_on_mark_propagate() {
        assert_eq!(mark_propagate_case(true), vec![0]);
    }

    #[test]
    fn h16_cdouble_prop_off_mark_propagate() {
        assert!(mark_propagate_case(false).is_empty());
    }
}

#[cfg(test)]
mod h16_int_bound_round_tests {
    use super::*;

    /// One column in [0, 10] with the row x <= 10; `integral` picks its type.
    /// Applies one bound change with `h16_int_bound_round = round` and returns the
    /// stored (lower, upper) bounds of the root domain.
    fn h16_int_bound_round_change(
        integral: bool,
        round: bool,
        boundtype: HighsBoundType,
        boundval: f64,
    ) -> (f64, f64) {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![10.0];
        lp.integrality_ = vec![if integral {
            HighsVarType::kInteger
        } else {
            HighsVarType::kContinuous
        }];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];

        let mut options = HighsOptions::default();
        options.h16_int_bound_round = round;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.new_data();
        ms.mipdata_.feastol = 1e-6;
        ms.setup_domain_propagation();

        let chg = HighsDomainChange {
            boundval,
            column: 0,
            boundtype,
        };
        crate::detach!(ms.mipdata_.domains[0], |d| d.change_bound(
            &mut ms,
            chg,
            HighsDomain_Reason::unspecified()
        ));
        (
            ms.mipdata_.domains[0].col_lower_[0],
            ms.mipdata_.domains[0].col_upper_[0],
        )
    }

    #[test]
    fn h16_int_bound_round_on_lower_ceils_integer_bound() {
        // ceil(2.0000001 - 1e-6) = 2 and ceil(1.9999999 - 1e-6) = 2 (1.16 changeBound).
        let (lo, _) = h16_int_bound_round_change(true, true, HighsBoundType::kLower, 2.0000001);
        assert_eq!(lo, 2.0);
        let (lo, _) = h16_int_bound_round_change(true, true, HighsBoundType::kLower, 1.9999999);
        assert_eq!(lo, 2.0);
    }

    #[test]
    fn h16_int_bound_round_on_upper_floors_integer_bound() {
        // floor(5.0000001 + 1e-6) = 5 and floor(4.9999999 + 1e-6) = 5 (1.16 changeBound).
        let (_, up) = h16_int_bound_round_change(true, true, HighsBoundType::kUpper, 5.0000001);
        assert_eq!(up, 5.0);
        let (_, up) = h16_int_bound_round_change(true, true, HighsBoundType::kUpper, 4.9999999);
        assert_eq!(up, 5.0);
    }

    #[test]
    fn h16_int_bound_round_on_continuous_unchanged() {
        // The rounding applies only to integral columns: a continuous bound is stored raw.
        let (lo, _) = h16_int_bound_round_change(false, true, HighsBoundType::kLower, 2.0000001);
        assert_eq!(lo, 2.0000001);
    }

    #[test]
    fn h16_int_bound_round_off_stores_raw_bound() {
        // 1.15.1 behaviour: the integer bounds are stored exactly as given.
        let (lo, _) = h16_int_bound_round_change(true, false, HighsBoundType::kLower, 2.0000001);
        assert_eq!(lo, 2.0000001);
        let (lo, _) = h16_int_bound_round_change(true, false, HighsBoundType::kLower, 1.9999999);
        assert_eq!(lo, 1.9999999);
        let (_, up) = h16_int_bound_round_change(true, false, HighsBoundType::kUpper, 5.0000001);
        assert_eq!(up, 5.0000001);
    }
}

#[cfg(test)]
mod h16_bound_range_inf_tests {
    use super::*;

    // One continuous column x0 in [0, inf), one row x0 <= 4 with coefficient 1; feastol 1e-6.
    fn inf_mip(h16_bound_range_inf: bool) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![kHighsInf];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![4.0];
        lp.integrality_ = vec![HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        ms.mipdata_.a_rstart_ = vec![0, 1];
        ms.mipdata_.a_rindex_ = vec![0];
        ms.mipdata_.a_rvalue_ = vec![1.0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.h16_bound_range_inf = h16_bound_range_inf;
        ms
    }

    #[test]
    fn h16_bound_range_inf_on_bound_range() {
        let tol = 1e-6;
        assert_eq!(
            bound_range(kHighsInf, 0.0, tol, HighsVarType::kContinuous, true, 0.3),
            kHighsInf
        );
        assert_eq!(
            bound_range(1.0, -kHighsInf, tol, HighsVarType::kContinuous, true, 0.3),
            kHighsInf
        );
        assert_eq!(
            bound_range(1.0, -kHighsInf, tol, HighsVarType::kInteger, true, 0.3),
            kHighsInf
        );
    }

    #[test]
    fn h16_bound_range_inf_off_bound_range() {
        let tol = 1e-6;
        // 1.15.1: inf - 0 = inf, and max(0.3 * inf, 1000 * tol) = inf, so inf - inf is NaN.
        assert!(bound_range(kHighsInf, 0.0, tol, HighsVarType::kContinuous, false, 0.3).is_nan());
        // Finite bounds: the flag changes nothing.
        for (upper, lower) in [(4.0, 1.0), (-2.0, -3.5), (0.0, 0.0)] {
            for var_type in [HighsVarType::kContinuous, HighsVarType::kInteger] {
                assert_eq!(
                    bound_range(upper, lower, tol, var_type, false, 0.3).to_bits(),
                    bound_range(upper, lower, tol, var_type, true, 0.3).to_bits()
                );
            }
        }
    }

    #[test]
    fn h16_bound_range_inf_on_row_capacity() {
        // Column 0 has no upper bound: the row's capacity threshold is infinite.
        let ms = inf_mip(true);
        let mut dom = HighsDomain::new(&ms);
        dom.col_upper_[0] = kHighsInf;
        dom.capacity_threshold_ = vec![0.0];
        dom.recompute_capacity_threshold(&ms, 0);
        assert_eq!(dom.capacity_threshold_[0], kHighsInf);
    }

    #[test]
    fn h16_bound_range_inf_off_row_capacity() {
        // 1.15.1: the NaN range is absorbed by cpp_max, leaving feastol.
        let ms = inf_mip(false);
        let mut dom = HighsDomain::new(&ms);
        dom.col_upper_[0] = kHighsInf;
        dom.capacity_threshold_ = vec![0.0];
        dom.recompute_capacity_threshold(&ms, 0);
        assert_eq!(dom.capacity_threshold_[0], 1e-6);
    }

    #[test]
    fn h16_bound_range_inf_on_update_threshold() {
        // Lower bound -inf: the new threshold for the upper-bound change is infinite.
        let ms = inf_mip(true);
        let mut dom = HighsDomain::new(&ms);
        dom.col_lower_[0] = -kHighsInf;
        dom.col_upper_[0] = 1.0;
        let mut threshold = 0.0;
        dom.update_threshold_ub_change(&ms, 0, 0.5, 1.0, &mut threshold);
        assert_eq!(threshold, kHighsInf);
    }

    #[test]
    fn h16_bound_range_inf_off_update_threshold() {
        let ms = inf_mip(false);
        let mut dom = HighsDomain::new(&ms);
        dom.col_lower_[0] = -kHighsInf;
        dom.col_upper_[0] = 1.0;
        let mut threshold = 0.0;
        dom.update_threshold_ub_change(&ms, 0, 0.5, 1.0, &mut threshold);
        assert_eq!(threshold, 1e-6);
    }
}

#[cfg(test)]
mod h16_row_redundant_tests {
    use super::*;

    // One continuous column in [0, inf) with coefficient 1 in one row; feastol 1e-6.
    fn row_mip(row_lower: f64, row_upper: f64, h16_row_redundant: bool) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![kHighsInf];
        lp.row_lower_ = vec![row_lower];
        lp.row_upper_ = vec![row_upper];
        lp.integrality_ = vec![HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        ms.mipdata_.a_rstart_ = vec![0, 1];
        ms.mipdata_.a_rindex_ = vec![0];
        ms.mipdata_.a_rvalue_ = vec![1.0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.h16_row_redundant = h16_row_redundant;
        ms
    }

    // Row upper 2.0 and d1 = round(2.0 + 1e-6). The double sum 2.0 + feastol is d1, so the 1.15.1
    // test max activity d1 <= d1 holds. d1 - 2.0 = 1e-6 + 1.4e-16 exactly (about 0.3 ulp(2.0)
    // above feastol), so the 1.16 CDouble test (activitymax - rowUpper <= feastol) fails.
    fn upper_d1() -> f64 {
        2.0 + 1e-6
    }

    #[test]
    fn h16_row_redundant_on_is_redundant_row() {
        // 1.16: the row is not redundant, the CDouble difference exceeds feastol.
        let ms = row_mip(-kHighsInf, 2.0, true);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymin_ = vec![HighsCDouble::from(0.0)];
        dom.activitymininf_ = vec![0];
        dom.activitymax_ = vec![HighsCDouble::from(upper_d1())];
        dom.activitymaxinf_ = vec![0];
        assert!(!dom.is_redundant_row(&ms, 0));
    }

    #[test]
    fn h16_row_redundant_off_is_redundant_row() {
        // 1.15.1: the double max activity d1 <= rowUpper + feastol holds, so the row is redundant.
        let ms = row_mip(-kHighsInf, 2.0, false);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymin_ = vec![HighsCDouble::from(0.0)];
        dom.activitymininf_ = vec![0];
        dom.activitymax_ = vec![HighsCDouble::from(upper_d1())];
        dom.activitymaxinf_ = vec![0];
        assert!(dom.is_redundant_row(&ms, 0));
    }

    #[test]
    fn h16_row_redundant_on_upper_infeasible_helper() {
        // isRowUpperInfeasible: activitymin - rowUpper > feastol, in CDouble.
        let ms = row_mip(-kHighsInf, 2.0, true);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymininf_ = vec![0];
        dom.activitymin_ = vec![HighsCDouble::from(upper_d1())];
        assert!(dom.is_row_upper_infeasible(&ms, 0));
        dom.activitymin_ = vec![HighsCDouble::from(2.0)];
        assert!(!dom.is_row_upper_infeasible(&ms, 0));
    }

    #[test]
    fn h16_row_redundant_on_lower_infeasible_helper() {
        // isRowLowerInfeasible: rowLower - activitymax > feastol, in CDouble.
        let ms = row_mip(2.0 + 1e-6, kHighsInf, true);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymaxinf_ = vec![0];
        dom.activitymax_ = vec![HighsCDouble::from(2.0)];
        assert!(dom.is_row_lower_infeasible(&ms, 0));
        dom.activitymax_ = vec![HighsCDouble::from(upper_d1())];
        assert!(!dom.is_row_lower_infeasible(&ms, 0));
    }

    fn mark_propagate_case(h16_row_redundant: bool) -> Vec<HighsInt> {
        // Lower side free, upper side as above with min activity -1 and a loose threshold.
        let ms = row_mip(-kHighsInf, 2.0, h16_row_redundant);
        let mut dom = HighsDomain::new(&ms);
        dom.activitymin_ = vec![HighsCDouble::from(-1.0)];
        dom.activitymininf_ = vec![0];
        dom.activitymax_ = vec![HighsCDouble::from(upper_d1())];
        dom.activitymaxinf_ = vec![0];
        dom.capacity_threshold_ = vec![10.0];
        dom.propagateflags_ = vec![0];
        dom.propagateinds_.clear();
        dom.mark_propagate(&ms, 0);
        dom.propagateinds_.clone()
    }

    #[test]
    fn h16_row_redundant_on_mark_propagate() {
        assert_eq!(mark_propagate_case(true), vec![0]);
    }

    #[test]
    fn h16_row_redundant_off_mark_propagate() {
        assert!(mark_propagate_case(false).is_empty());
    }
}
