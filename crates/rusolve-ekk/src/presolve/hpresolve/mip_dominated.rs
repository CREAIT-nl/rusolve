//! HiGHS source: presolve/HPresolve.cpp (dominatedColumns 1045-1454, enumerateSolutions 5223-5696)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! MIP-only methods. Read `presolve/hpresolve/mod.rs` first
//! (ownership, "MIP mode", containers, slices, bound arrays, lambdas).

use super::*;
use crate::detach;
use crate::prelude::*;

/// Captured locals of `HPresolve::dominatedColumns` (C++ lambdas capture them by reference).
struct DomCtx<'a> {
    /// C++ `postsolve_stack`
    ps: &'a mut HighsPostsolveStack,
    /// C++ `signatures`: (rows with finite lower side, rows with finite upper side) bit sets
    signatures: Vec<(u32, u32)>,
    /// C++ `numDomChecks` (counted, never read)
    #[allow(dead_code)]
    num_dom_checks: usize,
    /// C++ `numDomChecksPredBndAnalysis`
    num_dom_checks_pred_bnd_analysis: usize,
    /// C++ `numFixedCols`
    num_fixed_cols: HighsInt,
    /// C++ `numFixedColsPredBndAnalysis`
    num_fixed_cols_pred_bnd_analysis: HighsInt,
    /// C++ `numModifiedBndsPredBndAnalysis`
    num_modified_bnds_pred_bnd_analysis: HighsInt,
    /// HiGHS 1.16 #3472 (h16_dominated_cols) C++ `check_rowpositions`: row positions copied by
    /// `getRowPositions` in checkRow (1.16.0 line 148); inert when the flag is off
    check_rowpositions: Vec<HighsInt>,
}

/// C++ lambda `isBinary` of dominatedColumns (reads the presolved model).
fn dom_is_binary(this: &HPresolve, col: HighsInt) -> bool {
    let i = col as usize;
    this.model.integrality_[i] == HighsVarType::kInteger
        && this.model.col_lower_[i] == 0.0
        && this.model.col_upper_[i] == 1.0
}

/// C++ lambda `addSignature` of dominatedColumns.
fn dom_add_signature(
    dc: &mut DomCtx<'_>,
    row: HighsInt,
    col: HighsInt,
    row_lower_finite: u32,
    row_upper_finite: u32,
) {
    let row_hashed_pos = hash_i32(row) >> 59;
    debug_assert!(row_hashed_pos < 32);
    dc.signatures[col as usize].0 |= row_lower_finite << row_hashed_pos;
    dc.signatures[col as usize].1 |= row_upper_finite << row_hashed_pos;
}

/// C++ lambda `checkDominationNonZero` of dominatedColumns.
fn dom_check_domination_nonzero(this: &HPresolve, row: HighsInt, mut aj: f64, mut ak: f64) -> bool {
    let ri = row as usize;
    if this.is_ranged(row) {
        // the row is an equality or ranged row, therefore the coefficients must
        // be parallel, otherwise one of the inequalities given by the row rules
        // out domination
        if (aj - ak).abs() > this.options.small_matrix_value {
            return false;
        }
        return true;
    }

    // normalize row to a <= constraint
    if this.model.row_upper_[ri] == kHighsInf {
        aj = -aj;
        ak = -ak;
    }

    // the coefficient of the dominating column needs to be smaller than or
    // equal to the coefficient of the dominated column
    if aj > ak + this.options.small_matrix_value {
        return false;
    }
    true
}

/// C++ lambda `checkDomination` of dominatedColumns.
fn dom_check_domination(
    this: &mut HPresolve,
    dc: &mut DomCtx<'_>,
    scalj: HighsInt,
    j: HighsInt,
    scalk: HighsInt,
    k: HighsInt,
) -> bool {
    let (ji, ki) = (j as usize, k as usize);
    // increment counter
    dc.num_dom_checks += 1;

    // rule out domination from integers to continuous variables
    if this.model.integrality_[ji] == HighsVarType::kInteger
        && this.model.integrality_[ki] == HighsVarType::kContinuous
    {
        return false;
    }

    // check the signatures
    let mut sj_minus = dc.signatures[ji].0;
    let mut sj_plus = dc.signatures[ji].1;
    if scalj == -1 {
        std::mem::swap(&mut sj_plus, &mut sj_minus);
    }

    let mut sk_minus = dc.signatures[ki].0;
    let mut sk_plus = dc.signatures[ki].1;
    if scalk == -1 {
        std::mem::swap(&mut sk_plus, &mut sk_minus);
    }

    // the set of rows with a negative coefficient must be a superset of the
    // dominated columns set of rows with a negative coefficient
    if (!sj_minus & sk_minus) != 0 {
        return false;
    }

    // the set of rows with a positive coefficient must be a subset of the
    // dominated columns set of rows with a positive coefficient
    if (sj_plus & !sk_plus) != 0 {
        return false;
    }

    // the dominating columns cost must be smaller or equal to the dominated
    // columns cost
    if scalj as f64 * this.model.col_cost_[ji]
        > scalk as f64 * this.model.col_cost_[ki] + this.options.small_matrix_value
    {
        return false;
    }

    // finally check the column vectors
    let mut it = this.get_column_vector(j);
    while let Some(nz) = it.next_nz(&*this) {
        let row = nz.index();
        let ak_pos = this.find_nonzero(row, k);
        let ak_val = if ak_pos == -1 {
            0.0
        } else {
            this.a_value[ak_pos as usize]
        };
        if !dom_check_domination_nonzero(
            this,
            row,
            scalj as f64 * nz.value(),
            scalk as f64 * ak_val,
        ) {
            return false;
        }
    }

    // check row only occurring in the column vector of k
    let mut it = this.get_column_vector(k);
    while let Some(nz) = it.next_nz(&*this) {
        let row = nz.index();
        let aj_pos = this.find_nonzero(row, j);
        // only rows in which aj does not occur are left to check
        if aj_pos != -1 {
            continue;
        }
        if !dom_check_domination_nonzero(this, row, 0.0, scalk as f64 * nz.value()) {
            return false;
        }
    }

    true
}

/// C++ lambda `fixCol` of dominatedColumns.
fn dom_fix_col(
    this: &mut HPresolve,
    dc: &mut DomCtx<'_>,
    col: HighsInt,
    direction: HighsInt,
) -> HPresolveResult {
    dc.num_fixed_cols += 1;
    if direction > 0 {
        hpresolve_checked_call!(this.fix_col_to_upper(&mut *dc.ps, col));
    } else {
        hpresolve_checked_call!(this.fix_col_to_lower(&mut *dc.ps, col));
    }
    // remove row singletons and doubleton equations
    hpresolve_checked_call!(this.remove_row_singletons(&mut *dc.ps));
    hpresolve_checked_call!(this.remove_doubleton_equations(&mut *dc.ps));
    HPresolveResult::kOk
}

/// HiGHS 1.16 (#3411, #3285, #3290): HPresolveDominatedCols.cpp::dominatedColumns tightenBounds lambda (1.16.0 lines 211-280); gates h16_int_bound_rounding, h16_dominated_cols
/// C++ lambda `tightenBounds` of dominatedColumns.
#[allow(clippy::too_many_arguments)]
fn dom_tighten_bounds(
    this: &mut HPresolve,
    dc: &mut DomCtx<'_>,
    col: HighsInt,
    col_bound: f64,
    col_is_at_upper: bool,
    other_col: HighsInt,
    other_col_bound: f64,
    other_col_coeff_pattern: HighsInt,
) -> HPresolveResult {
    // HiGHS 1.16 #3411 / #3285 / #3290 (1.16.0 tightenBounds lambda): flags read once at the top
    let h16_int_bound_rounding = this.options.h16_int_bound_rounding;
    let h16_dominated_cols = this.options.h16_dominated_cols;
    let ci = col as usize;
    // bound should be finite
    debug_assert!(other_col_bound.abs() != kHighsInf);
    // return if variable is already fixed
    if this.model.col_lower_[ci] == this.model.col_upper_[ci] {
        return HPresolveResult::kOk;
    }
    // initialise bounds
    let mut lower_bound: f64;
    let mut upper_bound: f64;
    // predictive bound analysis, see Theorem 3 from Gamrath et al.'s paper
    if col_is_at_upper {
        // (i) x_j <= MINL^k_j(otherColBound)
        upper_bound = this.compute_implied_upper_bound(
            col,
            other_col,
            other_col_bound,
            other_col_coeff_pattern,
        );
        // (iii) x_j >= min{colBound, MAXL^k_j(otherColBound)}
        lower_bound = cpp_min(
            col_bound,
            this.compute_implied_lower_bound(
                col,
                other_col,
                other_col_bound,
                other_col_coeff_pattern,
            ),
        );
        if this.model.col_cost_[ci] <= 0.0 {
            // (v) if c_j <= 0, then x_j >= min{colBound,
            //                                  MINU^k_j(otherColBound)}
            let mut worst_case_upper = this.compute_worst_case_upper_bound(
                col,
                other_col,
                other_col_bound,
                other_col_coeff_pattern,
            );
            // HiGHS 1.16 #3411: 1.16.0 drops this rounding (computeBound rounds instead)
            if h16_int_bound_rounding {
            } else if this.model.integrality_[ci] != HighsVarType::kContinuous {
                worst_case_upper = (worst_case_upper + this.primal_feastol).floor();
            }
            lower_bound = cpp_max(lower_bound, cpp_min(col_bound, worst_case_upper));
        }
    } else {
        // (ii) x_k >= MAXL^j_k(otherColBound)
        lower_bound = this.compute_implied_lower_bound(
            col,
            other_col,
            other_col_bound,
            other_col_coeff_pattern,
        );
        // (iv) x_k <= max{colBound, MINL^j_k(otherColBound)}
        upper_bound = cpp_max(
            col_bound,
            this.compute_implied_upper_bound(
                col,
                other_col,
                other_col_bound,
                other_col_coeff_pattern,
            ),
        );
        if this.model.col_cost_[ci] >= 0.0 {
            // (vi) if c_k >= 0, then x_k <= max{colBound,
            //                                   MAXU^j_k(otherColBound)}
            let mut worst_case_lower = this.compute_worst_case_lower_bound(
                col,
                other_col,
                other_col_bound,
                other_col_coeff_pattern,
            );
            // HiGHS 1.16 #3411: 1.16.0 drops this rounding (computeBound rounds instead)
            if h16_int_bound_rounding {
            } else if this.model.integrality_[ci] != HighsVarType::kContinuous {
                worst_case_lower = (worst_case_lower - this.primal_feastol).ceil();
            }
            upper_bound = cpp_min(upper_bound, cpp_max(col_bound, worst_case_lower));
        }
    }
    // update bounds
    // HiGHS 1.16 #3285: 1.16.0 guards the update against an infinite bound (lines 259-260)
    let lower_update = if h16_dominated_cols {
        lower_bound < kHighsInf && lower_bound > this.model.col_lower_[ci] + this.primal_feastol
    } else {
        lower_bound > this.model.col_lower_[ci] + this.primal_feastol
    };
    if lower_update {
        // HiGHS 1.16 #3411: 1.16.0 drops this rounding (computeBound rounds instead)
        if h16_int_bound_rounding {
        } else if this.model.integrality_[ci] != HighsVarType::kContinuous {
            lower_bound = (lower_bound - this.primal_feastol).ceil();
        }
        if lower_bound == this.model.col_upper_[ci] {
            dc.num_fixed_cols_pred_bnd_analysis += 1;
            // HiGHS 1.16 #3290: 1.16.0 returns the fixCol result directly (line 263)
            if h16_dominated_cols {
                return dom_fix_col(this, dc, col, 1);
            }
            hpresolve_checked_call!(dom_fix_col(this, dc, col, 1));
        } else if this.model.integrality_[ci] != HighsVarType::kContinuous {
            dc.num_modified_bnds_pred_bnd_analysis += 1;
            hpresolve_checked_call!(this.change_col_lower(col, lower_bound));
        }
    }
    // HiGHS 1.16 #3285: 1.16.0 guards the update against an infinite bound (lines 269-270)
    let upper_update = if h16_dominated_cols {
        upper_bound > -kHighsInf && upper_bound < this.model.col_upper_[ci] - this.primal_feastol
    } else {
        upper_bound < this.model.col_upper_[ci] - this.primal_feastol
    };
    if upper_update {
        // HiGHS 1.16 #3411: 1.16.0 drops this rounding (computeBound rounds instead)
        if h16_int_bound_rounding {
        } else if this.model.integrality_[ci] != HighsVarType::kContinuous {
            upper_bound = (upper_bound + this.primal_feastol).floor();
        }
        if upper_bound == this.model.col_lower_[ci] {
            dc.num_fixed_cols_pred_bnd_analysis += 1;
            // HiGHS 1.16 #3290: 1.16.0 returns the fixCol result directly (line 273)
            if h16_dominated_cols {
                return dom_fix_col(this, dc, col, -1);
            }
            hpresolve_checked_call!(dom_fix_col(this, dc, col, -1));
        } else if this.model.integrality_[ci] != HighsVarType::kContinuous {
            dc.num_modified_bnds_pred_bnd_analysis += 1;
            hpresolve_checked_call!(this.change_col_upper(col, upper_bound));
        }
    }
    HPresolveResult::kOk
}

/// C++ lambda `varsHaveSameType` of dominatedColumns.
fn dom_vars_have_same_type(this: &HPresolve, col1: HighsInt, col2: HighsInt) -> bool {
    let i1 = this.model.integrality_[col1 as usize];
    let i2 = this.model.integrality_[col2 as usize];
    (i1 != HighsVarType::kContinuous && i2 != HighsVarType::kContinuous)
        || (i1 != HighsVarType::kInteger && i2 != HighsVarType::kInteger)
}

/// C++ lambda `checkCols` of dominatedColumns.
#[allow(clippy::too_many_arguments)]
fn dom_check_cols(
    this: &mut HPresolve,
    dc: &mut DomCtx<'_>,
    col: HighsInt,
    k: HighsInt,
    direction: HighsInt,
    multiplier: HighsInt,
    bound_implied: bool,
    has_cliques: bool,
    same_var_type: bool,
) -> HPresolveResult {
    let (ci, ki) = (col as usize, k as usize);
    // compute direction for dominated variable
    // (1: lower bound; -1: upper bound)
    let direction_k: HighsInt = multiplier * direction;
    // get bounds
    let dominating_bound = if direction > 0 {
        this.model.col_upper_[ci]
    } else {
        this.model.col_lower_[ci]
    };
    let dominated_bound = if direction_k > 0 {
        this.model.col_lower_[ki]
    } else {
        this.model.col_upper_[ki]
    };
    // check if bounds are finite
    // HiGHS quirk: the test is `direction * bound != kHighsInf` (not a sign-aware
    // finiteness test); it is kept as written.
    let is_dominating_bound_finite = direction as f64 * dominating_bound != kHighsInf;
    // HiGHS quirk: compares `direction_k * bound` with `-kHighsInf` (kept as written).
    let is_dominated_bound_finite = direction_k as f64 * dominated_bound != -kHighsInf;
    // check whether variable 'k' can potentially be fixed. check if there are
    // cliques in advance instead of directly searching for a common clique.
    let try_to_fix = is_dominated_bound_finite && (bound_implied || has_cliques);
    // check whether predictive bound analysis can be performed. both
    // variables need to have the same type.
    let try_to_strengthen_bounds =
        (is_dominating_bound_finite || is_dominated_bound_finite) && same_var_type;
    if try_to_fix || try_to_strengthen_bounds {
        // increment counter for number of domination checks due to predictive
        // bound analysis
        if !try_to_fix {
            dc.num_dom_checks_pred_bnd_analysis += 1;
        }
        // check for domination
        if dom_check_domination(this, dc, direction, col, direction_k, k) {
            // Re-check the implied bound condition since earlier fixings in
            // this dominatedColumns call may have changed the model state
            let current_bound_implied = if direction > 0 {
                this.is_upper_implied(col)
            } else {
                this.is_lower_implied(col)
            };
            if try_to_fix
                && (current_bound_implied
                    || this.mip_mut().mipdata_.cliquetable.have_common_clique(
                        HighsCliqueTable_CliqueVar::new(col, (direction > 0) as HighsInt),
                        HighsCliqueTable_CliqueVar::new(k, (direction_k > 0) as HighsInt),
                    ))
            {
                // direction =  1, multiplier =  1:
                // case (i)   ub(x_j) =  inf,  x_j >  x_k: set x_k = lb(x_k)
                // direction =  1, multiplier = -1:
                // case (ii)  ub(x_j) =  inf,  x_j > -x_k: set x_k = ub(x_k)
                // direction = -1, multiplier =  1:
                // case (iii) lb(x_j) = -inf, -x_j > -x_k: set x_k = ub(x_k)
                // direction = -1, multiplier = -1:
                // case (iv)  lb(x_j) = -inf, -x_j >  x_k: set x_k = lb(x_k)
                hpresolve_checked_call!(dom_fix_col(this, dc, k, -direction_k));
            } else if try_to_strengthen_bounds {
                // tighten bounds via predictive bound analysis, see Theorem 3
                // from Gamrath et al.'s paper
                if is_dominated_bound_finite {
                    hpresolve_checked_call!(dom_tighten_bounds(
                        this,
                        dc,
                        col,
                        dominating_bound,
                        direction > 0,
                        k,
                        dominated_bound,
                        direction * direction_k,
                    ));
                }
                if this.col_deleted[ci] == 0 && is_dominating_bound_finite {
                    hpresolve_checked_call!(dom_tighten_bounds(
                        this,
                        dc,
                        k,
                        dominated_bound,
                        direction_k < 0,
                        col,
                        dominating_bound,
                        direction * direction_k,
                    ));
                }
            }
        }
    }
    HPresolveResult::kOk
}

/// HiGHS 1.16 (#3472): HPresolveDominatedCols.cpp::dominatedColumns checkRow lambda (lines 351-396, 1.16.0), gate h16_dominated_cols
/// On path of `dom_check_row` when `h16_dominated_cols` is set. Iterates a copy of the row
/// positions (`getRowPositions`, 1.16.0 line 354) and skips stale positions (line 356), instead
/// of the stored-row slice of 1.15.1. `bound_implied` is the caller's value, as in 1.15.1; the
/// 1.16.0 per-row recompute of `boundImplied` (line 371) is #3218, deferred, not ported.
#[allow(clippy::too_many_arguments)]
fn dom_check_row_positions(
    this: &mut HPresolve,
    dc: &mut DomCtx<'_>,
    row: HighsInt,
    col: HighsInt,
    direction: HighsInt,
    best_val: f64,
    bound_implied: bool,
    has_cliques: bool,
) -> HPresolveResult {
    if this.row_deleted[row as usize] != 0 {
        return HPresolveResult::kOk;
    }
    let mut check_rowpositions = std::mem::take(&mut dc.check_rowpositions);
    this.get_row_positions(row, &mut check_rowpositions);
    for &pos in check_rowpositions.iter() {
        let p = pos as usize;
        if this.a_value[p] == 0.0 || this.a_row[p] != row {
            continue;
        }
        // get column index
        let k = this.a_col[p];

        // skip column 'col' and deleted columns
        if k == col || this.col_deleted[k as usize] != 0 {
            continue;
        }

        // get coefficient
        let ak = this.a_value[p];

        // check if variables have the same type
        let same_var_type = dom_vars_have_same_type(this, col, k);

        // skip checks if nothing to do
        if !bound_implied && !has_cliques && !same_var_type {
            continue;
        }

        // try to fix variables or strengthen bounds
        if dom_check_domination_nonzero(
            this,
            row,
            direction as f64 * best_val,
            direction as f64 * ak,
        ) {
            hpresolve_checked_call!(dom_check_cols(
                this,
                dc,
                col,
                k,
                direction,
                1,
                bound_implied,
                has_cliques,
                same_var_type,
            ));
            if this.col_deleted[col as usize] != 0 || this.row_deleted[row as usize] != 0 {
                break;
            }
        }

        if this.col_deleted[k as usize] == 0
            && this.a_value[p] != 0.0
            && this.a_row[p] == row
            && dom_check_domination_nonzero(
                this,
                row,
                direction as f64 * best_val,
                -direction as f64 * ak,
            )
        {
            hpresolve_checked_call!(dom_check_cols(
                this,
                dc,
                col,
                k,
                direction,
                -1,
                bound_implied,
                has_cliques,
                same_var_type,
            ));
            if this.col_deleted[col as usize] != 0 || this.row_deleted[row as usize] != 0 {
                break;
            }
        }
    }
    dc.check_rowpositions = check_rowpositions;
    HPresolveResult::kOk
}

/// C++ lambda `checkRow` of dominatedColumns.
#[allow(clippy::too_many_arguments)]
fn dom_check_row(
    this: &mut HPresolve,
    dc: &mut DomCtx<'_>,
    row: HighsInt,
    col: HighsInt,
    direction: HighsInt,
    best_val: f64,
    bound_implied: bool,
    has_cliques: bool,
) -> HPresolveResult {
    // HiGHS 1.16 #3472 (h16_dominated_cols): HPresolveDominatedCols.cpp lines 351-396 (1.16.0)
    if this.options.h16_dominated_cols {
        return dom_check_row_positions(
            this,
            dc,
            row,
            col,
            direction,
            best_val,
            bound_implied,
            has_cliques,
        );
    }
    this.store_row(row);
    let only_pred_bnd_analysis = !bound_implied && !has_cliques;
    let mut it = this.get_stored_row();
    while let Some(nz) = it.next_nz(&*this) {
        // get column index
        let k = nz.index();

        // skip column 'col' and deleted columns
        if k == col || this.col_deleted[k as usize] != 0 {
            continue;
        }

        // get coefficient
        let ak = nz.value();

        // check if variables have the same type
        let same_var_type = dom_vars_have_same_type(this, col, k);

        // skip checks if nothing to do
        if only_pred_bnd_analysis && !same_var_type {
            continue;
        }

        // try to fix variables or strengthen bounds
        // check already known non-zeros in respective columns in advance to
        // avoid (potentially slow) element-wise comparison if possible.
        if dom_check_domination_nonzero(
            this,
            row,
            direction as f64 * best_val,
            direction as f64 * ak,
        ) {
            hpresolve_checked_call!(dom_check_cols(
                this,
                dc,
                col,
                k,
                direction,
                1,
                bound_implied,
                has_cliques,
                same_var_type,
            ));
            if this.col_deleted[col as usize] != 0 {
                break;
            }
        }

        if this.col_deleted[k as usize] == 0
            && dom_check_domination_nonzero(
                this,
                row,
                direction as f64 * best_val,
                -direction as f64 * ak,
            )
        {
            hpresolve_checked_call!(dom_check_cols(
                this,
                dc,
                col,
                k,
                direction,
                -1,
                bound_implied,
                has_cliques,
                same_var_type,
            ));
            if this.col_deleted[col as usize] != 0 {
                break;
            }
        }
    }
    HPresolveResult::kOk
}

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::dominatedColumns (lines 1045-1454, 410 lines)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    pub fn dominated_columns(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::DOMINATED_COLS);
        // See section 6.4 "Dominated columns", Achterberg et al., Presolve Reductions
        // in Mixed Integer Programming, INFORMS Journal on Computing 32(2):473-506.
        // See also Gamrath, G., Koch, T., Martin, A. et al., Progress in presolving
        // for mixed integer programming, Math. Prog. Comp. 7, 367–398 (2015).

        let num_col = self.model.num_col_;
        // non-zero signatures for comparing columns
        let mut dc = DomCtx {
            ps: postsolve_stack,
            signatures: vec![(0u32, 0u32); num_col as usize],
            num_dom_checks: 0,
            num_dom_checks_pred_bnd_analysis: 0,
            num_fixed_cols: 0,
            num_fixed_cols_pred_bnd_analysis: 0,
            num_modified_bnds_pred_bnd_analysis: 0,
            check_rowpositions: Vec::new(),
        };

        for i in 0..self.a_value.len() {
            if self.a_value[i] == 0.0 {
                continue;
            }

            let row = self.a_row[i];
            let col = self.a_col[i];
            let row_lower_finite = (self.model.row_lower_[row as usize] != -kHighsInf) as u32;
            let row_upper_finite = (self.model.row_upper_[row as usize] != kHighsInf) as u32;

            if self.a_value[i] > 0.0 {
                dom_add_signature(&mut dc, row, col, row_lower_finite, row_upper_finite);
            } else {
                dom_add_signature(&mut dc, row, col, row_upper_finite, row_lower_finite);
            }
        }

        // count number of fixed columns and modified bounds
        let mut num_cols: HighsInt = 0;

        // parameters for predictive bound analysis
        const MAX_AVERAGE_NUM_DOM_CHECKS_PRED_BND_ANALYSIS: usize = 10000;
        const MIN_AVERAGE_NUM_REDS_PRED_BND_ANALYSIS: f64 = 1e-2;

        // perform predictive bound analysis?
        let mut allow_pred_bnd_analysis = true;

        for j in 0..num_col {
            let ju = j as usize;
            // skip deleted columns
            if self.col_deleted[ju] != 0 {
                continue;
            }

            // increment counter for number of columns
            num_cols += 1;

            // initialise
            let mut best_row_plus: HighsInt = -1;
            let mut best_row_plus_len: HighsInt = kHighsIInf;
            let mut aj_best_row_plus = 0.0;
            let mut best_row_minus: HighsInt = -1;
            let mut best_row_minus_len: HighsInt = kHighsIInf;
            let mut aj_best_row_minus = 0.0;

            let mut it = self.get_column_vector(j);
            while let Some(nz) = it.next_nz(&*self) {
                let row = nz.index();
                let ri = row as usize;
                let scale: HighsInt = if self.model.row_upper_[ri] != kHighsInf {
                    1
                } else {
                    -1
                };

                let val = scale as f64 * nz.value();
                if val > 0.0 && self.rowsize[ri] < best_row_plus_len {
                    best_row_plus = row;
                    best_row_plus_len = self.rowsize[ri];
                    aj_best_row_plus = nz.value();
                }

                if val < 0.0 && self.rowsize[ri] < best_row_minus_len {
                    best_row_minus = row;
                    best_row_minus_len = self.rowsize[ri];
                    aj_best_row_minus = nz.value();
                }
            }

            // check if bounds are implied or there are cliques
            let lower_implied = self.is_lower_implied(j);
            let upper_implied = self.is_upper_implied(j);
            let has_neg_cliques = dom_is_binary(self, j)
                && self.mip().mipdata_.cliquetable.num_cliques_3(j, false) > 0;
            let has_pos_cliques = dom_is_binary(self, j)
                && self.mip().mipdata_.cliquetable.num_cliques_3(j, true) > 0;

            // use row 'bestRowMinus'
            if best_row_minus != -1 && (allow_pred_bnd_analysis || lower_implied || has_neg_cliques)
            {
                hpresolve_checked_call!(dom_check_row(
                    self,
                    &mut dc,
                    best_row_minus,
                    j,
                    -1,
                    aj_best_row_minus,
                    lower_implied,
                    has_neg_cliques,
                ));
            }

            // use row 'bestRowPlus'
            if self.col_deleted[ju] == 0
                && best_row_plus != -1
                && (allow_pred_bnd_analysis || upper_implied || has_pos_cliques)
            {
                hpresolve_checked_call!(dom_check_row(
                    self,
                    &mut dc,
                    best_row_plus,
                    j,
                    1,
                    aj_best_row_plus,
                    upper_implied,
                    has_pos_cliques,
                ));
            }

            // do not use predictive bound analysis if it requires many domination
            // checks and only yields few fixings or improved bounds on average
            let average_num_dom_checks_pred_bnd_analysis =
                dc.num_dom_checks_pred_bnd_analysis / (num_cols as usize);
            let average_num_reds_pred_bnd_analysis = (dc.num_fixed_cols_pred_bnd_analysis
                + dc.num_modified_bnds_pred_bnd_analysis)
                as f64
                / (num_cols as f64);
            // rusolve mip_presolve_domcol_full: never switch the
            // predictive bound analysis off (HiGHS stops it after 300000 checks if it yields < 1e-2
            // reductions per column, which on set covering models happens before the first fixing)
            allow_pred_bnd_analysis = allow_pred_bnd_analysis
                && (self.options.mip_presolve_domcol_full
                    || dc.num_dom_checks_pred_bnd_analysis
                        <= 30 * MAX_AVERAGE_NUM_DOM_CHECKS_PRED_BND_ANALYSIS
                    || (average_num_dom_checks_pred_bnd_analysis
                        <= MAX_AVERAGE_NUM_DOM_CHECKS_PRED_BND_ANALYSIS
                        && average_num_reds_pred_bnd_analysis
                            >= MIN_AVERAGE_NUM_REDS_PRED_BND_ANALYSIS));
        }

        if dc.num_fixed_cols > 0 || dc.num_modified_bnds_pred_bnd_analysis > 0 {
            crate::highs_log_dev!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Fixed {} dominated columns and strengthened {} bounds\n",
                dc.num_fixed_cols,
                dc.num_modified_bnds_pred_bnd_analysis
            );
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::enumerateSolutions (lines 5223-5696, 474 lines)
    /// HiGHS 1.16 (#3529): HPresolveEnumeration.cpp::enumerateSolutions (lines 16-530)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::enumerateSolutions (lines 5901-6375)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    ///
    /// PORT NOTE: `mipsolver->profiling_` start/stop calls are omitted (profiling is a stub).
    /// The probing scores (`computeProbingScore`, pure, evaluated while no bound has changed)
    /// are computed up front for the integer columns, because the mutable phase runs with the
    /// solver lent to `with_mip_model!` and cannot call `self` methods.
    pub fn enumerate_solutions(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::ENUMERATE);
        // enumerate all solutions for pure binary constraints with a small number of
        // variables

        // prepare probing
        let mut first_call = false;
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 5901-5905 at 7b890c255
        let x1b = self.options.mip_presolve_clique_table;
        if x1b {
            self.mip_mut().mipdata_.cliquetable.set_presolve_flag(true);
        }
        let prepare_result = self.prepare_probing(postsolve_stack, &mut first_call);
        if prepare_result != HPresolveResult::kOk {
            if x1b {
                self.mip_mut().mipdata_.cliquetable.set_presolve_flag(false);
            }
            return prepare_result;
        }

        let num_col = self.model.num_col_ as usize;
        let probing_scores: Vec<(i64, HighsInt)> = (0..num_col)
            .map(|cj| {
                if self.model.integrality_[cj] != HighsVarType::kContinuous {
                    self.compute_probing_score(cj as HighsInt)
                } else {
                    (0, 0)
                }
            })
            .collect();
        let random_seed = self.options.random_seed;
        // HiGHS 1.16 #3529 (h16_enum_sig64): HPresolveEnumeration.cpp lines 31-51, 97-236 (1.16.0)
        let h16_enum_sig64 = self.options.h16_enum_sig64;

        let domain = &self.mip().mipdata_.domains[0];
        let col_lower = domain.col_lower_.clone();
        let col_upper = domain.col_upper_.clone();
        let mut st = EnumState {
            solutions: vec![vec![0; ENUM_MAX_NUM_SOLUTIONS]; ENUM_MAX_ROW_SIZE],
            vars: vec![0; ENUM_MAX_ROW_SIZE],
            branches: vec![(0, 0); ENUM_MAX_ROW_SIZE],
            worst_case_bounds: vec![0; num_col],
            worst_case_lower_bound: vec![kHighsInf; num_col],
            worst_case_upper_bound: vec![-kHighsInf; num_col],
            col_lower,
            col_upper,
            num_cliques_found: 0,
        };

        let result = with_mip_model!(self, ms => enum_run(&mut ms, &probing_scores, random_seed, h16_enum_sig64, &mut st));
        if result != HPresolveResult::kOk {
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 6155 at 7b890c255
            // (every non-OK `enum_run` result comes from `handleInfeasibility`)
            if x1b {
                self.mip_mut().mipdata_.cliquetable.set_presolve_flag(false);
            }
            return result;
        }

        // finalise probing
        let mut num_vars_fixed: HighsInt = 0;
        let mut num_bnds_tightened: HighsInt = 0;
        let mut num_vars_substituted: HighsInt = 0;
        let mut lifted_non_zeros: HighsInt = 0;
        let finalise_result = self.finalise_probing(
            postsolve_stack,
            first_call,
            &mut num_vars_fixed,
            &mut num_bnds_tightened,
            &mut num_vars_substituted,
            &mut lifted_non_zeros,
        );
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 6354-6361 at 7b890c255
        if x1b && finalise_result != HPresolveResult::kOk {
            self.mip_mut().mipdata_.cliquetable.set_presolve_flag(false);
        }
        hpresolve_checked_call!(finalise_result);

        if num_vars_fixed > 0 || num_bnds_tightened > 0 || num_vars_substituted > 0 {
            crate::highs_log_dev!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Enumeration presolve fixed {} columns, tightened {} bounds and performed {} substitutions\n",
                num_vars_fixed,
                num_bnds_tightened,
                num_vars_substituted
            );
        }

        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 6375 at 7b890c255
        if x1b {
            self.mip_mut().mipdata_.cliquetable.set_presolve_flag(false);
        }
        self.check_limits(postsolve_stack)
    }
}

/// C++ `maxRowSize`
const ENUM_MAX_ROW_SIZE: usize = 8;
/// C++ `maxNumRowsChecked`
const ENUM_MAX_NUM_ROWS_CHECKED: HighsInt = 400;
/// C++ `maxNumSolutions = 1 << maxRowSize`
const ENUM_MAX_NUM_SOLUTIONS: usize = 1 << ENUM_MAX_ROW_SIZE;
/// C++ `maxPercentageRowOverlap`
const ENUM_MAX_PERCENTAGE_ROW_OVERLAP: usize = 50;
/// C++ `maxNumFails`
const ENUM_MAX_NUM_FAILS: HighsInt = 6;

/// C++ `std::tuple<double, double, HighsInt, HighsInt, uint32_t> candidateRow`.
type EnumCandidateRow = (f64, f64, HighsInt, HighsInt, u32);

/// Captured mutable locals of `HPresolve::enumerateSolutions` (C++ lambdas capture by reference).
struct EnumState {
    /// C++ `solutions[i][sol]`
    solutions: Vec<Vec<HighsInt>>,
    /// C++ `vars`
    vars: Vec<HighsInt>,
    /// C++ `branches[i]`: (numDomainChanges, numChangedCols)
    branches: Vec<(usize, usize)>,
    /// C++ `worstCaseBounds`
    worst_case_bounds: Vec<HighsInt>,
    /// C++ `worstCaseLowerBound`
    worst_case_lower_bound: Vec<f64>,
    /// C++ `worstCaseUpperBound`
    worst_case_upper_bound: Vec<f64>,
    /// C++ `col_lower` (copy of the domain bounds at entry, updated per row)
    col_lower: Vec<f64>,
    /// C++ `col_upper`
    col_upper: Vec<f64>,
    /// C++ `numCliquesFound` (counted, never read)
    #[allow(dead_code)]
    num_cliques_found: HighsInt,
}

/// C++ lambda `getBinaryRow` of enumerateSolutions. Returns false (and leaves `numnzs`
/// as far as it got) when the row has a non-binary variable or too many variables.
fn enum_get_binary_row(
    ms: &HighsMipSolver,
    row: HighsInt,
    binvars: &mut [HighsInt],
    numnzs: &mut usize,
) -> bool {
    *numnzs = 0;
    let domain = &ms.mipdata_.domains[0];
    let start = ms.mipdata_.a_rstart_[row as usize];
    let end = ms.mipdata_.a_rstart_[row as usize + 1];
    for j in start..end {
        // get index
        let col = ms.mipdata_.a_rindex_[j as usize];
        // skip fixed variables
        if domain.is_fixed(col) {
            continue;
        }
        // skip row if there are non-binary variables or maximum number of
        // elements is reached
        if !domain.is_binary(ms, col) || *numnzs >= ENUM_MAX_ROW_SIZE {
            return false;
        }
        // store binary variable
        binvars[*numnzs] = col;
        *numnzs += 1;
    }
    if *numnzs == 0 {
        return false;
    }
    pdqsort_by(&mut binvars[..*numnzs], |a, b| a < b);
    true
}

/// C++ lambda `computeRowScore` of enumerateSolutions.
fn enum_compute_row_score(
    probing_scores: &[(i64, HighsInt)],
    binvars: &[HighsInt],
    numnzs: usize,
) -> (f64, f64) {
    let mut score: i64 = 0;
    let mut score2: HighsInt = 0;
    for &col in &binvars[..numnzs] {
        let probing_score = probing_scores[col as usize];
        score += probing_score.0;
        score2 += probing_score.1;
    }
    (score as f64 / numnzs as f64, score2 as f64 / numnzs as f64)
}

/// C++ lambda `computeRowSignature` of enumerateSolutions.
fn enum_compute_row_signature(binvars: &[HighsInt], numnzs: usize) -> u32 {
    let mut signature: u32 = 0;
    for &col in &binvars[..numnzs] {
        let col_hashed_pos = hash_i32(col) >> 59;
        debug_assert!(col_hashed_pos < 32);
        signature |= 1u32 << col_hashed_pos;
    }
    signature
}

/// C++ lambda `computeRowOverlap` of enumerateSolutions.
fn enum_compute_row_overlap(
    binvars: &[HighsInt],
    binvars2: &[HighsInt],
    numnzs: usize,
    numnzs2: usize,
) -> usize {
    let mut overlap = 0;
    let mut ir = 0;
    let mut ir2 = 0;
    while ir < numnzs && ir2 < numnzs2 {
        if binvars[ir] < binvars2[ir2] {
            ir += 1;
        } else if binvars[ir] > binvars2[ir2] {
            ir2 += 1;
        } else {
            ir += 1;
            ir2 += 1;
            overlap += 1;
        }
    }
    overlap
}

/// C++ lambda `compileRows` of enumerateSolutions.
fn enum_compile_rows(
    ms: &HighsMipSolver,
    probing_scores: &[(i64, HighsInt)],
    random_seed: HighsInt,
    rows: &mut Vec<EnumCandidateRow>,
) {
    let mut binvars = vec![0; ENUM_MAX_ROW_SIZE];
    let mut numnzs: usize = 0;
    let mut random = HighsRandom::new(random_seed as HighsUInt);
    for i in 0..ms.num_row() {
        // skip redundant rows
        if ms.mipdata_.domains[0].is_redundant_row(ms, i) {
            continue;
        }
        // skip non-binary rows
        if !enum_get_binary_row(ms, i, &mut binvars, &mut numnzs) {
            continue;
        }
        // compute row score
        let score = enum_compute_row_score(probing_scores, &binvars, numnzs);
        // add row to vector
        rows.push((
            -score.0,
            -score.1,
            random.integer(),
            i,
            enum_compute_row_signature(&binvars, numnzs),
        ));
    }
}

/// C++ lambda `removeSimilarRows` of enumerateSolutions.
fn enum_remove_similar_rows(ms: &HighsMipSolver, rows: &mut Vec<EnumCandidateRow>) {
    if rows.len() <= 1 {
        return;
    }
    let mut binvars = vec![0; ENUM_MAX_ROW_SIZE];
    let mut binvars2 = vec![0; ENUM_MAX_ROW_SIZE];
    let mut numnzs: usize = 0;
    let mut numnzs2: usize = 0;
    let mut num_rows_accepted: HighsInt = 0;
    let mut num_rows_removed: HighsInt = 0;
    let mut num_comparisons: usize = 0;
    for i in 0..rows.len() - 1 {
        // get row index and skip removed rows
        let r = rows[i].3;
        if r == -1 {
            continue;
        }
        // check if maximum number of rows is reached
        num_rows_accepted += 1;
        if num_rows_accepted >= ENUM_MAX_NUM_ROWS_CHECKED {
            break;
        }
        // get indices of binary variables in the row
        enum_get_binary_row(ms, r, &mut binvars, &mut numnzs);
        // initialise counters
        let mut num_rows_active: HighsInt = 0;
        let old_num_rows_removed = num_rows_removed;
        for ii in i + 1..rows.len() {
            // get row index and skip removed rows
            let r2 = rows[ii].3;
            if r2 == -1 {
                continue;
            }
            // compare signatures to see if there may be overlap
            num_rows_active += 1;
            if (rows[i].4 & rows[ii].4) == 0 {
                continue;
            }
            // get indices of binary variables in the row
            enum_get_binary_row(ms, r2, &mut binvars2, &mut numnzs2);
            // check if there is too much overlap (compute overlap coefficient)
            num_comparisons += 1;
            let overlap = enum_compute_row_overlap(&binvars, &binvars2, numnzs, numnzs2);
            if (100 * overlap) / cpp_min(numnzs, numnzs2) > ENUM_MAX_PERCENTAGE_ROW_OVERLAP {
                // mark row for removal
                num_rows_removed += 1;
                rows[ii].3 = -1;
            }
        }
        // stop iterating if at most one of the remaining rows is active (not
        // deleted)
        if num_rows_active - num_rows_removed + old_num_rows_removed <= 1 {
            break;
        }
    }
    let _ = num_comparisons;
    if num_rows_removed > 0 {
        rows.retain(|p| p.3 != -1);
    }
    if rows.len() > ENUM_MAX_NUM_ROWS_CHECKED as usize {
        rows.truncate(ENUM_MAX_NUM_ROWS_CHECKED as usize);
    }
}

/// HiGHS 1.16 (#3529): HPresolveEnumeration.cpp::enumerateSolutions `CandidateRow` (lines 31-51).
/// On-path candidate row (`h16_enum_sig64`): the 1.15.1 tuple plus the 64-bit signature, the
/// binary count and the coincidence count used by `mayBeSimilar`. The row index is unique, so
/// `less` compares (score, score2, tieBreaker, row) only, as the C++ `operator<` does.
#[derive(Clone)]
struct EnumCandidateRowH16 {
    /// C++ `score` (negated, so that the best rows come first)
    score: f64,
    /// C++ `score2`
    score2: f64,
    /// C++ `tieBreaker`
    tie_breaker: HighsInt,
    /// C++ `row` (-1 if the row is removed)
    row: HighsInt,
    /// C++ `signature` (uint64_t)
    signature: u64,
    /// C++ `numBinaries`
    num_binaries: HighsInt,
    /// C++ `numCoincidences`: binaries minus the number of distinct signature bits
    num_coincidences: HighsInt,
}

impl EnumCandidateRowH16 {
    /// C++ `CandidateRow::operator<` (lines 46-50).
    fn less(&self, other: &EnumCandidateRowH16) -> bool {
        (self.score, self.score2, self.tie_breaker, self.row)
            < (other.score, other.score2, other.tie_breaker, other.row)
    }
}

/// HiGHS 1.16 (#3529): HPresolveEnumeration.cpp::enumerateSolutions lambda `computeRowSignature`
/// (lines 97-107), 64-bit variant: bucket `hash(col) >> 58`.
/// `HighsHashHelpers::popcnt` (HighsHash.h lines 106/128/177) is `count_ones()` here.
fn enum_compute_row_signature_u64(binvars: &[HighsInt], numnzs: usize) -> u64 {
    let mut signature: u64 = 0;
    for &col in &binvars[..numnzs] {
        let col_hashed_pos = hash_i32(col) >> 58;
        debug_assert!(col_hashed_pos < 64);
        signature |= 1u64 << col_hashed_pos;
    }
    signature
}

/// HiGHS 1.16 (#3529): HPresolveEnumeration.cpp::enumerateSolutions lambda `requiredRowOverlap`
/// (lines 109-114). Kept as the C++ expression (`div_ceil` would give the same value).
#[allow(clippy::manual_div_ceil)]
fn enum_required_row_overlap(numnzs: usize) -> usize {
    ((ENUM_MAX_PERCENTAGE_ROW_OVERLAP + 1) * numnzs + 99) / 100
}

/// HiGHS 1.16 (#3529): HPresolveEnumeration.cpp::enumerateSolutions lambda `mayBeSimilar`
/// (lines 116-129). The signature check never rejects a pair whose real overlap reaches
/// `enum_required_row_overlap`: a shared column sets the same bit in both signatures, so the
/// common bits are at least the overlap minus the coincidences of either row.
fn enum_may_be_similar(row: &EnumCandidateRowH16, row2: &EnumCandidateRowH16) -> bool {
    let num_coincidences = row.num_coincidences.min(row2.num_coincidences);
    let num_common_bits = (row.signature & row2.signature).count_ones() as HighsInt;
    (num_common_bits + num_coincidences) as usize
        >= enum_required_row_overlap(row.num_binaries.min(row2.num_binaries) as usize)
}

/// HiGHS 1.16 (#3529): HPresolveEnumeration.cpp::enumerateSolutions lambda `compileRows`
/// (lines 152-172), on-path row type.
fn enum_compile_rows_h16(
    ms: &HighsMipSolver,
    probing_scores: &[(i64, HighsInt)],
    random_seed: HighsInt,
    rows: &mut Vec<EnumCandidateRowH16>,
) {
    let mut binvars = vec![0; ENUM_MAX_ROW_SIZE];
    let mut numnzs: usize = 0;
    let mut random = HighsRandom::new(random_seed as HighsUInt);
    for i in 0..ms.num_row() {
        // skip redundant rows
        if ms.mipdata_.domains[0].is_redundant_row(ms, i) {
            continue;
        }
        // skip non-binary rows
        if !enum_get_binary_row(ms, i, &mut binvars, &mut numnzs) {
            continue;
        }
        // compute row score
        let score = enum_compute_row_score(probing_scores, &binvars, numnzs);
        // compute signature and add row to vector
        let signature = enum_compute_row_signature_u64(&binvars, numnzs);
        let num_coincidences = numnzs as HighsInt - signature.count_ones() as HighsInt;
        rows.push(EnumCandidateRowH16 {
            score: -score.0,
            score2: -score.1,
            tie_breaker: random.integer(),
            row: i,
            signature,
            num_binaries: numnzs as HighsInt,
            num_coincidences,
        });
    }
}

/// HiGHS 1.16 (#3529): HPresolveEnumeration.cpp::enumerateSolutions lambda `removeSimilarRows`
/// (lines 175-225), on-path row type. The loop order and the overlap test are those of the
/// 1.15.1 `enum_remove_similar_rows`; only the pre-filter `mayBeSimilar` is new.
fn enum_remove_similar_rows_h16(ms: &HighsMipSolver, rows: &mut Vec<EnumCandidateRowH16>) {
    if rows.len() <= 1 {
        return;
    }
    let mut binvars = vec![0; ENUM_MAX_ROW_SIZE];
    let mut binvars2 = vec![0; ENUM_MAX_ROW_SIZE];
    let mut numnzs: usize = 0;
    let mut numnzs2: usize = 0;
    let mut num_rows_accepted: HighsInt = 0;
    let mut num_rows_removed: HighsInt = 0;
    let mut num_comparisons: usize = 0;
    for i in 0..rows.len() - 1 {
        // get row index and skip removed rows
        let r = rows[i].row;
        if r == -1 {
            continue;
        }
        // check if maximum number of rows is reached
        num_rows_accepted += 1;
        if num_rows_accepted >= ENUM_MAX_NUM_ROWS_CHECKED {
            break;
        }
        // get indices of binary variables in the row
        enum_get_binary_row(ms, r, &mut binvars, &mut numnzs);
        // initialise counters
        let mut num_rows_active: HighsInt = 0;
        let old_num_rows_removed = num_rows_removed;
        for ii in i + 1..rows.len() {
            // get row index and skip removed rows
            let r2 = rows[ii].row;
            if r2 == -1 {
                continue;
            }
            // compare signatures to see if there may be too much overlap
            num_rows_active += 1;
            if !enum_may_be_similar(&rows[i], &rows[ii]) {
                continue;
            }
            // get indices of binary variables in the row
            enum_get_binary_row(ms, r2, &mut binvars2, &mut numnzs2);
            // check if there is too much overlap (compute overlap coefficient)
            num_comparisons += 1;
            let overlap = enum_compute_row_overlap(&binvars, &binvars2, numnzs, numnzs2);
            if (100 * overlap) / cpp_min(numnzs, numnzs2) > ENUM_MAX_PERCENTAGE_ROW_OVERLAP {
                // mark row for removal
                num_rows_removed += 1;
                rows[ii].row = -1;
            }
        }
        // stop iterating if at most one of the remaining rows is active (not
        // deleted)
        if num_rows_active - num_rows_removed + old_num_rows_removed <= 1 {
            break;
        }
    }
    let _ = num_comparisons;
    if num_rows_removed > 0 {
        rows.retain(|p| p.row != -1);
    }
    if rows.len() > ENUM_MAX_NUM_ROWS_CHECKED as usize {
        rows.truncate(ENUM_MAX_NUM_ROWS_CHECKED as usize);
    }
}

/// C++ lambda `findBranchVar` of enumerateSolutions.
fn enum_find_branch_var(ms: &HighsMipSolver, st: &EnumState, num_vars: usize) -> HighsInt {
    // find variable for branching
    for i in 0..num_vars {
        if !ms.mipdata_.domains[0].is_fixed(st.vars[i]) {
            return st.vars[i];
        }
    }
    -1
}

/// C++ lambda `doBranch` of enumerateSolutions (just performs the initial lower branch).
fn enum_do_branch(
    ms: &mut HighsMipSolver,
    st: &mut EnumState,
    num_vars: usize,
    num_branches: &mut HighsInt,
) {
    // find variable for branching
    let branchvar = enum_find_branch_var(ms, st, num_vars);
    debug_assert!(branchvar >= 0);

    // branch downwards
    *num_branches += 1;
    st.branches[*num_branches as usize] = (
        ms.mipdata_.domains[0].get_domain_change_stack().len(),
        ms.mipdata_.domains[0].get_changed_cols().len(),
    );
    detach!(ms.mipdata_.domains[0], |d| d.change_bound_2(
        ms,
        HighsBoundType::kUpper,
        branchvar,
        0.0,
        HighsDomain_Reason::branching(),
    ));
}

/// C++ lambda `doBacktrack` of enumerateSolutions.
fn enum_do_backtrack(
    ms: &mut HighsMipSolver,
    st: &mut EnumState,
    num_branches: &mut HighsInt,
) -> bool {
    while *num_branches >= 0 {
        // get column index
        let (num_domain_changes, num_changed_cols) = st.branches[*num_branches as usize];
        let (col, bndtype) = {
            let domchg = &ms.mipdata_.domains[0].get_domain_change_stack()[num_domain_changes];
            (domchg.column, domchg.boundtype)
        };
        // backtrack
        detach!(ms.mipdata_.domains[0], |d| d.backtrack(ms));
        ms.mipdata_.domains[0].clear_changed_cols_2(num_changed_cols);
        if bndtype == HighsBoundType::kUpper {
            // branch upwards
            detach!(ms.mipdata_.domains[0], |d| d.change_bound_2(
                ms,
                HighsBoundType::kLower,
                col,
                1.0,
                HighsDomain_Reason::branching(),
            ));
            break;
        } else {
            // remove branch
            st.branches[*num_branches as usize] = (0, 0);
            *num_branches -= 1;
        }
    }
    // check if enumeration is complete
    *num_branches >= 0
}

/// C++ lambda `identicalVars` of enumerateSolutions.
fn enum_identical_vars(st: &EnumState, num_solutions: usize, index1: usize, index2: usize) -> bool {
    for sol in 0..num_solutions {
        if st.solutions[index1][sol] != st.solutions[index2][sol] {
            return false;
        }
    }
    true
}

/// C++ lambda `complementaryVars` of enumerateSolutions.
fn enum_complementary_vars(
    st: &EnumState,
    num_solutions: usize,
    index1: usize,
    index2: usize,
) -> bool {
    for sol in 0..num_solutions {
        if st.solutions[index1][sol] != 1 - st.solutions[index2][sol] {
            return false;
        }
    }
    true
}

/// C++ lambda `handleInfeasibility` of enumerateSolutions (profiling omitted).
fn enum_handle_infeasibility(infeasible: bool) -> HPresolveResult {
    if infeasible {
        return HPresolveResult::kPrimalInfeasible;
    }
    HPresolveResult::kOk
}

/// C++ lambda `removeWorstCaseBounds` of enumerateSolutions.
fn enum_remove_worst_case_bounds(
    st: &mut EnumState,
    pos: usize,
    num_worst_case_bounds: &mut usize,
) {
    let col = st.worst_case_bounds[pos] as usize;
    st.worst_case_lower_bound[col] = kHighsInf;
    st.worst_case_upper_bound[col] = -kHighsInf;
    st.worst_case_bounds[pos] = st.worst_case_bounds[*num_worst_case_bounds - 1];
    st.worst_case_bounds[*num_worst_case_bounds - 1] = 0;
    *num_worst_case_bounds -= 1;
}

/// C++ lambda `updateWorstCaseBounds` of enumerateSolutions.
fn enum_update_worst_case_bounds(ms: &HighsMipSolver, st: &mut EnumState, col: HighsInt) -> bool {
    let ci = col as usize;
    let domain = &ms.mipdata_.domains[0];
    // update worst-case bounds
    st.worst_case_lower_bound[ci] = cpp_min(st.worst_case_lower_bound[ci], domain.col_lower_[ci]);
    st.worst_case_upper_bound[ci] = cpp_max(st.worst_case_upper_bound[ci], domain.col_upper_[ci]);
    // check if worst-case bounds are not tighter than global bounds
    st.worst_case_lower_bound[ci] <= st.col_lower[ci]
        && st.worst_case_upper_bound[ci] >= st.col_upper[ci]
}

/// C++ lambda `handleSolution` of enumerateSolutions.
#[allow(clippy::too_many_arguments)]
fn enum_handle_solution(
    ms: &mut HighsMipSolver,
    st: &mut EnumState,
    num_vars: usize,
    num_solutions: &mut usize,
    num_worst_case_bounds: &mut usize,
    min_num_active_cols: &mut usize,
    max_num_active_cols: &mut usize,
) {
    // propagate
    detach!(ms.mipdata_.domains[0], |d| d.propagate(ms));
    if ms.mipdata_.domains[0].infeasible() {
        return;
    }
    // handling of worst-case bounds
    if *num_solutions == 0 {
        // initialize
        let changed: Vec<HighsInt> = ms.mipdata_.domains[0].get_changed_cols().to_vec();
        for col in changed {
            st.worst_case_bounds[*num_worst_case_bounds] = col;
            *num_worst_case_bounds += 1;
            enum_update_worst_case_bounds(ms, st, col);
        }
    } else {
        let mut i: usize = 0;
        while i < *num_worst_case_bounds {
            let col = st.worst_case_bounds[i];
            if !ms.mipdata_.domains[0].is_changed_col(col) {
                // no bound changes for this variable -> reset worst-case
                // bounds and remove variable
                enum_remove_worst_case_bounds(st, i, num_worst_case_bounds);
            } else {
                // update worst-case bounds
                if enum_update_worst_case_bounds(ms, st, col) {
                    enum_remove_worst_case_bounds(st, i, num_worst_case_bounds);
                } else {
                    i += 1;
                }
            }
        }
    }
    // store solution and compute minimum and maximum number of active
    // variables (i.e. those having solution value of 1)
    let mut num_active_cols: usize = 0;
    for i in 0..num_vars {
        let sol_value: HighsInt = if ms.mipdata_.domains[0].col_lower_[st.vars[i] as usize] == 0.0 {
            0
        } else {
            1
        };
        st.solutions[i][*num_solutions] = sol_value;
        if sol_value != 0 {
            num_active_cols += 1;
        }
    }
    *min_num_active_cols = cpp_min(*min_num_active_cols, num_active_cols);
    *max_num_active_cols = cpp_max(*max_num_active_cols, num_active_cols);
    *num_solutions += 1;
}

/// The mutable part of `HPresolve::enumerateSolutions`, run with the solver lent to
/// `with_mip_model!` (the presolved model is `ms.mipdata_.presolved_model`).
/// HiGHS 1.16 (#3529): HPresolveEnumeration.cpp::enumerateSolutions (lines 16-530, rows at 232-236, 386-390)
fn enum_run(
    ms: &mut HighsMipSolver,
    probing_scores: &[(i64, HighsInt)],
    random_seed: HighsInt,
    h16_enum_sig64: bool,
    st: &mut EnumState,
) -> HPresolveResult {
    // maximum size of a row and maximum number of rows that will be checked
    // compile rows and sort them
    let mut rows: Vec<EnumCandidateRow> = Vec::new();
    // HiGHS 1.16 #3529 (h16_enum_sig64): HPresolveEnumeration.cpp lines 232-236 (1.16.0)
    let mut rows_h16: Vec<EnumCandidateRowH16> = Vec::new();
    if h16_enum_sig64 {
        enum_compile_rows_h16(ms, probing_scores, random_seed, &mut rows_h16);
        pdqsort_by(&mut rows_h16, |a, b| a.less(b));

        // remove similar rows
        enum_remove_similar_rows_h16(ms, &mut rows_h16);
    } else {
        enum_compile_rows(ms, probing_scores, random_seed, &mut rows);
        pdqsort_by(&mut rows, |a, b| a < b);

        // remove similar rows
        enum_remove_similar_rows(ms, &mut rows);
    }

    // loop over candidate rows
    let num_candidate_rows = if h16_enum_sig64 {
        rows_h16.len()
    } else {
        rows.len()
    };
    let mut num_fails: HighsInt = 0;
    for k in 0..num_candidate_rows {
        // get row index
        let row = if h16_enum_sig64 {
            rows_h16[k].row
        } else {
            rows[k].3
        };
        // skip redundant rows
        if ms.mipdata_.domains[0].is_redundant_row(ms, row) {
            continue;
        }
        // check row
        let mut num_vars: usize = 0;
        if !enum_get_binary_row(ms, row, &mut st.vars, &mut num_vars) {
            continue;
        }

        // main loop
        let mut num_branches: HighsInt = -1;
        let mut num_worst_case_bounds: usize = 0;
        let mut num_solutions: usize = 0;
        let mut min_num_active_cols = num_vars;
        let mut max_num_active_cols: usize = 0;
        loop {
            let mut backtrack = ms.mipdata_.domains[0].infeasible();
            if !backtrack {
                backtrack = enum_find_branch_var(ms, st, num_vars) < 0;
                if backtrack {
                    enum_handle_solution(
                        ms,
                        st,
                        num_vars,
                        &mut num_solutions,
                        &mut num_worst_case_bounds,
                        &mut min_num_active_cols,
                        &mut max_num_active_cols,
                    );
                }
            }
            // branch or backtrack
            if !backtrack {
                enum_do_branch(ms, st, num_vars, &mut num_branches);
            } else if !enum_do_backtrack(ms, st, &mut num_branches) {
                break;
            }
        }

        // no solutions -> infeasible
        hpresolve_checked_call!(enum_handle_infeasibility(num_solutions == 0));

        // store current number of bound changes etc.
        let old_num_changed_cols = ms.mipdata_.domains[0].get_changed_cols().len();
        let old_num_cliques = ms.mipdata_.cliquetable.num_cliques();
        let old_num_substitutions = ms.mipdata_.cliquetable.get_substitutions_2().len();

        // check if all variables form a clique
        if max_num_active_cols == 1 || min_num_active_cols == num_vars - 1 {
            st.num_cliques_found += 1;
            let mut clique: Vec<HighsCliqueTable_CliqueVar> = (0..num_vars)
                .map(|i| {
                    HighsCliqueTable_CliqueVar::new(
                        st.vars[i],
                        (max_num_active_cols == 1) as HighsInt,
                    )
                })
                .collect();
            let equality = min_num_active_cols == max_num_active_cols;
            detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                ms,
                &mut clique,
                equality,
                kHighsIInf
            ));
            hpresolve_checked_call!(enum_handle_infeasibility(
                ms.mipdata_.domains[0].infeasible()
            ));
        }

        // analyse worst-case bounds
        for i in 0..num_worst_case_bounds {
            let col = st.worst_case_bounds[i];
            let ci = col as usize;
            if st.worst_case_lower_bound[ci] > ms.mipdata_.domains[0].col_lower_[ci] {
                // tighten lower bound
                let bound = st.worst_case_lower_bound[ci];
                detach!(ms.mipdata_.domains[0], |d| d.change_bound_2(
                    ms,
                    HighsBoundType::kLower,
                    col,
                    bound,
                    HighsDomain_Reason::unspecified(),
                ));
                hpresolve_checked_call!(enum_handle_infeasibility(
                    ms.mipdata_.domains[0].infeasible()
                ));
            }
            if st.worst_case_upper_bound[ci] < ms.mipdata_.domains[0].col_upper_[ci] {
                // tighten upper bound
                let bound = st.worst_case_upper_bound[ci];
                detach!(ms.mipdata_.domains[0], |d| d.change_bound_2(
                    ms,
                    HighsBoundType::kUpper,
                    col,
                    bound,
                    HighsDomain_Reason::unspecified(),
                ));
                hpresolve_checked_call!(enum_handle_infeasibility(
                    ms.mipdata_.domains[0].infeasible()
                ));
            }
            // clean up
            st.worst_case_lower_bound[ci] = kHighsInf;
            st.worst_case_upper_bound[ci] = -kHighsInf;
            st.worst_case_bounds[i] = 0;
        }

        for i in 0..num_vars - 1 {
            // get column index
            let col = st.vars[i];
            for ii in i + 1..num_vars {
                // get column index
                let col2 = st.vars[ii];
                // check if two binary variables take identical or complementary
                // values in all feasible solutions
                if enum_identical_vars(st, num_solutions, i, ii) {
                    // add clique (1 - x_1) + x_2 = 1 to clique table
                    st.num_cliques_found += 1;
                    let mut clique = [
                        HighsCliqueTable_CliqueVar::new(col, 0),
                        HighsCliqueTable_CliqueVar::new(col2, 1),
                    ];
                    detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                        ms,
                        &mut clique,
                        true,
                        kHighsIInf
                    ));
                    hpresolve_checked_call!(enum_handle_infeasibility(
                        ms.mipdata_.domains[0].infeasible()
                    ));
                } else if enum_complementary_vars(st, num_solutions, i, ii) {
                    // add clique (1 - x_1) + (1 - x_2) = 1 to clique table
                    st.num_cliques_found += 1;
                    let mut clique = [
                        HighsCliqueTable_CliqueVar::new(col, 0),
                        HighsCliqueTable_CliqueVar::new(col2, 0),
                    ];
                    detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                        ms,
                        &mut clique,
                        true,
                        kHighsIInf
                    ));
                    hpresolve_checked_call!(enum_handle_infeasibility(
                        ms.mipdata_.domains[0].infeasible()
                    ));
                }
            }
        }

        // update bounds
        let changed: Vec<HighsInt> = ms.mipdata_.domains[0].get_changed_cols().to_vec();
        let num_changed_cols = changed.len();
        for col in changed {
            let ci = col as usize;
            st.col_lower[ci] = ms.mipdata_.domains[0].col_lower_[ci];
            st.col_upper[ci] = ms.mipdata_.domains[0].col_upper_[ci];
        }

        // clear changed cols
        ms.mipdata_.domains[0].clear_changed_cols();

        // check if solution enumeration failed
        if num_changed_cols != old_num_changed_cols
            || ms.mipdata_.cliquetable.num_cliques() != old_num_cliques
            || ms.mipdata_.cliquetable.get_substitutions_2().len() != old_num_substitutions
        {
            num_fails = 0;
        } else {
            num_fails += 1;
            if num_fails > ENUM_MAX_NUM_FAILS {
                break;
            }
        }
    }

    HPresolveResult::kOk
}

#[cfg(test)]
mod h16_domcol_tests {
    //! HiGHS 1.16 dominated-column regression tests (upstream HiGHS 1.16.0 check/TestMipSolver.cpp models of
    //! #3411, #3290 and #3472). Each runs the public solver end to end.
    use crate::prelude::*;

    /// Solves `lp` with the given presolve setting and HiGHS 1.16 flags; returns the model status and
    /// the objective value.
    fn solve_with(
        lp: HighsLp,
        presolve: &str,
        int_bound_rounding: bool,
        dominated_cols: bool,
        mip_feastol: Option<f64>,
    ) -> (HighsModelStatus, f64) {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_int_bound_rounding", int_bound_rounding),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_dominated_cols", dominated_cols),
            HighsStatus::kOk
        );
        // The Rust option setter does not register mip_feasibility_tolerance (see highs_api_mip.rs
        // option_kind), so the field is set directly.
        if let Some(tol) = mip_feastol {
            highs.options_.mip_feasibility_tolerance = tol;
        }
        assert_eq!(
            highs.set_option_value_str("presolve", presolve),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        (
            highs.get_model_status(),
            highs.get_info().objective_function_value,
        )
    }

    /// Upstream issue-3389 model: tol = 2^-20 and coeff = 2 / tol, so all data is exact. With
    /// y = 1 the rows give 1 <= x <= 2 - tol / 2, and x = 2 - tol / 2 is integral within the MIP
    /// feasibility tolerance.
    fn issue_3389_lp(tol: f64) -> HighsLp {
        let coeff = 2.0 / tol;
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 4;
        lp.sense_ = ObjSense::kMaximize;
        lp.col_cost_ = vec![1.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![3.0, 3.0];
        lp.row_lower_ = vec![-0.5 * coeff, -kHighsInf, 1.0, -kHighsInf];
        lp.row_upper_ = vec![kHighsInf, 0.5 * coeff - 1.0, kHighsInf, 1.0];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kInteger];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 2, 6];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 2, 3];
        lp.a_matrix_.value_ = vec![
            -coeff,
            -coeff,
            1.5 * coeff - 1.0,
            1.5 * coeff - 1.0,
            1.0,
            1.0,
        ];
        lp
    }

    #[test]
    fn h16_domcol_int_bound_rounding_issue_3389() {
        let tol = 1.0 / 1048576.0;
        // Presolve on with the flag off is the 1.15.1 path (it reports Infeasible on this model)
        // and is not asserted here.
        for (presolve, flag) in [("off", false), ("off", true), ("on", true)] {
            let (status, objective) =
                solve_with(issue_3389_lp(tol), presolve, flag, false, Some(tol));
            assert_eq!(
                status,
                HighsModelStatus::kOptimal,
                "presolve {presolve}, flag {flag}"
            );
            assert!(
                (objective - 1.0).abs() < 1e-8 || (objective - (2.0 - 0.5 * tol)).abs() < 1e-8,
                "presolve {presolve}, flag {flag}: objective {objective}"
            );
        }
    }

    /// Upstream issue-3389-binding-row model: with y = 1 row 0 gives x <= 2 - tol / 2 and row 1
    /// gives x <= 2 - tol / 8; each worst-case bound is rounded with its own tolerance.
    fn issue_3389_binding_row_lp(tol: f64) -> HighsLp {
        let coeff = 2.0 / tol;
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 4;
        lp.sense_ = ObjSense::kMaximize;
        lp.col_cost_ = vec![1.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![3.0, 3.0];
        lp.row_lower_ = vec![
            -(1.0 - 0.5 * tol),
            -coeff * (1.0 - 0.125 * tol),
            1.0,
            -kHighsInf,
        ];
        lp.row_upper_ = vec![kHighsInf, kHighsInf, kHighsInf, 1.0];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kInteger];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 2, 6];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 2, 3];
        lp.a_matrix_.value_ = vec![-1.0, -coeff, 1.0, coeff, 1.0, 1.0];
        lp
    }

    #[test]
    fn h16_domcol_int_bound_rounding_issue_3389_binding_row() {
        let tol = 1.0 / 1048576.0;
        for (presolve, flag) in [("off", false), ("off", true), ("on", true)] {
            let (status, objective) = solve_with(
                issue_3389_binding_row_lp(tol),
                presolve,
                flag,
                false,
                Some(tol),
            );
            assert_eq!(
                status,
                HighsModelStatus::kOptimal,
                "presolve {presolve}, flag {flag}"
            );
            assert!(
                (objective - 1.0).abs() < 1e-8 || (objective - (2.0 - 0.5 * tol)).abs() < 1e-8,
                "presolve {presolve}, flag {flag}: objective {objective}"
            );
        }
    }

    /// Upstream dominated-column-double-fixing model (12 binary columns, 10 rows).
    fn dominated_double_fixing_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 12;
        lp.num_row_ = 10;
        lp.col_cost_ = vec![
            -6.0, -4.0, 4.0, 8.0, 1.0, 5.0, -3.0, -3.0, -10.0, 4.0, 0.0, 0.0,
        ];
        lp.col_lower_ = vec![0.0; 12];
        lp.col_upper_ = vec![1.0; 12];
        lp.integrality_ = vec![HighsVarType::kInteger; 12];
        lp.row_lower_ = vec![
            -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf, -2.0, -kHighsInf,
            -kHighsInf, -1.0, -kHighsInf,
        ];
        lp.row_upper_ = vec![3.0, -2.0, 3.0, 1.0, 4.0, -2.0, -1.0, 4.0, -1.0, 1.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 12;
        lp.a_matrix_.num_row_ = 10;
        lp.a_matrix_.start_ = vec![0, 2, 8, 11, 13, 14, 17, 20, 23, 27, 30, 33, 38];
        lp.a_matrix_.index_ = vec![
            5, 9, 1, 3, 4, 5, 6, 8, 1, 4, 5, 5, 8, 2, 1, 6, 7, 2, 7, 8, 1, 4, 8, 1, 3, 4, 9, 4, 5,
            6, 1, 8, 9, 4, 6, 7, 8, 9,
        ];
        lp.a_matrix_.value_ = vec![
            -4.0, 1.0, -1.0, -4.0, -3.0, 1.0, -1.0, 3.0, 4.0, 1.0, -4.0, -1.0, 1.0, -1.0, -2.0,
            1.0, -1.0, -2.0, -1.0, -4.0, -1.0, -4.0, -2.0, -3.0, 1.0, -1.0, -1.0, 4.0, 4.0, -1.0,
            1.0, 4.0, -1.0, -2.0, -1.0, 1.0, 3.0, 1.0,
        ];
        lp
    }

    #[test]
    fn h16_domcol_dominated_double_fixing_issue_3290() {
        // Presolve off never runs dominatedColumns, so both flag values must prove infeasibility.
        for flag in [false, true] {
            let (status, _) = solve_with(dominated_double_fixing_lp(), "off", false, flag, None);
            assert_eq!(
                status,
                HighsModelStatus::kInfeasible,
                "presolve off, flag {flag}"
            );
        }
        // Presolve on with the new dominated-column path (#3472) must prove infeasibility too.
        let (status, _) = solve_with(dominated_double_fixing_lp(), "on", false, true, None);
        assert_eq!(
            status,
            HighsModelStatus::kInfeasible,
            "presolve on, flag on"
        );
    }

    /// Upstream dominated_columns_container_overflow MPS model (12 binary columns, 13 rows,
    /// objective constant 2 from the RHS on the objective row).
    fn container_overflow_lp() -> HighsLp {
        // Columns: t1..t5, d1_1, d1_2, d1_3, d1_5, d2_1, d2_3, d2_5. Each entry is (row, value);
        // rows: 0 total_disk, 1 total_gcu, 2 table_disk, 3 table_gcu, 4 req2_served,
        // 5 req1_served, 6 avail_2_3, 7 avail_1_3, 8 avail_2_1, 9 avail_1_2, 10 avail_1_1,
        // 11 avail_2_5, 12 avail_1_5.
        let columns: [(f64, &[(usize, f64)]); 12] = [
            (
                1.0,
                &[
                    (0, 200.0),
                    (1, 200.0),
                    (2, 200.0),
                    (3, 200.0),
                    (8, -1.0),
                    (10, -1.0),
                ],
            ),
            (
                1.0,
                &[(0, 100.0), (1, 100.0), (2, 100.0), (3, 100.0), (9, -1.0)],
            ),
            (-1.0, &[(0, 200.0), (1, 200.0), (6, -1.0), (7, -1.0)]),
            (-1.0, &[(0, 100.0), (1, 100.0)]),
            (0.0, &[(11, -1.0), (12, -1.0)]),
            (2000.0, &[(5, 1.0), (10, 1.0)]),
            (0.0, &[(5, 1.0), (9, 1.0)]),
            (2000.0, &[(5, 1.0), (7, 1.0)]),
            (20000.0, &[(5, 1.0), (12, 1.0)]),
            (0.0, &[(4, 1.0), (8, 1.0)]),
            (0.0, &[(4, 1.0), (6, 1.0)]),
            (40000.0, &[(4, 1.0), (11, 1.0)]),
        ];
        let mut lp = HighsLp::default();
        lp.num_col_ = 12;
        lp.num_row_ = 13;
        lp.offset_ = 2.0;
        lp.col_cost_ = columns.iter().map(|col| col.0).collect();
        lp.col_lower_ = vec![0.0; 12];
        lp.col_upper_ = vec![1.0; 12];
        lp.integrality_ = vec![HighsVarType::kInteger; 12];
        lp.row_lower_ = vec![-kHighsInf; 13];
        lp.row_upper_ = vec![
            300.000001, 300.000001, 100.000001, 100.000001, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0,
        ];
        lp.row_lower_[4] = 1.0;
        lp.row_lower_[5] = 1.0;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 12;
        lp.a_matrix_.num_row_ = 13;
        // start_ already holds the leading 0 from HighsSparseMatrix::default().
        for (_, entries) in columns.iter() {
            for &(row, value) in entries.iter() {
                lp.a_matrix_.index_.push(row as HighsInt);
                lp.a_matrix_.value_.push(value);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        lp
    }

    #[test]
    fn h16_domcol_container_overflow_issue_3472() {
        // Upstream TestPresolve.cpp dominated-cols-container-overflow (default options, presolve
        // on): Optimal with objective 2. Presolve on with the flag off is the 1.15.1 path (Rust
        // copies the stored row, so it also reaches objective 2) and is not asserted here.
        let (status, objective) = solve_with(container_overflow_lp(), "on", false, true, None);
        assert_eq!(status, HighsModelStatus::kOptimal, "presolve on, flag on");
        assert!(
            (objective - 2.0).abs() < 1e-8,
            "presolve on, flag on: objective {objective}"
        );
    }
}

#[cfg(test)]
mod h16_enum_tests {
    use super::*;
    use crate::lp_data::hconst::kHighsInf;
    use crate::lp_data::highs_callback::HighsCallback;
    use crate::lp_data::highs_lp::HighsLp;
    use crate::lp_data::highs_options::HighsOptions;
    use crate::lp_data::hstruct::HighsSolution;
    use crate::mip::mip_solver::HighsMipSolver;
    use crate::prelude::*;

    /// Pure-binary MIP `min c.x` s.t. `sum_{j in S_r} x_j <= rhs_r`, all columns binary.
    fn enum_mip_lp(num_col: HighsInt, costs: &[f64], rows: &[(Vec<HighsInt>, f64)]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = num_col;
        lp.num_row_ = rows.len() as HighsInt;
        lp.col_cost_ = costs.to_vec();
        lp.col_lower_ = vec![0.0; num_col as usize];
        lp.col_upper_ = vec![1.0; num_col as usize];
        lp.integrality_ = vec![HighsVarType::kInteger; num_col as usize];
        lp.row_lower_ = vec![-kHighsInf; rows.len()];
        lp.row_upper_ = rows.iter().map(|(_, rhs)| *rhs).collect();
        // column-wise matrix: rows are visited in increasing order, so each column's
        // entries are sorted by row index
        let mut col_entries: Vec<Vec<HighsInt>> = vec![Vec::new(); num_col as usize];
        for (r, (cols, _)) in rows.iter().enumerate() {
            for &j in cols {
                col_entries[j as usize].push(r as HighsInt);
            }
        }
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = num_col;
        lp.a_matrix_.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.start_ = vec![0];
        for entries in &col_entries {
            for &r in entries {
                lp.a_matrix_.index_.push(r);
                lp.a_matrix_.value_.push(1.0);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        lp
    }

    /// Ten binaries and fourteen heavily overlapping rows (windows of a 3-step walk mod 10),
    /// each row of size s allowing at most s - 2 ones. Returns the LP and its brute-force optimum.
    fn enum_similar_rows_model() -> (HighsLp, f64) {
        let num_col: HighsInt = 10;
        let costs: Vec<f64> = (0..num_col).map(|j| -(1.0 + 0.37 * j as f64)).collect();
        let rows: Vec<(Vec<HighsInt>, f64)> = (0..14)
            .map(|r| {
                let size = 4 + (r % 3) as HighsInt;
                let cols: Vec<HighsInt> = (0..size)
                    .map(|k| (2 * r as HighsInt + 3 * k) % num_col)
                    .collect();
                (cols, (size - 2) as f64)
            })
            .collect();
        let mut best = f64::INFINITY;
        for mask in 0u32..(1u32 << num_col) {
            let bit = |j: HighsInt| (mask >> j) & 1;
            let feasible = rows.iter().all(|(cols, rhs)| {
                let s: HighsInt = cols.iter().map(|&j| bit(j) as HighsInt).sum();
                (s as f64) <= *rhs
            });
            if !feasible {
                continue;
            }
            let obj: f64 = (0..num_col)
                .map(|j| costs[j as usize] * bit(j) as f64)
                .sum();
            best = best.min(obj);
        }
        (enum_mip_lp(num_col, &costs, &rows), best)
    }

    /// Solves the model through the MIP solver (presolve on, so enumeration runs) with the
    /// `h16_enum_sig64` option set to `on`, and returns the optimal objective.
    fn enum_solve_objective(lp: &HighsLp, on: bool) -> f64 {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.log_options.output_flag = false;
        options.threads = 1;
        options.presolve = "on".to_string();
        options.mip_heuristic_run_feasibility_jump = false;
        options.mip_heuristic_run_fpr = false;
        options.mip_heuristic_run_local_mip = false;
        options.h16_enum_sig64 = on;
        let solution = HighsSolution::default();
        let mut callback = HighsCallback::default();
        let mut ms = HighsMipSolver::new(&mut callback, &options, lp, &solution, false, 0);
        ms.run();
        ms.solution_objective_
    }

    #[test]
    fn h16_enum_sig64_signature_and_coincidences() {
        // Two columns in the same 64-bucket (hash >> 58) give one signature bit and one
        // coincidence; a single column sets the bit at hash >> 58, not at the 1.15.1 >> 59.
        let mut seen: Vec<Option<HighsInt>> = vec![None; 64];
        let mut pair: Option<(HighsInt, HighsInt)> = None;
        for col in 0..200 {
            let bucket = (hash_i32(col) >> 58) as usize;
            if let Some(other) = seen[bucket] {
                pair = Some((other, col));
                break;
            }
            seen[bucket] = Some(col);
        }
        let (c1, c2) = pair.expect("64 buckets: a collision among 200 columns");
        let sig_pair = enum_compute_row_signature_u64(&[c1, c2], 2);
        assert_eq!(sig_pair.count_ones(), 1);
        assert_eq!(2 - sig_pair.count_ones() as HighsInt, 1);
        let sig_single = enum_compute_row_signature_u64(&[c1], 1);
        assert_eq!(sig_single, 1u64 << (hash_i32(c1) >> 58));
        assert_ne!(sig_single, 0);
    }

    #[test]
    fn h16_enum_may_be_similar_never_rejects() {
        assert_eq!(enum_required_row_overlap(1), 1);
        assert_eq!(enum_required_row_overlap(4), 3);
        assert_eq!(enum_required_row_overlap(8), 5);

        // Deterministic sets of 1..=8 columns out of 40; a pair that the 1.15.1 overlap test
        // removes (100 * overlap / min > 50) must never be rejected by the signature filter.
        let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 33
        };
        let mut sets: Vec<Vec<HighsInt>> = Vec::new();
        for _ in 0..60 {
            let size = 1 + (next() % 8) as usize;
            let mut cols: Vec<HighsInt> = Vec::new();
            while cols.len() < size {
                let col = (next() % 40) as HighsInt;
                if !cols.contains(&col) {
                    cols.push(col);
                }
            }
            cols.sort();
            sets.push(cols);
        }
        let rows: Vec<EnumCandidateRowH16> = sets
            .iter()
            .enumerate()
            .map(|(i, cols)| {
                let signature = enum_compute_row_signature_u64(cols, cols.len());
                EnumCandidateRowH16 {
                    score: 0.0,
                    score2: 0.0,
                    tie_breaker: 0,
                    row: i as HighsInt,
                    signature,
                    num_binaries: cols.len() as HighsInt,
                    num_coincidences: cols.len() as HighsInt - signature.count_ones() as HighsInt,
                }
            })
            .collect();
        let mut filtered_pairs = 0usize;
        let mut removable_pairs = 0usize;
        for a in 0..sets.len() {
            for b in a + 1..sets.len() {
                let (sa, sb) = (&sets[a], &sets[b]);
                let overlap = enum_compute_row_overlap(sa, sb, sa.len(), sb.len());
                let removable =
                    (100 * overlap) / sa.len().min(sb.len()) > ENUM_MAX_PERCENTAGE_ROW_OVERLAP;
                let may = enum_may_be_similar(&rows[a], &rows[b]);
                if removable {
                    removable_pairs += 1;
                    assert!(may, "signature filter rejected a removable pair {a},{b}");
                }
                if !may {
                    filtered_pairs += 1;
                }
            }
        }
        assert!(removable_pairs > 0, "no removable pairs in the fixture");
        assert!(filtered_pairs > 0, "the signature filter rejected nothing");
    }

    #[test]
    fn h16_enum_on_same_optimum() {
        // With h16_enum_sig64 on, enumeration runs the 1.16 row filter; the optimum must be the
        // brute-force optimum of the pure-binary model.
        let (lp, optimum) = enum_similar_rows_model();
        let obj = enum_solve_objective(&lp, true);
        assert!(
            (obj - optimum).abs() < 1e-6,
            "on: {obj} vs brute force {optimum}"
        );
    }

    #[test]
    fn h16_enum_off_same_optimum() {
        // Option off: the 1.15.1 path, same model, same brute-force optimum.
        let (lp, optimum) = enum_similar_rows_model();
        let obj = enum_solve_objective(&lp, false);
        assert!(
            (obj - optimum).abs() < 1e-6,
            "off: {obj} vs brute force {optimum}"
        );
    }
}
