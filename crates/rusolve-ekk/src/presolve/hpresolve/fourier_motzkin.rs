//! HiGHS source: presolve/HPresolveFourierMotzkin.cpp (1.16.0, lines 14-861)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiGHS #3255 (`mip_fm_presolve`, default off): Fourier-Motzkin
//! elimination of continuous columns, following Zhang, Ploskas and Sahinidis, "A novel
//! linear optimization presolve technique based on Fourier-Motzkin elimination", MPC 18
//! (2026) 345-378. Read `presolve/hpresolve/mod.rs` first.
//!
//! Translation notes:
//! * The C++ lambdas are `fm_*` methods (those that read the presolve state) or free
//!   functions (those that do not); captured workspace vectors are passed explicitly.
//! * `std::unordered_map` (`rowAncestry`, `rowOriginals`) -> `BTreeMap`. Only lookups,
//!   inserts and erases are order-free; `fourierMotzkinBlockFinalise` iterates
//!   `rowAncestry`, which fixes the order of each parent row's descendant list. That
//!   order only changes the summation order of the dual recovery in postsolve;
//!   presolve output does not depend on it.
//! * `std::set` -> `BTreeSet`.
//! * `startPresolveRuleLog` / `stopPresolveRuleLog` are not called: the port's presolve
//!   log has no slot for the 1.16 rule (as for the zero-cost-singleton rule). The
//!   `analysis_.logging_on_` restore in `finalise` is kept.
//! * `getRowVector(r)` stored in `FmeRowData` is materialised (fme.rs module doc).

use super::*;
use crate::prelude::*;
use crate::presolve::highs_postsolve_stack::{
    FmeAncestryEntry, FmeBlockStep, FmeNewRow, FmeRowData, FmeRowHeader, FmeStepHeader,
};

// max. absolute coefficient
const MAX_COEF: f64 = 1e3;

// max. number of consecutive failures (while trying to build the heap)
const MAX_NUM_FAILS: HighsInt = 100;
// max. size of the heap
const MAX_HEAP_SIZE: HighsInt = 10000;

// sentinel row indices for variable bounds and objective row
const K_UPPER_BOUND_ROW: HighsInt = -2;
const K_LOWER_BOUND_ROW: HighsInt = -3;
const K_OBJECTIVE_ROW: HighsInt = -4;

/// HiGHS 1.16: HPresolveFourierMotzkin.cpp `struct Heap::Entry` (lines 39-43)
#[derive(Clone, Copy, Debug)]
pub(crate) struct FmHeapEntry {
    col: HighsInt,
    ne_red: i64,
    mr_red: i64,
}

/// HiGHS 1.16: HPresolveFourierMotzkin.cpp `struct Heap` (lines 38-130): an indexed
/// max-heap ordered by (neRed, mrRed).
#[derive(Clone, Debug, Default)]
pub(crate) struct FmHeap {
    entries: Vec<FmHeapEntry>,
    pos: Vec<HighsInt>,
}

impl FmHeap {
    pub(crate) fn empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub(crate) fn size(&self) -> HighsInt {
        self.entries.len() as HighsInt
    }

    pub(crate) fn top(&self) -> HighsInt {
        self.entries[0].col
    }

    pub(crate) fn contains(&self, col: HighsInt) -> bool {
        self.pos[col as usize] != -1
    }

    pub(crate) fn reset(&mut self, num_col: HighsInt, reserve_size: HighsInt) {
        self.entries.clear();
        self.entries.reserve(reserve_size as usize);
        self.pos.clear();
        self.pos.resize(num_col as usize, -1);
    }

    pub(crate) fn push(&mut self, col: HighsInt, ne_red: i64, mr_red: i64) {
        self.pos[col as usize] = self.size();
        self.entries.push(FmHeapEntry {
            col,
            ne_red,
            mr_red,
        });
    }

    pub(crate) fn insert(&mut self, col: HighsInt, ne_red: i64, mr_red: i64) {
        self.push(col, ne_red, mr_red);
        self.sift_up(self.pos[col as usize]);
    }

    pub(crate) fn remove(&mut self, col: HighsInt) {
        let p = self.pos[col as usize];
        if p == -1 {
            return;
        }
        self.swap(p, self.size() - 1);
        self.pos[col as usize] = -1;
        self.entries.pop();
        self.sift_up(p);
        self.sift_down(p);
    }

    pub(crate) fn update(&mut self, col: HighsInt, ne_red: i64, mr_red: i64) {
        let p = self.pos[col as usize];
        if p == -1 {
            return;
        }
        self.entries[p as usize].ne_red = ne_red;
        self.entries[p as usize].mr_red = mr_red;
        self.sift_up(p);
        self.sift_down(p);
    }

    pub(crate) fn heapify(&mut self) {
        let mut i = self.size() / 2 - 1;
        while i >= 0 {
            self.sift_down(i);
            i -= 1;
        }
    }

    fn better(&self, i: HighsInt, j: HighsInt) -> bool {
        let (ei, ej) = (&self.entries[i as usize], &self.entries[j as usize]);
        if ei.ne_red != ej.ne_red {
            return ei.ne_red > ej.ne_red;
        }
        ei.mr_red > ej.mr_red
    }

    fn swap(&mut self, i: HighsInt, j: HighsInt) {
        if i == j {
            return;
        }
        self.entries.swap(i as usize, j as usize);
        self.pos[self.entries[i as usize].col as usize] = i;
        self.pos[self.entries[j as usize].col as usize] = j;
    }

    fn sift_up(&mut self, mut i: HighsInt) {
        if i >= self.size() {
            return;
        }
        while i > 0 {
            let parent = (i - 1) / 2;
            if !self.better(i, parent) {
                break;
            }
            self.swap(i, parent);
            i = parent;
        }
    }

    fn sift_down(&mut self, mut i: HighsInt) {
        let n = self.size();
        if i >= n {
            return;
        }
        loop {
            let mut best = i;
            let left = 2 * i + 1;
            let right = 2 * i + 2;
            if left < n && self.better(left, best) {
                best = left;
            }
            if right < n && self.better(right, best) {
                best = right;
            }
            if best == i {
                break;
            }
            self.swap(i, best);
            i = best;
        }
    }
}

/// HiGHS 1.16: HPresolveFourierMotzkin.cpp `struct newRowEntry` (lines 132-135)
#[derive(Clone, Copy, Debug)]
struct NewRowEntry {
    col: HighsInt,
    val: HighsCDouble,
}

/// HiGHS 1.16: HPresolveFourierMotzkin.cpp `struct newRow` (lines 137-144)
#[derive(Clone, Debug)]
struct NewRow {
    entries: Vec<NewRowEntry>,
    upper: f64,
    plus_index: HighsInt,
    minus_index: HighsInt,
    plus_scale: f64,
    minus_scale: f64,
}

/// HiGHS 1.16: HPresolveFourierMotzkin.cpp `struct NewRowOrigin` (lines 146-151)
#[derive(Clone, Copy, Debug)]
struct NewRowOrigin {
    plus_row: HighsInt,
    minus_row: HighsInt,
    plus_scale: f64,
    minus_scale: f64,
}

/// `acceptCoef` lambda (lines 159-162)
fn accept_coef(val: f64) -> bool {
    let absval = val.abs();
    absval == 0.0 || (absval >= 1.0 / MAX_COEF && absval <= MAX_COEF)
}

/// `isReduction` lambda (lines 380-382)
fn is_reduction(ne_red: i64, mr_red: i64) -> bool {
    ne_red > 0 || (ne_red == 0 && mr_red > 0)
}

type RowOriginals = BTreeMap<HighsInt, BTreeSet<HighsInt>>;
type RowAncestry = BTreeMap<HighsInt, Vec<FmeAncestryEntry>>;

/// `insertOriginals` lambda (lines 384-399)
fn insert_originals(
    rows: &mut BTreeSet<HighsInt>,
    originals: &RowOriginals,
    row: HighsInt,
    col: HighsInt,
) {
    if row == K_UPPER_BOUND_ROW {
        rows.insert(-(2 * col + 1));
    } else if row == K_LOWER_BOUND_ROW {
        rows.insert(-(2 * col + 2));
    } else if let Some(set) = originals.get(&row) {
        rows.extend(set.iter().copied());
    } else {
        rows.insert(row);
    }
}

/// `mergeOriginals` lambda (lines 401-408)
fn merge_originals(
    rows: &mut BTreeSet<HighsInt>,
    originals: &RowOriginals,
    plus_row: HighsInt,
    minus_row: HighsInt,
    col: HighsInt,
) {
    rows.clear();
    insert_originals(rows, originals, plus_row, col);
    insert_originals(rows, originals, minus_row, col);
}

/// `cernikovRedundant` lambda (lines 410-417)
pub(crate) fn cernikov_redundant(
    rows: &mut BTreeSet<HighsInt>,
    originals: &RowOriginals,
    plus_row: HighsInt,
    minus_row: HighsInt,
    col: HighsInt,
    num_cols_elim: HighsInt,
) -> bool {
    merge_originals(rows, originals, plus_row, minus_row, col);
    rows.len() as HighsInt > num_cols_elim + 2
}

/// `findRowIndex` lambda (lines 534-539)
fn find_row_index(row: HighsInt, rows: &[FmeRowData]) -> HighsInt {
    for (i, r) in rows.iter().enumerate() {
        if r.header.row == row {
            return i as HighsInt;
        }
    }
    -1
}

/// `inheritAncestry` lambda (lines 551-566)
fn inherit_ancestry(
    row_ancestry: &mut RowAncestry,
    new_model_row: HighsInt,
    parent_row: HighsInt,
    parent_row_index: HighsInt,
    step_index: HighsInt,
    scale: f64,
    is_minus: bool,
) {
    if parent_row < 0 {
        return;
    }
    if let Some(parent) = row_ancestry.get(&parent_row).cloned() {
        for a in &parent {
            row_ancestry
                .entry(new_model_row)
                .or_default()
                .push(FmeAncestryEntry {
                    step: a.step,
                    parent_row_index: a.parent_row_index,
                    scale: a.scale * scale,
                    is_minus: a.is_minus,
                });
        }
    }
    if parent_row_index >= 0 {
        row_ancestry
            .entry(new_model_row)
            .or_default()
            .push(FmeAncestryEntry {
                step: step_index,
                parent_row_index,
                scale,
                is_minus,
            });
    }
}

impl HPresolve {
    /// `isCandidate` lambda (lines 164-175)
    fn fm_is_candidate(&self, col: HighsInt) -> bool {
        let cu = col as usize;
        if self.col_deleted[cu] != 0 {
            return false;
        }
        if self.colsize[cu] == 0 {
            return false;
        }
        if col == self.model.fme_obj_col_ {
            return false;
        }
        if self.model.integrality_[cu] != HighsVarType::kContinuous {
            return false;
        }
        if !accept_coef(self.model.col_cost_[cu]) {
            return false;
        }
        if !self.options.presolve_fm_obj_reformulation && self.model.col_cost_[cu] != 0.0 {
            return false;
        }
        for nz in self.get_column_vector(col).iter(self) {
            if self.is_equation(nz.index()) || !accept_coef(nz.value()) {
                return false;
            }
        }
        true
    }

    /// `computeCandidates` lambda (lines 177-184)
    fn fm_compute_candidates(&self, candidates: &mut Vec<HighsInt>) -> bool {
        candidates.clear();
        for col in 0..self.model.num_col_ {
            if self.fm_is_candidate(col) {
                candidates.push(col);
            }
        }
        let colsize = &self.colsize;
        pdqsort_by(&mut candidates[..], |&a, &b| {
            colsize[a as usize] < colsize[b as usize]
        });
        !candidates.is_empty()
    }

    /// `checkRows` lambda (lines 186-238)
    fn fm_check_rows(
        &self,
        col: HighsInt,
        obj_row_cols: &[HighsInt],
        i_plus: &mut Vec<HighsInt>,
        i_minus: &mut Vec<HighsInt>,
        ne_plus: &mut i64,
        ne_minus: &mut i64,
    ) {
        *ne_plus = 0;
        *ne_minus = 0;
        i_plus.clear();
        i_minus.clear();
        for nz in self.get_column_vector(col).iter(self) {
            let row = nz.index();
            let r = row as usize;
            if self.row_deleted[r] != 0 {
                continue;
            }

            if self.is_ranged(row) {
                i_plus.push(row);
                *ne_plus += self.rowsize[r] as i64;
                i_minus.push(row);
                *ne_minus += self.rowsize[r] as i64;
            } else {
                let direction: HighsInt = if self.model.row_upper_[r] != kHighsInf {
                    1
                } else {
                    -1
                };

                if direction as f64 * nz.value() > 0.0 {
                    i_plus.push(row);
                    *ne_plus += self.rowsize[r] as i64;
                } else {
                    i_minus.push(row);
                    *ne_minus += self.rowsize[r] as i64;
                }
            }
        }

        // include finite variable bounds as singleton rows
        let cu = col as usize;
        if self.model.col_upper_[cu] != kHighsInf {
            i_plus.push(K_UPPER_BOUND_ROW);
            *ne_plus += 1;
        }
        if self.model.col_lower_[cu] != -kHighsInf {
            i_minus.push(K_LOWER_BOUND_ROW);
            *ne_minus += 1;
        }

        // simulate the objective constraint row for candidates with nonzero
        // cost when the reformulation has not yet been performed
        if !obj_row_cols.is_empty() && self.model.col_cost_[cu] != 0.0 {
            let obj_row_size = obj_row_cols.len() as i64;
            if self.model.col_cost_[cu] > 0.0 {
                i_plus.push(K_OBJECTIVE_ROW);
                *ne_plus += obj_row_size;
            } else {
                i_minus.push(K_OBJECTIVE_ROW);
                *ne_minus += obj_row_size;
            }
        }
    }

    /// `collectAffectedCols` lambda (lines 240-262)
    fn fm_collect_affected_cols(
        &self,
        col: HighsInt,
        set: &[HighsInt],
        obj_row_cols: &[HighsInt],
        mark: &mut [HighsInt],
        other_mark: &[HighsInt],
        affected_cols: &mut Vec<HighsInt>,
    ) {
        for &row in set {
            if row == K_OBJECTIVE_ROW {
                for &k in obj_row_cols {
                    if k == col {
                        continue;
                    }
                    let ku = k as usize;
                    if mark[ku] == 0 && other_mark[ku] == 0 {
                        affected_cols.push(k);
                    }
                    mark[ku] += 1;
                }
            } else {
                if row < 0 {
                    continue;
                }
                for nz in self.get_row_vector(row).iter(self) {
                    let k = nz.index();
                    if k == col {
                        continue;
                    }
                    let ku = k as usize;
                    if mark[ku] == 0 && other_mark[ku] == 0 {
                        affected_cols.push(k);
                    }
                    mark[ku] += 1;
                }
            }
        }
    }

    /// `checkNonZeros` lambda (lines 264-311)
    #[allow(clippy::too_many_arguments)]
    fn fm_check_non_zeros(
        &self,
        col: HighsInt,
        obj_row_cols: &[HighsInt],
        i_plus: &mut Vec<HighsInt>,
        i_minus: &mut Vec<HighsInt>,
        p_plus: &mut [HighsInt],
        p_minus: &mut [HighsInt],
        affected_cols: &mut Vec<HighsInt>,
        ne_red: &mut i64,
        mr_red: &mut i64,
    ) -> bool {
        // initialise
        *ne_red = 0;
        *mr_red = 0;

        // check rows
        let mut ne_plus: i64 = 0;
        let mut ne_minus: i64 = 0;
        self.fm_check_rows(
            col,
            obj_row_cols,
            i_plus,
            i_minus,
            &mut ne_plus,
            &mut ne_minus,
        );

        if i_plus.is_empty() || i_minus.is_empty() {
            // other presolve reductions may handle this case (e.g., implied free
            // column substitution)
            i_plus.clear();
            i_minus.clear();
            return false;
        }

        // take into account other variables present in the rows
        self.fm_collect_affected_cols(col, i_plus, obj_row_cols, p_plus, p_minus, affected_cols);
        self.fm_collect_affected_cols(col, i_minus, obj_row_cols, p_minus, p_plus, affected_cols);

        // compute correction term
        let mut correction: i64 = 0;
        for &k in affected_cols.iter() {
            let ku = k as usize;
            correction += p_plus[ku] as i64 * p_minus[ku] as i64;
            p_plus[ku] = 0;
            p_minus[ku] = 0;
        }

        let m_plus = i_plus.len() as i64;
        let m_minus = i_minus.len() as i64;
        let ne_old = ne_plus + ne_minus;
        // note that we subtract the entries for column 'col' since these are
        // eliminated
        let ne_new = m_plus * (ne_minus - m_minus) + m_minus * (ne_plus - m_plus) - correction;
        *ne_red = ne_old - ne_new;
        *mr_red = m_plus + m_minus - m_plus * m_minus;
        true
    }

    /// `checkNewRow` lambda (lines 313-346)
    fn fm_check_new_row(&self, nr: &NewRow, is_redundant: &mut bool) -> HPresolveResult {
        let mut implied_lower = HighsCDouble::from(0.0);
        let mut implied_upper = HighsCDouble::from(0.0);
        let mut lower_finite = true;
        let mut upper_finite = true;
        *is_redundant = false;
        for e in &nr.entries {
            let lb = self.model.col_lower_[e.col as usize];
            let ub = self.model.col_upper_[e.col as usize];
            if e.val > 0.0 {
                lower_finite = lower_finite && lb != -kHighsInf;
                if lower_finite {
                    implied_lower += e.val * lb;
                }
                upper_finite = upper_finite && ub != kHighsInf;
                if upper_finite {
                    implied_upper += e.val * ub;
                }
            } else {
                lower_finite = lower_finite && ub != kHighsInf;
                if lower_finite {
                    implied_lower += e.val * ub;
                }
                upper_finite = upper_finite && lb != -kHighsInf;
                if upper_finite {
                    implied_upper += e.val * lb;
                }
            }
            if !lower_finite && !upper_finite {
                return HPresolveResult::kOk;
            }
        }

        let lower = if lower_finite {
            f64::from(implied_lower)
        } else {
            -kHighsInf
        };
        let upper = if upper_finite {
            f64::from(implied_upper)
        } else {
            kHighsInf
        };

        // check for infeasibility
        if lower > nr.upper + self.primal_feastol {
            return HPresolveResult::kPrimalInfeasible;
        }

        // check for redundancy
        *is_redundant = upper <= nr.upper + self.primal_feastol;

        HPresolveResult::kOk
    }

    /// `getRowData` lambda (lines 348-362): returns (absCoef, direction, bound).
    /// `findNonzero` splays the row, as in C++.
    fn fm_get_row_data(
        &mut self,
        row: HighsInt,
        col: HighsInt,
        multiplier: HighsInt,
    ) -> (f64, HighsInt, f64) {
        if row < 0 {
            // artificial lower / upper bound row
            let bound = if multiplier > 0 {
                self.model.col_upper_[col as usize]
            } else {
                -self.model.col_lower_[col as usize]
            };
            (1.0, 1, bound)
        } else {
            let p_pos = self.find_nonzero(row, col);
            debug_assert!(p_pos != -1);
            let aval = self.a_value[p_pos as usize];
            let direction: HighsInt = if multiplier as f64 * aval > 0.0 {
                1
            } else {
                -1
            };
            let abs_coef = aval.abs();
            let bound = if direction > 0 {
                self.model.row_upper_[row as usize]
            } else {
                -self.model.row_lower_[row as usize]
            };
            (abs_coef, direction, bound)
        }
    }

    /// `collectRowEntries` lambda (lines 364-378)
    fn fm_collect_row_entries(
        &self,
        row: HighsInt,
        col: HighsInt,
        scale: f64,
        new_row_entries: &mut Vec<NewRowEntry>,
        new_row_mark: &mut [HighsInt],
    ) {
        if row < 0 {
            return;
        }
        for nz in self.get_row_vector(row).iter(self) {
            if nz.index() == col {
                continue;
            }
            let val = scale * nz.value();
            let k = nz.index() as usize;
            if new_row_mark[k] == -1 {
                new_row_mark[k] = new_row_entries.len() as HighsInt;
                new_row_entries.push(NewRowEntry {
                    col: nz.index(),
                    val: HighsCDouble::from(val),
                });
            } else {
                new_row_entries[new_row_mark[k] as usize].val += val;
            }
        }
    }

    /// `reformulateObjective` lambda (lines 419-491): min c^T x + offset becomes min z
    /// with c^T x - z <= -offset, so that FM can eliminate continuous columns with
    /// nonzero cost.
    fn fm_reformulate_objective(&mut self, postsolve_stack: &mut HighsPostsolveStack) {
        if self.model.fme_obj_col_ != -1 {
            debug_assert!(self.col_deleted[self.model.fme_obj_col_ as usize] == 0);
            return;
        }

        let z_col = self.model.num_col_;
        self.model.num_col_ += 1;
        self.model.a_matrix_.num_col_ += 1;

        // extend model vectors
        self.model.col_cost_.push(1.0);
        self.model.col_lower_.push(-kHighsInf);
        self.model.col_upper_.push(kHighsInf);
        self.model.integrality_.push(HighsVarType::kContinuous);
        let last_start = *self.model.a_matrix_.start_.last().unwrap();
        self.model.a_matrix_.start_.push(last_start);
        if !self.model.col_names_.is_empty() {
            self.model.col_names_.push("fme_obj_z".to_string());
        }

        // extend presolve vectors
        self.colhead.push(-1);
        self.colsize.push(0);
        self.col_deleted.push(0);
        self.impl_col_lower.push(-kHighsInf);
        self.impl_col_upper.push(kHighsInf);
        self.col_lower_source.push(-1);
        self.col_upper_source.push(-1);
        self.impl_row_dual_source_by_col.push(BTreeSet::new());
        self.changed_col_flag.push(1);
        self.num_probes.push(0);

        // update implied bound structures (setBoundArrays is a no-op in the port)
        self.implied_row_bounds.set_bound_arrays();
        self.implied_dual_row_bounds
            .set_num_sums(self.model.num_col_);

        // register in postsolve stack
        postsolve_stack.append_col_to_model();

        // build the objective constraint row: c^T x - z <= -offset
        let offset = self.model.offset_;
        let mut obj_indices: Vec<HighsInt> = Vec::new();
        let mut obj_values: Vec<f64> = Vec::new();
        for j in 0..z_col {
            let ju = j as usize;
            if self.col_deleted[ju] == 0 && self.model.col_cost_[ju] != 0.0 {
                obj_indices.push(j);
                obj_values.push(self.model.col_cost_[ju]);
            }
        }
        obj_indices.push(z_col);
        obj_values.push(-1.0);

        // zero out original costs and offset
        for j in 0..z_col as usize {
            self.model.col_cost_[j] = 0.0;
        }
        self.model.offset_ = 0.0;

        // add the constraint row to the matrix (C++ ignores the result)
        let _ = self.add_to_matrix_row(
            postsolve_stack,
            -kHighsInf,
            -offset,
            &obj_indices,
            &obj_values,
        );

        // register reduction so getReducedPrimalSolution can compute z
        let cost_entries: Vec<Nonzero> = obj_indices
            .iter()
            .zip(obj_values.iter())
            .map(|(&i, &v)| Nonzero::new(i, v))
            .collect();
        postsolve_stack.fourier_motzkin_obj_col(z_col, offset, &cost_entries);

        self.model.fme_obj_col_ = z_col;

        self.shrink_problem(postsolve_stack);
    }

    /// `collectCandidatesAndBuildHeap` lambda (lines 493-531)
    #[allow(clippy::too_many_arguments)]
    fn fm_collect_candidates_and_build_heap(
        &self,
        candidates: &mut Vec<HighsInt>,
        heap: &mut FmHeap,
        i_plus: &mut Vec<HighsInt>,
        i_minus: &mut Vec<HighsInt>,
        p_plus: &mut Vec<HighsInt>,
        p_minus: &mut Vec<HighsInt>,
        affected_cols: &mut Vec<HighsInt>,
        obj_row_cols: &[HighsInt],
    ) -> bool {
        // compute candidates
        if !self.fm_compute_candidates(candidates) {
            return false;
        }
        // set up data structures for heap
        let num_col = self.model.num_col_ as usize;
        heap.reset(self.model.num_col_, candidates.len() as HighsInt);
        p_plus.clear();
        p_plus.resize(num_col, 0);
        p_minus.clear();
        p_minus.resize(num_col, 0);
        i_plus.reserve((self.model.num_row_ as usize).saturating_sub(i_plus.len()));
        i_minus.reserve((self.model.num_row_ as usize).saturating_sub(i_minus.len()));
        affected_cols.reserve(num_col.saturating_sub(affected_cols.len()));
        // inspect candidates (with limits)
        let mut num_fails: HighsInt = 0;
        for &col in candidates.iter() {
            let mut ne_red: i64 = 0;
            let mut mr_red: i64 = 0;
            let elim_candidate = self.fm_check_non_zeros(
                col,
                obj_row_cols,
                i_plus,
                i_minus,
                p_plus,
                p_minus,
                affected_cols,
                &mut ne_red,
                &mut mr_red,
            );
            affected_cols.clear();
            if !elim_candidate || !is_reduction(ne_red, mr_red) {
                // count number of failures
                num_fails += 1;
                if num_fails > MAX_NUM_FAILS {
                    break;
                }
                continue;
            }
            // add to heap
            num_fails = 0;
            heap.push(col, ne_red, mr_red);
            if heap.size() >= MAX_HEAP_SIZE {
                break;
            }
        }
        if heap.empty() {
            return false;
        }
        // heapify
        heap.heapify();
        true
    }

    /// `collectRows` lambda (lines 541-549)
    fn fm_collect_rows(&self, rows: &[HighsInt]) -> Vec<FmeRowData> {
        let mut result: Vec<FmeRowData> = Vec::new();
        for &r in rows {
            if r < 0 {
                continue;
            }
            result.push(FmeRowData {
                header: FmeRowHeader {
                    row: r,
                    row_lower: self.model.row_lower_[r as usize],
                    row_upper: self.model.row_upper_[r as usize],
                },
                row_vec: self.get_row_vector(r).iter(self).collect(),
            });
        }
        result
    }

    /// `printLog` lambda (lines 568-576)
    fn fm_print_log(&self, cols_removed: HighsInt, rows_removed: HighsInt, rows_added: HighsInt) {
        crate::highs_log_dev!(
            &self.options.log_options,
            HighsLogType::kInfo,
            "Fourier-Motzkin ({} objective reformulation) added {} rows and eliminated {} rows and {} columns\n",
            if self.options.presolve_fm_obj_reformulation {
                "with"
            } else {
                "without"
            },
            rows_added,
            rows_removed,
            cols_removed
        );
    }

    /// `finalise` lambda (lines 153-157)
    fn fm_finalise(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        logging_on: bool,
    ) -> HPresolveResult {
        self.analysis_.logging_on_ = logging_on;
        self.check_limits(postsolve_stack)
    }

    /// HiGHS 1.16: HPresolveFourierMotzkin.cpp::HPresolve::fourierMotzkin (lines 14-861)
    pub fn fourier_motzkin(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        num_cols_eliminated: &mut HighsInt,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::FOURIER_MOTZKIN);
        let logging_on = self.analysis_.logging_on_;

        // workspace vectors
        let mut candidates: Vec<HighsInt> = Vec::new();
        let mut i_plus: Vec<HighsInt> = Vec::new();
        let mut i_minus: Vec<HighsInt> = Vec::new();
        let mut p_plus: Vec<HighsInt> = Vec::new();
        let mut p_minus: Vec<HighsInt> = Vec::new();
        let mut affected_cols: Vec<HighsInt> = Vec::new();

        // indexed max-heap
        let mut heap = FmHeap::default();

        // precompute the objective row: columns with nonzero cost
        // used to simulate the objective constraint in checkRows before
        // reformulation actually happens
        let mut obj_row_cols: Vec<HighsInt> = Vec::new();
        if self.options.presolve_fm_obj_reformulation && self.model.fme_obj_col_ == -1 {
            for j in 0..self.model.num_col_ {
                let ju = j as usize;
                if self.col_deleted[ju] == 0 && self.model.col_cost_[ju] != 0.0 {
                    obj_row_cols.push(j);
                }
            }
        }

        // compute candidates and build initial heap
        if !self.fm_collect_candidates_and_build_heap(
            &mut candidates,
            &mut heap,
            &mut i_plus,
            &mut i_minus,
            &mut p_plus,
            &mut p_minus,
            &mut affected_cols,
            &obj_row_cols,
        ) {
            return self.fm_finalise(postsolve_stack, logging_on);
        }

        // vectors for computing new row entries
        let mut new_row_entries: Vec<NewRowEntry> = Vec::new();
        let mut new_row_mark: Vec<HighsInt> = vec![-1; self.model.num_col_ as usize];

        // vector for storing new rows
        let mut new_rows: Vec<NewRow> = Vec::new();

        // workspace for filtering new rows
        let mut row_lower: Vec<f64> = Vec::new();
        let mut row_upper: Vec<f64> = Vec::new();
        let mut row_indices: Vec<Vec<HighsInt>> = Vec::new();
        let mut row_values: Vec<Vec<f64>> = Vec::new();
        let mut new_row_origins: Vec<NewRowOrigin> = Vec::new();

        // vector for saving affected candidates
        let mut save_affected_cols: Vec<HighsInt> = Vec::new();

        // counters for numbers of eliminations
        *num_cols_eliminated = 0;
        let mut num_cols_eliminated_block: HighsInt = 0;
        let mut num_rows_eliminated: HighsInt = 0;
        let mut num_rows_added: HighsInt = 0;

        // FM block data for postsolve
        let mut block_steps: Vec<FmeBlockStep> = Vec::new();

        // surviving row to its ancestry (which parent rows it descends from)
        let mut row_ancestry: RowAncestry = BTreeMap::new();

        // distinct original parent rows for each derived row (Cernikov check)
        let mut row_originals: RowOriginals = BTreeMap::new();
        let mut merged_originals: BTreeSet<HighsInt> = BTreeSet::new();

        // main loop: eliminate variables from heap
        while !heap.empty() {
            let col = heap.top();
            heap.remove(col);

            // if this candidate has nonzero cost and objective has not yet been
            // reformulated, perform the reformulation now and rebuild the heap
            if self.model.fme_obj_col_ == -1 && self.model.col_cost_[col as usize] != 0.0 {
                // finalise any in-progress FM block before reformulating, since
                // reformulateObjective pushes other reductions onto the data stack
                if !block_steps.is_empty() {
                    postsolve_stack.fourier_motzkin_block_finalise(&block_steps, &row_ancestry);
                    self.fm_print_log(
                        num_cols_eliminated_block,
                        num_rows_eliminated,
                        num_rows_added,
                    );
                    block_steps.clear();
                    row_ancestry.clear();
                    row_originals.clear();
                    num_cols_eliminated_block = 0;
                    num_rows_eliminated = 0;
                    num_rows_added = 0;
                }
                // reformulate objective
                self.fm_reformulate_objective(postsolve_stack);
                // clear vector for objective and resize marker
                obj_row_cols.clear();
                new_row_mark.resize(self.model.num_col_ as usize, -1);
                // re-compute candidates and re-build heap
                if !self.fm_collect_candidates_and_build_heap(
                    &mut candidates,
                    &mut heap,
                    &mut i_plus,
                    &mut i_minus,
                    &mut p_plus,
                    &mut p_minus,
                    &mut affected_cols,
                    &obj_row_cols,
                ) {
                    return self.fm_finalise(postsolve_stack, logging_on);
                }
                continue;
            }

            // compute affected columns
            let mut ne_red: i64 = 0;
            let mut mr_red: i64 = 0;
            let elim_candidate = self.fm_check_non_zeros(
                col,
                &obj_row_cols,
                &mut i_plus,
                &mut i_minus,
                &mut p_plus,
                &mut p_minus,
                &mut affected_cols,
                &mut ne_red,
                &mut mr_red,
            );

            // heap data should be up-to-date
            debug_assert!(elim_candidate && is_reduction(ne_red, mr_red));
            let _ = elim_candidate;

            let step_idx = block_steps.len() as HighsInt;

            // perform elimination: generate new rows
            new_rows.clear();
            for &p_row in i_plus.iter() {
                let (p_coef_abs, p_direction, p_bound) = self.fm_get_row_data(p_row, col, 1);

                for &m_row in i_minus.iter() {
                    let (m_coef_abs, m_direction, m_bound) = self.fm_get_row_data(m_row, col, -1);

                    // scale factor to preserve violation tolerances (see section 4.3):
                    let s = (p_coef_abs * m_coef_abs) / (p_coef_abs + m_coef_abs);
                    let p_scale = s / p_coef_abs;
                    let m_scale = s / m_coef_abs;

                    // collect row entries
                    self.fm_collect_row_entries(
                        p_row,
                        col,
                        p_direction as f64 * p_scale,
                        &mut new_row_entries,
                        &mut new_row_mark,
                    );
                    self.fm_collect_row_entries(
                        m_row,
                        col,
                        m_direction as f64 * m_scale,
                        &mut new_row_entries,
                        &mut new_row_mark,
                    );

                    // reset marker before removing near-zeros
                    for e in &new_row_entries {
                        new_row_mark[e.col as usize] = -1;
                    }

                    // remove near-zero entries
                    let small_matrix_value = self.options.small_matrix_value;
                    new_row_entries.retain(|e| !(e.val.abs() <= small_matrix_value));

                    // store new row
                    let new_upper = f64::from(
                        HighsCDouble::from(p_scale) * p_bound
                            + HighsCDouble::from(m_scale) * m_bound,
                    );
                    new_rows.push(NewRow {
                        entries: new_row_entries.clone(),
                        upper: new_upper,
                        plus_index: p_row,
                        minus_index: m_row,
                        plus_scale: p_direction as f64 * p_scale,
                        minus_scale: m_direction as f64 * m_scale,
                    });

                    // clear vector
                    new_row_entries.clear();
                }
            }

            // add new rows, filtering out redundant ones
            row_lower.clear();
            row_upper.clear();
            row_indices.clear();
            row_values.clear();
            new_row_origins.clear();

            for nr in &new_rows {
                let mut redundant = false;
                hpresolve_checked_call!(self.fm_check_new_row(nr, &mut redundant));
                if redundant {
                    continue;
                }

                // Cernikov redundancy check
                if cernikov_redundant(
                    &mut merged_originals,
                    &row_originals,
                    nr.plus_index,
                    nr.minus_index,
                    col,
                    *num_cols_eliminated,
                ) {
                    continue;
                }

                let mut indices: Vec<HighsInt> = Vec::with_capacity(nr.entries.len());
                let mut values: Vec<f64> = Vec::with_capacity(nr.entries.len());
                for e in &nr.entries {
                    indices.push(e.col);
                    values.push(f64::from(e.val));
                }
                row_lower.push(-kHighsInf);
                row_upper.push(nr.upper);
                row_indices.push(indices);
                row_values.push(values);
                new_row_origins.push(NewRowOrigin {
                    plus_row: nr.plus_index,
                    minus_row: nr.minus_index,
                    plus_scale: nr.plus_scale,
                    minus_scale: nr.minus_scale,
                });
            }

            // serialize row data for postsolve before addToMatrix invalidates slices
            let plus_rows = self.fm_collect_rows(&i_plus);
            let minus_rows = self.fm_collect_rows(&i_minus);

            // push row data for this elimination step onto the postsolve stack
            postsolve_stack.fourier_motzkin_block_push_step(col, &plus_rows, &minus_rows);

            // save block metadata
            debug_assert!(self.model.col_cost_[col as usize] == 0.0);
            block_steps.push(FmeBlockStep {
                header: FmeStepHeader {
                    col_lower: self.model.col_lower_[col as usize],
                    col_upper: self.model.col_upper_[col as usize],
                    col,
                    num_plus: plus_rows.len() as HighsInt,
                    num_minus: minus_rows.len() as HighsInt,
                },
                new_rows: Vec::new(),
            });

            // add new rows to matrix
            let first_new_row = self.model.num_row_;
            if !self.add_to_matrix_rows(
                postsolve_stack,
                &row_lower,
                &row_upper,
                &row_indices,
                &row_values,
            ) {
                return self.fm_finalise(postsolve_stack, logging_on);
            }
            num_rows_added += row_indices.len() as HighsInt;

            // build FmeNewRow data and ancestry for this step
            let mut step_new_rows: Vec<FmeNewRow> = Vec::with_capacity(new_row_origins.len());
            for (k, origin) in new_row_origins.iter().enumerate() {
                let new_model_row = first_new_row + k as HighsInt;
                let p_idx = find_row_index(origin.plus_row, &plus_rows);
                let m_idx = find_row_index(origin.minus_row, &minus_rows);
                inherit_ancestry(
                    &mut row_ancestry,
                    new_model_row,
                    origin.plus_row,
                    p_idx,
                    step_idx,
                    origin.plus_scale,
                    false,
                );
                inherit_ancestry(
                    &mut row_ancestry,
                    new_model_row,
                    origin.minus_row,
                    m_idx,
                    step_idx,
                    origin.minus_scale,
                    true,
                );
                merge_originals(
                    &mut merged_originals,
                    &row_originals,
                    origin.plus_row,
                    origin.minus_row,
                    col,
                );
                row_originals.insert(new_model_row, merged_originals.clone());
                step_new_rows.push(FmeNewRow {
                    row: new_model_row,
                    plus_parent_idx: p_idx,
                    minus_parent_idx: m_idx,
                });
            }
            block_steps.last_mut().unwrap().new_rows = step_new_rows;

            // mark column as deleted
            self.mark_col_deleted(col);
            num_cols_eliminated_block += 1;
            *num_cols_eliminated += 1;

            // remove old rows containing col (skip bound rows)
            for &rp in i_plus.iter() {
                if rp < 0 {
                    continue;
                }
                row_ancestry.remove(&rp);
                row_originals.remove(&rp);
                self.remove_row(rp);
                num_rows_eliminated += 1;
            }
            for &rm in i_minus.iter() {
                if rm < 0 {
                    continue;
                }
                row_ancestry.remove(&rm);
                row_originals.remove(&rm);
                if self.row_deleted[rm as usize] != 0 {
                    continue;
                }
                self.remove_row(rm);
                num_rows_eliminated += 1;
            }

            // update affected candidates in the heap
            std::mem::swap(&mut save_affected_cols, &mut affected_cols);
            for &k in save_affected_cols.iter() {
                // check if variable is a candidate
                let is_candidate_col = self.fm_is_candidate(k);
                // skip variable if it is not on the heap and no candidate
                if !heap.contains(k) && !is_candidate_col {
                    continue;
                }
                // check column non-zeros
                let mut ne: i64 = 0;
                let mut mr: i64 = 0;
                let elim_candidate = is_candidate_col
                    && self.fm_check_non_zeros(
                        k,
                        &obj_row_cols,
                        &mut i_plus,
                        &mut i_minus,
                        &mut p_plus,
                        &mut p_minus,
                        &mut affected_cols,
                        &mut ne,
                        &mut mr,
                    );
                affected_cols.clear();
                if !elim_candidate || !is_reduction(ne, mr) {
                    // no candidate or not beneficial -> remove from heap
                    heap.remove(k);
                } else if !heap.contains(k) {
                    // new candidate -> insert into heap
                    heap.insert(k, ne, mr);
                } else {
                    // update heap
                    heap.update(k, ne, mr);
                }
            }
            save_affected_cols.clear();

            if self.check_limits(postsolve_stack) != HPresolveResult::kOk {
                break;
            }
        }

        if num_cols_eliminated_block > 0 {
            // finalize the FM block
            postsolve_stack.fourier_motzkin_block_finalise(&block_steps, &row_ancestry);

            // log message
            self.fm_print_log(
                num_cols_eliminated_block,
                num_rows_eliminated,
                num_rows_added,
            );
        }

        self.fm_finalise(postsolve_stack, logging_on)
    }
}

#[cfg(test)]
mod x1b_fm_tests {
    use super::*;
    use crate::lp_data::highs_options::{kIpmString, kSimplexString};
    use crate::util::highs_timer::HighsTimer;

    /// The heap pops by (neRed, mrRed) descending and keeps `pos` consistent.
    #[test]
    fn x1b_fm_a_heap_order() {
        let mut heap = FmHeap::default();
        heap.reset(6, 6);
        heap.push(0, 1, 5);
        heap.push(1, 3, 0);
        heap.push(2, 3, 2);
        heap.push(3, -1, 9);
        heap.push(4, 0, 1);
        heap.heapify();
        assert_eq!(heap.top(), 2);
        heap.update(4, 7, 0);
        heap.insert(5, 3, 1);
        heap.remove(3);
        assert!(!heap.contains(3));
        let mut order = Vec::new();
        while !heap.empty() {
            let col = heap.top();
            heap.remove(col);
            order.push(col);
        }
        assert_eq!(order, vec![4, 2, 5, 1, 0]);
        assert!((0..6).all(|col| !heap.contains(col)));
    }

    /// A new LP has no FM objective column.
    #[test]
    fn x1b_fm_a_highs_lp_default_fme() {
        assert_eq!(HighsLp::default().fme_obj_col_, -1);
    }

    /// Cernikov's rule counts distinct original rows (bound rows encoded per
    /// column) against the number of eliminated columns plus 2.
    #[test]
    fn x1b_fm_b_cernikov_redundant() {
        let mut originals: RowOriginals = BTreeMap::new();
        let mut merged = BTreeSet::new();
        // two original rows: 2 <= 0 + 2
        assert!(!cernikov_redundant(&mut merged, &originals, 0, 1, 7, 0));
        assert_eq!(merged.iter().copied().collect::<Vec<_>>(), vec![0, 1]);
        // bound rows of column 7: -(2*7+1) and -(2*7+2)
        assert!(!cernikov_redundant(
            &mut merged,
            &originals,
            K_UPPER_BOUND_ROW,
            K_LOWER_BOUND_ROW,
            7,
            0
        ));
        assert_eq!(merged.iter().copied().collect::<Vec<_>>(), vec![-16, -15]);
        // derived rows 10 and 11 from {0,1} and {1,2,3}: 4 originals > 1 + 2
        originals.insert(10, [0, 1].into_iter().collect());
        originals.insert(11, [1, 2, 3].into_iter().collect());
        assert!(cernikov_redundant(&mut merged, &originals, 10, 11, 7, 1));
        assert!(!cernikov_redundant(&mut merged, &originals, 10, 11, 7, 2));
    }

    /// The LP of Zhang, Ploskas and Sahinidis (TestPresolveRules.cpp lines 680-696,
    /// 1.16.0), stored column-wise.
    fn paper_lp(col_upper0: f64, row_upper: [f64; 3], cost: [f64; 4]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 3;
        lp.col_cost_ = cost.to_vec();
        lp.col_lower_ = vec![0.0; 4];
        lp.col_upper_ = vec![kHighsInf; 4];
        lp.col_upper_[0] = col_upper0;
        lp.row_lower_ = vec![-kHighsInf; 3];
        lp.row_upper_ = row_upper.to_vec();
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 1, 4, 6, 9];
        lp.a_matrix_.index_ = vec![0, 0, 1, 2, 1, 2, 0, 1, 2];
        lp.a_matrix_.value_ = vec![-1.0, 1.0, 2.0, 3.0, 1.0, -1.0, -1.0, 2.0, 3.0];
        lp
    }

    fn fm_options(obj_reformulation: bool) -> HighsOptions {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_presolve_row_append = true;
        options.mip_fm_presolve = true;
        options.presolve_fm_obj_reformulation = obj_reformulation;
        options
    }

    fn fm_presolve(lp: HighsLp, options: &HighsOptions) -> (HPresolve, HighsPostsolveStack) {
        let mut postsolve_stack = HighsPostsolveStack::default();
        postsolve_stack.initialize_index_maps(lp.num_row_, lp.num_col_);
        let mut timer = HighsTimer::default();
        timer.start();
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, options, options.presolve_reduction_limit, Some(timer)));
        (presolve, postsolve_stack)
    }

    /// The objective reformulation appends one continuous free column z with
    /// cost 1, zeroes the other costs and adds the row c^T x - z <= -offset.
    #[test]
    fn x1b_fm_c_reformulate_objective() {
        let mut lp = paper_lp(40.0, [-30.0, 50.0, 40.0], [1.0, 0.0, 3.0, 4.0]);
        lp.offset_ = 2.5;
        let (mut presolve, mut stack) = fm_presolve(lp, &fm_options(true));
        presolve.fm_reformulate_objective(&mut stack);
        let m = &presolve.model;
        assert_eq!(m.num_col_, 5);
        assert_eq!(m.fme_obj_col_, 4);
        assert_eq!(m.col_cost_, vec![0.0, 0.0, 0.0, 0.0, 1.0]);
        assert_eq!(m.offset_, 0.0);
        assert_eq!(m.col_lower_[4], -kHighsInf);
        assert_eq!(m.col_upper_[4], kHighsInf);
        assert_eq!(m.num_row_, 4);
        assert_eq!(m.row_upper_[3], -2.5);
        assert_eq!(m.row_lower_[3], -kHighsInf);
        // objective row entries: x0 (1), x2 (3), x3 (4), z (-1)
        let mut row: Vec<(HighsInt, f64)> = presolve
            .get_row_vector(3)
            .iter(&presolve)
            .map(|nz| (nz.index(), nz.value()))
            .collect();
        row.sort_by_key(|e| e.0);
        assert_eq!(row, vec![(0, 1.0), (2, 3.0), (3, 4.0), (4, -1.0)]);
        // a second call is a no-op
        presolve.fm_reformulate_objective(&mut stack);
        assert_eq!(presolve.model.num_col_, 5);
    }

    /// On the paper LP without costs, FM eliminates columns (direct call).
    #[test]
    fn x1b_fm_d_direct_eliminates() {
        let lp = paper_lp(40.0, [-30.0, 50.0, 40.0], [0.0; 4]);
        let (mut presolve, mut stack) = fm_presolve(lp, &fm_options(false));
        let mut num_cols_eliminated: HighsInt = 0;
        let result = presolve.fourier_motzkin(&mut stack, &mut num_cols_eliminated);
        assert_eq!(result, HPresolveResult::kOk);
        assert!(num_cols_eliminated > 0);
        assert_eq!(presolve.num_deleted_cols, num_cols_eliminated);
    }

    /// Solves `lp` with presolve `use_presolve` and `solver`, checking optimality and zero
    /// infeasibilities as solveAndCheck (TestPresolveRules.cpp lines 1596-1643) does.
    fn solve_and_check(
        lp: &HighsLp,
        obj_reformulation: bool,
        solver: &str,
        use_presolve: bool,
    ) -> (f64, HighsLp) {
        let mut h = Highs::default();
        h.options_.set_all_h16(true);
        h.options_.output_flag = false;
        h.options_.mip_presolve_row_append = true;
        h.options_.mip_fm_presolve = true;
        h.options_.presolve_fm_obj_reformulation = obj_reformulation;
        h.options_.presolve_rule_test = kH16PresolveRuleFourierMotzkin;
        h.set_option_value_str("presolve", if use_presolve { "on" } else { "off" });
        h.set_option_value_str(
            "run_crossover",
            if solver == kIpmString { "off" } else { "on" },
        );
        h.set_option_value_str("solver", solver);
        assert_eq!(h.pass_model(lp.clone()), HighsStatus::kOk);
        h.run();
        assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
        assert_eq!(h.get_info().num_primal_infeasibilities, 0);
        assert_eq!(h.get_info().num_dual_infeasibilities, 0);
        if use_presolve {
            // The presolved model size, as recorded in the run data
            let plp = h.presolve_.get_reduced_problem();
            let run_data = h.get_run_data();
            assert_eq!(run_data.presolved_model_num_col, plp.num_col_);
            assert_eq!(run_data.presolved_model_num_row, plp.num_row_);
            assert_eq!(run_data.presolved_model_num_nz, plp.a_matrix_.num_nz());
            // Ensure that any basis postsolve is correct (IPM runs without crossover)
            if solver != kIpmString {
                assert_eq!(run_data.num_simplex_iterations_after_postsolve, 0);
            }
        } else {
            // TestRunData.cpp lines 198-204: no presolve leaves the counts invalid
            let run_data = h.get_run_data();
            assert!(run_data.valid);
            assert_eq!(run_data.presolved_model_num_col, kHighsIllegalIntMeasure);
            assert_eq!(
                run_data.num_simplex_iterations_after_postsolve,
                kHighsIllegalIntMeasure
            );
        }
        (
            h.get_info().objective_function_value,
            h.presolve_.get_reduced_problem().clone(),
        )
    }

    /// presolveOffOn (TestPresolveRules.cpp lines 1645-1659) plus an objective comparison.
    fn presolve_off_on(lp: &HighsLp, obj_reformulation: bool) -> Vec<HighsLp> {
        let (obj_off, _) = solve_and_check(lp, obj_reformulation, kSimplexString, false);
        let mut presolved = Vec::new();
        for solver in [kSimplexString, kIpmString] {
            let (obj_on, plp) = solve_and_check(lp, obj_reformulation, solver, true);
            assert!(
                (obj_on - obj_off).abs() <= 1e-6 * (1.0 + obj_off.abs()),
                "{solver}: {obj_on} vs {obj_off}"
            );
            presolved.push(plp);
        }
        presolved
    }

    /// Port of TestPresolveRules.cpp "test-fourier-motzkin" (1.16.0 lines
    /// 656-715).
    #[test]
    fn x1b_fm_d_test_fourier_motzkin() {
        // FM example from paper
        presolve_off_on(&paper_lp(40.0, [-30.0, 50.0, 40.0], [0.0; 4]), false);
        // FM example from paper - tightened
        presolve_off_on(&paper_lp(5.0, [-30.0, 75.0, 50.0], [0.0; 4]), false);
        // FM example from paper - tightened and with costs: objective reformulation
        let presolved = presolve_off_on(
            &paper_lp(5.0, [-30.0, 75.0, 50.0], [1.0, 2.0, 3.0, 4.0]),
            true,
        );
        for plp in presolved {
            assert_eq!(plp.num_col_, 1);
            assert_eq!(plp.num_row_, 6);
            assert_eq!(plp.a_matrix_.start_[plp.num_col_ as usize], 6);
        }
    }

    /// With `mip_fm_presolve` off the rule test does nothing.
    #[test]
    fn x1b_fm_w_rule_test_off_is_noop() {
        let lp = paper_lp(40.0, [-30.0, 50.0, 40.0], [0.0; 4]);
        let mut options = fm_options(false);
        options.mip_fm_presolve = false;
        options.presolve_rule_test = kH16PresolveRuleFourierMotzkin;
        let (mut presolve, mut stack) = fm_presolve(lp, &options);
        assert_eq!(
            presolve.presolve_rule_test(&mut stack),
            HPresolveResult::kOk
        );
        assert_eq!(presolve.num_deleted_cols, 0);
        assert_eq!(presolve.model.num_row_, 3);
    }
}
