//! HiGHS source: presolve/HighsPostsolveStack.cpp::undoFourierMotzkinBlock (1.16.0 lines
//! 1372-1673)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiGHS #3255 (`mip_fm_presolve`, default off): undo of one
//! Fourier-Motzkin block: primal values (Algorithm 3), row and column duals (Algorithm
//! 4), then a basis from the duals. The C++ lambdas are nested `fn`s / closures with
//! the captured state passed explicitly.

use super::*;
use crate::prelude::*;

impl HighsPostsolveStack {
    /// HiGHS 1.16: HighsPostsolveStack.cpp::undoFourierMotzkinBlock (lines 1372-1673)
    pub fn undo_fourier_motzkin_block(
        steps: &[FmeStepData],
        options: &HighsOptions,
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        let tol = options.primal_feasibility_tolerance;
        let dual_tol = options.dual_feasibility_tolerance;

        let num_steps = steps.len();

        // primal postsolve (Algorithm 3): process in reverse elimination order
        fn tighten_bounds(
            rows: &[FmeParentRow],
            col_value: &[f64],
            lower_bound: &mut f64,
            upper_bound: &mut f64,
        ) {
            for row in rows {
                let aij = row.coef;
                let mut sum = HighsCDouble::from(0.0);
                for nz in &row.entries {
                    sum += HighsCDouble::from(nz.value) * col_value[nz.index as usize];
                }
                let direction: HighsInt = if aij > 0.0 { 1 } else { -1 };
                let rhs_upper = if direction > 0 {
                    row.header.row_upper
                } else {
                    row.header.row_lower
                };
                let rhs_lower = if direction > 0 {
                    row.header.row_lower
                } else {
                    row.header.row_upper
                };
                if direction as f64 * rhs_upper != kHighsInf {
                    let bound = f64::from((rhs_upper - sum) / aij);
                    *upper_bound = upper_bound.min(bound);
                }
                if direction as f64 * rhs_lower != -kHighsInf {
                    let bound = f64::from((rhs_lower - sum) / aij);
                    *lower_bound = lower_bound.max(bound);
                }
            }
        }

        for s in (0..num_steps).rev() {
            let step = &steps[s];
            let col = step.header.col as usize;
            let mut lower = step.header.col_lower;
            let mut upper = step.header.col_upper;

            tighten_bounds(&step.plus_rows, &solution.col_value, &mut lower, &mut upper);
            tighten_bounds(
                &step.minus_rows,
                &solution.col_value,
                &mut lower,
                &mut upper,
            );

            if lower <= tol && upper >= -tol {
                solution.col_value[col] = 0.0;
            } else if lower > 0.0 {
                solution.col_value[col] = lower;
            } else {
                solution.col_value[col] = upper;
            }
        }

        if !solution.dual_valid {
            return;
        }

        // dual postsolve (Algorithm 4): process in reverse elimination order
        fn recover_dual(rows: &[FmeParentRow], row_dual: &mut [f64]) {
            for row in rows {
                let mut dual = HighsCDouble::from(0.0);
                for desc in &row.descendants {
                    dual += HighsCDouble::from(row_dual[desc.row as usize]) * desc.scale_factor;
                }
                row_dual[row.header.row as usize] += f64::from(dual);
            }
        }

        for s in (0..num_steps).rev() {
            let step = &steps[s];
            let col = step.header.col as usize;

            // u_i = sum_{k in K^j_i} lambda_k * scaleFactor
            recover_dual(&step.plus_rows, &mut solution.row_dual);
            recover_dual(&step.minus_rows, &mut solution.row_dual);

            // col_dual = -sum a_{ij} * row_dual[i] (cost is zero after reformulation)
            let mut col_dual = HighsCDouble::from(0.0);
            let mut visited: Vec<u8> = vec![0; solution.row_dual.len()];
            for row in &step.plus_rows {
                col_dual -=
                    HighsCDouble::from(row.coef) * solution.row_dual[row.header.row as usize];
                visited[row.header.row as usize] = 1;
            }
            for row in &step.minus_rows {
                if visited[row.header.row as usize] != 0 {
                    continue;
                }
                col_dual -=
                    HighsCDouble::from(row.coef) * solution.row_dual[row.header.row as usize];
            }
            solution.col_dual[col] = f64::from(col_dual);
        }

        // basis postsolve: use dual solution to determine basis status
        if !basis.valid {
            return;
        }

        // pre-compute lower and upper slacks for each row
        let compute_slacks = |solution: &HighsSolution,
                              col: usize,
                              rows: &[FmeParentRow],
                              lower_slacks: &mut Vec<f64>,
                              upper_slacks: &mut Vec<f64>| {
            let n = rows.len();
            lower_slacks.resize(n, 0.0);
            upper_slacks.resize(n, 0.0);
            for r in 0..n {
                let mut activity = HighsCDouble::from(rows[r].coef) * solution.col_value[col];
                for nz in &rows[r].entries {
                    activity +=
                        HighsCDouble::from(nz.value) * solution.col_value[nz.index as usize];
                }
                let act = f64::from(activity);
                lower_slacks[r] = if rows[r].header.row_lower != -kHighsInf {
                    act - rows[r].header.row_lower
                } else {
                    kHighsInf
                };
                upper_slacks[r] = if rows[r].header.row_upper != kHighsInf {
                    rows[r].header.row_upper - act
                } else {
                    kHighsInf
                };
            }
        };

        // row must be basic if it has zero dual and activity strictly
        // between bounds (complementary slackness)
        let row_must_be_basic =
            |solution: &HighsSolution, row: HighsInt, lower_slack: f64, upper_slack: f64| {
                solution.row_dual[row as usize].abs() <= dual_tol
                    && lower_slack > tol
                    && upper_slack > tol
            };

        // assign row as basic
        fn assign_basic_row_status(
            solution: &HighsSolution,
            basis: &mut HighsBasis,
            dual_tol: f64,
            row: HighsInt,
            basic_assigned: &mut HighsInt,
        ) -> bool {
            let r = row as usize;
            if basis.row_status[r] == HighsBasisStatus::kBasic
                || solution.row_dual[r].abs() > dual_tol
            {
                return false;
            }
            basis.row_status[r] = HighsBasisStatus::kBasic;
            *basic_assigned += 1;
            true
        }

        // assign row as non-basic
        fn assign_non_basic_row_status(
            solution: &HighsSolution,
            basis: &mut HighsBasis,
            dual_tol: f64,
            row: HighsInt,
            lower_slack: f64,
            upper_slack: f64,
        ) {
            let r = row as usize;
            if solution.row_dual[r] > dual_tol {
                basis.row_status[r] = HighsBasisStatus::kLower;
            } else if solution.row_dual[r] < -dual_tol {
                basis.row_status[r] = HighsBasisStatus::kUpper;
            } else {
                basis.row_status[r] = if upper_slack < lower_slack {
                    HighsBasisStatus::kUpper
                } else {
                    HighsBasisStatus::kLower
                };
            }
        }

        // assign row status
        #[allow(clippy::too_many_arguments)]
        fn assign_row_status(
            solution: &HighsSolution,
            basis: &mut HighsBasis,
            dual_tol: f64,
            row: HighsInt,
            lower_slack: f64,
            upper_slack: f64,
            basic_assigned: &mut HighsInt,
            force_non_basic: bool,
        ) {
            if force_non_basic
                || !assign_basic_row_status(solution, basis, dual_tol, row, basic_assigned)
            {
                assign_non_basic_row_status(
                    solution,
                    basis,
                    dual_tol,
                    row,
                    lower_slack,
                    upper_slack,
                );
            }
        }

        // collect a single candidate for basic assignment
        #[allow(clippy::too_many_arguments)]
        fn collect_candidate(
            solution: &HighsSolution,
            basis: &mut HighsBasis,
            dual_tol: f64,
            row: HighsInt,
            forced_non_basic: bool,
            lower_slack: f64,
            upper_slack: f64,
            entries: &[Nonzero],
            parent_index: HighsInt,
            is_minus: bool,
            candidates: &mut Vec<(HighsInt, HighsInt, bool)>,
        ) {
            if basis.row_status[row as usize] == HighsBasisStatus::kBasic {
                return;
            }
            if forced_non_basic || solution.row_dual[row as usize].abs() > dual_tol {
                assign_non_basic_row_status(
                    solution,
                    basis,
                    dual_tol,
                    row,
                    lower_slack,
                    upper_slack,
                );
                return;
            }
            let mut non_basic_count: HighsInt = 0;
            for nz in entries {
                if basis.col_status[nz.index as usize] != HighsBasisStatus::kBasic {
                    non_basic_count += 1;
                }
            }
            candidates.push((non_basic_count, parent_index, is_minus));
        }

        for s in (0..num_steps).rev() {
            let step = &steps[s];
            let col = step.header.col as usize;
            let num_plus = step.header.num_plus as usize;
            let num_minus = step.header.num_minus as usize;

            // compute slacks
            let mut plus_lower_slack: Vec<f64> = Vec::new();
            let mut plus_upper_slack: Vec<f64> = Vec::new();
            let mut minus_lower_slack: Vec<f64> = Vec::new();
            let mut minus_upper_slack: Vec<f64> = Vec::new();
            compute_slacks(
                solution,
                col,
                &step.plus_rows,
                &mut plus_lower_slack,
                &mut plus_upper_slack,
            );
            compute_slacks(
                solution,
                col,
                &step.minus_rows,
                &mut minus_lower_slack,
                &mut minus_upper_slack,
            );

            // non-basic propagation: if a generated row is non-basic (with nonzero
            // dual), both its parents are forced non-basic. mark them so the greedy
            // passes skip them. only force if the parent doesn't must-be-basic.
            let mut forced_non_basic_plus: Vec<bool> = vec![false; num_plus];
            let mut forced_non_basic_minus: Vec<bool> = vec![false; num_minus];
            for nr in &step.new_rows {
                // get indices of parent rows
                let p = nr.plus_parent_idx;
                let m = nr.minus_parent_idx;
                // skip basic rows (zero dual) and degenerate non-basic rows (with zero
                // dual)
                if p < 0 || m < 0 || solution.row_dual[nr.row as usize].abs() <= dual_tol {
                    continue;
                }
                let (p, m) = (p as usize, m as usize);
                // mark rows that do not have to be basic
                if !row_must_be_basic(
                    solution,
                    step.plus_rows[p].header.row,
                    plus_lower_slack[p],
                    plus_upper_slack[p],
                ) {
                    forced_non_basic_plus[p] = true;
                }
                if !row_must_be_basic(
                    solution,
                    step.minus_rows[m].header.row,
                    minus_lower_slack[m],
                    minus_upper_slack[m],
                ) {
                    forced_non_basic_minus[m] = true;
                }
            }

            // mark ranged rows (appearing in both plus and minus sets)
            let mut is_minus_row_ranged: Vec<bool> = vec![false; num_minus];
            let mut num_ranged: HighsInt = 0;
            for m in 0..num_minus {
                for p in 0..num_plus {
                    if step.minus_rows[m].header.row == step.plus_rows[p].header.row {
                        is_minus_row_ranged[m] = true;
                        num_ranged += 1;
                        break;
                    }
                }
            }

            // count number of basic new rows
            let num_new_rows = step.new_rows.len() as HighsInt;
            let mut num_basic_new_rows: HighsInt = 0;
            for nr in &step.new_rows {
                if basis.row_status[nr.row as usize] == HighsBasisStatus::kBasic {
                    num_basic_new_rows += 1;
                }
            }

            // how many basic variables are needed?
            let basic_needed =
                (num_plus as HighsInt + num_minus as HighsInt - num_ranged - num_new_rows)
                    + num_basic_new_rows;
            let mut basic_assigned: HighsInt = 0;

            // determine col status
            let col_must_be_basic = solution.col_value[col] > step.header.col_lower + tol
                && solution.col_value[col] < step.header.col_upper - tol;
            let col_can_be_basic = col_must_be_basic || solution.col_dual[col].abs() <= dual_tol;

            // pass 1: assign all must-be-basic (col and rows)
            if col_must_be_basic {
                basis.col_status[col] = HighsBasisStatus::kBasic;
                basic_assigned += 1;
            }
            for p in 0..num_plus {
                if forced_non_basic_plus[p] {
                    continue;
                }
                if row_must_be_basic(
                    solution,
                    step.plus_rows[p].header.row,
                    plus_lower_slack[p],
                    plus_upper_slack[p],
                ) {
                    assign_basic_row_status(
                        solution,
                        basis,
                        dual_tol,
                        step.plus_rows[p].header.row,
                        &mut basic_assigned,
                    );
                }
            }
            for m in 0..num_minus {
                if is_minus_row_ranged[m] || forced_non_basic_minus[m] {
                    continue;
                }
                if row_must_be_basic(
                    solution,
                    step.minus_rows[m].header.row,
                    minus_lower_slack[m],
                    minus_upper_slack[m],
                ) {
                    assign_basic_row_status(
                        solution,
                        basis,
                        dual_tol,
                        step.minus_rows[m].header.row,
                        &mut basic_assigned,
                    );
                }
            }

            // pass 2: assign can-be-basic col (if not already assigned)
            if !col_must_be_basic {
                if col_can_be_basic && basic_assigned < basic_needed {
                    basis.col_status[col] = HighsBasisStatus::kBasic;
                    basic_assigned += 1;
                } else if solution.col_value[col] <= step.header.col_lower + tol {
                    basis.col_status[col] = HighsBasisStatus::kLower;
                } else {
                    basis.col_status[col] = HighsBasisStatus::kUpper;
                }
            }

            // pass 3: assign can-be-basic rows, sorted by non-basic support count
            // to reduce risk of rank deficiency in degenerate cases
            let mut candidates: Vec<(HighsInt, HighsInt, bool)> = Vec::new();
            for p in 0..num_plus {
                collect_candidate(
                    solution,
                    basis,
                    dual_tol,
                    step.plus_rows[p].header.row,
                    forced_non_basic_plus[p],
                    plus_lower_slack[p],
                    plus_upper_slack[p],
                    &step.plus_rows[p].entries,
                    p as HighsInt,
                    false,
                    &mut candidates,
                );
            }
            for m in 0..num_minus {
                if is_minus_row_ranged[m] {
                    continue;
                }
                collect_candidate(
                    solution,
                    basis,
                    dual_tol,
                    step.minus_rows[m].header.row,
                    forced_non_basic_minus[m],
                    minus_lower_slack[m],
                    minus_upper_slack[m],
                    &step.minus_rows[m].entries,
                    m as HighsInt,
                    true,
                    &mut candidates,
                );
            }
            // sort descending by non-basic support count
            pdqsort_by(&mut candidates[..], |a, b| a.0 > b.0);
            for cand in &candidates {
                let parent_index = cand.1 as usize;
                let force = basic_assigned >= basic_needed;
                if cand.2 {
                    assign_row_status(
                        solution,
                        basis,
                        dual_tol,
                        step.minus_rows[parent_index].header.row,
                        minus_lower_slack[parent_index],
                        minus_upper_slack[parent_index],
                        &mut basic_assigned,
                        force,
                    );
                } else {
                    assign_row_status(
                        solution,
                        basis,
                        dual_tol,
                        step.plus_rows[parent_index].header.row,
                        plus_lower_slack[parent_index],
                        plus_upper_slack[parent_index],
                        &mut basic_assigned,
                        force,
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod x1b_fm_tests {
    use super::*;

    fn row(r: HighsInt, lower: f64, upper: f64, nz: &[(HighsInt, f64)]) -> FmeRowData {
        FmeRowData {
            header: FmeRowHeader {
                row: r,
                row_lower: lower,
                row_upper: upper,
            },
            row_vec: nz
                .iter()
                .enumerate()
                .map(|(p, &(i, v))| HighsSliceNonzero::new(i, v, p as HighsInt))
                .collect(),
        }
    }

    /// Eliminate x (col 0, bounds [0, 10]) from {x + y <= 2 (row 0), -x + y <= 0 (row 1)};
    /// the appended row 2 is their sum, 2y <= 2.
    fn eliminate_x() -> HighsPostsolveStack {
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(2, 2);
        stack.append_rows_to_model(1);
        let plus = [row(0, -kHighsInf, 2.0, &[(0, 1.0), (1, 1.0)])];
        let minus = [row(1, -kHighsInf, 0.0, &[(0, -1.0), (1, 1.0)])];
        stack.fourier_motzkin_block_push_step(0, &plus, &minus);
        let step = FmeBlockStep {
            header: FmeStepHeader {
                col_lower: 0.0,
                col_upper: 10.0,
                col: 0,
                num_plus: 1,
                num_minus: 1,
            },
            new_rows: vec![FmeNewRow {
                row: 2,
                plus_parent_idx: 0,
                minus_parent_idx: 0,
            }],
        };
        let mut ancestry: BTreeMap<HighsInt, Vec<FmeAncestryEntry>> = BTreeMap::new();
        ancestry.insert(
            2,
            vec![
                FmeAncestryEntry {
                    step: 0,
                    parent_row_index: 0,
                    scale: 1.0,
                    is_minus: false,
                },
                FmeAncestryEntry {
                    step: 0,
                    parent_row_index: 0,
                    scale: 1.0,
                    is_minus: true,
                },
            ],
        );
        stack.fourier_motzkin_block_finalise(&[step], &ancestry);
        stack
    }

    #[test]
    fn x1b_fm_e_push_pop_block() {
        let mut stack = eliminate_x();
        stack.reduction_values.reset_position();
        let steps = HighsPostsolveStack::pop_fourier_motzkin_block(&mut stack.reduction_values);
        assert_eq!(steps.len(), 1);
        let s = &steps[0];
        assert_eq!(s.header.col, 0);
        assert_eq!((s.header.col_lower, s.header.col_upper), (0.0, 10.0));
        assert_eq!((s.header.num_plus, s.header.num_minus), (1, 1));
        assert_eq!(s.plus_rows.len(), 1);
        assert_eq!(s.minus_rows.len(), 1);
        let p = &s.plus_rows[0];
        assert_eq!(p.header.row, 0);
        assert_eq!(p.header.row_upper, 2.0);
        assert_eq!(p.coef, 1.0);
        assert_eq!(p.entries.len(), 1);
        assert_eq!((p.entries[0].index, p.entries[0].value), (1, 1.0));
        assert_eq!(p.descendants.len(), 1);
        assert_eq!(
            (p.descendants[0].row, p.descendants[0].scale_factor),
            (2, 1.0)
        );
        let m = &s.minus_rows[0];
        assert_eq!(m.header.row, 1);
        assert_eq!(m.coef, -1.0);
        assert_eq!(m.descendants.len(), 1);
        assert_eq!(s.new_rows.len(), 1);
        assert_eq!(s.new_rows[0].row, 2);
    }

    #[test]
    fn x1b_fm_f_undo_primal_dual() {
        let mut stack = eliminate_x();
        // Presolved model: col y (orig 1) and the appended row 2.
        stack.compress_index_maps(&[-1, -1, 0], &[-1, 0]);
        let mut solution = HighsSolution::default();
        solution.value_valid = true;
        solution.dual_valid = true;
        solution.col_value = vec![1.0];
        solution.col_dual = vec![0.0];
        solution.row_value = vec![2.0];
        solution.row_dual = vec![-1.0];
        let mut basis = HighsBasis::default();
        basis.valid = true;
        basis.col_status = vec![HighsBasisStatus::kBasic];
        basis.row_status = vec![HighsBasisStatus::kUpper];
        let mut options = HighsOptions::default();
        options.mip_presolve_row_append = true;
        options.mip_fm_presolve = true;
        stack.undo(&options, &mut solution, &mut basis, -1, false);

        // x in [1, 1] from x <= 2 - y and x >= y.
        assert_eq!(solution.col_value, vec![1.0, 1.0]);
        // Parent duals inherit the descendant dual (scale 1); x has zero reduced cost.
        assert_eq!(solution.row_dual, vec![-1.0, -1.0]);
        assert_eq!(solution.col_dual[0], 0.0);
        // Back in the original 2x2 space; x basic, both parent rows nonbasic at upper.
        assert_eq!(basis.col_status[0], HighsBasisStatus::kBasic);
        assert_eq!(
            basis.row_status,
            vec![HighsBasisStatus::kUpper, HighsBasisStatus::kUpper]
        );
        let num_basic = basis
            .col_status
            .iter()
            .chain(basis.row_status.iter())
            .filter(|&&s| s == HighsBasisStatus::kBasic)
            .count();
        assert_eq!(num_basic, 2);
    }
}
