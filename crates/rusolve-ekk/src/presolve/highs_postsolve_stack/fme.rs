//! HiGHS source: presolve/HighsPostsolveStack.h (1.16.0 lines 62-103, 125-136, 321-336,
//! 613-728), presolve/HighsPostsolveStack.cpp (1.16.0 lines 1297-1370)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiGHS #3255 (`mip_fm_presolve`, default off): the postsolve records of
//! the Fourier-Motzkin presolve. A block of eliminations is pushed step by step
//! (`fourier_motzkin_block_push_step`, one call per eliminated column) and closed by
//! `fourier_motzkin_block_finalise`, which records the descendants of every parent row,
//! the new-row origins, the step headers and the step count, then the reduction. The
//! block undo is in `undo_fme.rs`.
//!
//! PORT NOTES:
//! * `FmeRowData<RowStorageFormat>` holds a `HighsMatrixSlice` view of the presolve
//!   matrix. Here `row_vec` is the row's nonzeros collected in slice order when the
//!   data is built; the only caller (HPresolve FM main loop) pushes the step right
//!   after collecting the rows, with no matrix change in between, so the entries are
//!   the ones the C++ slice yields.
//! * `fourierMotzkinBlockFinalise` iterates `std::unordered_map<HighsInt,
//!   std::vector<FmeAncestryEntry>> rowAncestry`; the port takes a `BTreeMap` (row
//!   order). The order only decides the order of the descendants of each parent row,
//!   which the dual undo sums.

use super::*;
use crate::prelude::*;

/// HiGHS 1.16: HighsPostsolveStack.h::FmeRowHeader (lines 62-66)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FmeRowHeader {
    pub row: HighsInt,
    pub row_lower: f64,
    pub row_upper: f64,
}
impl_highs_data_stack_item!(FmeRowHeader {
    row: HighsInt,
    row_lower: f64,
    row_upper: f64
});

/// HiGHS 1.16: HighsPostsolveStack.h::FmeRowData (lines 68-72); see the module doc.
#[derive(Clone, Debug, Default)]
pub struct FmeRowData {
    pub header: FmeRowHeader,
    pub row_vec: Vec<HighsSliceNonzero>,
}

/// HiGHS 1.16: HighsPostsolveStack.h::FmeStepHeader (lines 74-80)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FmeStepHeader {
    pub col_lower: f64,
    pub col_upper: f64,
    pub col: HighsInt,
    pub num_plus: HighsInt,
    pub num_minus: HighsInt,
}
impl_highs_data_stack_item!(FmeStepHeader {
    col_lower: f64,
    col_upper: f64,
    col: HighsInt,
    num_plus: HighsInt,
    num_minus: HighsInt
});

/// HiGHS 1.16: HighsPostsolveStack.h::FmeDescendant (lines 82-85)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FmeDescendant {
    pub row: HighsInt,
    pub scale_factor: f64,
}
impl_highs_data_stack_item!(FmeDescendant {
    row: HighsInt,
    scale_factor: f64
});

/// HiGHS 1.16: HighsPostsolveStack.h::FmeNewRow (lines 87-91)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FmeNewRow {
    pub row: HighsInt,
    pub plus_parent_idx: HighsInt,
    pub minus_parent_idx: HighsInt,
}
impl_highs_data_stack_item!(FmeNewRow {
    row: HighsInt,
    plus_parent_idx: HighsInt,
    minus_parent_idx: HighsInt
});

/// HiGHS 1.16: HighsPostsolveStack.h::FmeAncestryEntry (lines 93-98)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FmeAncestryEntry {
    pub step: HighsInt,
    pub parent_row_index: HighsInt,
    pub scale: f64,
    pub is_minus: bool,
}

/// HiGHS 1.16: HighsPostsolveStack.h::FmeBlockStep (lines 100-103)
#[derive(Clone, Debug, Default)]
pub struct FmeBlockStep {
    pub header: FmeStepHeader,
    pub new_rows: Vec<FmeNewRow>,
}

/// HiGHS 1.16: HighsPostsolveStack.h::FourierMotzkinObjCol (lines 125-136)
#[derive(Clone, Copy, Debug, Default)]
pub struct FourierMotzkinObjCol {
    pub offset: f64,
    pub col: HighsInt,
}
impl_highs_data_stack_item!(FourierMotzkinObjCol {
    offset: f64,
    col: HighsInt
});

/// HiGHS 1.16: HighsPostsolveStack.h::FmeParentRow (lines 321-326)
#[derive(Clone, Debug, Default)]
pub struct FmeParentRow {
    pub header: FmeRowHeader,
    pub coef: f64,
    pub entries: Vec<Nonzero>,
    pub descendants: Vec<FmeDescendant>,
}

/// HiGHS 1.16: HighsPostsolveStack.h::FmeStepData (lines 328-333)
#[derive(Clone, Debug, Default)]
pub struct FmeStepData {
    pub header: FmeStepHeader,
    pub plus_rows: Vec<FmeParentRow>,
    pub minus_rows: Vec<FmeParentRow>,
    pub new_rows: Vec<FmeNewRow>,
}

impl FourierMotzkinObjCol {
    /// HiGHS 1.16: HighsPostsolveStack.cpp::FourierMotzkinObjCol::transformToPresolvedSpace
    /// (lines 1297-1304). Not called by 1.16.0; kept for the port's
    /// completeness.
    pub fn transform_to_presolved_space(&self, cost_entries: &[Nonzero], primal_sol: &mut [f64]) {
        let mut val = self.offset;
        for entry in cost_entries {
            val += entry.value * primal_sol[entry.index as usize];
        }
        primal_sol[self.col as usize] = val;
    }

    /// HiGHS 1.16: HighsPostsolveStack.cpp::FourierMotzkinObjCol::undo (lines 1306-1313)
    pub fn undo(&self, cost_entries: &[Nonzero], solution: &mut HighsSolution) {
        if !solution.dual_valid {
            return;
        }
        let z_dual = solution.col_dual[self.col as usize];
        for entry in cost_entries {
            solution.col_dual[entry.index as usize] += entry.value * z_dual;
        }
        solution.col_dual[self.col as usize] = 0.0;
    }
}

impl HighsPostsolveStack {
    /// The plus or minus half of `fourierMotzkinBlockPushStep` (HighsPostsolveStack.h
    /// lines 617-635 / 637-655, identical code).
    fn fme_push_rows(&mut self, col: HighsInt, rows: &[FmeRowData]) {
        let mut headers: Vec<FmeRowHeader> = Vec::with_capacity(rows.len());
        let mut coefs: Vec<f64> = Vec::with_capacity(rows.len());
        for rd in rows {
            let mut translated: Vec<Nonzero> = Vec::new();
            let mut coef = 0.0;
            for nz in &rd.row_vec {
                if nz.index() == col {
                    coef = nz.value();
                } else {
                    translated.push(Nonzero::new(
                        self.orig_col_index[nz.index() as usize],
                        nz.value(),
                    ));
                }
            }
            self.reduction_values.push_vec(&translated);
            coefs.push(coef);
            headers.push(FmeRowHeader {
                row: self.orig_row_index[rd.header.row as usize],
                row_lower: rd.header.row_lower,
                row_upper: rd.header.row_upper,
            });
        }
        self.reduction_values.push_vec(&coefs);
        self.reduction_values.push_vec(&headers);
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::fourierMotzkinBlockPushStep (lines 612-657)
    pub fn fourier_motzkin_block_push_step(
        &mut self,
        col: HighsInt,
        plus_rows: &[FmeRowData],
        minus_rows: &[FmeRowData],
    ) {
        self.fme_push_rows(col, plus_rows);
        self.fme_push_rows(col, minus_rows);
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::fourierMotzkinBlockFinalise (lines 659-716);
    /// `rowAncestry` is a `BTreeMap` (module doc).
    pub fn fourier_motzkin_block_finalise(
        &mut self,
        block_steps: &[FmeBlockStep],
        row_ancestry: &BTreeMap<HighsInt, Vec<FmeAncestryEntry>>,
    ) {
        let num_steps = block_steps.len();

        let mut plus_descendants_all: Vec<Vec<Vec<FmeDescendant>>> = vec![Vec::new(); num_steps];
        let mut minus_descendants_all: Vec<Vec<Vec<FmeDescendant>>> = vec![Vec::new(); num_steps];
        for s in 0..num_steps {
            plus_descendants_all[s].resize(block_steps[s].header.num_plus as usize, Vec::new());
            minus_descendants_all[s].resize(block_steps[s].header.num_minus as usize, Vec::new());
        }
        for (&row, entries) in row_ancestry {
            let orig_row = self.orig_row_index[row as usize];
            for a in entries {
                let d = FmeDescendant {
                    row: orig_row,
                    scale_factor: a.scale,
                };
                if a.is_minus {
                    minus_descendants_all[a.step as usize][a.parent_row_index as usize].push(d);
                } else {
                    plus_descendants_all[a.step as usize][a.parent_row_index as usize].push(d);
                }
            }
        }

        for s in 0..num_steps {
            debug_assert_eq!(
                plus_descendants_all[s].len(),
                block_steps[s].header.num_plus as usize
            );
            for p in 0..block_steps[s].header.num_plus as usize {
                self.reduction_values.push_vec(&plus_descendants_all[s][p]);
            }
            debug_assert_eq!(
                minus_descendants_all[s].len(),
                block_steps[s].header.num_minus as usize
            );
            for m in 0..block_steps[s].header.num_minus as usize {
                self.reduction_values.push_vec(&minus_descendants_all[s][m]);
            }
        }

        for step in block_steps {
            let translated: Vec<FmeNewRow> = step
                .new_rows
                .iter()
                .map(|nr| FmeNewRow {
                    row: self.orig_row_index[nr.row as usize],
                    plus_parent_idx: nr.plus_parent_idx,
                    minus_parent_idx: nr.minus_parent_idx,
                })
                .collect();
            self.reduction_values.push_vec(&translated);
        }

        for step in block_steps {
            let header = FmeStepHeader {
                col_lower: step.header.col_lower,
                col_upper: step.header.col_upper,
                col: self.orig_col_index[step.header.col as usize],
                num_plus: step.header.num_plus,
                num_minus: step.header.num_minus,
            };
            self.reduction_values.push(&header);
        }

        self.reduction_values.push(&(num_steps as HighsInt));
        self.reduction_added(ReductionType::kFourierMotzkinBlock);
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::fourierMotzkinObjCol (lines 718-728)
    pub fn fourier_motzkin_obj_col(
        &mut self,
        col: HighsInt,
        offset: f64,
        cost_entries: &[Nonzero],
    ) {
        self.reduction_values.push(&FourierMotzkinObjCol {
            offset,
            col: self.orig_col_index[col as usize],
        });
        let translated_entries: Vec<Nonzero> = cost_entries
            .iter()
            .filter(|entry| entry.index != col)
            .map(|entry| Nonzero::new(self.orig_col_index[entry.index as usize], entry.value))
            .collect();
        self.reduction_values.push_vec(&translated_entries);
        self.reduction_added(ReductionType::kFourierMotzkinObjCol);
    }

    /// HiGHS 1.16: HighsPostsolveStack.cpp::popFourierMotzkinBlock (lines 1315-1370)
    pub fn pop_fourier_motzkin_block(stack: &mut HighsDataStack) -> Vec<FmeStepData> {
        let mut num_steps: HighsInt = 0;
        stack.pop(&mut num_steps);

        let mut steps: Vec<FmeStepData> = vec![FmeStepData::default(); num_steps as usize];

        // step headers
        for s in (0..num_steps as usize).rev() {
            stack.pop(&mut steps[s].header);
        }

        // new row origins
        for s in (0..num_steps as usize).rev() {
            stack.pop_vec(&mut steps[s].new_rows);
        }

        // descendants
        for s in (0..num_steps as usize).rev() {
            let num_minus = steps[s].header.num_minus as usize;
            steps[s]
                .minus_rows
                .resize(num_minus, FmeParentRow::default());
            for m in (0..num_minus).rev() {
                stack.pop_vec(&mut steps[s].minus_rows[m].descendants);
            }
            let num_plus = steps[s].header.num_plus as usize;
            steps[s].plus_rows.resize(num_plus, FmeParentRow::default());
            for p in (0..num_plus).rev() {
                stack.pop_vec(&mut steps[s].plus_rows[p].descendants);
            }
        }

        // row data
        for s in (0..num_steps as usize).rev() {
            // minus row data
            let mut minus_headers: Vec<FmeRowHeader> = Vec::new();
            let mut minus_coefs: Vec<f64> = Vec::new();
            stack.pop_vec(&mut minus_headers);
            stack.pop_vec(&mut minus_coefs);
            let num_minus = minus_coefs.len();
            for r in 0..num_minus {
                steps[s].minus_rows[r].header = minus_headers[r];
                steps[s].minus_rows[r].coef = minus_coefs[r];
            }
            for r in (0..num_minus).rev() {
                stack.pop_vec(&mut steps[s].minus_rows[r].entries);
            }

            // plus row data
            let mut plus_headers: Vec<FmeRowHeader> = Vec::new();
            let mut plus_coefs: Vec<f64> = Vec::new();
            stack.pop_vec(&mut plus_headers);
            stack.pop_vec(&mut plus_coefs);
            let num_plus = plus_coefs.len();
            for r in 0..num_plus {
                steps[s].plus_rows[r].header = plus_headers[r];
                steps[s].plus_rows[r].coef = plus_coefs[r];
            }
            for r in (0..num_plus).rev() {
                stack.pop_vec(&mut steps[s].plus_rows[r].entries);
            }
        }

        steps
    }
}
